import pysam
import gzip
from pysam.libcalignedsegment import SAM_FLAGS
import numpy as np
import time
import numba
import pickle

import sys
import argparse
import zarr
from pyfaidx import Fasta
import enum
import re

BAM_CMATCH, BAM_CINS, BAM_CDEL, BAM_CREF_SKIP, BAM_CSOFT_CLIP, BAM_CHARD_CLIP, BAM_CPAD, BAM_CEQUAL, BAM_CDIFF, BAM_CBACK = range(10)

SAM_QNAME = 0x00000001
SAM_FLAG  = 0x00000002
SAM_RNAME = 0x00000004
SAM_POS   = 0x00000008
SAM_MAPQ  = 0x00000010
SAM_CIGAR = 0x00000020
SAM_RNEXT = 0x00000040
SAM_PNEXT = 0x00000080
SAM_TLEN  = 0x00000100
SAM_SEQ   = 0x00000200
SAM_QUAL  = 0x00000400
SAM_AUX   = 0x00000800
SAM_RGAUX = 0x00001000

class Event(enum.IntEnum):

    # single-read 1-position
    SOFT_LEFT_FWD = enum.auto()
    SOFT_RIGHT_FWD = enum.auto()
    SOFT_LEFT_REV = enum.auto()
    SOFT_RIGHT_REV = enum.auto()

    INS_FWD = enum.auto()
    INS_REV = enum.auto()

    # single-read 2-position
    DEL_FWD = enum.auto()
    DEL_REV = enum.auto()
    # read-pair 2-position
    SPLIT_PLUS_PLUS = enum.auto()
    SPLIT_MINUS_MINUS = enum.auto()
    SPLIT_MINUS_PLUS = enum.auto()
    SPLIT_PLUS_MINUS = enum.auto()

    DISCORDANT_PLUS_MINUS = enum.auto()
    DISCORDANT_PLUS_PLUS = enum.auto()
    DISCORDANT_MINUS_MINUS = enum.auto()
    DISCORDANT_MINUS_PLUS = enum.auto()

    MATE_UNMAPPED = enum.auto()

    def __str__(self):
        return str(self.value)

CONSUME_REF = {BAM_CMATCH, BAM_CDEL, BAM_CREF_SKIP, BAM_CEQUAL, BAM_CDIFF}

def sa_end(pos, cigar, patt=re.compile("(\d+)(\w)")):
    # given a cigar as a string and a position. find the end
    # used to find end of supp alignment given in sa tag.
    result = pos
    for (l, op) in patt.findall(cigar):
        if op in "MDN=X":
            result += int(l)
    return result

def add_split_event(aln, sa, li, min_mapping_quality, self_left, self_right):
    # given an SA from aln, add the appropriate event to li.
    if aln.mapping_quality < min_mapping_quality: return
    (rname, pos, strand,cigar,mapq, nm) = sa.split(',')
    sa_left = cigar[len(cigar) - 1] != 'S'

    mapq = int(mapq)
    if mapq < min_mapping_quality: return
    pos = int(pos)

    lookup = {
              ('-', True): Event.SPLIT_MINUS_MINUS,
              ('+', True): Event.SPLIT_PLUS_MINUS,
              ('-', False): Event.SPLIT_MINUS_PLUS,
              ('+', False): Event.SPLIT_PLUS_PLUS,

              (True, '-'): Event.SPLIT_MINUS_MINUS,
              (True, '+'): Event.SPLIT_MINUS_PUS,
              (False, '-'): Event.SPLIT_PLUS_MINUS,
              (False, '+'): Event.SPLIT_PLUS_PLUS,
    }

    # NOTE!!! not handling interchromosomals for now.
    if rname != aln.reference_name: return
    
    # NOTE: current code treats this:
    #
    # primary:      supp:
    # 100M50S->     <- 100S50M
    # which is a *GOOD* DEL signal as Event.SPLIT_PLUS_MINUS
    #
    # and this:
    #
    # supp:        primary:
    # 50S100M->    <-50M100S
    # which is *NOT* a DEL signals as Event.SPLIT_PLUS_MINUS
    # the same (though positions are different)

    # split is to left of primary.
    if pos < aln.reference_start:
        if sa_left:
            li.append((rname, pos, 
                       aln.reference_name, (aln.reference_start if self_left else aln.reference_end),
                       lookup[(strand, aln.is_reverse)]))
        else:
            # this, with self_left should be most common
            li.append((rname, sa_end(pos, cigar),
                       aln.reference_name, aln.reference_start if self_left else aln.reference_end,
                       lookup[(strand, aln.is_reverse)]))
    else: # split is right of primary
        if sa_left:
            # this with self right should be most common
            li.append((aln.reference_name, aln.reference_start if self_left else aln.reference_end,
                       rname, pos,
                       lookup[(aln.is_reverse, strand)]))
        else:
            li.append((aln.reference_name, aln.reference_start if self_left else aln.reference_end,
                       rname, sa_end(pos, cigar),
                       lookup[(aln.is_reverse, strand)]))


def proper_pair(aln):
    # TODO: Luca can adjust this as needed.
    return aln.is_proper_pair

def add_events(a, b, li, min_clip_len, min_cigar_event_length=10, min_mapping_quality=10):

    for aln in (a, b):
        offset = aln.reference_start
        if aln.mapping_quality < min_mapping_quality: continue
        # we store if this read is left or right clipped here to avoid double
        # iteration over cigar.
        self_left = False
        self_right = False

        for op, length in aln.cigartuples:
            self_left = offset == aln.reference_start and op == BAM_CSOFT_CLIP
            self_right = op == BAM_CSOFT_CLIP

            if op == BAM_CDEL and length >= min_cigar_event_length:
                li.append((aln.reference_name, offset, aln.reference_name, offset + length,
                    Event.DEL_REV if aln.is_reverse else Event.DEL_FWD))
            if op in CONSUME_REF:
                offset += length
        sa_tag = None
        try:
            sa_tag = aln.get_tag("SA")
        except KeyError:
            continue

        for sa in sa_tag.strip(';').split(";"):
            add_split_event(aln, sa, li, min_mapping_quality, self_left, self_right)
            break # only add first split read event.

    if proper_pair(a): return

    if a.mapping_quality < min_mapping_quality and b.mapping_quality < min_mapping_quality: return
    # TODO: Luca, what if one read has low MQ?

    lookup = {(True, True): Event.DISCORDANT_MINUS_MINUS,
              (True, False): Event.DISCORDANT_MINUS_PLUS,
              (False, False): Event.DISCORDANT_PLUS_PLUS,
              (False, True): Event.DISCORDANT_PLUS_MINUS}

    if a.is_unmapped and not b.is_unmapped:
        li.append((b.reference_name, best_position(b, "left"), b.reference_name, best_position(b, "right"), Event.MATE_UNMAPPED))
        return
    if b.is_unmapped and not a.is_unmapped:
        li.append((a.reference_name, best_position(a, "left"), a.reference_name, best_position(a, "right"), Event.MATE_UNMAPPED))
        return


    # TODO: how to decide which positions to use for discordant?
    # we know that a < b because it came first in the bam
    li.append((a.reference_name, best_position(a, "left"), b.reference_name, best_position(b, "right"), lookup[(a.is_reverse, b.is_reverse)]))

def best_position(aln, side):
    cigar = aln.cigartuples
    if cigar is not None and len(cigar) > 1:
        if cigar[0][0] == BAM_CSOFT_CLIP:
            return aln.reference_start
        if cigar[-1][0] == BAM_CSOFT_CLIP:
            return aln.reference_end

    if side == "right": return aln.reference_start

    return aln.reference_end

def write_depths(depths, outdir):
    for chrom, array in depths.items():
        if len(array) > 10000000:
            print(f"[sv-channels] writing: {chrom}", file=sys.stderr)
        d = zarr.open(f"{outdir}/depths.{chrom}.bin", mode='w',
                      shape=array.shape, chunks=(10000,), dtype='i4')
        d[:] = array

@numba.jit(nopython=True)
def _set_depth(arr, s, e):
    # this is the fastest way I've found to do depth.
    arr[s:e] += 1

def set_depth(aln, depths, chrom_lengths, outdir):
    if aln.flag & (SAM_FLAGS.FUNMAP): return

    if aln.reference_name not in depths:
        write_depths(depths, outdir)
        depths.clear()
        depths[aln.reference_name] = np.zeros((chrom_lengths[aln.reference_name],), dtype='i4')

    arr = depths[aln.reference_name]
    for s, e in aln.get_blocks():
        _set_depth(arr, s, e)

def soft_and_ins(aln, li, min_event_len, min_mapping_quality=15):
    cigar = aln.cigartuples
    if cigar is None or len(cigar) < 2: return
    if aln.mapping_quality < min_mapping_quality: return
    # check each end of read
    offset = aln.reference_start
    for i, (op, length) in enumerate(cigar):
        if op == BAM_CSOFT_CLIP and length >= min_event_len:
            if i == 0: # softclip at left end of read
                event = Event.SOFT_LEFT_REV if aln.is_reverse else Event.SOFT_LEFT_FWD
            else: # softclip at right end of read
                event = Event.SOFT_RIGHT_REV if aln.is_reverse else Event.SOFT_RIGHT_FWD
            li.append((aln.reference_name, offset, event))

        if op == BAM_CINS and length >= min_event_len:
            li.append((aln.reference_name, offset, Event.INS_REV if aln.is_reverse else Event.INS_FWD))

        if op in CONSUME_REF:
            offset += length

def iterate(bam, fai, outdir="sv-channels", min_clip_len=16,
        min_mapping_quality=10):

    depths = {}

    t0 = time.time()

    # keep reads in here. only do anything with pairs. if a read is already in
    # here we have a pair since we are skipping supplementary and secondary.
    pairs = {} 

    chrom_lengths = dict(zip(bam.references, bam.lengths))

    # (chr, pos, Event)
    softs = []
    # (chrA, posA, chrB, posB, Event)
    events = []
    processed_pairs = 0

    # instead of incrementing numpy array for every cigar event.
    # add here, then when it gets large-enough we can do something faster.
    for i, b in enumerate(bam): # TODO add option to iterate over VCF of putative SV sites.

        if i == 2000000 or (i % 20000000 == 0 and i > 0):
            print(f"[sv-channels] i:{i} ({b.reference_name}:{b.reference_start}) processed-pairs:{processed_pairs} len pairs:{len(pairs)} events:{len(events)} reads/second:{i/(time.time() - t0):.0f}", file=sys.stderr)

        if b.flag & (SAM_FLAGS.FQCFAIL | SAM_FLAGS.FDUP
                | SAM_FLAGS.FSECONDARY | SAM_FLAGS.FSUPPLEMENTARY): continue 

        set_depth(b, depths, chrom_lengths, outdir)

        if b.query_name in pairs:
            processed_pairs += 1
            a = pairs.pop(b.query_name)
            soft_and_ins(a, softs, min_clip_len, min_mapping_quality)
            soft_and_ins(b, softs, min_clip_len, min_mapping_quality)
            add_events(a, b, events, min_clip_len)
        else:
            pairs[b.query_name] = b

    print(f"[sv-channels] processed-pairs:{processed_pairs} len pairs:{len(pairs)} events:{len(events)}", file=sys.stderr)
    write_depths(depths, outdir)

    return dict(soft_and_insertions=softs, events=events, depths=depths)

def write_text(li, path):
    li.sort()
    fh = gzip.open(path, mode="wt")
    for row in li:
        print("\t".join(str(x) for x in row), file=fh)
    fh.close()

def main():
    p = argparse.ArgumentParser()
    p.add_argument("-o", "--out-dir", help="sample-specific output directory",
            default="sv-channels")
    p.add_argument("reference")
    p.add_argument("bam")

    a = p.parse_args()

    reqd = SAM_QNAME | SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR | SAM_RNEXT | SAM_PNEXT | SAM_TLEN | SAM_AUX

    bam = pysam.AlignmentFile(a.bam, "r", threads=2,
            format_options=[("required_fields=%d" % reqd).encode()],
            reference_filename=a.reference)

    fai = Fasta(a.reference)

    d = iterate(bam, fai, a.out_dir)
    write_text(d["soft_and_insertions"], f"{a.out_dir}/sv-channels.soft_and_insertions.txt.gz")
    write_text(d["events"], f"{a.out_dir}/sv-channels.events2d.txt.gz")


if __name__ == "__main__":

    main()
