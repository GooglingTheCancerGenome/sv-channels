import pysam
import gzip
from pysam.libcalignedsegment import SAM_FLAGS
import numpy as np
import time
import pickle
import array

import sys
import os
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

    # single-read 1-position.
    # these do not include reads that are soft-clipped on both ends.
    SOFT_LEFT_FWD = enum.auto()
    SOFT_RIGHT_FWD = enum.auto()
    SOFT_LEFT_REV = enum.auto()
    SOFT_RIGHT_REV = enum.auto() # 4

    INS_FWD = enum.auto()
    INS_REV = enum.auto() # 6

    # lots of variants with high NM are usually bad alignments
    HIGH_NM = enum.auto()

    # soft-clipped on both ends usually noise
    SOFT_BOTH = enum.auto()

    # single-read 2-position
    DEL_FWD = enum.auto()
    DEL_REV = enum.auto() # 10

    # read-pair 2-position
    SPLIT_LOW_QUALITY       = enum.auto()
    SPLIT_INTER_CHROMOSOMAL = enum.auto()

    MATE_UNMAPPED = enum.auto() # 13

    # 2 position orphanable
    # NOTE: orphanable events must be at end of this enum
    SPLIT_PLUS_PLUS = enum.auto()
    SPLIT_MINUS_MINUS = enum.auto()
    SPLIT_MINUS_PLUS = enum.auto()
    SPLIT_PLUS_MINUS = enum.auto() # 17

    DISCORDANT_PLUS_MINUS = enum.auto() # 18
    DISCORDANT_PLUS_PLUS = enum.auto() # 19
    DISCORDANT_MINUS_MINUS = enum.auto() # 20
    DISCORDANT_MINUS_PLUS = enum.auto() # 21

    def __str__(self):
        return str(self.value) + "\t" + self.name

orphanable_events = {
    Event.SPLIT_PLUS_PLUS,
    Event.SPLIT_PLUS_MINUS,
    Event.SPLIT_MINUS_MINUS,
    Event.SPLIT_MINUS_PLUS,
    Event.DISCORDANT_PLUS_MINUS,
    Event.DISCORDANT_PLUS_PLUS,
    Event.DISCORDANT_MINUS_MINUS,
    Event.DISCORDANT_MINUS_PLUS,
}


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
    (rname, pos, strand,cigar,mapq, nm) = sa.split(',')
    sa_left = cigar[len(cigar) - 1] != 'S'

    mapq = int(mapq)
    if min(mapq, aln.mapping_quality) == 0: return
    pos = int(pos)

    interchromosomal = rname != aln.reference_name

    lookup = {('-', True): Event.SPLIT_MINUS_MINUS,
              ('+', True): Event.SPLIT_PLUS_MINUS,
              ('-', False): Event.SPLIT_MINUS_PLUS,
              ('+', False): Event.SPLIT_PLUS_PLUS,

              (True, '-'): Event.SPLIT_MINUS_MINUS,
              (True, '+'): Event.SPLIT_MINUS_PLUS,
              (False, '-'): Event.SPLIT_PLUS_MINUS,
              (False, '+'): Event.SPLIT_PLUS_PLUS,
           }

    if pos < aln.reference_start:
        if sa_left:
            li.append((rname, pos,
                       aln.reference_name, (aln.reference_start if self_left else aln.reference_end),
                       lookup[(strand, aln.is_reverse)], aln.qname))
        else:
            li.append((rname, sa_end(pos, cigar) - 1,
                       aln.reference_name, aln.reference_start if self_left else aln.reference_end - 1,
                       lookup[(strand, aln.is_reverse)], aln.qname))
    else:
        if sa_left:
            li.append((aln.reference_name, aln.reference_start - 1 if self_left else aln.reference_end,
                       rname, pos - 1,
                       lookup[(aln.is_reverse, strand)], aln.qname))
        else:
            li.append((aln.reference_name, aln.reference_start if self_left else aln.reference_end,
                       rname, sa_end(pos, cigar),
                       lookup[(aln.is_reverse, strand)], aln.qname))
    if 0 < min(mapq, aln.mapping_quality) < min_mapping_quality:
        t = list(li[-1])
        t[4] = Event.SPLIT_LOW_QUALITY
        li[-1] = tuple(t)
    elif interchromosomal:
        t = list(li[-1])
        t[4] = Event.SPLIT_INTER_CHROMOSOMAL
        li[-1] = tuple(t)


def proper_pair(aln, max_insert_size):
    if (aln.flag & SAM_FLAGS.FPROPER_PAIR) == 0: return False
    return max_insert_size is None or abs(aln.template_length) <= max_insert_size

def add_events(a, b, li, min_clip_len, min_cigar_event_length=10, min_mapping_quality=10, max_insert_size=None):

    for aln in (a, b):
        offset = aln.reference_start
        # we store if this read is left or right clipped here to avoid double
        # iteration over cigar.
        self_left = False
        self_right = False

        if aln.mapping_quality >= min_mapping_quality:
            chrom = aln.reference_name
            for i, (op, length) in enumerate(aln.cigartuples):
                if i == 0:
                    self_left = op == BAM_CSOFT_CLIP
                self_right = op == BAM_CSOFT_CLIP

                if op == BAM_CDEL and length >= min_cigar_event_length:
                    li.append((chrom, offset, chrom, offset + length,
                        Event.DEL_REV if aln.is_reverse else Event.DEL_FWD, aln.qname))
                if op in CONSUME_REF:
                    offset += length

        if aln.has_tag("SA"):
            sa_tag = aln.get_tag("SA")

            for sa in sa_tag.strip(';').split(";"):
                add_split_event(aln, sa, li, min_mapping_quality, self_left, self_right)
                #break # only add first split read event.

    if proper_pair(a, max_insert_size): return
    #if aln.is_proper_pair: print(abs(aln.template_length))

    if a.mapping_quality < min_mapping_quality and b.mapping_quality < min_mapping_quality: return
    if proper_pair(a, max_insert_size): return

    lookup = {(True, True): Event.DISCORDANT_MINUS_MINUS,
              (True, False): Event.DISCORDANT_MINUS_PLUS,
              (False, False): Event.DISCORDANT_PLUS_PLUS,
              (False, True): Event.DISCORDANT_PLUS_MINUS}

    if a.is_unmapped and not b.is_unmapped:
        # "right" and "left" don't make sense here, but "right" gives start and "left" gives end, which works.
        chrom = b.reference_name
        li.append((chrom, best_position(b, "right"), chrom, best_position(b, "left"), Event.MATE_UNMAPPED, a.qname))
        return
    if b.is_unmapped and not a.is_unmapped:
        chrom = a.reference_name
        li.append((chrom, best_position(a, "right"), chrom, best_position(a, "left"), Event.MATE_UNMAPPED, a.qname))
        return

    # we know that a < b because it came first in the bam
    li.append((a.reference_name, best_position(a, "left"), b.reference_name, best_position(b, "right"), lookup[(a.is_reverse, b.is_reverse)], a.qname))

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
        # convert back to numpy array for cumsum
        array = np.asarray(array, dtype='i4')
        np.cumsum(array, out=array)
        if len(array) > 10000000:
            print(f"[sv-channels] writing: {chrom}", file=sys.stderr)
        d = zarr.open(f"{outdir}/depths.{chrom}.bin", mode='w',
                      shape=array.shape, chunks=(10000,), dtype='i4')
        d[:] = array
        del d

def set_depth(aln, depths, chrom_lengths, outdir, min_mapq):
    if aln.mapping_quality < min_mapq: return
    if aln.flag & SAM_FLAGS.FUNMAP: return

    chrom = aln.reference_name
    if chrom not in depths:
        write_depths(depths, outdir)
        depths.clear()
        depths[chrom] = np.zeros((chrom_lengths[chrom],), dtype='i4')
        # faster to set individual values in array.array
        b = array.array('i')
        b.frombytes(depths[chrom].tobytes())
        depths[chrom] = b

    arr = depths[chrom]
    for s, e in aln.get_blocks():
        arr[s] += 1
        arr[e-1] -= 1

def soft_and_ins(aln, li, events2d, min_event_len, min_mapping_quality=15, high_nm=11):
    cigar = aln.cigartuples

    chrom = aln.reference_name

    if aln.mapping_quality >= min_mapping_quality:
        try:
            nm = aln.get_tag("NM")
            if nm >= high_nm:
                for op, length in cigar:
                    if op == BAM_CINS or op == BAM_CDEL:
                        nm -= length
                if nm >= high_nm:
                    li.append((chrom, int((aln.reference_start + aln.reference_end) / 2), Event.HIGH_NM, aln.qname))
        except KeyError:
            pass

    if cigar is None or len(cigar) < 2: return
    if cigar[0][0] == BAM_CSOFT_CLIP and cigar[-1][0] == BAM_CSOFT_CLIP:
        events2d.append((chrom, aln.reference_start, chrom, aln.reference_end, Event.SOFT_BOTH, aln.qname))

    if aln.mapping_quality < min_mapping_quality: return
    # check each end of read
    offset = aln.reference_start
    for op, length in cigar:
        if op == BAM_CSOFT_CLIP and length >= min_event_len:
            event = None
            if offset == aln.reference_start:
                event = Event.SOFT_LEFT_REV if aln.is_reverse else Event.SOFT_LEFT_FWD
            else:
                event = Event.SOFT_RIGHT_REV if aln.is_reverse else Event.SOFT_RIGHT_FWD
            li.append((chrom, offset, event, aln.qname))

        if op == BAM_CINS and length >= min_event_len:
            li.append((chrom, offset, Event.INS_REV if aln.is_reverse else Event.INS_FWD, aln.qname))

        if op in CONSUME_REF:
            offset += length

def chop(li, n, look_back):
    if look_back == 0: return
    for i, v in enumerate(li[-look_back:], start=max(0, len(li) - look_back)):
       if len(v) == n: continue
       assert li[i] == v
       li[i] = v[:n]

def iterate(bam, fai, outdir="sv-channels", min_clip_len=14,
        min_mapping_quality=10, max_insert_size=None, debug=False, high_nm=11):

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

    fail_flags = (SAM_FLAGS.FQCFAIL | SAM_FLAGS.FDUP | SAM_FLAGS.FSECONDARY | SAM_FLAGS.FSUPPLEMENTARY)

    # this is stupid, but helps performance
    # by default we store read name and other diagnostic info in softs and
    # events. if --debug is not used, then we use `chop` to remove that.
    # this tells chop how far back to look so it doesn't need to iterate over
    # the entire list.
    chop_mod = 400000

    # instead of incrementing numpy array for every cigar event.
    # add here, then when it gets large-enough we can do something faster.
    for i, b in enumerate(bam): # TODO add option to iterate over VCF of putative SV sites.

        if i == 2000000 or (i % 10000000 == 0 and i > 0):
            print(f"[sv-channels] i:{i} ({b.reference_name}:{b.reference_start}) processed-pairs:{processed_pairs} len pairs:{len(pairs)} events:{len(events)} reads/second:{i/(time.time() - t0):.0f}", file=sys.stderr)
            sys.stderr.flush()

        if not debug:
            # removing the debugging stuff, otherwise memory ballons.
            if len(softs) % chop_mod == 0:
                chop(softs, 3, chop_mod)
            if len(events) % chop_mod == 0:
                chop(events, 5, chop_mod)

        if b.flag & fail_flags: continue

        set_depth(b, depths, chrom_lengths, outdir, min_mapping_quality)

        if b.query_name in pairs:
            processed_pairs += 1
            a = pairs.pop(b.query_name)
            soft_and_ins(a, softs, events, min_clip_len, min_mapping_quality, high_nm)
            soft_and_ins(b, softs, events, min_clip_len, min_mapping_quality, high_nm)
            add_events(a, b, events, min_clip_len, min_cigar_event_length=10, max_insert_size=max_insert_size)
        else:
            # we dont use sequence or base-qualities so set them to empty to reduce memory.
            b.query_sequence = None
            b.query_qualities = None
            pairs[b.query_name] = b

    print(f"[sv-channels] processed-pairs:{processed_pairs} len pairs:{len(pairs)} events:{len(events)}", file=sys.stderr)
    write_depths(depths, outdir)
    if not debug:
        chop(softs, 3, chop_mod)
        chop(events, 5, chop_mod)
    write_text(softs, f"{outdir}/sv-channels.soft_and_insertions.txt.gz")
    write_text(events, f"{outdir}/sv-channels.events2d.txt.gz")

def find_max_insert_size(bam_path, reference, p=0.985):
    assert p < 1 and p > 0, "expected value between 0 and 1 in find_max_insert_size"
    reqd = SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_RNEXT | SAM_PNEXT
    bam = pysam.AlignmentFile(bam_path, "r", threads=2,
            format_options=[("required_fields=%d" % reqd).encode()],
            reference_filename=reference)
    fail_flags = (SAM_FLAGS.FREAD2 | SAM_FLAGS.FQCFAIL | SAM_FLAGS.FDUP | SAM_FLAGS.FSECONDARY | SAM_FLAGS.FSUPPLEMENTARY | SAM_FLAGS.FUNMAP | SAM_FLAGS.FMUNMAP)

    isizes = []
    skip = 400000
    n = 300000
    for aln in bam:
        flag = aln.flag
        if aln.mapping_quality < 30 or (flag & fail_flags): continue
        if aln.reference_id != aln.next_reference_id: continue
        if (flag & SAM_FLAGS.FREVERSE) == (flag & SAM_FLAGS.FMREVERSE): continue
        isizes.append(abs(aln.template_length))
        if len(isizes) - skip > n: break
    if len(isizes) - skip > n: isizes = isizes[skip:]
    isizes.sort()
    return isizes[max(0, int(len(isizes) * p - 1))]

def write_text(li, path, debug=False):
    li.sort()
    fh = gzip.open(path, mode="wt")
    if debug:
        for row in li:
            print("\t".join(str(x) for x in row), file=fh)
    else:
        for row in li:
            row = list(row)
            if isinstance(row[-1], str): # read-name
                row = row[:-1]
            row[-1] = row[-1].value
            print("\t".join(str(x) for x in row), file=fh)
    fh.close()

def main():
    p = argparse.ArgumentParser()
    p.add_argument("-o", "--out-dir", help="sample-specific output directory",
            default="sv-channels")
    p.add_argument("--min-mapping-quality", help="skip reads with mapping quality below this value (default=%(default)s)", type=int, default=10)
    p.add_argument("--debug", help="save debug information in output txt.gz files", action="store_true", default=False)
    p.add_argument("reference")
    p.add_argument("bam")

    a = p.parse_args()
    assert os.path.isfile(a.bam), "[svchannels] a file (or link) is required for the [bam] argument"

    reqd = SAM_QNAME | SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR | SAM_RNEXT | SAM_PNEXT | SAM_TLEN | SAM_AUX
    p = 0.99

    max_insert = find_max_insert_size(a.bam, a.reference, p)
    print(f"[sv-channels] calculated insert size of {p * 100:.2f}th percentile as: {max_insert}", file=sys.stderr)

    bam = pysam.AlignmentFile(a.bam, "r", threads=2,
            format_options=[("required_fields=%d" % reqd).encode()],
            reference_filename=a.reference)

    fai = Fasta(a.reference)

    iterate(bam, fai, a.out_dir, min_mapping_quality=a.min_mapping_quality, max_insert_size=max_insert, debug=a.debug, high_nm=10)
    print("[sv-channels] done with iteration", file=sys.stderr)


if __name__ == "__main__":

    main()
