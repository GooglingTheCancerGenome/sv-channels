import sys
import os
import re
import time
from pathlib import Path
from pyfaidx import Fasta
from pysam import VariantFile

import numpy as np
np.set_printoptions(threshold=5000)

import zarr
from matplotlib import pyplot as plt

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.insert(0, os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from svchannels.extract_signals import Event, orphanable_events

DISCORDANTS = [Event.DISCORDANT_PLUS_MINUS, Event.DISCORDANT_PLUS_PLUS, Event.DISCORDANT_MINUS_MINUS, Event.DISCORDANT_MINUS_PLUS]


def breakpoints(bedpe):

    gt_lookup = {
            (1, 1): "1/1",
            (0, 1): "0/1",
            (0, 0): "0/0",
            (None, None): "./."
    }

    if bedpe.endswith((".vcf", ".vcf.gz", ".bcf")):
        for v in VariantFile(bedpe):
            svt = ""
            try:
                svt = v.info["SVTYPE"]
            except KeyError:
                pass
            samples = v.samples
            if len(samples) == 1:
                try:
                    svt += "\t" + gt_lookup[samples[0]['GT']]
                except:
                    svt += "\t./."
            yield (v.chrom, v.start, v.stop, svt)
    else:
        for line in xopen(bedpe):
            toks = line.strip().split()
            # if toks[6] != 'DEL': continue
            if toks[0] != toks[3]: continue

            chrom = toks[0]
            l, r = int(toks[1]), int(toks[4])
            yield (chrom, l, r, toks[6])

def find_signals_for_event(chrom, apos, bpos, signals2d_a, signals2d_b, signals1d, disc_a, disc_b, expand=250):
    # signals2d_a is sorted by A psoition, so we can use search-sorted
    # we can use that to find things unique to a as well as things shared with a and b.
    # but we need signals2d_b to find things unique to b.
    result = {}

    assert apos <= bpos, ("svchannels: expected apos < bpos, got:", apos, bpos)

    # only need to search apos since we propagate b, a to a, b before this...
    # but...
    aidxl = np.searchsorted(signals2d_a["a_pos"], apos - expand, side="left")
    aidxr = np.searchsorted(signals2d_a["a_pos"], apos + expand, side="right")

    subset = signals2d_a[aidxl: aidxr]

    # same chrom and with expected bounds on both ends.
    selection = (subset["b_pos"] >= (bpos - expand)) & (subset["b_pos"] <= (bpos + expand)) & (
            subset["a_chrom"] == chrom) & (subset["b_chrom"] == chrom)
    # slower
    #selection = subset.eval('(b_pos >= (@bpos - @expand)) & (b_pos <= (@bpos + @expand)) & (a_chrom == @chrom) & (b_chrom == @chrom)')
    result['shared'] = subset[selection]
    result['left-only'] = subset[~selection][["a_chrom", "a_pos", "event"]]
    result['left-only'].columns = ("chrom", "pos", "event")
    result['left-only'] = result['left-only'].to_dict(orient='records')

    # now repeat for B side.
    bidxl = np.searchsorted(signals2d_b["b_pos"], bpos - expand, side="left")
    bidxr = np.searchsorted(signals2d_b["b_pos"], bpos + expand, side="right")

    subsetb = signals2d_b[bidxl: bidxr]
    selectionb = (subsetb["a_pos"] >= apos - expand) & (subsetb["a_pos"] <= apos + expand) & (
            subsetb["a_chrom"] == chrom) & (subsetb["b_chrom"] == chrom)
    result['right-only'] = subsetb[~selectionb][["b_chrom", "b_pos", "event"]]
    result['right-only'].columns = ("chrom", "pos", "event")
    result['right-only'] = result['right-only'].to_dict(orient='records')

    # 1d search ...
    idxl = np.searchsorted(signals1d["pos"], apos - expand, side="left")
    idxr = np.searchsorted(signals1d["pos"], apos + expand, side="right")

    bidxl = np.searchsorted(signals1d["pos"], bpos - expand, side="left")
    bidxr = np.searchsorted(signals1d["pos"], bpos + expand, side="right")

    bidxl = max(bidxl, idxl)
    bidxr = max(bidxr, bidxl)

    result['left-only'].extend(signals1d[idxl:idxr].to_dict(orient='records'))
    result['right-only'].extend(signals1d[bidxl:bidxr].to_dict(orient='records'))

    # DISCORDANTS
    # now get discordants. we must divide by 8. everthing else is the same.
    apos >>= 3
    bpos >>= 3
    aidxl = np.searchsorted(disc_a["a_pos"], apos - expand, side="left")
    aidxr = np.searchsorted(disc_a["a_pos"], apos + expand, side="right")
    subset = disc_a[aidxl: aidxr]

    selection = (subset["b_pos"] >= (bpos - expand)) & (subset["b_pos"] <= (bpos + expand)) & (
            subset["a_chrom"] == chrom) & (subset["b_chrom"] == chrom)

    result['shared-discordants'] = subset[selection].copy(deep=True)

    lo = subset[~selection][["a_chrom", "a_pos", "event"]].copy(deep=True)
    lo.columns = ("chrom", "pos", "event")
    result['left-only-discordants'] = lo.to_dict(orient='records')

    # discordants b-size
    bidxl = np.searchsorted(disc_b["b_pos"], bpos - expand, side="left")
    bidxr = np.searchsorted(disc_b["b_pos"], bpos + expand, side="right")

    subsetb = disc_b[bidxl: bidxr]
    selectionb = (subsetb["a_pos"] >= apos - expand) & (subsetb["a_pos"] <= apos + expand) & (
            subsetb["a_chrom"] == chrom) & (subsetb["b_chrom"] == chrom)

    ro = subsetb[~selectionb][["b_chrom", "b_pos", "event"]].copy(deep=True)
    ro.columns = ("chrom", "pos", "event")
    result['right-only-discordants'] = ro.to_dict(orient='records')

    return result


def fill_arr(channels, events, posns, opos, offset, shift=0):
    pmax = channels.shape[1]

    for i in range(len(events)):
        # shift is by 3 for discordants to get back into array coords
        pos = posns[i] - (opos >> shift) + offset
        if pos < pmax and pos >= 0:
            channels[events[i], pos] += 1

def fill_orphan_dicts(channels, dicts, opos, offset, shift=0):
    pmax = channels.shape[1]
    m = max(Event) - 1
    for d in dicts:

        pos = d['pos'] - (opos >> shift) + offset
        if pos < pmax and pos >= 0:

            if Event(d['event']) in orphanable_events:
                channels[d['event'] + len(orphanable_events), pos] += 1
            else:
                channels[d['event'], pos] += 1

ONE_HOT = [np.array(l, dtype='c') for l in "ACTGN"]
N_onehot = len(ONE_HOT) # ACTGN

def generate_channels_for_event(chrom, apos, bpos, signals1d, signals2d_a, signals2d_b,
                                edisc_a, edisc_b,
                                expand, gap, depths, fasta):
    r = find_signals_for_event(chrom, apos, bpos, signals2d_a, signals2d_b, signals1d, edisc_a, edisc_b, expand=expand)

    ChannelShape = (max(Event) + 1 + len(orphanable_events) + N_onehot, 4 * expand + gap)
    channels = np.zeros(ChannelShape, dtype=np.int32)
    # TODO: handle apos - expand < 0
    #channels[0, 2*expand: 2 * expand + gap] = 0 # -1
    channels[0, 0:2 * expand] = depths[apos - expand:apos + expand]
    channels[0, 2 * expand + gap:] = depths[bpos - expand:bpos + expand]

    channels[-N_onehot:, 0: 2 * expand] = onehot(fasta, chrom, apos - expand, apos + expand)
    channels[-N_onehot:, 2 * expand + gap:] = onehot(fasta, chrom, bpos - expand, bpos + expand)

    fill_arr(channels, np.asarray(r["shared"]["event"]), np.asarray(r["shared"]["a_pos"]), apos, expand)
    fill_arr(channels, np.asarray(r["shared"]["event"]), np.asarray(r["shared"]["b_pos"]), bpos, 3 * expand + gap)
    fill_arr(channels, np.asarray(r["shared-discordants"]["event"]), np.asarray(r["shared-discordants"]["a_pos"]), apos, expand, shift=3)
    fill_arr(channels, np.asarray(r["shared-discordants"]["event"]), np.asarray(r["shared-discordants"]["b_pos"]), bpos, 3 * expand + gap, shift=3)

    if len(r["left-only"]) > 0:
        fill_orphan_dicts(channels, r['left-only'], apos, expand)
    if len(r["left-only-discordants"]) > 0:
        fill_orphan_dicts(channels, r['left-only-discordants'], apos, expand, shift=3)

    if len(r["right-only"]) > 0:
        fill_orphan_dicts(channels, r['right-only'], bpos, 3 * expand + gap)
    if len(r["right-only-discordants"]) > 0:
        fill_orphan_dicts(channels, r['right-only-discordants'], bpos, 3 * expand + gap, shift=3)

    """
    for i, row in enumerate(channels):
        if i == 0:
            print("depth:", row)
        else:
            if i > max(Event):
                print("ORPHAN:", Event(i-len(orphanable_events)), row)
            else:
                print(Event(i), row)
    """
    return channels

def xopen(filepath):
    import io
    with open(filepath, 'rb') as test_f:
        if test_f.read(2) == b'\x1f\x8b':
            import gzip
            return io.TextIOWrapper(gzip.open(filepath, 'rb'))
        else:
            return io.TextIOWrapper(open(filepath, 'rb'))


def plot_event(chan, event, expand, gap):
    chrom, l, r = event
    depth = chan[0]
    mat = chan[1:].astype('float')
    counts = mat.sum(axis=1).astype(int)
    ys, xs = np.where(mat > 0)
    vals = mat[ys, xs]
    fig, axes = plt.subplots(2, 1, figsize=(12, 12), sharex=True, gridspec_kw={'height_ratios': [1, 5]})
    axes[1].scatter(xs, ys, c=vals, s = 36, marker="s")
    axes[1].set_xlabel("relative position")
    axes[1].set_ylabel("channel")
    axes[1].set_yticks(range(max(Event) + len(orphanable_events) + len(ONE_HOT)))
    axes[1].set_yticklabels(
               [f"{Event(i + 1).name} ({counts[i]})" for i in range(max(Event))] +
               [f"orphan: {Event(i).name} ({counts[i + len(orphanable_events) - 1]})" for i in range(max(Event) + 1 - len(orphanable_events), max(Event) + 1) ] +
               [o for o in ONE_HOT])
    axes[1].set_ylim(-1, max(Event) + len(orphanable_events) + len(ONE_HOT))
    axes[1].set_xlim(0, len(depth))
    axes[0].plot(depth)
    axes[0].axvline(x=expand, zorder=-1, ls='--', c='gray')
    axes[0].axvline(x=3*expand+gap, zorder=-1, ls='--', c='gray')
    axes[0].axvspan(2*expand, 2*expand + gap, color='lightgray', edgecolor=None)
    axes[1].axvspan(2*expand, 2*expand + gap, color='lightgray', edgecolor=None)
    axes[0].set_ylabel("depth")

    fig.suptitle(f"{chrom}:{l}-{r}")

    plt.tight_layout()
    plt.show()


def onehot(fa, chrom, start, end, ATGC=ONE_HOT):
    fr = np.array(fa.get_seq(chrom, start + 1, end), dtype='c') # get_seq expects 1-based coords
    oh = np.zeros((len(ATGC), end - start), dtype=np.int8)
    for i, l in enumerate(ATGC):
        oh[i] = (fr == l)
    return oh


def main(args=sys.argv[1:]):
    import argparse
    import pandas as pd
    p = argparse.ArgumentParser()
    p.add_argument("directory", help="sample-specific directory created by svchannels extract for input")
    p.add_argument("output", help="sample-specific output directory")
    p.add_argument("bedpe", help="bedpe or VCF of structural variants")
    p.add_argument(
        '--expand',
        type=int,
        default=62,
        help="Specify expansion around breakpoints (default: %(default)s) for creating the channel. This will automatically extend 8x the specified value for discordant reads")
    p.add_argument(
        '--gap',
        type=int,
        default=8,
        help="Specify width of gap (default: %(default)s)")
    p.add_argument("--reference", help="reference fasta file", required=True)
    p.add_argument("--shift",
                   type=int,
                   default=10,
                   help="number of base pairs to shift the windows outwards and inwards",
                   required=True)

    a = p.parse_args(args)

    od = Path(a.output)
    for parent in reversed(od.parents):
        parent.mkdir(mode=0o774, exist_ok=True)
    od.mkdir(mode=0o774, exist_ok=True)

    depths_by_chrom = zarr.open_group(os.path.join(a.directory, "depths.zarr.zip"))

    # we need positions sorted by a and by b so we do that here.
    e2d_a = pd.read_table(f"{a.directory}/sv-channels.events2d.txt.gz", compression="gzip",
                        usecols=list(range(5)),
                        names=["a_chrom", "a_pos", "b_chrom", "b_pos", "event"],
                        dtype={"a_chrom": str, "b_chrom": str, "a_pos": np.int32, "b_pos": np.int32})

    edisc_a = e2d_a.iloc[e2d_a.event.isin(DISCORDANTS).values].copy(deep=True)
    # we compress discordant position by 8
    edisc_a.a_pos = edisc_a.a_pos.values >> 3
    edisc_a.b_pos = edisc_a.b_pos.values >> 3

    e2d_a = e2d_a.iloc[~e2d_a.event.isin(DISCORDANTS).values]

    e2d_a.sort_values("a_pos", inplace=True)
    edisc_a.sort_values("a_pos", inplace=True)

    e2d_b = e2d_a.copy(deep=True)
    e2d_b.sort_values("b_pos", inplace=True)

    edisc_b = edisc_a.copy(deep=True)
    edisc_b.sort_values("b_pos", inplace=True)

    fasta = Fasta(a.reference, as_raw=True)

    e1d = pd.read_table(f"{a.directory}/sv-channels.soft_and_insertions.txt.gz", compression="gzip",
                        usecols=list(range(3)),
                        names=["chrom", "pos", "event"],
                        dtype={"chrom": str, "pos": np.int32, "event": np.int8})
    print(f"[svchannels] read {len(e2d_a)} 2D events and {len(e1d)} 1d events", file=sys.stderr)
    t0 = time.time()
    expand = a.expand
    gap = a.gap
    shift = a.shift

    # write the lines considered to file
    fh_svp = open(os.path.join(a.output, 'sv_positions.bedpe'), 'w')

    dups_toks = set()
    events = [] # we iterate over events first so we can size the zarr array.

    for chrom, l, r, svt in breakpoints(a.bedpe):

        if l < expand: continue
        if r < expand: continue
        if r < l:
            l, r = r, l

        clen = len(depths_by_chrom[chrom])
        if r >= clen - expand: continue
        # here, sv_chan is shape (n_channels, 2 * expand + gap) can accumulate these and send to learner.
        toks_id = f'{chrom}_{l}_{r}'
        if toks_id in dups_toks: continue
        dups_toks.add(toks_id)
        # no shift
        fh_svp.write(f'{chrom}\t{l}\t{l+1}\t{chrom}\t{r}\t{r+1}\t{svt}\n')
        events.append((chrom, l, r))
        # shift outwards
        shift_out1 = np.random.randint(low=0, high=shift + 1)
        shift_out2 = np.random.randint(low=0, high=shift + 1)
        fh_svp.write(f'{chrom}\t{l - shift_out1}\t{l - shift_out2 + 1}\t{chrom}\t{r + shift_out1}\t{r + shift_out2 + 1}\t{svt}\n')
        events.append((chrom, l - shift_out1, r + shift_out2))
        # shift inwards
        shift_in1 = np.random.randint(low=0, high=shift + 1)
        shift_in2 = np.random.randint(low=0, high=shift + 1)
        fh_svp.write(f'{chrom}\t{l + shift_in1}\t{l + shift_in2 + 1}\t{chrom}\t{r - shift_in1}\t{r - shift_in2 + 1}\t{svt}\n')
        events.append((chrom, l + shift_in1, r - shift_in2))

    fh_svp.close()

    print(f"[svchannels] creating array", file=sys.stderr)
    ChannelShape = (len(events), max(Event) + 1 + len(orphanable_events) + N_onehot, 4 * expand + gap)
    #Z = zarr.open(a.directory + 'sv_chan.zarr', mode='w', shape=ChannelShape,
    #        chunks=(16, ChannelShape[1], ChannelShape[2]))
    #chunks=(16, ChannelShape[1], ChannelShape[2])
    chunks=(1, ChannelShape[1], ChannelShape[2])
    Z = zarr.zeros(ChannelShape, chunks=chunks, dtype='i4')
    print(f"[svchannels] shape of sv_chan array will be {Z.shape}", file=sys.stderr)
    t1 = time.time()
    n_vars = 10000
    for i, event in enumerate(events):

        ch = generate_channels_for_event(event[0], event[1], event[2], e1d, e2d_a, e2d_b,
                                             edisc_a, edisc_b,
                                             expand, gap, depths_by_chrom[event[0]], fasta)
        Z[i] = ch #.reshape((1, Z.shape[1], Z.shape[2]))
        if i % n_vars == 0 and i > 0:
            td = time.time() - t1
            print(f"{i}/{Z.shape}. in {td:.0f} seconds ({n_vars/td:.1f} SVs/second \n{Z.info}", file=sys.stderr)
            sys.stderr.flush()
            t1 = time.time()
        #plot_event(ch, event, expand, gap)

    zarr.save_array(os.path.join(a.output, 'channels.zarr.zip'), Z)
    t = time.time() - t0
    print(f"[svchannels] generated channels for {len(events)} SVs in {t:.1f} seconds ({len(events)/t:.0f} SVs/second)", file=sys.stderr)


if __name__ == "__main__":
    main()
