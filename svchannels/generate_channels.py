import sys
import os
import re
import time
from pyfaidx import Fasta


import numpy as np
np.set_printoptions(threshold=5000)

import zarr

from matplotlib import pyplot as plt

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.insert(0, os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from svchannels.extract_signals import Event, orphanable_events


def find_signals_for_event(chrom, apos, bpos, signals2d_a, signals2d_b, signals1d, expand=250):
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
    return result


def fill_arr(channels, events, posns, opos, offset):
    pmax = channels.shape[1]

    for i in range(len(events)):
        pos = posns[i] - opos + offset
        if pos < pmax and pos >= 0:
            channels[events[i], pos] += 1

def fill_orphan_dicts(channels, dicts, opos, offset):
    pmax = channels.shape[1]
    m = max(Event) - 1
    for d in dicts:

        pos = d['pos'] - opos + offset
        if pos < pmax and pos >= 0:

            if Event(d['event']) in orphanable_events:
                channels[d['event'] + len(orphanable_events), pos] += 1
            else:
                channels[d['event'], pos] += 1

ONE_HOT = [np.array(l, dtype='c') for l in "ACTGN"]
N_onehot = len(ONE_HOT) # ACTGN

def generate_channels_for_event(chrom, apos, bpos, signals1d, signals2d_a, signals2d_b, expand, gap, depths, fasta):
    r = find_signals_for_event(chrom, apos, bpos, signals2d_a, signals2d_b, signals1d, expand=expand)

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

    if len(r["left-only"]) > 0:
        fill_orphan_dicts(channels, r['left-only'], apos, expand)
    if len(r["right-only"]) > 0:
        # TODO this is broken from find_signals (not getting because not sorted on b
        fill_orphan_dicts(channels, r['right-only'], bpos, 3 * expand + gap)

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


def read_bins(directory):
    import glob
    bins = glob.glob(directory + "/depths*.bin")
    bin_patt = f"^{directory.rstrip('/')}\/depths\.(.+)\.bin$"
    bin_patt = re.compile(bin_patt)
    depths_by_chrom = {}
    for bin in bins:
        try:
            chrom = bin_patt.match(bin).groups(1)[0]
            depths_by_chrom[chrom] = zarr.open(bin, mode='r')
        except AttributeError:
            raise ValueError(f"unable to find chrom in {bin}")
    return depths_by_chrom


def plot_event(chan, toks, expand, gap):
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

    fig.suptitle(f"{toks[0]}:{toks[1]}-{toks[4]}")

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
    p.add_argument("directory", help="sample-specific directory created by svchannels extract")
    p.add_argument("bedpe")
    p.add_argument(
        '--expand',
        type=int,
        default=500,
        help="Specify width of single windows (default: %(default)s)")
    p.add_argument("--prefix", help="output prefix", default="")
    p.add_argument(
        '--gap',
        type=int,
        default=10,
        help="Specify width of gap (default: %(default)s)")
    p.add_argument("--reference", help="reference fasta file", required=True)

    a = p.parse_args(args)

    depths_by_chrom = read_bins(a.directory)
    e2d_a = pd.read_table(f"{a.directory}/sv-channels.events2d.txt.gz", compression="gzip",
                        usecols=list(range(5)),
                        names=["a_chrom", "a_pos", "b_chrom", "b_pos", "event"],
                        dtype={"a_chrom": str, "b_chrom": str})
    e2d_a.sort_values("a_pos", inplace=True)
    e2d_b = e2d_a.copy(deep=True)
    e2d_b.sort_values("b_pos", inplace=True)

    fasta = Fasta(a.reference, as_raw=True)

    e1d = pd.read_table(f"{a.directory}/sv-channels.soft_and_insertions.txt.gz", compression="gzip",
                        usecols=list(range(3)),
                        names=["chrom", "pos", "event"],
                        dtype={"chrom": str})
    print(f"[svchannels] read {len(e2d_a)} 2D events and {len(e1d)} 1d events", file=sys.stderr)
    t0 = time.time()
    expand = a.expand
    gap = a.gap


    # write the lines considered to file
    file_object = open(a.prefix + 'sv_positions.bedpe', 'w')

    dups_toks = set()
    events = [] # we iterate over events first so we can size the zarr array.

    for line in xopen(a.bedpe):

        toks = line.strip().split()
        # if toks[6] != 'DEL': continue
        if toks[0] != toks[3]: continue
        if int(toks[1]) < expand: continue
        if int(toks[4]) < expand: continue
        if int(toks[4]) < int(toks[1]):
            toks[4], toks[1] = toks[1], toks[4]

        clen = len(depths_by_chrom[toks[0]])
        if int(toks[4]) >= clen - expand: continue
        # here, sv_chan is shape (n_channels, 2 * expand + gap) can accumulate these and send to learner.
        toks_id = str(int(toks[1])) + '_' + str(int(toks[4]))
        if toks_id in dups_toks: continue
        dups_toks.add(toks_id)
        file_object.write(line)

        events.append((toks[0], int(toks[1]), int(toks[4])))

    file_object.close()

    print(f"[svchannels] creating array", file=sys.stderr)
    ChannelShape = (0, max(Event) + 1 + len(orphanable_events) + N_onehot, 4 * expand + gap)
    Z = zarr.open(a.prefix + 'sv_chan.zarr', mode='w', shape=ChannelShape,
            chunks=(1, 1, 4 * expand + gap))
    print(f"shape of sv_chan array will be {(len(events), Z.shape[0], Z.shape[1])}", file=sys.stderr)
    for i, event in enumerate(events):

        ch = generate_channels_for_event(event[0], event[1], event[2], e1d, e2d_a, e2d_b,
                                             expand, gap, depths_by_chrom[event[0]], fasta)
        Z.append(ch.reshape((1, Z.shape[1], Z.shape[2])))
        if i % 100 == 0:
            print(f"{i}/{Z.shape}", file=sys.stderr)
            sys.stderr.flush()

    t = time.time() - t0
    print(f"generated channels for {len(events)} SVs in {t:.1f} seconds ({n/t:.0f} SVs/second)", file=sys.stderr)


if __name__ == "__main__":
    import pandas as pd

    if False:
        outdir = "sv-channels"
        chrom = "1"
        depths = zarr.open(f"{outdir}/depths.{chrom}.bin", mode='r')

        signals2d = np.array([
            ("chr1", 11100, "chr3", 11400, Event(1)),
            ("chr1", 11200, "chr1", 11300, Event(2)),
            ("chr1", 11200, "chr1", 11400, Event(3)),
            ("chr1", 11200, "chr2", 11400, Event(4)),
            ("chr1", 11400, "chr1", 11200, Event(5)),
            ("chr1", 11500, "chr1", 11900, Event(6))
        ], dtype=[('a_chrom', 'S8'), ('a_pos', np.int32),
                  ('b_chrom', 'S8'), ('b_pos', np.int32),
                  ('event', 'i2')
                  ])
        signals2d = pd.DataFrame(signals2d)

        signals1d = np.array([
            ("chr1", 11100, Event(3)),
            ("chr1", 11400, Event(3)),
        ], dtype=[('chrom', 'S8'), ('pos', np.int32), ('event', 'i2')])
        signals1d = pd.DataFrame(signals1d)

        apos, bpos = (11100, 11400)
        print(f"querying {apos}, {bpos}")

        # print("all")
        # print(signals2d.tolist())

        expand = 100
        gap = 10

        # r = find_signals_for_event("chr1", 100, 400, signals2d, signals1d, expand=expand)
        # for k, v in r.items():
        #    print(k)
        #    print(v.tolist())

        print(generate_channels_for_event("chr1", apos, bpos, signals1d, signals2d, expand, gap, depths))

    main()
