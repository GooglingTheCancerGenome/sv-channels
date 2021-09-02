import sys
import os
import re
import time

import numpy as np
import numba
import zarr

PACKAGE_PARENT = '..'
SCRIPT_DIR = os.path.dirname(os.path.realpath(os.path.join(os.getcwd(), os.path.expanduser(__file__))))
sys.path.insert(0, os.path.normpath(os.path.join(SCRIPT_DIR, PACKAGE_PARENT)))

from svchannels.extract_signals import Event

def find_signals_for_event(apos, bpos, signals2d, signals1d, expand=250):
    # note: only works on same chromosome as implemented.
    result = {}

    assert apos <= bpos, ("svchannels: expected apos < bpos, got:", apos, bpos)

    # only need to search apos since we propagate b, a to a, b before this...
    # but...
    aidxl = np.searchsorted(signals2d["a_pos"], apos - expand, side="left")
    aidxr = np.searchsorted(signals2d["a_pos"], apos + expand, side="right")

    # need b-side for events like splits from b that go to another chrom.
    bidxl = np.searchsorted(signals2d["b_pos"], bpos - expand, side="left")
    bidxr = np.searchsorted(signals2d["b_pos"], bpos + expand, side="right")

    # make sure we don't grab same values more than once.
    bidxl = max(bidxl, aidxr)
    bidxr = max(bidxr, bidxl)

    subset = signals2d[aidxl: aidxr]

    # same chrom and with expected bounds on both ends.
    selection = (subset["b_pos"] >= (bpos - expand)) & (subset["b_pos"] <= (bpos + expand)) & (subset["a_chrom"] == subset["b_chrom"])
    result['shared'] = subset[selection]
    result['left-only'] = subset[~selection][["a_chrom", "a_pos", "event"]] #.tolist()

    subsetb = signals2d[bidxl: bidxr]
    selectionb = (subsetb["a_pos"] >= apos - expand) & (subsetb["a_pos"] <= apos + expand) & (subsetb["a_chrom"] == subsetb["b_chrom"])
    result['right-only'] =  subsetb[~selectionb][["a_chrom", "a_pos", "event"]] #.tolist()

    
    # 1d search ...
    idxl = np.searchsorted(signals1d["pos"], apos - expand, side="left")
    idxr = np.searchsorted(signals1d["pos"], apos + expand, side="right")

    bidxl = np.searchsorted(signals1d["pos"], bpos - expand, side="left")
    bidxr = np.searchsorted(signals1d["pos"], bpos + expand, side="right")

    bidxl = max(bidxl, idxl)
    bidxr = max(bidxr, bidxl)

    result['left-only'].append(signals1d[idxl:idxr])
    result['right-only'].append(signals1d[bidxl:bidxr])

    #result['left-only'] = np.array(result['left-only'], dtype=[('chrom', 'S8'), ('pos', np.int32), ('event', 'i2')])
    #result['right-only'] = np.array(result['right-only'], dtype=[('chrom', 'S8'), ('pos', np.int32), ('event', 'i2')])
    result['left-only'].columns = ['chrom', 'pos', 'event'] 
    result['right-only'].columns = ['chrom', 'pos', 'event'] 

    return result

#@numba.jit(nopython=True)
def fill_arr(channels, events, posns, opos, offset):

    pmax = channels.shape[1]
   
    for i in range(len(events)):
        pos = posns[i] - opos + offset
        if pos < pmax and pos >= 0:
          channels[events[i], pos] += 1


def generate_channels_for_event(apos, bpos, signals1d, signals2d, expand, gap, depths):
    r = find_signals_for_event(apos, bpos, signals2d, signals1d, expand=expand)

    channels = np.zeros((max(Event) + 1, 4 * expand + gap), dtype=np.int32)
    # TODO: handle apos - expand < 0
    channels[0, 0:2*expand] = depths[apos - expand:apos + expand]
    channels[0, 2*expand+gap:] = depths[bpos - expand:bpos + expand]

    fill_arr(channels, np.asarray(r["shared"]["event"]), np.asarray(r["shared"]["a_pos"]), apos, expand)
    fill_arr(channels, np.asarray(r["shared"]["event"]), np.asarray(r["shared"]["b_pos"]), bpos, 3 * expand + gap)

    if len(r["left-only"]) > 0:
        fill_arr(channels, np.asarray(r["left-only"]["event"]),  np.asarray(r["left-only"]["pos"]), apos, expand)
    if len(r["right-only"]) > 0:
        fill_arr(channels, np.asarray(r["right-only"]["event"]), np.asarray(r["right-only"]["pos"]), bpos, 3 * expand + gap)

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

def main(args=sys.argv[1:]):
    import argparse
    import pandas as pd
    p = argparse.ArgumentParser()
    p.add_argument("directory", help="sample-specific directory created by svchannels extract")
    p.add_argument("bedpe")

    a = p.parse_args(args)

    depths_by_chrom = read_bins(a.directory)
    e2d = pd.read_table(f"{a.directory}/sv-channels.events2d.txt.gz", compression="gzip",
                        names=["a_chrom", "a_pos", "b_chrom", "b_pos", "event"],
                        dtype={"a_chrom": str, "b_chrom": str})

    e1d = pd.read_table(f"{a.directory}/sv-channels.soft_and_insertions.txt.gz", compression="gzip",
                        names=["chrom", "pos", "event"],
                        dtype={"chrom": str})
    print(f"[svchannels] read {len(e2d)} 2D events and {len(e1d)} 1d events", file=sys.stderr)
    t0 = time.time()
    n = 0
    expand = 250
    gap = 10

    for line in xopen(a.bedpe):
        toks = line.strip().split()
        if toks[0] != toks[3]: continue
        if int(toks[1]) < expand: continue
        if int(toks[4]) < expand: continue
        if int(toks[4]) < int(toks[1]):
            toks[4], toks[1] = toks[1], toks[4]
        n += 1
        # here, sv_chan is shape (n_channels, 2 * expand + gap) can accumulate these and send to learner.
        sv_chan = generate_channels_for_event(int(toks[1]), int(toks[4]), e1d, e2d, expand, gap, depths_by_chrom[toks[0]])
 
    t = time.time() - t0
    print(f"generated {n} channels in {t:.1f} seconds ({n/t:.0f} SVs/second)", file=sys.stderr)
    

    


if __name__ == "__main__":
    import pandas as pd

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

    #print("all")
    #print(signals2d.tolist())

    expand = 100
    gap = 10

    #r = find_signals_for_event(100, 400, signals2d, signals1d, expand=expand)
    #for k, v in r.items():
    #    print(k)
    #    print(v.tolist())

    #print(generate_channels_for_event(apos, bpos, signals1d, signals2d, expand, gap, depths))

    main()
