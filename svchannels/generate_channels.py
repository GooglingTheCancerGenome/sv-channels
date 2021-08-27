import numpy as np
import numba
import zarr

import sys
import os

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
    result['left-only'] = subset[~selection][["a_chrom", "a_pos", "event"]].tolist()

    subsetb = signals2d[bidxl: bidxr]
    selectionb = (subsetb["a_pos"] >= apos - expand) & (subsetb["a_pos"] <= apos + expand) & (subsetb["a_chrom"] == subsetb["b_chrom"])
    result['right-only'] =  subsetb[~selectionb][["a_chrom", "a_pos", "event"]].tolist()

    
    # 1d search ...
    idxl = np.searchsorted(signals1d["pos"], apos - expand, side="left")
    idxr = np.searchsorted(signals1d["pos"], apos + expand, side="right")

    bidxl = np.searchsorted(signals1d["pos"], bpos - expand, side="left")
    bidxr = np.searchsorted(signals1d["pos"], bpos + expand, side="right")

    bidxl = max(bidxl, idxl)
    bidxr = max(bidxr, bidxl)

    result['left-only'].extend(signals1d[idxl:idxr])
    result['right-only'].extend(signals1d[bidxl:bidxr])
    result['left-only'] = np.array(result['left-only'], dtype=[('chrom', 'S8'), ('pos', np.int32), ('event', 'i2')])
    result['right-only'] = np.array(result['right-only'], dtype=[('chrom', 'S8'), ('pos', np.int32), ('event', 'i2')])

    return result

@numba.jit(nopython=True)
def fill_arr(channels, events, posns, opos, offset):

    pmax = channels.shape[1]
   
    for i in range(len(events)):
        pos = posns[i] - opos + offset
        if pos < pmax and pos >= 0:
          channels[events[i], pos] += 1


def generate_channels_for_event(apos, bpos, signals1d, signals2d, expand, gap, depths):
    r = find_signals_for_event(apos, bpos, signals2d, signals1d, expand=expand)

    channels = np.zeros((max(Event), 4 * expand + gap), dtype=np.int32)
    # TODO: handle apos - expand < 0
    channels[0, 0:2*expand] = depths[apos - expand:apos + expand]
    channels[0, 2*expand+gap:] = depths[bpos - expand:bpos + expand]

    fill_arr(channels, np.asarray(r["shared"]["event"]), np.asarray(r["shared"]["a_pos"]), apos, expand)
    fill_arr(channels, np.asarray(r["shared"]["event"]), np.asarray(r["shared"]["b_pos"]), bpos, 3 * expand + gap)

    fill_arr(channels, r["left-only"]["event"],  r["left-only"]["pos"], apos, expand)
    fill_arr(channels, r["right-only"]["event"], r["right-only"]["pos"], bpos, 3 * expand + gap)

    return channels

if __name__ == "__main__":

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

    signals1d = np.array([
        ("chr1", 11100, Event(3)),
        ("chr1", 11400, Event(3)),
        ], dtype=[('chrom', 'S8'), ('pos', np.int32), ('event', 'i2')])

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

    print(generate_channels_for_event(apos, bpos, signals1d, signals2d, expand, gap, depths))
