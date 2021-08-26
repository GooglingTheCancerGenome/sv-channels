import numpy as np
import numba

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
    aidxl = np.searchsorted(signals2d.a_pos, apos - expand, side="left")
    aidxr = np.searchsorted(signals2d.a_pos, apos + expand, side="right")

    # need b-side for events like splits from b that go to another chrom.
    bidxl = np.searchsorted(signals2d.b_pos, bpos - expand, side="left")
    bidxr = np.searchsorted(signals2d.b_pos, bpos + expand, side="right")

    # make sure we don't grab same values more than once.
    bidxl = max(bidxl, aidxr)
    bidxr = max(bidxr, bidxl)

    subset = signals2d[aidxl: aidxr]

    # same chrom and with expected bounds on both ends.
    selection = (subset.b_pos >= (bpos - expand)) & (subset.b_pos <= (bpos + expand)) & (subset.a_chrom == subset.b_chrom)
    result['shared'] = subset[selection]
    result['left-only'] = subset[~selection][["a_chrom", "a_pos", "event"]].tolist()

    subsetb = signals2d[bidxl: bidxr]
    selectionb = (subsetb.a_pos >= apos - expand) & (subsetb.a_pos <= apos + expand) & (subsetb.a_chrom == subsetb.b_chrom)
    result['right-only'] =  subsetb[~selectionb][["a_chrom", "a_pos", "event"]].tolist()

    
    # 1d search ...
    idxl = np.searchsorted(signals1d.pos, apos - expand, side="left")
    idxr = np.searchsorted(signals1d.pos, apos + expand, side="right")

    bidxl = np.searchsorted(signals1d.pos, bpos - expand, side="left")
    bidxr = np.searchsorted(signals1d.pos, bpos + expand, side="right")

    bidxl = max(bidxl, idxl)
    bidxr = max(bidxr, bidxl)

    result['left-only'].extend(signals1d[idxl:idxr])
    result['right-only'].extend(signals1d[bidxl:bidxr])

    return result

@numba.jit(nopython=True)
def fill_2d(channels, events, posns, pos, offset):
   
    for i in range(len(events)):
        pos = posns[i] - pos + offset
        print(posns[i], pos, offset)
        channels[events[i], pos] += 1


def generate_channels_for_event(r, apos, bpos, expand, gap):

    channels = np.zeros((max(Event), 4 * expand + gap), dtype=np.int32)

    fill_2d(channels, np.asarray(r["shared"].event), np.asarray(r["shared"].a_pos), apos, expand)
    fill_2d(channels, np.asarray(r["shared"].event), np.asarray(r["shared"].b_pos), bpos, 3 * expand + gap)


    return channels

if __name__ == "__main__":

    signals2d = np.array([
        ("chr1", 100, "chr3", 400, Event(1)),
        ("chr1", 200, "chr1", 300, Event(2)),
        ("chr1", 200, "chr1", 400, Event(3)),
        ("chr1", 200, "chr2", 400, Event(4)),
        ("chr1", 400, "chr1", 200, Event(5)),
        ("chr1", 500, "chr1", 900, Event(6))
        ], dtype=[('a_chrom', 'S8'), ('a_pos', np.int32),
                  ('b_chrom', 'S8'), ('b_pos', np.int32),
                  ('event', 'i2')
                  ]).view(np.recarray)

    signals1d = np.array([
        ("chr1", 100, Event(3)),
        ("chr1", 400, Event(3)),
        ], dtype=[('chrom', 'S8'), ('pos', np.int32), ('event', 'i2')]).view(np.recarray)


    apos, bpos = (100, 400)
    print(f"querying {apos}, {bpos}")

    print("all")
    print(signals2d.tolist())

    expand = 100
    gap = 10

    r = find_signals_for_event(100, 400, signals2d, signals1d, expand=expand)
    for k, v in r.items():
        print(k)
        if hasattr(v, "tolist"):
            print(v.tolist())
        else:
            print(v)

    print(generate_channels_for_event(r, apos, bpos, expand, gap))
