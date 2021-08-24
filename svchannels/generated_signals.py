import numpy as np
import numba

def make_array(expand, gap, dtype=np.int32):
    return np.zeros(4 * expand + gap, dtype=dtype)

def generate_signals(apos, bpos, signals2d, signals1d, expand=250):
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
    print(subset, bpos)

    # same chrom and with expected bounds on both ends.
    selection = (subset.b_pos >= (bpos - expand)) & (subset.b_pos <= (bpos + expand)) & (subset.a_chrom == subset.b_chrom)
    result['shared'] = subset[selection]
    result['left-only'] = subset[~selection]

    subsetb = signals2d[bidxl: bidxr]
    selectionb = (subsetb.a_pos >= apos - expand) & (subsetb.a_pos <= apos + expand) & (subsetb.a_chrom == subsetb.b_chrom)
    result['right-only'] =  subsetb[~selectionb]

    return result


if __name__ == "__main__":

    signals2d = np.array([
        ("chr1", 200, "chr1", 400),
        ("chr1", 200, "chr2", 400),
        ("chr1", 400, "chr1", 200),
        ("chr1", 500, "chr1", 900)
        ], dtype=[('a_chrom', 'S8'), ('a_pos', np.int32),
                  ('b_chrom', 'S8'), ('b_pos', np.int32)
                  ]).view(np.recarray)
    print(signals2d.a_chrom)

    print(generate_signals(100, 400, signals2d, []))





    

