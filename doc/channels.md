sv-channels extracts signals from a BAM or CRAM file and outputs them in the form of *channels*, 1-dimensional arrays
where the index in the array reflects relative genomic position and the value indicates the number of pieces of evidence
for that channel at a given position. These channels are used to train a deep neural network to filter structural variants.

The channels used in `sv-channels` are somewhat self-documented in the Enum used by the code:

```python
class Event(enum.IntEnum):
    # NOTE that python enums start at 1. The 0th track is depth.

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
```

These are the signals extracted. When considering a specific SV, the final 8 channels (starting at `SPLIT_PLUS_PLUS` can be either 
paired or "orphaned". For example, a paired split read would be when the ends (splits) of the split read are close to the break-points
defined by the SV. An orphaned split would occur when, for example, the left-end fell within the window, but that split went to another
location in the genome that did not support the event.

Reads are extracted according to the --min-mapping-quality parameter given on the command-line. Reads with a lower mapping-quality are not
used expect for `SPLIT_LOW_QUALITY`.
