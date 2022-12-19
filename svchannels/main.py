import sys
from . import extract_signals
from . import generate_channels
from . import label
from . import score
from . import train

commands = {
        'extract-signals': ("extract SV signals from a BAM/CRAM", extract_signals.main),
        'generate-channels': ("generate matricies (channels) from signals and a set of SVs", generate_channels.main),
        'label': ("Label SVs", label.main),
        'score': ("score a set of variants (channels) given a trained model", score.main),
        'train': ("train a model from a set of variants (channels)", train.main)
    }


def main(args=sys.argv[1:] if len(sys.argv) > 0 else []):
    if len(args) == 0 or not args[0] in commands:
        print("svchannels:", file=sys.stderr)
        if len(args) > 0 and args[0] not in ('-h', '--help'):
            print(f"Command '{args[0]}' not found.\n Available commands are:", file=sys.stderr)
        for cmd in commands:
            print(f'\t{cmd}: {commands[cmd][0]}', file=sys.stderr)
        sys.exit(1)

    commands[args[0]][1](args[1:])


if __name__ == "__main__":
    main()
