import sys
from . import extract_signals
from . import generate_channels

commands = {
        'extract-signals': extract_signals.main,
        'generate-channels': generate_channels.main,
    }

def main(args=sys.argv[1:] if len(sys.argv) > 0 else []):
    if len(args) == 0 or not args[0] in commands:
        print(f"Command '{args[0]}' not found.\n Available commands are:", file=sys.stderr)
        for cmd in commands:
            print(f'\t{cmd}', file=sys.stderr)
        sys.exit(1)

    commands[args[0]](args[1:])

if __name__ == "__main__":
    main()
