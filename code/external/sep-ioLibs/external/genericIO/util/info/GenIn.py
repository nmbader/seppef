#!/usr/bin/env python3
import argparse
import genericIO
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Print info about files')
    parser.add_argument('files', metavar='Files', type=str, nargs='+',
                        help='A list of files to give more info')
    parser.add_argument("--io", type=str,choices=[@GEN_IO_TYPES@], help='IO type. Defaults to defaultIO')
    args = parser.parse_args()

    io=genericIO.defaultIO;

    if args.io:
        io=genericIO.io(args.io)

    for f in args.files:
        print(genericIO.defaultIO.getRegFile(f))




