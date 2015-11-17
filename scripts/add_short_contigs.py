#!/usr/bin/env python3

#
# Copyright 2015 Dovetail Genomics LLC
#
#


"""
This script adds short contigs back into the scaffold fasta. 

Limitations:
- It does not check that you are using the correct length cutoff. If
the length cutoff for HiRise is changed, then the cutoff for this 
script must be the same. Otherwise, you may end up with missing or
duplicated sequence.

- It does not generate scaffold names in the same way as HiRise scaffolds.
Instead, it keeps the same headers as the broken.fa file. 

- It does not change the line length of the broken.fa file. If the HiRise
scaffolds were written with a different line length, then the short
contigs will have a different line length than the scaffolds.
"""

import argparse


def main(broken=None, scaffolds=None, length=1000, output=None):
    """Takes as input a broken.fa file and a HiRise scaffold file and 
    outputs the scaffolds with the short contigs appended to the end."""

    # Open output and write HiRise scaffolds
    with open(output, 'w') as out_handle:
        with open(scaffolds) as scaff_handle:
            for line in scaff_handle:
                out_handle.write(line)

        # Writes small contigs to output
        for header, seq, string_seq in get_seqs(broken):
            if len(seq) < length:
                out_handle.write(header)
                out_handle.write(string_seq)
                out_handle.write("\n")

            
def get_seqs(fasta):
    """Takes as input a fasta file and yields the header line,
    the sequence stripped on newlines, and the sequence with new
    lines as they were in the fasta file."""

    with open(fasta) as fasta_handle:
        header = None
        for line in fasta_handle:
            if line.startswith(">"):
                if header:
                    yield header, ''.join(seq), '\n'.join(seq)
                seq = []
                header = line
            else:
                seq.append(line.rstrip())

        # Yield final sequence
        yield header, ''.join(seq), '\n'.join(seq)
            
def get_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--broken", help="The broken.fa file.")
    parser.add_argument("-s", "--scaffolds", help="The HiRise scaffold fasta file.", default="/dev/stdin")
    parser.add_argument("-l", "--length", help="The maximum length of contigs to include", default=1000)
    parser.add_argument("-o", "--output", help="Output file", default="/dev/stdout")
    parser.add_argument("-d", "--debug", help="Run pdb", default=False, action="store_true")
    return parser.parse_args()

if __name__ == "__main__":
    args = get_args()

    if args.debug:
        import pdb
        pdb.set_trace()

    main(broken=args.broken, scaffolds=args.scaffolds,
         length=args.length, output=args.output)
