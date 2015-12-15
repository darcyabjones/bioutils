#! /usr/bin/env python3

from __future__ import print_function

program = "bioutils convert"
version = "0.1.0"
author = "Darcy Jones"
date = "1 October 2015"
email = "darcy.ab.jones@gmail.com"
short_blurb = (
    "Converts sequences between file formats."
    )
license = (
    '{program}-{version}\n'
    '{short_blurb}\n\n'
    'Copyright (C) {date},  {author}'
    '\n\n'
    'This program is free software: you can redistribute it and/or modify '
    'it under the terms of the GNU General Public License as published by '
    'the Free Software Foundation, either version 3 of the License, or '
    '(at your option) any later version.'
    '\n\n'
    'This program is distributed in the hope that it will be useful, '
    'but WITHOUT ANY WARRANTY; without even the implied warranty of '
    'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the '
    'GNU General Public License for more details.'
    '\n\n'
    'You should have received a copy of the GNU General Public License '
    'along with this program. If not, see <http://www.gnu.org/licenses/>.'
    )

license = license.format(**locals())

############################# Import all modules #############################

import sys
import argparse

from Bio import SeqIO
from Bio.Alphabet import generic_dna
from Bio.Alphabet import generic_rna
from Bio.Alphabet import generic_protein

################################## Functions #################################

def inhandler(fp, mode='r'):
    if fp == sys.stdin or fp == '-':
        return sys.stdin
    elif isinstance(fp, str):
        return open(fp, mode)
    else:
        return fp

def outhandler(fp, mode='w'):
    if fp == sys.stdout or fp == '-':
        return sys.stdout
    elif isinstance(fp, str):
        return open(fp, mode)
    else:
        return fp

alphabets = {
    'dna': generic_dna,
    'rna': generic_rna,
    'aa': generic_protein,
    }

#################################### Main ####################################

def main(infile, outfile, fromfmt='fasta', tofmt='fasta', alphabet='dna'):
    with inhandler(infile) as inhandle, outhandler(outfile) as outhandle:
        SeqIO.convert(
            inhandle,
            fromfmt,
            outhandle,
            tofmt,
            alphabet=alphabets[alphabet]
            )
    return

if __name__== '__main__':
    arg_parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=license,
        epilog=(
            'Example usage:\n'
            '$ %(prog)s -i my_fasta.fna -o my_fasta.gbk --from fasta --to genbank --alphabet dna\n'
            )
        )
    arg_parser.add_argument(
        "-i", "--infile",
        default=sys.stdin,
        type=argparse.FileType('r'),
        help=(
            "Path to the input file."
            "Default is stdin."
            )
        )
    arg_parser.add_argument(
        "-o", "--outfile",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help=(
            "Path to write output to."
            "Default is stdout."
            )
        )
    arg_parser.add_argument(
        "-f", "--from",
        dest='fromfmt',
        default='fasta',
        choices=[
            'fasta', 'embl', 'genbank', 'pir',
            'seqxml', 'swiss', 'tab', 'uniprot-xml'
            ],
        help="The format of the input sequence. Default is 'fasta'."
        )
    arg_parser.add_argument(
        "-t", "--to",
        dest='tofmt',
        default='fasta',
        choices=[
            'fasta', 'embl', 'genbank', 'pir',
            'seqxml', 'swiss', 'tab', 'uniprot-xml'
            ],
        help="The format of the output sequence. Default is 'fasta'."
        )
    arg_parser.add_argument(
        "-a", "--alphabet",
        choices=['dna', 'rna', 'aa'],
        default='dna',
        help="",
        )
    arg_parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s {}'.format(version),
        )


    args = arg_parser.parse_args()

    main(**args.__dict__)
