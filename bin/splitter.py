#! /usr/bin/env python3

from __future__ import print_function

program = "bioutils splitter"
version = "0.1.0"
author = "Darcy Jones"
date = "1 October 2015"
email = "darcy.ab.jones@gmail.com"
short_blurb = ('')
license = (
    '{program}-{version}\n'
    '{short_blurb}\n'
    '\n'
    'Copyright (C) {date},  {author}\n'
    '\n'
    'This program is free software: you can redistribute it and/or modify '
    'it under the terms of the GNU General Public License as published by '
    'the Free Software Foundation, either version 3 of the License, or '
    '(at your option) any later version.\n'
    '\n'
    'This program is distributed in the hope that it will be useful, '
    'but WITHOUT ANY WARRANTY; without even the implied warranty of '
    'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the '
    'GNU General Public License for more details.\n'
    '\n'
    'You should have received a copy of the GNU General Public License '
    'along with this program. If not, see <http://www.gnu.org/licenses/>.'
    )

license = license.format(**locals())

############################ Import all modules ##############################

import os
from os.path import split as psplit
from os.path import splitext as splitext
import re
import argparse
import sys
from collections import defaultdict

################################# Functions ##################################

fmt_to_ext = {
    'fasta': '.fa',
    'embl': '.embl',
    'genbank': '.gbk',
    'pir': '.pir',
    'seqxml': '.xml',
    'swiss': '.sw',
    'tab': '.tsv',
    'uniprot-xml': '.xml',
    }

def inhandler(fp, mode='r'):
    if fp == sys.stdin or fp == '-':
        return sys.stdin
    elif isinstance(fp, str):
        return open(fp, mode)
    else:
        return fp

def main(infile, prefix, num_records, fmt):
    num_records = int(num_records)

    chunk_num = 0
    filepaths = list()

    if prefix is None and infile == sys.stdin:
        prefix = 'split'
    elif prefix is None:
        prefix = splitext(infile.name)[0]

    if infile == sys.stdin:
        ext = fmt_to_ext[fmt]
    else:
        ext = splitext(infile.name)[1]

    with inhandler(infile) as inhandle:
        sequences = SeqIO.parse(inhandle, format=fmt)

        record_num = 0
        new_sequences = list()

        for sequence in sequences:
            if record_num >= num_records:
                fp = prefix + '-' + str(chunk_num) + ext
                SeqIO.write(new_sequences, fp, format=fmt)
                new_sequences = list()
                chunk_num += 1
                record_num = 0
            else:
                new_sequences.append(sequences)

        if len(new_sequences) > 0:
            fp = prefix + '-' + str(chunk_num) + ext
            SeqIO.write(new_sequences, fp, format=fmt)

    return

############################ Argument Handling ###############################

if __name__== '__main__':
    arg_parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=license,
        )
    arg_parser.add_argument(
        "-i", "--infile",
        default=sys.stdin,
        type=argparse.FileType('r'),
        help=(
            "Path to the input file. "
            "Default is stdin."
            )
        )
    arg_parser.add_argument(
        "-p","--prefix",
        default=None,
        help=""
        )
    arg_parser.add_argument(
        "-n",
        dest='num_records',
        default=100,
        type=int,
        help=""
        )
    arg_parser.add_argument(
        "-f", "--format",
        dest='fmt',
        default='fasta',
        choices=[
            'fasta', 'embl', 'genbank', 'pir',
            'seqxml', 'swiss', 'tab', 'uniprot-xml'
            ],
        help="The format of the sequences. Default is 'fasta'."
        )
    arg_parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s {}'.format(version),
        )

    args = arg_parser.parse_args()

    main(**args.__dict__)
