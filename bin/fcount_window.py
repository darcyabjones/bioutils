#! /usr/bin/env python3

from __future__ import print_function

import os
from os.path import split as psplit
from os.path import splitext as splitext
import re
import argparse
import sys
from collections import defaultdict

from Bio import SeqIO
import gffutils

program = "extract_features"
version = "0.1.0"
author = "Darcy Jones"
date = "1 October 2015"
email = "darcy.ab.jones@gmail.com"
short_blurb = (
    'Extracts GFF3 features from a reference fasta file.'
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


bedgraph_header = "track type=bedGraph\n"
bedgraph_row = "{seqid}\t{start}\t{end}\t{value}\n"

wig_header = "track type=wiggle\n"
wig_chrom_header = "fixedStep chrom={seqid} start=1 step={step} span={size}\n"
wig_row = "{value}\n"

def main(gff,
         outfile,
         format="bedgraph",
         window_size=10000,
         window_step=100,
         fasta=None,
         ftype='gene'):

    if fasta is not None:
        seqs = SeqIO.parse(fasta, format="fasta")
        lengths = dict()
        for seq in seqs:
            lengths[seq.id] = len(seq.seq)
    else:
        lengths = None

    # Parse GFF annotations.


    gff_db = gffutils.create_db(
        gff.name,
        dbfn=':memory:',
        id_spec=None,
        force=True,
        merge_strategy="create_unique"
        )

    if lengths is None:
        scafs = gff_db.execute("""
            SELECT seqid, MAX(end)
            FROM features GROUP BY seqid;
            """)
        lengths = {s["seqid"]: s["MAX(end)"] for s in scafs.fetchall()}

    out_rows = defaultdict(list)

    outfile.write(bedgraph_header)

    for scaf, length in lengths.items():
        i = 0
        while i < length:
            start = i
            end = i + window_size
            if end >= length:
                end = length

            features = gff_db.region(
                seqid=scaf,
                start=start + 1,
                end=end,
                featuretype=ftype,
                completely_within=True
                )

            count = sum(1 for _ in features)
            row = {
                "seqid": scaf,
                "start": start,
                "end": end,
                "value": count
                }

            outfile.write(bedgraph_row.format(**row))

            if end >= length:
                break

            i += window_step
    return

"########################### Argument Handling ##############################"

if __name__ == '__main__':
    arg_parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=license,
        )
    arg_parser.add_argument(
        "-f", "--fasta",
        default=None,
        type=argparse.FileType('r'),
        help=(
            "Path to the input file."
            "Default is stdin.\n"
            "Note that only one of infile and gff can take input from stdin."
            ),
        )
    arg_parser.add_argument(
        "-g", "--gff",
        default=sys.stdin,
        type=argparse.FileType('r'),
        help=(
            "Path to the GFF file."
            "Default is stdin.\n"
            ),
        )
    arg_parser.add_argument(
        "-o", "--outfile",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help=(
            "Path to write output to."
            "Default is stdout."
            ),
        )
    arg_parser.add_argument(
        "-r", "--format",
        default="bedgraph",
        choices=["bedgraph", "wig", "bigwig"],
        help="",
        )
    arg_parser.add_argument(
        "-t", "--type",
        dest='ftype',
        default='gene',
        help=(
            "The type of features to count."
            "Default is 'gene'."
            ),
        )
    arg_parser.add_argument(
        "-s", "--size",
        default=10000,
        type=int,
        help="Window size",
        dest="window_size"
        )
    arg_parser.add_argument(
        "-p", "--step",
        default=100,
        type=int,
        help="Window step",
        dest="window_step"
        )
    arg_parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s {}'.format(version),
        )

    args = arg_parser.parse_args()

    main(**args.__dict__)
