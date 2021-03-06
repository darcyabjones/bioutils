#! /usr/bin/env python3

from __future__ import print_function

program = "bioutils deduplicate"
version = "0.1.0"
author = "Darcy Jones"
date = "1 October 2015"
email = "darcy.ab.jones@gmail.com"
short_blurb = (
    "Remove duplicate sequences from a sequence file."
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
import re
from collections import defaultdict
from difflib import SequenceMatcher

from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid

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

PUNCTUATION = re.compile(':|;|\||_|\s|\*|\.|-|,|\[|\]|\{|\}|\(|\)|\'|\"')
def junk(x):
    if PUNCTUATION.match(x) is not None:
        return True
    else:
        return False

def get_close_matches(word, possibilities, cutoff=0.6, junkfn=None):
    matches = list()
    for poss in possibilities:
        match = SequenceMatcher(junkfn, word, poss)
        if match.ratio() > cutoff:
            matches.append(poss)
    return matches

def get_unique(seqids, cutoff=0.6, junkfn=None):
    if len(seqids) < 2:
        return seqids

    matches = get_close_matches(
        seqids[0],
        seqids[1:],
        cutoff=cutoff,
        junkfn=junkfn,
        )

    seqids = seqids[:1] + [s for s in seqids[1:] if s not in matches]
    if len(seqids) > 2:
        return seqids[0] + get_unique(seqids[1:], cutoff=cutoff, junkfn=junkfn)
    else:
        return seqids


#################################### Main ####################################

def main(infile, outfile, match_id=True, cutoff=0.6, fmt='fasta'):
    with inhandler(infile) as inhandle:
        checksums = defaultdict(list)
        new_seqs = dict()
        seqs = SeqIO.parse(inhandle, format=fmt)
        for seq in seqs:
            new_seqs[seq.id] = seq
            sgid = seguid(seq.seq)
            checksums[sgid].append(seq.id)
        seqs = new_seqs

    keep = list()
    for sgid, seqids in checksums.items():
        if match_id:
            keep.extend(get_unique(
                seqids,
                cutoff=cutoff,
                junkfn=junk,
                ))
        else:
            keep.append(seqids[0])

    with outhandler(outfile) as outhandle:
        new_seqs = list()
        for seqid in keep:
            new_seqs.append(seqs[seqid])
        SeqIO.write(new_seqs, outfile, format=fmt)
    return


if __name__== '__main__':
    arg_parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=license,
        epilog=(
            'Example usage:\n'
            '$ %(prog)s -i my_fasta.fna -o my_fasta.faa\n'
            )
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
        "-o", "--outfile",
        default=sys.stdout,
        type=argparse.FileType('w'),
        help=(
            "Path to write output to. "
            "Default is stdout."
            )
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
        "-m", "--no-id",
        dest='match_id',
        default=True,
        action='store_false',
        help="",
        )
    arg_parser.add_argument(
        "-c", "--cutoff",
        default=0.6,
        type=float,
        help="",
        )
    arg_parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s {}'.format(version),
        )


    args = arg_parser.parse_args()

    main(**args.__dict__)
