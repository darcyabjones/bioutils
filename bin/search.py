#! /usr/bin/env python3

from __future__ import print_function

program = "seq_grep"
version = "0.1.0"
author = "Darcy Jones"
date = "1 October 2015"
email = "darcy.ab.jones@gmail.com"
short_blurb = (
    "Retrieve sequences from a sequence file based on regular expression matches."
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

from Bio import SeqIO

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

def re_checker(field, regexes):
    if field == 'id':
        def any_re(record):
            for regex in regexes:
                if regex.search(record.id) is not None:
                    return True
            return False
    elif field == 'description':
        def any_re(record):
            for regex in regexes:
                if regex.search(record.description) is not None:
                    return True
            return False
    elif field == 'seq':
        def any_re(record):
            for regex in regexes:
                if regex.search(str(record.seq)) is not None:
                    return True
            return False
    return any_re

#################################### Main ####################################

def main(
        infile,
        outfile,
        regex,
        field='seq',
        ignorecase=False,
        locale=False,
        dotall=False,
        multiline=False,
        debug=False,
        unicode_=False,
        fmt='fasta'
        ):
    regex_flags=0
    if ignorecase:
        regex_flags |= re.IGNORECASE
    if locale:
        regex_flags |= re.LOCALE
    if dotall:
        regex_flags |= re.DOTALL
    if multiline:
        regex_flags |= re.MULTILINE
    if debug:
        regex_flags |= re.DEBUG
    if unicode_:
        regex_flags |= re.UNICODE

    regexes = list()
    for r in regex:
        regexes.append(re.compile(r, flags=regex_flags))

    checker = re_checker(field, regexes)

    with inhandler(infile) as inhandle, outhandler(outfile) as outhandle:
        out_seqs = list()
        seqs = SeqIO.parse(inhandle, format=fmt)
        for seq in seqs:
            if checker(seq):
                out_seqs.append(seq)

        SeqIO.write(out_seqs, outhandle, format=fmt)
    return

if __name__== '__main__':
    arg_parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=license,
        epilog=(
            'Example usage:\n'
            '$ %(prog)s -i my_fasta.faa -o my_fasta.faa "\*\S"\n'
            '\n'
            'Would find all sequences with internal stops.\n'
            )
        )
    arg_parser.add_argument(
        "regex",
        nargs='+',
        help="The regular expressions to search for."
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
        "-e", "--field",
        default='seq',
        choices=['id', 'description', 'name'],
        help="The field to look for the pattern in. Default is 'seq' (ie the sequence)."
        )
    arg_parser.add_argument(
        "--ignorecase",
        default=False,
        action='store_true',
        help=(
            'Perform case-insensitive matching; expressions like [A-Z] '
            'will match lowercase letters, too. '
            'This is not affected by the current locale.'
            ),
        )
    arg_parser.add_argument(
        "--locale",
        default=False,
        action='store_true',
        help=r'Make \w, \W, \b, \B, \s and \S dependent on the current locale.',
        )
    arg_parser.add_argument(
        "--dotall",
        default=False,
        action='store_true',
        help=(
            "Make the '.' special character match any character at all, "
            "including a newline; without this flag, '.' will match "
            'anything except a newline.'
            ),
        )
    arg_parser.add_argument(
        "--multiline",
        default=False,
        action='store_true',
        help=(
            "When specified, the pattern character '^' matches at the "
            "beginning of the string and at the beginning of each line "
            "(immediately following each newline); and the pattern character "
            "'$' matches at the end of the string and at the end of each line "
            "(immediately preceding each newline).\n"
            "By default, '^' matches only at the beginning of the string, "
            "and '$' only at the end of the string and immediately before the "
            "newline (if any) at the end of the string."
            ),
        )
    arg_parser.add_argument(
        "--unicode",
        dest='unicode_',
        default=False,
        action='store_true',
        help=(
            r'Make \w, \W, \b, \B, \d, \D, \s and \S dependent on the '
            r'Unicode character properties database.'
            ),
        )
    arg_parser.add_argument(
        "--debug",
        default=False,
        action='store_true',
        help='Display debug information about compiled expression.',
        )
    arg_parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s {}'.format(version),
        )


    args = arg_parser.parse_args()

    main(**args.__dict__)
