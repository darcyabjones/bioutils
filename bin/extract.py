#! /usr/bin/env python3

from __future__ import print_function

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

############################ Import all modules ##############################

import os
from os.path import split as psplit
from os.path import splitext as splitext
import re
import argparse
import sys
from collections import defaultdict

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature
from Bio.SeqFeature import CompoundLocation
from BCBio import GFF

################################## Classes ###################################

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


def subfeatures(feature):
    """ Return a location object from GFF output.

    The BCBio GFF parser adds exons, mRNA, and CDS's to features as
    sub_features which has since been depreciated in Biopython in favour
    of the CompoundLocation object. This function returns a Location
    object that we can use to extract sequences.

    Keyword arguments:
    feature -- A SeqFeature record with the _sub_features attribute.

    Returns:
    An ExactLocation or CompundLocation object.
    """
    new_features = list()

    if not hasattr(feature, "_sub_features"):
        return [feature]
    if feature._sub_features is None:
        return [feature]

    sub_features = feature._sub_features
    feature._sub_features = []
    new_features.append(feature)

    sub_features_cds = [f for f in sub_features if f.type.lower() == 'cds']
    sub_features_exons = [f for f in sub_features if f.type.lower() == 'exon']
    sub_features_mrna = [f for f in sub_features if f.type.lower() == 'mrna']
    sub_features_others = [f for f in sub_features
                           if f.type.lower() not in {'cds', 'exon', 'mrna'}]
    new_features.extend(sub_features_others)

    strand = feature.strand
    if len(sub_features_exons) > 0:
        sub_features_exons.sort(key=lambda l: l.location.start)
        locations = [f.location for f in sub_features_exons]
        if strand == -1:
            """ When calling CompoundLocation.extract() the sequences are
            extracted in the order that they are encountered. For features on
            the - strand, we need to reverse this order. """
            locations.reverse()
        qualifiers = sub_features_exons[0].qualifiers
        sub_feature = SeqFeature(
            id=sub_features_exons[0].id,
            type="CDS",
            strand=strand,
            qualifiers=qualifiers,
            location=CompoundLocation(locations)
            )
        new_features.append(sub_feature)

    if len(sub_features_cds) > 0:
        if len(sub_features_cds) > 1:
            sub_features_cds.sort(key=lambda l: l.location.start)
            locations = [f.location for f in sub_features_cds]
            if strand == -1:
                locations.reverse()
            locations = CompoundLocation(locations)
        else:
            # One CDS returns an ExactLocation
            locations = sub_features_cds[0].location

        qualifiers = sub_features_cds[0].qualifiers
        sub_feature = SeqFeature(
            id=sub_features_cds[0].id,
            type="CDS",
            strand=strand,
            qualifiers=qualifiers,
            location=locations,
            )
        new_features.append(sub_feature)

    if len(sub_features_mrna) > 0:
        for mrna in sub_features_mrna:
            new_features.extend(subfeatures(mrna))

    return new_features


def main(infile, gff, outfile, ftype='CDS', use_phase=False, translate=False):
    with inhandler(infile) as handle:
        ref_seq = SeqIO.to_dict(SeqIO.parse(handle, format="fasta"))
    # Parse GFF annotations.

    with inhandler(gff) as handle:
        genome_with_features = GFF.parse(
            handle,
            base_dict=ref_seq
            )
        """ bcbio-gff codes exons, mRNA etc as subfeatures which is now
        depreciated in biopython, this code fixes that issue. """
        new_genome_with_features = list()
        for scaffold in genome_with_features:
            new_features = list()
            for feature in scaffold.features:
                gene_features = subfeatures(feature)
                new_features.extend(gene_features)
            scaffold.features = new_features
            new_genome_with_features.append(scaffold)
        """ Genome with features doesn't have scaffolds without any gff
        features. Here I update the existing records in genome with the
        new ones containing features. """
        ref_seq.update(SeqIO.to_dict(new_genome_with_features))

    sequences = list()
    for scaffold, sequence in ref_seq.items():
        for feature in sequence.features:
            if feature.type != ftype:
                continue
            start = feature.location.start
            end = feature.location.end
            phase = int(feature.qualifiers['phase'][0])
            strand = feature.location.strand

            if use_phase:
                fseq = feature.extract(sequence)[phase:]
            else:
                fseq = feature.extract(sequence)

            fseq.id = feature.id
            fseq.name = feature.id

            strand = '-' if strand == -1 else '+'
            fseq.description = "{}:{}-{}[{}]{}".format(
                scaffold,
                start,
                end,
                strand,
                phase
                )
            if translate:
                tseq = fseq.seq.translate()
                fseq.seq = tseq
            sequences.append(fseq)

    with outhandler(outfile) as handle:
        SeqIO.write(sequences, handle, 'fasta')
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
            "Note that only one of 'infile' and 'gff' can take input from stdin."
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
        "-t", "--type",
        dest='ftype',
        default='CDS',
        help=(
            "The type of features to extract sequences from."
            "Default is 'CDS'."
            ),
        )
    arg_parser.add_argument(
        "-p", "--translate",
        default=False,
        action='store_true',
        help="Translate sequences."
        )
    arg_parser.add_argument(
        "-a", "--phase",
        dest='use_phase',
        default=False,
        action='store_true',
        help="Use phase information."
        )
    arg_parser.add_argument(
        '--version',
        action='version',
        version='%(prog)s {}'.format(version),
        )

    args = arg_parser.parse_args()

    main(**args.__dict__)
