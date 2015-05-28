#!/usr/local/bin/env python
#python 2.7.5 requires biopython

########### dispomFormat ############

#Writes a new fasta file with non-wrapped sequences and \r newline characters appropriate for dispom

#Version 1. Darcy Jones, January 2014.
#Contact Darcy Jones, darcy.ab.jones@gmail.com

###Import modules
from Bio import SeqIO;
import argparse;
import sys;

###Code

def main(infile, outfile, trim=None):
    sequences=SeqIO.parse(infile, 'fasta');
    #   Read in standard fasta file.

    seq_dict={};
    len_dict={};
    deleted_list=[];

    for seq in sequences:
        if 'N' not in seq.seq[0:trim]:
        #   If N is not in the sequence. If trim == None this will exclude any sequences with N's in them, if trim is specified this will exclude any sequences with Ns in the sequence beyond the trim threshold.
            seq_dict[seq.id]=seq;
            seq_length=len(seq.seq);
            if seq_length in seq_dict:
                len_dict[seq_length].append(seq.id);
            else:
                len_dict[seq_length]=[seq.id];
        else:
            deleted_list.append(seq.id);

    # Delete short sequences and remove N's if trim is specified.
    if trim==None:
        del len_dict[max(len_dict.keys())];
        for len_ in len_dict:
            for id_ in len_dict[len_]:
                del seq_dict[id_];
                deleted_list.append(id_);
    else:
        shortest=max(len_dict.keys()); # keeps track of the shortest allowed length to keep within parameters.
        for len_ in len_dict:
            if len_<trim: # If the length of the sequence is less than the minimum length that youve set, exclude all sequences of that length.
                for id_ in len_dict[len_]:
                    del seq_dict[id_];
                    deleted_list.append(id_);
            elif len_<shortest: 
                shortest=len_;

        # Find and trim sequences with N's.
        for seq in seq_dict:
            firstN=seq_dict[seq].seq.lower().find('n');
            if firstN < shortest and firstN != -1: # If an N occurs earlier than the previous shortest value, the position becomes the new shortest value.
                shortest=firstN;

        trim=shortest # Set trim to the shortest allowed length.

    if outfile==sys.stdout:
        outfile_handle=outfile;
    else:
        outfile_handle=open(outfile, 'w');
    for seq in seq_dict:
        record=">{}\r{}\r".format(seq_dict[seq].id, seq_dict[seq].seq[0:trim]);
        outfile_handle.write(record);
    if outfile!=sys.stdout:
        outfile_handle.close;
    print(deleted_list);

if __name__=='__main__':
    ###Argument Handling
    arg_parser = argparse.ArgumentParser(description='Writes a new fasta file with non-wrapped sequences and \\r newline characters. Appropriate for dispom input.');
    arg_parser.add_argument("infile", help="Directory to line wrapped input file. To use stdin (for piped input) enter '-'");
    arg_parser.add_argument("outfile", help="Directory of new fasta to be written in non-wrapped format. To use stdout (for piped output) enter '-'");
    arg_parser.add_argument("-t", "--trim", type=int, default=None, help='Enter an integer, trims sequences up to this integer untill there are no N\'s and all sequences are the same length. Excludes sequences that contain N\'s or are still short after this trim threshold is met. Default = trimming not allowed.')
    args = arg_parser.parse_args();

    infile=args.infile;
    outfile=args.outfile;
    if infile=='-':
        infile=sys.stdin;
    if outfile=='-':
        outfile=sys.stdout;
    trim=args.trim;

    main(infile, outfile, trim);