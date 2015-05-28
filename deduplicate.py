#!/usr/bin/env python3

#Version 1. Darcy Jones, January 2014.
#Contact Darcy Jones, darcy.ab.jones@gmail.com

###Import modules

import sys;
import argparse;
import os;
import re;
from difflib import SequenceMatcher, get_close_matches;
from Bio import SeqIO
from Bio.SeqUtils.CheckSum import seguid; #required for sequence comparison


###Variable Definitions

###Function definitions.

def checkProfile(id_, sgid):
	# Returns T/F. True means it is unique
	decision=True; 
	if sgid in profile_dict:
		possible_sequences = get_close_matches(id_, profile_dict[sgid], cutoff=0.1);
		if len(possible_sequences)>0:
			id_comparison = SequenceMatcher(None, id_, possible_sequences[0]);
			print(id_, possible_sequences, id_comparison.ratio(),'\n')
			if id_comparison.ratio() > LOS:
				decision=False; # Too similar
			else:
				decision=True;
		else:
			decision=True;
	else:
		decision=True;
	return decision;

def checkInput(id_, sgid):
	# Returns T/F. True means it is unique
	decision=True; 
	if sgid in seq_dict:
		possible_sequences = get_close_matches(id_, seq_dict[sgid], cutoff=0.1);
		if len(possible_sequences)>0:
			id_comparison = SequenceMatcher(None, id_, possible_sequences[0]);
			print(id_, possible_sequences, id_comparison.ratio(),'\n')
			if id_comparison.ratio() > LOS:
				decision=False; # Too similar
			else:
				decision=True;
		else:
			decision=True;
	else:
		decision=True;
	return decision;



###Code

def main(input_path, profile_path = None, output_path = '-'):
	print(input_path)
	print(profile_path)
	print(output_path)

	if profile_path =="-":
		profile_path = sys.stdin;

	global profile_dict;
	profile_dict = {};

	if profile_path != None:
		profile = SeqIO.parse(profile_path, 'fasta');
		for sequence in profile:
			full_id = "{}{}".format(sequence.id, sequence.description)
			sgid = seguid(sequence.seq);
			id_ = re.sub(':|;|\||_|\s|\*|\.|-|,|\[|\]|\{|\}|\(|\)|\'|\"','',full_id)
			if sgid in profile_dict:
				profile_dict[seguid(sequence.seq)].add(id_);
			else:
				profile_dict[seguid(sequence.seq)]={id_};

	global seq_dict, LOS;
	seq_dict=dict()
	LOS=0.6

	write_list=[]
	for file_ in input_path:
		if os.path.isfile(file_):
			these_seqs = SeqIO.parse(file_, 'fasta');
			for seq in these_seqs:
				full_id = "{}{}".format(seq.id, seq.description)
				id_=re.sub(':|;|\||_|\s|\*|\.|-|,|\[|\]|\{|\}|\(|\)|\'|\"','',full_id);
				sgid=seguid(seq.seq)
				if checkProfile(id_, sgid) == True:
					if checkInput(id_, sgid) == True:
						write_list.append(seq);
						if sgid in seq_dict:
							seq_dict[sgid].add(id_);
						else:
							seq_dict[sgid] = {id_};

	if output_path == '-':
		output_path = sys.stdout;

	SeqIO.write(write_list, output_path, 'fasta')

if __name__== '__main__':

	###Argument Handling
	arg_parser=argparse.ArgumentParser(description='description');
	arg_parser.add_argument("-f", '--input_path', nargs='*', help="directory to input file(s) or directory, enter '-' for stdin (default). You may indicate multiple files by separating with a space eg. -f one.faa two.faa etc.");
	arg_parser.add_argument("-p", '--profile_path', default=None, help="directory to profile file");
	arg_parser.add_argument("-o", '--output_path', default="-", help="directory to output file, enter '-' for stdout (default)");
	args = arg_parser.parse_args();

	input_path=args.input_path;
	profile_path=args.profile_path;
	output_path=args.output_path;

	main(input_path, profile_path, output_path);
