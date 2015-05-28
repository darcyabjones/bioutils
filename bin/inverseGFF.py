#!/usr/local/bin/env python
#python 2.7.5 requires biopython

########### inverseGFF ############

#A tool that takes GFF file with N-block boundaries and an assembled genome and reconstructs contigs.


#Version 1. Darcy Jones, December 2013.
#Contact Darcy Jones, darcy.ab.jones@gmail.com

###Import modules
import sys;
import argparse;
from Bio import SeqIO;

###Function definitions.

def gffList(gff_handle_):
#	Reads through all lines in .gff file, extracts seqid, start and stop positions, returns list of lists with format ['seqid'[(start, stop)...]].
	temp_dict={};
	for line in gff_handle_:
		record=line.rstrip('\n').split('\t');
		#	removes newline characters, splits into list by tabs.
		if len(record)>1:
		#	Excludes first line.
			if record[0] not in temp_dict:
			#	first dict entry for a scaffold.
				temp_dict[record[0]]=[(int(record[3]), int(record[4]))];
			else:
			# continues adding (start, stop) tuples for each scaffold (seqid).
				temp_dict[record[0]].append((int(record[3]), int(record[4])));
	
	for scaffold in temp_dict:
		temp_dict[scaffold].sort(key=lambda tup: tup[0]); # Sort in place by start position, small to big
	return temp_dict;

def sequenceBoundaries(gff_list_, output_dir_, seq_index_):
	scaf_blocks_={};
	for scaffold in seq_index_:
		if scaffold in gff_list_:
			inverse_nn_=gff_list_[scaffold];
			temp_list=[];
			try:
				scaffold_length=len(seq_index_[scaffold].seq);
			except(KeyError):
				print('Error retrieving the length of the sequence {}'.format(scaffold));
				raise;
			if inverse_nn_[0][0]>1:
			# if start position for the first N block is greater than 1.
				block=(1, inverse_nn_[0][0]-1);
				temp_list.append(block);
				text='{}\t{}\t{}\n'.format(scaffold, block[0], block[1]);
				output_dir_.write(text);
			i=1; #1 because we've already looked at 0.
			while i<(len(gff_list_[scaffold])):
			#	for all N blocks except the first and last.
				block=(inverse_nn_[i-1][1], inverse_nn_[i][0]-1);
				#	block is the sequence between the last position of the previous N block and the first position of the current N block
				temp_list.append(block);
				text='{}\t{}\t{}\n'.format(scaffold, block[0], block[1]);
				output_dir_.write(text);
				i+=1;
			#For the last N-block
			block=(inverse_nn_[len(inverse_nn_)-1][1], scaffold_length);
			temp_list.append(block);
			text='{}\t{}\t{}\n'.format(scaffold, block[0], block[1]);
			output_dir_.write(text);
			scaf_blocks_[scaffold]=temp_list;
			#	Adds all of the sequence blocks to a new dictionary to keep.
		else:
			scaffold_length=len(seq_index_[scaffold].seq);
			scaf_blocks_[scaffold]=[(1, scaffold_length)];
	return scaf_blocks_;

def getSequences(scaf_blocks_, seq_index_, output_fasta_):
	from Bio.SeqRecord import SeqRecord;
	for scaffold in scaf_blocks_:
		seq_block=seq_index_[scaffold].seq;
		i=1
		for block in scaf_blocks_[scaffold]:
			new_seq_block=seq_block[block[0]-1:block[1]];
			new_record=SeqRecord(new_seq_block, '{}_{}'.format(seq_index_[scaffold].id, i), '', '');
			i+=1;
			SeqIO.write(new_record, output_fasta_, 'fasta');

###Code.

def main(input_dir, output_dir, fasta_dir, fasta_dir_output=None):
	if input_dir=='-':
		input_dir=sys.stdin;
	if output_dir=='-':
		output_dir=sys.stdout;

	if not isinstance(input_dir, file):
		with open(input_dir, 'rU') as gff_handle:
			gff_list=gffList(gff_handle);
	else:
		gff_list=gffList(input_dir);
	
	seq_index=SeqIO.index(fasta_dir, 'fasta');
	
	if not isinstance(output_dir, file):
		with open(output_dir, 'w') as output_handle:
			scaf_blocks=sequenceBoundaries(gff_list, output_handle, seq_index);
	else:
		scaf_blocks=sequenceBoundaries(gff_list, output_dir, seq_index);

	if fasta_dir_output != None:
		with open(fasta_dir_output, 'w') as outfasta_handle:
			getSequences(scaf_blocks, seq_index, outfasta_handle);



if __name__== '__main__':
	###Argument handling.

	arg_parser = argparse.ArgumentParser(description='');
	arg_parser.add_argument("file_dir", help="Directory to .gff file. To use stdin (for piped input) enter '-'");
	arg_parser.add_argument("fasta_dir", help="Directory to fasta file containing sequences with blocks.");
	arg_parser.add_argument("file_dir_output", help="Directory/name of output file to view blocks in. To use stdout (for piped output) enter '-'.");
	arg_parser.add_argument("-o","--fasta_dir_output", default=None, help="Directory for new fasta file to be written to. If run without this flag, inverseGFF will not generate a new fasta");
	args = arg_parser.parse_args();

	###Variable Definitions

	input_dir=args.file_dir;
	output_dir=args.file_dir_output;

	fasta_dir=args.fasta_dir;
	fasta_dir_output=args.fasta_dir_output;

	main(input_dir, output_dir, fasta_dir, fasta_dir_output);

