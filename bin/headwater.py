#!/usr/local/bin/env python
#python 2.7.5 requires biopython

########### Headwater ############

#A tool that takes a GFF file and a genome file, and returns a file with x nucleotides upstream of each gene.


#Version 1. Darcy Jones, December 2013.
#Contact Darcy Jones, darcy.ab.jones@gmail.com

###Import modules
import sys;
import argparse;
from Bio import SeqIO;
from Bio.Seq import Seq;
from Bio.SeqRecord import SeqRecord;
import os;

###Function definitions.

def gffList(gff_handle_):
#   Reads through all lines in .gff file, extracts seqid, start and stop positions, returns dictionary keyed by seqid, with value = list of dictionaries with format={id:'id', name:'name', start:'start', stop:'stop'};
    temp_dict={};
    i=0
    for line in gff_handle_:
        record=line.rstrip('\n').split('\t');
        #   Removes newline characters, splits into list by tabs.
        if len(record)>1:
        #   Excludes first line.
            attributes_list=record[8].split(';')
            attributes={key.lower():value for key, value in [attr.split('=') for attr in attributes_list if '=' in attr]};
            #   Column 9 in the GFF file contains ';' separated attributes of the format: key=val. The above script first strips by ';' then by '=' to create a dictionary with the same key:val setup. keys will always be in lower case.
            if len(attributes)>0:
                if 'id' in attributes:
                    seq_id=attributes['id'];
                elif 'name' in attributes:
                    seq_id=attributes['id'];
                else:
                    seq_id=record[0]+'gffLine'+repr(i);

                if 'name' in attributes:
                    seq_name=attributes['name']
                else:
                    seq_name=seq_id;

                seq_desc=' '.join([attributes[attr] for attr in attributes if attr not in {'name', 'id'}])
            else:
                seq_id=record[0]+'gffLine'+repr(i);
                seq_name=seq_id;
                seq_desc='';
            if record[0] not in temp_dict:
                # First dict entry for a scaffold. if seqid not in temp_dict
                #'strand' returns TRUE if not negative strand, incuding all other characters.
                temp_dict[record[0]]=[{'id':seq_id,'name':seq_name,'description':seq_desc,'start':int(record[3]),'stop':int(record[4]),'strand':record[6]!='-'}];
                #   record[0] is seqid, record[3] is start, record[4] is stop.
            else:
                # Continues adding (start, stop) tuples for each scaffold (seqid).
                #'strand' returns TRUE if not negative strand, incuding all other characters.
                temp_dict[record[0]].append({'id':seq_id,'name':seq_name,'description':seq_desc,'start':int(record[3]),'stop':int(record[4]),'strand':record[6]!='-'});
                # Record[0] is seqid, record[3] is start, record[4] is stop.
        i+=1;
    for sequence in temp_dict:
        temp_dict[sequence].sort(key=lambda dic: dic['start']); # Sort in place by start position, small to big

    return temp_dict;

def extendSeq(gff_dict_, number_, last_gene):
    if last_gene==True:
        for seqid in gff_dict_: #For each scaffold
            i=0;        
            while i < len(gff_dict_[seqid]): #For each gene
                if gff_dict_[seqid][i]['strand']: #If strand equals TRUE i.e. not negative
                    if i<1:
                        promoterPos=0;
                    else:
                        promoterPos=gff_dict_[seqid][i-1]['stop'];
                    if number_!=None:
                        if promoterPos<(gff_dict_[seqid][i]['start']-int(number_)):
                            promoterPos=(gff_dict_[seqid][i]['start']-int(number_))
                else:
                    if i+1>=len(gff_dict_[seqid]):
                        promoterPos=None; # Takes everything until the end of the scaffold.
                    else:
                        promoterPos=gff_dict_[seqid][i+1]['start'];
                    if number_!=None:
                        if promoterPos>(gff_dict_[seqid][i]['stop']+int(number_)) or promoterPos==None:
                            promoterPos=(gff_dict_[seqid][i]['stop']+int(number_));
                gff_dict_[seqid][i]['promoterPos']=promoterPos;
                i+=1;           
    else:
        for seqid in gff_dict_: #For each scaffold
            i=0;
            while i < len(gff_dict_[seqid]): #For each gene
                if gff_dict_[seqid][i]['strand']: 
                    promoterPos=gff_dict_[seqid][i]['start'] - (int(number_)+1)
                else:
                    promoterPos = gff_dict_[seqid][i]['stop'] + int(number_)
                if promoterPos < 0:
                    promoterPos=0;
                gff_dict_[seqid][i]['promoterPos']=promoterPos;
                i+=1;
    return gff_dict_;

def formatFinder(file_):
#   Attempts to find sequence format by file extension.
    extToFormat_={
    '.abi':'abi','.ab':'abi','.ace':'ace','.ebl':'embl',
    '.emb':'embl','.dat':'embl','exclude':'exclude','.fasta':'fasta',
    '.fast':'fasta','.seq':'fasta','.fa':'fasta','.fsa':'fasta',
    '.nt':'fasta','.aa':'fasta','.faa':'fasta',
    '.fastq':'fastq','.fq':'fastq','.gb':'genbank','.gbank':'genbank',
    '.genbank':'genbank','.ig':'ig','.imgt':'imgt',
    '.xml':'seqxml','.sff':'sff',
    '.sp':'swiss', '.tab':'tab','.qual':'qual' 
    };

    f_ext=os.path.splitext(file_);
    #   creates a list of strings, separated by '.'. Used to find 
    #   file extensions.
    #   if there is at least one string separated by a '.' in the
    #   filename.

    if f_ext[1] in extToFormat_:
    #   standard extensions
        format = extToFormat_[f_ext[1]];
        #   uses outer most extension to predict input format
    else:
        format ='exclude';
    return format;

###Code.

def main(gff_input, seq_input, number=500, seq_format=None, seq_output='-', verbose=0, last_gene=False):
    if verbose>0:
        print('######### Beginning upstreamGetter ###########');

    if seq_format==None:
        seq_format=formatFinder(seq_input);
    if seq_format=='exclude':
        sys.exit('The format for the sequence file {} is not supported or has a non-standard extension. Please try running again and specify the -f (format) flag.'.format(seq_input));

    if verbose>1:
        print('Using parameters');
        print('\tGFF file: {}'.format(gff_input));
        print('\tSequence file: {}'.format(seq_input));
        print('\tRead/Write format: {}'.format(seq_format));
        print('\tNumber: {}'.format(number));
        print('\tWriting sequences to: {}'.format(seq_output));

    if verbose>0:
        print('### Reading GFF data ###');
    with open(gff_input, 'rU') as gff_handle:
        gff_dict=gffList(gff_handle);
    if verbose>1:
        print('Found {} seqids in GFF'.format(len(gff_dict)));

    if verbose>0:
        print('### Finding upstream region to be written ###');
    gff_dict=extendSeq(gff_dict, number, last_gene);

    seq_index=SeqIO.index(seq_input, seq_format);

    if verbose>0:
        print('### Extracting upstream regions ###');
    seq_write=[];
    no_genes=[];
    for seq in seq_index:
        scaffold=seq_index[seq]
        if seq in gff_dict:
            if verbose>1:
                print('Extracting upstream regions for genes in sequence {}'.format(seq))
            for gene in gff_dict[seq]:
                if verbose >1:
                    if gene['strand']:
                        upstreamRange='{}-{}'.format(gene['promoterPos'], gene['start'])
                    else:
                        if gene['promoterPos']==None:
                            upstreamRange='{}-END'.format(gene['stop'])
                        else:
                            upstreamRange='{}-{}'.format(gene['stop'], gene['promoterPos'])
                    print("\tGene: {}, Upstream sequence Range: {}, Gene sequence Range: {}-{}".format(gene['id'], upstreamRange, gene['start'], gene['stop']));
                if gene['strand']: #Strand validate TRUE if NOT negative stranded
                    record=SeqRecord(scaffold.seq[gene['promoterPos']:gene['start']-1], id=gene['id'], name=gene['name'], description=gene['description'])
                else:
                    record=SeqRecord(scaffold.seq[gene['stop']:gene['promoterPos']], id=gene['id'], name=gene['name'], description=gene['description'])
                    record.seq=record.seq.reverse_complement();
                seq_write.append(record);
        else:
            no_genes.append(seq);

    if verbose>1 and len(no_genes)>0:
        print('No GFF data available for {} sequences'.format(len(no_genes)));
    if verbose>0:
        print('### Writing new sequences ###');
    count=SeqIO.write(seq_write, seq_output, seq_format);

    if verbose>0:
        print('Successfully wrote {} sequences to {}'.format(count, seq_output));

    if verbose>0:
        print('######### End upstreamGetter ###########');

if __name__== '__main__':
    ###Argument handling.
    arg_parser = argparse.ArgumentParser(description='');
    arg_parser.add_argument("gff_dir", help="Directory to .gff file. To use stdin (for piped input) enter '-'");
    arg_parser.add_argument("seq_dir", help="Directory to sequence file containing genome");
    arg_parser.add_argument("-n","--number", default=None, help="Number of upstream nucleotides to include in the output. Default=500");
    arg_parser.add_argument("-l","--last_gene", action='store_true', default=False, help="Instead of extracting a number of upstream nucleotides, extract all nucleotides from end of last gene to beginning of current gene. Number is ignored if this is flagged.");
    arg_parser.add_argument("-f","--seq_format", default=None, help="Sequence format to read and write the sequences in. By default the program will attempt to find the format.");
    arg_parser.add_argument("-o","--seq_output", default='-', help="Directory for new sequence file to be written to. Default is '-' (stdout)");
    arg_parser.add_argument("-v", '--verbose', action="count", help="Toggle counter. -v gives limited running feedback. -vv gives lots of running feedback. Default (no -v flag) does not give any running feedback. When using stdout output, verbose is automatically disabled even if flagged.");
    args = arg_parser.parse_args();

    ###Variable Definitions
    verbose=args.verbose;
    gff_input=args.gff_dir;
    if gff_input=='-':
        gff_input=sys.stdin;
    seq_output=args.seq_output;
    if seq_output=='-' or seq_output==sys.stdout:
        seq_output=sys.stdout;
        verbose=0;
    seq_input=args.seq_dir;
    number=args.number;
    seq_format=args.seq_format;
    last_gene=args.last_gene;
    if not last_gene and number == None:
        number=500;

    main(gff_input, seq_input, number, seq_format, seq_output, verbose, last_gene);

