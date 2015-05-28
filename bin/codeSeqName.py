#!/usr/local/bin/env python
#python 2.7.5 requires biopython

############### codeSeqName.py ##################

#A sequence formatting tool that can take sequence file(s) or strings 
# and code the sequence names and information, putting out a new sequence
# file with coded names, and a tab delimited file mapping the original name to the coded name.
# You can then run the code again using any file containing the coded names and the map file to restore the names.

#Version 1. Darcy Jones, December 2013.
#Contact Darcy Jones, darcy.ab.jones@gmail.com

###Import modules
from Bio import SeqIO, AlignIO;
import sys;
import argparse;
import os;
import random;
from collections import defaultdict;

def codemapNameCheck(suffix_):
    if os.path.isfile('{}.codemap'.format(suffix_)):
        if isint(suffix_[-1]) and suffix_[-2]=='_':
            i=int(suffix_[-1])+1;
            suffix_=suffix_[-1]+str(i);
        else:
            suffix_=suffix_+'_1';
    return suffix_;

def getDir(dir_, extToFormat_, permitted_):
#   Takes a directory and outputs a list of files. 
    file_list_=[];
    try:
        file_list_=[os.path.join(dir_, f) for f in os.listdir(dir_) 
        if os.path.isfile(os.path.join(dir_, f)) and f[0]!='.' 
        and extToFormat_[os.path.splitext(f)[1]] in permitted_];
        #   Returns a list of .phy files found within 
        #   the folder specified by the input flag.
        #   Excludes hidden files eg .ds-store and the specified output file.
    except:
        error='Error. A problem occurred trying to open the directory \
        {}. Check that you have permission to use this directory \
        and/or that the directory has the correct filesystem \
        format'.format(dir_);
        print(error);
    return file_list_;

def formatFinder(file_, extToFormat_, permitted_):
#   Attempts to find sequence format by file extension.
    f_ext=os.path.splitext(file_);
    #   creates a list of strings, separated by '.'. Used to find 
    #   file extensions.
    #   if there is at least one string separated by a '.' in the
    #   filename.

    if extToFormat_[f_ext[1]] in permitted_: 
    #   standard extensions
        format = extToFormat_[f_ext[1]];
        #   uses outer most extension to predict input format
    else:
        format ='exclude';
    return format;

def isint(s):
    try:
        int(s);
        return True;
    except ValueError:
        return False;

def listHandler(input_dir_, output_dir_, format_, suffix_, extToFormat_, formatToExt, permitted_):
#   Takes the list of input directories, checks that they exist, and that they have a valid format, and handles folders in the list, by adding each supported file to the final dictionary of input files.
    file_list_=defaultdict(list);
    i=0;
    j=0;
    while j<len(input_dir_):
    #   for j in input_dir:
        if i>=len(output_dir_):
        # if there are still more input files but the specified list of output files is finished, start using default output.
            output_dir_.append(None);
        elif output_dir_[i]==sys.stdout:
            pass;
        elif output_dir_[i]=='-':
            output_dir_[i]=sys.stdout
        elif output_dir_[i][-1]=='/':
        #   If you've specified a folder as output.
            try:
                os.mkdir(output_dir_[i]); 
                #   Attempts to create the specified output directory
            except OSError:
            #   Raised if the directory already exists.
                pass;

        if len(output_dir_)==1:
        #   If there is only one specified output path and it is not the default output type, extend the list so that all files are written to that output, eg. write to folder, stdout, or merge to one file..
            if output_dir_[0] == '+' or output_dir_[0] == None:
                output_dir_=['+' for input_ in input_dir_];
            elif output_dir_[0] == sys.stdout or output_dir_[0] == '-':
                output_dir_=[output_dir_[0] for input_ in input_dir_];
            elif os.path.isdir(output_dir_[0]):
                output_dir_=[output_dir_[0] for input_ in input_dir_];

        if input_dir_[j]=='-':
            input_dir_[j]=sys.stdin;

        if output_dir_[i] == '+' or output_dir_[i] == None:
        #   if output_dir is default (write to same directory as input with same name, but different format);
            if input_dir_[j] == sys.stdin:
                if format_ == None:
                    _format_=formatFinder(input_dir_[j], extToFormat_, permitted_);
                    #   find input format.
                else:
                    _format_=format_;
                file_list_['stdin'+j+suffix_+formatToExt[_format_]].append({'in_file':input_dir_[j],'format':_format_});
                    #   NB in this case for stdin, the output is written as stdin1_suffix.format, stdin2_suffix.format etc. 
            
            elif os.path.isdir(input_dir_[j]):
            #   If the input item is a directory get all permitted items in the directory (getDir()). 
                for _file_ in getDir(input_dir_[j], extToFormat_, permitted_):
                #   Loop through the files found in the directory.
                    if format_ == None:
                        _format_=formatFinder(_file_, extToFormat_, permitted_);
                        #   find input format.
                    else:
                        _format_=format_;
                    file_list_[os.path.splitext(_file_)[0]+suffix_+formatToExt[_format_]].append({'in_file':_file_,'format':_format_});
            
            elif os.path.isfile(input_dir_[j]):
            #   If the input item is a file.
                if format_ == None:
                    _format_=formatFinder(input_dir_[j], extToFormat_, permitted_);
                    #   find input format.
                else:
                    _format_=format_;
                file_list_[os.path.splitext(input_dir_[j])[0]+suffix_+formatToExt[_format_]].append({'in_file':input_dir_[j],'format':_format_});

        elif output_dir_[i] == sys.stdout:
            if input_dir_[j] == sys.stdin:
                if format_ == None:
                    _format_=formatFinder(input_dir_[j], extToFormat_, permitted_);
                    #   find input format.
                else:
                    _format_=format_;
                file_list_[output_dir_[i]].append({'in_file':input_dir_[j],'format':_format_});

            elif os.path.isdir(input_dir_[j]):
                for _file_ in getDir(input_dir_[j], extToFormat_, permitted_):
                    if format_ == None:
                        _format_=formatFinder(_file_, extToFormat_, permitted_);
                        #   find input format.
                    else:
                        _format_=format_;
                    file_list_[output_dir_[i]].append({'in_file':_file_,'format':_format_});

            elif os.path.isfile(input_dir_[j]):
                in_format_=formatFinder(input_dir_[j], extToFormat_, permitted_);
                for out_format_ in output_format_:
                    file_list_[output_dir_[i]].append({'in_file':input_dir_[j],'format':in_format_});

        else:
        #   if the output is to be written to a directory or a file.
            if input_dir_[j] == sys.stdin and os.path.isdir(output_dir_[i]):
            #   If the input is from stdin and output is to a directory.
                if format_ == None:
                    _format_=formatFinder(input_dir_[j], extToFormat_, permitted_);
                    #   find input format.
                else:
                    _format_=format_;
                file_list_[os.path.join(output_dir_[i], 'stdin'+j+suffix_+formatToExt[_format_])].append({'in_file':input_dir_[j],'format':_format_});
                    #   NB in this case for stdin, the output is written as stdin1.format, stdin2.format etc. 

            elif input_dir_[j] == sys.stdin and (output_dir_[i] != None or output_dir_[i] != '+'):
            #   If input is from stdin and output is to a file
                if format_ == None:
                    _format_=formatFinder(input_dir_[j], extToFormat_, permitted_);
                    #   find input format.
                else:
                    _format_=format_;
                file_list_[output_dir_[i]+formatToExt[_format_]].append({'in_file':input_dir_[j],'format':_format_});

            elif os.path.isfile(input_dir_[j]) and os.path.isdir(output_dir_[i]):
            #   If input is from a file and output is to a directory. 
                if format_ == None:
                    _format_=formatFinder(input_dir_[j], extToFormat_, permitted_);
                    #   find input format.
                else:
                    _format_=format_;
                file_list_[os.path.join(output_dir_[i], os.path.splitext(os.path.split(input_dir_[j])[1])[0]+suffix_+formatToExt[_format_])].append({'in_file':input_dir_[j],'format':_format_});
            
            elif os.path.isdir(input_dir_[j]) and os.path.isdir(output_dir_[i]):
            #   If input is batch from a directory and output is to a directory.
                for _file_ in getDir(input_dir_[j], extToFormat_, permitted_):
                #   Get all of the supported files from the input directory and loop through them.
                    if format_ == None:
                        _format_=formatFinder(_file_, extToFormat_, permitted_);
                        #   find input format.
                    else:
                        _format_=format_;
                    file_list_[os.path.join(output_dir_[i], os.path.splitext(os.path.split(_file_)[1])[0]+suffix_+formatToExt[_format_])].append({'in_file':_file_,'format':_format_});
                    #   Nb output file name take the input file, splits off the original extension and path information and tacks on the new path and extension.

            elif os.path.isfile(input_dir_[j]) and output_dir_[i] != None:
            #   If the input is a file and output is to a file.
                if format_ == None:
                    _format_=formatFinder(input_dir_[j], extToFormat_, permitted_);
                    #   find input format.
                else:
                    _format_=format_;
                file_list_[output_dir_[i]+formatToExt[_format_]].append({'in_file':input_dir_[j],'format':_format_});
            elif os.path.isdir(input_dir_[j]) and output_dir_[i] != None:
            #   If the input is a directory and output is to a file.
                for _file_ in getDir(input_dir_[j], extToFormat_, permitted_):
                    if format_ == None:
                        _format_=formatFinder(_file_, extToFormat_, permitted_);
                        #   find input format.
                    else:
                        _format_=format_;
                    file_list_[output_dir_[i]+formatToExt[_format_]].append({'in_file':_file_,'format':_format_});
        j+=1;
        i+=1;
    return file_list_;

def codeGenerator(length, seen_set, chars='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'):
#   Cred for this lovely bit, http://stackoverflow.com/questions/2257441/python-random-string-generation-with-upper-case-letters-and-digits
    code= ''.join(random.choice(chars) for x in range(length));
    if code in seen_set:
    #   If we've already used the code for a different sequence run the generator again.
        codeGenerator(length, seen_set);
    return code;

def taxonFinder(record_):
#   Query NCBI database for taxonomic information, return everythin it can use dots for what it can't.
    taxon_={'super_kingdom':'.','kingdom':'.', 'phylum':'.', 'class':'.','sub_class':'.','order':'.','sub_order':'.','family':'.','sub_family':'.','genus':'.','sub_genus':'.','species':'.','sub_species':'.'};
    return taxon_

### Format Dictionaries
extToFormat={
    '.abi':'abi','.ab':'abi','.ace':'ace','.aln':'clustal','.ebl':'embl',
    '.emb':'embl','.dat':'embl','exclude':'exclude','e':'e','.fasta':'fasta',
    '.fast':'fasta','.seq':'fasta','.fa':'fasta','.fsa':'fasta',
    '.nt':'fasta','.aa':'fasta','.fasta_aln':'fasta','.faa':'fasta',
    '.fastq':'fastq','.fq':'fastq','.gb':'genbank','.gbank':'genbank',
    '.genbank':'genbank','.ig':'ig','.imgt':'imgt','.nex':'nexus',
    '.nxs':'nexus','.phd':'phred','.phy':'phylip','.phylip':'phylip',
    '.phlp':'phylip','.phyl':'phylip','.ph':'phylip','.pir':'pir',
    '.xml':'seqxml','.sff':'sff','.stk':'stockholm','.sth':'stockholm',
    '.sp':'swiss', '.tab':'tab','.qual':'qual' 
};

formatToExt={
    'abi':'.abi','abi':'.ab','ace':'.ace','clustal':'.aln','embl':'.ebl',
    'exclude':'exclude','e':'e','fasta':'.fa','fastq':'.fastq','genbank':'.gb',
    'ig':'.ig','imgt':'.imgt','nexus':'.nex','paup':'.nex','phred':'.phd',
    'phylip':'.phy','pir':'.pir','seqxml':'.xml','sff':'.sff',
    'stockholm':'.stk','swiss':'.sp','tab':'.tab','qual':'.qual',
    'uniprot-xml':'.xml','seqxml':'.xml','fastq-solexa':'.fastq',
    'fastq-illumina':'.fastq', 'fastq-sanger':'.fastq'
};

permitted_formats={
    'clustal', 'embl', 'fasta', 'fastq', 'fastq-sanger', 'fastq-solexa', 
    'fastq-illumina', 'genbank','imgt', 'nexus', 'phd', 'phylip', 
    'phylip-sequential', 'phylip-relaxed' 'seqxml', 'sff',
    'stockholm', 'tab, qual', 'exclude'
};

### Code

def main(input_dir, output_dir, format=None, code_map_dir=None, suffix=None, truncate=9):
    if code_map_dir == None:
        if suffix==None:
            suffix='_coded';
        seen=set();
        file_dict=listHandler(input_dir, output_dir, format, suffix, extToFormat, formatToExt, permitted_formats);
        suffix=codemapNameCheck(suffix);
        with open('{}.codemap'.format(suffix), 'w') as map_file:
        #   Open the file to write the code map to.
            for out_file_name in file_dict:
                format_=file_dict[out_file_name][0]['format'];
                out_sequences=[];
                if format_ != 'exclude':
                    with open(out_file_name, 'w') as out_file:
                        for in_file_name in file_dict[out_file_name]:
                            sequences=SeqIO.parse(in_file_name['in_file'], in_file_name['format']);
                            for sequence in sequences:
                                if format != 'exclude':
                                    code_name=codeGenerator(truncate, seen);
                                    current_id=sequence.id;
                                    current_name=sequence.name;
                                    current_description=sequence.description;
                                    taxon=taxonFinder(sequence);
                                    map_file.write('{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n'.format(code_name, current_id, current_name, current_description, taxon['super_kingdom'], taxon['kingdom'], taxon['phylum'], taxon['class'], taxon['sub_class'], taxon['order'], taxon['sub_order'], taxon['family'], taxon['sub_family'], taxon['genus'], taxon['sub_genus'], taxon['species'], taxon['sub_species']));
                                    #   Writes line to tab delimited file.
                                    sequence.id=code_name;
                                    sequence.name='';
                                    sequence.description='';
                                    #   Overwrites the sequence descriptions with new coded name.
                                    seen.add(code_name);
                                    out_sequences.append(sequence);
                        SeqIO.write(out_sequences, out_file, format_);

    elif os.path.isfile(code_map_dir):
        if suffix==None:
            suffix='_decoded';
        map_dict=dict();
        with open(code_map_dir, 'r') as map_file:
            for line in map_file:
                record=line.rstrip('\n').split('\t');
                map_dict[record[0]]={'id':record[1], 'name':record[2], 'description':record[3]};
        file_dict= {output_dir[0]:[{'in_file':input_dir[0]}]}#listHandler(input_dir, output_dir, format, suffix, extToFormat, formatToExt, permitted_formats);
        for out_file_name in file_dict:
            with open(out_file_name, 'w') as out_file_:
                for in_file_name in file_dict[out_file_name]:
                    with open(in_file_name['in_file'], 'rU') as in_file_:
                        for line in in_file_:
                            for key, val in map_dict.iteritems():
                                line=line.replace(key, val['id']);
                            out_file_.write(line);

    else:
        print('Error, the code_map_dir specified ({}) is not a file'.format(code_map_dir));

if __name__ == '__main__':

    ### Argument handling
    arg_parser=argparse.ArgumentParser(description='A sequence formatting tool that can take sequence file(s) or strings and code the sequence names and information, putting out a new sequence file with coded names, and a tab delimited file mapping the original name to the coded name. You can then run the code again using any file containing the coded names and the map file to restore the names.');
    arg_parser.add_argument("-f", '--input_file_dir', nargs='*', default=['-'], 
        help="directory to input file, default is (-) stdin."
    );
    #directory to input file(s) and directories, default is (-) stdin. You may indicate multiple files or directories by separating with a space eg. -i one.fa dir1/ two.fa etc. If your input is a folder of files all files in that directory with a supported extension will be used.

    arg_parser.add_argument("-o", '--output_file_dir', nargs='*', 
        default=['-'], help="directory to output file, default is (-) stdout."
    );
    # directory to output file(s) or directories, default is (-) stdout. You may indicate multiple output directories by separating with a space eg. -i one.faa two.pir etc

    arg_parser.add_argument("-g", '--format', default=None, help="Format to read write the sequence file as when creating the coded file. By default codeSeqName will Attempt to find the format automatically from the file extension.", choices=['clustal', 'embl', 'fasta', 'fastq', 'fastq-sanger', 'fastq-solexa', 'fastq-illumina', 'genbank','imgt', 'nexus', 'phd', 'phylip', 'phylip-sequential', 'phylip-relaxed', 'seqxml', 'sff','stockholm', 'tab, qual', 'exclude']
    );

    arg_parser.add_argument("-c", '--code_map_dir', default=None, help="Directory to code map file. Only one code map file is permitted at a time. If this argument is ommitted, codeSeqName will create a new code file etc."
    );

    arg_parser.add_argument("-s", '--suffix', default=None, help="Specify a short extension to be added to the end of the sequence file name to help distinguish the coded or decoded file from the original. This name will also be the name of the code map file. default =filename_coded.fa, _coded.codemap"
    );

    arg_parser.add_argument(
        "-t", '--truncate', type=int, default=9, help="Enter an integer.\
        Keeps coded sequence names to at or below x alphanumeric characters. Default= 9"
    );

    args = arg_parser.parse_args();

    ###Variables
    input_dir=args.input_file_dir; 
    #   Directory to input file(s)

    output_dir=args.output_file_dir; 
    #   Directory to output file(s)

    format=args.format;

    code_map_dir=args.code_map_dir; 
    #   Directory to .codemap file

    suffix=args.suffix; 
    #   Suffix to add to end of seq file, and to be name for code map file

    truncate=args.truncate;
    #   Stores an integer to truncate sequence name length to.

    main(input_dir, output_dir, format, code_map_dir, suffix, truncate);
