#!/usr/bin/env python

# view synth_riboseq_005_filtered_sequence_transcript_Aligned.out.bam | awk '{print $1"\t"$3"\t"$4}' | awk -F"|" '{printf $1"\t"$2; for(i=8;i<=NF-1;i++){if($i~/CDS:/) {gsub("CDS:","",$i); gsub("-","\t",$i); printf "\t"$i;}} printf "\t"$NF"\n"}'

import sys
import re
import string
from Bio import SeqIO
import pysam
from BCBio.GFF import GFFExaminer
from translation import *

# Give a transcript reference fasta file, and return a dictionary from transcript ID to CDS range.
def get_gencode_cds_range(transcript_ref_file):
	"""
	transcript.fa header:
	0 transcript-id|
	1 gene-id|
	2 Havana-gene-id (if the gene contains manually annotated transcripts, '-' otherwise)|
	3 Havana-transcript-id (if this transcript was manually annotated, '-' otherwise)|
	4 transcript-name|
	5 gene-name|
	6 sequence-length|
	7 5'-UTR (3'-UTR if reverse strand) location in the transcript|
	8 CDS location in the transcript|
	9 3'-UTR (5'-UTR if reverse strand) location in the transcript
	"""
	f_ref = open(transcript_ref_file)
	trans2cds = {}
	for line in f_ref:
		if line.startswith('>'):
			line = line.lstrip('>').rstrip().rstrip('|')
			linearray = line.split('|')
			transID = linearray[0]
			for a in linearray:
				if a.startswith("CDS:"):
					start, stop = a.lstrip("CDS:").split('-')
					trans2cds[transID] = [string.atoi(start), string.atoi(stop)]
	f_ref.close()
	return trans2cds

# Give an gencode gtf annotation, and return a dictionary from transcript ID to all exon genome absolute positions.
def transcript2exons(gencode_annotation_file):
	f_annotation = open(gencode_annotation_file)
	trans2exons = {}
	for line in f_annotation:
		if line.startswith("#"):
			continue
		linearray = line.strip().split('\t')
		if linearray[2] == 'CDS' and (linearray[8].find("exon_number") and linearray[8].find("transcript_id") != -1):
			description = linearray[8].strip(';').split('; ')
			description = [description[i].replace('"','').split(' ') for i in range(0,len(description))]
			for a in description:
				aname, aval = a
				if aname == "transcript_id":
					transID = aval
				if aname == "exon_number":
					exonID = aval
			if trans2exons.has_key(transID):
				trans2exons[transID].append([linearray[0], linearray[6], string.atoi(exonID), string.atoi(linearray[3]), string.atoi(linearray[4])])
			else:
				trans2exons[transID] = [[linearray[0], linearray[6], string.atoi(exonID), string.atoi(linearray[3]), string.atoi(linearray[4])]]
	f_annotation.close()
	return trans2exons

# Check whether the information from exons are correct.
def checkexons(trans2exons, trans2cds):
	for key, value in trans2cds.items():
		if not trans2exons.has_key(key):
			print "Error: " + key
			continue
		len_from_exons = sum([trans2exons[key][i][4]-trans2exons[key][i][3]+1 for i in range(0, len(trans2exons[key]))])
		len_from_CDS = trans2cds[key][1]-trans2cds[key][0]+1
		if len_from_CDS - len_from_exons <> 0 and len_from_CDS - len_from_exons <> 3:
			print "Error: " + key
			continue

# 
def get_absolute_position(align_result_file, trans2exons, offset, output_file):
	f_alignfile = open(align_result_file)
	f_out = open(output_file, 'w')
	for line in f_alignfile:
		linearray = line.strip().split('\t')
		this_readname = linearray[0]
		this_transcriptname = linearray[1]
		this_genename = linearray[2]
		this_CDSstart = string.atoi(linearray[3])
		this_CDSstop = string.atoi(linearray[4])
		this_mapstart = string.atoi(linearray[5])
		this_readlen = string.atoi(linearray[6])
		this_exon_pos = trans2exons[this_transcriptname]
		this_ribosome_pos = this_mapstart + offset[this_readlen]
		if this_ribosome_pos < this_CDSstart:
			if this_ribosome_pos < this_CDSstart-3:
				continue
			else:
				this_ribosome_pos = this_CDSstart
		if this_ribosome_pos > this_CDSstop:
			if this_ribosome_pos > this_CDSstop+3:
				continue
			else:
				this_ribosome_pos = this_CDSstop
		this_ribosome_pos_in_CDS = this_ribosome_pos - this_CDSstart + 1
		pos_in_exon = this_ribosome_pos_in_CDS
		for i in range(0, len(this_exon_pos)):
			exon_info = this_exon_pos[i]
			exon_length = exon_info[4] - exon_info[3] + 1
			if pos_in_exon > exon_length:
				pos_in_exon = pos_in_exon - exon_length
			else:
				if exon_info[1] == '+':
					absolute_pos = exon_info[3] + pos_in_exon - 1
				if exon_info[1] == '-':
					absolute_pos = exon_info[4] - pos_in_exon + 1
				break
		f_out.write('\t'.join(map(str, [this_readname, this_exon_pos[0][0], this_exon_pos[0][1], absolute_pos, this_ribosome_pos_in_CDS, this_transcriptname, this_genename])) + '\n')
	f_out.close()

	
# Main
path_annotation = sys.argv[1]
trans2exons = transcript2exons(path_annotation)

path_ref = sys.argv[2]
trans2cds = get_gencode_cds_range(path_ref)

checkexons(trans2exons, trans2cds)

offset = {25:12,26:12,27:12,28:12,29:12,30:12,31:13,32:13,33:13,34:14,35:14,36:14}

path_in = sys.argv[3]
path_out = sys.argv[4]
get_absolute_position(path_in, trans2exons, offset, path_out)
