#!/usr/bin/env python
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def get_mRNA_abd(fn):
	print "reading in mRNA TPM..."
	tid2abd = {}
	abd_file = open(fn)
	line = abd_file.readline()
	for line in abd_file:
		words = line.rstrip().split()
		if(int(words[1])>0):
			tid2abd[words[0]] = int(words[1])
	abd_file.close()
	return tid2abd

def append_abd_to_ref_fa(tid2abd, ifname, ofname):
	print "appending abundance to reference fasta"
	ifile = open(ifname, "rU")
	ofile = open(ofname, 'w')
	for rec in SeqIO.parse(ifile, "fasta"):
		tid = rec.id
		if tid in tid2abd:
			header = "{0}${1}".format(tid, tid2abd[tid])
			new_rec = SeqRecord(rec.seq, id=header, description=rec.description)
			SeqIO.write(new_rec,ofile, "fasta")
	ifile.close()
	ofile.close()

if __name__ == "__main__":
	if len(sys.argv)!=4:
		print "Usage: python make_ref_fa_with_abd.py abd_fname transcript_fasta new_fasta"
		exit(1)
	sfn = sys.argv[1]
	ifname = sys.argv[2]
	ofname = sys.argv[3]
	tid2abd = get_mRNA_abd(sfn)
	append_abd_to_ref_fa(tid2abd, ifname, ofname)
	print "done."
