#!/usr/bin/env python
import sys
from Bio import SeqIO

if len(sys.argv)!=3:
    print "Usage: python select_rrna_from_ncrna.py rrna_fa ofile"
    exit(1)

rrna_fa = sys.argv[1]
out_fa = sys.argv[2]

ofile = open(out_fa, 'w')

j = 0
rrna_file = open(rrna_fa)
for rec in SeqIO.parse(rrna_file, "fasta"):
    theader = rec.description
    words = theader.split()
    for w in words:
        if w.startswith("gene_biotype:"):
            gbt = w.lstrip("gene_biotype:")
            if gbt == "rRNA":
                SeqIO.write(rec, ofile, "fasta")
                j += 1
rrna_file.close()
print "{0} rRNA seqs included".format(j)
ofile.close()

