import sys
from Bio import SeqIO

if len(sys.argv)!=3:
    print "Usage: python find_transcript_synonyms.py path_transcripts ofile"
    exit(1)

def build_Seq_to_transIDs(path_transcripts):
 Seqfile = open(path_transcripts, 'rU')
 Seq2trans = {}
 for rec in SeqIO.parse(Seqfile, "fasta"):
  transID = rec.id.split('|')[0]
  for i in rec.id.split('|'):
   if i.startswith('CDS:'):
    start, stop = map(int, i.strip('CDS:').split('-'))
  Seq_CDS = str(rec.seq)[start-1:stop]
  if Seq2trans.has_key(Seq_CDS):
   Seq2trans[Seq_CDS].extend([transID])
  else:
   Seq2trans[Seq_CDS] = [transID]
 Seqfile.close()
 return Seq2trans

# Main
path_transcripts = sys.argv[1] #"/data/hcui/temp/20170506/Playground/ref/Processed/gencode.v26.pc_transcripts_filter.fa"
Seq2trans = build_Seq_to_transIDs(path_transcripts)

path_out = sys.argv[2] #"/data/hcui/temp/20170506/Playground/ref/Processed/synonym.txt"
file_out = open(path_out, 'w')
for key, value in Seq2trans.items():
 outline = '\t'.join(value) + '\n'
 file_out.write(outline)
file_out.close()
