from Bio import SeqIO
import math

# Create the range of values for the histogram, and a general dict 
# With the two methods used for predicting genes 

ran = range(0,500000,10)
mgmrange = {}.fromkeys(ran, 0)
prorange = {}.fromkeys(ran, 0)
predictor  = {'mgm': mgmrange, 'pro': prorange}

# Open up the fasta, and check the length. Saving the value in the dict 
with open("data/raw/gc/BBMO.v2_95id_clust.fasta", "r") as handle:
    for record in SeqIO.parse(handle, "fasta"):
      length = len(record.seq)
      leng_floor = math.floor(length * 0.1) * 10
      if "mgm" in record.id:
        predictor["mgm"][leng_floor] +=1
      else:
        predictor["pro"][leng_floor] +=1
        

# Saving results in a simple txt 
outfile = open('data/helper_seqstats.txt', 'w')
 
# Sort out results 
for typ in predictor.keys():
  for key in sorted(predictor[typ].keys()):
    outfile.write("{0}\t{1}\t{2}\n".format(typ, key, predictor[typ][key]))
 

