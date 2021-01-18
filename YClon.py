import os
import time
import sys
import platform
from alive_progress import alive_bar
import pandas as pd
import numpy as np
from collections import defaultdict
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.cluster import AgglomerativeClustering
import hdbscan
import sparse_dot_topn.sparse_dot_topn as ct


def directory_path(file_path):
  OS = platform.system()
  temp =  file_path.split(os.sep)
  file_path = file_path.replace(temp[len(temp)-1],"")
  return(file_path)

def build_kmers_tf_idf(sequence, ksize=3): 
    ngrams = zip(*[sequence[i:] for i in range(ksize)])
    return [''.join(ngram) for ngram in ngrams]


filename =  sys.argv[1]
method = "AHAM"
thr = 0.09
sequence_column = "cdr3"
vcolumn = "v_call"
jcolumn = "j_call"
seqID = "sequence_id"
separator = "\t"

for x in range(1,len(sys.argv)):
  if sys.argv[x].find("--input") != -1:
    filename = sys.argv[x+1]
  elif sys.argv[x].find("--method") != -1: 
    method = sys.argv[x+1]
  elif sys.argv[x].find("--thr") != -1:
    thr = float(sys.argv[x+1])
  elif sys.argv[x].find("--sequence") != -1:
    sequence_column = sys.argv[x+1]
  elif sys.argv[x].find("--v_gene") != -1:
    vcolumn = sys.argv[x+1]
  elif sys.argv[x].find("--j_gene") != -1:
    jcolumn = sys.argv[x+1]
  elif sys.argv[x].find("--seq_id") != -1:
    seqID = sys.argv[x+1]
  elif sys.argv[x].find("--sep") != -1:
    separator = sys.argv[x+1]





start_time = time.time()



path = directory_path(filename)


f = open(filename, 'r')
i=0
fail = 0
clonotypes = {}
file_size = 0
print("Processing file")
for x in f:
  file_size +=1
  if x.find(seqID) == -1:
    x = x.strip()
    data = x.split(separator)
    try:
      cdr3len = str(len(data[junc_indx]))
      vGene = data[vGene_indx].split('*') #include all v gene alleles
      jGene = data[jGene_indx].split('*') #include all v gene alleles
      key = vGene[0]+","+jGene[0]+","+cdr3len
      clonotypes.setdefault(key, [])
      proclonotype = [data[seq_id_indx],data[junc_indx]]
      clonotypes[key].append(proclonotype)
    except:
      fail +=1
  else:
    x = x.strip()
    data = x.split(separator)

    seq_id_indx = data.index(seqID)
    try:
      junc_indx = data.index(sequence_column)
    except:
      print("There is no column named cdr3")
    vGene_indx = data.index(vcolumn)
    jGene_indx = data.index(jcolumn)
    colunas = [data[seq_id_indx],data[junc_indx]]
    # print(colunas)






  
f.close()
temp_filename = path+"YClon_temp.txt"
size = len(clonotypes)
temp = open(temp_filename, 'w')
count = 0
total_clust = 0
unico_pq_VJLen = 0
maior = 0
a = 0
with alive_bar(len(clonotypes), title="Clonotyping") as bar: 
  for key in clonotypes: #each key is the combination of V gene, J gene and the length of cdr3, the values are the sequence ID and cdr3 sequence
    bar()
    if len(clonotypes[key]) > 1:
      teste = pd.DataFrame(clonotypes[key])
      teste.columns = colunas

      junc_seq = teste[sequence_column]
      seq_id = teste[seqID]

      vectorizer = TfidfVectorizer(min_df=1, analyzer=build_kmers_tf_idf)
      tf_idf_matrix = vectorizer.fit_transform(junc_seq)  
      dist = 1 - cosine_similarity(tf_idf_matrix)
      if method ==  "AHAM":
        clusterer = AgglomerativeClustering(distance_threshold = thr, n_clusters= None, affinity="precomputed",linkage='average')
      elif method == "HDBSCAN":
        clusterer = hdbscan.HDBSCAN(min_cluster_size=3, cluster_selection_epsilon= thr, metric='precomputed')
      cluster_teste = clusterer.fit(dist)
      clone = cluster_teste.labels_
      # if cluster_teste.labels_.max() != -1:
      maior = count + cluster_teste.labels_.max() + 1
      for i in range(0, len(clonotypes[key])):
        if clone[i] != -1:
          clone_id = count + int(clone[i]) +1
        else:
          clone_id = maior + 1
          maior += 1
          total_clust += 1
        
        temp.write(str(seq_id[i])+","+str(clone_id)+"\n")	
      count = maior
      total_clust += cluster_teste.labels_.max() + 1	
    else:
      unico_pq_VJLen +=1
      count += 1
      total_clust += 1
      teste = pd.DataFrame(clonotypes[key])
      teste.columns = colunas
      seq_id = teste[seqID]
      for i in range(0, len(clonotypes[key])):
        temp.write(str(seq_id[i])+","+str(count)+"\n")





		

in_temp = open(temp_filename, 'r')
in_airr = open(filename, 'r')
filename_temp = filename.split(".")
out_filename = filename_temp[0]+"_YClon_clonotyped."+filename_temp[1]
out = open(out_filename, 'w+')
clonotipo = {}
maior ={}
maximo = 0


for x in in_temp:
	data = x.split(",")
	clonotipo[data[0]] = data[1]
	if data[1] not in maior:
		maior[data[1]] = list(data[0])
	else:
		maior[data[1]].append(data[0])
for i in maior:
	if maximo < len(maior[i]):
		maximo = len(maior[i])
# print(maximo)


# with alive_bar(file_size, title="Writing output file") as bar: 
with alive_bar(file_size, title="Writing output file") as bar:
  for x in in_airr:
    bar()
    if x.find(seqID) == -1:
      data = x.split(separator)
      if data[0] in clonotipo:
        for i in range(0, len(data)):
          out.write(data[i].strip()+separator)
        out.write(clonotipo[data[0]])
    else:
      out.write(x.strip()+separator+"clone_id\n")



os.remove(temp_filename)




current_time = time.time()
elapsed_time = current_time - start_time


print("The work was completed in: " + "%.3f" % int(elapsed_time))
