#!/usr/bin/env python3
import sys

version = "1.4 - May 2nd 2023"
try:
	if sys.argv.index("--version") != -1 :
		print("Yclon version "+version)
		exit()
except:
	pass

import os
import time
start_time = time.time()
import platform
from alive_progress import alive_bar
import pandas as pd
import numpy as np
from collections import defaultdict
# from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.cluster import AgglomerativeClustering


def directory_path(file_path):
	OS = platform.system()
	temp =	file_path.split(os.sep)
	file_path = file_path.replace(temp[len(temp)-1],"")
	return(file_path)

def simpson_di(data):
	#from https://gist.github.com/martinjc/f227b447791df8c90568


    """ Given a hash { 'species': count } , returns the Simpson Diversity Index
    
    >>> simpson_di({'a': 10, 'b': 20, 'c': 30,})
    0.3888888888888889
    """

    def p(n, N):
        """ Relative abundance """
        if n == 0:
            return 0
        else:
            return float(n)/N

    N = sum(data.values())
    
    return sum(p(n, N)**2 for n in data.values() if n != 0)


def shannon_di(data):
    #from https://gist.github.com/audy/783125

    """ Given a hash { 'species': count } , returns the SDI
    
    >>> sdi({'a': 10, 'b': 20, 'c': 30,})
    1.0114042647073518"""
    
    from math import log as ln
    
    def p(n, N):
        """ Relative abundance """
        if n == 0:
            return 0
        else:
            return (float(n)/N) * ln(float(n)/N)
            
    N = sum(data.values())
    
    return -sum(p(n, N) for n in data.values() if n != 0)

separator = "\t"

if("rarefy" in sys.argv):
	from rare import rarefy

	clones = {}
	filename = sys.argv[2]
	print(filename)
	
	try:
		repertoire = pd.read_csv(filename,usecols = ['clone_id','clone_seq_count'], sep=separator)
	except:
		print("Please provide the path to a file with clone_id and clone_seq_count columns")
		exit()
	repertoire = repertoire.drop_duplicates(subset=["clone_id"])

	# print(repertoire["clone_id"].to_numpy())
	# exit()
	# print(repertoire.T[repertoire.T[0]=="clone_seq_count"].to_numpy)
	repertoire = repertoire.T
	a = []
	for x in range(len(repertoire.columns)):
		a.append("col"+str(x))
	
	repertoire.columns = a
	import random
	# path = 'https://github.com/Auerilas/ecopy-data/raw/master/data/BCI.csv'
	# downloadPath = path.format()
	# loaded = pd.read_csv(downloadPath)
	# print(sum(repertoire.loc["clone_seq_count"]))
	repertoire.loc["clone_seq_count"] = [random.randint(0,2) for x in range(1,24000)]
	repertoire.loc["clone_id"] = [random.randint(0,2) for x in range(1,24000)]
	# print(repertoire)
	# print(sum(loaded.iloc[1]))
	print(rarefy(repertoire.to_numpy(),'rarefy'))
	exit()

if("diversity" in sys.argv):
	clones = {}
	filename = sys.argv[2]
	print(filename)
	try:
		repertoire = pd.read_csv(filename,usecols = ['clone_id','clone_seq_count'], sep=separator)
	except:
		print("Please provide the path to a file with clone_id and clone_seq_count columns")
		exit()
	repertoire = repertoire.drop_duplicates(subset=["clone_id"])
	for x in repertoire["clone_id"]:
		clones[x]=repertoire[repertoire["clone_id"]==x]["clone_seq_count"].values[0]
		# exit()
	# clones = repertoire.set_index('clone_id')
	# clones = repertoire.to_dict('records')
	# print(clones)
	print("---------------------------------------------")
	print("DIVERSITY REPORT")
	print("\n")
	print("Simpson diversity:"+str(1-simpson_di(clones)))
	print("\n")
	print("Shannon diverity: "+str(shannon_di(clones)))
	print("\n")
	print("Shannon eveness: "+str(shannon_di(clones)/len(clones)))
	print("---------------------------------------------")
	print("\n")
	exit()









def most_common(lst):
    return max(set(lst), key=lst.count)

filename = ""

if (any(".tsv" in i for i in sys.argv) == True) and (any("--input" in i for i in sys.argv) == False):
	filename =	sys.argv[1]
elif any("--input" in i for i in sys.argv) == False:
	print("Please provide a file path")
	exit()
else:
	filename = sys.argv[sys.argv.index("--input")+1]

print(filename)

method = "AHAM"
thr = 0.09
sequence_column = "cdr3"
vcolumn = "v_call"
jcolumn = "j_call"
seqID = "sequence_id"
memory_usage = "default"
path = directory_path(filename)
filename_temp = filename.split(".")
out_filename = filename_temp[0]+"_YClon_clonotyped."+filename_temp[1]
ksize = 3
short_output = False
every_in_the_folder = False


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
	elif sys.argv[x].find("--kmer_length") != -1:
		ksize = int(sys.argv[x+1])
	elif sys.argv[x].find("--low_memory") != -1:
		memory_usage = "low"
	elif sys.argv[x].find("--dir_out") != -1:
		path = sys.argv[x+1]
		out_filename = os.path.join(path,os.path.basename(out_filename))
	elif sys.argv[x].find("--out") != -1:
		out_filename = os.path.join(path,sys.argv[x+1]+"_YClon_clonotyped."+filename_temp[1])
		#path_and_name = os.path.join(path,sys.argv[x+1])
	elif sys.argv[x].find("--short_output") != -1:
		short_output = True
		out_small_name = out_filename.replace("_YClon_clonotyped.","_YClon_clonotyped_only_essential_columns.")
		# os.path.join(path,sys.argv[x+1]+"_YClon_clonotyped_only_essential_columns."+filename_temp[1])
		# if sys.argv[x+1] != "--":
		# 	out_small_name = sys.argv[x+1]
	elif sys.argv[x].find("--folder") != -1:
		every_in_the_folder = True
		filename = sys.argv[x+1]


def build_kmers_tf_idf(sequence, ksize=ksize): 
		ngrams = zip(*[sequence[i:] for i in range(ksize)])
		return [''.join(ngram) for ngram in ngrams]


def clonotyping(filename, thr, sequence_column, vcolumn, jcolumn, seqID, separator, memory_usage, short_output):
	f = open(filename, 'r')
	x = f.readline().strip()
	head = x.split(separator)

	# for i in range(len(head)):
	# 	head[i] = head[i]

	number_of_columns = len(head)

	try:
		seq_id_indx = head.index(seqID)
	except:
		print("\nWARNING\nThere is no column named "+seqID+"\n")
		exit()

	try:
		junc_indx = head.index(sequence_column)
	except:
		print("\nWARNING\nThere is no column named "+sequence_column+"\n")
		exit()

	try:
		vGene_indx = head.index(vcolumn)
	except:
		print("\nWARNING\nThere is no column named "+vcolumn+"\n")
		exit()

	try:
		jGene_indx = head.index(jcolumn)
	except:
		print("\nWARNING\nThere is no column named "+jcolumn+"\n")
		exit()

	colunas = [head[seq_id_indx],head[junc_indx]]


	i=0
	fail = 0
	clonotypes = {}
	file_size = 0
	print("Processing file")
	for x in f:
		file_size +=1
		# x = x.strip()
		data = list(x.split(separator))
		# print(len(data))
		# print(data)
		# exit()
		if len(data)!= number_of_columns:
			continue
		try:
			if len(data[junc_indx]) < 4:
				continue
			cdr3len = str(len(data[junc_indx]))
			vGene = data[vGene_indx].split(',')[0].split('*') #include all v gene alleles
			jGene = data[jGene_indx].split(',')[0].split('*') #include all v gene alleles
			# print(vGene[0]+" "+jGene[0]+"\n")
			if (jGene != "") and (vGene != "") and (cdr3len != 0):
				key = vGene[0]+","+jGene[0]+","+cdr3len
			clonotypes.setdefault(key, [])
			proclonotype = [data[seq_id_indx],data[junc_indx]]
			clonotypes[key].append(proclonotype)
		except:
			fail +=1
			



		
	f.close()
	temp_filename = path+"YClon_temp.txt"
	size = len(clonotypes)
	temp = open(temp_filename, 'w')
	count = 0
	total_clust = 0
	unico_pq_VJLen = 0
	maior = 0
	a = 0
	if memory_usage == "low":
		with alive_bar(len(clonotypes), title="Clonotyping") as bar: 
			for key in clonotypes: #each key is the combination of V gene, J gene and the length of cdr3, the values are the sequence ID and cdr3 sequence
				bar()
				# if bar.current() == 3336:
				#	 print(clonotypes[key])
				if len(clonotypes[key]) > 1:
					pre_clone = pd.DataFrame(clonotypes[key])
					pre_clone.columns = colunas
					unique_cdr3 = pre_clone.drop_duplicates(sequence_column)
				else:
					unique_cdr3 = clonotypes[key]
				if len(unique_cdr3) > 1:
					# pre_clone = pd.DataFrame(clonotypes[key])
					# pre_clone.columns = colunas
					# unique_cdr3 = pre_clone.drop_duplicates(sequence_column)


					# junc_seq = pre_clone[sequence_column]
					seq_id = pre_clone[seqID]

					junc_seq = unique_cdr3[sequence_column]
					seq_unique = unique_cdr3[seqID]


					# vectorizer = TfidfVectorizer(min_df=1, analyzer=build_kmers_tf_idf)
					vectorizer = CountVectorizer(min_df=1, analyzer=build_kmers_tf_idf)
					tf_idf_matrix = vectorizer.fit_transform(junc_seq)	
					dist = 1 - cosine_similarity(tf_idf_matrix)
					if method ==	"AHAM":
						clusterer = AgglomerativeClustering(distance_threshold = thr, n_clusters= None, metric="precomputed",linkage='average')
					elif method == "HDBSCAN":
						clusterer = hdbscan.HDBSCAN(min_cluster_size=3, cluster_selection_epsilon= thr, metric='precomputed')
					cluster_pre_clone = clusterer.fit(dist)
					clone = cluster_pre_clone.labels_
					# if cluster_pre_clone.labels_.max() != -1:
					maior = count + cluster_pre_clone.labels_.max() + 1
					# print("clone_size: "+str(len(clone)))
					# print("junc_seq: "+str(len(junc_seq)))
					# print(junc_seq.index)
					junc_seq = junc_seq.reset_index(drop=True) 
					# print(junc_seq[0])
					for i in range(0, len(clone)):
						if clone[i] != -1:
							clone_id = count + int(clone[i]) +1
						else:
							clone_id = maior + 1
							maior += 1
							total_clust += 1
						
						# print(junc_seq[i])

						all_again = pre_clone.loc[pre_clone[sequence_column] == junc_seq[i]][seqID]
						# print(all_again)
						for k in all_again:
							temp.write(k+","+str(clone_id)+"\n")	
					count = maior
					total_clust += cluster_pre_clone.labels_.max() + 1	
				else:
					unico_pq_VJLen +=1
					count += 1
					total_clust += 1
					pre_clone = pd.DataFrame(clonotypes[key])
					pre_clone.columns = colunas
					seq_id = pre_clone[seqID]
					for i in range(0, len(clonotypes[key])):
						temp.write(str(seq_id[i])+","+str(count)+"\n")
	else:
		with alive_bar(len(clonotypes), title="Clonotyping") as bar: 
			for key in clonotypes: #each key is the combination of V gene, J gene and the length of cdr3, the values are the sequence ID and cdr3 sequence
				bar()
				if len(clonotypes[key]) > 1:
					teste = pd.DataFrame(clonotypes[key])
					teste.columns = colunas

					junc_seq = teste[sequence_column]
					seq_id = teste[seqID]


					# vectorizer = TfidfVectorizer(min_df=1, analyzer=build_kmers_tf_idf)
					vectorizer = CountVectorizer(min_df=1, analyzer=build_kmers_tf_idf)
					tf_idf_matrix = vectorizer.fit_transform(junc_seq)	
					dist = 1 - cosine_similarity(tf_idf_matrix)
					if method ==	"AHAM":
						clusterer = AgglomerativeClustering(distance_threshold = thr, n_clusters= None, metric="precomputed",linkage='average')
					elif method == "HDBSCAN":
						clusterer = hdbscan.HDBSCAN(min_cluster_size=3, cluster_selection_epsilon= thr, metric='precomputed')
					cluster_teste = clusterer.fit(dist)
					clone = cluster_teste.labels_

					# if cluster_teste.labels_.max() != -1:
					maior = count + cluster_teste.labels_.max() + 1

					for i in range(0, len(clonotypes[key])):
						# print("vai escrever")
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

	pre_clone = []


			
	temp.close()
	in_temp = open(temp_filename, 'r')
	in_airr = open(filename, 'r')
	filename_temp = filename.split(".")
	out = open(out_filename, 'w+')
	clonotipo = {}
	maior ={}
	maximo = 0

	
	for x in in_temp:
		data = x.strip().split(",")
		clonotipo[','.join(data[0:-1])] = data[-1].strip()
		if data[-1] not in maior:
			maior[data[-1]] = []
			maior[data[-1]].append(','.join(data[0:-1]))
		else:
			maior[data[-1]].append(','.join(data[0:-1]))

	most_common_cdr3 = {}
	most_common_seq_id = {}
	seq_list = []
 
	with alive_bar(file_size+1, title="Writing output file") as bar:
		for x in in_airr:
			bar()
			if x.find(seqID) == -1:
				data = x.split(separator)

				if data[seq_id_indx] in clonotipo:
					for i in range(0, len(data)):
						out.write(data[i].strip()+separator)
					out.write(clonotipo[data[seq_id_indx]]+"\n")
					if clonotipo[data[seq_id_indx]] not in most_common_cdr3:
						most_common_cdr3[clonotipo[data[seq_id_indx]]] = []
						most_common_seq_id[clonotipo[data[seq_id_indx]]] = []
						most_common_cdr3[clonotipo[data[seq_id_indx]]].append(data[junc_indx].strip())
						most_common_seq_id[clonotipo[data[seq_id_indx]]].append(data[seq_id_indx])
					else:
						most_common_cdr3[clonotipo[data[seq_id_indx]]].append(data[junc_indx].strip())
						most_common_seq_id[clonotipo[data[seq_id_indx]]].append(data[seq_id_indx])
			else:
				data = x.split(separator)
				seq_id_indx = data.index(seqID)
				out.write(x.strip()+separator+"clone_id\n")


	if short_output == True:
		out_small = open(out_small_name, 'w+')

		in_airr = open(filename, 'r')


		#with alive_bar(file_size+1, title="Writing short output file") as bar:
		for x in in_airr:
			bar()
			if x.find(seqID) == -1:
				data = x.strip().split(separator)
				if data[seq_id_indx] in clonotipo:
					out_small.write(data[seq_id_indx]+separator+data[vGene_indx]+separator+data[jGene_indx]+separator+data[junc_indx]+separator)
					out_small.write(clonotipo[data[seq_id_indx]])
			else:
				data = x.split(separator)
				seq_id_indx = data.index(seqID)
				out_small.write(seqID+separator+vcolumn+separator+jcolumn+separator+sequence_column+separator+"clone_id\n")
				#out_small.write(x.strip()+separator+"clone_idn")


	os.remove(temp_filename)
	out.close()
	temp_filename = path+"YClon_temp.txt"
	temp = open(temp_filename, 'w')
	out = open(out_filename, 'r')


	out_report = open(out_filename.split("_")[0]+"_YClon_clone_report.tsv","w")
	out_report.write("sequence_id\tseq_count\tmost_common_cdr3\tclone_id\n")
	# print(most_common_seq_id)
	# print(maior)

	#write the seq_count of each clone in row

	for x in out:
		if x.find(seqID) != -1:
			temp.write(x.strip()+separator+"clone_seq_count\n")
			i = x.strip().split(separator).index("clone_id")
		else:
			temp.write(x.strip()+separator+str(len(maior[x.strip().split(separator)[i]]))+"\n")
	temp.close()
	out.close()
	os.rename(temp_filename,out_filename)
	for i in most_common_cdr3:
	  # print(i)
	  cdr3 = most_common(most_common_cdr3[i])
	  # print(most_common_seq_id[i])
	  out_report.write(most_common_seq_id[i][most_common_cdr3[i].index(cdr3)]+"\t"+str(len(maior[i]))+"\t"+cdr3+"\t"+i+"\n")

	current_time = time.time()
	elapsed_time = current_time - start_time

	if short_output == True:
		os.remove(out_filename)
	print("The work was completed in: " + "%.3f" % int(elapsed_time))


path = directory_path(filename)

if every_in_the_folder == False:
	clonotyping(filename, thr, sequence_column, vcolumn, jcolumn, seqID, separator, memory_usage, short_output)
else:
	files = os.listdir(path)
	for x in files:
		clonotyping(path+"/"+x, thr, sequence_column, vcolumn, jcolumn, seqID, separator, memory_usage, short_output)



