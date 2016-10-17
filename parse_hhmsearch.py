
# -*- coding: utf-8 -*-

##########################################################################################
##########################################################################################
##
##								Library
##
##########################################################################################
##########################################################################################

import sys, os
import glob
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
#get_ipython().magic('matplotlib notebook')
from matplotlib.backends.backend_pdf import PdfPages


##########################################################################################
##########################################################################################
##
##								Function
##
##########################################################################################
##########################################################################################


##########################################################################################
##########################################################################################

def read_hhr(hhr_file, dataframe):

	"""
	Function that read a .hhr and return write it in the pandas.Dataframe  
	
	:param hhr_file: the .hhr file that hhsearch write  
	:type: str  
	:param dataframe: the dataframe that we want to add the identities  
	:type: pandas.Dataframe  
	"""

	name=''
	first_line=True
	with open(hhr_file, 'r') as hhr :
		for line in hhr :
			if first_line:
				name2=dict_annotation[line.split()[1].split("_NB-0")[0].split(".")[0]]
				first_line=False
			if '>' in line :
				name= dict_annotation[line[1:].split(".")[0].split("_NB-0")[0]]
			elif name!='' :
				Probab, E_value, Score, Aligned_cols, Identities, Similarity, Sum_probs= line.split()
				if float(E_value[8:]) <= 1e-3 :
					if dataframe.loc[name, name2] == 0 :
						if name == name2 :
							dataframe.loc[name, name2] = 0
						else :
							dataframe.loc[name, name2] = int(Identities[11:-1])
				name=''
	return dataframe
				




##########################################################################################
##########################################################################################

def read_hhr_pvalue(hhr_file, dataframe):

	"""
	Function that read a .hhr and return write it in the pandas.Dataframe and select the p-value instead of e-value  
	
	:param hhr_file: the .hhr file that hhsearch write  
	:type: str  
	:param dataframe: the dataframe that we want to add the identities  
	:type: pandas.Dataframe
	"""

	tabular=False
	name=''
	first_line=True
	dict_p_value = {}
	with open(hhr_file, 'r') as hhr :
		for line in hhr :
			if first_line:
				name2=dict_annotation[line.split()[1].split("_NB-0")[0].split(".")[0]]
				first_line=False
			if "No Hit" in line :
				tabular = True
			elif "\n" == line :
				tabular = False
			elif tabular :
				name1 = dict_annotation[line.split()[1].split(".")[0].split("_NB-0")[0]]
				if (name1,name2) not in dict_p_value :
					dict_p_value[(name1,name2)] = float(line.split()[4])
			elif '>' in line :
				name= dict_annotation[line[1:].split(".")[0].split("_NB-0")[0]]
			elif name!='' :
				Probab, E_value, Score, Aligned_cols, Identities, Similarity, Sum_probs= line.split()
				if dict_p_value[(name,name2)] <= 1e-3 :
					if dataframe.loc[name, name2] == 0 :
						if name == name2 :
							dataframe.loc[name, name2] = 0
						else :
							dataframe.loc[name, name2] = int(Identities[11:-1])
				name=''

	return dataframe
				
##########################################################################################
##########################################################################################

def create_adjency_matrix(dict_annotation, file_hhr):

	"""
	Function that create the adjency matrix (create a identity dataframe)

	:param dict_annotation: dictionnary with the iformation from the annotation table 
	"""
	length = len(info_for_dict[:,1])
	identity_df = pd.DataFrame(data=np.zeros((length, length), dtype=int), index=info_for_dict[:,1] ,columns=info_for_dict[:,1])
	for fileHHR in files_hhr :
		identity_df = read_hhr(fileHHR, identity_df)
	
	plt.figure()
	plt.imshow(identity_df, interpolation='nearest')
	plt.title("Identity matrix with all profiles")
	pdf.savefig()
	plt.close()
	
	return identity_df

##########################################################################################
##########################################################################################
##
##								Main
##
##########################################################################################
##########################################################################################
# Analyse avec juste T2SS, T4P, Archaellum, Tad


files_hhr = glob.glob("/Users/rdenise/Documents/Analysis_similarity/result_hhsearch/*")



info_for_dict = np.loadtxt("/Users/rdenise/Documents/Analysis_similarity/annotation_list.txt", dtype="string", delimiter=";")
dict_annotation = {line[0]:line[1] for line in info_for_dict}


with PdfPages("All_graph.pdf") as pdf :




graph = nx.from_numpy_matrix(identity_df.values)
graph = nx.relabel_nodes(graph, dict(enumerate(identity_df.columns)))
outdeg = graph.degree()
to_remove = [n for n in outdeg if outdeg[n] == 0]
print(to_remove)
graph.remove_nodes_from(to_remove)
nx.write_gml(graph, '/Users/rdenise/Documents/Analysis_similarity/network.gml')



nx.draw_networkx(graph)



length = len(info_for_dict[:,1])
identity_df_pvalue = pd.DataFrame(data=np.zeros((length, length), dtype=int), index=info_for_dict[:,1] ,columns=info_for_dict[:,1])
for fileHHR in files_hhr :
	identity_df_pvalue = read_hhr_pvalue(fileHHR, identity_df_pvalue)
plt.imshow(identity_df_pvalue, interpolation='nearest')



graph = nx.from_numpy_matrix(identity_df_pvalue.values)
graph = nx.relabel_nodes(graph, dict(enumerate(identity_df_pvalue.columns)))
outdeg = graph.degree()
to_remove = [n for n in outdeg if outdeg[n] == 0]
print(to_remove)
graph.remove_nodes_from(to_remove)
nx.write_gml(graph, '/Users/rdenise/Documents/Analysis_similarity/network_pvalue.gml')



nx.draw_networkx(graph)



files_hhr_less = glob.glob("/Users/rdenise/Documents/Analysis_similarity/result_hhsearch_less_T4SS/*")



info_for_dict2 = np.loadtxt("/Users/rdenise/Documents/Analysis_similarity/annotation_list2.txt", dtype="string", delimiter=";")
dict_annotation2 = {line[0]:line[1] for line in info_for_dict2}



length = len(info_for_dict2[:,1])
identity_df2 = pd.DataFrame(data=np.zeros((length, length), dtype=int), index=info_for_dict2[:,1] ,columns=info_for_dict2[:,1])
for fileHHR in files_hhr_less :
	identity_df2 = read_hhr(fileHHR, identity_df2)
plt.imshow(identity_df2, interpolation='nearest')



graph = nx.from_numpy_matrix(identity_df2.values)
graph = nx.relabel_nodes(graph, dict(enumerate(identity_df2.columns)))
outdeg = graph.degree()
to_remove = [n for n in outdeg if outdeg[n] == 0]
print(to_remove)
graph.remove_nodes_from(to_remove)
nx.write_gml(graph, '/Users/rdenise/Documents/Analysis_similarity/network2.gml')



nx.draw_networkx(graph)



length = len(info_for_dict2[:,1])
identity_df_pvalue2 = pd.DataFrame(data=np.zeros((length, length), dtype=int), index=info_for_dict2[:,1] ,columns=info_for_dict2[:,1])
for fileHHR in files_hhr_less :
	identity_df_pvalue2 = read_hhr_pvalue(fileHHR, identity_df_pvalue2)
plt.imshow(identity_df_pvalue2, interpolation='nearest')



graph = nx.from_numpy_matrix(identity_df_pvalue2.values)
graph = nx.relabel_nodes(graph, dict(enumerate(identity_df_pvalue2.columns)))
outdeg = graph.degree()
to_remove = [n for n in outdeg if outdeg[n] == 0]
print(to_remove)
graph.remove_nodes_from(to_remove)
nx.write_gml(graph, '/Users/rdenise/Documents/Analysis_similarity/network_pvalue2.gml')



nx.draw_networkx(graph)





