
# -*- coding: utf-8 -*-

##########################################################################################
##########################################################################################
##
##								Library
##
##########################################################################################
##########################################################################################

import argparse
from textwrap import dedent
import sys, os
import time
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


def create_folder(mypath):

	"""
	Created the folder that I need to store my result if it doesn't exist
	:param mypath: path where I want the folder (write at the end of the path)
	:type: string
	:return: Nothing
	"""

	try:
		os.makedirs(mypath)
	except OSError:
		pass

	return

##########################################################################################
##########################################################################################

def create_dict_length_profile(hhr_files, dict_annotation):

	"""
	Function that read a .hhr and return write it in the pandas.Dataframe and select the p-value instead of e-value

	:param hhr_file: the .hhr files that hhsearch write
	:type: str
	:param dict_annotation: dict that contain the name of the profiles before the first "." and the name transformed.
	:type: dict
	:return: the dictionary of the length of all the profiles
	:rtype: dict
	"""

	dict_length = {}

	for hhr_file in hhr_files :
		with open(hhr_file, 'r') as hhr :
			for line in hhr :
				if "Query" in line:
					name = dict_annotation[line.split()[1].split("_NB-0")[0].split(".")[0][:30]]
				elif "Match_columns" in line :
					dict_length[name] = float(line.split()[1])
				else :
					break

	return dict_length

##########################################################################################
##########################################################################################

def read_hhr_evalue(hhr_file, dataframe, dict_annotation, dict_length, MIN_COVERAGE):

	"""
	Function that read a .hhr and return write it in the pandas.Dataframe

	:param hhr_file: the .hhr file that hhsearch write
	:type: str
	:param dataframe: the dataframe that we want to add the identities
	:type: pandas.Dataframe
	:param dict_annotation: dict that contain the name of the profiles before the first "." and the name transformed.
	:type: dict
	:param dict_length: the dict of the length of each profiles
	:type: dict
	:param MIN_COVERAGE: The minimum coverage wanted by the user
	:type: float
	:return: the dataframe with the identities added
	:rtype: pandas.Dataframe
	"""

	name=''
	first_line=True
	with open(hhr_file, 'r') as hhr :
		for line in hhr :
			if first_line:
				name2 = dict_annotation[line.split()[1].split("_NB-0")[0].split(".")[0][:30]]
				first_line=False
			if '>' in line :
				name = dict_annotation[line.rstrip()[1:].split(".")[0].split("_NB-0")[0][:30]]
			elif name!='' :
				Probab, E_value, Score, Aligned_cols, Identities, Similarity, Sum_probs= line.split()
				minimum_length = min(dict_length[name],dict_length[name2])
				coverage = float(Aligned_cols.split("=")[1])/minimum_length
				if float(E_value[8:]) <= 1e-3 and coverage > MIN_COVERAGE:
					if dataframe.loc[name, name2] == 0 :
						if name == name2 :
							dataframe.loc[name, name2] = 0
						else :
							dataframe.loc[name, name2] = int(Identities[11:-1])
				name=''
	return dataframe

##########################################################################################
##########################################################################################

def read_hhr_pvalue(hhr_file, dataframe, dict_annotation, dict_length, MIN_COVERAGE):

	"""
	Function that read a .hhr and return write it in the pandas.Dataframe and select the p-value instead of e-value

	:param hhr_file: the .hhr file that hhsearch write
	:type: str
	:param dataframe: the dataframe that we want to add the identities
	:type: pandas.Dataframe
	:param dict_annotation: dict that contain the name of the profiles before the first "." and the name transformed.
	:type: dict
	:param dict_length: the dict of the length of each profiles
	:type: dict
	:param MIN_COVERAGE: The minimum coverage wanted by the user
	:type: float
	:return: the dataframe with the identities added
	:rtype: pandas.Dataframe
	"""

	tabular=False
	name=''
	first_line=True
	dict_p_value = {}
	with open(hhr_file, 'r') as hhr :
		for line in hhr :
			if first_line:
				name2 = dict_annotation[line.split()[1].split("_NB-0")[0].split(".")[0][:30]]
				first_line=False
			if "No Hit" in line :
				tabular = True
			elif "\n" == line :
				tabular = False
			elif tabular :
				# NOTE Si nom trop long probleme
				name1 = dict_annotation[line.split()[1].split(".")[0].split("_NB-0")[0][:30]]
				if (name1,name2) not in dict_p_value :
					dict_p_value[(name1,name2)] = float(line.split()[4])
			elif '>' in line :
				name = dict_annotation[line.rstrip()[1:].split(".")[0].split("_NB-0")[0][:30]]
			elif name!='' :
				Probab, E_value, Score, Aligned_cols, Identities, Similarity, Sum_probs = line.split()
				minimum_length = min(dict_length[name],dict_length[name2])
				coverage = float(Aligned_cols.split("=")[1])/minimum_length
				if dict_p_value[(name,name2)] <= 1e-3 and coverage > MIN_COVERAGE:
					if dataframe.loc[name, name2] == 0 :
						if name == name2 :
							dataframe.loc[name, name2] = 0
						else :
							dataframe.loc[name, name2] = int(Identities[11:-1])
				name=''

	return dataframe

##########################################################################################
##########################################################################################

def create_adjency_matrix(dict_annotation, fileshhr, value, dict_length, coverage_min):

	"""
	Function that create the adjency matrix (create a identity dataframe)

	:param dict_annotation: dict that contain the name of the profiles before the first "." and the name transformed.
	:type: dict
	:param fileshhr: the list of all the .hhr files.
	:type: list of str
	:param value: The name of the value (Pvalue or Evalue you want to use)
	:type: str
	:param dict_length: the dict of the length of each profiles
	:type: dict
	:param coverage_min: The minimum coverage wanted by the user
	:type: float
	:return: a dataframe with the information of the percentage of identity for all the profiles
	:rtype: pandas.Dataframe
	"""

	length = len(info_for_dict[:,1])
	identity_df = pd.DataFrame(data=np.zeros((length, length), dtype=int), index=info_for_dict[:,1] ,columns=info_for_dict[:,1])

	print("#################")
	print("## Read HHR files")
	print("#################")

	progression=0

	if value == "Evalue" :
		for fileHHR in fileshhr :
			progression+=1
			sys.stdout.write("{:.2f}% : {}/{} sequences\r".format(progression/float(length)*100, progression, length))
			sys.stdout.flush()

			identity_df = read_hhr_evalue(fileHHR, identity_df, dict_annotation, dict_length, coverage_min)
	elif value == "Pvalue" :
		for fileHHR in fileshhr :

			progression+=1
			sys.stdout.write("{:.2f}% : {}/{} sequences\r".format(progression/float(length)*100, progression, length))
			sys.stdout.flush()

			identity_df = read_hhr_pvalue(fileHHR, identity_df, dict_annotation, dict_length, coverage_min)

	print()
	print("Done !")

	plt.figure()
	plt.imshow(identity_df, interpolation='nearest')
	plt.title("Identity matrix with all profiles using {}".format(value))
	pdf.savefig()
	plt.close()

	return identity_df

##########################################################################################
##########################################################################################

def make_graph(identity_df, PATH_TO_RESULTS, value, coverage):

	"""
	Create the graph with the indentity dataframe that containt only the node with at least a number of degree of one.

	:param identity_df: the dataframe with the identities added
	:type: pandas.Dataframe
	:param: path to the result folder
	:type: str
	:param value: The name of the value (Pvalue or Evalue you want to use)
	:type: str
	:param coverage: The minimum coverage wanted by the user
	:type: float
	:return: Nothing
	"""

	graph = nx.from_numpy_matrix(identity_df.values)
	graph = nx.relabel_nodes(graph, dict(enumerate(identity_df.columns)))
	outdeg = graph.degree()
	to_remove = [n for n in outdeg if outdeg[n] == 0]

	print("#######################")
	print("## All the node removed")
	print("#######################")

	print(to_remove)

	graph.remove_nodes_from(to_remove)
	nx.write_gml(graph, os.path.join(PATH_TO_RESULTS,'network_{}_{}.gml'.format(value, coverage)))

	plt.figure()
	nx.draw_networkx(graph)
	plt.title("Graph of the identity between the profiles using {} and a coverage of {}".format(value, coverage))
	pdf.savefig()
	plt.close()

	return

##########################################################################################
##########################################################################################
##
##								Main
##
##########################################################################################
##########################################################################################

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
     description=dedent("""


-----------------------------------------

  _   _ _   _ ____     ____                 _
 | | | | | | |  _ \   / ___|_ __ __ _ _ __ | |__
 | |_| | |_| | |_) | | |  _| '__/ _` | '_ \| '_ \\
 |  _  |  _  |  _ <  | |_| | | | (_| | |_) | | | |
 |_| |_|_| |_|_| \_\  \____|_|  \__,_| .__/|_| |_|
                                     |_|


-----------------------------------------
""") )

general_option = parser.add_argument_group(title = "General input dataset options")
general_option.add_argument("-hhr",'--hhrfolder',
 							required=True,
							metavar="<path>",
							dest="HHRFolder",
							help="Path to the HHR result folder")
general_option.add_argument("-a",'--annotationFile',
 							required=True,
							metavar="<file>",
							dest="annotFile",
							help="File with the information about the annotations")
general_option.add_argument("-o",'--output',
 							default=None,
							dest="output",
							metavar='<path>',
							help="Using <path> for output files (default: HHRFolder directory)")
general_option.add_argument("-c",'--coverage_min',
 							default=None,
							dest="coverage_min",
							type=float,
							metavar='<COVERAGE_MIN>',
							help="The minimum coverage for the alignmenent in the results (between 0 and 1)")


args = parser.parse_args()


files_hhr = glob.glob(os.path.join(args.HHRFolder,"*"))

if args.output :
	PATH_TO_RESULTS=os.path.join(args.output, "results")
else :
	PATH_TO_RESULTS=(args.HHRFolder, "results")

create_folder(PATH_TO_RESULTS)

info_for_dict = np.genfromtxt(args.annotFile, dtype=str, delimiter=";")
DICT_ANNOTATION = {line[0][:30]:line[1] for line in info_for_dict}

DICT_LENGTH = create_dict_length_profile(files_hhr, DICT_ANNOTATION)

if args.coverage_min :
	COVERAGE = args.coverage_min
else :
	COVERAGE = 0.8

with PdfPages(os.path.join(PATH_TO_RESULTS,"All_graph.pdf")) as pdf :
		# NOTE using Evalue
		identity_DF=create_adjency_matrix(DICT_ANNOTATION, files_hhr, "Evalue", DICT_LENGTH, COVERAGE)
		make_graph(identity_DF, PATH_TO_RESULTS, "Evalue", COVERAGE)

		# NOTE using Pvalue
		identity_DF=create_adjency_matrix(DICT_ANNOTATION, files_hhr, "Pvalue", DICT_LENGTH, COVERAGE)
		make_graph(identity_DF, PATH_TO_RESULTS, "Pvalue", COVERAGE)
