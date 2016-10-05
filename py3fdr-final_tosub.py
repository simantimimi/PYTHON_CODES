'''
TITLE:  GENE ENRICHMENT WITH HYPERGEOMETRIC CALCULATION and MULTIPLE TESTING AND P-VALUE CORRECTION

Team:Computational biology
Contributors: Simanti Bhattacharya, Amit Das and Manimaran Nellai Chinnaih; 
Date Created: 29-09-2016;
Code Version:1
Python Version: Python 3.4

'''
#======================================================================
'''
##PREREQUISITES:
DATABASE FORMAT: GOID|TERM|ASPECT\tGENE1\tGENE2\tGENE3,...\tGENE'N' (where \t stands for tab separated)
INPUT GENE SET FORMAT: GENE1,GENE2,GENE3,...,GENE'N' and NAME SHOULD BE DISEASENAME.txt OR DRUGNAME.txt LIKE THIS

**CHANGE in PATH\TO\Anaconda\Lib\site-packages\statsmodels\stats\multitest.py line number 260 -->from return reject_, pvals_corrected_, alphacSidak, alphacBonf --> to--> return pvals_corrected_
'''
#======================================================================
import operator
#import statsmodels
#from statsmodels import *
import statsmodels.sandbox.stats.multicomp
from scipy.stats import hypergeom
import scipy.stats as stats
import csv
import numpy as np
import sys
import time
import os
start_time=time.time()

#======================================================================

#### %*%*%*% MODULE5: USEFULL WHEN GENEID AS INPUT %*%*%*%  #######

#======================================================================
'''def id2nam(x):
	f2 = open("geneids2_names.txt",'r')
	for i in f2:
		i = i.strip()
		gdet = i.split("\t")
		if x == gdet[0]:
			return(gdet[1])
	f2.close()'''

#======================================================================

#### %*%*%*% MODULE3: HYPERGEOMETRIC CALCULATION %*%*%*%  #######

#======================================================================


def hypergeometry(d, i, out):
	inputgenes = inputgene(i) 						### CALLING MODULE1; NOTE: i==sys.argv[2]==inputgene sets
	f = open(d,'r') 								### NOTE: d==data==sys.argv[1]==database
	datas = f.readlines()
	inputlen = len(inputgenes) 						### THIS IS 'S' OR THE INPUT NUMBER OR SAMPLE 
	#print(inputlen)
	#print(inputgenes)
	useroutputfile = out 							### out==sys.argv[3]==user defined output file in .csv
	global M
	M = universegene(d)								### CALLING MODULE2; NOTE: d==data==sys.argv[1]==database; GETS 'P'/ THE COUNT OF TOTAL POPULATION/UNIVERSE SET
	with open(useroutputfile,'w',newline='') as outputfile:
		global a
		global head
		a = csv.writer(outputfile,delimiter='$')
		head = ("FIELD","GOID|GO_TERM|ASPECT","Universe set","No. of genes","No. of input genes","No. of overlapping genes","Gene Names","p-value","Adj_BH") ### Add "ADJ_Bon" heading for bonferroni and row[0] at 138th position.
		a.writerow(head)
		b=[]
		rv2=[]
		for i in range(0,len(datas)):
			match = []
			paths = datas[i].strip().split("\t")		### paths CONTAINS DATABASE INFORMATION, SO paths[0] FIRST COL OF DATABASE
			#print (paths)##WORKING
			pathgenes = paths[1:len(paths)]				
			pathlen = len(pathgenes)					### THIS IS 'p' OR SUCCESS FROM POPULATION
			if len(pathgenes) > 1:
				pathgenes[-1] = pathgenes[-1].strip("\n")     
			else:  
				pathgenes[0] = pathgenes[0].strip("\n")
			#print(pathlen)
			#print(pathgenes)
			match = list(set(pathgenes).intersection(inputgenes))       
			matchlen  = len(match)  					### THIS IS 's' OR SUCCESS FROM SAMPLE/INPUT
			#print(matchlen)
			#print(match)
			#print(M)
			#print(inputlen)
			hpd = stats.hypergeom(M,pathlen,inputlen)   ### USING P,p,S
			pval = np.float64()
			pval = hpd.pmf(matchlen)						### HYPERGEOMETRIC P-VALUE WITH hpd AND 's'
			#print(pval) ##WORKING
			name = paths[0]
			ff1=sys.argv[1].split(".")[0] # name1 		##CHANGE ACCORDING TO YOUR INPUT
			#name = paths[0]+'|'+paths[1]
			t2 = []
			tname = []
			for l in match:
				t2.append(l) 
				#tname.append(id2nam(l))
			str1 = " ".join(t2)
			#str2 = " ".join(tname)  
			if matchlen > 0:  								### GETTING THOSE GOID ROWS FOR WHICH THERE IS AN OVERLAP BETWEEN INPUT AND ANNOTATED GENE SETS
				b.append(pval)
				rv1 = (ff1+"$"+name+"$"+str(M)+"$"+str(pathlen)+"$"+str(inputlen)+"$"+str(matchlen)+"$"+str1+"$"+str(pval))
				rv2.append(rv1)
		#print (b)#LIST OF PVALS
		#padj=statsmodels.sandbox.stats.multicomp.multipletests(b, alpha=0.05, method='hs', is_sorted=False, returnsorted=False)
		sorted(b)
		padj_fdr=statsmodels.sandbox.stats.multicomp.multipletests(b, alpha=0.01, method='fdr_bh')    ### ALTERNATIVE METHODS: 'fdr_bh', 'fdr_i', 'fdr_p', 'fdri', 'fdrp' for Benjamini-Hochberg
		padj_b=statsmodels.sandbox.stats.multicomp.multipletests(b, alpha=0.01, method='b')    ### ALTERNATIVE METHODS: 'b', 'bonf', 'bonferroni' for bonferroni
		for v1,v2,v3 in zip (rv2,padj_fdr,padj_b): 
			v4 = (v3,v2,v1)
			a.writerow (v4)
		#print (outputfile)
		outputfile.close()
		mysort(useroutputfile) 							### THIS FUNCTION IS NOT WORKING DO MANUAL SORTING

#======================================================================

#### %*%*%*% MODULE4: SORTING FINAL RESULT %*%*%*%  #######

#======================================================================

def mysort(useroutputfile):
	with open(useroutputfile,'r') as data:
		data1 = csv.reader(data,delimiter='$')
		headers = next(data1)
		global sortedlist
		sortedlist = sorted(data1,key=lambda x:float(x[1]))    # 0 specifies according to first column we want to sort
		#print (sortedlist)
		data.close()
	with open(useroutputfile,'w',newline='') as wtsor:
		filewriter=csv.writer(wtsor,delimiter='$') 
		filewriter.writerow(head)
		for row in sortedlist:
			new =(row[2].split("$")[0],row[2].split("$")[1],row[2].split("$")[2],row[2].split("$")[3],row[2].split("$")[4],row[2].split("$")[5],row[2].split("$")[6],row[2].split("$")[7],row[1]) #### you can add row[0] i.e. Bonferroni p value and change head accordingly
			#print (new)
			filewriter.writerow(new)

#======================================================================
			
# %*%*%*% MODULE2: PROCESSING DATABASE SET: GIVES THE COUNT OF TOTAL POPULATION/UNIVERSE SET/ P %*%*%*%  #

#======================================================================

def universegene(myd):
	all_gene = []
	f1 = open(myd,'r')
	c1 = f1.readlines()
	#print(c1) 
	for i in range(0,len(c1)):
		c1[i] = c1[i].strip("\n")
		ing = c1[i].split("\t")
		tst = set(all_gene + ing[1:]) ## TAKING UP THE GENE COLUMN THAT STARTS AT 2ND COLUMN TO LAST
		all_gene = list(tst)
	c2 = len(all_gene) ## P OR THE TOTAL POPULATION
	#print(c2)
	return c2 
	
#======================================================================

#### %*%*%*% MODULE1: PROCESSING INPUT SET %*%*%*%  #######

#======================================================================

def inputgene(i):
	fin = open(i,'r')     
	for j in fin:
		j = j.strip("\n")
		userinput  = j.split(',')
		for k in range(0,len(userinput)):
			userinput[k] = userinput[k].strip()
#			print (userinput)
	return userinput					## RETURNS USERINPUT IN NEWLINE

#======================================================================

#### %*%*%*% FOR RUNNING THE CODE MULTIPLE TIMES on DIFFERENT INPUT FILES %*%*%*%  #######

#======================================================================

'''data = ("YOURDATABASE.txt") ## in the form of field (GOID/GOTERM)\t(tab separated genes) 
for root, dir, files in os.walk("C:\\Users\\path \\to \\your \\input \\files"):
	for file in files:
		if file.endswith('.txt'):  ##change according to your input files extension
			inputg = file
			name1=inputg.split(".")[0]
			print (name1)
			output=inputg.split(".")[0]+"_out.csv"
			hypergeometry(data,inputg,output)
			print (output.split(".")[0]+"\n**********Done************")
			#print(len(i1))'''
#======================================================================

#### %*%*%*% FOR SINGLE RUN WITH SYS.ARGVs %*%*%*%  #######

#======================================================================

data = ("db.txt") #sys.argv[1] 						## DATABASE
inputg = sys.argv[1]					## INPUT GENESET
output = sys.argv[1].split(".")[0]+"_out.csv"#sys.argv[1].split(".")[0]+"_out.csv"					##CONSUMING MORE TIME ## OUTPUT; FORMAT: NAME.CSV BUT AUTOMATICALLY TAKES INPUT FILE NAME AS OUTFILE NAME_OUT.CSV
#output = sys.argv[3]					## OUTPUT; FORMAT: NAME.CSV
hypergeometry(data,inputg,output)

print ("time taken: ---%s seconds---"%(time.time()-start_time))