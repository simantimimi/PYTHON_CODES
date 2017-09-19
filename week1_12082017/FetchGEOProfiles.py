import urllib2
from urllib.request import urlopen

#__author__ = 'prabhasa.ravikirthi'

from Bio import Entrez
import xml.etree.ElementTree as ET
import csv
import time

gene_list=['duchenne muscular dystrophy']


outfile=csv.writer(open("try.csv",'w'),delimiter=",")
outfile.writerow(["GDS accesion number","Gene Symbol","Experiment Title","Organism","Phenotype"])
Entrez.email = "simanti.bhattacharya@gmail.com"
for query_gene in gene_list:
    #AND Homo sapiens[Organism]
    search_handle = Entrez.esearch(db="geoprofiles"
                               ,term=query_gene
                               #,term=query_gene+"[Gene Symbol] AND up down genes[filter]"
                               # ,term="ptgs2[Gene Symbol] AND Homo sapiens[Organism]"
                               ,usehistory=True
                               # ,retmax=2
                                )
    #print(search_handle)
    records=Entrez.read(search_handle,validate=True)
    # print records
    #print(records)
    if records['Count']!='0':
        # print records["IdList"]
        # q=','.join(records["IdList"][1:3])
        # print q
        # data_handle=Entrez.esummary(db="gds",id=72368286,query_key=records["QueryKey"],WebEnv=records["WebEnv"])
        # data=Entrez.read(data_handle)

        # for d in data["DocumentSummarySet"]["DocumentSummary"]:
        #     print d["GDS"], d["groups"]
            # print d

        profile_handle=Entrez.esummary(db="geoprofiles"
                               # ,id='94412796,94400592')
                               # ,retmode='xml')
                               # ,retmax=2
                               ,query_key=records["QueryKey"],WebEnv=records["WebEnv"])
        profile=Entrez.read(profile_handle)
        # print profile
        print("Number of records fetched for ", query_gene, " : ",
              len(profile["DocumentSummarySet"]["DocumentSummary"]))

        genep=[]
        for p in profile["DocumentSummarySet"]["DocumentSummary"]:
            info=[p["GDS"],p["geneName"],p["title"],p["taxon"]]
            while True:
                try:
                    state_handle=Entrez.esummary(db="gds",id=p["GDS"])
                    state=Entrez.read(state_handle)
                    break
                except urllib2.URLError:
                    print("passed")
                    pass
                    # state_handle=Entrez.esummary(db="gds",id=p["GDS"])
                    # state=Entrez.read(state_handle)
            time.sleep(1)
            if len(state)>1: print(p["GDS"], len(state))
            info.append(state[0]["subsetInfo"])
            # print p["GDS"]
            if info not in genep:
                #print info
                genep.append(info)
        print("number of unique gene-disease associations for", query_gene, " : ", len(genep))

        outfile=csv.writer(open("E:\\geo\GEOprofiles_ALL_03.12.2015.csv",'ab'),delimiter=",")
        # outfile.writerow(["GDS accesion number","Gene Symbol","Experiment Title","Organims","Phenotype"])
        for g in genep:
            g[2]=g[2].encode('ascii','ignore')
            print(g)
            outfile.writerow(g)
    else:
        print("Number of records fetched for ", query_gene, " : 0")
