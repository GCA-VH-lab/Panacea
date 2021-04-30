__author__		= "Chayan Kumar Saha"
__copyright__	= "MIT License: Copyright (c) 2020 Chayan Kumar Saha"
__email__		= "chayan.sust7@gmail.com"

import argparse
import os, sys, os.path

usage= '''  Description: Filter prdicted toxins based on Abundance '''

parser = argparse.ArgumentParser(description=usage)
parser.add_argument("-i", "--input", required=True, help="*_TreeOrderOperonTA.txt")
parser.add_argument("-d", "--desc", required=True, help="*_outdesc.txt")
parser.add_argument("-l", "--list", required=True, help="List of input that was used for TAGs.py")
parser.add_argument("-t", "--threshold", help="threshold here is a number of different species that have same TA cluster, default=2")
parser.add_argument("-v","--version", action="version", version='%(prog)s 1.0.0')
args = parser.parse_args()
parser.parse_args()


if args.threshold:
	thresh=args.threshold
else:
	thresh="2"


taxaGCF={}
serialTaxa={}
i=0
with open(args.list, 'r') as lIn:
	for line in lIn:
		Line=line.rstrip().split('\t')
		i+=1
		taxa=Line[1]+'#'+str(i)+'_'+Line[2]
		taxaGCF[taxa]=Line[0]
		serialTaxa[str(i)]=taxa

#print(taxaGCF)#'WP_090558408.1#2281_Bacillus_subtilis_Ia1a': 'GCF_900095345.1'
#print(len(taxaGCF))#2281


def cognateFind(item1,item2):
	item1List=item1.split('|')
	item1List.remove(item2)
	return ''.join(map(str,item1List))


def speciesCheck(item1, item2):
	query=item1.split('#')[0]
	species=item1.split('#')[1][item1.split('#')[1].index('_')+1:]
	TaPairList=item2.split('|')
	IndexQuery=TaPairList.index(query)
	clustIndex=''
	if IndexQuery==0:
		clustIndex=1
	if IndexQuery==1:
		clustIndex=0
	ClustAccession=TaPairList[int(clustIndex)]
	return ClustAccession+'\t'+species

speciesDict={}
LineList=[]
cognateSet=set()
#j=0
with open(args.input, 'r') as TAin:
	for line in TAin:
		if line[0]!='#':
			Line=line.rstrip().split('\t')
			LineList.append(Line)
			cognateSet.add(cognateFind(Line[3],Line[1])+'\t'+Line[0].split('#')[1].split('_')[0])
			speciesDict[speciesCheck(Line[0],Line[3]).split('\t')[0]+'\t'+taxaGCF[Line[0]]]=taxaGCF[Line[0]]


def getTA(item):
	query=item[1]
	queryOrder=item[3].split('|').index(query)
	TaOrder=''
	if queryOrder==0:
		TaOrder=1
	elif queryOrder==1:
		TaOrder=0
	TAfam=item[5].split('|')[int(TaOrder)]
	occuranceNum=(item[6].split('/')[0])
	return str(TAfam)+'\t'+occuranceNum

occuranceSet=set()
TADict={}#TAGroup(key)=Clus\tOccur
for items in LineList:
	TADict[items[5]]=getTA(items)
	occuranceSet.add(int(getTA(items).split('\t')[1]))


TopList=sorted(occuranceSet,reverse=True)


clusterName=set()
for items in TopList:
	for ids in TADict:
		if int(TADict[ids].split('\t')[1])==items:
			clusterName.add(int(TADict[ids].split('\t')[0]))

TopCluster=sorted(clusterName)

clusList=[]

with open(args.desc, 'r') as desIn:
	for line in desIn:
		if '\t' in line:
			Line=line.rstrip().split('\t')
			#print(Line[0])
			for cluster in TopCluster:
				#print(cluster)
				if Line[0].split('(')[0]+'('==str(cluster)+'(':
					clusterACC=str(cluster)+'#'+Line[1]
					clusList.append(clusterACC)


clustDict={}
for item in TopCluster:
	Clus_accList=[]
	for element in clusList:
		if element.split('#')[0]==str(item):
			Clus_accList.append(element.split('#')[1])
	clustDict[str(item)]=Clus_accList


speciesClusDict={}
for ids in clustDict:
	speciesClustSet=set()
	for items in clustDict[ids]:
		for keys in speciesDict:
			if keys.split('\t')[0]==items:
				speciesClustSet.add(speciesDict[keys])
	speciesClusDict[ids]=len(speciesClustSet)

print('#Cluster', 'Discarded_predicted_PanT', 'Taxa', sep='\t')
for ids in speciesClusDict:
	if int(speciesClusDict[ids])<int(thresh):
		for items in clustDict[ids]:
			for cogitems in cognateSet:
				if items==cogitems.split('\t')[0]:
					print(ids, items, serialTaxa[cogitems.split('\t')[1]], sep='\t')
