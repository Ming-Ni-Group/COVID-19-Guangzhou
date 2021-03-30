#!/usr/local/bin/python
# -*- coding: utf-8 -*-  
import os,sys,linecache
import numpy as np



def readCt(Person_SIDF):
	sampleCtdic = {}
	Person_SIDFls = open(Person_SIDF).readlines()
	for Person_SIDFl in Person_SIDFls:
		if "Patient" not in Person_SIDFl and Person_SIDFl != "\n":
			sample_SID = Person_SIDFl.split("\t")[0]
			Ct = Person_SIDFl.split("\t")[4].split("\n")[0]
			sampleCtdic[sample_SID] = Ct		
	return sampleCtdic


def fileLst(depthpath):
	filenameLst=[]
	for filename in os.listdir(depthpath):
		filenameLst.append(depthpath + "/" + filename)
	return filenameLst

def outHead():	
	headLst = []
	headLst.append("sample")
	headLst.append("Ct")
	#headLst.append("flag")
	headLst.append("Cov10X")
	#headLst.append("perct_20X")
	#headLst.append("perct_50X")
	#headLst.append("mean")
	#headLst.append("median")
	#headLst.append("dp_1_4")
	#headLst.append("dp_1_2")	
	#headLst.append("dp_3_4")
	head = "\t".join(headLst)
	return head



def CovStat(samp,defF,sample_Ctdic):
	genomLen = 29903
	posiLst = []
	dayu_10Xcount = 0
	dayu_20Xcount = 0
	dayu_50Xcount = 0
	defFls = open(defF).readlines()
	for defFl in defFls:
		posi = defFl.split("\t")[1]
		dep = int(defFl.split("\t")[2])
		posiLst.append(dep)
		if dep > 10:
			dayu_10Xcount += 1
	dayu_10X_Percent = round(float(dayu_10Xcount)*100/genomLen,2)

	a =	posiLst

	mean = round(np.mean(a),2)
	median = round(np.median(a),2)
	dp_1_4=int(np.percentile(a,25))
	dp_1_2=int(np.percentile(a,50))
	dp_3_4=int(np.percentile(a,75))
	OUTLst = []
	
	OUTLst.append(str(samp))
	OUTLst.append(str(sample_Ctdic[samp]))
	OUTLst.append(str(dayu_10X_Percent))
	#OUTLst.append(str(dayu_20X_Percent))
	#OUTLst.append(str(dayu_50X_Percent))
	#OUTLst.append(str(mean))
	#OUTLst.append(str(median))
	#OUTLst.append(str(dp_1_4))
	#OUTLst.append(str(dp_1_2))	
	#OUTLst.append(str(dp_3_4))

	out = "\t".join(OUTLst)
	return out




