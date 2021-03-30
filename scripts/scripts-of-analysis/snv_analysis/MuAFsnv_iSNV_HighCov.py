#!/usr/local/bin/python
# -*- coding: utf-8 -*-  
import sys,os,linecache
from optparse import OptionParser
import numpy as np

import MuAFsnv_iSNV_HighCov_module  as shared_snvMOD
#1


# The script is used to view the MuAF distribution  of synonymous and nonsynonymous annotations  of singleton and shared iSNV sites from the iSNVtable.
################################################################################################################################
#
# For help as a standalone program type: python MuAFsnv_iSNV_Ct29.py  -h
#
# Examples:
#	   python MuAFsnv_iSNV_Ct29.py  -i  ../../../data/iSNV_with_SNP.all.txt  -t  ../../../data/sampCtValuefile.txt  -e  ../../../data/allsampMutatAnno.txt  -o ../../../result/Fig4c-single_share.snvMuAF.txt
#		
#	
################################################################################################################################



def main():
	flagLst = ["share-","single-"]
	if (os.path.exists(boxplot_outF)):
		os.remove(boxplot_outF)
	boxplot_head = shared_snvMOD.boxplot_outHead()
	boxplot_outFO = open(boxplot_outF,'a')
	boxplot_outFO.write(boxplot_head + "\n")
	boxplot_outFO.close()
	seqIDDic = shared_snvMOD.sampleID(ID_seqIDF)[0]
	L1L2_samps = shared_snvMOD.sampleID(ID_seqIDF)[1]
	print L1L2_samps
	sampleCtdic = shared_snvMOD.sampCt(sampsCtF,minCtValue,maxCtValue)
	for flag in flagLst:	
		singleXins = shared_snvMOD.outAllSampFStat(allSampsAnno,sampleCtdic,flag,seqIDDic,L1L2_samps)
		sampsLst = sampleCtdic.keys()
		PosiPosi = singleXins[0]
		print PosiPosi
		PosiDict = singleXins[1]
		shareNumDic = singleXins[2]
		
		effDic = singleXins[3]
		lieDic = shared_snvMOD.iSNVlieDic(iSNV,seqIDDic)
			
		boxout = shared_snvMOD.boxplot_syn_miss(sampsLst,lieDic,iSNV,effDic,shareNumDic,lowFreq,upFreq,PosiPosi,flag)
		boxplot_outFO = open(boxplot_outF,'a')
		boxplot_outFO.write(boxout + "\n")
		boxplot_outFO.close()



if __name__ == '__main__':
	################################################################################################################################
	# Parameters
	################################################################################################################################
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i","--inputfile",
	                  dest = "inputfile",
	                  default = "",
	                  metavar = "file",
	                  help = "Path to input iSNV Table [required]")
	parser.add_option("-t","--sampCtfile",
	                  dest = "sampCtfile",
	                  default = "",
	                  metavar = "file",
	                  help = "Path to Ct value Table [required]")
	parser.add_option("-e","--effsnp",
	                  dest = "effsnp",
	                  default = "",
	                  metavar = "file",
	                  help = "Path to All Samp snpeff stat file [required]")

	parser.add_option("-s","--ID_seqIDF",
	                  dest = "ID_seqIDF",
	                  default = "",
	                  metavar = "file",
	                  help = "Path to All Samp seqID and ID [required]")

	parser.add_option("-o","--boxplot_outF",
	                  dest = "boxplot_outF",
	                  default = "",
	                  metavar = "file",
	                  help = "Name of output file  [required]")

	parser.add_option("-m","--minMuAF ",
	                  dest = "minMuAF",
	                  default = "0.05",
	                  metavar = "",
	                  help = "Minimum MuAF for the iSNV cite to be counted.[optional]")

	parser.add_option("-M","--maxMuAF",
	                  dest = "maxMuAF",
	                  default = "0.95",
	                  metavar = "",
	                  help = "Maximum MuAF for the iSNV cite to be counted.[optional]")

	parser.add_option("-c","--minCtValue",
	                  dest = "minCtValue",
	                  default = "10",
	                  metavar = "",
	                  help = "Minimum Ct Value  for the samples to be counted,Ct value >= minCtValue.[optional]")

	parser.add_option("-C","--maxCtValue",
	                  dest = "maxCtValue",
	                  default = "29",
	                  metavar = "",
	                  help = "Maximum Ct Value  for the samples to be counted,Ct value < maxCtValue.[optional]")



	(options,args) = parser.parse_args()

	iSNV           = os.path.abspath(options.inputfile)
	sampsCtF       = os.path.abspath(options.sampCtfile)	
	allSampsAnno   = os.path.abspath(options.effsnp)
	ID_seqIDF      = os.path.abspath(options.ID_seqIDF)
	boxplot_outF   = options.boxplot_outF
	lowFreq        = options.minMuAF
	upFreq         = options.maxMuAF
	minCtValue     = options.minCtValue
	maxCtValue     = options.maxCtValue

	
	
 


	main()
