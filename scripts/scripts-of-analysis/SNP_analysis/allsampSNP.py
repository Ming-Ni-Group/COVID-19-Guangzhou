#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys,os,time,linecache
from optparse import OptionParser
import allsampSNP_module  as geneM


# This script is used to count the base substitution types for singleton and shared SNP sites. 
################################################################################################################################
#
# Can be used as a standalone program or functions can be called from other Python scripts
# For help as a standalone program type: python allsampSNP.py  -h  or python allsampSNP.py  -help
#
# Examples:
#	    python allsampSNP.py -i  data/iSNV_with_SNP.all-Q20.txt   -t  data/samples-Ct-Value  -e  data/allsampMutatAnno.txt  -g share -o result/
#		python allsampSNP.py -i data/iSNV_with_SNP.all-Q20.txt   -t  data/samples-Ct-Value  -e  data/allsampMutatAnno.txt  -g single -o result/
#		python allsampSNP.py -i data/iSNV_with_SNP.all-Q20.txt   -t  data/samples-Ct-Value  -e  data/allsampMutatAnno.txt  -g shareAndsingle -o result/
#
#		 
# Notes:
#    -  Everything is calculated against the reference. 
#	 -  SNP :Freq >= 0.95
#	
################################################################################################################################




def main():
	startTime = time.time()

	SAMPIDs = geneM.sample_ID(iSNVTable)
	sampleCtdic = geneM.sampCt(sampsCtF)
	print "samples Num:--------------" + str(len(sampleCtdic))
	print 
	print 


	singleXins = geneM.outAllSampFStat(outAllSampF,sampleCtdic,single_shareFlag)
	singlePosL = singleXins[0]
	singlePosD = singleXins[1]

	geneLenD = geneM.RegionLen()
	Statouts = geneM.AllSampFStat(outAllSampF,singlePosL,singlePosD)
	Htu_snvDic = Statouts


	H_snv_outF = outputpath + "/Sfig1-" + single_shareFlag + "--SNP-effStat.txt" 

	Fexit = os.path.exists(H_snv_outF)
	if Fexit == True:
		os.remove(H_snv_outF)
	geneM.BaseChgnOut(Htu_snvDic,geneLenD,H_snv_outF)


	endTime = time.time()
	print 
	sys.stdout.write("Total time taken: "+str(endTime-startTime)+" seconds\n")


if __name__ == "__main__":
	################################################################################################################################
	# Parameters
	################################################################################################################################
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i","--inputfile",
	                  dest = "inputfile",
	                  default = "",
	                  metavar = "path",
	                  help = "Path to input iSNV Table [required]")
	parser.add_option("-t","--sampCtfile",
	                  dest = "sampCtfile",
	                  default = "",
	                  metavar = "path",
	                  help = "Input file of samples' Ct values [required]")

	parser.add_option("-e","--effsnp",
	                  dest = "effsnp",
	                  default = "",
	                  metavar = "path",
	                  help = "Path to All Samp snpeff stat file [required]")

	parser.add_option("-g","--single_shareFlag",
	                  dest = "single_shareFlag",
	                  default = "",
	                  action = "store",
	                  type = "string",
	                  help = "Flag of snv-sample type : single | share | shareAndsingle  [required]")

	parser.add_option("-o","--outputpath",
	                  dest = "outputpath",
	                  default = "",
	                  metavar = "path",
	                  help = "Path to save output files in [required]")

	(options,args) = parser.parse_args()

	iSNVTable       = os.path.abspath(options.inputfile)
	sampsCtF       = os.path.abspath(options.sampCtfile)	
	outAllSampF      = os.path.abspath(options.effsnp)
	single_shareFlag = options.single_shareFlag
	outputpath       = os.path.abspath(options.outputpath)
	

	main()












