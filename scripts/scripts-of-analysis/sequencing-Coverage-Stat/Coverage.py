#!/usr/local/bin/python
# -*- coding: utf-8 -*-  
import os,sys,linecache,time
from optparse import OptionParser
import numpy as np
import Coverage_module

#The script is to count the genomic coverage of each sample 
################################################################################################################################
#
# For help as a standalone program type: python Ct_Coverage.py   -h
#
# Examples:
#	python Coverage.py   -i  data/person200/  -c  data/sampCtValuefile.txt  -o  result/CoverageVSCt.txt
#	#parameters:
#		-i  Path to depth files generated from 'samtools depth' commond for all samples.
#		-c  Input file including the  Ct value of each sample
#		-o  Output file ,Statistics of Genome coverage of loci with output depth >= 10X for all samples
#		
################################################################################################################################


def main():
################################################################################################################################
# Parameters
################################################################################################################################
	usage = "usage:python  %prog [option]"
	parser = OptionParser(usage=usage)

	parser.add_option("-i","--inputpath",
	                  dest = "inputpath",
	                  default = "",
	                  metavar = "path",
	                  help = "Path to depth files generated from 'samtools depth' commond for all samples.  [required]")

	parser.add_option("-o","--outputfile",
	                  dest = "outputfile",
	                  default = "",
	                  metavar = "file",
	                  help = "output file,Statistics of Genome coverage of loci with output depth >= 10X for all samples [required]")
	parser.add_option("-c","--inputCtfile",
	                  dest = "inputCtfile",
	                  default = "",
	                  metavar = "file",
	                  help = "Input file of Ct value [required]")
	(options,args) = parser.parse_args()
	depthpath       = os.path.abspath(options.inputpath)
	sampsCtF       = os.path.abspath(options.inputCtfile)
	outCoverageStatF = os.path.abspath(options.outputfile)

	startTime = time.time()
	#sampsCtF = "/home/liuhj/liuhj2_44/nCov_202004/Fig/Ct_Coverage/sampsCt.txt"
	sample_Ctdic = Coverage_module.readCt(sampsCtF)
	print "Start"
	print "input path :  " + depthpath
	print "output file   : " + outCoverageStatF
	head=Coverage_module.outHead()
	if (os.path.exists(outCoverageStatF)):
		os.remove(outCoverageStatF)
	defFO = open(outCoverageStatF,'a')
	defFO.write(head + "\n")
	defFO.close()

	filenameLst = Coverage_module.fileLst(depthpath)
	for depthF in filenameLst:
		samp = depthF.split("/")[-1].split(".depth.txt")[0]
		print  "Calculating :   " + samp
		sampCovOut = Coverage_module.CovStat(samp,depthF,sample_Ctdic)
		defFO = open(outCoverageStatF,'a')
		defFO.write(sampCovOut + "\n")
		defFO.close()

	endTime = time.time()
	print "Total samples: " + str(len(filenameLst))
	sys.stdout.write("Total time taken: "+str(endTime-startTime)+" seconds\n")



if __name__ == "__main__":
	main()

