#!/usr/bin/python
# -*- coding: utf-8 -*-

def sample_ID(iSNVTable):
	iSNVTablels = open(iSNVTable,'r').readlines()
	for iSNVTablel in iSNVTablels:
		if "pos	" in iSNVTablel:
			sampIDs = iSNVTablel.split("\n")[0].split("\t")[2:-1]			
	return sampIDs 


def RegionLen():
	refGeneDic = {"orf1a":"266..13468","orf1b":"13468..21555","S":"21563..25384"\
	,"ORF3a":"25393..26220","E":"26245..26472","M":"26523..27191",\
	"ORF6":"27202..27387","ORF7a":"27394..27759","ORF8":"27894..28259",\
	"N":"28274..29533","ORF10":"29558..29674",}
	geneLenD = {}
	codingLen = 0
	for refGene in refGeneDic:
		geneIDStrt = int(refGeneDic[refGene].split("..")[0])
		geneIDEnd = int(refGeneDic[refGene].split("..")[1])
		genelen = geneIDEnd - geneIDStrt + 1
		geneLenD[refGene] = genelen
		codingLen += genelen

	NoncodingLen = 29903-codingLen
	geneLenD['Coding'] = codingLen
	geneLenD['Non-Coding'] = NoncodingLen
	return geneLenD




def sampCt(sampsCtF):
	sampleCtdic = {}
	sampsCtFls = open(sampsCtF).readlines()
	for Person_SIDFl in sampsCtFls:
		if "Patient" not in Person_SIDFl and Person_SIDFl != "\n":
			sample_SID = Person_SIDFl.split("\t")[0]
			Ct = Person_SIDFl.split("\t")[4].split("\n")[0]
			
			sampleCtdic[sample_SID] = Ct
	return sampleCtdic
	


def outAllSampFStat(outAllSampF,sampleCt25dic,single_shareFlag):	
	Ct25PosiSampDic = {}
	
	outAllSampFls = open(outAllSampF).readlines()
	for outAllSampFl in outAllSampFls:
		if "sample" not in outAllSampFl and outAllSampFl != "\n":
			linekey = outAllSampFl.split("\n")[0].split("\t")
			sample = linekey[0].split("/")[-1]
			if sample in sampleCt25dic.keys():
				posi = linekey[1]
				Freq = linekey[17]
				if int(posi) < 29782 and int(posi) > 86 and int(posi) != 12436:
					if Freq != "NO" and Freq != "NA" and float(Freq) >= 95:						
						if posi not in Ct25PosiSampDic:				
							Ct25PosiSampDic[posi] = [sample]
						else :					
							Ct25PosiSampDic[posi].append(sample)

	Ct25singlePosL = []
	Ct25singlePosD = Ct25PosiSampDic
	for Ct25Posi in Ct25PosiSampDic.keys():
		if single_shareFlag == "single":			
			if len(Ct25PosiSampDic[Ct25Posi]) == 1:
				Ct25singlePosL.append(int(Ct25Posi))
			else:
				Ct25singlePosD.pop(Ct25Posi)
		if single_shareFlag == "share":
			if len(Ct25PosiSampDic[Ct25Posi]) > 1:
				Ct25singlePosL.append(int(Ct25Posi))
			else:
				Ct25singlePosD.pop(Ct25Posi)
		if single_shareFlag == "shareAndsingle":
			if len(Ct25PosiSampDic[Ct25Posi]) >= 1:
				Ct25singlePosL.append(int(Ct25Posi))
			else:
				Ct25singlePosD.pop(Ct25Posi)

	Ct25singlePosL.sort()
	return Ct25singlePosL,Ct25singlePosD



def Htu(Htu_snvDic,Annotation,HGVS_c):
	if Annotation in ['missense_variant','stop_gained']:
		AnnFlag = "Nonsyn"
	elif Annotation in ['synonymous_variant']:
		AnnFlag = "Syn"
	else:
		AnnFlag = "Non-Coding"
	if HGVS_c not in Htu_snvDic:
		Htu_snvDic[HGVS_c] = {}
		if AnnFlag not in Htu_snvDic[HGVS_c]:
			Htu_snvDic[HGVS_c][AnnFlag] = 1
	else:
		if AnnFlag not in Htu_snvDic[HGVS_c]:
			Htu_snvDic[HGVS_c][AnnFlag] = 1
		else:
			Htu_snvDic[HGVS_c][AnnFlag] += 1
	return Htu_snvDic						


		
def AllSampFStat(outAllSampF,Ct25singlePosL,Ct25singlePosD):
	Htu_snvDic = {}
	outAllSampFls = open(outAllSampF).readlines()
	for outAllSampFl in outAllSampFls:
		if "sample" not in outAllSampFl and outAllSampFl != "\n":			
			linekey = outAllSampFl.split("\n")[0].split("\t")
			sample = linekey[0].split("/")[-1]
			posi = linekey[1]

			if int(posi) in Ct25singlePosL and sample in Ct25singlePosD[posi] :	 			
				ref = linekey[2]
				alle = linekey[3]
				if ref == 'T':
					ref = "U"
				if alle == 'T':
					alle = "U"
				Annotation = linekey[4]
				Gene_Name = linekey[5]
				Feature_Type = linekey[6]
				Transcript_BioType = linekey[7]
				rank = linekey[8]
				rankall = linekey[9]
				HGVS_c_1 = linekey[10]
				HGVS_c = ref + ">" + alle
				HGVS_p = linekey[11]
				cDNA_pos = linekey[12]
				cDNA_posall = linekey[13]
				AA_pos = linekey[14]
				AA_posall = linekey[15]
				charRank = linekey[16]
				Freq = linekey[17]
				
				if Gene_Name == "orf1ab":
					if int(posi) >= 266 and int(posi) <= 13468:
						Gene_Name = "orf1a"
					if int(posi) >= 13468 and int(posi) <= 21555:
						Gene_Name = "orf1b"

				if float(Freq) >= 95:   
					Htu_snvDic = Htu(Htu_snvDic,Annotation,HGVS_c)

	return Htu_snvDic


def BaseChgnOut(Htu_Dic,geneLenD,H_outF):  ##isnv/kb
	H_outFO = open(H_outF,'a')
	H_outFO.write("baseChange" +  "\t" + "eff" +  "\t" + "Count"+ "\n")
	H_outFO.close()
	outLst = []
	for Htu_key in Htu_Dic.keys():			
		for Charposi in ['Nonsyn','Syn',"Non-Coding"]:
			lineLst = []
			lineLst.append(Htu_key)
			lineLst.append(Charposi)
			if Charposi in  Htu_Dic[Htu_key]:
				CharposiCout = Htu_Dic[Htu_key][Charposi]
			else:
				CharposiCout = 0
			lineLst.append(str(CharposiCout))				
			#lineLst.append(str(round(float(CharposiCout)/geneLenD[Charposi],3)))
			outline = "\t".join(lineLst)
			outLst.append(outline)
	out = "\n".join(outLst)
	H_outFO = open(H_outF,'a')
	H_outFO.write(out+ "\n")
	H_outFO.close()



def RplotsingleData(outAllSampF,Ct25singlePosL,Ct25singlePosD,RplotsingleDataF,Flag):	
	outAllSampFls = open(outAllSampF).readlines()
	for outAllSampFl in outAllSampFls:
		if "sample" not in outAllSampFl and outAllSampFl != "\n":
			linekey = outAllSampFl.split("\n")[0].split("\t")
			sample = linekey[0]
			posi = linekey[1]
			Freq = linekey[17]
			ref = linekey[2]
			alle = linekey[3]
			cChange=  ref + ">" + alle

			if int(posi) in Ct25singlePosL and  sample in Ct25singlePosD[posi]:								
				outsinglesnpeffFO = open(RplotsingleDataF,'a')
				outsinglesnpeffFO.write(posi + "\t" + sample + "\t" + Flag  + "\t" + Freq  +"\t" + cChange + "\n")
				outsinglesnpeffFO.close()
				









