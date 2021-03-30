# -*- coding: utf-8 -*-
"""
Created on Sun Apr 12 15:49:57 2020

@author: Jiayong Zhong

get nCOV19 iSNVs from MAFFT alignments

# Command version
python nCoVSNV.py <alignmfile> <Postionfile> <refHeader>

"""
import sys
from Bio import SeqIO
import pandas as pd
#sys.argv[1] alignmfile
#sys.argv[2] Positionfile
#sys.argv[3] refHeader="EPI_ISL_402125"
alignmfile, Positionfile, refHeader = sys.argv[1], sys.argv[2], sys.argv[3]
# Loading fasta 
alignmDic = SeqIO.to_dict( SeqIO.parse(alignmfile, "fasta") )
# load position
with open(Positionfile, "r") as fileID:   
    positions = fileID.readlines()
positions = [int(i.rstrip("\n")) for i in positions]
psite = ["p.%d"%i for i in positions]
## alignment file check: all of the  alignment seq are the same length
seqlen=[len(seqstr) for seqstr in alignmDic.values()]
if max(seqlen)==min(seqlen):
    print("Alignment check!")

# refheader 
refseq = str(alignmDic[refHeader].seq)

# refseq position
refposition, ALLalignposition = [], []
gap_P=-1
P=0 # first position 0
for i in range(len(refseq)):
    if refseq[i] != "-":
        refposition.append(P)
        ALLalignposition.append(i)
        P +=1
    else:
        refposition.append(gap_P)
               
# find position function
def FindPostion(find_P, alignmDic, seqkeys, refposition):
    find_P = find_P-1 # 0 zero 
    p_index = refposition.index(find_P)
    # keys orders
    sk = list( seqkeys.values() )[0]
    # get position in other seq
    p_base = []
    for skey in sk:
        p_base.append( str.upper( alignmDic[skey][p_index] ) )
    return(p_base)

# loop find
P_base_list = []
seqkeys = {"Header":list(alignmDic.keys())}
### specific postion in postion
for find_P in positions:
    p_base = FindPostion(find_P, alignmDic, seqkeys, refposition)
    P_base_list.append(p_base)

# export
export_p_base_df = dict(zip(psite, P_base_list) )
## merge and export
export_p_base_df = pd.DataFrame( dict(seqkeys, **export_p_base_df) ) 

# all postions of samples
AllDic = dict()
allkeys = [sk for sk in alignmDic.keys()]
allkeys.insert(0, refHeader)
allpositionlist = []
for skey in allkeys:
    seqstr = str.upper( str(alignmDic[skey].seq) )
    seqstr = [seqstr[i] for i in ALLalignposition]
    allpositionlist.append(seqstr)
AllDic = dict(zip(allkeys, allpositionlist))
Allpf = pd.DataFrame(AllDic)
### export
exportfile = alignmfile.replace(".fa", "piSNVs.csv")
export_p_base_df.to_csv(exportfile, index=False)
allexportfile = alignmfile.replace(".fa", "_all_iSNVs.csv")
Allpf.to_csv(allexportfile, index=False)
