#!/usr/bin/env python
#coding=utf-8
import sys,collections,gzip,pickle
import os
from tqdm import tqdm
import pandas
from multiprocessing import Pool as ProcessPool
import itertools
from ete4 import Tree
import networkx
from networkx.algorithms.components.connected import connected_components
import SupplementSrc
import numpy as np
from Bio import SeqIO
import argparse
def GetSyntenic(MafFile,Prefix):
	Syntenicd={}
	#Syntenicd=collections.defaultdict(list)
	Synl=[]
	i=0
	MafHandle=open(MafFile,"rb")
	MagicNumber=MafHandle.read(2)
	if MagicNumber == b'\x1f\x8b':
		MafHandle=gzip.open(MafFile)
	for x in MafHandle:
		x=x.decode().rstrip()
		if x.startswith("#") or (not x):
			continue
		l=x.split()
		if l[0]=="a":
			if Synl:
				BlockName=Prefix+"_block%s"%i
				#print(BlockName,len(Synl))
				i+=1
				for si in Synl:
					sp=si[1].split(".")[0] #sp里面不能有"."
					ch=".".join(si[1].split(".")[1:])
					strand=si[4]
					seq=si[-1]
					if strand=="+":
						start=int(si[2])
						end=int(si[2])+int(si[3])
					elif strand=="-":
						start=int(si[5])-int(si[2])-int(si[3]) #$6-$3-$4
						end=int(si[5])-int(si[2]) #$6-$3
					else:
						print("strand err")
						sys.exit()
					if not sp in Syntenicd:
						Syntenicd[sp]=collections.defaultdict(list)
					Syntenicd[sp][BlockName].append((ch,start,end,strand,seq))
			Synl=[]
		if l[0]=="s":
			Synl.append(l)
	BlockName=Prefix+"_block%s"%i
	#print(BlockName,len(Synl))
	for si in Synl:
		sp=si[1].split(".")[0]
		ch=".".join(si[1].split(".")[1:])
		strand=si[4]
		seq=si[-1]
		if strand=="+":
			start=int(si[2])
			end=int(si[2])+int(si[3])
		elif strand=="-":
			start=int(si[5])-int(si[2])-int(si[3]) #$6-$3-$4
			end=int(si[5])-int(si[2]) #$6-$3
		else:
			print("strand err")
			sys.exit()
		if not sp in Syntenicd:
			Syntenicd[sp]=collections.defaultdict(list)
		Syntenicd[sp][BlockName].append((ch,start,end,strand,seq))
	for sp in Syntenicd:
		Fr=open(Prefix+"_%s_synt.bed"%sp,"w")
		for BlockName in Syntenicd[sp]:
			for ch,start,end,strand,seq in Syntenicd[sp][BlockName]:
				Fr.write(ch+"\t"+str(start)+"\t"+str(end)+"\t"+BlockName+"\t.\t"+strand+"\n")
		Fr.close()
	return Syntenicd
def GetSyntenicGene(Syntenicd,Fd,Prefix): #need gene.bed
	BlockGened={}
	Blockd={}
	BlockList=[]
	BinfoList=[]
	for sp in Syntenicd.keys():
		BlockGened[sp]=collections.defaultdict(list)
		lines=os.popen("bedtools intersect -loj -a %s_%s_synt.bed -b %s"%(Prefix,sp,Fd[sp]["gene"])).readlines()
		for x in lines:
			x=x.rstrip()
			l=x.split("\t")
			try:
				Blockd[l[3]]
			except KeyError:
				Blockd[l[3]]={}
			Blockd[l[3]][sp]=[]
			if l[-2]!="-1":
				BlockGened[sp][l[3]].append((int(l[7]),int(l[8]),l[9],l[-1])) #start,end,name,strand
	for sp in Syntenicd:
		for BlockName in tqdm(Syntenicd[sp],total=len(Syntenicd[sp].keys())):
			for ch,Bstart,Bend,strand,seq in Syntenicd[sp][BlockName]:
				BinfoList.append([BlockName,sp,ch,Bstart,Bend,strand,seq])
				for gStart,gEnd,gName,gStrand in BlockGened[sp][BlockName]:
					RelpStart=gStart-Bstart
					if RelpStart<0:
						RelpStart=0
					RelpEnd=gEnd-Bstart
					if RelpEnd>=Bend-Bstart+1:
						RelpEnd=Bend-Bstart
					print(BlockName,ch,Bstart,Bend)
					print(gName,gStart,gEnd)
					print(RelpStart,RelpEnd)
					print("seq",len(seq))
					trp={}
					t=0
					for ni,a in enumerate(seq):
						if a!="-":
							trp[t]=ni
							t+=1
					print(trp)
					RelpStart=trp[RelpStart]
					if RelpEnd!=0:
						RelpEnd=trp[RelpEnd-1]+1
					else:
						RelpEnd=trp[RelpEnd]+1
					print(RelpStart,RelpEnd)
					BlockList.append((BlockName,len(seq),sp,gName,RelpStart,RelpEnd,gStrand))
					Blockd[BlockName][sp].append((RelpStart,RelpEnd,gName,gStrand))
	with open("%s_Blockd.pickle"%Prefix,"wb") as f:
		pickle.dump(Blockd,f)
	DfBlock=pandas.DataFrame(BlockList,columns=['BlockName','BlockLength','Species','GeneName','Start','End','Strand'])
	DfBlock.to_csv("%s_blockgene.xls"%Prefix,sep="\t",index=False)
	DfBlockInfo = pandas.DataFrame(BinfoList, columns=['BlockName','Species','Chr','Start','End','Strand','Seq'])
	DfBlockInfo.to_csv("%s_blockinfo.xls"%Prefix,sep="\t",index=False)
	return Blockd
def GetOrthoCluster(Syntenicd,Blockd,AllGened,GenePairsScored,Blastd,ProteinSeqd,GeneToTransd,RefSp,Fd,Prefix,NoUnannotated):
	OrthoCluster={}
	GeneToBlock=collections.defaultdict(list)
	AlreadGene={}
	UsedGene=collections.defaultdict(list)
	for BlockName in tqdm(Blockd,total=len(Blockd.keys())):
		#OrthoCluster[BlockName]={}
		#Splist=list(Blockd[BlockName].keys())
		for RefStart,RefEnd,RefName,RefStrand in Blockd[BlockName][RefSp]:
			try:
				OrthoCluster[RefName]
			except KeyError:
				OrthoCluster[RefName]={}
			GeneToBlock[RefName].append(BlockName)
			UsedGene[RefSp].append(RefName)
			a=pandas.Interval(RefStart,RefEnd)
			for sp in Blockd[BlockName]:
				if sp==RefSp:
					continue
				try:
					OrthoCluster[RefName][sp]	
				except KeyError:
					OrthoCluster[RefName][sp]=[]
				for gStart,gEnd,gName,gStrand in Blockd[BlockName][sp]:
					#GeneToBlock[gName].append(BlockName,gStart,gEnd)
					b=pandas.Interval(gStart,gEnd)
					if a.overlaps(b):
						UsedGene[sp].append(gName)
						try:
							Bitscore,Pident,Nident,Qlen,_,_=Blastd[(RefName,gName)]
						except KeyError:
							Bitscore,Pident,Nident,Qlen=0,0,0,1
						try:
							gBs,gRefName=AlreadGene[gName]
							flag=1
						except KeyError:
							flag=0
						if flag==1:
							if gRefName!=RefName:
								if Bitscore>gBs:
									if OrthoCluster[RefName][sp]==[]:
										OrthoCluster[RefName][sp]=[Bitscore,Pident,Nident,Qlen,gName]			
										AlreadGene[gName]=[Bitscore,RefName]
										OrthoCluster[gRefName][sp]=[]
									else:
										rBs,_,_,_,_=OrthoCluster[RefName][sp]
										if Bitscore>rBs:
											OrthoCluster[RefName][sp]=[Bitscore,Pident,Nident,Qlen,gName]			
											AlreadGene[gName]=[Bitscore,RefName]
											OrthoCluster[gRefName][sp]=[]
						else:
							if OrthoCluster[RefName][sp]==[]:
								OrthoCluster[RefName][sp]=[Bitscore,Pident,Nident,Qlen,gName]
								AlreadGene[gName]=[Bitscore,RefName]
							else:
								rBs,_,_,_,rgName=OrthoCluster[RefName][sp]
								if Bitscore>rBs:
									OrthoCluster[RefName][sp]=[Bitscore,Pident,Nident,Qlen,gName]
									AlreadGene[gName]=[Bitscore,RefName]
									del AlreadGene[rgName]
									UsedGene[sp].remove(rgName)

	with open("%s_OrthoCluster.pickle"%Prefix,"wb") as f:
		pickle.dump(OrthoCluster,f)
	Trans={}
	UnUsedGene=collections.defaultdict(list)
	for sp in AllGened:
		UnUsedGene[sp]=list(set(AllGened[sp])-set(UsedGene[sp]))
	for RefName in OrthoCluster:
		for sp in OrthoCluster[RefName]:
			if OrthoCluster[RefName][sp]:
				Bitscore,Pident,Nident,Qlen,gName=OrthoCluster[RefName][sp]
				if GenePairsScored[RefName] and Bitscore<GenePairsScored[RefName]:
					UnUsedGene[sp].append(gName)
					print("syntenic unused",gName,RefName)
	with open("%s_UnUsedGene.pickle"%Prefix,"wb") as f:
		pickle.dump(UnUsedGene,f)
	os.system("mkdir %s_blast_out"%Prefix)
	SpList=list(Syntenicd.keys())
	OrthoParams=[]
	for RefName in tqdm(OrthoCluster,total=len(OrthoCluster.keys())):
		#print("Trans",sys.getsizeof(Trans))
		#print("OrthoParams",sys.getsizeof(OrthoParams))
		for sp in SpList: #OrthoCluster[RefName]:
			if sp==RefSp:
				continue
			for g in UnUsedGene[sp]:
				try:
					Bitscore,Pident,Nident,Qlen,_,_=Blastd[(RefName,g)]
				except KeyError:
					Bitscore,Pident,Nident,Qlen=0,0,0,1
				try:
					BestBitscore,BestPident,BestNident,BestQlen,Bg=Trans[(sp,RefName)]
					if Bitscore>BestBitscore:
						Trans[(sp,RefName)]=(Bitscore,Pident,Nident,Qlen,g)
				except KeyError:
					Trans[(sp,RefName)]=(0,0,0,1,"")
				#Trans[sp][RefName].append((Bitscore,Pident,Nident,Qlen,g))
				try:
					Bitscore,Pident,Nident,Qlen,_,_=Blastd[(g,RefName)]
				except KeyError:
					Bitscore,Pident,Nident,Qlen=0,0,0,1
				try:
					BestBitscore,BestPident,BestNident,BestQlen,Bg=Trans[(sp,g)]
					if Bitscore>BestBitscore:
						Trans[(sp,g)]=(Bitscore,Pident,Nident,Qlen,RefName)
				except KeyError:
					Trans[(sp,g)]=(0,0,0,1,"")
				#Trans[sp][g].append((Bitscore,Pident,Nident,Qlen,RefName))
			if not sp in OrthoCluster[RefName]:
				continue
			if not NoUnannotated:
				if OrthoCluster[RefName][sp]==[]:
					Fr=open("%s_blast_out/"%(Prefix)+Prefix+"_"+RefName+"_protein.faa","w")
					for ProteinId in GeneToTransd[RefName]:
						Fr.write(">"+ProteinId+"\n")
						Fr.write(ProteinSeqd[RefSp][ProteinId]+"\n")
					Fr.close()
					Fr=open("%s_blast_out/"%(Prefix)+Prefix+"_"+RefName+"_"+sp+"_region.bed","w")
					for BlockName in GeneToBlock[RefName]:
						t=0
						for ch,Bstart,Bend,strand,seq in Syntenicd[sp][BlockName]:
							RegionName=BlockName+"_"+sp+"_region%s"%t
							Fr.write(ch+"\t"+str(Bstart)+"\t"+str(Bend)+"\t"+RegionName+"\t.\t"+strand+"\n")
							t+=1
					Fr.close()
					OrthoParams.append((RefName,sp,Fd[sp]["genome"],Prefix))
	del Syntenicd,Blockd,ProteinSeqd
	TransProtein={}
	for RefName in tqdm(OrthoCluster,total=len(OrthoCluster.keys())):
		for sp in SpList:
			if sp==RefSp:
				continue
			gBs,gPi,gNi,gLen,g1=Trans[(sp,RefName)]
			if g1=="":
				print(sp,RefName)
				TransProtein[(RefName,sp)]=[0,0,0,1,None]
				continue
			rBs,rPi,rNi,rLen,r1=Trans[(sp,g1)]
			if r1==RefName and r1!="":
				TransProtein[(RefName,sp)]=[gBs,gPi,gNi,gLen,g1]
			else:
				TransProtein[(RefName,sp)]=[0,0,0,1,None]
	with open("%s_TransProtein.pickle"%Prefix,"wb") as f:
		pickle.dump(TransProtein,f) #translocation
	OrthoDNA={}
	if not NoUnannotated:
		Fr=open("GetBlastOut_log.txt","w")
		pool = ProcessPool(50)
		OrthoDNAResults=pool.starmap(GetBlastOut,OrthoParams)
		pool.close()
		pool.join()
		for r in OrthoDNAResults:
			RefName,sp,ri=r
			if ri:
				OrthoDNA[(RefName,sp)]=ri #sorted(RegionBlast)[-1]
			else:
				Fr.write(RefName+"\t"+sp+"\n")
		Fr.close()
		with open("%s_OrthoDNA.pickle"%Prefix,"wb") as f:
			pickle.dump(OrthoDNA,f) #unannot
	return OrthoCluster,OrthoDNA,TransProtein
def GetBlastOut(RefName,sp,GenomeFastaFile,Prefix): #need genome
	Cmd1="bedtools sort -i %s_blast_out/%s|bedtools merge -d 1000 -i -|awk '{t+=1;print $0\"\t%s\"t}'|bedtools getfasta -name -fi %s -bed - >%s_blast_out/%s"%(Prefix,Prefix+"_"+RefName+"_"+sp+"_region.bed",sp+"_region",GenomeFastaFile,Prefix,Prefix+"_"+RefName+"_"+sp+"_region.fasta") # -d need change
	Cmd2="miniprot -t 14 -I --trans %s_blast_out/%s %s_blast_out/%s >%s_blast_out/%s_miniprot.out"%(Prefix,Prefix+"_"+RefName+"_"+sp+"_region.fasta",Prefix,Prefix+"_"+RefName+"_protein.faa",Prefix,Prefix+"_"+RefName+"_"+sp)
	SupplementSrc.GetRightSystem(Cmd1)
	SupplementSrc.GetRightSystem(Cmd2)
	Fr=open(Prefix+"_blast_out/"+Prefix+"_"+RefName+"_"+sp+"_region_predict.faa","w")
	t=0
	for x in open("%s_blast_out/%s_miniprot.out"%(Prefix,Prefix+"_"+RefName+"_"+sp)):
		x=x.rstrip()
		if x.startswith("##STA"):
			pep=x.split()[1]
			Fr.write(">"+ProteinName+".pep%s\n"%t)
			Fr.write(pep+"\n")
			t+=1
		else:
			l=x.split("\t")
			name,_,ch,coor=l[5].split(":")
			print(name,ch,coor)
			OriginalStar,OriginalEnd=coor.split("-")
			Start=int(OriginalStar)+int(l[7])
			End=int(OriginalStar)+int(l[8])
			ProteinName=name+"::"+ch+":"+str(Start)+"-"+str(End)
			print(ProteinName)
	Fr.close()
	if t!=0:
		Cmd1="makeblastdb -in %s_blast_out/%s -out %s_blast_out/%s -dbtype prot -title %s"%(Prefix,Prefix+"_"+RefName+"_"+sp+"_region_predict.faa",Prefix,Prefix+"_"+RefName+"_"+sp+"_region_predict",Prefix+"_"+RefName+"_"+sp+"_region_predict")
		Cmd2="blastp -db %s_blast_out/%s -out %s_blast_out/%s -query %s_blast_out/%s -num_threads 10 -evalue 10 -outfmt \"6 qseqid qlen qstart qend sseqid slen sstart send qcovs bitscore evalue pident nident\""%(Prefix,Prefix+"_"+RefName+"_"+sp+"_region_predict",Prefix,Prefix+"_"+RefName+"_"+sp+"_region.out",Prefix,Prefix+"_"+RefName+"_protein.faa")
		SupplementSrc.GetRightSystem(Cmd1)
		SupplementSrc.GetRightSystem(Cmd2)
		RegionBlastTemp=collections.defaultdict(list)
		for x in open("%s_blast_out/"%(Prefix)+Prefix+"_"+RefName+"_"+sp+"_region.out"):
			x=x.rstrip()
			l=x.split("\t")
			RegionBlastTemp[(l[0],l[4])].append((float(l[-4]),float(l[-1])/float(l[1])))
		RegionBlast=[]
		for pair in RegionBlastTemp:
			a=sorted(RegionBlastTemp[pair])[-1]
			RegionBlast.append((a[0],a[1],pair[1])) #need change
		#print(RefName,sp)
		#print(RegionBlast)
		#print(sorted(RegionBlast)[-1][1:])
		if RegionBlast:
			return RefName,sp,sorted(RegionBlast)[-1][1:]
		else:
			return RefName,sp,None #blastp 比对不出结果
	else:
		return RefName,sp,None #miniprot 比对不出结果
def FilterOrthoCluster(GenePairsScored,GeneDnaPairsIdentify,OrthoCluster,OrthoDNA,TransProtein,RefSp,Prefix,NoUnannotated):
	FinalOrthoCluster={}
	for RefName in OrthoCluster:
		FinalOrthoCluster[RefName]={}
		for sp in OrthoCluster[RefName]:
			if OrthoCluster[RefName][sp]:
				Bitscore,Pident,Nident,Qlen,gName=OrthoCluster[RefName][sp]
				print(RefName,GenePairsScored[RefName],gName,Bitscore)
				if GenePairsScored[RefName] and Bitscore>=GenePairsScored[RefName]:
					FinalOrthoCluster[RefName][sp]=gName
				else:
					Bitscore,Pident,Nident,Qlen,gName=TransProtein[(RefName,sp)]
					#if Nident/Qlen*100>=TransIdentify:
					if GenePairsScored[RefName] and Bitscore>=GenePairsScored[RefName]:
						FinalOrthoCluster[RefName][sp]=gName+":translocation"
					else:
						FinalOrthoCluster[RefName][sp]=None
			else:
				FinalOrthoCluster[RefName][sp]=None
	for RefName,sp in TransProtein:
		if sp==RefSp:
			print("sp error")
			sys.exit()
		if (not sp in OrthoCluster[RefName]) or (OrthoCluster[RefName][sp]==[]):
			Bitscore,Pident,Nident,Qlen,gName=TransProtein[(RefName,sp)]
			#if Nident/Qlen*100>=TransIdentify: # may need change
			if GenePairsScored[RefName] and Bitscore>=GenePairsScored[RefName]:
				FinalOrthoCluster[RefName][sp]=gName+":translocation"
	if not NoUnannotated:
		for RefName,sp in OrthoDNA:
			Ident,gName=OrthoDNA[(RefName,sp)]
			if Ident>=GeneDnaPairsIdentify: # may need change
				if not FinalOrthoCluster[RefName][sp]:
					FinalOrthoCluster[RefName][sp]=gName
	df=pandas.DataFrame.from_dict(FinalOrthoCluster,orient="index").reset_index().rename(columns={"index":RefSp})
	print(df)
	print(len(OrthoCluster.keys()))
	print(len(FinalOrthoCluster.keys()))
	d=collections.defaultdict(list)
	for RefName in FinalOrthoCluster:
		if not FinalOrthoCluster[RefName]:
			d[RefSp].append(RefName)
	dfTemp=pandas.DataFrame.from_dict(d)
	print(dfTemp)
	df=pandas.concat([df,dfTemp])
	print(df)
	df["Group"]=["OrthoGroup%s"%i for i,RefName in enumerate(FinalOrthoCluster)]
	df.set_index("Group",inplace=True)
	df.to_csv(Prefix+"_OrthoGroupCluster.xls",sep="\t",index_label="Group")
	return df.fillna(np.nan)
def GetTreeFile(FaaFile):
	pr=FaaFile.split(".faa")[0]
	os.system("mafft --anysymbol --thread 10 %s >%s_align.faa"%(FaaFile,pr)) #L-INS-i (Probably most accurate, very slow),default:FFT-NS-2 fast
	os.system("trimal -keepheader -in %s_align.faa -out %s_align_trim.faa -automated1"%(pr,pr))
	if os.path.exists("%s_align_trim.faa"%pr):
		os.system("fasttree -quote %s_align_trim.faa >%s.tree"%(pr,pr))
	else:
		print(FaaFile)
		os.system("fasttree -quote %s_align.faa >%s.tree"%(pr,pr))
	#os.system("fasttree -quote %s_align.faa >%s.tree"%(pr,pr))
	TreeHanlde = Tree(open("%s.tree"%(pr)))
	TreeHanlde.set_outgroup(TreeHanlde.get_farthest_leaf()[0]) #rooted tree
	TreeFile="%s_deal.tree"%(pr)
	TreeHanlde.write(outfile=TreeFile)
	return TreeFile
def GetRepeatAndParalogGene(OrthoClusterDf,RefSp,Blastd,GenePairsScored,ProteinSeqd,GeneToTransd,AllGened,TransToGened,Inflation,Prefix):
	RefGeneList=list(OrthoClusterDf[RefSp])
	al=[]
	Merged={}
	AllGeneSetd={}
	ParalogList=[]
	for sp in AllGened:
		AllGeneSetd[sp]=set(AllGened[sp])
	Fr=open(Prefix+"_paralog_mcl.txt","w") #paralog
	for g1,g2 in itertools.permutations(RefGeneList,2):
		try:
			if GenePairsScored[g1] and Blastd[(g1,g2)][0]>=GenePairsScored[g1]:
				Fr.write(g1+"\t"+g2+"\t"+str(Blastd[(g1,g2)][0])+"\n") #paralog
				ParalogList.append(g1)
				ParalogList.append(g2)
		except KeyError:
			continue
	Fr.close()
	os.system("mcl %s --abc -I %s -o %s -te 50 -v all"%(Prefix+"_paralog_mcl.txt",Inflation,Prefix+"_paralog_mcl.out"))#-I越大，分的cluster越细
	OrthoClusterDf.reset_index(drop=True,inplace=True)
	OrthoClusterDf.set_index(RefSp,inplace=True)
	os.mkdir("%s_tree"%Prefix)
	RepeatList=[]
	FaaList=[]
	ParalogList=set(ParalogList)
	for i,x in enumerate(open(Prefix+"_paralog_mcl.out")):
		x=x.rstrip()
		pl=x.split("\t")
		Fr=open("%s_tree/%s_paralog_cluster%s.faa"%(Prefix,RefSp,i),"w")
		t=0
		for g in pl:
			transl=[]
			for gi in GeneToTransd[g]:
				transl.append((len(ProteinSeqd[RefSp][gi]),gi))
			gi=sorted(transl)[-1][-1]
			Fr.write(">"+gi+"\n")
			Fr.write(ProteinSeqd[RefSp][gi]+"\n")
			gd=OrthoClusterDf.loc[g,:].to_dict()
			for sp in gd:
				if pandas.isna(gd[sp]):
					continue
				if "translocation" in gd[sp]:
					gd[sp]=gd[sp].split(":")[0]
				if gd[sp] in AllGeneSetd[sp]:
					transl=[]
					for gi in GeneToTransd[gd[sp]]:
						transl.append((len(ProteinSeqd[sp][gi]),gi))
					gi=sorted(transl)[-1][-1]
					Fr.write(">"+gi+"\n")
					Fr.write(ProteinSeqd[sp][gi]+"\n")
				else:
					pseq=SeqIO.index("%s_blast_out/%s_%s_%s_region_predict.faa"%(Prefix,Prefix,g,sp),"fasta")
					Fr.write(">"+gd[sp]+"_"+str(t)+"\n")
					Fr.write(str(pseq[gd[sp]].seq)+"\n")
					t+=1
		Fr.close()
		FaaList.append("%s_tree/%s_paralog_cluster%s.faa"%(Prefix,RefSp,i))
	pool = ProcessPool(50) 
	TreeFileList=pool.map(GetTreeFile,FaaList)
	pool.close()
	pool.join()
	for TreeFile in TreeFileList:
		#print(TreeFile)
		if os.path.getsize(TreeFile)==0:
			print(0,TreeFile)
			continue
		Rl=GetRepeatGene(TreeFile,AllGeneSetd,TransToGened,Blastd,RefSp)
		#TreeHanlde = Tree(open(TreeFile))
		#pl=[]
		#for ti in TreeHanlde.leaf_names():
		#	if (not "region" in ti) and (TransToGened[ti] in AllGeneSetd[RefSp]):
		#		pl.append(TransToGened[ti])
		#pl=set(pl)
		#print(TreeHanlde)
		for ri in Rl:
			print("Repeat",ri)
			#pl=list(set(pl)-set(ri))
			RepeatList.append(ri)
			for rii in ri:
				ParalogList.remove(rii)
		#if len(pl)>=1:
			#ParalogList.append(pl)
			#print("Paralog",pl)
		#else:
		#	print("Confusing paralog",pl)
	with open("%s_RepeatList.txt"%Prefix,"w") as f:
		for r in RepeatList:
			f.write("\t".join(r)+"\n")
	#with open("%s_ParalogList.txt"%Prefix,"w") as f:
	#	for r in ParalogList:
	#		f.write("\t".join(r)+"\n")
	return RepeatList,ParalogList
def GetRepeatGene(TreeFile,AllGeneSetd,TransToGened,Blastd,RefSp):
	Rl=[]
	TreeHanlde = Tree(open(TreeFile))
	for node in TreeHanlde.traverse("postorder"):
		if node.is_leaf:
			continue
		tl=list(node.leaf_names())
		spl=[]
		gl=[]
		refgl=[]
		for t in tl:
			if t is None:
				continue
			if "region" in t:
				g=t
			else:
				g=TransToGened[t]
			gl.append(g)
			for sp in AllGeneSetd:
				if g in AllGeneSetd[sp] or g.startswith(sp):
					spl.append(sp)
					if sp==RefSp:
						refgl.append(g)
					break
		if len(set(spl))==1 and spl[0]==RefSp and len(set(gl))!=1:
			Rl.append(list(set(gl)))
			print(gl)
			print(node)
		elif spl.count(RefSp)>1:
			Rel=[]
			for ni in list(node.traverse("postorder"))[:-1]:
				if ni.dist==0:
					Rel.append(ni.name)
			if len(Rel)==len(list(node.traverse("postorder"))[:-1]):
				Rl.append(list(set(refgl)))
				print("dist==0")
				print(refgl)
				print(node)
	G = networkx.Graph()
	for ri in Rl:
		print(ri)
		G.add_nodes_from(ri)
		for g1,g2 in itertools.combinations(ri,2):
			try:
				Blastd[(g1,g2)]
				G.add_edge(g1,g2)
			except KeyError:
				try:
					Blastd[(g2,g1)]
					G.add_edge(g2,g1)
				except KeyError:
					continue
		#G.add_edges_from(list(itertools.combinations(ri,2)))
	#print(len(list(itertools.chain.from_iterable(Rl))))
	Rl=connected_components(G)
	Rl=list(Rl)
	#print(len(list(itertools.chain.from_iterable([list(ri) for ri in Rl]))))
	return Rl
def GetGeneClassify(Df,ChimerasGenesList,RepeatGeneList,ParalogGeneList,OutGroupSpList,RefSp,Prefix):
	ChimerasGenesList=set(ChimerasGenesList)
	DfTemp=Df[~Df[RefSp].isin(ChimerasGenesList)] #remove Chimeras
	#for ri in RepeatGeneList:
	#	ri2=list(set(ri)-ChimerasGenesList)
	#	if len(ri)!=len(ri2):
	#		DfTemp=DfTemp[~DfTemp[RefSp].isin(ri2)] #remove Repeat about Chimera (v8)
	#for pi in ParalogGeneList:
	#	pi2=list(set(pi)-ChimerasGenesList)
	#	if len(pi)!=len(pi2):
	#		DfTemp=DfTemp[~DfTemp[RefSp].isin(pi2)] #remove Paralog about Chimera (v8)
	print("DfTemp shape",DfTemp.shape)
	#DfMerge=pandas.concat([DfTemp,DfTemp2])# final merge cluster
	#DfTemp.to_csv("merge.txt",sep="\t")
	DfTemp=DfTemp.copy()
	for sp in OutGroupSpList:
		DfTemp=DfTemp[DfTemp[sp].isna()]
	DfNewGene=DfTemp
	DfNewGene.to_csv("newgene.txt",sep="\t")
	#Pl=[]
	Dl=[]
	Dl+=ParalogGeneList
	#for pi in ParalogGeneList:
		#a=DfNewGene[DfNewGene[RefSp].isin(pi)]
		#if not a.empty:
			#a=a.apply(lambda x: np.NaN if len(list(set(x)))==1 and pandas.isna(list(set(x))[0]) else ",".join(x.dropna()))
			#Pl.append(a.to_frame().T)
			#Pl.append(a)
		#Dl+=pi
	DfNewGeneParalog=DfNewGene[DfNewGene[RefSp].isin(ParalogGeneList)].copy()
    #DfNewGeneParalog=pandas.concat(Pl)
	Rl=[]
	for ri in RepeatGeneList:
		a=DfNewGene[DfNewGene[RefSp].isin(ri)]
		if not a.empty:
			Rl.append(a)
		Dl+=ri
	DfNewGeneRepeat=pandas.concat(Rl)
	DfNewGeneOrphan=DfNewGene[~DfNewGene[RefSp].isin(Dl)]
	#DfSpeciOnly.to_csv(Prefix+"_SpeciesSpecificGenesOnly.xls",sep="\t",index_label="Group")
	#DfSpeciRepeat.to_csv(Prefix+"_SpeciesSpecificGenesRepeat.xls",sep="\t",index_label="Group")
	DfNewGeneParalog.to_csv(Prefix+"_NewGenesParalog.xls",sep="\t",index_label="Group")
	DfNewGeneRepeat.to_csv(Prefix+"_NewGenesRepeat.xls",sep="\t",index_label="Group")
	DfNewGeneOrphan.to_csv(Prefix+"_NewGenesOrphan.xls",sep="\t",index_label="Group")
def PrepareData(Target,Prefix,Args):
	Fd={}
	pl=[]
	tl=[]
	OutGroupSpList=[]
	for x in open(Target): #target:sp\tgenome.fasta\tprotein.faa\tgene.bed\tgene_to_protein.txt/\t(in/out)
		x=x.rstrip()
		l=x.split("\t")
		sp=l[0]
		Fd[sp]={}
		Fd[sp]["genome"]=l[1]
		Fd[sp]["protein"]=l[2]
		Fd[sp]["gene"]=l[3]
		pl.append((sp,l[2]))
		tl.append((sp,l[4]))
		if l[-1]=="out" and not Args.outgroupfile: #outgroup
			OutGroupSpList.append(sp)
	if Args.outgroupfile:
		OutGroupSpList=[sp.rstrip() for sp in open(Args.outgroupfile)]
	TransToGened={}
	GeneToTransd=collections.defaultdict(list)
	AllGened={}
	ProteinSeqd={}
	for sp,f in tl:
		AllGened[sp]=[]
		for x in open(f):
			x=x.rstrip()
			l=x.split("\t")
			if len(l)<2:
				continue
			gene=l[0]
			protein=l[1]
			if not Args.genetotranspickle:
				GeneToTransd[gene].append(protein)
			if not Args.transtogenepickle:
				TransToGened[protein]=gene
			if not Args.allgenepickle:
				AllGened[sp].append(gene)
	if not Args.proteinpickle:
		for sp,f in pl:
			if f.endswith(".gz"):
				faa=gzip.open(f,"rt")
			else:
				faa=open(f,"r")
			d={}
			for record in SeqIO.parse(faa,"fasta"):
				d[record.id]=str(record.seq)
			ProteinSeqd[sp]=d
			faa.close()
		with open("%s_ProteinSeqd.pickle"%Prefix,"wb") as f:
			pickle.dump(ProteinSeqd,f)
	else:
		with open("%s_ProteinSeqd.pickle"%Prefix,"rb") as f:
			ProteinSeqd=pickle.load(f)
	if not Args.blastout:
		pfl=[]
		gzpfl=[]
		for sp,f in pl:
			if f.endswith(".gz"):
				gzpfl.append(f)
			else:
				pfl.append(f)
		if os.path.exists("%s_protein.faa"%Prefix):
			os.remove("%s_protein.faa"%Prefix)
		if pfl:
			os.system("cat %s >>%s_protein.faa"%(" ".join(pfl),Prefix))
		if gzpfl:
			os.system("zcat %s >>%s_protein.faa"%(" ".join(gzpfl),Prefix))#All protein names need to be different
		os.system("makeblastdb -in %s_protein.faa -out %s_protein -dbtype prot -title %s_protein"%(Prefix,Prefix,Prefix))
		os.system("blastp -seg yes -db %s_protein -query %s_protein.faa -out %s_protein_blast.out -num_threads 50 -evalue 0.001 -outfmt \"6 qseqid qlen qstart qend sseqid slen sstart send qcovs bitscore evalue pident nident\""%(Prefix,Prefix,Prefix))
		BlastOut="%s_protein_blast.out"%Prefix
	else:
		BlastOut=Args.blastout
	if not Args.blastpickle:
		Blastd=GetBlastd("%s_protein_blast.out"%Prefix,TransToGened)
		with open("%s_Blastd.pickle"%Prefix,"wb") as f:
			pickle.dump(Blastd,f)
	else:
		with open(Args.blastpickle,"rb") as f:
			Blastd=pickle.load(f)
	if not Args.genetotranspickle:
		with open("%s_GeneToTransd.pickle"%Prefix,"wb") as f:
			pickle.dump(GeneToTransd,f)
	else:
		with open(Args.genetotranspickle,"rb") as f:
			GeneToTransd=pickle.load(f)
	if not Args.transtogenepickle:
		with open("%s_TransToGened.pickle"%Prefix,"wb") as f:
			pickle.dump(TransToGened,f)
	else:
		with open(Args.transtogenepickle,"rb") as f:
			TransToGened=pickle.load(f)
	if not Args.allgenepickle:
		with open("%s_AllGened.pickle"%Prefix,"wb") as f:
			pickle.dump(AllGened,f)
	else:
		with open(Args.allgenepickle,"rb") as f:
			AllGened=pickle.load(f)
	return Fd,TransToGened,GeneToTransd,AllGened,ProteinSeqd,Blastd,BlastOut,OutGroupSpList
def GetBlastd(OriginalBlastOut,TransToGened):
	length=int(os.popen("wc -l %s"%OriginalBlastOut).readlines()[0].split()[0])
	d=collections.defaultdict(list)
	for x in tqdm(open(OriginalBlastOut),total=length):
		x=x.rstrip()
		l=x.split("\t")
		d[(TransToGened[l[0]],TransToGened[l[4]])].append((float(l[-4]),float(l[-2]),float(l[-1]),float(l[1]),l[0],l[4])) #bitscore,pident,nident,qlen,qseqid,sseqid
	Blastd={}
	for n in d:
		Blastd[n]=sorted(d[n])[-1]
	return Blastd
def CheckTarget(Target,Args):
	if not os.path.exists(Target):
		print("Target file does not exist:", Target)
		sys.exit(1)		
	spl=[]
	pl=[]
	gl=[]
	tl=[]
	for x in open(Target):
		x=x.rstrip()
		l=x.split("\t")
		if len(l)!=6:
			print("Target file format error, each line should contain 6 fields: species, genome, protein, gene bed, gene_to_protein, out/in group")
			sys.exit(1)
		spl.append(l[0])
		if not os.path.exists(l[1]):
			print("Genome fasta file does not exist:", l[1])
			sys.exit(1)
		if not os.path.exists(l[2]):
			print("Protein fasta file does not exist:", l[2])
			sys.exit(1)
		if not os.path.exists(l[3]):
			print("Gene bed file does not exist:", l[3])
			sys.exit(1)
		if not os.path.exists(l[4]):
			print("Gene to protein mapping file does not exist:", l[4])
			sys.exit(1)
		if l[-1] not in ["in","out"]:
			print("Last field should be 'in' or 'out' for outgroup specification.")
			sys.exit(1)
		tl.append(l[-1])
		for record in SeqIO.parse(l[2],"fasta"):
			pl.append(record.id)
		gl+=[xi.rstrip().split("\t")[3] for xi in open(l[3])]
	if len(spl)!=len(set(spl)):
		print("Species names in the target file should be unique.")
		sys.exit(1)
	if len(pl)!=len(set(pl)):
		print("Protein IDs in the target file should be unique.")
		sys.exit(1)
	if len(gl)!=len(set(gl)):
		print("Gene IDs in the target file should be unique.")
		sys.exit(1)
	if not "in" in tl or not "out" in tl:
		print("Target file should specify at least one species as 'in' and one as 'out' for outgroup.")
		sys.exit(1)
	if Args.refspecies not in spl:
		print("Reference species (%s) is not in the target file."%Args.refspecies)
		sys.exit(1)
	if Args.outgroupfile:
		for sp in open(Args.outgroupfile,"r"):
			sp=sp.rstrip()
			if sp not in spl:
				print("Outgroup species %s is not in the target file."%sp)
				sys.exit(1)
def GetParser():
	parser = argparse.ArgumentParser(prog='SNGFinder',description='Identifying new genes based on the syntenic method.')
	parser.add_argument('--prefix', required=True,help='Result file prefix.')
	parser.add_argument('--maf', required=True,help="Whole-genome multiple sequence alignment file in maf format.")
	parser.add_argument("--target",required=True,help="A target file.")
	parser.add_argument("--blastpickle",required=False,help="The pickle file that stores the best alignment results of blast. if this file is provided, skip the GetBlastd function.")
	parser.add_argument("--proteinpickle",required=False,help="The pickle file that stores the protein sequence. if this file is provided, it is not necessary to read the protein sequence from the FASTA file.")
	parser.add_argument("--genetotranspickle",required=False,help="The pickle file that stores the correspondence between genes and transcripts.")
	parser.add_argument("--transtogenepickle",required=False,help="The pickle file that stores the corresponding between transcripts and genes.")
	parser.add_argument("--allgenepickle",required=False,help="The pickle file that stores all gene names for each species.")
	parser.add_argument("--blastout",required=False,help="all vs. all blast results.")
	parser.add_argument("--syntenicpickle",required=False,help="The pickle file that stores the syntenic information between species. If this file is provided, skip the GetSyntenic function.")
	parser.add_argument("--blockpickle",required=False,help="The pickle file that stores the syntenic gene block information between species. If this file is provided, skip the GetSyntenicGene function.")
	parser.add_argument("--thpickle",required=False,help="The pickle file that stores the threshold for judging ortholog or paralog. If this file is provided, skip the GetRBH function.")
	parser.add_argument("--chimericlistfile",required=False,help="The file that stores the list of chimeric genes. If this file is provided, skip the GetChimeras function.")
	parser.add_argument("--refspecies",required=True,help="Reference species name, which should be consistent with the species name in the target file.")
	parser.add_argument("--outgroupfile",required=False,help="The file that stores the outgroup species, one species per line. If this file is provided, the outgroup species will be determined based on this file instead of the target file.")
	parser.add_argument("--dnaidentify",type=float,default=0.5,help="The threshold for judging ortholog based on DNA sequence identity. If the identity of the best blast hit in the syntenic region is greater than or equal to this threshold, it will be considered as an ortholog.")
	parser.add_argument("--inflation",type=float,default=3,help="The inflation parameter for MCL clustering of paralog genes. The larger the value, the finer the clustering.")
	parser.add_argument("--overlapth",type=float,default=0.05,help="Maximum overlap between blast hits allowed when judging chimeric genes.")
	parser.add_argument("--leftpercentile",type=float,default=20,help="The left percentile of the blast hits bit score distribution used to judge chimeric genes. If the blast hits bit score is less than this percentile, it will be removed.")
	parser.add_argument("--rightpercentile",type=float,default=100,help="The right percentile of the blast hits bit score distribution used to judge chimeric genes. If the blast hits bit score is greater than this percentile, it will be removed.")
	parser.add_argument("--nochimeras",action='store_true',help="If this option is set, the chimeric genes will not be identified.")
	parser.add_argument("--nounannotated",action='store_true',help="If this option is set, the unannotated genes will not be identified.")
	return parser
if __name__ == '__main__':
	Parser=GetParser()
	Args = Parser.parse_args()
	Prefix=Args.prefix
	MafFile=Args.maf
	Target=Args.target
	RefSp=Args.refspecies
	GeneDnaPairsIdentify=Args.dnaidentify #0.5
	NoChimeras=Args.nochimeras
	NoUnannotated=Args.nounannotated
	Inflation=Args.inflation #3
	OverlapTh=Args.overlapth
	LeftPercentile=Args.leftpercentile
	RightPercentile=Args.rightpercentile
	#OutGroupSpList=[x.rstrip() for x in open(OutGroupFile)]
	#print("GeneDnaPairsIdentify type",type(GeneDnaPairsIdentify))
	CheckTarget(Target, Args)
	Fd,TransToGened,GeneToTransd,AllGened,ProteinSeqd,Blastd,BlastOut,OutGroupSpList=PrepareData(Target,Prefix,Args)
	if not Args.syntenicpickle:
		Syntenicd=GetSyntenic(MafFile,Prefix)
		with open('%s_Syntenicd.pickle'%Prefix,'wb') as f:
			pickle.dump(Syntenicd, f)
	else:
		with open(Args.syntenicpickle,"rb") as f:
			Syntenicd=pickle.load(f)
	if not Args.blockpickle:
		Blockd=GetSyntenicGene(Syntenicd,Fd,Prefix)
	else:
		with open(Args.blockpickle,"rb") as f:
			Blockd=pickle.load(f)
	SupplementSrc.AllGened=AllGened
	SupplementSrc.Blastd=Blastd
	if not Args.thpickle:
		GenePairsScored=SupplementSrc.GetRBH(RefSp)
		with open("%s_BitScore_threshold.pickle"%Prefix,"wb") as f:
			pickle.dump(GenePairsScored,f)
	else:
		with open(Args.thpickle,"rb") as f:
			GenePairsScored=pickle.load(f)
	OrthoCluster,OrthoDNA,TransProtein=GetOrthoCluster(Syntenicd,Blockd,AllGened,GenePairsScored,Blastd,ProteinSeqd,GeneToTransd,RefSp,Fd,Prefix,NoUnannotated)
	ChimerasGenesList=[]
	if not NoChimeras:
		if not Args.chimericlistfile:
			ChimerasGened,ChimerasTransd,OriginalChimerasd=SupplementSrc.GetChimeras(BlastOut,RefSp,TransToGened,OverlapTh,LeftPercentile,RightPercentile)#change fixed parameters
			pandas.DataFrame.from_dict(ChimerasTransd,orient="index").to_csv("%s_ChimerasTrans.xls"%Prefix,sep="\t",index_label=RefSp)
			with open("%s_ChimerasTransd.pickle"%Prefix,"wb") as f:
				pickle.dump(ChimerasTransd,f)
			with open("%s_OriginalChimerasd.pickle"%Prefix,"wb") as f:
				pickle.dump(OriginalChimerasd,f)
			Fr=open("%s_ChimerasGenesList.txt"%Prefix,"w") #Check
			for gl in ChimerasGened:
				ChimerasGenesList+=gl.split(",")
				Fr.write(gl+"\n")
			Fr.close()
		else:
			ChimerasGenesList=[]
			for gl in open(Args.chimericlistfile):
				gl=gl.rstrip()
				ChimerasGenesList+=gl.split(",")
	Df=FilterOrthoCluster(GenePairsScored,GeneDnaPairsIdentify,OrthoCluster,OrthoDNA,TransProtein,RefSp,Prefix,NoUnannotated)
	RepeatGeneList,ParalogGeneList=GetRepeatAndParalogGene(Df.copy(),RefSp,Blastd,GenePairsScored,ProteinSeqd,GeneToTransd,AllGened,TransToGened,Inflation,Prefix)
	GetGeneClassify(Df.copy(),ChimerasGenesList,RepeatGeneList,ParalogGeneList,OutGroupSpList,RefSp,Prefix)
