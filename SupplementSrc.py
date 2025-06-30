#!/usr/bin/env python
#coding=utf-8
#20%过滤，同一转录本多个片段求和后计算 (v2)
#有5% overlap，按照转录本计算 (v3又重新进行了优化)
import collections,pickle,sys,os
from tqdm import tqdm
from multiprocessing import Pool as ProcessPool
import pandas
import numpy as np
global AllGened,Blastd
def GetRightSystem(Cmd):
	ExitStatus=1
	Attempts = 0
	while ExitStatus!=0:
		ExitStatus=os.system(Cmd)
		Attempts+=1
		if Attempts==3: #3 need  change
			print("GetRightSystem break: %s"%Cmd)
			break
	return ExitStatus
def GetRefd(RefSp,RefGene):
	#print(RefGene)
	Refid={}
	for sp in AllGened:
		if sp==RefSp:
			continue
		sl=[]
		for TargetGene in AllGened[sp]:
			try:
				Bitscore,Pident,Nident,Qlen,_,_=Blastd[(RefGene,TargetGene)]
				sl.append((float(Bitscore),TargetGene))
			except KeyError:
				pass
		if sl:
			Refid[sp]=sorted(sl)[-1]
	return RefGene,Refid
def GetTard(RefSp,TargetGene):
	#print(TargetGene)
	sl=[]
	for RefGene in AllGened[RefSp]:
		try:
			Bitscore,Pident,Nident,Qlen,_,_=Blastd[(TargetGene,RefGene)]
			sl.append((float(Bitscore),RefGene))
		except KeyError:
			pass
	if sl:
		return TargetGene,sorted(sl)[-1]
	else:
		return TargetGene,None
def GetRBH(RefSp):
	Refd={}
	TardTemp={}
	Params=[(RefSp,RefGene) for RefGene in AllGened[RefSp]]
	pool = ProcessPool(50)
	R=pool.starmap(GetRefd,Params)
	pool.close()
	pool.join()
	for ri in R:
		RefGene,Refid=ri
		Refd[RefGene]={}
		for sp in Refid:
			Refd[RefGene][sp]=Refid[sp]
	Tard={}
	for sp in AllGened:
		if sp==RefSp:
			continue
		Tard[sp]={}
		Params=[(RefSp,TargetGene) for TargetGene in AllGened[sp]]
		pool=ProcessPool(50)
		R=pool.starmap(GetTard,Params)
		pool.close()
		pool.join()
		for ri in R:
			TargetGene,Item=ri
			Tard[sp][TargetGene]=Item
	Thd={}
	Lastl=[]
	for RefGene in Refd:
		sl=[]
		for sp in Refd[RefGene]:
			Bitscore,TargetGene=Refd[RefGene][sp]
			if Tard[sp][TargetGene]:
				BitscoreTemp,RefGeneTemp=Tard[sp][TargetGene]
				if RefGeneTemp==RefGene:
					sl.append((Bitscore,sp,TargetGene))
		print(RefGene,sorted(sl))
		if sl:
			Thd[RefGene]=sorted(sl)[0][0]
		else:
			for sp in Refd[RefGene]:
				Bitscore,TargetGene=Refd[RefGene][sp]
				sl.append((Bitscore,sp,TargetGene))
			print("NORBH",RefGene,sorted(sl))
			if sl:
				Thd[RefGene]=sorted(sl)[0][0]
			else:
				Thd[RefGene]=None
				Lastl.append(RefGene)
	sd={}
	for g in Lastl:
		sd[g]=[]
		for RefGene in AllGened[RefSp]:
			try:
				Bitscore,Pident,Nident,Qlen,_,_=Blastd[(g,RefGene)]
				sd[g].append((float(Bitscore),RefGene))
			except KeyError:
				try:
					Bitscore,Pident,Nident,Qlen,_,_=Blastd[(RefGene,g)]
					sd[g].append((float(Bitscore),RefGene))
				except KeyError:
					continue
	for g in sd:
		if sd[g]:
			Thd[g]=sorted(sd[g])[0][0]
	return Thd
def GetOverlap(a, b):
	return max(0, min(a[1], b[1]) - max(a[0], b[0]))
def GetChimeras(BlastOut,RefSp,TransToGened,OverlapTh,LeftPercentile,RightPercentile):
	d={}
	rd={}
	#ChiList=[]
	for sp in AllGened:
		d[sp]=set(AllGened[sp])
	for sp in d:
		if sp==RefSp:
			continue
		rd[(RefSp,sp)]=collections.defaultdict(list)
		rd[(sp,RefSp)]=collections.defaultdict(list)
	#Fd,Lengthd=FilterBlast(BlastOut,TransToGened,d[RefSp]) #FilterBlast
	Fd,Lengthd=FilterBlast2(BlastOut,TransToGened,d,RefSp,LeftPercentile,RightPercentile) #Change,for bitscore
	for pair in Fd:
		spair=()
		if TransToGened[pair[0]] in d[RefSp]:
			for sp in d:
				if sp==RefSp:
					continue
				if TransToGened[pair[1]] in d[sp]:
					spair=(RefSp,sp)
					break
		if TransToGened[pair[1]] in d[RefSp]:
			for sp in d:
				if sp==RefSp:
					continue
				if TransToGened[pair[0]] in d[sp]:
					spair=(sp,RefSp)
					break
		if spair:
			for ri in Fd[pair]:
				rd[spair][pair[0]].append((ri[0],ri[1],ri[2],pair[1]))
	FinalR=collections.defaultdict(list)
	for sp1,sp2 in rd:
		#print(sp1,sp2)
		for gene in rd[(sp1,sp2)]:
			#print(gene)
			pl=rd[(sp1,sp2)][gene] ##start,end,score,id
			flag=0
			Tempalignd=collections.defaultdict(list)
			for pi in pl:
				Tempalignd[pi[-1]].append(pi)
			Alignl=[]
			for ti in Tempalignd:
				score=0
				tl=[]
				for pi in sorted(Tempalignd[ti],key=lambda x: x[2],reverse=True):
					if tl==[]:
						score+=pi[2]
						tl.append(pi)
					else:
						for tpi in tl:
							if GetOverlap(pi,tpi)!=0:
								break
						else:
							tl.append(pi)
							score+=pi[2]
				Alignl.append((score,ti,Tempalignd[ti]))
			while flag==0 and Alignl:
				Alignl = sorted(Alignl,reverse=True)
				#print("Alignl",Alignl) 
				cl=GetInterval2(Alignl,Lengthd[gene],OverlapTh,TransToGened) #Change,for transcripts
				#cl=GetInterval(pl,Lengthd[gene],OverlapTh) #overlap hit<=5%
				#print("GetInterval",cl)
				cl2=MergeInterval(cl) #final no overlap interval
				#print("MergeInterval",cl2)
				cov=0
				gi=[]
				gtd=collections.defaultdict(list)
				for ci in cl2:
					cov+=abs(ci[1]-ci[0])+1
					for gii in set(ci[2]):
						gi.append(TransToGened[gii])
						gtd[TransToGened[gii]].append(gii)
				tflag=0
				for gti in gtd:
					if len(set(gtd[gti]))>1:
						print("two or more same transcripts",cl2)
						tflag=1
				cov=cov/Lengthd[gene]
				print("transcript,coverge,support,tflag",(gene,cov,set(gi),tflag))
				if cov>=0.8 and len(set(gi))>=2 and tflag!=1: #0.8 need change
					cl=FurtherFilterCl(cl,Lengthd[gene],TransToGened) #change,进一步去除多余片段
					cl2=MergeInterval(cl) #change
					print("chimeras",cl2)
					FinalR[gene].append((cl2,(sp1,sp2)))
					flag=1
				else:
					Alignl=Alignl[1:]
			#if cl:
			#	sys.exit()
	FinalRToGene={} #同一基因不同转录本的结果会被合并到一起
	FinalRToTrans={}
	for g in FinalR:
		if TransToGened[g] in d[RefSp]:
			if not TransToGened[g] in FinalRToGene:
				FinalRToGene[TransToGened[g]]={}
			if not g in FinalRToTrans:
				FinalRToTrans[g]={}
			ifRef=1
		else:
			ifRef=0
		for fi in FinalR[g]:
			if fi[1][0]==RefSp:
				sp=fi[1][1]
			elif fi[1][1]==RefSp:
				sp=fi[1][0]
			if ifRef==0:
				#RefGene=""
				#print("fi[0]",fi[0])
				RefGeneTempList=[]
				RefTransTempList=[]
				for gi in fi[0]:
					#RefGene+=",".join(set([TransToGened[gii] for gii in set(gi[-1])]))+","
					for gii in set(gi[-1]):
						RefGeneTempList.append(TransToGened[gii])
						RefTransTempList.append(gii)
				RefGene=",".join(sorted(set(RefGeneTempList)))
				RefTrans=",".join(sorted(set(RefTransTempList)))
				if not RefGene in FinalRToGene:
					FinalRToGene[RefGene]={}
				if not RefTrans in FinalRToTrans:
					FinalRToTrans[RefTrans]={}
				FinalRToGene[RefGene][sp]=TransToGened[g]
				FinalRToTrans[RefTrans][sp]=g
			else:
				#ti=""
				tiGeneList=[]
				tiTransList=[]
				#print("fi[0]",fi[0])
				for gi in fi[0]:
					#ti+=",".join([TransToGened[gii] for gii in set(gi[-1])])+","
					for gii in set(gi[-1]):
						tiGeneList.append(TransToGened[gii])
						tiTransList.append(gii)
				tiGene=",".join(sorted(set(tiGeneList)))
				tiTrans=",".join(sorted(set(tiTransList)))
				#print(g,sp,ti)
				FinalRToGene[TransToGened[g]][sp]=tiGene
				FinalRToTrans[g][sp]=tiTrans
	
	return FinalRToGene,FinalRToTrans,FinalR
def FurtherFilterCl(cl,TransLength,TransToGened): #start,end,score,id
	cl=sorted(cl,key=lambda ll:ll[2],reverse=True)
	flag=0
	n=0
	while flag==0:
		gi=[]
		cov=0
		Cltemp=MergeInterval(cl[:2+n])
		for ci in Cltemp:
			cov+=abs(ci[1]-ci[0])+1
			for gii in set(ci[2]):
				gi.append(TransToGened[gii])
		cov=cov/TransLength
		if cov>=0.8 and len(set(gi))>=2: #need change
			flag=1
		else:
			n+=1
	return cl[:2+n]
def FilterBlast(BlastOut,TransToGened,RefGeneSet):
	length=int(os.popen("wc -l %s"%BlastOut).readlines()[0].split()[0])
	Lengthd={}
	d=collections.defaultdict(list)
	for x in tqdm(open(BlastOut),total=length):
		x=x.rstrip()
		l=x.split("\t")
		Lengthd[l[0]]=float(l[1])
		if float(l[8])>=20 and float(l[8])<80 and (TransToGened[l[0]] in RefGeneSet or TransToGened[l[4]] in RefGeneSet): #l[8] qcov,80 need change
			d[(l[0],l[4])].append((int(l[2]),int(l[3]),float(l[-4])))
	return d,Lengthd
def FilterBlast2(BlastOut,TransToGened,SpGeneSet,RefSp,LeftPercentile,RightPercentile):
	length=int(os.popen("wc -l %s"%BlastOut).readlines()[0].split()[0])
	Lengthd={}
	d={}
	for sp in SpGeneSet:
		if sp==RefSp:
			continue
		d[(RefSp,sp)]=collections.defaultdict(list)
	for x in tqdm(open(BlastOut),total=length):
		x=x.rstrip()
		l=x.split("\t")
		Lengthd[l[0]]=float(l[1])
		if (TransToGened[l[0]] in SpGeneSet[RefSp]) and (not TransToGened[l[4]] in SpGeneSet[RefSp]):
			TarSp=""
			for sp in SpGeneSet:
				if TransToGened[l[4]] in SpGeneSet[sp] and sp!=RefSp:
					TarSp=sp
					break
			if TarSp=="":
				print("SpGeneSet error",l[4])
				sys.exit()
			d[(RefSp,TarSp)][(l[0],l[4])].append((int(l[2]),int(l[3]),float(l[-4])))
		elif (TransToGened[l[4]] in SpGeneSet[RefSp]) and (not TransToGened[l[0]] in SpGeneSet[RefSp]):
			TarSp=""
			for sp in SpGeneSet:
				if TransToGened[l[0]] in SpGeneSet[sp] and sp!=RefSp:
					TarSp=sp
					break
			if TarSp=="":
				print("SpGeneSet error",l[0])
				sys.exit()
			d[(RefSp,TarSp)][(l[0],l[4])].append((int(l[2]),int(l[3]),float(l[-4])))
	Fd={}
	for SpPair in d:
		sl=[]
		for gPair in d[SpPair]:
			for b in d[SpPair][gPair]:
				sl.append(b[-1])
		mi=np.percentile(sl,LeftPercentile)  #5,change 
		ma=np.percentile(sl,RightPercentile) #95,change 
		#print(sl)
		print(SpPair)
		print("bitscore threshold",mi,ma)
		for gPair in d[SpPair]:
			for b in d[SpPair][gPair]:
				if b[-1]>=mi and b[-1]<=ma: #change
					Fd[gPair]=d[SpPair][gPair]
					break
	return Fd,Lengthd
def GetInterval2(Alignl,TransLength,OverlapTh,TransToGened):# (score,ti,pl)
	rl=[]
	nl=[]
	for n,ali in enumerate(Alignl):
		pl=ali[2]
		if n==0:
			rl+=pl
			nl.append(TransToGened[ali[1]])
			continue
		if TransToGened[ali[1]] in nl:
			continue
		flag=0
		for ri in rl:
			for pi in pl:
				Olen=GetOverlap(ri,pi)
				if Olen/TransLength>OverlapTh or Olen/(abs(pi[1]-pi[0])+1)>0.9 or Olen/(abs(ri[1]-ri[0])+1)>0.9: #0.05,change,小片段包含在另一个片段内的去除
					flag=1
					break
			if flag==1:
				break
		if flag!=1:
			rl+=pl
			nl.append(TransToGened[ali[1]])
	return rl
def GetInterval(Pl,GeneLength,OverlapTh): #OverlapTh=0.05
	#Pl = sorted(Pl, key=lambda x: x[2],reverse=True)
	print("sorted",Pl)
	rl=[]
	ol=[]
	for n,pi in enumerate(Pl):
		if n==0:
			rl.append(pi)
			continue
		for ri in rl:
			if ri[-1]!=pi[-1]:
				Olen=GetOverlap(ri,pi)
				#print("overlap",Olen,ri,pi)
				if Olen/GeneLength>OverlapTh: #need change,0.05
					ol.append(pi[-1])
					break
		else:
			rl.append(pi)
	#url=set(url)
	#print("url",url)
	#rl=[]
	#rl=list(set(Pl)-url)
	Frl=[]
	for pi in rl:
		if not pi[-1] in ol:
			Frl.append(pi)
	return Frl
def MergeInterval(intervals):
	if len(intervals) == 1:
		return [[intervals[0][0],intervals[0][1],[intervals[0][3]]]]
	if len(intervals) == 0:
		return intervals
	intervals.sort(key=lambda x:x[0])
	result = [[intervals[0][0],intervals[0][1],[intervals[0][3]]]]
	for interval in intervals[1:]:
		if interval[0] <= result[-1][1]:
			result[-1][1] = max(result[-1][1], interval[1])
			result[-1][2].append(interval[3])
		else:
			result.append([interval[0],interval[1],[interval[3]]])
	return result
if __name__ == '__main__':
	RefSp="Oryza_japonica"
	Prefix="OryzaSp11_ref_Oryza_japonica"
	#global Blastd,AllGened
	with open("OryzaSp11_Blastd.pickle","rb") as f:
		Blastd=pickle.load(f)
	with open("OryzaSp11_AllGened.pickle","rb") as f:
		AllGened=pickle.load(f)
	#with open("OryzaSp11_GeneToTransd.pickle","rb") as f:
	#	GeneToTransd=pickle.load(f)
	Thd=GetRBH(RefSp)
	with open("%s_BitScore_threshold.pickle"%Prefix,"wb") as f:
		pickle.dump(Thd,f)
	#print(Thd)
	#with open("OryzaSp11_TransToGened.pickle","rb") as f:
	#	TransToGened=pickle.load(f)
	#BlastOut="all_protein_out.txt"
	#FinalRToGene,FinalR=GetChimeras(BlastOut,RefSp,TransToGened)
	#df=pandas.DataFrame.from_dict(FinalRToGene,orient="index")
	#df.to_csv("test_chimeras.xls",sep="\t",index_label=RefSp)
	#print(Chimerasd)
