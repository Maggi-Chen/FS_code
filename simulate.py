allt=open('../goodgene_gencode.v39.transcripts.fa','r').read().split('>')[1:]

protein=[c for c in allt if 'protein_coding' in c and 'ENSG0' not in c.split('|')[5] ]


f1=open('highexpress.fa','w')
f2=open('mediumexpress.fa','w')
f3=open('lowexpress.fa','w')


import random 

# assign normal transcripts to low/medium/high expression groups
for c in protein:
	name=c.split('\n')[0]
	seq=''.join(c.split('\n')[1:])
	if len(seq)<100:
		continue
	randi=random.randint(1,3)
	if randi==1:
		f1.write('>'+name+'\n'+seq+'\n')
	if randi==2:
		f2.write('>'+name+'\n'+seq+'\n')
	if randi==3:
		f3.write('>'+name+'\n'+seq+'\n')

done=0

fused_transcript=[]


# pick 300 random genes & pair as gene fusions
goodfuse=[c for c in protein if 2000<=len(''.join(c.split('\n')[1:]))<=50000]
geneid=list(set([c.split('|')[1].split('.')[0] for c in goodfuse]))
print len(geneid)
nonoverlap=open('../nonoverlapgene','r').read().split('\n')[:-1]
geneid=[c for c in geneid if c in nonoverlap]
print len(geneid)

fuseid=random.sample(geneid,400)
random.shuffle(fuseid)

f=open('fused_transcript_all.fa','w')


g1=open('fused_transcript_high.fa','w')
g2=open('fused_transcript_medium.fa','w')
g3=open('fused_transcript_low.fa','w')


# simulate 100 gene fusions with breakpoint in exon
for i in range(100):	
	gene1=fuseid[2*i]
	gene2=fuseid[2*i+1]

	trans1=[c for c in goodfuse if c.split('|')[1].split('.')[0] ==gene1]
	trans1=trans1[random.randint(0,len(trans1)-1)]
	seq1=''.join(trans1.split('\n')[1:])
	trans2=[c for c in goodfuse if c.split('|')[1].split('.')[0] ==gene2]
	trans2=trans2[random.randint(0,len(trans2)-1)]
	seq2=''.join(trans2.split('\n')[1:])

	bppos1=random.randint(1000,len(seq1)-1000)
	bppos2=random.randint(1000,len(seq2)-1000)
	leftright=random.randint(1,2)
	if leftright==1:
		leriinfo='left,right'
		newseq=seq1[:bppos1]+seq2[bppos2:]
	if leftright==2:
		leriinfo='left,left'
		newseq=seq1[:bppos1]+seq2[:bppos2]
	if leftright==3:
		leriinfo='right,right'
		newseq=seq1[bppos1:]+seq2[bppos2:]
	if leftright==4:
		leriinfo='right,left'
		newseq=seq1[bppos1:]+seq2[:bppos2]
	newname='fusedgene_'+str(done)+'_'+leriinfo+','+str(bppos1)+','+str(bppos2)+'...'+trans1.split('\n')[0]+'...'+trans2.split('\n')[0]
	f.write('>'+newname+'\n'+newseq+'\n')

	highlow=random.randint(1,3)
	if highlow==1:
		f1.write('>'+newname+'\n'+newseq+'\n')
		g1.write('>'+newname+'\n'+newseq+'\n')
	if highlow==2:
		f2.write('>'+newname+'\n'+newseq+'\n')
		g2.write('>'+newname+'\n'+newseq+'\n')
	if highlow==3:
		f3.write('>'+newname+'\n'+newseq+'\n')
		g3.write('>'+newname+'\n'+newseq+'\n')

refseq=open('/data/user/maggic/svstudy/data/reference/hg38.fa','r').read().split('>')[1:-1]
ref={}

dd={}
dd['A']='T'
dd['T']='A'
dd['G']='C'
dd['C']='G'

def getrev(seq):
	new=''
	for c in seq:
		new=dd[c]+new
	return new


for c in refseq:
	c=c.split('\n')[:-1]
	chrom=c[0].split(' ')[0]
	if '_' in chrom:
		continue
	seq=''.join(c[1:])
	ref[chrom]=seq
	if chrom=='chr1':
		print 'chr1 len' ,len(seq)

introngene=fuseid[200:]
print len(introngene)
gtf=open('/data/project/chonglab/Maggic/python/defusion/Homo_sapiens.GRCh38.104.chrname.gtf','r').read().split('\n')[5:]

gtfinfo=[c for c in gtf if 'gene_id "' in c and  c.split('gene_id "')[1].split('"')[0] in introngene]


doneintronfused=0


# simulate gene fusions with breakpoint in introns
for i in range(100,200):
	gene1=fuseid[2*i]
	gene2=fuseid[2*i+1]

	gtf1=[c for c in gtfinfo if  c.split('gene_id "')[1].split('"')[0] ==gene1]
	print gene1,len(gtf1)
	traid=[c for c in gtf1 if c.split('\t')[2]=='transcript'][0].split('transcript_id "')[1].split('"')[0]
	gtf1=[c for c in gtf1 if 'transcript_id "' in c and c.split('transcript_id "')[1].split('"')[0] == traid and c.split('\t')[2]=='exon']
	chromseq=ref[gtf1[0].split('\t')[0]]
	transcriptseq=''
	strand=gtf1[0].split('\t')[6]
	for c in gtf1:
		if  strand=='-':
			transcriptseq+=getrev(chromseq[int(c.split('\t')[3]):int(c.split('\t')[4])])
		else:
			transcriptseq+=chromseq[int(c.split('\t')[3]):int(c.split('\t')[4])]
	if len(transcriptseq)<500:
		print 'transcript length <500'
		continue
	bp1=random.randint(500,len(transcriptseq))

	fusedseq=transcriptseq[:bp1]


	gtf2=[c for c in gtfinfo if 'gene_id "'+gene2+'";' in c and c.split('\t')[2]!='transcript']
	gene2pos=[c for c in gtf2 if c.split('\t')[2]=='gene'][0]
	allexon=[c for c in gtf2 if c.split('\t')[2]!='gene']
	allexonpos=[]
	for c in allexon:
		allexonpos+=[[int(c.split('\t')[3]),int(c.split('\t')[4])]]

	testround=0
	goodout=0

	while testround<=5:
		try:
			bp1=random.randint(int(gene2pos.split('\t')[3])+500,int(gene2pos.split('\t')[4])-1000)
			bp2=random.randint(bp1+500,min(int(gene2pos.split('\t')[4])-500,bp1+5000))
			ifintron=0
			print bp1,bp2
			for d in allexonpos:
				if d[0]<=bp1<=d[1] or d[0]<=bp2<=d[1]:
					ifintron==1;break
			if ifintron==0:
				goodout=1
				break
			else:
				testround+=1
		except:
			testround+=1


	if goodout==0:
		print 'failed find intron bp '
		aa=input('key');continue

	intronicseq=ref[gene2pos.split('\t')[0]][bp1:bp2]
	if gene2pos.split('\t')[6]=='-':
		intronicseq=getrev(intronicseq)
	fusedseq=fusedseq+intronicseq
	newname=gene1+'_'+gene2+'_'+str(bp1)+'_'+strand+'_'+str(bp1)+'_'+str(bp2)+'_'+gene2pos.split('\t')[6]
	f.write('>'+newname+'\n'+fusedseq+'\n')
	highlow=random.randint(1,3)
	if highlow==1:
		f1.write('>'+newname+'\n'+fusedseq+'\n')
		g1.write('>'+newname+'\n'+fusedseq+'\n')
	if highlow==2:
		f2.write('>'+newname+'\n'+fusedseq+'\n')
		g2.write('>'+newname+'\n'+fusedseq+'\n')
	if highlow==3:
		f3.write('>'+newname+'\n'+fusedseq+'\n')
		g3.write('>'+newname+'\n'+fusedseq+'\n')
	doneintronfused+=1
	if doneintronfused>=50:
		break
		

f1.close()
f2.close()
f3.close()
g1.close()
g2.close()
g3.close()


f.close()

















