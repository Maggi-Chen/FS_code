import sys
datatype='pacbio'
# ground truth gene fusion list
vcf=open('genome/genefusion_breakpoint','r').read().split('\n')[:-1]
# fusionseeker gene fusion list
pac=open('simulation_rep1/confident_genefusion_polished.txt','r').read().split('\n')[:-1]


goodchrom=['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']

# select gene fusions in autosomes and sex chromosomes
new=[c for c in pac if c.split('\t')[4] in goodchrom and c.split('\t')[6] in goodchrom]

pac=[]
# select unique gene pairs
for c in new:
	ifdone=0
	for d in pac:
		if (c.split('\t')[1]==d.split('\t')[1] and c.split('\t')[2]==d.split('\t')[2] ) or (c.split('\t')[1]==d.split('\t')[2] and c.split('\t')[2]==d.split('\t')[1]):
			ifdone=1
			break
	if ifdone==0:
		pac+=[c]

detected=[]
correct=[]

for c in vcf:
	gene1=c.split('\t')[0]
	gene2=c.split('\t')[1]

	for d in pac:
		# two gene fusions are considered as match if two fused genes are identical 
		if (d.split('\t')[1]==gene1 and d.split('\t')[2]==gene2) or (d.split('\t')[1]==gene2 and d.split('\t')[2]==gene1):
			detected+=[c]
			correct+=[d]
			break

# print number of ground truth, true positive, reported gene fusions
print len(vcf),len(detected),len(pac)


# write results to txt files
f=open('detected_'+datatype,'w')
for c in detected:
	f.write(c+'\n')
f.close()
f=open(datatype+'_results','w')
false=[c for c in pac if c not in correct]
correct=list(set(correct))

for c in correct:
	f.write('TP\t'+c+'\n')
for c in false:
	f.write('FP\t'+c+'\n')
f.close()

missed=[c for c in vcf if c not in detected]
f=open('missed_'+datatype,'w')
for c in missed:
	f.write(c+'\n')
f.close()

