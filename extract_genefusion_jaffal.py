samplelist=open('filelist','r').read().split('\n')[:-1]

for samplename in samplelist:
	a=open(samplename+'/jaffa_results.csv','r').read().split('\n')[1:-1]
	# keep only HighConfidence and LowConfidence events, filter PotentialTransSplicing
	a=[c for c in a if 'Confidence' in c]
	genefusion=[]
	for c in a:
		c=c.split('"')
		genefusion+=[c[3].split(':')[0]+'\t'+c[3].split(':')[1]+'\t'+c[5]+'\t'+c[6].split(',')[1]+'\t'+c[9]+'\t'+c[10].split(',')[1]+'\t'+c[12].split(',')[3]+'\t'+c[15]+'\t'+c[13]]

	f=open(samplename+'/jaffa_gene_fusion','w')
	for c in genefusion:
		f.write(c+'\n')
	f.close()


