allf=open('filelist','r').read().split('\n')[:-1]

for filename in allf:
	a=open(filename,'r').read().split('\n')[:-1]
	a=[c for c in a if 'SumGF' in c]

	b=[]
	for c in a:
		c=c.split(' ')
		b+=[c[0].split(':')[0].split('\t')[1]+'\t'+c[0].split(':')[1]+'\t'+c[1]+'\t'+c[2].split(':')[0]+'\t'+c[2].split(':')[1]+'\t'+c[3].split(':')[0]+'\t'+c[3].split(':')[1]]

	f=open('genefusion_'+filename,'w')
	for c in b:
		f.write(c+'\n')
	f.close()
