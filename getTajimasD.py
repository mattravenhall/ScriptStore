print('Importing libraries',flush=True)
import pandas as pd
import numpy as np
import multiprocessing

tpedFile = './'
tfamFile = './'
annotFile = './'

numThreads = 20     # Number of threads to use when multi-threading Tajima's D calculation - NB. limited to 52 due to current naming convention
minSNPs = 4

print('Reading in snp calls file',flush=True)
dat = pd.read_table(tpedFile,header=None,sep=' ')

print('Reading in sample IDs file',flush=True)
tfam = pd.read_table(tfamFile,header=None,sep=' ')
IDs = [x for x in pd.concat([tfam.ix[:,0]+'_A',tfam.ix[:,0]+'_B'],axis=1).stack().values] # interleave column 1_A and column 2_B

print('Reading in annotation file',flush=True)
annot = pd.read_table(annotFile,header=1,sep=' ')

dat.columns = ['Chromosome','snpID','-','Location']+IDs
dat['Gene'] = 'Intergenic'

print('Annotating SNPs',flush=True)
newGenes = []
for x in range(0,dat.shape[0]):
    y = (annot['Gene'].ix[(annot['Chromosome'] == str(dat['Chromosome'].ix[x])) & (annot['Start'] < dat['Location'].ix[x]) & (annot['End'] > dat['Location'].ix[x])].str.cat())
    if y == '':
        newGenes.append('Intergenic')
    else:
        newGenes.append(y)

dat['Gene'] = newGenes
dat['-'] = 'BORDER'

print('Removing intergenic regions',flush=True)
dat <- dat.ix[(dat['Gene'] != 'Intergenic') & (dat['Gene'] != '')]

# I think that 0=missing, 1=ref, 2=alt
# So, remove all SNPs with 0s.
print('Removing missing SNPs',flush=True)
dat = dat.ix[['Intergenic' not in dat.ix[x].values for x in range(0,dat.shape[0])]]

if dat.shape[0] == 0:
    sys.exit('Error: No complete case SNPs present')

# And make 1 to 0, and 2 to 1.
print('Properly formatting calls',flush=True)
tmp = dat['Gene']
dat = pd.concat([dat.ix[:,0:3],dat.ix[:,4:(dat.shape[1]-1)]],axis=1)
dat['Gene'] = tmp

# Split calls from IDs
print('Prepping for Tajimas D',flush=True)
calls = dat.ix[:,4:(dat.shape[1]-1)]
info = dat.ix[:,[0,1,2,3,'Gene']]
dat.to_csv('ReadyForTajimasD.csv')

#print('Reading in data',flush=True)
#dat = pd.read_csv('ReadyForTajimasD.csv')
print('Subsetting to valid genes',flush=True)
datt = dat.ix[dat.Gene.isin(dat.Gene.value_counts().reset_index(name='count').query('count > 4')['index'])]
dattt = datt.iloc[:,4:(dat.shape[1]-1)]

def iterGenes(subset,subname):
    indx = 0
    PWdiff = {}
    dT = {}
    dW = {}
    print('Iterating genes',flush=True)
    for eachGene in subset.Gene.unique():
        indx += 1
        if indx % 50 == 0:
            print('Progress ('+subname+'): '+str(indx)+' of '+str(subset.Gene.unique().size)+'genes',flush=True)
        subsett = subset.iloc[:,4:(dat.shape[1]-1)]
        currGene = subsett.ix[subset.Gene == eachGene]
        n = currGene.shape[1]
        PWdiffTotal = 0
        for i1 in range(0,n):
            for i2 in range(0,n):
                if i2 > i1:
                    PWdiffTotal += sum(currGene.ix[:,i1] != currGene.ix[:,i2])
        PWdiff[eachGene] = PWdiffTotal
        dT[eachGene] = PWdiff[eachGene]/(n*(n-1)/2)
        segSites = currGene.sum(axis=1) / currGene.shape[1]
        segSites = sum(dattt.head().sum(axis=1) != 1 | 0)
        harmonic = sum([1/x for x in range(1,n)])
        dW[eachGene] = segSites/harmonic

    PWdiff = pd.DataFrame.from_dict(PWdiff,orient='index')
    dT = pd.DataFrame.from_dict(dT,orient='index')
    dW = pd.DataFrame.from_dict(dW,orient='index')
    TajD = pd.concat([PWdiff,dT,dW],axis=1)
    TajD.columns = ['PWdiff','dT','dW']
    TajD.to_csv(subname+'.csv')

# Create subsets and their subprocesses
r1 = 0					# gene subset range index 1
cindx = 0 				# subprocess index
childName='ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'	# for subprocess file names
jobs = []				# for holding jobs for wait check

for x in range(0,numThreads):
    r2 = r1 + (len(dat.Gene.unique()) // numThreads)
    if x == (numThreads - 1):
        r2 += len(dat.Gene.unique()) % numThreads
    childSubset = datt.ix[datt.Gene.isin(dat.Gene[r1:r2])]
    r1 = r2
    p = multiprocessing.Process(target=iterGenes,args=(childSubset,'sub'+childName[cindx],))
    cindx += 1
    jobs.append(p)
    p.start()

# Wait for all subprocesses to finish
for job in jobs:
    job.join()

# Read in all subset files
for each in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
    if each == 'A':
        mergedSubs = pd.read_csv('subA.csv')
    else:
        tmpSubCSV = pd.read_csv('sub'+each+'.csv')
        mergedSubs = pd.concat([mergedSubs,tmpSubCSV])
# Consider adding a system call to remove subfiles, or store as variables instead.

print("Calculating Tajima's D",flush=True)
mergedSubs['D'] = (mergedSubs['dT'] - mergedSubs['dW']) / (np.std(mergedSubs['dT'] - mergedSubs['dW']))

print('Writing output',flush=True)
mergedSubs.to_csv('TajimasD_py.csv')
