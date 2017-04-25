# vcfToSNPs.py
# For pulling raw SNP calls from a vcf file.

import pandas as pd

vcf = 'in.vcf'
nChromosomes = 22
outName = 'SNP_output'

for chrom in range(1,nChromosomes+1):
    # Read in vcf, via a generator, filtering for a given chromosome
    iter_csv = pd.read_csv(vcf,iterator=True,chunksize=1000,header=None,comment='#',sep='\t',dtype='str')
    df = pd.concat([chunk[chunk.ix[:,0] == str(chrom)] for chunk in iter_csv])

    # SNP matrix; rows as SNPs, columns as samples
    for ind in range(9,df.shape[1]):
        if ind == 9:
            SNPS = pd.DataFrame(pd.DataFrame(df.ix[:,ind].str.split(':').tolist()).ix[:,0].str.split('/').tolist())
        else:
            newCols = pd.DataFrame(pd.DataFrame(df.ix[:,ind].str.split(':').tolist()).ix[:,0].str.split('/').tolist())
            SNPS = pd.concat([SNPS,newCols],axis=1)
    SNPS.to_csv(outName+chrom+'.snps',sep=' ',header=False,index=False)

    # map, rows as SNPs, columns as chr-ID-bp-bp
    Map = pd.DataFrame(data={'A':df.ix[:,0],'B':'ID_'+df.ix[:,1],'C':df.ix[:,1],'D':df.ix[:,1]})
    Map.to_csv(outName+chrom+'.map',sep=' ',header=False,index=False)
