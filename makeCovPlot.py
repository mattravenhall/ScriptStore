#!/home/matt/anaconda3/bin/python3

#################################################################################
# Process Input #
#################

# Import all the things please
import os
import sys
import subprocess
import pandas as pd
import matplotlib.pyplot as plt

covStore = '/path/to/full/cov/files/'
outDirectory = './'
threads = 5
title = ''
filename='test'
neighbours=[]
figureX = 15
figureY = 2

if len(sys.argv) == 1 or sys.argv[1].lower() == 'help':
	print('Coverage Plotter\n'+
	'Required:\n'+
	'\t--region=chromosome:basepair-basepair\n'+
	'\t--samples=ID,ID,ID,ID\n'+
	'\t--genes=basepair-basepair,basepair-basepair\n'+
	'\t--SVs=basepair-basepair,basepair-basepair\n'+
	'Optional:\n'+
	'\t--threads=int\n'+
	'\t--covStore=/path/to/dir\n'+
	'\t--title=string\n'+
	'\t--filename=string\n'+
	'\t--figureX=int [15]\n'+
	'\t--figureY=int [2]\n'+
	'\t--neighbours=basepair-basepair,basepair-basepair\n')
	sys.exit()

def checkArg(valueGiven, givenArg, isBool=False):
    try:
        outValue = valueGiven.split('=')[1]
    except IndexError:
        print('\033[91mError: No argument supplied to '+givenArg+'.\033[0m\n',flush=True)
        sys.exit()
    if outValue == '':
        print('\033[91mError: No argument supplied to '+givenArg+'.\033[0m\n',flush=True)
        sys.exit()
    if isBool and (outValue.lower() in ['t','true']):
        return(True)
    elif isBool and (outValue.lower() in ['f','false']):
        return(False)
    elif isBool:
        print('\033[91mError: Incorrect argument "'+outValue+'" supplied to '+givenArg+'.\033[0m\n',flush=True)
        sys.exit()
    else:
        return(outValue)

for arg in sys.argv[1:]:
	if arg.startswith('--region'):
		chromosome, region = checkArg(arg, '--region').split(':')
		regionStart, regionEnd = region.split('-')
	elif arg.startswith('--samples'):
		samples = checkArg(arg, '--samples').split(',')
	elif arg.startswith('--genes'):
		genes = checkArg(arg, '--genes').split(',')
	elif arg.startswith('--neighbours'):
		neighbours = checkArg(arg, '--neighbours').split(',')
	elif arg.startswith('--SVs'):
		SVs = checkArg(arg, '--SVs').split(',')
	elif arg.startswith('--threads'):
		threads = int(checkArg(arg, '--threads'))
	elif arg.startswith('--covStore'):
		covStore = checkArg(arg, '--covStore')
		if not os.path.isdir(covStore):
			print('Given covStore directory does not exist')
			sys.exit()
	elif arg.startswith('--outDir'):
		outDirectory = checkArg(arg, '--outDir')
		if not os.path.isdir(outDirectory):
			print('Given outDir directory does not exist')
			sys.exit()
	elif arg.startswith('--title'):
		title = checkArg(arg, '--title')
	elif arg.startswith('--filename'):
		filename = checkArg(arg, '--filename')
	elif arg.startswith('--figureX'):
		figureX = int(checkArg(arg, '--figureX'))
	elif arg.startswith('--figureY'):
		figureY = int(checkArg(arg, '--figureY'))
	else:
		print('\033[91mArgument {0} not recognised.\033[0m\n'.format(arg))
		sys.exit()

missingVar = []
for var in ['region', 'samples', 'genes', 'SVs']:
	if var not in locals():
		missingVar.append(var)
if len(missingVar) > 0:
	print('Missing required variables: {0}'.format(missingVar))
	sys.exit()

print('Running Coverage Pipeline For:')
print(' Region: {0} ({1:,} bp)'.format(region,int(regionEnd)-int(regionStart)))
print(' Sample/s: {0}'.format(samples))
print(' Gene/s: {0}'.format(genes))
print(' SV/s: {0}'.format(SVs))
print(' Neighbouring Genes: {0}'.format(neighbours))
print(' Threads: {0}'.format(threads))
print(' Path To Coverage: {0}'.format(covStore))
print(' Figure Title: {0}'.format(chromosome+':'+region if title == '' else title))
print(' Output Filename: {0}.png'.format(filename))
print('')

#################################################################################
# Grab Coverage #
#################

print('Grabbing Coverage Subsets')
commands = ["awk -F $'\t' '{ if ($1 == \""+chromosome+"\" && $2 > "+regionStart+" && $2 < "+regionEnd+") print $0 }' <(zcat -dc "+covStore+sample+".coverage.gz) | gzip > "+sample+'_'+chromosome+':'+region+".cov.gz" for sample in samples]
processes = [subprocess.Popen(['/bin/bash','-c',cmd]) for cmd in commands]
print(' Waiting for {0} subprocess(es).'.format(len(processes)))
for p in processes: p.wait()
print('')

#################################################################################
# Plot Coverage #
#################

print('Creating Coverage Plot')
fig, ax = plt.subplots(nrows=len(samples), ncols=1,figsize=(figureX,figureY*len(samples)))

if len(samples) == 1:
	sample = samples[0]
	covFile = pd.read_csv('{0}_{1}:{2}-{3}.cov.gz'.format(sample,chromosome,regionStart,regionEnd), compression='gzip', header=None, sep='\t')

	ax.plot(covFile[1], covFile[2], lw=0.5)
	ax.set_title(chromosome+':'+region if title == '' else title)
	ax.set_xlabel('Position (bp)')
	ax.set_ylabel(sample)

	for s in SVs:
		start, end = s.split('-')
		ax.axvspan(start, end, color='red', alpha=0.1)
	for g in genes:
		start, end = g.split('-')
		ax.axvspan(start, end, color='green', alpha=0.2)
	for n in neighbours:
		start, end = n.split('-')
		ax.axvspan(start, end, color='grey', alpha=0.2)
	ax.axhline(y=0,color='black',lw=0.5)

else:
	for i, sample in enumerate(samples):
		covFile = pd.read_csv('{0}_{1}:{2}-{3}.cov.gz'.format(sample,chromosome,regionStart,regionEnd), compression='gzip', header=None, sep='\t')

		ax[i].plot(covFile[1], covFile[2], lw=0.5)
		if i == 0: ax[i].set_title(chromosome+':'+region if title == '' else title)
		ax[i].set_xlabel('Position (bp)')
		ax[i].set_ylabel(sample)

		for s in SVs:
			start, end = s.split('-')
			ax[i].axvspan(start, end, color='red', alpha=0.1)
		for g in genes:
			start, end = g.split('-')
			ax[i].axvspan(start, end, color='green', alpha=0.2)
		for n in neighbours:
			start, end = n.split('-')
			ax[i].axvspan(start, end, color='grey', alpha=0.2)
		ax[i].axhline(y=0,color='black',lw=0.5)

fig.savefig(filename+'.png', bbox_inches='tight')
plt.close(fig)

print('')
print('Done!')
