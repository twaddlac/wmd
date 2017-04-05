#!/usr/local/bin/python2.7

import sys
import subprocess
import argparse
import re
import cogent
from cogent.db.ncbi import EFetch
import cogent.parse.genbank
import matplotlib
import pysam
from Bio import SeqIO
import os
import time
import ggplot
import math

def main():

	global args
	global refName
	global refNamePrefix

	parser = argparse.ArgumentParser(description='Execute modular SNV analysis pipeline.')
	parser.add_argument('--name', '-n', metavar='name', nargs='+', help='Sample names (exculding .fastq)', required=True)
	parser.add_argument('--mappers', '-m', metavar='map1', nargs='+',help='Mapping algorithms to execute (default: novoalign)', default='novoalign', choices=['bowtie2','novoalign','bwa'])
	#parser.add_argument('--variant_callers', '-v', metavar='vc1', nargs='+', help='Variant calling algorithms to execute (default gconf)', default='gconf',\
	#choices=['samtools','gconf','gatk','SOAPSnp'])
	parser.add_argument('--reference','-r', metavar='ref', nargs=1, help='Fasta reference (gi, accession number, or fasta file). (default: hg18)', default='hg18')
	#parser.add_argument('--in', '-i', metavar='seq', nargs=1, help='Query read file', required=True)
	parser.add_argument('--threads', '-t', metavar='t', nargs=1, help='Number of threads to run for each program (default: 8)', default=8, type=int)
	parser.add_argument('--all', '-a', help='Run all mappers and variant callers', action='store_true')
	parser.add_argument('--mincov','-c',help='Minimum coverage to call variant',default=4, type=int)
	parser.add_argument('--paired','-p',help='Samples are paired-end/mate-pair. Prefix MUST be the same between both samples with corresponding *.1 *.2 appendices',action='store_true')
	parser.add_argument('--minqual','-q',help='Minimum mean quality score to filter reads (default: 12)', default=12, type=int)
	parser.add_argument('--illumina', help='Illumina quality encoding', action='store_true')
	parser.add_argument('--verbose', '-V' ,help='Verbose output', action='store_true')
	args = parser.parse_args()

	if os.path.isfile('wmd-gconf-bam-files.txt'):
	        os.remove('wmd-gconf-bam-files.txt')

	if '.fasta' in args.reference[0] or '.fa' in args.reference[0]:
		refName = args.reference[0]
	else:
		refName = args.reference[0]+'.fasta'
                print 'Downloading reference from NCBI'
                ef = EFetch(id=args.reference[0], rettype='fasta')
                lines = ef.read().splitlines()
                with open(refName, 'w') as f:
                	for line in lines:
                		f.write(line+'\n')
                print 'Download successful'

	print 'Indexing reference'
	refNamePrefix = refName.split('.fa')
	bowtie2index = ['bowtie2-build','-q',refName,refName]
        novoindex = ['novoindex',refNamePrefix[0]+'.ndx',refName]
        bwaindex = ['bwa','index','-a','bwtsw',refName]
	samtoolsindex = ['samtools','faidx',refName]
	indexCommands = []

	if 'bowtie2' in args.mappers or args.all\
	and not os.path.exists(refName+'.1.bt2')\
	and not os.path.exists(refName+'.2.bt2')\
	and not os.path.exists(refName+'.3.bt2')\
	and not os.path.exists(refName+'.4.bt2')\
	and not os.path.exists(refName+'.rev.1.bt2')\
	and not os.path.exists(refName+'.rev.2.bt2'): 		indexCommands.append(bowtie2index)

	if 'novoalign' in args.mappers or args.all\
	and not os.path.exists(refNamePrefix[0]+'.ndx'): 	indexCommands.append(novoindex)

	if 'bwa' in args.mappers or args.all\
	and not os.path.exists(refName+'.pac')\
	and not os.path.exists(refName+'.ann')\
	and not os.path.exists(refName+'.bwt')\
	and not os.path.exists(refName+'.amb')\
	and not os.path.exists(refName+'.pac'):			indexCommands.append(bwaindex)
	
	if not os.path.exists(refName+'.fai'):			indexCommands.append(samtoolsindex)
	
	processes=[]
	for index in indexCommands:
		processes.append(subprocess.Popen(index))

	exit_codes =[p.wait() for p in processes]
	
        print 'Done indexing reference'

	if not os.path.exists(refNamePrefix[0]+'.dict'):
		print 'Making sequence dictionary'
		subprocess.call(['java','-jar','/opt/wmd/picard-tools-1.104/CreateSequenceDictionary.jar','R=',refName,'O=',refNamePrefix[0]+'.dict','QUIET=','true'])

	for sampleName in args.name:

	        # echo "echo blah" | qsub

	        #qsub -hold_jid job1, job2 ... to wait for all of the jobs to be completed
        	#Popen.wait = checks to see if the Popen process is done
	        # save all bam files and run gconf on them at once
	        print 'Starting\t',sampleName

	        print 'Quality filtereing reads'

		prinseq = ['prinseq-lite.pl','-fastq',sampleName+'.fastq','-min_len','50','-min_qual_mean',str(args.minqual),'-out_good',sampleName+'.filtered',\
                          '-out_bad','null']

		if args.paired:
			prinseq.insert(3, '-fastq2')
			prinseq.insert(4, sampleName+'.2.fastq')

		#if args.illumina:
			#prinseq.insert(1,'-phred64')

	        with open(sampleName+'.filtered.prinseq-stats.tsv','w') as stats:

			subprocess.call(prinseq,stderr=stats)
	                #if args.paired:

        	                #subprocess.call(['prinseq-lite.pl','-fastq',sampleName+'.1.fastq','fastq2',sampleName+'.2.fastq','-min_len','50','-min_qual_mean',str(args.minqual),\
				#'-out_good',sampleName+'.filtered','-out_bad','null'],stderr=stats)

                	#else:
                        	#subprocess.call(['prinseq-lite.pl','-fastq',sampleName+'.fastq','-min_len','50','-min_qual_mean',str(args.minqual),'-out_good',sampleName+'.filtered',\
				#'-out_bad','null'],stderr=stats)
	        #subprocess.call(['fastq_quality_filter','-Q33', '-q 20', '-p 75', '-i',i+'.fastq', '-o',i+'.filtered.fastq'])
        	print 'Done filtering'

	        if not args.paired:
	                print 'Calculating read quality statistics'
        	        if not os.path.exists(sampleName+'_filtered_fastqc'):
                	        os.makedirs(sampleName+'_filtered_fastqc')
	                subprocess.call(['fastqc','-t', str(args.threads), '-o', sampleName+'_filtered_fastqc', sampleName+'.filtered.fastq'])

        	        print 'Done calculating read quality statistics'

	        sampleName = sampleName+'.filtered'

               	#subprocess.call(['ln','-s',refName,sampleName+'.ref.fasta'])
		#subprocess.call(['ln','-s',refName+'.fai',sampleName+'.ref.fasta.fai'])

	        with open(sampleName+'.log','w') as log:
        	        log.write('Starting\t'+sampleName+'\n')

	        sample = WMD(sampleName)

        	if 'bowtie2' in args.mappers or args.all:
                	#bowtie2(i)
	                sample.runBowtie2()
        	        gatkPreProcessing(sampleName,'bowtie2')
                	#gconf(i,'bowtie2')
	                MapStats(sampleName,'bowtie2')
        	if 'novoalign' in args.mappers or args.all:
                	#novoalign(i)
	                sample.runNovoalign()
        	        gatkPreProcessing(sampleName,'novoalign')
                	#gconf(i,'novoalign')
	                MapStats(sampleName,'novoalign')
        	if 'bwa' in args.mappers or args.all:
                	#bwa(i)
	                sample.runBwa()
        	        gatkPreProcessing(sampleName,'bwa')
                	#gconf(i,'bwa')
	                MapStats(sampleName,'bwa')

        	#if '.fasta' in args.reference[0] or '.fa' in args.reference[0]:
                #	subprocess.call(['unlink',sampleName+'.ref.fasta'])

	gconf('wmd-gconf-bam-files.txt')

	return



class WMD(object):

	def __init__(self, name):
		self.name = name
		#index commands
		#self.bowtie2index = ['bowtie2-build','-q',self.name+'.ref.fasta',self.name+'.ref.fasta']
		#self.novoindex = ['novoindex',self.name+'.ref.ndx',self.name+'.ref.fasta']
		#self.bwaindex = ['bwa','index','-a','bwtsw',name+'.ref.fasta']

		#mapping commands
		self.bowtie2 = ['bowtie2','-p',str(args.threads),'--very-sensitive-local',refName,self.name+'.fastq','-S',self.name+'.bowtie2.sam']
		self.bowtie2pe = ['bowtie2','-p',str(args.threads),'--very-sensitive-local',refName,'-1',self.name+'.1.fastq','-2',self.name+'.2.fastq','-S',self.name+'.bowtie2.sam']
		#self.novoalign = ['novoalign','-o','SAM','-d',self.name+'.ref.ndx','-f',self.name+'.fastq']
		#self.novoalignpe = ['novoalign','-o','SAM','-d',self.name+'.ref.ndx','-f',self.name+'.1.fastq',self.name+'.2.fastq']
		self.bwa = ['bwa','mem','-t',str(args.threads),refName,name+'.fastq']
		self.bwape = ['bwa','mem','-t',str(args.threads),refName,name+'.1.fastq',name+'.2.fastq']
		
		if args.illumina:
			#self.bowtie2.insert(1,'--phred64')
			#self.bowtie2pe.insert(1,'--phred64')
			self.bowtie2pe.indert(1,'--rf')
			#self.novoalign.insert(1,'-F')
			#self.novoalign.insert(2,'ILMFQ')
			self.novoalign.insert(3,'-i')
			self.novoalign.insert(4,'MP')
			#self.novoalignpe.insert(1,'-F')
                        #self.novoalignpe.insert(2,'ILMFQ')
                        self.novoalignpe.insert(3,'-i')
                        self.novoalignpe.insert(4,'MP')

	#def indexIt(self,command):
		#print 'Indexing '+self.name+' '+command[0]
                #with open(self.name+'.log','a') as log:
                #        subprocess.call(command,stderr=log, stdout=log)
		#print 'Finished indexing '+self.name+' '+command[0]


	def mapIt(self,command):
		print 'Aligning '+self.name+' with '+command[0]
		stdoutSam = ('bwa')
		if any(sam in command[0] for sam in stdoutSam): # if mapper outputs sam to STDOUT
			with open(self.name+'.'+command[0]+'.sam','w') as outfile,open(self.name+'.log','a') as log:
				subprocess.call(command,stderr=log, stdout=outfile)
		else:
			with open(self.name+'.log','a') as log:
				subprocess.call(command,stderr=log, stdout=log)
		print 'Finished aligning '+self.name+' with '+command[0]

	
	def runBowtie2(self):
		name = self.name
		#self.indexIt(self.bowtie2index)
		if args.paired:
			self.mapIt(self.bowtie2pe)
		else:
			self.mapIt(self.bowtie2)
		samToBam(name,'bowtie2')
		#subprocess.call(['rm',name+'.ref.fasta.1.bt2',name+'.ref.fasta.2.bt2',name+'.ref.fasta.3.bt2',name+'.ref.fasta.4.bt2',name+'.ref.fasta.rev.1.bt2',name+'.ref.fasta.rev.2.bt2'])
	        return


	def runNovoalign(self):
		name = self.name
		#self.indexIt(self.novoindex)
		wc = ""
		samFiles = []
		processess = []
                if args.paired:
			wc = subprocess.Popen(['wc','-l',name+'.1.fastq'],stdout=subprocess.PIPE)
			out = wc.communicate()
			wc.wait()
			lines = str(math.trunc(math.ceil(int(out[0].split(' ')[0])/int(args.threads))))
			first = subprocess.Popen(['split','-d','-l',lines,name+'.1.fastq',name+'.1.fastq.'])
			second = subprocess.Popen(['split','-d','-l',lines,name+'.2.fastq',name+'.2.fastq.'])
			first.wait()
			second.wait()
			for i in range(0,int(args.threads)):
				if int(len(str(i))) == 1:
                                        suffix = '0'+str(i)
                                else:
                                        suffix = str(i)
				with open(self.name+'.novoalign.'+str(i),'w') as sam:
					samFiles.append(self.name+'.novoalign.'+str(i))
					processess.append(subprocess.Popen(['novoalign','-o','SAM','-d',refNamePrefix[0]+'.ndx','-f',self.name+'.1.fastq.'+suffix,self.name+'.2.fastq.'+suffix], stdout=sam))
                        #self.mapIt(self.novoalignpe)      
                else:
			wc = subprocess.Popen(['wc','-l',name+'.fastq'],stdout=subprocess.PIPE)
			out = wc.communicate()
			wc.wait()
                        lines = str(math.trunc(math.ceil(int(out[0].split(' ')[0])/int(args.threads))))
                        subprocess.call(['split','-d','-l',lines,name+'.fastq',name+'.fastq.'])
                        for i in range(0,int(args.threads)):
				if int(len(str(i))) == 1:
					suffix = '0'+str(i)
				else:
					suffix = str(i)
                                with open(self.name+'.novoalign.'+str(i),'w') as sam:
					samFiles.append(self.name+'.novoalign.'+str(i))
                                        processess.append(subprocess.Popen(['novoalign','-o','SAM','-d',refNamePrefix[0]+'.ndx','-f',self.name+'.fastq.'+suffix], stdout=sam))
                exit_codes = [p.wait() for p in processess]

		subprocess.call('rm *fastq.*', shell=True)

		mergeCommand = ['java','-jar','/opt/wmd/picard-tools-1.104/MergeSamFiles.jar']
		for i in samFiles:
			mergeCommand.append('I=')
			mergeCommand.append(i)
		mergeCommand.extend(['O=',name+'.novoalign.sam','USE_THREADING=','true'])
		subprocess.call(mergeCommand)

		for sam in samFiles:
			subprocess.call(['rm',sam])

			#self.mapIt(self.novoalign)
		samToBam(name,'novoalign')
		#subprocess.call(['rm',name+'.ref.ndx'])
       		print 'Done with novoalign'
	        return


	def runBwa(self):
		name = self.name
		#self.indexIt(self.bwaindex)
		if args.paired:
			self.mapIt(self.bwape)
		else:
			self.mapIt(self.bwa)
		samToBam(name,'bwa')
		#subprocess.call(['rm',name+'.ref.fasta.amb',name+'.ref.fasta.ann',name+'.ref.fasta.bwt',name+'.ref.fasta.pac',name+'.ref.fasta.sa'])
		print 'Done with BWA'
		return


def samToBam(name,mapper):

	print 'Converting '+name+'-'+mapper+' SAM to BAM'	

	with open(name+'.log','a') as log:
		log.write('******sam2bam******')	
		with open(name+'.'+mapper+'.temp', 'w') as outfile:
        	        subprocess.call(['samtools','view','-bS',name+'.'+mapper+'.sam'], stdout=outfile,stderr=log)
        subprocess.call(['samtools','sort',name+'.'+mapper+'.temp',name+'.'+mapper])
        subprocess.call(['samtools','index',name+'.'+mapper+'.bam'])

	subprocess.call(['rm',name+'.'+mapper+'.temp',name+'.'+mapper+'.sam'])

	print 'Finished converting '+name+'-'+mapper+' SAM to BAM'
	
	return

def gconf(fof):
	print 'Starting GConf'

	subprocess.call(['gc.0.12.pl','--fof',fof,'-r',refName,'--samtools','/opt/wmd/samtools-0.1.19/samtools','-w','-z','-l','/opt/wmd/gc_pvalue_lookup.transposed.tsv',\
	'-c',str(args.mincov),'-t',str(args.threads)])
    
	print 'Finished GConf'
    
	return

def gatkPreProcessing(name,mapper):

	print 'Starting GATK pre-processing'
	with open(name+'.log','a') as log, open('wmd-gconf-bam-files.txt','a') as bams:
		log.write('******GATK Preprocessing******')

		print 'Adding ReadGroups'
		subprocess.call(['java','-jar','/opt/wmd/picard-tools-1.104/AddOrReplaceReadGroups.jar','I=',name+'.'+mapper+'.bam','O=',name+'.'+mapper+'.rg.bam','LB=','Ion',\
		'PL=','Ion','PU=','Ion','SM=',name],stdout=log,stderr=log)

		print 'Indexing bam'
		subprocess.call(['samtools','index',name+'.'+mapper+'.rg.bam'])

		print 'Finding indels'
		subprocess.call(['java','-jar','/opt/wmd/GenomeAnalysisTK-2.7-4-g6f46d11/GenomeAnalysisTK.jar','-T','RealignerTargetCreator','-R',refName,\
		'-I',name+'.'+mapper+'.rg.bam','-o',name+'.'+mapper+'.rg.intervals'],stdout=log,stderr=log)

		print 'Realigning indels'
		subprocess.call(['java','-jar','/opt/wmd/GenomeAnalysisTK-2.7-4-g6f46d11/GenomeAnalysisTK.jar','-T','IndelRealigner','-R',refName,\
		'-I',name+'.'+mapper+'.rg.bam','-targetIntervals',name+'.'+mapper+'.rg.intervals','-o',name+'.'+mapper+'.rg.realigned.bam'],stdout=log, stderr=log)

		bams.write(name+'.'+mapper+'.rg.realigned.bam\n')

		subprocess.call(['rm',name+'.'+mapper+'.bam',name+'.'+mapper+'.bam.bai',name+'.'+mapper+'.rg.intervals',name+'.'+mapper+'.rg.bam',name+'.'+mapper+'.rg.bam.bai'])

	return

def MapStats(name,mapper):
    
    print 'Starting MapStats'

    try:   
	bamfile = pysam.Samfile(name+'.'+mapper+'.rg.realigned.bam', 'rb')
    except IOError:
	print "Bam file didn't open"
   
    mapped = bamfile.mapped
    unmapped = bamfile.unmapped
    totalreadbases = 0
    refs = bamfile.references
    reflens = bamfile.lengths
    totalrefbases = 0
    totalcov = 0
    refbasescovered = 0
    readbasesmapped = 0
    uss = {} # unique start sites hash
    usscount = 0
    
    for seq in SeqIO.parse(name+".fastq","fastq"):
        totalreadbases += len(seq.seq)

    for i in range(0,len(reflens)):
        totalrefbases += reflens[i]
        for p in bamfile.pileup(refs[i],1,reflens[i]):
            if p.n > 0:
                refbasescovered += 1
            totalcov += p.n
        
        for alignedread in bamfile.fetch(refs[i],1,reflens[i]):
            #print "tid",alignedread.tid,"\taend",alignedread.aend,"\talen",alignedread.alen,"\tpos",alignedread.pos,"\tpositions",alignedread.positions
            uss[refs[i]+str(alignedread.aend)] = 1
            uss[refs[i]+str((alignedread.pos - 1))] = 1
            
    for key in uss.keys():
        usscount += 1
   
    with open(name+'.'+mapper+'.rg.realigned.map-stats.tsv', 'w') as stats:
	#out = "Percent Reads Mapped = %s\nPercent Ref Bases Covered = %s\nPercent Read Bases Mapped = %s\nAverage Coverage = %s\nPercent Unique Start Sites = %s" % \
        #(((mapped/(mapped + unmapped)) * 100), ((refbasescovered/totalrefbases)*100),((totalcov/totalreadbases)*100),((totalcov/totalrefbases)*100),\
        #((usscount/totalrefbases)*100))
	stats.write("Percent Reads Mapped = %s\nPercent Ref Bases Covered = %s\nPercent Read Bases Mapped = %s\nAverage Coverage = %s\nPercent Unique Start Sites = %s" % \
	(((mapped/(mapped + unmapped)) * 100), ((refbasescovered/totalrefbases)*100),((totalcov/totalreadbases)*100),((totalcov/totalrefbases)*100),\
	((usscount/totalrefbases)*100)))
	#stats.write(out)

    
    
    #for p in bamfile.pileup():
    #    print p.pos,"\t",p.n
    
    bamfile.close()

    print "Done calculating MapStats"

    return

if __name__ == '__main__':
	main()
