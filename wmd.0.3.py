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
#import ggplot
import math
import smtplib
import jinja2

def main():

	global args
	global refName
	global refNamePrefix
	global mappingStatsTable
	
	readStatsTable = []
	filteredStatsTable = []
	mappingStatsTable = {}
	gconfStatsTable = {}

	templateLoader = jinja2.FileSystemLoader( searchpath="/opt/wmd/" )
	templateEnv = jinja2.Environment( loader=templateLoader )

	TEMPLATE_FILE = 'wmd-template.html'
	template = templateEnv.get_template(TEMPLATE_FILE)


	# recognize zipped files
	# zip filtered fastq files
	# correct % read bases mapped
	# calculate percent valid unique
	# output results to table
	# add MrFast, crack, soap
	# enable output of each subroutine to be passed to whatever subsequent subroutine
	# find longest common substring for paired end
	# get valid pairs
	parser = argparse.ArgumentParser(description='Execute modular SNV analysis pipeline.')
	parser.add_argument('--name', '-n', metavar='name', nargs='+', help='Sample names (exculding .fastq)', required=True)
	parser.add_argument('--mappers', '-m', metavar='map1', nargs='+',help='Mapping algorithms to execute (default: novoalign)', default='novoalign', choices=['bowtie2','novoalign','bwa'])
	parser.add_argument('--variant_callers', '-v', metavar='vc1', nargs='+', help='Variant calling algorithms to execute (default gconf)', default='gconf',\
	choices=['samtools','gconf','gatk','varscan'])
	parser.add_argument('--reference','-r', metavar='ref', nargs=1, help='Fasta reference (gi, accession number, or fasta file).', required=True)
	parser.add_argument('--threads', '-t', metavar='t', nargs=1, help='Number of threads to run for each program (default: 8)', default=8, type=int)
	parser.add_argument('--all', '-a', help='Run all mappers and variant callers', action='store_true')
	parser.add_argument('--mincov','-c',help='Minimum coverage to call variant',default=4, type=int)
	parser.add_argument('--paired','-p',help='Samples are paired-end/mate-pair. Prefix MUST be the same between both samples with corresponding *.1 *.2 appendices',action='store_true')
	parser.add_argument('--minqual','-q',help='Minimum mean quality score to filter reads (default: 12)', default=12, type=int)
	parser.add_argument('--illumina', help='Illumina quality encoding', action='store_true')
	parser.add_argument('--percent',help='Percent of reads to use per fastq file', type=float, default=100)
	parser.add_argument('--delete_indices', help='Delete reference indices. (default:off)', action='store_true')
	parser.add_argument('--email', help='Send email to addresses when run has completed.', type=str, nargs='+', metavar='user@mail.com')
	parser.add_argument('--verbose', '-V' ,help='Verbose output', action='store_true')
	args = parser.parse_args()

	# use glob to find paired end files

	if args.all:
		args.mappers = ['bowtie2','novoalign','bwa']

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

	if 'bowtie2' in args.mappers\
	and not os.path.isfile(refName+'.1.bt2')\
        and not os.path.isfile(refName+'.2.bt2')\
        and not os.path.isfile(refName+'.3.bt2')\
        and not os.path.isfile(refName+'.4.bt2')\
        and not os.path.isfile(refName+'.rev.1.bt2')\
        and not os.path.isfile(refName+'.rev.2.bt2'):		indexCommands.append(bowtie2index)

	if 'novoalign' in args.mappers\
	and not os.path.isfile(refNamePrefix[0]+'.ndx'): 	indexCommands.append(novoindex)

	if 'bwa' in args.mappers\
	and not os.path.isfile(refName+'.pac')\
	and not os.path.isfile(refName+'.ann')\
	and not os.path.isfile(refName+'.bwt')\
	and not os.path.isfile(refName+'.amb')\
	and not os.path.isfile(refName+'.pac'):			indexCommands.append(bwaindex)
	
	if not os.path.isfile(refName+'.fai'):			indexCommands.append(samtoolsindex)
	
	processes=[]
	for index in indexCommands:
		processes.append(subprocess.Popen(index))

	exit_codes =[p.wait() for p in processes]
	
        print 'Done indexing reference'

	if not os.path.isfile(refNamePrefix[0]+'.dict'):
		print 'Making sequence dictionary'
		subprocess.call(['java','-jar','/opt/wmd/picard-tools-1.104/CreateSequenceDictionary.jar','R=',refName,'O=',refNamePrefix[0]+'.dict','QUIET=','true'])

	#remove trailing .fastq
	
	temp = [x.split('.fastq')[0] for x in args.name]
	args.name = temp

	for sampleName in args.name:

	        # echo "echo blah" | qsub

	        #qsub -hold_jid job1, job2 ... to wait for all of the jobs to be completed
        	#Popen.wait = checks to see if the Popen process is done
	        # save all bam files and run gconf on them at once

		#prefix = sampleName.split('.fastq')
		#sampleName = prefix[0]
		
	        print 'Starting\t',sampleName

	        print 'Quality filtereing reads'

		prinseq = ['prinseq-lite.pl']

		if args.paired:
			prinseq.extend(['-fastq',sampleName+'.1.fastq','-fastq2',sampleName+'.2.fastq'])
		else:
			prinseq.extend(['-fastq',sampleName+'.fastq'])

		prinseq.extend(['-min_len','50','-min_qual_mean',str(args.minqual),'-out_good',sampleName+'.filtered','-out_bad','null'])

		#if args.illumina:
			#prinseq.insert(1,'-phred64')
		
	        #with open(sampleName+'.filtered.prinseq-stats.tsv','w') as stats:

		prinseqProc = subprocess.Popen(prinseq,stderr=subprocess.PIPE)
		prinseqOut = prinseqProc.communicate()
		data = prinseqOut[1].split('\n\t')
		data.pop(0)
		#data.remove('Sequences filtered by specified parameters:')
		#print data
		filtered = False
		for entry in data:
			row = entry.split(':')

			if 'Sequences filtered by specified parameters:' in entry:
				filtered = True

			elif filtered:
			
				filteredStatsTable.append([row[0].strip(),row[1].strip()])

			else:
				readStatsTable.append([row[0].strip(),row[1].strip()])

		#print readStatsTable

		if args.paired:
			subprocess.call(['mv',sampleName+'.filtered_1.fastq',sampleName+'.1.filtered.fastq'])
			subprocess.call(['mv',sampleName+'.filtered_2.fastq',sampleName+'.2.filtered.fastq'])
	                #if args.paired:
		#return
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

	        #sampleName = sampleName+'.filtered'

               	#subprocess.call(['ln','-s',refName,sampleName+'.ref.fasta'])
		#subprocess.call(['ln','-s',refName+'.fai',sampleName+'.ref.fasta.fai'])

	        with open(sampleName+'.log','w') as log:
        	        log.write('Starting\t'+sampleName+'\n')

	        sample = WMD(sampleName)

        	if 'bowtie2' in args.mappers:
                	#bowtie2(i)
	                sample.runBowtie2()
        	        #gatkPreProcessing(sampleName,'bowtie2')
                	#gconf(i,'bowtie2')
	                #MapStats(sampleName,'bowtie2')
        	if 'novoalign' in args.mappers:
                	#novoalign(i)
	                sample.runNovoalign()
        	        #gatkPreProcessing(sampleName,'novoalign')
                	#gconf(i,'novoalign')
	                #MapStats(sampleName,'novoalign')
        	if 'bwa' in args.mappers:
                	#bwa(i)
	                sample.runBwa()
        	        #gatkPreProcessing(sampleName,'bwa')
                	#gconf(i,'bwa')
	                #MapStats(sampleName,'bwa')

        	#if '.fasta' in args.reference[0] or '.fa' in args.reference[0]:
                #	subprocess.call(['unlink',sampleName+'.ref.fasta'])

	gconf('wmd-gconf-bam-files.txt')

	templateVars = {'statsTitle' : 'WMD Read Statistics',
			'readStatsTable' : readStatsTable,
			'filteredTitle' : 'Sequences filtered by specified parameters',
			'filteredStatsTable' : filteredStatsTable
		       }

	outputHTML = template.render(templateVars)

	with open('wmd-stats.html','w') as html:
		html.write(outputHTML)

	if args.email:

		message = """From: WMD <localhost@vandelay.genetics.pitt.edu>
			     Subject: WMD Run Finished
	
			     Your WMD run has completed. DO NOT REPLY.
		          """

		s = smtplib.SMTP('localhost')
		s.sendmail('localhost@vandelay.genetics.pitt.edu',args.email,message)
		s.quit()

	return



class WMD(object):

	def __init__(self, name):
		self.name = name+'.filtered'
		self.first = name+'.1.filtered'
		self.second = name+'.2.filtered'
		#index commands
		#self.bowtie2index = ['bowtie2-build','-q',self.name+'.ref.fasta',self.name+'.ref.fasta']
		#self.novoindex = ['novoindex',self.name+'.ref.ndx',self.name+'.ref.fasta']
		#self.bwaindex = ['bwa','index','-a','bwtsw',name+'.ref.fasta']

		#mapping commands
		self.bowtie2 = ['bowtie2','-p',str(args.threads),'--very-sensitive-local',refName,self.name+'.fastq','-S',self.name+'.bowtie2.sam']
		self.bowtie2pe = ['bowtie2','-p',str(args.threads),'--very-sensitive-local',refName,'-1',self.first+'.fastq','-2',self.second+'.fastq','-S',self.name+'.bowtie2.sam']
		#self.novoalign = ['novoalign','-o','SAM','-d',self.name+'.ref.ndx','-f',self.name+'.fastq']
		#self.novoalignpe = ['novoalign','-o','SAM','-d',self.name+'.ref.ndx','-f',self.name+'.1.fastq',self.name+'.2.fastq']
		self.bwa = ['bwa','mem','-t',str(args.threads),refName,name+'.fastq']
		self.bwape = ['bwa','mem','-t',str(args.threads),refName,self.first+'.fastq',self.second+'.fastq']
	
		#variant calling commands
		self.varscan = ['java','-jar','/opt/wmd/varscan/VarScan.v2.3.7.jar','--min-coverage',str(args.mincov),'--output-vcf','1',name+'.'+mapper+'.mpileup.tsv']
		self.gatk = ['java','-jar','/opt/wmd/gatk/GenomeAnalysisTK.jar','-T','HaplotypeCaller','--nt',str(args.threads),'-R',refName,'-I',name+'.'+mapper
	
		if args.illumina:
			#self.bowtie2.insert(1,'--phred64')
			#self.bowtie2pe.insert(1,'--phred64')
			self.bowtie2pe.insert(1,'--rf')
			#self.novoalign.insert(1,'-F')
			#self.novoalign.insert(2,'ILMFQ')
			#self.novoalign.insert(3,'-i')
			#self.novoalign.insert(4,'MP')
			#self.novoalignpe.insert(1,'-F')
                        #self.novoalignpe.insert(2,'ILMFQ')
                        #self.novoalignpe.insert(3,'-i')
                        #self.novoalignpe.insert(4,'MP')

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
		gatkPreProcessing(name,'bowtie2')
		MapStats(name,'bowtie2')
		if args.delete_indices:	subprocess.call(['rm',name+'.ref.fasta.1.bt2',name+'.ref.fasta.2.bt2',name+'.ref.fasta.3.bt2',name+'.ref.fasta.4.bt2',name+'.ref.fasta.rev.1.bt2',name+'.ref.fasta.rev.2.bt2'])
	        return


	def runNovoalign(self):
		name = self.name
		first = self.first
		second = self.second
		#self.indexIt(self.novoindex)
		wc = ""
		samFiles = []
		processes = []
		# convert lines to # reads, divide by threads, multiply by 4 (to get # reads)
                if args.paired:
			wc = subprocess.Popen(['wc','-l',first+'.fastq'],stdout=subprocess.PIPE)
			out = wc.communicate()
			wc.wait()
			reads = int(out[0].split(' ')[0])/4
			readsPerFile = math.trunc(math.ceil(reads/int(args.threads)))
			linesPerFile = str(readsPerFile * 4)
			#cut = str(linesPerFile*4)
			#print out[0].split(' ')[0],reads,readsPerFile,linesPerFile
			#lines = str(math.trunc(math.ceil(int(out[0].split(' ')[0])/int(args.threads))))
			firstFastq = subprocess.Popen(['split','-d','-l',linesPerFile,first+'.fastq',first+'.fastq.'])
			secondFastq = subprocess.Popen(['split','-d','-l',linesPerFile,second+'.fastq',second+'.fastq.'])
			firstFastq.wait()
			secondFastq.wait()
			for i in range(0,int(args.threads)):
				if int(len(str(i))) == 1:
                                        suffix = '0'+str(i)
                                else:
                                        suffix = str(i)
				with open(name+'.novoalign.'+str(i),'w') as sam:
					samFiles.append(name+'.novoalign.'+str(i))
					#p = subprocess.Popen(['novoalign','-o','SAM','-d',refNamePrefix[0]+'.ndx','-f',first+'.fastq.'+suffix,second+'.fastq.'+suffix], stdout=sam)
					p = subprocess.Popen('novoalign -r Random -e 1 -o SAM -d '+refNamePrefix[0]+'.ndx -f '+first+'.fastq.'+suffix+' '+second+'.fastq.'+suffix+' | perl -pe \'s/\t\t/\t/g\'', stdout=sam, shell=True)
					#processes.append(subprocess.Popen(['novoalign','-o','SAM','-d',refNamePrefix[0]+'.ndx','-f',first+'.fastq.'+suffix,second+'.fastq.'+suffix], stdout=sam))
					processes.append(p)
                        #self.mapIt(self.novoalignpe)      
                else:
			wc = subprocess.Popen(['wc','-l',name+'.fastq'],stdout=subprocess.PIPE)
			out = wc.communicate()
			wc.wait()
			reads = int(out[0].split(' ')[0])/4
                        readsPerFile = math.trunc(math.ceil(reads/int(args.threads)))
                        linesPerFile = str(readsPerFile * 4)
                        #lines = str(math.trunc(math.ceil(int(out[0].split(' ')[0])/int(args.threads))))
                        subprocess.call(['split','-d','-l',linesPerFile,name+'.fastq',name+'.fastq.'])
                        for i in range(0,int(args.threads)):
				if int(len(str(i))) == 1:
					suffix = '0'+str(i)
				else:
					suffix = str(i)
                                with open(name+'.novoalign.'+str(i),'w') as sam:
					samFiles.append(name+'.novoalign.'+str(i))
                                        #processes.append(subprocess.Popen(['novoalign','-o','SAM','-d',refNamePrefix[0]+'.ndx','-f',name+'.fastq.'+suffix], stdout=sam))
					p = subprocess.Popen('novoalign -r Random -e 1 -o SAM -d '+refNamePrefix[0]+'.ndx -f '+name+'.fastq.'+suffix+' | perl -pe \'s/\t\t/\t/g\'', stdout=sam, shell=True)
					processes.append(p)
                exit_codes = [p.wait() for p in processes]

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
		gatkPreProcessing(name,'novoalign')
		MapStats(name,'novoalign')
		if args.delete_indices:	subprocess.call(['rm',name+'.ref.ndx'])
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
		gatkPreProcessing(name,'bwa')
		MapStats(name,'bwa')
		if args.delete_indices:	subprocess.call(['rm',name+'.ref.fasta.amb',name+'.ref.fasta.ann',name+'.ref.fasta.bwt',name+'.ref.fasta.pac',name+'.ref.fasta.sa'])
		print 'Done with BWA'
		return


	def runVarscan(self):
		name = self.name
		

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
		subprocess.call(['java','-jar','/opt/wmd/gatk/GenomeAnalysisTK.jar','-T','RealignerTargetCreator','-R',refName,\
		'-I',name+'.'+mapper+'.rg.bam','-o',name+'.'+mapper+'.rg.intervals'],stdout=log,stderr=log)

		print 'Realigning indels'
		subprocess.call(['java','-jar','/opt/wmd/gatk/GenomeAnalysisTK.jar','-T','IndelRealigner','-R',refName,\
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
   
	mapped = float(bamfile.mapped)
        unmapped = float(bamfile.unmapped)
        totalreadbases = float()
        refs = bamfile.references
        reflens = bamfile.lengths
        totalrefbases = float()
        totalcov = float()
        refbasescovered = float()
        readbasesmapped = float()
        uss = {} # unique start sites hash
        usscount = float()
	uniqueMappedReads = {}


        fastqFiles = []
	#temp = args.name
        if args.paired:
		fastqFiles.append(name+'.1.filtered')
		fastqFiles.append(name+'.2.filtered')
	else:
		fastqFiles.append(name)	

	for fastq in fastqFiles:
		for seq in SeqIO.parse(fastq+".fastq","fastq"):
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
		stats.write("Percent Reads Mapped = %.2f\nPercent Ref Bases Covered = %.2f\nPercent Read Bases Mapped = %.2f\nAverage Coverage = %.2f\nPercent Unique Start Sites = %.2f" % \
		((mapped/(mapped+unmapped))*100,(refbasescovered/totalrefbases)*100,(totalcov/totalreadbases)*100,totalcov/totalrefbases,(usscount/totalrefbases)*100))
		#stats.write(out)

    
    
    	#for p in bamfile.pileup():
    	#    print p.pos,"\t",p.n
    
    	bamfile.close()

    	print "Done calculating MapStats"

    	return

if __name__ == '__main__':
	main()
