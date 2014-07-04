#!/usr/bin/python

import numpy as np
import sys
import os
import pysam
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pymethyl import MethIndex
import cPickle

motifIndex = {'+':{'CG':1,'CHG':2,'CHH':3}, '-':{'CG':4,'CHG':5,'CHH':6}}
methTypes = ('CG','CHG','CHH')

class methExperiment:
	def __init__(self):
		self.samples = {}
		self.gff = {}
		self.fa = "Unset"
		self.emissions = {}
		self.windowSize = 500.0
	def addSample(self, name, samples = []):
		self.samples[name] = []
		print "Added sample: %s" % (name)
		for i in samples: self.addReplicate(name, i)
	def addReplicate(self,name,replicate):
		if not checkFiles([replicate]):
			self.samples[name].append(replicate)
			print "Added %s to %s" % (replicate, name)
	def addGFF(self, name, gffFile):
		#psedudogene,chromosome,gene,transposable_element,rest is NC
		if not checkFiles([gffFile]):
			if gffFile not in self.gff:
				self.gff[gffFile] = name
		else:
			return 1
	def setReference(self, fa):
		if not checkFiles([fa]):
			self.fa = fa
			self.FA = pysam.Fastafile(self.fa)
	def printExperiment(self): # prints the experiment layout
		print "\nGFF Files:"
		for i in range(len(self.gff)-1):
			k = self.gff.keys()[i]
			print "  |-  [ %s ]\t %s" % (self.gff[k], k)
		k = self.gff.keys()[-1]
		print "   -  [ %s ]\t %s" % (self.gff[k], k)
		print "\nReference: %s\n" % (self.fa)
		for sample in self.samples:
			print "%s:" % (sample)
			numReps = len(self.samples[sample])
			for i in range(numReps-1):
				print "  |- %s" % (self.samples[sample][i])
			print "   - %s\n" % (self.samples[sample][numReps-1])
	def makeMethBins(self, minRegionLen = 700): # creates a dictinary of methylation frequencies per feature
		print "Making methylation bins"
		self.dataDict = {}
		for sampleName in self.samples:
			self.dataDict[sampleName] = {}
			methIndexes = []
			for mf in self.samples[sampleName]:
				methIndexes.append(MethIndex(mf,self.fa))
			for gffFile in self.gff:
				self.dataDict[sampleName][gffFile] = {}
				for methType in methTypes:
					self.dataDict[sampleName][gffFile][methType] = []
				count = 0
				for struct in gff_gen(gffFile, minRegionLen):
					count += 1
					methylBins = geneBin(struct, self.FA, methIndexes)
					for i in range(3):
						methType = methTypes[i]
						if np.isfinite(methylBins[i]):
							self.dataDict[sampleName][gffFile][methType].append(methylBins[i])
				print "%s: %i" % (gffFile,count)
	def plotMethBins(self,outputType='eps'): # prints the frequencies by type from makeMethBins as boxplots
		for sampleName in self.samples:
			fig = plt.figure(figsize=(11,7))
			count = 1
			for methType in methTypes:
				dataList = []
				names = []
				for gff in self.gff:
					dataList.append(self.dataDict[sampleName][gff][methType])
					names.append(self.gff[gff])
				makeBoxplot("31"+str(count), dataList, methType, names, fig)
				count += 1
			plt.xticks(range(len(names)+1)[1:],names)
			plt.subplots_adjust(hspace=0.1, top=0.95, bottom=0.05, left=0.07, right=0.96)
			plt.suptitle("Methylation Probabilities by Region")
			plt.savefig(sampleName+"_probabilities."+outputType)
	def makeEmissions(self): # generate emission probabilies from makeMethBins
		for sampleName in self.samples:
			self.emissions[sampleName] = {}
			for gff in self.gff:
				self.emissions[sampleName][gff] = {}
				print gff
				for methType in methTypes:
					Data = self.dataDict[sampleName][gff][methType]
					hist, edges = np.histogram(Data, bins=10, range=(0,1)) # edges are half open [0, 1)
					self.emissions[sampleName][gff][methType] = hist
		cPickle.dump(self.emissions, open("emissions.pickle",'wb'))
	def readEmissions(self):
		if not checkFiles(['emissions.pickle']):
			# loads self.emissions[sampleName][gff][methType] array of 10 hist values.
			self.emissions = cPickle.load(open("emissions.pickle",'rb'))
		else:
			print "Make emissions first"
	def printEmissions(self):
		print "\nEmission Probabilities"
		for sampleName in self.samples:
			for gff in self.gff:
				print gff
				for methType in methTypes:
					string = "%s:\t" % (methType)
					hist = self.emissions[sampleName][gff][methType]
					for i in hist:
						string += "%d\t" % (i)
					print string + "sum: %d" % (sum(hist))
		print "[0\t"+")\t".join(map(str,np.arange(0.1,1.1,0.1)))+']'
	def makeTransitions(self): # makes transition probabilities
		gffStats = {}
		print "\nGFF Stats:"
		for gffFile in self.gff:
			if self.gff[gffFile] != "NC":
				gffStats[gffFile] = simpleStats(gffFile)
		self.transMatrix = makeTransMatrix(gffStats, self.gff, self.fa+'.fai', self.windowSize)

def getChromSizes(inFai): # returns a dicationary of chromosome sizes
	chromLens = {}
	for line in open(inFai,'r'):	#Chr1	30427671	6	79	80
		tmp = line.split('\t')
		chromLens[tmp[0]] = int(tmp[1])
	return chromLens

def simpleStats(inFile): # returns the count and average size of a feature type
	print inFile
	count = 0
	total = 0
	for line in open(inFile,'r'):
		#Chr1	TAIR10	transposable_element	11897	11976
		count += 1
		tmp = line.split('\t')
		start = int(tmp[3])
		end = int(tmp[4])
		total += end-start+1
	print "\tOccurrences: %d" % (count)
	print "\tAverage Size: %.1f" % (total/float(count))
	return (count, total/float(count))

def makeTransMatrix(gffStats, gff, inFai, windowSize): #generate transition probabilities
	chromSizes = getChromSizes(inFai)
	genomeBases = sum(chromSizes.values())
	#self.gff[gffFile] = name
	#NC	TG	G	PG	TE
	letterIndex = gff.values()
	NC, G, TE = map(letterIndex.index, ["NC","G","TE"])
	transMatrix = np.zeros((5,5))
	basesNC = genomeBases-np.sum(map(lambda x: x[0]*x[1], gffStats.values()))
	windowsNC = basesNC/windowSize
	regionCounts = sum(map(lambda x: x[0], gffStats.values()))
	transMatrix[NC][NC] = windowsNC/float(windowsNC+regionCounts)
	te2gene = 4178.0/31189.0 #make this dynamic
	gene2te = 4178.0/28775.0 #make this dynamic
	for gffFile,v in gffStats.iteritems():
		curLetter = gff[gffFile]
		curIndex = letterIndex.index(curLetter)
		count, avg = v
		windowsPerRegion = avg/windowSize
		toRegion = count/float(windowsNC+regionCounts)
		transMatrix[NC][curIndex] = toRegion
		if curLetter == "G":
			stayProb = 1.0/windowsPerRegion
			leaveProb = 1-(gene2te+stayProb)
			transMatrix[curIndex][TE] = gene2te
		elif curLetter == "TE":
			stayProb = 1.0/windowsPerRegion
			leaveProb = 1-(te2gene+stayProb)
			transMatrix[curIndex][G] = te2gene
		else:
			stayProb = 1.0/windowsPerRegion
			leaveProb = (1.0-stayProb)
		transMatrix[curIndex][NC] = leaveProb
		transMatrix[curIndex][curIndex] = stayProb
	print "\nTransition Matrix"
	print "\t"+'\t'.join(letterIndex)+'\tSum'
	for i in range(len(letterIndex)):
		print letterIndex[i]+'\t'+'\t'.join(map(lambda x: str(round(x,2)), transMatrix[i]))+'\t'+str(round(sum(transMatrix[i]),2))
	return transMatrix
		
def checkFiles(fileList): # check to see if any files in a list don't exist
	if not fileList: return 0
	ret = 0
	for f in fileList:
		if not os.path.isfile(f):
			ret = 1
			print "Could not find: %s" % (f)
	return ret

def geneBin(geneStruct, FA, methReps):
	gChr, gStrand, gStart, gEnd, gID = geneStruct
	geneLen = gEnd-gStart+1
	Y = np.zeros((3,geneLen))			#methyl values
	C = np.zeros((3,geneLen))			#count values
	for rep in methReps:
		methA, contA = rep.fetch(chrom=gChr,start=gStart,end=gEnd)
		methNA = np.array(methA)
		contNA = np.array(contA) #make np for slicing
		rNot1 = methNA != -1
		for methIndex in range(1,len(methTypes)+1): #1-indexed
			for i in range(0,4,3): #handle both strands
				methUse = contNA == methIndex+i
				locUse = np.logical_and(rNot1, methUse)
				C[methIndex-1, locUse] += 1
				Y[methIndex-1, locUse] += methNA[locUse]
			#if gStrand == '+':
			#	methUse = contNA == methIndex+3
			#else:
			#	methUse = contNA == methIndex
			#locUse = np.logical_and(rNot1, methUse)
			#C[methIndex-1, locUse] += 1
			#Y[methIndex-1, locUse] += methNA[locUse]
	outBins = []
	for i in range(3):
		cgt1 = C[i,:] > 1
		Y[i,cgt1] = Y[i,cgt1]/C[i,cgt1] #average across replicates
		outBins.append(meanFunc(Y[i,:],C[i,:]))
	return outBins #(CG, CHG, CHH)

def meanFunc(Y,C): # returns the means of bins
	gZero = Y[C > 0]
	if len(gZero): return np.mean(gZero)
	return np.nan

def gff_gen(inFile, minSize):
	#Chr1    TAIR10  transposable_element_gene       433031  433819  .       -       .       ID=AT1G02228;Note=transposable_element_gene;Name=AT1G02228;Derives_from=AT1TE01405
	IF = open(inFile,'r')
	for line in IF:
		tmp = line.split('\t')
		chrom = tmp[0]
		strand = tmp[6]
		start, end = map(int, tmp[3:5])
		idNum = tmp[8].split(';')[0]
		idNum = idNum.split('=')[1]
		if end-start > minSize:
			yield (chrom, strand, start, end, idNum)
	IF.close()

def makeBoxplot(sp, dataList, methType, names, fig):
	ax = fig.add_subplot(int(sp))
	ax.boxplot(dataList, sym='')
	ax.set_ylabel(methType)
	if sp != "313":
		plt.setp(ax.get_xticklabels(), visible=False)
		plt.setp(ax.get_xticklines(), visible=False)
