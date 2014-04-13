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

#H Dictionaries
hD = {'A':'H', 'T':'H', 'C':'H', 'G':'G', 'N':'N', 'R':'N', 'Y':'N', 'K':'N', 'M':'N', 'S':'N', 'W':'N', 'B':'N', 'D':'N', 'H':'N', 'V':'N', 'X':'N'}
compDict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N', 'R':'N', 'Y':'N', 'K':'N', 'M':'N', 'S':'N', 'W':'N', 'B':'N', 'D':'N', 'H':'N', 'V':'N', 'X':'N'}
motifIndex = {'+':{'CG':1,'CHG':2,'CHH':3}, '-':{'CG':4,'CHG':5,'CHH':6}}
methTypes = ('CG','CHG','CHH')

class methExperiment:
	def __init__(self):
		self.samples = {}
		self.gff = {}
		self.fa = "Unset"
		self.probs = {}
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
	def printExperiment(self):
		print "\nGFF File: %s\nReference File: %s\n" %(self.gff,self.fa)
		for sample in self.samples:
			print "%s:" % (sample)
			numReps = len(self.samples[sample])
			for i in range(numReps-1):
				print "  |- %s" % (self.samples[sample][i])
			print "   - %s" % (self.samples[sample][numReps-1])
	def makeMethBins(self):
		self.dataDict = {}
		minRegionLen = 700
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

	def plotData(self):
		for sampleName in self.samples:
			fig = plt.figure(figsize=(11,7))
			count = 1
			for methType in ("CG","CHH","CHG"):
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
			plt.savefig(sampleName+"_probabilities.eps")
			plt.savefig(sampleName+"_probabilities.svg")

	def printProbs(self):
		for sampleName in self.samples:
			self.probs[sampleName] = {}
			for gff in self.gff:
				self.probs[sampleName][gff] = {}
				print gff
				for methType in methTypes:
					Data = self.dataDict[sampleName][gff][methType]
					# edges are half open [0, 1)
					hist, edges = np.histogram(Data, bins=10, range=(0,1))
					self.probs[sampleName][gff][methType] = hist
					string = "%s:\t" % (methType)
					for i in hist:
						string += "%d\t" % (i)
					string += "sum: %d" % (sum(hist))
					print string
			print edges

	def writeProbs(self):
		for sampleName in self.samples:
			cPickle.dump(self.probs, open(sampleName+"_probs.pickle",'wb'))
	def readProbs(self):
		for sampleName in self.samples:
			cPickle.load(open(sampleName+"_probs.pickle",'rb'))
		
def checkFiles(fileList):
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
	outBins = []
	for i in range(3):
		cgt1 = C[i,:] > 1
		Y[i,cgt1] = Y[i,cgt1]/C[i,cgt1] #average across replicates
		outBins.append(meanFunc(Y[i,:],C[i,:]))
	return outBins #(CG, CHG, CHH)

def meanFunc(Y,C): #getting weird output values
	gZero = Y[C > 0]
	if len(gZero):
		return np.mean(gZero)
	else:
		return np.nan

def gff_gen(inFile, minSize):
	#Chr1    TAIR10  transposable_element_gene       433031  433819  .       -       .       ID=AT1G02228;Note=transposable_element_gene;Name=AT1G02228;Derives_from=AT1TE01405
	IF = open(inFile,'r')
	for line in IF:
		tmp = line.split('\t')
		chrom = tmp[0]
		strand = tmp[6]
		start = int(tmp[3])
		end = int(tmp[4])
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

if __name__ == "__main__":
	main()
