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
		self.probs = {}
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

	def plotData(self,outputType='eps'):
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
			plt.savefig(sampleName+"_probabilities."+outputType)

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
		# loads self.probs[sampleName][gff][methType] array of 10 hist values.
		for sampleName in self.samples:
			cPickle.load(open(sampleName+"_probs.pickle",'rb'))
	def makeGffStats(self):
		gffStats = {}
		for gffFile in self.gff:
			if self.gff[gffFile] != "NC":
				gffStats[gffFile] = simpleStats(gffFile)
		self.transMatrix = simpleAnalyze(gffStats, self.gff, self.fa+'.fai', self.windowSize)

def getChromSizes(inFai):
	chromLens = {}
	for line in open(inFai,'r'):
		#Chr1	30427671	6	79	80
		tmp = line.split('\t')
		chromLens[tmp[0]] = int(tmp[1])
	return chromLens

def simpleStats(inFile):
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

def calcNC(genomeBases, gffStats, gff):
	total = genomeBases
	for gffFile in gff:
			total -= gffStats[gffFile][0]*gffStats[gffFile][1]
	return total

def simpleAnalyze(gffStats, gff, inFai, windowSize):
	chromDict = getChromSizes(inFai)
	genomeBases = sum(chromDict.values())
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
	print map(sum, transMatrix)
	print "\t"+'\t'.join(letterIndex)
	for i in range(len(letterIndex)):
		print letterIndex[i]+'\t'+'\t'.join(map(lambda x: str(round(x,2)), transMatrix[i]))
	return transMatrix
		
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

def regionFinder(E,TP,fa,chromLens,states): #get this integrated
	#TP[from][to]
	# Emissions ##############
	#[ 0.   0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9  1. ]
	#[0 1)
	#("NC","G","TG","TE","PG")

	#get emission probabilities parsed like this
	#pass in state name order
	#E[state][meth][count][index]

	EP = parseEP(E) #EP[state][meth][bin]

	for sampleName in self.samples:
		methIndexes = []
		for mf in self.samples[sampleName]:
			methIndexes.append(MethIndex(mf,fa))
			useVoting(chromLens, experiment, fastaFile, TP, EP, "Chr1", states, 50)
			sys.stderr.write("Starting sample "+samNum+"\n")
			#window(chromLens, experiment, fastaFile, TP, EP, "Chr1", states)

def window(chromLens, experiment, fastaFile, TP, EP, chrom, states):
	starts, methSeq = parseChrom(fastaFile, experiment, chrom, 1, chromLens)
	statePath = viterbi(TP, EP, methSeq, states, fastaFile, chrom)
	printPath(starts, statePath, states, chromLens, chrom)

def useVoting(chromLens, experiment, fastaFile, TP, EP, chrom, states, slideSize):
	voteArray = makeVoteArray(chromLens["Chr1"])
	for start in xrange(1,windowSize,slideSize):
		starts, methSeq = parseChrom(fastaFile, experiment, "Chr1", start, chromLens)
		statePath = viterbi(TP, EP, methSeq, states, fastaFile, "Chr1")
		del methSeq
		vote(voteArray, statePath, starts)
		print "%d/%d" % (start, windowSize)
	#plotVotes(states, chromLens, chrom, voteArray)
	#maxVotes = np.argmax(voteArray, axis=1)
	TEvotes = voteArray[:,3] >= 5
	maxVotes = np.zeros(len(TEvotes))
	maxVotes[TEvotes] = 3
	#print TEs
	printFeatures(maxVotes, 3, states)

def printFeatures(maxVotes, targetIndex, states):
	features = []
	start = 0
	end = 0
	for index in xrange(len(maxVotes)):
		if maxVotes[index] == 3:
			if start == 0:
				start = index+1
			else:
				end = index+1
		else:
			if end != 0:
				features.append((start,end))
				start,end = 0,0
	for i in features:
		print "Chr1\t%s\t%d\t%d" % (states[targetIndex], i[0], i[1])

def plotVotes(states, chromLens, chrom, voteArray):
	plt.figure(0)
	for state in xrange(len(states)):
		X = np.arange(1,chromLens[chrom]+1,dtype="float")
		sums = np.sum(voteArray,axis=1,dtype="float")
		Y = voteArray[:,state]/sums
		Y[0] = 0
		Y[-1] = 0
		plt.fill(X,Y,alpha=0.3,label=states[state])
	plt.legend()
	plt.show()

def vote(voteArray, statePath, starts):
	for index in xrange(len(statePath)):
		state = statePath[index]
		startIndex = starts[index]-1
		endIndex = startIndex+windowSize
		voteArray[startIndex:endIndex,state] += 1

def makeVoteArray(size):
	return np.zeros((size,5), dtype=int)

def printPath(starts, statePath, states, chromLens, chrom):
	for i in xrange(len(statePath)):
		if statePath[i] == 3:
			start = starts[i]
			end = start+windowSize-1
			print "%s\t%d\t%d" % (chrom, start, end)

def viterbi(TP, EP, methSeq, states, fastaFile, chrom):
	#TP[from][to]
	startProb = [0.0, -np.inf, -np.inf, -np.inf, -np.inf]
	TP = np.log2(TP)
	numStates = len(states)
	#print states
	pathMat = np.empty([len(states), len(methSeq)], dtype=int)
	probMat = np.empty([len(states), len(methSeq)])
	pathMat[:,0] = -1
	for i in xrange(numStates):
		probMat[i,0] = startProb[i]+calcEmission(EP[i], methSeq[0])
	for i in xrange(1,len(methSeq)):
		if sum(np.isnan(methSeq[i])) == 3:
				pMax = np.amax(probMat[:,i-1])
				pMaxIndex = probMat[:,i-1].argmax()
				for x in xrange(numStates):
					pathMat[x,i] = pMaxIndex
					probMat[x,i] = -np.inf
				probMat[pMaxIndex,i] = pMax
		else:
			for x in xrange(numStates):
				priors = []
				for y in xrange(numStates):
					priors.append(TP[y][x]+probMat[y,i-1])
				pMax = max(priors)
				pMaxIndex = priors.index(pMax)
				pathMat[x,i]= pMaxIndex
				probMat[x,i] = pMax+calcEmission(EP[x], methSeq[i])
	maxProb = np.amax(probMat[:,-1])
	maxIndex = probMat[:,-1].argmax()
	path = [maxIndex]
	for i in xrange(len(methSeq)-1,0,-1): # backtrace path
		path.append(pathMat[path[-1],i])
	#print "plog: %.3f" % (maxProb)
	path = path[::-1]
	return path

def calcEmission(stateEP, data): #do something with this
	#CG, CHH, CHG
	total = 0.0
	for i in xrange(len(data)):
		if data[i] == 1:
			total += stateEP[i][9]
		elif np.isnan(data[i]):
			pass
		else:
			total += stateEP[i][int(data[i])]
	return total
	
def parseChrom(fastaFile, experiment, chrom, start, chromLens):
	FA = pysam.Fastafile(fastaFile)
	chromLen = chromLens[chrom]
	starts = np.arange(start, chromLen+1, windowSize, dtype=int)
	if starts[-1]+windowSize-1 > chromLen:
		starts = starts[:-1]
	methSeq = range(len(starts))
	for index in xrange(len(starts)):
		start = starts[index]	#1-indexed
		end = start+windowSize-1
		methFreqs = getFreqs(FA, experiment, chrom, start, end)
		methSeq[index] = methFreqs[:]
	return (starts, methSeq)

def getFreqs(FA, experiment, chrom, start, end):
	output = []
	Y = np.zeros(3,end-start+1)	
	C = np.zeros(3,end-start+1, dtype=int)
	for methReplicate in experiment:
		methA, contA = methReplicate.fetch(chrom=chrom, start=start, end=end)
		methNA = np.array(methA)
		contNA = np.array(contA)
		for methIndex in xrange(3):
			methLocs = contNA == methIndex+1 || contNA == methIndex+3
			Y[methIndex,methLocs] += methNA[methLocs]
			C[methIndex,methLocs] += contNA[methLocs]
	for methIndex in xrange(3):
		if sum(C[methIndex]) == 0:
			output.append(np.nan)
		else:
			output.append(np.mean(Y[methIndex,C[methIndex,C[methIndex] > 0]]))
	return tuple(output)
	
def parseEP(E):
	#E[state][meth][count][index]
	return map(parseState, E)
def parseState(EPstate):
	return map(parseMeth, EPstate)
def parseMeth(EPmeth):
	# added a pseudo-count so regions to reduce stringency
	total = float(EPmeth[1]+len(EPmeth[0]))
	return map(lambda x: np.log2((x+1.0)/total), EPmeth[0])
