import numpy as np
from team import getChromSizes
from team import methTypes
from pymethyl import MethIndex 
import sys
import pysam

windowSize = 500

def regionFinder(E,TP,fa,gff,samples): #get this integrated
	states = gff.values()
	chromSizes = getChromSizes(fa+'.fai')
	#TP[from][to]
	# Emissions ##############
	#[ 0.   0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9  1. ]
	#[0 1)
	#("NC","G","TG","TE","PG")

	EP = parseEP(E,gff) #EP[sample][state][meth][index]

	for sampleName in samples:
		methFiles = []
		for mf in samples[sampleName]:
			methFiles.append(MethIndex(mf,fa))
		sys.stderr.write("Starting sample "+sampleName+"\n")
		useVoting(chromSizes, methFiles, fa, TP, EP[sampleName], "Chr1", states, 50)
		#window(chromLens, methFiles, fastaFile, TP, EP, "Chr1", states)

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
	teIndex = states.index("TE")
	for i in xrange(len(statePath)):
		if statePath[i] == teIndex:
			start = starts[i]
			end = start+windowSize-1
			print "%s\t%d\t%d" % (chrom, start, end)

def viterbi(TP, EP, methSeq, states, fastaFile, chrom):
	stateNC = states.index("NC")
	#TP[from][to]
	startProb = [-np.inf, -np.inf, -np.inf, -np.inf, -np.inf]
	startProb[stateNC] = 1.0
	TP = np.log2(TP)
	numStates = len(states)
	#print states
	pathMat = np.empty([len(states), len(methSeq)], dtype=int)
	probMat = np.empty([len(states), len(methSeq)])
	pathMat[:,0] = -1
	for i in xrange(numStates):
		probMat[i,0] = startProb[i]+calcEmission(EP[i], methSeq[0])
	for i in xrange(1,len(methSeq)):
		if sum(np.isnan(methSeq[i])) == 3: # no methylation. stay in current state
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
	Y = np.zeros((3,end-start+1)	)
	C = np.zeros((3,end-start+1), dtype=int)
	for methReplicate in experiment:
		methA, contA = methReplicate.fetch(chrom=chrom, start=start, end=end)
		methNA = np.array(methA)
		contNA = np.array(contA)
		for methIndex in xrange(3):
			methLocs = np.logical_or(contNA == methIndex+1, contNA == methIndex+3)
			Y[methIndex,methLocs] += methNA[methLocs]
			C[methIndex,methLocs] += contNA[methLocs]
	for methIndex in xrange(3):
		if sum(C[methIndex]) == 0:
			output.append(np.nan)
		else:
			output.append(np.mean(Y[methIndex,C[methIndex,C[methIndex] > 0]]))
	return tuple(output)
	
def parseEP(E, gff, pseudoCount=1):
	#self.emissions[sampleName][gff][methType] array of 10 hist values.
	#EP[sample][state][meth][index]
	EP = {}
	for sample in E:
		EP[sample] = []
		for gffFile in gff:
			methData = []
			for methType in methTypes:
				total = sum(E[sample][gffFile][methType])
				methData.append(map(lambda x: np.log2((x+float(pseudoCount))/float(total)), E[sample][gffFile][methType]))
			EP[sample].append(methData)
	return EP
