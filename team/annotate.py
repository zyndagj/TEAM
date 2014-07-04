import numpy as np
from team import getChromSizes
from team import methTypes
from pymethyl import MethIndex 
import sys
import pysam

#firstBases = 1000000

def regionFinder(E,TP,fa,gff,samples, windowSize, slideSize): #get this integrated
	states = gff.values()
	chromSizes = getChromSizes(fa+'.fai')
	#TP[from][to]
	# Emissions ##############
	#[ 0.   0.1  0.2  0.3  0.4  0.5  0.6  0.7  0.8  0.9  1. ]
	#[0 1)
	#("NC","G","TG","TE","PG")

	EP = parseEP(E,gff) #EP[sample][state][meth][index]
	TP = np.log2(TP)

	for sampleName in samples:
		methFiles = []
		for mf in samples[sampleName]:
			methFiles.append(MethIndex(mf,fa))
		sys.stderr.write("Starting sample "+sampleName+"\n")
		for chrom in chromSizes.values():
			useVoting(chromSizes, methFiles, TP, EP, chrom, states, windowSize, slideSize, "TE", sampleName)
			#window(chromLens, methFiles, fastaFile, TP, EP, "Chr1", states)

#def window(chromLens, experiment, fastaFile, TP, EP, chrom, states):
#	starts, methSeq = parseChrom(experiment, chrom, 1, chromLens, windowSize)
#	statePath = viterbi(TP, EP, methSeq, states, fastaFile, chrom)
#	printPath(starts, statePath, states, chromLens, chrom, windowSize)

def useVoting(chromLens, experiment, TP, EP, chrom, states, windowSize, slideSize, targetState, sampleName, voteThresh=5):
	#voteArray = makeVoteArray(firstBases)
	voteArray = makeVoteArray(chromLens[chrom])
	for start in xrange(1,windowSize,slideSize):
		starts, methSeq = parseChrom(experiment, chrom, start, chromLens, windowSize)
		statePath = viterbi(TP, EP, methSeq, states, chrom)
		vote(voteArray, statePath, starts, windowSize)
		print "Finished %d/%d" % (start, windowSize)
		del methSeq
		del starts
		del statePath
	#plotVotes(states, chromLens, chrom, voteArray)
	maxVotes = np.argmax(voteArray, axis=1)
	writeFeatures(maxVotes == states.index(targetState), targetState, chrom, sampleName) #print TEs
	#targetVotes = voteArray[:,states.index(targetState)] >= voteThresh
	#printFeatures(targetVotes, targetState, chrom) #print TEs

def writeFeatures(targetVotes, targetState, chrom, sampleName):
	features = []
	start = 0
	end = 0
	for index in xrange(len(targetVotes)):
		if targetVotes[index]:
			if start == 0:
				start = index+1
			else:
				end = index+1
		else:
			if end != 0:
				features.append((start,end))
				start,end = 0,0
	OF = open(sampleName+'.gff','w')
	fCount = 0
	for i in features:
		OF.write("%s\tTEAM\t%s\t%i\t%i\t.\t.\t.\tID=%s_%i\n" % (chrom, targetState, i[0], i[1], sampleName, fCount))
		fCount += 1
	OF.close()

def plotVotes(states, chromLens, chrom, voteArray):
	plt.figure(0)
	for state in xrange(len(states)):
		X = np.arange(1,chromLens[chrom]+1,dtype=float)
		sums = np.sum(voteArray,axis=1,dtype=float)
		Y = voteArray[:,state]/sums
		Y[0] = 0
		Y[-1] = 0
		plt.fill(X,Y,alpha=0.3,label=states[state])
	plt.legend()
	plt.show()

def vote(voteArray, statePath, starts, windowSize):
	for index in xrange(len(statePath)):
		startIndex = starts[index]-1
		endIndex = startIndex+windowSize
		voteArray[startIndex:endIndex,statePath[index]] += 1

def makeVoteArray(size):
	return np.zeros((size,5), dtype=int)

def printPath(starts, statePath, states, chromLens, chrom, windowSize):
	teIndex = states.index("TE")
	for i in xrange(len(statePath)):
		if statePath[i] == teIndex:
			start = starts[i]
			end = start+windowSize-1
			print "%s\t%d\t%d" % (chrom, start, end)

def viterbi(TP, EP, methSeq, states, chrom):
	stateNC = states.index("NC")
	#TP[from][to]
	startProb = [-np.inf, -np.inf, -np.inf, -np.inf, -np.inf]
	startProb[stateNC] = 0.0
	numStates = len(states)
	pathMat = np.empty([numStates, len(methSeq)], dtype=int)
	probMat = np.empty([numStates, len(methSeq)], dtype=float)
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
				del priors
	maxProb = np.amax(probMat[:,-1])
	maxIndex = probMat[:,-1].argmax()
	path = [maxIndex]
	for i in xrange(len(methSeq)-1,0,-1): # backtrace path
		path.append(pathMat[path[-1],i])
	#print "plog: %.3f" % (maxProb)
	del pathMat
	del probMat
	retPath = path[::-1]
	del path
	return retPath

def calcEmission(stateEP, data): #do something with this
	total = 0.0
	for i in xrange(len(data)):
		if data[i] == 1:
			total += stateEP[i][9]
		elif np.isnan(data[i]):
			pass
		else:
			total += stateEP[i][int(data[i]*10.0)]
	return total

def parseChrom(experiment, chrom, start, chromLens, windowSize):
	#chromLen = firstBases
	chromLen = chromLens[chrom]
	starts = np.arange(start, chromLen+1, windowSize, dtype=int)
	if starts[-1]+windowSize-1 > chromLen:
		starts = starts[:-1]
	bigY, bigC = getFreqs(experiment, chrom, start, chromLen)
	methSeq = range(len(starts))#methSeq = np.empty((len(starts),3), dtype=float)
	cgt1 = bigC > 1
	bigY[cgt1] = bigY[cgt1]/bigC[cgt1] #average y's
	output = range(3)
	for index in xrange(len(starts)):
		sIndex = index*windowSize
		eIndex = sIndex+windowSize
		for methIndex in xrange(3):
			if sum(bigC[methIndex,sIndex:eIndex]) == 0:
				output[methIndex]=np.nan#methSeq[index, methIndex] = np.nan
			else:
				cg0 = bigC[methIndex,sIndex:eIndex] > 0
				#methSeq[index, methIndex] = np.mean(bigY[methIndex,sIndex:eIndex][cg0])
				output[methIndex]=np.mean(bigY[methIndex,sIndex:eIndex][cg0])
				del cg0
		methSeq[index] = tuple(output)
	del bigY
	del bigC
	retObj = (tuple(starts), tuple(methSeq))
	del starts
	del methSeq
	return retObj

def getFreqs(experiment, chrom, start, end):
	Y = np.zeros((3,end-start+1), dtype=float)
	C = np.zeros((3,end-start+1), dtype=int)
	for methReplicate in experiment:
		methA, contA = methReplicate.fetch(chrom=chrom, start=start, end=end)
		methNA = np.array(methA, dtype=float)
		del methA
		contNA = np.array(contA, dtype=int)
		del contA
		for methIndex in xrange(3):
			methLocs = np.logical_or(contNA == methIndex+1, contNA == methIndex+4)
			Y[methIndex,methLocs] += methNA[methLocs]
			C[methIndex,methLocs] += 1
			del methLocs
		del methNA
		del contNA
	retObj = (Y, C)
	del Y
	del C
	return retObj

def parseEP(E, gff, pseudoCount=1):
	#self.emissions[sampleName][gff][methType] array of 10 hist values.
	#EP[sample][state][meth][index]
	EP = []
	for sample in E:
		for gffFile in gff:
			methData = []
			for methType in methTypes:
				total = sum(E[sample][gffFile][methType])
				methData.append(map(lambda x: np.log2((x+float(pseudoCount))/float(total)), E[sample][gffFile][methType]))
			EP.append(methData)
	return EP
