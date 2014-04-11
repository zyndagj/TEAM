#!/usr/bin/python

import sys
import os.path
import struct
from cFetch import cFetch

class MethIndex:
	"""
	Creates and queries methylation index for quick random access.
	"""

	def __init__(self, methFile, faFile):
		"""
		Index constructor
		
		Arguments
		=================================
		methFile		- Methylation input file from MethylCoder
		faFile			- Fasta file of reference (index required)
		"""
		self.methFile = methFile
		self.methBin = methFile + '.bin'
		self.methBinIndex = self.methBin + '.idx'
		self.FA = faFile
		self.FAI = self.FA + '.fai'
		self.chromDict = {}
		self.seekDict = {}
		self.readFAI()
		if not os.path.exists(self.methBin) or not os.path.exists(self.methBinIndex):
			self.makeIndex()
		else:
			self.readBinIndex()

	def fetch(self, chrom, start = -1, end = -1, minCov = 5, maxCov = 200):
		"""
		Fetches a region of methylation ratios.
		
		If no start is given, 1-end is returned.
		If no end is given, start- is returned.

		If depth < minCov or depth > maxCov, then -1 returned
		at that position.
		
		Agurments
		=================================
		chrom   - Chromosome
		start   - Start of region (1-indexed)
		end	- End of region (1-indexed)
		minCov	- Minimum coverage needed (Default: 5)
		maxCov	- Maximum coverage allowed (Default: 200)
		"""
		#think about
		# http://stackoverflow.com/questions/8461798/how-can-i-struct-unpack-many-numbers-at-once
		if chrom not in self.seekDict:
			print 'Not a real chromosome'
			return []
		if end < start and end != -1:
			print 'Bad coordinates'
			return []
		if start == -1:
			seekStart = self.seekDict[chrom]
			if end == -1:
				countEnd = self.chromDict[chrom]
			else:
				countEnd = end
		else:
			seekStart = self.seekDict[chrom] + (start - 1) * 6
			if end == -1:
				countEnd = self.chromDict[chrom]-start+1
			else:
				countEnd = min(end-start+1,self.chromDict[chrom]-start+1)
		# returns (outMeth, outContext)
		return cFetch(self.methBin, seekStart, countEnd, minCov, maxCov)

	def readFAI(self):
		"""
		Reads the sizes of each chromosome from the FASTA index
		"""
		#FAI Format  http://www.biostars.org/p/1495/
                #chrName chrLen chrSeek lineBases lineLen
                #Chr1    30427671        6       79      80
                #Line len is bases+\n
		if not os.path.exists(self.FAI):
			sys.exit('Please index the fasta file')
		IFAI = open(self.FAI, 'r')
		for line in IFAI:
			tmp = line.rstrip('\n').split()
			chrom = tmp[0]
			chromLen = int(tmp[1])
			self.chromDict[chrom] = chromLen

		IFAI.close()

	def makeIndex(self):
		"""
		Makes the .bin and .bin.idx methylation index files
		"""
		#chr     pos     strand  context ratio   eff_CT_count    C_count CT_count        rev_G_count     rev_GA_count    CI_lower   CI_upper
                #Chr1    15      +       AACCC   1.000   1.00    1       1       8       8       0.207   1.000
		IM = open(self.methFile, 'r')
		OI = open(self.methBin, 'wb')
		curPos = 1
		line = IM.readline()
		order = []
		### Grab initial chromsome name
		if line[:3] != 'chr':
			curChrom = line.split()[0]
			IM.seek(0)
		else:
			curChrom = IM.readline().split()[0]
			IM.seek(0)
			IM.readline()
		order.append(curChrom)
		### Parse data
		for line in IM:
			tmp = line.split('\t')
			chrom = tmp[0]
			pos = int(tmp[1])
			contextNum = self.parseContext(tmp[3])
			C = int(tmp[6])
			CorT = int(tmp[7])
			if pos < curPos:
				self.fillChrom(curChrom, curPos, OI)
				curPos = 1
				curChrom = chrom
				order.append(curChrom)
			while curPos < pos:
				self.writeBlank(OI)
				curPos += 1
			self.writeData(OI, C, CorT, contextNum)
			curPos += 1
		self.fillChrom(curChrom, curPos, OI)
		IM.close()
		OI.close()
		self.makeBinIndex(order)

	def getInitChrom(IM, line):
		if line[:3] != 'chr':
			curChrom = line.split()[0]
			IM.seek(0)
		else:
			curChrom = IM.readline().split()[0]
			IM.seek(0)
			IM.readline()
		return curChrom
		
	def parseContext(self, context):
		"""
		Parses the methylation motif

		 1	2	3		 4	5	6
		+CG	CHG	CHH		-CG	CHG	CHH
		"""
		if context[2] == "C":
			if context[3] == "G":
				return 1
			elif context[4] == "G":
				return 2
			else:
				return 3
		else:
			if context[1] == "C":
				return 4
			elif context[0] == "C":
				return 5
			else:
				return 6
	
	def fillChrom(self, curChrom, curPos, F):
		"""
		Writes data for rest of chromosome
		"""
		curLim = self.chromDict[curChrom]
		for i in xrange(curLim - curPos + 1):
			self.writeBlank(F)

	def writeBlank(self, F):
		"""
		Writes two blank values (65535) and one 0 to .bin file.
		65535 is the largest USHORT.
		"""
		F.write('\xff\xff\xff\xff\x00\x00')

	def writeData(self, F, C, CorT, contextNum):
		"""
		Writes three unsigned shorts to .bin file.
		"""
		F.write(struct.pack('HHH', C, CorT, contextNum))

	def makeBinIndex(self, chromOrder):
		"""
		Makes the bin index based on the order the chromosomes
		were written.
		"""
		location = 0
		OBI = open(self.methBinIndex, 'w')
		for chrom in chromOrder:
			size = self.chromDict[chrom]
			OBI.write(chrom + '\t')
			OBI.write(str(location) + '\n')
			self.seekDict[chrom] = location
			location = location + size * 6
		OBI.close()

	def readBinIndex(self):
		"""
		Loads the bin index into a seek dictionary.
		"""
		IF = open(self.methBinIndex, 'r')
		for line in IF:
			tmp = line.rstrip('\n').split('\t')
			self.seekDict[tmp[0]] = int(tmp[1])

		IF.close()
