#!/usr/bin/python3

#Questions: Are seq 1 and seq 2 the same length?

from which_pyqt import PYQT_VER
if PYQT_VER == 'PYQT5':
	from PyQt5.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT4':
	from PyQt4.QtCore import QLineF, QPointF
elif PYQT_VER == 'PYQT6':
	from PyQt6.QtCore import QLineF, QPointF
else:
	raise Exception('Unsupported Version of PyQt: {}'.format(PYQT_VER))

import random
from collections import namedtuple
from enum import Enum

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
# NOTMATCH = 0
INDEL = 5
SUB = 1

DEBUG = False

class direction(Enum):
	LEFT = 1
	TOP = 2
	DIAGONAL = 3
	NONE = 4

class GeneSequencing:

	def __init__( self ):
		pass
	def printDict(self, scoreDict):
		for key in scoreDict:
			print(key, scoreDict[key])
	def compareC(self, char1, char2, direction):
		if char1 == char2 and direction == direction.DIAGONAL:
			return MATCH # Most uptimal case where the characters are the same and we are comparing diagonally
		elif (direction == direction.DIAGONAL):
			return SUB # If the characters are different but still comparing diagonally
		else: #If we are comparing left or top, so its always an indel
			return INDEL
	def is_valid(self, rowIndex, colIndex):
		return abs(rowIndex - colIndex) <= MAXINDELS
	
	def start_index(self, rowIndex):
		return max(0, rowIndex - MAXINDELS)
	def end_index(self, rowIndex, length):
		return min(length, rowIndex + MAXINDELS + 1) #Plus one because of indexing, and inclusivity	
	def backtrack(self, seq1, seq2, directionDict):
		alignment1 = []
		alignment2 = []
		rowIndex, colIndex = len(seq1), len(seq2)

		while rowIndex > 0 and colIndex > 0:
			currentDirection = directionDict[(rowIndex, colIndex)]

			if currentDirection == direction.DIAGONAL:
				alignment1.append(seq1[rowIndex-1])
				alignment2.append(seq2[colIndex-1])
				rowIndex -= 1
				colIndex -= 1
			elif currentDirection == direction.LEFT:
				alignment1.append(seq1[rowIndex-1])
				alignment2.append('-')
				rowIndex -= 1
			elif currentDirection == direction.TOP:
				alignment1.append('-')
				alignment2.append(seq2[colIndex-1])
				colIndex -= 1
			elif currentDirection == direction.NONE:
				alignment1.append(seq1[rowIndex-1])
				alignment2.append(seq2[colIndex-1])
				rowIndex -= 1
				colIndex -= 1

		# Reverse the alignments and join them into strings
		alignment1 = ''.join(reversed(alignment1))
		alignment2 = ''.join(reversed(alignment2))

		return alignment1, alignment2

# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
# how many base pairs to use in computing the alignment

	def align(self, seq1: str, seq2: str, banded, align_length):
		self.banded = banded
		self.MaxCharactersToAlign = MAXINDELS
		
		seq1 = seq1[:align_length]
		seq2 = seq2[:align_length]

		#Check to see if the sequences are the same first
		# if(seq1 == seq2):
		# 	cost = -3 * len(seq1)
		# 	return {'align_cost': cost, 'seqi_first100':'', 'seqj_first100':''}

		#Check to see if answer is in bound aka valid
		if (banded and abs(len(seq1) - len(seq2)) > MAXINDELS):
			return {'align_cost': float('inf'), 'Not in banded segment':'', 'Not in banded segment':''}
			
		#MY_STORAGE section
		scoreDict = {} #Make a new dictionary that stores the score of each cell from the grid values
		directionDict = {} #Make a new dictionary that stores the direction of each cell from the grid values
		# ROW FIRST, the COLUMN 

		# Initialize the zeroth row and column
		for i in range(len(seq1) + 1):
			scoreDict[(i, 0)] = i * INDEL
			directionDict[(i, 0)] = direction.TOP

		for j in range(len(seq2) + 1):
			scoreDict[(0, j)] = j * INDEL
			directionDict[(0, j)] = direction.LEFT

		myScore : int = self.compareC(seq1[0], seq2[0], direction.DIAGONAL) 
		scoreDict[(1,1)] = myScore
		directionDict[(1,1)] = direction.NONE

		if (banded):
			# Start on the first row and run down and fill it out, then do the next row
			for rowIndex in range(len(seq1)):
				for colIndex in range(self.start_index(rowIndex), self.end_index(rowIndex, len(seq2))):
					# if DEBUG:
					# 	print("i: ", rowIndex, "j: ", colIndex)
					bestScore : float = float('inf')
					if (self.is_valid(rowIndex-1,colIndex)): #We can try pulling from the left
						compareScore : int = self.compareC(seq1[rowIndex], seq2[rowIndex], direction.LEFT) #Compare the two characters
						finalScore = scoreDict[(rowIndex,colIndex+1)] + compareScore
						if (finalScore < bestScore): #If this is the best score, store it
							bestScore = finalScore
							scoreDict[(rowIndex+1,colIndex+1)] = finalScore
							directionDict[(rowIndex+1,colIndex+1)] = direction.LEFT

					if (self.is_valid(rowIndex,colIndex -1)): #We can try pulling from the top
						compareScore : int = self.compareC(seq1[rowIndex], seq2[rowIndex], direction.TOP)
						finalScore = scoreDict[(rowIndex+1,colIndex)] + compareScore
						if (finalScore < bestScore): #If this is the best score, store it
							bestScore = finalScore
							scoreDict[(rowIndex+1,colIndex+1)] = finalScore
							directionDict[(rowIndex+1,colIndex+1)] = direction.TOP

					#We can try pulling from the diagonal
					compareScore : int = self.compareC(seq1[rowIndex], seq2[rowIndex], direction.DIAGONAL)
					finalScore = scoreDict[(rowIndex,colIndex)] + compareScore
					if (finalScore < bestScore):
						bestScore = finalScore
						scoreDict[(rowIndex+1,colIndex+1)] = finalScore
						directionDict[(rowIndex+1,colIndex+1)] = direction.DIAGONAL
		else: #Not banded
			# Start on the first row and run down and fill it out, then do the next row
			for rowIndex in range(len(seq1)):
				for colIndex in range(len(seq2)):
					# if DEBUG:
					# 	print("i: ", rowIndex, "j: ", colIndex)
					bestScore : float = float('inf')
					#We can try pulling from the left
					compareScore : int = self.compareC(seq1[rowIndex], seq2[rowIndex], direction.LEFT) #Compare the two characters
					finalScore = scoreDict[(rowIndex,colIndex+1)] + compareScore
					if (finalScore < bestScore): #If this is the best score, store it
						bestScore = finalScore
						scoreDict[(rowIndex+1,colIndex+1)] = finalScore
						directionDict[(rowIndex+1,colIndex+1)] = direction.LEFT

					#We can try pulling from the top
					compareScore : int = self.compareC(seq1[rowIndex], seq2[rowIndex], direction.TOP)
					finalScore = scoreDict[(rowIndex+1,colIndex)] + compareScore
					if (finalScore < bestScore): #If this is the best score, store it
						bestScore = finalScore
						scoreDict[(rowIndex+1,colIndex+1)] = finalScore
						directionDict[(rowIndex+1,colIndex+1)] = direction.TOP

					#We can try pulling from the diagonal
					compareScore : int = self.compareC(seq1[rowIndex], seq2[rowIndex], direction.DIAGONAL)
					finalScore = scoreDict[(rowIndex,colIndex)] + compareScore
					if (finalScore < bestScore):
						bestScore = finalScore
						scoreDict[(rowIndex+1,colIndex+1)] = finalScore
						directionDict[(rowIndex+1,colIndex+1)] = direction.DIAGONAL

		
		score:int = scoreDict[(len(seq1),len(seq2))]

		returnVal = self.backtrack(seq1, seq2, directionDict)
		alignment1 = returnVal[0][:100]
		alignment2 = returnVal[1][:100]

		# alignment1 = 'abc-easy  DEBUG:({} chars,align_len={}{})'.format(
		# 	len(seq1), align_length, ',BANDED' if banded else '')
		# alignment2 = 'as-123--  DEBUG:({} chars,align_len={}{})'.format(
		# 	len(seq2), align_length, ',BANDED' if banded else '')
###################################################################################################

		return {'align_cost':score, 'seqi_first100':alignment1, 'seqj_first100':alignment2}
