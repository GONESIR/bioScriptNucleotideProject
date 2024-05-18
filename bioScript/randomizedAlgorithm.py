#!/usr/bin/env python3
from collections import Counter
import argparse
import re
import random
parser = argparse.ArgumentParser()

required = parser.add_argument_group('Required arguments')
required.add_argument('-i', '--input', type=str, help='Path to input file', required=True)
args = parser.parse_args()
file_in = open(args.input)
lines = file_in.readlines()

data = []

for line in lines:
    data.append(line.rstrip("\n"))

intValues = data[0].split(" ")

def randomizedMotifSearch(dna, k, t, max_iterations=1000):
    best_motifs = []
    for i in range(t):
        best_motifs.append(dna[i][:k])
    best_score = score(best_motifs)
    iteration = 0
    while iteration < max_iterations:
        motifs = runningRandomizedMotifs(dna, k, t)
        current_score = score(motifs)
        if current_score < best_score:
            best_motifs = motifs
            best_score = current_score    
        iteration += 1
    return best_motifs


def runningRandomizedMotifs(dna,k,t):
    best_motifs = randomSelection(dna,k,t)
    while True:
        profile = matrixGenerator(best_motifs)
        motifs = createMotifs(profile, dna)
        if score(motifs) < score(best_motifs):
            best_motifs = motifs
        else:
            return best_motifs

def randomSelection(dna,k,t):
    motifs = []
    for i in dna:
        point = random.randint(0, ((len(i)-k)))
        motif = i[point:point + k]
        motifs.append(motif)
    return motifs

def matrixGenerator(motifs):
    profile =[]
    for i in range(len(motifs[0])): # Laplace Rule
        baseA = 1 
        baseC = 1 
        baseG = 1
        baseT = 1
        for motif in motifs:
            if motif[i] == 'A':
                baseA += 1
            elif motif[i] == 'C':
                baseC += 1
            elif motif[i] == 'G':
                baseG += 1
            elif motif[i] == 'T':
                baseT += 1
        profile.append([baseA / (len(motifs) + 4), baseC/ (len(motifs) + 4), \
                        baseG / (len(motifs) + 4), \
                        baseT / (len(motifs) + 4)])  # Laplace Rule                       
    return profile

def score(motifs):
    consensus = findConsesus(motifs)
    score = 0
    for motif in motifs:
        score += hemmingDistance(consensus, motif)
    return score

def findConsesus(motifs):
    consensus = ""
    for i in range(len(motifs[0])):
        baseCounts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        for motif in motifs:
            baseCounts[motif[i]] += 1
        consensusBase = max(baseCounts, key=baseCounts.get)
        consensus += consensusBase
    return consensus

def createMotifs(profile, dna):
    motifs = []
    for i in dna:
        motifs.append(profileMostProbableKmer(i, profile, len(profile)))
    return motifs

def profileMostProbableKmer(dna_string, profile, k):
    bestPattern = dna_string[0:0 + k]
    bestProbability = 0
    for i in range(len(dna_string) - k + 1):
        stringKmer = dna_string[i:i + k]
        newProbablity = calculateProbability(profile, stringKmer)
        if newProbablity > bestProbability:
            bestPattern = stringKmer
            bestProbability = newProbablity
    return bestPattern

def calculateProbability(profile, string):
    probablity = 1
    for i in range(0, len(string)):
        if string[i] == 'A':
            probablity = probablity * profile[i][0]
        elif string[i] == 'C':
            probablity = probablity * profile[i][1]
        elif string[i] == 'G':
            probablity = probablity * profile[i][2]
        elif string[i] == 'T':
            probablity = probablity * profile[i][3]
    return probablity

def hemmingDistance(str1, str2):
    distance = 0
    for i in range(len(str1)):
        if str1[i] != str2[i]:
            distance += 1
    return distance


t = min(int(intValues[1]), len(data) - 1)
finalSequences = randomizedMotifSearch(data[1:], int(intValues[0]), int(intValues[1]))
for finalSequence in finalSequences:
    print(finalSequence)