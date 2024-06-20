#!/usr/bin/env python3
from collections import Counter
import argparse
import re
import random
parser = argparse.ArgumentParser()

required = parser.add_argument_group('Required arguments')
required.add_argument('-i', '--input', type=str, help='Path to input file', required=True)
args = parser.parse_args()
fileIn = open(args.input)
spectrum = list(map(int, fileIn.readline().strip().split()))

def expand(candidatePeptides, aminoAcidandMasses):
    newPeptides = []
    for peptide in candidatePeptides:
        for aminoAcid in aminoAcidandMasses:
            newPeptides.append(peptide + [aminoAcid])
    return newPeptides


def Mass(peptide, aminoAcidandMasses):
    return sum(aminoAcidandMasses[aminoAcid] for aminoAcid in peptide)

def cycloSpectrum(candidatePeptide, aminoAcidandMasses):
    peptideLength = len(candidatePeptide)
    spectrum = [0]
    longerPeptide = candidatePeptide + candidatePeptide
    for k in range(1, peptideLength):
        for i in range(peptideLength):
            subpeptide = longerPeptide[i:i + k]
            spectrum.append(Mass(subpeptide, aminoAcidandMasses))
    spectrum.append(Mass(candidatePeptide, aminoAcidandMasses))
    return sorted(spectrum)

def consistencyCheck(candidatePeptide, spectrum, aminoAcidandMasses):
    peptideLength = len(candidatePeptide)
    newSpectrum = [0]
    for i in range(peptideLength):
        for j in range(i + 1, peptideLength + 1):
            subpeptide = candidatePeptide[i:j]
            mass = Mass(subpeptide, aminoAcidandMasses)
            newSpectrum.append(mass)
    for mass in newSpectrum:
        if mass not in spectrum:
            return False        
    if Mass(candidatePeptide, aminoAcidandMasses) > parentMass(spectrum):
        return False
    else:
        return True
    
def parentMass(spectrum):
    return max(spectrum)

def CyclopeptideSequencing(spectrum):
    finalPeptides = set()
    aminoAcidandMasses = {
        'G': 57, 'A': 71, 'S': 87, 'P': 97, 'V': 99, 'T': 101, 'C': 103, 'I': 113, 'L': 114,
        'N': 115, 'D': 128, 'K': 129, 'E': 131, 'M': 137, 'H': 147, 'F': 156, 'R': 163, 'Y': 186
    }
    candidatePeptides = [[]]
    matchingPeptides = {''}
    while candidatePeptides:
        newPeptides = []
        for candidatePeptide in candidatePeptides:
            massPeptide = Mass(candidatePeptide, aminoAcidandMasses)
            if massPeptide == parentMass(spectrum):
                if cycloSpectrum(candidatePeptide, aminoAcidandMasses) == spectrum:
                    matchingPeptides.add("-".join(map(str, [aminoAcidandMasses[aminoAcid] for aminoAcid in candidatePeptide])))
            elif not not consistencyCheck(candidatePeptide, spectrum, aminoAcidandMasses):
                newPeptides.extend(expand([candidatePeptide], aminoAcidandMasses))
        candidatePeptides = newPeptides
    finalPeptides.update(matchingPeptides)
    return finalPeptides

displayedPeptides = []    
matchedPeptides = CyclopeptideSequencing(spectrum)
for peptide in matchedPeptides:
    formattedPeptide  = '-'.join(peptide.split('-'))
    displayedPeptides.append(formattedPeptide)

print(' '.join(displayedPeptides))