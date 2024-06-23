#!/usr/bin/env python3
from collections import Counter
import argparse
parser = argparse.ArgumentParser()
required = parser.add_argument_group('Required arguments')
required.add_argument('-i', '--input', type=str, help='Path to input file', required=True)
args = parser.parse_args()
fileIn = open(args.input)
sampleData = fileIn.read().splitlines()

def genomeToEdges(genome):
    listOfEdges = []
    for chromosome in genome:
        nodes = genomeToCycleRepresentation(chromosome)
        listOfEdges.extend([(nodes[j], nodes[j + 1]) for j in range(len(nodes) - 1)])
        listOfEdges.append([nodes[-1], nodes[0]])
    return listOfEdges

def findConnectedEdge(current, edges):
    for edge in edges:
        if current[0] in edge or current[1] in edge:
            return edge
    return None

def genomeToCycleRepresentation(chromosome):
    listOfNodes = []
    for genomicSegment in chromosome:
        if genomicSegment < 0:
            listOfNodes.extend([-2 * genomicSegment, -2 * genomicSegment - 1])      
        else:
            listOfNodes.extend([2 * genomicSegment - 1, 2 * genomicSegment])
           
    return listOfNodes

def twoBreakDistance(genome1, genome2):
    edges = genomeToEdges(genome1 + genome2)
    blocks = {node for edge in edges for node in edge}
    cycles = []

    while edges:
        currentEdge = edges.pop(0)
        cycle = [currentEdge]

        while currentEdge := findConnectedEdge(currentEdge, edges):
            cycle.append(currentEdge)
            edges.remove(currentEdge)

        cycles.append(cycle)
    distance =  ((len(blocks) - 2 * len(cycles)) / 2)-1   

    return distance

p = [list(map(int, chrom.strip().split())) for chrom in sampleData[0][1:-1].split(")(")]
q = [list(map(int, chrom.strip().split())) for chrom in sampleData[1][1:-1].split(")(")]


result = twoBreakDistance(p, q)
print(int(result))