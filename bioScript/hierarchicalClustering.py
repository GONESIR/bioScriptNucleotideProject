#!/usr/bin/env python3
from collections import Counter
import argparse
parser = argparse.ArgumentParser()
required = parser.add_argument_group('Required arguments')
required.add_argument('-i', '--input', type=str, help='Path to input file', required=True)
args = parser.parse_args()
fileIn = open(args.input)
n = int(fileIn.readline().strip())
distanceMatrix = [list(map(float, fileIn.readline().split())) for item in range(n)]

def findClosestClusters(clusters, distanceMatrix):
    minDistance = float('inf')
    closestI, closestJ = -1, -1

    for i in range(len(clusters) - 1):
        for j in range(i + 1, len(clusters)):
            distance = computeAverageDistance(clusters[i], clusters[j], distanceMatrix)
            if distance < minDistance:
                minDistance = distance
                closestI, closestJ = i, j

    return closestI, closestJ

def mergeClusters(cluster1, cluster2):
    return cluster1 + cluster2

def computeAverageDistance(cluster1, cluster2, distance_matrix):
    return sum(distance_matrix[i][j] for i in cluster1 for j in cluster2) / (len(cluster1) * len(cluster2))


def hierarchicalClustering(distanceMatrix):
    newListOfClusters = []
    n = len(distanceMatrix)
    clusters = [[i] for i in range(n)]
    
    while len(clusters) > 1:
        i, j = findClosestClusters(clusters, distanceMatrix)
        newCluster = mergeClusters(clusters[i], clusters[j])
        clusters = [item for item in clusters if item not in [clusters[i], clusters[j]]]
        clusters.append(newCluster)
        newListOfClusters.append(newCluster)

    return newListOfClusters

resultantClusters = hierarchicalClustering(distanceMatrix)

for cluster in resultantClusters:
    print(' '.join(map(lambda x: str(x + 1), cluster)))