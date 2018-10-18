import random,sys,math
import matplotlib.pyplot as plt
import numpy as np
from ProteinEnergy import structureEnergy

def transform(x):
	temp = []
	for i in x:
		if i >= -3.0 and i < -1.0:
			temp.append('L')
		elif i >= -1.0 and i < 1.0:
			temp.append('C')
		elif i >= 1.0 and i <= 3.0:
			temp.append('R')
	return temp

np.random.seed(3611)

aminoAcid = "HHPPHPPHPPHPPHPPHPPHPPHH"

dimension = len(aminoAcid) - 2
population = []

populationSize = len(aminoAcid)
maxIterations = 2000

F = 0.8
Cr = 0.7

lowerBound = -3
upperBound = 3

bestVector = np.zeros(dimension)
bestFitness = 2		#Maximum value can be 1

for i in range(populationSize):
	vector = np.random.uniform(low = lowerBound, high = upperBound, size = dimension)
	fitness = structureEnergy(transform(vector),aminoAcid)
	population.append(vector)
	if fitness < bestFitness:
		bestFitness = fitness
		bestVector = vector

print("%d %d" %(0,bestFitness))


for i in range(maxIterations):
	for j in range(populationSize):
		indexes = [idx for idx in range(populationSize) if idx != j]
		
		r1, r2, r3 = np.random.choice(indexes,3,replace = False)
		x1 = population[r1]
		x2 = population[r2]
		x3 = population[r3]
		
		donorVector = x1 + F*(x2 - x3)
		donorVector = np.clip(donorVector,lowerBound,upperBound)
		
		#Binomial CrossOver
		crossPoints = np.random.rand(dimension) < Cr
		trialVector = np.where(crossPoints, donorVector, population[j])
		
		trialVectorFitness = structureEnergy(transform(trialVector),aminoAcid)

		if trialVectorFitness < structureEnergy(transform(population[j]),aminoAcid):
			population[j] = trialVector

			if trialVectorFitness < bestFitness:
				bestVector = trialVector
				bestFitness = trialVectorFitness
		#print(''.join(transform(population[j])))
	print("%d %d" %(i + 1,bestFitness))