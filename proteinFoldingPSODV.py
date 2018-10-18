import random,sys,math
from operator import add,sub
import matplotlib.pyplot as plt
import numpy as np
from ProteinEnergy import structureEnergy

class particle:
	def __init__(self, positionVector, velocityVector, fitness):
		self.positionVector = positionVector
		self.velocityVector = velocityVector
		self.fitness = fitness


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

np.random.seed(344470)

aminoAcid = "HHHHHHHHHHHHPHPHPPHHPPHHPPHPPHHPPHHPPHPPHHPPHHPPHPHPHHHHHHHHHHHH"
dimension = len(aminoAcid) - 2

posMin = -3
posMax = 3

vMax = 2
vMin = -2

#Parameters of PSO

maxIterations = 12000

populationSize = 50

inertiaWeight = 0.7

dampingFactor = 1.0

c2 = 2

CR = 0.7

beta = 0.8

#Initialisation

population = []

gbestPositionVector = np.zeros(dimension)
gbestCost = sys.maxsize

for i in range(populationSize):
	vector = np.random.uniform(low = posMin, high = posMax, size = dimension)
	fitness = structureEnergy(transform(vector),aminoAcid)
	population.append(particle(vector,np.zeros(dimension),fitness))

	if fitness < gbestCost:
		gbestPositionVector = vector
		gbestCost = fitness

gbestCostAtEachIteration = []

gbestCostAtEachIteration.append(gbestCost)

print("%d %d" % (0,gbestCost))

#Main Loop starts

counter = 0

for it in range(maxIterations):
	counter += 1
	if counter == 300:
		counter = 0
		for pi in range(populationSize):
			population[pi].positionVector = np.random.uniform(low = posMin, high = posMax, size = dimension)
			population[pi].fitness = structureEnergy(transform(population[pi].positionVector),aminoAcid)
			if population[pi].fitness < gbestCost:
				gbestCost = population[pi].fitness
				gbestPositionVector = population[pi].positionVector

	for i in range(populationSize):
		
		indexes = [idx for idx in range(populationSize) if idx != i]
		
		r1, r2 = np.random.choice(indexes,2,replace = False)
		x1 = population[r1].positionVector
		x2 = population[r2].positionVector
		
		differenceVector = beta*(x1- x2)
		differenceVector = np.clip(differenceVector,posMin,posMax)

		socialInfluence = gbestPositionVector - population[i].positionVector
		socialInfluence = c2*random.random()*socialInfluence

		for j in range(dimension):
			if random.random() < CR:
				population[i].velocityVector[j] = inertiaWeight*population[i].velocityVector[j] + differenceVector[j] + socialInfluence[j]

		population[i].velocityVector = np.clip(population[i].velocityVector , vMin, vMax )

		trialVector = population[i].positionVector + population[i].velocityVector

		np.clip(trialVector , posMin, posMax, out = trialVector)


		fitness = structureEnergy(transform(trialVector),aminoAcid)
		if fitness  < population[i].fitness:
			population[i].positionVector = trialVector


		#print(population[i].positionVector)
		#print(''.join(transform(population[i].positionVector)))

		if fitness < gbestCost:
			gbestPositionVector = population[i].positionVector
			gbestCost = fitness
			counter = 0

	gbestCostAtEachIteration.append(gbestCost)

	inertiaWeight *= dampingFactor

	print("%d %d" % (it + 1,gbestCost))	