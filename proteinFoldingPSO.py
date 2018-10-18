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
		self.pbestPositionVector = positionVector
		self.pbestCost = fitness


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

np.random.seed(211209)

aminoAcid = "PPHHHPHHHHHHHHPPPHHHHHHHHHHPHPPPHHHHHHHHHHHHPPPPHHHHHHPHHPHP"
dimension = len(aminoAcid) - 2

posMin = -3
posMax = 3

vMax = 2
vMin = -2

#Parameters of PSO

maxIterations = 4000

populationSize = len(aminoAcid)

inertiaWeight = 0.7

dampingFactor = 1.0

c1 = 2
c2 = 2

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

for it in range(maxIterations):
	for i in range(populationSize):
		
		personalInfluence = population[i].pbestPositionVector - population[i].positionVector
		personalInfluence = c1*random.random()*personalInfluence

		socialInfluence = gbestPositionVector - population[i].positionVector
		socialInfluence = c2*random.random()*socialInfluence

		population[i].velocityVector = inertiaWeight * population[i].velocityVector + personalInfluence + socialInfluence

		population[i].velocityVector = np.clip(population[i].velocityVector , vMin, vMax )

		population[i].positionVector = population[i].positionVector + population[i].velocityVector

		np.clip(population[i].positionVector , posMin, posMax, out = population[i].positionVector)

		fitness = structureEnergy(transform(population[i].positionVector),aminoAcid)

		#print(population[i].positionVector)
		#print(''.join(transform(population[i].positionVector)))

		if fitness < population[i].pbestCost:
			population[i].pbestPositionVector = population[i].positionVector
			population[i].pbestCost = fitness

			if fitness < gbestCost:
				gbestPositionVector = population[i].positionVector
				gbestCost = fitness

	gbestCostAtEachIteration.append(gbestCost)

	inertiaWeight *= dampingFactor

	print("%d %d" % (it + 1,gbestCost))	


#Plots
'''
plt.xlabel("Iteration Number")
plt.ylabel("Fitness")
plt.plot([i for i in range(maxIterations + 1)],gbestCostAtEachIteration)
plt.grid()
plt.show()
'''