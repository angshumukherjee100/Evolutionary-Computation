import random,sys,math
from operator import add,sub
import matplotlib.pyplot as plt
import numpy as np

class particle:
	def __init__(self, positionVector, velocityVector, fitness):
		self.positionVector = positionVector
		self.velocityVector = velocityVector
		self.fitness = fitness
		self.pbestPositionVector = positionVector
		self.pbestCost = fitness


def sphere(x):
	z = 0.0
	for xi in x:
		z += xi**2
	return z


random.seed()

dimension = 100

posMin = -100
posMax = 100

vMax = 20
vMin = -20

#Parameters of PSO

maxIterations = 10000

populationSize = 20

inertiaWeight = 1

dampingFactor = 1

c1 = 2
c2 = 2

#Initialisation

population = []

gbestPositionVector = np.zeros(dimension)
gbestCost = sys.maxsize

for i in range(populationSize):
	vector = np.random.uniform(low = posMin, high = posMax, size = dimension)
	fitness = sphere(vector)
	population.append(particle(vector,np.zeros(dimension),fitness))

	if fitness < gbestCost:
		gbestPositionVector = vector
		gbestCost = fitness

gbestCostAtEachIteration = []

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

		np.clip(population[i].positionVector , posMin, posMax, out = population[i].positionVector )

		fitness = sphere(population[i].positionVector)

		if fitness < population[i].pbestCost:
			population[i].pbestPositionVector = population[i].positionVector
			population[i].pbestCost = fitness

			if fitness < gbestCost:
				gbestPositionVector = population[i].positionVector
				gbestCost = fitness

	gbestCostAtEachIteration.append(gbestCost)

	inertiaWeight *= dampingFactor

	print("%d %e" % (it + 1,gbestCost))			

#Plots
plt.xlabel("Iteration Number")
plt.ylabel("Fitness")

plt.semilogy([i for i in range(1,maxIterations + 1)],gbestCostAtEachIteration)
plt.grid()
plt.show()