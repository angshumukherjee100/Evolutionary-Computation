import random,sys,math
from operator import add,sub
import matplotlib.pyplot as plt 


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

posMin = -600
posMax = 600

vMax = 50
vMin = -50

#Parameters of PSO

maxIterations = 100

populationSize = 20

inertiaWeight = 1

dampingFactor = 0.99

c1 = 2
c2 = 2

#Initialisation

population = []

gbestPositionVector = [0]*dimension
gbestCost = sys.maxsize

for i in range(populationSize):
	vector = [random.uniform(posMin,posMax) for j in range(dimension)]
	fitness = sphere(vector)
	population.append(particle(vector,[0]*dimension,fitness))

	if fitness < gbestCost:
		gbestPositionVector = list(vector)
		gbestCost = fitness

gbestCostAtEachIteration = []


#Main Loop starts

for it in range(maxIterations):
	for i in range(populationSize):

		func1 = lambda x : c1*random.random()*x
		func2 = lambda x : c2*random.random()*x
		
		personalInfluence = list(map(sub, population[i].pbestPositionVector, population[i].positionVector))
		personalInfluence = list(map(func1, personalInfluence))

		socialInfluence = list(map(sub, gbestPositionVector, population[i].positionVector))
		socialInfluence = list(map(func2, socialInfluence))

		population[i].velocityVector = list(map(lambda x: inertiaWeight * x, population[i].velocityVector))

		population[i].velocityVector = list(map(add,population[i].velocityVector,personalInfluence))
		population[i].velocityVector = list(map(add,population[i].velocityVector,socialInfluence))


		for velocityComponent in range(len(population[i].velocityVector)):
			if population[i].velocityVector[velocityComponent] > vMax:
				population[i].velocityVector[velocityComponent] = vMax
			if population[i].velocityVector[velocityComponent] < vMin:
				population[i].velocityVector[velocityComponent] = vMin


		population[i].positionVector = list(map(add,population[i].positionVector,population[i].velocityVector))

		for positionComponent in range(len(population[i].positionVector)):
			if population[i].positionVector[positionComponent] > posMax:
				population[i].positionVector[positionComponent] = posMax
			if population[i].positionVector[positionComponent] < posMin:
				population[i].positionVector[positionComponent] = posMin

		fitness = sphere(population[i].positionVector)

		if fitness < population[i].pbestCost:
			population[i].pbestPositionVector = list(population[i].positionVector)
			population[i].pbestCost = fitness

			if fitness < gbestCost:
				gbestPositionVector = list(population[i].positionVector)
				gbestCost = fitness

	gbestCostAtEachIteration.append(gbestCost)

	inertiaWeight *= dampingFactor

	print("%d %e" % (it + 1,gbestCost))			

#Plots
plt.xlabel("Iteration Number")
plt.ylabel("Fitness")

#plt.plot([i for i in range(1,maxIterations + 1)],gbestCostAtEachIteration)
plt.semilogy([i for i in range(1,maxIterations + 1)],gbestCostAtEachIteration)
plt.grid()
plt.show()