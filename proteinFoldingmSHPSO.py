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

def sample(center,radius):
	transformedCenter = transform(center)
	#print("Original :",end=' ')
	#print(''.join(transformedCenter))
	indexes = np.arange(len(center))
	hammingDistance = round(radius)
	if hammingDistance == 0:
		return np.array([])
	idxlist = np.random.choice(indexes, hammingDistance, replace = False)
	for i in idxlist:
		if transformedCenter[i] == 'L':
			toss = random.random()
			if toss < 0.5:
				transformedCenter[i] = 'R'
			else:
				transformedCenter[i] = 'C'
		elif transformedCenter[i] == 'C':
			toss = random.random()
			if toss < 0.5:
				transformedCenter[i] = 'R'
			else:
				transformedCenter[i] = 'L'
		else:
			toss = random.random()
			if toss < 0.5:
				transformedCenter[i] = 'C'
			else:
				transformedCenter[i] = 'L'
	#print(''.join(transformedCenter))
	returnCenter = center
	for i in idxlist:
		if(transformedCenter[i] == 'L'):
			returnCenter[i] = np.random.uniform(low = -3.0 , high = -1.0)
		elif(transformedCenter[i] == 'C'):
			returnCenter[i] = np.random.uniform(low = -1.0 , high = 1.0)
		else:
			returnCenter[i] = np.random.uniform(low = 1.0 , high = 3.0)

	return returnCenter

class particle:
	def __init__(self, positionVector, velocityVector, fitness,radius):
		self.positionVector = positionVector
		self.velocityVector = velocityVector
		self.fitness = fitness
		self.pbestPositionVector = positionVector
		self.pbestCost = fitness
		self.radius = radius

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

np.random.seed(12345)

#HPHPPHHPHPPHPHHPPHPH

aminoAcid = "HPHPPHHPHPPHPHHPPHPH"

dimension = len(aminoAcid) - 2
posMin = -3
posMax = 3

vMax = 6
vMin = -6

#Parameters of PSO

maxIterations = 2000
populationSize = len(aminoAcid)

inertiaWeight = 0.7

hypersphereRadius = len(aminoAcid)/2

dampingFactor = 1

c1 = 1.49618
c2 = 1.49618
c3 = 1.01

#Step 1

population = []

gbestPositionVector = np.zeros(dimension)
gbestCost = sys.maxsize
gbestParticle = -1

for i in range(populationSize):
	vector = np.random.uniform(low = posMin, high = posMax, size = dimension)
	fitness = structureEnergy(transform(vector),aminoAcid)
	population.append(particle(vector,np.zeros(dimension),fitness,hypersphereRadius))

	if fitness < gbestCost:
		gbestPositionVector = vector
		gbestCost = fitness
		gbestParticle = i

gbestCostAtEachIteration = []

print("%d %d" % (0,gbestCost))

gbestCostAtEachIteration.append(gbestCost)

#Step 2

for it in range(maxIterations):

	NGbestVector = population[gbestParticle].positionVector
	
	for i in range(populationSize):
		hbestVector = population[i].positionVector
		hbestCost = population[i].fitness

		pointsCount = 0

		points = []

		flag = True

		while pointsCount != 10:
			point = sample(population[i].positionVector , population[i].radius)
			if point != np.array([]):
				points.append(point)
				pointsCount += 1
			else:
				flag = False
				break

		population[i].radius = hypersphereRadius*(1 - (it/maxIterations))
		
		if flag == True:
			for j in range(len(points)):
				pointFitness = structureEnergy(transform(points[j]),aminoAcid)
				print(pointFitness)

				if pointFitness < hbestCost:
					hbestVector = points[j]
					hbestCost = pointFitness

			if hbestCost < population[i].pbestCost:
				population[i].pbestCost = hbestCost
				population[i].pbestPositionVector = hbestVector
		
			if i == gbestParticle:
				NGbestVector = hbestVector



#Step 3

	for i in range(populationSize):
		
		personalInfluence = population[i].pbestPositionVector - population[i].positionVector
		personalInfluence = c1*random.random()*personalInfluence

		socialInfluence = gbestPositionVector - population[i].positionVector
		socialInfluence = c2*random.random()*socialInfluence

		NGbestInfluence = NGbestVector - population[i].positionVector
		NGbestInfluence = c3*random.random()*NGbestInfluence

		population[i].velocityVector = inertiaWeight * population[i].velocityVector + personalInfluence + socialInfluence + NGbestInfluence

		population[i].velocityVector = np.clip(population[i].velocityVector , vMin, vMax )

		population[i].positionVector = population[i].positionVector + population[i].velocityVector

		population[i].positionVector = np.clip(population[i].positionVector , posMin, posMax)

		for dim in range(dimension):
			if population[i].positionVector[dim] == posMax or population[i].positionVector[dim] == posMin:
				population[i].positionVector[dim] = np.random.uniform(posMin,posMax)

		#print(''.join(transform(population[i].positionVector)))
		#print(population[i].positionVector)

		fitness = structureEnergy(transform(population[i].positionVector),aminoAcid)

		if fitness < population[i].pbestCost:
			population[i].pbestPositionVector = population[i].positionVector
			population[i].pbestCost = fitness

			if fitness < gbestCost:
				gbestPositionVector = population[i].positionVector
				gbestCost = fitness
				gbestParticle = i

	gbestCostAtEachIteration.append(gbestCost)

	inertiaWeight *= dampingFactor

	print("%d %d" % (it + 1,gbestCost))	