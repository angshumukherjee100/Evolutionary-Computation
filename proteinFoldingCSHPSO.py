import random,sys,math
import matplotlib.pyplot as plt 
import numpy as np
from scipy.special import gammainc
from ProteinEnergy import structureEnergy
from operator import itemgetter


class particle:
	def __init__(self, positionVector, velocityVector, fitness,radius):
		self.positionVector = positionVector
		self.velocityVector = velocityVector
		self.fitness = fitness
		self.pbestPositionVector = positionVector
		self.pbestFitness = fitness
		self.radius = radius

def sample(center,radius):
	r = radius
	ndim = center.size
	x = np.random.normal(size = (1,ndim))
	ssq = np.sum(x**2, axis = 1)
	fr = r*gammainc(ndim/2,ssq/2)**(1/ndim)/np.sqrt(ssq)
	frtiled = np.tile(fr.reshape(1,1),(1,ndim))
	p = center + np.multiply(x,frtiled)
	return p

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

def euclideanDistance(x,y):
	d = 0.0
	dim = len(x)
	for i in range(dim):
		d += (x[i]-y[i])**2
	return math.sqrt(d)


np.random.seed(23656)

aminoAcid = "HPHPPHHPHPPHPHHPPHPH"

dimension = len(aminoAcid) - 2
population = []

populationSize = 20
maxIterations = 2000

F = 0.8

posMin = -3
posMax = 3

vMax = 6
vMin = -6

inertiaWeight = 0.7

hypersphereRadius = 7

c1 = 1.49618
c2 = 1.49618
c3 = 1.01

gbestPositionVector = np.zeros(dimension)
gbestFitness = sys.maxsize
gbestParticle = -1

for i in range(populationSize):
	vector = np.random.uniform(low = posMin, high = posMax, size = dimension)
	fitness = structureEnergy(transform(vector),aminoAcid)
	population.append(particle(vector,np.zeros(dimension),fitness,hypersphereRadius))

	if fitness < gbestFitness:
		gbestPositionVector = vector
		gbestFitness = fitness
		gbestParticle = i

gbestFitnessAtEachIteration = []

print("%d %d" % (0,gbestFitness))

gbestFitnessAtEachIteration.append(gbestFitness)

for it in range(maxIterations):

	distances = []

	for i in range(populationSize):
		for j in range(i+1,populationSize):
			distances.append((i,j,euclideanDistance(population[i].positionVector,population[j].positionVector)))

	distances = sorted(distances, key = itemgetter(2))

	dict = {}
	for i in range(populationSize):
		dict[i] = -1
	#0 is for DE and 1 is for SHPSO in the dictionary
	for i in range(len(distances)):
		if dict[distances[i][0]] == -1:
			dict[distances[i][0]] = 0
		if dict[distances[i][1]] == -1:
			dict[distances[i][1]] = 1
	#Swarm Splitted

	#DE starts
	indexes = [idx for idx in range(len(dict)) if dict[idx] == 0]

	for j in indexes:
		tempIndexList = indexes.copy()
		tempIndexList.remove(j)
		r1, r2, r3 = np.random.choice(tempIndexList, 3, replace = False)
		x1 = population[r1].positionVector
		x2 = population[r2].positionVector
		x3 = population[r3].positionVector
		
		donorVector = x1 + F*(x2 - x3)
		donorVector = np.clip(donorVector,posMin,posMax)
		
		#Binomial CrossOver
		Cr = (0.9*(population[j].fitness/gbestFitness)) +  0.1
		crossPoints = np.random.rand(dimension) < Cr
		trialVector = np.where(crossPoints, donorVector, population[j].positionVector)
		
		trialVectorFitness = structureEnergy(transform(trialVector),aminoAcid)

		if trialVectorFitness < structureEnergy(transform(population[j].positionVector),aminoAcid):
			population[j].pbestPositionVector = trialVector
			population[j].pbestFitness = trialVectorFitness
	
	#Swarm merged

	NGbestVector = population[gbestParticle].positionVector
	
	for i in range(populationSize):
		hbestVector = population[i].positionVector
		hbestCost = population[i].fitness

		pointsCount = 0

		points = []

		while pointsCount != 10:
			point = sample(population[i].pbestPositionVector , population[i].radius)

			flag = True

			for vectorComponent in point[0]:
				if vectorComponent < posMin or vectorComponent > posMax:
					flag = False
					break

			if flag == True:
				points.append(point[0])
				pointsCount += 1

		population[i].radius = hypersphereRadius*(1 - (it/maxIterations))
		
		for j in range(len(points)):
			pointFitness = structureEnergy(transform(points[j]),aminoAcid)

			if pointFitness < hbestCost:
				hbestVector = points[j]
				hbestCost = pointFitness

		if hbestCost < population[i].pbestFitness:
			population[i].pbestFitness = hbestCost
			population[i].pbestPositionVector = hbestVector
		
		if i == gbestParticle:
			NGbestVector = hbestVector

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

		fitness = structureEnergy(transform(population[i].positionVector),aminoAcid)

		if fitness < population[i].pbestFitness:
			population[i].pbestPositionVector = population[i].positionVector
			population[i].pbestFitness = fitness

			if fitness < gbestFitness:
				gbestPositionVector = population[i].positionVector
				gbestFitness = fitness
				gbestParticle = i

	gbestFitnessAtEachIteration.append(gbestFitness)

	print("%d %d" % (it + 1,gbestFitness))
