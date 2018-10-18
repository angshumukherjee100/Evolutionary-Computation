import random
import numpy as np
import math
from ProteinEnergy import structureEnergy


def cauchyRandomWalk(x):
    temp = x + np.random.standard_cauchy(len(x))
    return temp

def transform(x):
    temp = []
    for i in x:
        if i < -3.0:
            temp.append('L')
        elif i >= -3.0 and i < -1.0:
            temp.append('L')
        elif i >= -1.0 and i < 1.0:
            temp.append('C')
        elif i >= 1.0 and i <= 3.0:
            temp.append('R')
        else:
            temp.append('R')
    return temp

aminoAcid = "HPHPPHHPHPPHPHHPPHPH"

lb = -3.0
ub = 3.0

F = 0.9
Cr = 0.7

dim = len(aminoAcid) - 2

SearchAgents_no = 20

Max_iter = 2000

np.random.seed()
    
Alpha_pos=np.zeros(dim)
Alpha_score=float("inf")
    
Beta_pos=np.zeros(dim)
Beta_score=float("inf")
    
Delta_pos=np.zeros(dim)
Delta_score=float("inf")
    
#Initialize the positions of search agents
population = []

for i in range(SearchAgents_no):
    vector = np.random.uniform(low = lb, high = ub, size = dim)
    fitness = structureEnergy(transform(vector),aminoAcid)
    
    if fitness < Alpha_score:
        Delta_score = Beta_score
        Delta_pos = Beta_pos
        Beta_score = Alpha_score
        Beta_pos = Alpha_pos
        Alpha_score = fitness 
        Alpha_pos = vector
            
    if (fitness > Alpha_score and fitness < Beta_score ):
        Delta_score = Beta_score
        Delta_pos = Beta_pos
        Beta_score = fitness
        Beta_pos = vector
            
    if (fitness > Alpha_score and fitness > Beta_score and fitness < Delta_score): 
        Delta_score = fitness # Update delta
        Delta_pos = vector

    population.append(vector)

print(str(0)+' '+ str(Alpha_score))
  
Convergence_curve = np.zeros(Max_iter + 1)

Convergence_curve[0] = Alpha_score
  
# Main loop
for l in range(0,Max_iter):
                   
    a = 2 - l*((2)/Max_iter)

    tempPos = cauchyRandomWalk(Alpha_pos)
    tempPos = np.clip(tempPos, lb, ub)
    tempfitness = structureEnergy(transform(tempPos),aminoAcid)        
    if tempfitness < Alpha_score:
        Delta_score = Beta_score
        Delta_pos = Beta_pos
        Beta_score = Alpha_score
        Beta_pos = Alpha_pos
        Alpha_pos = tempPos
        Alpha_score = tempfitness

    tempPos = cauchyRandomWalk(Beta_pos)
    tempPos = np.clip(tempPos, lb, ub)
    tempfitness = structureEnergy(transform(tempPos),aminoAcid)        
    if tempfitness <  Beta_score and tempfitness > Alpha_score:
        Delta_score = Beta_score
        Delta_pos = Beta_pos
        Beta_pos = tempPos
        Beta_score = tempfitness
    elif tempfitness < Alpha_score:
            Beta_pos = Alpha_pos
            Beta_score = Alpha_score
            Alpha_score = tempfitness
            Alpha_pos = tempPos
        
    tempPos = cauchyRandomWalk(Delta_pos)
    tempPos = np.clip(tempPos, lb, ub)
    tempfitness = structureEnergy(transform(tempPos),aminoAcid)        
    if tempfitness < Delta_score and tempfitness > Beta_score:
        Delta_pos = tempPos
        Delta_score = tempfitness
    elif tempfitness < Beta_score and tempfitness > Alpha_score:
        Delta_score = Beta_score
        Delta_pos = Beta_pos
        Beta_pos = tempPos
        Beta_score = tempfitness
    elif tempfitness < Alpha_score:
        Delta_score = Beta_score
        Delta_pos = Beta_pos
        Beta_pos = Alpha_pos
        Beta_score = Alpha_score


    for i in range(0,SearchAgents_no):
        if random.random() < 0.9:
            indexes = [idx for idx in range(SearchAgents_no) if idx != i]
            r1, r2, r3, r4, r5 = np.random.choice(indexes,5,replace = False)
            x1 = population[r1]
            x2 = population[r2]
            x3 = population[r3]
            x4 = population[r4]
            x5 = population[r5]
		
            donorVector = x1 + F*(x2 - x3) + F*(x4 - x5)
            donorVector = np.clip(donorVector,lb,ub)
		
            #Binomial CrossOver
            crossPoints = np.random.rand(dim) < Cr
            trialVector = np.where(crossPoints, donorVector, population[i])
            trialVectorFitness = structureEnergy(transform(trialVector),aminoAcid)
            if trialVectorFitness < structureEnergy(transform(population[i]),aminoAcid):
                population[i] = trialVector
        else:

            for j in range (0,dim):     
                           
                r1 = random.random()
                r2 = random.random()
                
                A1 = 2*a*r1 - a
                C1 = 2*r2
               
                D_alpha = abs(C1*Alpha_pos[j] - population[i][j])
                X1 = Alpha_pos[j] - A1*D_alpha

                r1 = random.random()
                r2 = random.random()
                
                A2 = 2*a*r1 - a
                C2 = 2*r2
                
                D_beta = abs(C2*Beta_pos[j] - population[i][j])
                X2 = Beta_pos[j] - A2*D_beta       
                
                r1 = random.random()
                r2 = random.random() 
                
                A3 = 2*a*r1 - a
                C3 = 2*r2
                
                D_delta = abs(C3*Delta_pos[j] - population[i][j])
                X3 = Delta_pos[j] - A3*D_delta;              
                
                population[i][j] = (X1+X2+X3)/3

        population[i] = np.clip(population[i], lb, ub)

        # Calculate objective function for each search agent
        fitness = structureEnergy(transform(population[i]),aminoAcid)
            
        # Update Alpha, Beta, and Delta
        if fitness < Alpha_score :
            Delta_score = Beta_score
            Delta_pos = Beta_pos
            Beta_score = Alpha_score
            Beta_pos = Alpha_pos
            Alpha_score = fitness # Update alpha
            Alpha_pos = population[i]
            
        if (fitness > Alpha_score and fitness < Beta_score):
            Delta_score = Beta_score
            Delta_pos = Beta_pos
            Beta_score = fitness  # Update beta
            Beta_pos = population[i]
            
        if (fitness > Alpha_score and fitness > Beta_score and fitness < Delta_score): 
            Delta_score = fitness # Update delta
            Delta_pos = population[i]
                 
        
    Convergence_curve[l + 1] = Alpha_score

    print(str(l + 1)+' '+ str(Alpha_score))