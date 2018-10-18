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

dim = len(aminoAcid) - 2

SearchAgents_no = 20

Max_iter = 4000

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
        for j in range (0,dim):     
                           
            r1 = random.random() # r1 is a random number in [0,1]
            r2 = random.random() # r2 is a random number in [0,1]
                
            A1 = 2*a*r1 - a; # Equation (3.3)
            C1 = 2*r2; # Equation (3.4)
               
            D_alpha = abs(C1*Alpha_pos[j] - population[i][j]); # Equation (3.5)-part 1
            X1 = Alpha_pos[j] - A1*D_alpha; # Equation (3.6)-part 1
                           
            r1 = random.random()
            r2 = random.random()
                
            A2 = 2*a*r1 - a; # Equation (3.3)
            C2 = 2*r2; # Equation (3.4)
                
            D_beta = abs(C2*Beta_pos[j] - population[i][j]); # Equation (3.5)-part 2
            X2 = Beta_pos[j] - A2*D_beta; # Equation (3.6)-part 2       
                
            r1 = random.random()
            r2 = random.random() 
                
            A3 = 2*a*r1 - a; # Equation (3.3)
            C3 = 2*r2; # Equation (3.4)
                
            D_delta = abs(C3*Delta_pos[j] - population[i][j]); # Equation (3.5)-part 3
            X3 = Delta_pos[j] - A3*D_delta; # Equation (3.5)-part 3             
                
            population[i][j] = (X1+X2+X3)/3  # Equation (3.7)

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