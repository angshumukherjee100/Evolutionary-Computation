def structureEnergy(orientation , proteinMolecule):
	#initialise grid

	numberOfMolecules = len(proteinMolecule)

	grid = []
	positionsOfAminoAcids = {}
	aminoAcidPositions = {}

	for i in range(2*numberOfMolecules):
		tempRow = []
		for j in range(2*numberOfMolecules):
			tempRow.append('.')
		grid.append(tempRow)

	#place amino acids on the grid

	i = numberOfMolecules
	j = numberOfMolecules

	grid[i][j] = proteinMolecule[0]
	positionsOfAminoAcids[(i,j)] = (0, proteinMolecule[0])
	aminoAcidPositions[0] = (i,j)

	j += 1
	grid[i][j] = proteinMolecule[1]
	positionsOfAminoAcids[(i,j)] = (1, proteinMolecule[1])
	aminoAcidPositions[1] = (i,j)

	prevDir = 'R'

	for pos in range(len(orientation)):
		
		if orientation[pos] == 'C':
			if prevDir == 'R':
				j += 1
			elif prevDir == 'L':
				j -= 1
			elif prevDir == 'U':
				i -= 1
			elif prevDir == 'D':
				i += 1
		elif orientation[pos] == 'R':
			if prevDir == 'R':
				i += 1
				prevDir = 'D'
			elif prevDir == 'L':
				i -= 1
				prevDir = 'U'
			elif prevDir == 'U':
				j += 1
				prevDir = 'R'
			elif prevDir == 'D':
				j -= 1
				prevDir = 'L'
		elif orientation[pos] == 'L':
			if prevDir == 'R':
				i -= 1
				prevDir = 'U'
			elif prevDir == 'L':
				i += 1
				prevDir = 'D'
			elif prevDir == 'U':
				j -= 1
				prevDir = 'L'
			elif prevDir == 'D':
				j += 1
				prevDir = 'R'
		if grid[i][j] == '.':
			grid[i][j] = proteinMolecule[pos + 2]
			positionsOfAminoAcids[(i,j)] = (pos + 2, proteinMolecule[pos + 2])
			aminoAcidPositions[pos + 2] = (i,j)
		else:
			return 1.0

		
	#compute energy
	#for row in grid:
	#		print(row)

	i = numberOfMolecules
	j = numberOfMolecules

	freeEnergy = 0

	#print(aminoAcidPositions)
	#print(positionsOfAminoAcids)

	for pos in range(numberOfMolecules):
		#check neighbours
		#print(pos)
		i = aminoAcidPositions[pos][0]
		j = aminoAcidPositions[pos][1]	
		if proteinMolecule[pos] == 'H':

			if (i + 1, j) in positionsOfAminoAcids:
				if positionsOfAminoAcids[(i + 1, j)][1] == 'H' and not(positionsOfAminoAcids[(i + 1, j)][0] == (pos + 1) or positionsOfAminoAcids[(i + 1, j)][0] == (pos - 1)):     
					freeEnergy -= 1
					
			if (i - 1, j) in positionsOfAminoAcids:
				if positionsOfAminoAcids[(i - 1, j)][1] == 'H' and not(positionsOfAminoAcids[(i - 1, j)][0] == (pos + 1) or positionsOfAminoAcids[(i - 1, j)][0] == (pos - 1)):     
					freeEnergy -= 1
					
			if (i , j + 1) in positionsOfAminoAcids:
				if positionsOfAminoAcids[(i, j + 1)][1] == 'H' and not(positionsOfAminoAcids[(i, j + 1)][0] == (pos + 1) or positionsOfAminoAcids[(i, j + 1)][0] == (pos - 1)):     
					freeEnergy -= 1
					
			if (i , j - 1) in positionsOfAminoAcids:
				if positionsOfAminoAcids[(i, j - 1)][1] == 'H' and not(positionsOfAminoAcids[(i, j - 1)][0] == (pos + 1) or positionsOfAminoAcids[(i, j - 1)][0] == (pos - 1)):     
					freeEnergy -= 1
				
		
	freeEnergy = freeEnergy / 2.0
	#print(freeEnergy)
	return freeEnergy

#structureEnergy('CCCCCCCCCCCCCCCCCC','HPHPPHHPHPPHPHHPPHPH')