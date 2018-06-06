def sphere(x):
	z = 0.0
	for xi in x:
		z += xi**2
	return z

def rosenbrock(x):
	z = 0.0
	for i in range(len(x) - 1):
		z += 100*((x[i]**2 - x[i + 1])**2) + (x[i] - 1)**2
	return z

def ackley(x):
	z = 20.0 + math.e
	z1 = 0.0
	z2 = 0.0
	
	for xi in x:
		z1 += xi**2
		z2 += math.cos(2*math.pi*xi)

	z1 = z1/len(x)
	z2 = z2/len(x)

	z1 = math.sqrt(z1)
	z1 = (-0.2)*z1
	z1 = (-20)*math.exp(z1)

	z2 =  (-1)*math.exp(z2)

	z = z + z1 + z2

	return z


def schwefel(x):
	alpha = 418.9829
	fitness = 0
	for i in range(len(x)):
		fitness -= x[i]*math.sin(math.sqrt(math.fabs(x[i])))

	return float(fitness) +  alpha*len(x)


def rastrigin(x):
	z = 0.0
	for xi in x:
		z += xi**2 - 10*math.cos(2*math.pi*xi) + 10.0

	return z

def griewank(x):
	z = 1.0
	z1 = 0.0
	z2 = 0.0

	for i in range(len(x)):
		z1 += (x[i]**2)/4000.0
		z2 *= math.cos((x[i])/math.sqrt(i + 1))

	z += z1 - z2

	return z