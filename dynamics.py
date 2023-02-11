###Jonathan Curtis
###Model A Langevin dynamics of XY model 
###12/09/2022

import numpy as np
from matplotlib import pyplot as plt 
import time 

def genThetas(L,T,ntimes,J=1.,dt = .05):
	"""Generates a sequence of angles for a given system size, temperature, number of time steps, Josephson coupling (default is 1.) and time step size (default is 5% J)"""
	### Uses periodic boundary conditions
	### defaul time step is 5% of J

	thetas = np.zeros((ntimes,L,L))

	for nt in range(1,ntimes):
		for k in range(L*L):
			x = k//L
			y = k % L
			thetas[nt,x,y] = thetas[nt-1,x,y] - J*dt*(
				np.sin( thetas[nt-1,x,y] - thetas[nt-1,(x+1)//L,y]) 
				+np.sin( thetas[nt-1,x,y] - thetas[nt-1,x-1,y]) 
				+np.sin( thetas[nt-1,x,y] - thetas[nt-1,x,(y+1)//L]) 
				+np.sin( thetas[nt-1,x,y] - thetas[nt-1,x,y-1])
				)

		thetas[nt,:,:] +=  np.random.normal(0.,2.*T*dt,size=(L,L))
		### We mod to range(-pi,pi)
		thetas[nt,:,:] = thetas[nt,:,:] % (2.*np.pi) - np.pi*np.ones((L,L))


	return thetas

def calcOP(thetas):
	"""Calculates the spatial average of order parameter over space"""
	s = thetas.shape
	ntimes = s[0]
	OP = np.zeros(ntimes,dtype=complex)
	OP = np.mean(np.exp(1.j*thetas),axis=(-1,-2))

	return OP


def main():

	L = 250### Lattice size -- LxL lattice
	J = 1.### Phase stiffness in units of Kelvin
	T = .7*.89*J ### Literature says BKT transition is at approximately .89 J 

	#dt = .05### Time step (must be very small) 
	nburn = 5000### Time steps we burn initially to equilibrate
	ntimes = 3000### Number of times steps we calculate and measure for

	ti = time.time()

	thetas = genThetas(L,T,nburn+ntimes)

	tf = time.time()
	print("Elapsed total time: ",tf-ti,"s")

	np.save("thetas.npy",thetas)



if __name__ == "__main__":
	main()

