###Jonathan Curtis
###Model A Langevin dynamics of XY model 
###12/09/2022

import numpy as np
from matplotlib import pyplot as plt 
import time 


def main():

	L = 250### Lattice size -- LxL lattice
	J = 1.### Phase stiffness in units of Kelvin
	T = .7*.89*J ### Literature says BKT transition is at approximately .89 J 

	#dt = .05### Time step (must be very small) 
	nburn = 5000### Time steps we burn initially to equilibrate
	ntimes = 3000### Number of times steps we calculate and measure for

	tpts = 500 ###We sample correlation functions for t > nburn and t < tpts 
	xpts = 100

	thetas = np.load("thetas.npy")

	Gtx = np.zeros((tpts,xpts),dtype=complex)

	samplepts = 200 ### We randomly select this many points to sample correlation function

	### We pick a spacetime separation of nt and nx
	for nt in range(tpts):
		for nx in range(xpts):

			samples = np.zeros(samplepts,dtype=complex)

			for k in range(samplepts):
				pt = np.random.randint(L,size=2)

				samples[k] = np.exp( 1.j* ( thetas[(nburn+nt),(pt[0] + nx)//L,pt[1]]-thetas[nburn,pt[0],pt[1]] ) ) 

			Gtx[nt,nx] = np.mean(samples)

	plt.plot(np.abs(Gtx[0,:])**2)
	plt.plot(np.abs(Gtx[100,:])**2)
	plt.plot(np.abs(Gtx[200,:])**2)
	plt.plot(np.abs(Gtx[400,:])**2)
	plt.show()


	plt.plot(np.abs(Gtx[:,0])**2)
	plt.plot(np.abs(Gtx[:,50])**2)
	plt.plot(np.abs(Gtx[:,90])**2)
	plt.show()

	plt.imshow(np.abs(Gtx)**2 )
	plt.colorbar()
	plt.show()




if __name__ == "__main__":
	main()

