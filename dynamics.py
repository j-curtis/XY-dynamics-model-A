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

	L = 200### Lattice size -- LxL lattice
	J = 1.### Phase stiffness in units of Kelvin
	nTs = 5
	Ts = np.linspace(0.*J,3.*J,nTs) ### Literature says BKT transition is at approximately .89 J 

	#dt = .05### Time step (must be very small) 
	nburn = 5000### Time steps we burn initially to equilibrate
	ntimes = 5000### Number of times steps we calculate and measure for

	ti = time.time()

#	opMean = np.zeros(nTs) ### Statistical average of local order parameter
	GxMean = np.zeros((nTs,L)) ### Statistical average of order parameter correlation function

	for n in range(nTs):

		t1 = time.time()
		thetas = genThetas(L,Ts[n],nburn+ntimes)
	#	opMean[n] = np.mean( np.exp(1.j*thetas[nburn:,0,0]) )
		GxMean[n,:] = np.mean( np.exp(1.j*(thetas[nburn:,:,0] - np.outer(thetas[nburn:,0,0],np.ones(L)) ) ), axis=0 )
		t2 = time.time()
		print("Temperature : ",n+1,"/",nTs,", Run time: ",t2-t1,"s")


	tf = time.time()
	print("Elapsed total time: ",tf-ti,"s")


	for n in range(nTs):
		plt.plot(np.abs(GxMean[n,:])**2)
	plt.show()





if __name__ == "__main__":
	main()

