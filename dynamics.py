###Jonathan Curtis
###Model A Langevin dynamics of XY model 
###12/09/2022

import numpy as np
from matplotlib import pyplot as plt 
import time 

def genThetas(L,T,ntimes,J=1.,dt = .05):
	"""Generates a sequence of angles for a given system size, temperature, number of time steps, Josephson coupling (default is 1.) and time step size (default is 5% J)"""
	### Returns thetas and elapsed time in seconds
	### Uses periodic boundary conditions
	### defaul time step is 5% of J

	thetas = np.zeros((ntimes,L,L))

	t0 = time.time()
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

	t1 = time.time()

	return thetas, t1-t0

def calcOP2(thetas):
	"""Calculates the system RMS of order parameter over time"""
	s = thetas.shape
	ntimes = s[0]
	OP = np.zeros(ntimes,dtype=complex)
	OP = np.mean(np.exp(1.j*thetas),axis=(-1,-2))

	return OP

L = 200### Lattice size -- LxL lattice

J = 1.### Phase stiffness in units of Kelvin
#temps = np.linspace(0.,np.pi*J,10) ### Temperatures we study, we perform a sweep
#ntemps = len(temps)
T = 3.5*J ### Literature says BKT transition is at approximately .89 J 

dt = .05### Time step (must be very small) 
ntimes = 1000### Number of times steps we calculate

thetas,et = genThetas(L,T,ntimes)

print("Elapsed time: ",et,"s")

op = calcOP2(thetas)

plt.plot(np.abs(op)**2)
plt.show()

quit()

OP = np.zeros(ntimes,dtype=complex)
OP = np.mean(np.exp(1.j*thetas),axis=(-1,-2))

plt.plot(np.abs(OP)**2)
plt.show()

quit()


thetas = np.zeros((ntimes,L,L))
#thetas[0,:,:] = ( np.random.default_rng().uniform(-np.pi,np.pi,L*L) ).reshape((L,L))

t0 = time.time()


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

tf = time.time()
print("Elapsed time: ",tf-t0,"s")

OP = np.zeros(ntimes,dtype=complex)
OP = np.mean(np.exp(1.j*thetas),axis=(-1,-2))

plt.plot(np.abs(OP)**2)
plt.show()



