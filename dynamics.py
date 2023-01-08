###Jonathan Curtis
###Model A Langevin dynamics of XY model 
###12/09/2022

import numpy as np
from matplotlib import pyplot as plt 

L = 20### Lattice size -- LxL lattice

J = 1.### Phase stiffness in units of Kelvin
Gamma = 1.### Phase diffusivity (unitless) (can be taken to one, setting the size of time step)
T = 0.2### Temperature in units of Kelvin

dt = .0005### Time step (must be very small)
ntimes = 1000### Number of times steps we calculate

thetas = np.zeros((ntimes,L,L))
#thetas[0,:,:] = ( np.random.default_rng().uniform(-np.pi,np.pi,L*L) ).reshape((L,L))


### Equations of motion
### Gamma d theta_j /dt = J sum_nn sin(theta_j - theta_nn) + eta_j
### < eta_j(t_n)eta_j(t_m) > = 2T Gamma  delta(t_n-t_m)

for nt in range(1,ntimes):
	for k in range(L*L):
		x = k//L
		y = k % L
		thetas[nt,x,y] = thetas[nt-1,x,y] - J*dt*(
			np.sin( thetas[nt-1,x,y] - thetas[nt-1,(x+1)//L,y]) 
			+np.sin( thetas[nt-1,x,y] - thetas[nt-1,x-1,y]) 
			+np.sin( thetas[nt-1,x,y] - thetas[nt-1,x,(y+1)//L]) 
			+np.sin( thetas[nt-1,x,y] - thetas[nt-1,x,y-1])
			)/Gamma
		thetas[nt,x,y] += dt*( np.random.default_rng().normal(0.,2.*T*Gamma/dt)) /Gamma


### This will calculate the correlation function < e^{itheta(0) - i theta(x)}> for each time, summed over sites in the lattice

Gx = np.zeros((ntimes,L),dtype=complex)
for nt in rangE(ntimes):
	for nx in range(L):
		Gx[nt,nx] = np.mean( np.exp( 1.j*thetas[nt,0,:] - 1.j*thetas[nt,nx,:] ) )


nt = 0
while nt < ntimes:
	plt.plot( np.abs(Gx[nt,:])) 
	plt.show()


	nt += 100





