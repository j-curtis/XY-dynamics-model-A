###Jonathan Curtis
###Model A Langevin dynamics of XY model 
###12/09/2022

import numpy as np
from matplotlib import pyplot as plt 

L = 40### Lattice size -- LxL lattice

J = 1.### Phase stiffness in units of Kelvin
T = .7### Temperature in units of Kelvin

dt = .001### Time step (must be very small) 
ntimes = 5000### Number of times steps we calculate

thetas = np.zeros((ntimes,L,L))
#thetas[0,:,:] = ( np.random.default_rng().uniform(-np.pi,np.pi,L*L) ).reshape((L,L))


### Equations of motion
### Gamma d theta_j /dt = J sum_nn sin(theta_j - theta_nn) + eta_j
### < eta_j(t_n)eta_j(t_m) > = 2T Gamma  delta(t_n-t_m)

for nt in range(1,ntimes):
	for k in range(L*L):
		x = k//L
		y = k % L
		###thetas[nt,x,y] = thetas[nt-1,x,y] - J*dt*(
		###	np.sin( thetas[nt-1,x,y] - thetas[nt-1,(x+1)//L,y]) 
		###	+np.sin( thetas[nt-1,x,y] - thetas[nt-1,x-1,y]) 
		###	+np.sin( thetas[nt-1,x,y] - thetas[nt-1,x,(y+1)//L]) 
		###	+np.sin( thetas[nt-1,x,y] - thetas[nt-1,x,y-1])
		###	)
		thetas[nt,x,y] = thetas[nt-1,x,y]- J * dt *np.sin(thetas[nt-1,x,y])
		thetas[nt,x,y] +=  np.random.default_rng().normal(0.,2.*T*dt)


### This will calculate the correlation function < e^{itheta(0) - i theta(x)}> for each time, summed over sites in the lattice

Gx = np.zeros((ntimes,L+1),dtype=complex)
for nt in range(ntimes):
	for nx in range(L+1):
		Gx[nt,nx] = np.mean( np.exp( 1.j*(thetas[nt,nx % L,:] - thetas[nt,0,:]) ) )


nt = 0
while nt < ntimes:
	plt.plot( np.abs(Gx[nt,:])) 
	plt.show()


	nt += 500

plt.plot(np.abs(np.mean(Gx,axis=0)))
plt.show()



