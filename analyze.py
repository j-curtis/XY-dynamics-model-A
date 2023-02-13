###Jonathan Curtis
###Model A Langevin dynamics of XY model 
###12/09/2022

import numpy as np
from matplotlib import pyplot as plt 
import time 

def csv2npy(fname):
	data = np.genfromtxt(fname+".csv",delimiter=" ")
	nt = data.shape[0]
	L = int(np.sqrt(data.shape[1]))
	data = data.reshape(nt,L,L)
	np.save(fname+".npy",data)


def plot(fname):

	L = 500### Lattice size -- LxL lattice
	J = 1.### Phase stiffness in units of Kelvin
	T = .5*J ### Literature says BKT transition is at approximately .89 J 
	#ntimes = 5000
	dt = .05

	data = np.load(fname+".npy")

	Gxmean = np.zeros(L,dtype=complex)
	Gxvar = np.zeros(L,dtype = complex)

	for n in range(L):
		Gxmean[n] = np.mean(np.exp(1.j*(data[:,n,0] - data[:,0,0] ) )  ,axis=0)
		Gxvar[n] = np.var(np.exp(1.j*(data[:,n,0] - data[:,0,0] ) )  ,axis=0)

	#Gx = np.mean(np.exp(1.j*(data[:,:,0] - np.outer(data[:,0,0],np.ones(L) ) ) ) ,axis=0)

	plt.plot(np.abs(Gxmean)**2)
	plt.show()

def main():
	fname = "thetas"
	csv2npy(fname)
	plot(fname)

if __name__ == "__main__":
	main()

