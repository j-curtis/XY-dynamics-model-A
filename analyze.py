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


def plotthetas(fname):

	#L = 200### Lattice size -- LxL lattice
	J = 1.### Phase stiffness in units of Kelvin
	T = 1.*J ### Literature says BKT transition is at approximately .89 J 
	#ntimes = 20000
	dt = .05*J

	data = np.load(fname+".npy")
	nt = data.shape[0]
	L = data.shape[1]

	Gxmean = np.zeros(L,dtype=complex)
	#Gxvar = np.zeros(L,dtype = complex)

	Gxmean = np.mean(np.exp(1.j*(data[:,0,:] - np.outer(data[:,0,0],np.ones(L) ) ) )  ,axis=0)

	#plt.plot(np.abs(Gxmean))
	#plt.show()

	Gtmean = np.mean(np.exp(1.j*(data[:,:,:] - np.tensordot(np.ones(nt),data[0,:,:],axes=0) ) )  ,axis=(-1,-2))
	plt.plot(np.abs(Gtmean))
	plt.show()

def plotMs(fname):
	L = 300
	J = 1.
	dt = .05 

	data = np.genfromtxt(fname+".csv",delimiter=", ")
	ntemps = data.shape[0]
	temps = data[:,0]
	m = data[:,1]**2 + data[:,2]**2
	plt.plot(temps,m)
	plt.show()

def plotvorticity(fname):
	J = 1.### Phase stiffness in units of Kelvin
	T = 1.*J ### Literature says BKT transition is at approximately .89 J 
	#ntimes = 20000
	dt = .05*J

	data = np.load(fname+".npy")
	nt = data.shape[0]
	L = data.shape[1]

	for i in range(30):
		plt.imshow(data[i*20,:30,:30])
		plt.show()


def computeFFT(fname):
	"""Computes the time-dependence spatial FFT of the vorticity in order to extract the momentum space autocorrelation functions"""
	data = np.load(fname+".npy")
	nt = data.shape[0]
	L = data.shape[1]

	### We now need to take the FFT for each time

	data_kspace = np.mean(np.fft.fft2(data,axes=(-1,-2)),axis=0)

	plt.imshow(data[0,:,:])
	plt.colorbar()
	plt.show()
	plt.imshow(np.abs(data_kspace)**2)
	plt.colorbar()
	plt.show()

	

def main():
	fname = "vorticity_10000s_T=1.0J"
	#csv2npy(fname)
	computeFFT(fname)

if __name__ == "__main__":
	main()







