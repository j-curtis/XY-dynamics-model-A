###Jonathan Curtis
###Processes and analyzes csv files to extract fourier transform noise profiles
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



def noise_vs_z(data):
	"""Computes teh dynamical noise vs z for the nv center"""
	L = data.shape[1]

	dataFFT = np.fft.rfftn(data,s=data.shape)

	zs = np.range(L)

	### filter function for a particular z is given by 
	### e^(-qz)/q for each voricity v(q), we then square and sum over q 

def genNVFilter(z,L):
	"""Generates momentum space filter for probe at depth z for momentum grid of LxL system"""
	###We first need to generate an appropriate momentum space mesh
	###According to the numpy conventions, for RFFT with output size LxL the momenta ranges as 0, 2pi/L, 2*2pi/L, ... (L-1)*2pi/L
	qs = np.linspace(0.,2.*np.pi,L)
	filterfunc = np.zeros((L,L))
	for j in range(L):
		for k in range(L):
			q = np.sqrt(qs[j]**2 + qs[k]**2)
			
			filterfunc[j,k] = np.exp(-z*q)/(2.*q + .0001) ### For q -> 0 we shift slightly, it should be ultimately suppressed anyways

	return filterfunc

def NVNoise(dataFFT,z,nts,Nsample):
	"""Compute the NV noise <Bz(t) Bz(0)> for a given depth z"""
	L = dataFFT.shape[1]

	filtfunc = genNVFilter(z,L)


	###Now we pick N time points t)i to average over (we compute for t_i -> t_i+nt) to obtain an ensemble average
	noise = np.zeros(nts)

	for N in range(Nsample):
		t0 = np.random.randint(0,dataFFT.shape[0]-nts)
		for nt in range(nts):
			noise[nt] += np.sum(dataFFT[t0,:,:]*filtfunc*filtfunc*np.conjugate(dataFFT[t0+nt,:,:]))

	return np.real(noise)

def correlationFFT(data,nts,Nsample):
	"""Compute the vorticity autocorrelation function < n(q,t) n^*(q,0) > """
	L = data.shape[1]
	correlation = np.zeros((nts,L,L),dtype=complex)
	meanFFT = np.zeros((L,L),dtype=complex)

	for N in range(Nsample):
		t0 = np.random.randint(0,data.shape[0]-nts)
		meanFFT += np.fft.fftn(data[t0,:,:],axes=[0,1])/float(Nsample)

	for N in range(Nsample):
		t0 = np.random.randint(0,data.shape[0]-nts)
		dataFFT0 = (np.fft.fftn(data[t0,:,:],axes=[0,1]) - meanFFT)
		correlation[0,:,:] += dataFFT0*np.conjugate(dataFFT0)/float(Nsample)

		for nt in range(1,nts):
			dataFFT = (np.fft.fftn(data[t0+nt,:,:],axes=[0,1]) - meanFFT)
			correlation[nt,:,:] += dataFFT*np.conjugate(dataFFT0)/float(Nsample)

	return correlation


def main():

	labels = ["0.8","0.9","1.0","3.0"]
	fnames = ["vorticity_L=350_t=10000s_T="+l+"J" for l in labels]
	ntemps = len(labels)
	temps = np.array([float(l) for l in labels])

	L = 350
	Lsave = 50
	nts = 50
	Nsample = 300
	z0 = 20
	zf = 100.
	nzs = 50

	zpts = np.linspace(z0,zf,nzs)

	FFTs = np.zeros((ntemps,nts,Lsave,Lsave),dtype=complex)
	noise = np.zeros((ntemps,nzs,nts))
	for i in range(ntemps):
		FFTs[i,:,:,:] = np.load(fnames[i]+"_correlationFFT.npy")

		for j in range(nzs):
			filterfunc = (genNVFilter(zpts[j],L)[:Lsave,:Lsave])**2

			noise[i,j,:] = np.sum(np.tensordot(np.ones(nts),filterfunc[1:,1:],axes=0)*FFTs[i,:,1:,1:], axis=(1,2) )

		plt.plot(noise[i,10,:])
		plt.xlabel(r'$t$')
		plt.yscale('log')
		plt.show()
		
		plt.plot(zpts,noise[i,:,10])
		plt.xlabel(r'$z [\xi_c]$')
		plt.yscale('log')
		plt.show()



if __name__ == "__main__":
	main()










