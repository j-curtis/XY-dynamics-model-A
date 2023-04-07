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

	return data



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

	prefix = "./L=100/vorticity_L=100_t=20000_T="
	labels = ["0.60","0.70","0.75","0.80","0.85","0.90","0.95","1.00","1.05","1.10","1.15","1.20","1.50"]

	### We downfold the time traces into correlation functions by sampling Nsample times and recording the time-delay correlation function for nts points
	### Ideally this would have an independent sampling from at least nts x Nsample points but we will make due with what we have and slightly oversample instead
	ntemps = len(labels)

	nts = 200
	Nsample = 200

	for j in range(ntemps):
		t0 = time.time()
		fname = prefix+labels[j]

		print("Loading data T="+labels[j])
		data = np.load(fname+".npy")

		correlation = correlationFFT(data,nts,Nsample)
		
		np.save(fname+"_correlation.npy",correlation)

		tf = time.time()
		print(str(tf-t0)+"s for T ="+labels[j]+"")

	nzpts = 50
	zpts = np.linspace(1.,50.,nzpts)

	L = 100

	filters = np.zeros((nzpts,L,L))
	for j in range(nzpts):
		filters[j,:,:] = genNVFilter(zpts[j],L)


	NVdata = np.zeros((ntemps,nzpts,nts))

	for j in range(ntemps):
		t0 = time.time()
		fname = prefix+labels[j]

		print("Loading correlation function T="+labels[j])
		data = np.load(fname+"_correlation.npy")

		for k in range(nzpts):
			NVdata[j,k,:] = np.sum(np.tensordot(np.ones(nts),filters[k,:,:],axes=0)*data,axis=(1,2)) 
	
		tf = time.time()
		print(str(tf-t0)+"s for T ="+labels[j])

	for j in range(ntemps):

		plt.plot(NVdata[j,10,:]/NVdata[j,10,0],label=r'$T = $'+labels[j])
		plt.title(r'z/$\xi_c$ = '+str(zpts[10]))
		plt.xlabel("t")

	plt.legend()
	plt.show()

	for j in range(ntemps):

		plt.plot(NVdata[j,20,:]/NVdata[j,20,0],label=r'$T = $'+labels[j])
		plt.title(r'z/$\xi_c$ = '+str(zpts[20]))
		plt.xlabel("t")

	plt.legend()
	plt.show()

	for j in range(ntemps):

		plt.plot(NVdata[j,40,:]/NVdata[j,40,0],label=r'$T = $'+labels[j])
		plt.title(r'z/$\xi_c$ = '+str(zpts[40]))
		plt.xlabel("t")

	plt.legend()
	plt.show()

	for j in range(ntemps):

		plt.plot(NVdata[j,:,0]/NVdata[j,0,0],label=r'$T = $'+labels[j])
		plt.title(r't = '+str(0))
		plt.xlabel("z")

	plt.legend()
	plt.show()

	for j in range(ntemps):

		plt.plot(NVdata[j,:,50]/NVdata[j,0,50],label=r'$T = $'+labels[j])
		plt.title(r't = '+str(50))
		plt.xlabel("z")

	plt.legend()
	plt.show()

	for j in range(ntemps):

		plt.plot(NVdata[j,:,100]/NVdata[j,0,100],label=r'$T = $'+labels[j])
		plt.title(r't = '+str(100))
		plt.xlabel("z")

	plt.legend()
	plt.show()

	for j in range(ntemps):

		plt.plot(NVdata[j,:,150]/NVdata[j,0,150],label=r'$T = $'+labels[j])
		plt.title(r't = '+str(150))
		plt.xlabel("z")

	plt.legend()
	plt.show()

	for j in range(ntemps):

		plt.plot(NVdata[j,:,-1]/NVdata[j,0,-1],label=r'$T = $'+labels[j])
		plt.title(r't = '+str(-1))
		plt.xlabel("z")

	plt.legend()
	plt.show()


if __name__ == "__main__":
	main()










