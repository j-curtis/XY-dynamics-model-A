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



def main():

	labels = ["0.3","0.8","0.9","1.0","3.0"]
	fnames = ["vorticity_L=350_t=10000s_T="+l+"J" for l in labels]
	ntemps = len(labels)
	temps = np.array([float(l) for l in labels])

	L = 350

	staticFFTRMS = np.zeros((ntemps,L,L))

	for i in range(ntemps):
		data = np.load(fnames[i]+".npy")
		fft = np.fft.fftn(data[-1000:,:,:],axes=[1,2])
		meanFFT = np.mean(fft,axis=0)
		staticFFTRMS[i,:,:] = np.real(np.mean(np.abs(fft)**2, axis= 0 ) - np.abs(meanFFT)**2)

	for i in range(ntemps):
		plt.plot(staticFFTRMS[i,0,:],label=r'$T/J = $'+labels[i])
	plt.legend()
	plt.show()

	#np.save("staticFFT.npy",staticFFT)
	
	quit()

	#np.save(fname+"_FFT.npy",dataFFT)


	nts = 1000
	Nsample  = 30

	plt.plot(NVNoise(dataFFT,z=.1,nts = nts,Nsample=Nsample),label=r'z = 0.1')
	plt.plot(NVNoise(dataFFT,z=.5,nts = nts,Nsample=Nsample),label=r'z = 0.5')
	plt.plot(NVNoise(dataFFT,z=1.,nts = nts,Nsample=Nsample),label=r'z = 1.0')
	plt.plot(NVNoise(dataFFT,z=5.,nts = nts,Nsample=Nsample),label=r'z = 5.0')
	plt.yscale('log')
	plt.xscale('log')
	plt.legend()
	plt.show()

	quit()

	for nt in range(20):
		dataFFT = np.fft.fftn(data[nt,:,:],s=(L,L))
		filt = np.real(np.fft.ifftn(filterfunc*dataFFT))

	
		plt.imshow(filt)
		plt.colorbar()
		plt.show()




if __name__ == "__main__":
	main()