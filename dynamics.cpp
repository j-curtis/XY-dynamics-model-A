#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <fstream>
#include <ctime>
#include <complex>

using namespace std;

const double pi = 3.14159265358979323846;

double* OPav(size_t L, size_t ntimes, double T, double J){

	double dt = .05/J;

	//Initialize RNG
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  	std::default_random_engine generator (seed);
  	std::normal_distribution<double> normal(0.0,sqrt(2.*T*dt));

  	//Array of theta(x,t) values
   	std::vector<std::vector<std::vector<double> > > thetas(ntimes, std::vector<std::vector<double> >(L, std::vector<double>(L) ) );

   	//We first run a burn run to thermalize the system
	for(int nt = 1; nt < ntimes; ++nt){
		for(int nx = 0; nx < L; nx++){
			for(int ny = 0; ny <L; ny++){

				//Implement the dynamics -- first the dE/dtheta terms
				thetas[nt][nx][ny] = thetas[nt-1][nx][ny] - J*dt*(
				 sin( thetas[nt-1][nx][ny] - thetas[nt-1][(nx+1)%L][ny] )
					+sin( thetas[nt-1][nx][ny] - thetas[nt-1][(nx-1)%L][ny] ) 
					+sin( thetas[nt-1][nx][ny] - thetas[nt-1][nx][(ny+1)%L] )
					+sin( thetas[nt-1][nx][ny] - thetas[nt-1][nx][(ny-1)%L] ) );

				//Now add a random noise term
				thetas[nt][nx][ny] += normal(generator);

				//Now bring back to interval [-pi,pi]
				thetas[nt][nx][ny] = std::fmod(thetas[nt][nx][ny] , 2.*pi);

			}
		}
	}



	//Now we calculate the space-time average of the order parameter 

	static double sum[2];
	sum[0] = 0.;
	sum[1] = 0.;

	for(int nt = 0; nt < ntimes; ++nt){
		sum[0] += (cos(thetas[nt][0][0]))/(double(ntimes));
		sum[1] += (sin(thetas[nt][0][0]))/(double(ntimes));

	}

	return sum;
}


int main() {

	const size_t L = 300;
	const size_t ntimes = 5000;


	double J = 1.;
	const size_t ntemps = 12;
	double temps[ntemps];
	for(int i =0; i < ntemps; i++){
		temps[i] = 1.8*double(i+1)*J/double(ntemps);
	}

	double dt = .05/J;
   	int t0 = std::time(NULL);

	double mag[ntemps][2]; 
	for(int i =0; i < ntemps; ++i){
		double * tmp = OPav(L,ntimes,temps[i],J);
		mag[i][0] = tmp[0];
		mag[i][1] = tmp[1];
	}
	int tf = std::time(NULL);

	std::cout<<(tf-t0)<<" seconds elapsed"<<std::endl;

	std::ofstream outfile;
	outfile.open("./M_vs_T.csv");
	for(int i =0; i < ntemps; i++){
		outfile<<temps[i]<<", "<<mag[i][0]<<", "<<mag[i][1]<<std::endl;
	}
	outfile.close();
	/*

	std::default_random_engine generator;
  	std::normal_distribution<double> normal(0.0,2.*T*dt);

   	std::vector<std::vector<std::vector<double> > > thetas(ntimes, std::vector<std::vector<double> >(L, std::vector<double>(L) ) );

   	int t0 = std::time(NULL);

   	//We first run a burn run to thermalize the system
	for(int nt = 1; nt < ntimes; ++nt){
		for(int nx = 0; nx < L; nx++){
			for(int ny = 0; ny <L; ny++){

				//Implement the dynamics -- first the dE/dtheta terms
				thetas[nt][nx][ny] = thetas[nt-1][nx][ny] - J*dt*(
				 sin( thetas[nt-1][nx][ny] - thetas[nt-1][(nx+1)%L][ny] )
					+sin( thetas[nt-1][nx][ny] - thetas[nt-1][(nx-1)%L][ny] ) 
					+sin( thetas[nt-1][nx][ny] - thetas[nt-1][nx][(ny+1)%L] )
					+sin( thetas[nt-1][nx][ny] - thetas[nt-1][nx][(ny-1)%L] ) );

				//Now add a random noise term
				thetas[nt][nx][ny] += normal(generator);

				//Now bring back to interval [-pi,pi]
				thetas[nt][nx][ny] = std::fmod(thetas[nt][nx][ny] , 2.*pi) - pi;

			}
		}
	}

	//Initialize the array back with the thermalized state
	for(int nx = 0; nx < L; nx++){
		for(int ny = 0; ny <L; ny++){
			thetas[0][nx][ny] = thetas[ntimes-1][nx][ny];
		}
	}

	//Repeat; sample dynamics from this batch
	for(int nt = 1; nt < ntimes; ++nt){
		for(int nx = 0; nx < L; nx++){
			for(int ny = 0; ny <L; ny++){

				//Implement the dynamics -- first the dE/dtheta terms
				thetas[nt][nx][ny] = thetas[nt-1][nx][ny] - J*dt*(
				 sin( thetas[nt-1][nx][ny] - thetas[nt-1][(nx+1)%L][ny] )
					+sin( thetas[nt-1][nx][ny] - thetas[nt-1][(nx-1)%L][ny] ) 
					+sin( thetas[nt-1][nx][ny] - thetas[nt-1][nx][(ny+1)%L] )
					+sin( thetas[nt-1][nx][ny] - thetas[nt-1][nx][(ny-1)%L] ) );

				//Now add a random noise term
				thetas[nt][nx][ny] += normal(generator);

				//Now bring back to interval [-pi,pi]
				thetas[nt][nx][ny] = std::fmod(thetas[nt][nx][ny] , 2.*pi) - pi;

			}
		}
	}

	int tf = std::time(NULL);

	std::cout<<(tf-t0)<<" seconds elapsed"<<std::endl;

	//We only output sampled thetas from every few steps to increase ergodicity and efficiency
	int samplestep = 500;

	std::ofstream outfile;
	outfile.open("./thetas.csv");
	for(int nt = 0; nt < ntimes; nt += samplestep){
		for(int nx =0; nx < L; nx++){
			for(int ny = 0; ny < L; ny++){
				outfile<<" "<<thetas[nt][nx][ny];
			}
		}
		outfile<<std::endl;
	}
	outfile.close();
	*/

	return 0;
}

