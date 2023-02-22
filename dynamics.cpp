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

	const size_t L = 200;
	const size_t ntimes = 5000;

	//Measure energy in units of J
	double J = 1.;
	double dt = .05/J;

	//Different temperatures
	/*
	const size_t ntemps = 10;
	double temps[ntemps];
	for(int i =0; i < ntemps; i++){
		temps[i] = (.5  + double(i+1)/double(ntemps))*J;
	}
	*/
	//Set temperature
	double T = 4.*J;

   	int t0 = std::time(NULL);
   	/*
   	//Run loop to collect M(T)
 	double mag[ntemps][2]; 
	for(int i =0; i < ntemps; ++i){
		double * tmp = OPav(L,ntimes,temps[i],J);
		mag[i][0] = tmp[0];
		mag[i][1] = tmp[1];
	}
	*/
	
	//Loop to generate full spacetime data for theta(x,t)

	//Initialize RNG
	unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  	std::default_random_engine generator (seed);
  	std::normal_distribution<double> normal(0.0,sqrt(2.*T*dt));

  	//Array of theta(x,t) values
   	std::vector<std::vector<std::vector<double> > > thetas(ntimes, std::vector<std::vector<double> >(L, std::vector<double>(L) ) );
	
	//Burn loop
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

	//Initialize the array back with the thermalized state
	for(int nx = 0; nx < L; nx++){
		for(int ny = 0; ny <L; ny++){
			thetas[0][nx][ny] = thetas[ntimes-1][nx][ny];
		}
	}

	//Sample loop
	//Now we collect the data and store as a csv file
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

	/*
	//We now compute the correlation function < e^i(theta(x,t)-theta(0))>
	const size_t ntd = 2000; //Number of time steps we evaluate the correlation delay function for 
   	std::vector<std::vector<double> > costx(ntd, std::vector<double>(L) );

   	for(int nt =0; nt < ntd; nt++){
   		for(int nx = 0; nx < L; nx++){
   			for(int i =0; i < L ; i ++){
   				costx[nt][nx] += cos(thetas[nt][nx][i] - thetas[nt][nx][0])/double(L);
   			}
   		}
   	}
   	*/


	
	//We now compute the time-traces of the vorticity
	//We define the vorticity for a site r such that v(r) = (theta(r+x)-theta(r))mod 2pi + (theta(r+x+y)-theta(r+x))mod 2pi  +(theta(r+y)-theta(r+x+y))mod 2pi  +(theta(r)-theta(r+y))mod 2pi 
   	std::vector<std::vector<std::vector<double> > > vort(ntimes, std::vector<std::vector<double> >(L, std::vector<double>(L) ) );

   	for(int t =0; t < ntimes; t++){
   		for(int x = 0; x < L; x++){
   			for(int y =0; y< L ; y++){
   				vort[t][x][y] = std::fmod(thetas[t][(x+1)%L][y] - thetas[t][x][y],2.*pi);
   				vort[t][x][y] += std::fmod(thetas[t][(x+1)%L][(y+1)%L] - thetas[t][(x+1)%L][y],2.*pi);
   				vort[t][x][y] += std::fmod(thetas[t][x][(y+1)%L] - thetas[t][(x+1)%L][(y+1)%L],2.*pi);
   				vort[t][x][y] += std::fmod(thetas[t][x][y] - thetas[t][x][(y+1)%L],2.*pi);

   			}
   		}
   	}

	//Save vorticity to csv file (large data set)
	std::ofstream outfile;
	outfile.open("./vorticity.csv");
	for(int t = 0; t < ntimes; t++){
		for(int x =0; x < L; x++){
			for(int y = 0; y < L; y++){
				outfile<<" "<<vort[t][x][y];
			}
		}
		outfile<<std::endl;
	}
	outfile.close();

   	/*
	//Save thermalized theta(t,x) to csv file (large data set)
	std::ofstream outfile;
	outfile.open("./thetas.csv");
	for(int nt = 0; nt < ntimes; nt++){
		for(int nx =0; nx < L; nx++){
			for(int ny = 0; ny < L; ny++){
				outfile<<" "<<thetas[nt][nx][ny];
			}
		}
		outfile<<std::endl;
	}
	outfile.close();
	*/
	
	//Save M(T) to csv file
	/*
	std::ofstream outfile;
	outfile.open("./M_vs_T.csv");
	for(int i =0; i < ntemps; i++){
		outfile<<temps[i]<<", "<<mag[i][0]<<", "<<mag[i][1]<<std::endl;
	}
	outfile.close();
	*/
	int tf = std::time(NULL);
	std::cout<<(tf-t0)<<" seconds elapsed"<<std::endl;

	
	return 0;
}

