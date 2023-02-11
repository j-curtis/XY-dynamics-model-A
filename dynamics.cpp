#include <iostream>
#include <cmath>
#include <random>
#include <vector>
#include <fstream>

using namespace std;

const double pi = 3.14159265358979323846;

int main() {

	const size_t L = 300;
	const size_t ntimes = 3000;


	double J = 1.;
	double T = 2.2;
	double dt = .05/J;


	std::default_random_engine generator;
  	std::normal_distribution<double> normal(0.0,2.*T*dt);

   	std::vector<std::vector<std::vector<double> > > thetas(ntimes, std::vector<std::vector<double> >(L, std::vector<double>(L) ) );

	//Zero out the first time slice entry for thetas
	for(int nx = 0; nx < L; nx++){
		for(int ny = 0; ny <L; ++ny){
				thetas[0][nx][ny] = 0.;
		}
	}

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
	
	std::ofstream outfile;
	outfile.open("./thetas.csv");
	for(int nt =0; nt < ntimes; nt++){
		for(int nx =0; nx < L; nx++){
			for(int ny = 0; ny < L; ny++){
				outfile<<" "<<thetas[nt][nx][ny];
			}
		}
		outfile<<std::endl;
	}
	outfile.close();

	return 0;
}

