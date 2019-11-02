#include "random.h"

void RandomNumbers::uniform_double (std::vector<double> &tableau_uni, double lower, double upper){
	std::uniform_real_distribution <> duniform(lower, upper);
	for (double n=0; n<tableau_uni.size(); ++n){
		tableau_uni [n] = duniform(rng);
	}
}

void RandomNumbers::normal(std::vector<double> &tableau_norm, double mean, double sd){
	std::normal_distribution <> dnormal (mean, sd);
	for (double n=0; n<tableau_norm.size(); ++n){
		tableau_norm[n] = dnormal(rng);
	}
}

void RandomNumbers::poisson(std::vector<int> &tableau_pois, double mean){
	std::poisson_distribution <> dpoisson(mean);
	for (double n=0; n<tableau_pois.size(); ++n){
		tableau_pois[n] = dpoisson(rng);
	}
}

double RandomNumbers::uniform_double (double lower, double upper){
	std::uniform_real_distribution <> duniform(lower, upper);
	return duniform(rng);
}

double RandomNumbers::normal(double mean, double sd){
	std::normal_distribution <> dnormal (mean, sd);
	return dnormal(rng);
}

int RandomNumbers::poisson(double mean){
	std::poisson_distribution <> dpoisson(mean);
	return dpoisson (rng);
}


RandomNumbers::RandomNumbers(unsigned long int s) 
{
    if (s == 0) {
        std::random_device rd;
        s = rd();
    }
    rng = std::mt19937(s);
}
	



