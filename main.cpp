#include <ctime>
#include <cmath>
#include <cstdio>

#include <fstream>
#include <iostream>

#include <typeinfo>
#include <vector>
#include <string>
#include <random>

#include <sys/stat.h>

#include "particle.h"
#include "vose.h"
#include "utils.h"

using namespace std;

template <typename T>
int sgn(T val)
{
	return (T(0) < val) - (val < T(0));
}

int main()
{
	//Open input file
	ifstream infile;

	infile.open("input_params.txt", ifstream::in);
	if (infile.fail())             //is it ok?
	{
		cout << "Input file did not open please check it\n";
		return 1;
	}

	//Determine how many particles there are
	unsigned int part_num = 0;
	//Get the timestep and duration
	double dt = 1, duration = 1;
	double stdev = 1;
	double mass = 1;
	int modulus = 1;
	double target_temp = 1;
	int Ng = 4;
	double length = 10;

	string input_type;
	string value;

	while ((infile >> input_type >> value))
	{
		if (input_type == "PART_NUM")
		{
			part_num = stoi(value);
		}
		else if (input_type == "STAN_DEV")
		{
			stdev = stod(value);
		}
		else if (input_type == "TIMESTEP")
		{
			dt = stod(value);
		}
		else if (input_type == "DURATION")
		{
			duration = stod(value);
		}
		else if (input_type == "LENGTH")
		{
			length = stod(value);
		}
		else if (input_type == "MODULUS")
		{
			modulus = stoi(value);
		}
		else if (input_type == "TEMP")
		{
			target_temp = stod(value);
		}
		else if (input_type == "GRIDSIZE")
		{
			Ng = stoi(value);
		}
	}
	printf("PART_NUM: %i\nSTDEV: %f\nMODULUS: %d\n", part_num, stdev, modulus);
	printf("TIMESTEP: %f\nDURATION: %f\ncell LENGTH:%f\n", dt, duration, length);
	printf("GRIDSIZE%d\n", Ng);
	infile.close();

	if (!dir_check("data/positions"))
	{
		const int dir_err = mkdir("./data/positions", S_IRWXU);
		if (-1 == dir_err)
		{
			printf("Error creating directory!\n");
			exit(1);
		}
	}
	else
	{
		system("rm -rf data/positions/*");
	}
	
	if (!dir_check("data/field"))
	{
		const int dir_err = mkdir("./data/field", S_IRWXU);
		if (-1 == dir_err)
		{
			printf("Error creating directory!\n");
			exit(1);
		}
	}
	else
	{
		system("rm -rf data/field/*");
	}

	if (!dir_check("data/charge"))
	{
		const int dir_err = mkdir("./data/charge", S_IRWXU);
		if (-1 == dir_err)
		{
			printf("Error creating directory!\n");
			exit(1);
		}
	}
	else
	{
		system("rm -rf data/charge/*");
	}

	if (!dir_check("data/potential"))
	{
		const int dir_err = mkdir("./data/potential", S_IRWXU);
		if (-1 == dir_err)
		{
			printf("Error creating directory!\n");
			exit(1);
		}
	}
	else
	{
		system("rm -rf data/potential/*");
	}
	
	//double temp;
	double* charge = new double[Ng];
	double* potential = new double[Ng];
	double* EField = new double[Ng];

	for (int i = 0; i < Ng; i++)
	{
		charge[i] = 0;
		potential[i] = 0;
		EField[i] = 0;
	}

	vector<particle> particles = initialize_system(part_num, mass, length);

	write_state(particles, charge, potential, EField, Ng, 0);

	//Forward Euler to start Leapfrog
	double dx = length / Ng;
	calc_charge(Ng, dx, charge, particles);
	calc_potential(Ng, dx, potential, charge);
	calc_EField(Ng, dx, EField, potential);
	

	for (auto &part: particles)
	{
		for (int i = 0; i < Ng; i++)
		{
			if (fabs(i * dx - part.pos[0]) <= dx / 2)
			{
				part.mom[0] += dt / 2 * EField[i];
				continue;         // To ensure only the closest EField is used.
			}
		}
	}

	write_state(particles, charge, potential, EField, Ng, 1);


	//Run the simulation
	unsigned int i = 2;
	//double scale = 1;
	for (double t = 0; t <= duration; t += dt)
	{
		update_state(particles, Ng, dt, length, charge, potential, EField);

		/*
		   // Thermostat for first half
		   temp = temperature(particles);
		   if (t<0.5*duration)
		   {
		    scale = sqrt(1+(target_temp/temp - 1)/10);
		    scale_mom(particles, scale);
		   }*/

		if ((i % modulus) == 0)
		{
			printf("i: %d\n", i);
			write_state(particles, charge, potential, EField, Ng, i);
		}
		i++;
	}

	delete[] charge;
	delete[] potential;
	delete[] EField;

	return 0;
}
