#include "utils.h"

bool dir_check(char const *pathname)
{
	struct stat info;

	if ( stat( pathname, &info ) != 0 )
	{
		printf( "cannot access %s\n", pathname );
		return false;
	}
	else if ( info.st_mode & S_IFDIR )
	{   // S_ISDIR() doesn't exist on my windows
		printf( "%s exists\n", pathname );
		return true;
	}
	else
	{
		printf( "%s is no directory\n", pathname );
		return false;
	}
}

vector<particle> initialize_system(int part_num, double mass, double length)
{
	// Initialize particles
	int grid_size = 1000;

	//Create vector to store Particles
	vector<particle> particles;
	for (unsigned int i = 0; i < part_num; i++)
	{
		particles.push_back(particle(i, mass));
	}

	auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
	mt19937 mt(seed);   // Seeds Mersenne Twister with Device RNG


	//Initialize positions
	{
		/*
		//Initialize particles on a cubic grid
		double thing = pow(part_num + 1, 1. / float(dim));

		int lattice_num = ceil(thing);
		double spacing = length / lattice_num;
		printf("lattice_num: %d spacing: %f\n", lattice_num, spacing);

		
		int part_counter = 0;
		for (int i = 0; i < lattice_num; i++)
		{
			for (int j = 0; j < lattice_num; j++)
			{
				for (int k = 0; k < lattice_num; k++)
				{
					if (part_counter < part_num)
					{
						particles[part_counter].pos[0] = (i + .5) * spacing;
						particles[part_counter].pos[1] = (j + .5) * spacing;
						particles[part_counter].pos[2] = (k + .5) * spacing;
						part_counter++;
					}
				}
			}
		}
		*/

		// Define probability density for perturbation
		double *density = new double[grid_size];
		//double arg_scale = 1/grid_size;
		double norm = 0;
		for (unsigned int i = 0; i < grid_size; i++)
		{
			density[i] = 1;
			norm += density[i];
		}

		//Normalize position probability density
		for (unsigned int i = 0; i < grid_size; i++)
		{
			density[i] /= norm;
		}

		// Sample indices
		vose * uniform_sampler = new vose(density, grid_size, mt);
		delete[] density;

		unsigned int j = 0;
		for (auto &part : particles)
		{
			for (int i = 0; i < 3; i++)
			{
				j = uniform_sampler->alias_method();
				//part.pos[i] += (float(j) / float(grid_size) - 1 / 2.) * spacing / 2;
				part.pos[i] = float(j)/float(grid_size)*length;
			}
		}

		delete uniform_sampler;
	}

	//Initialize momenta
	{
		// Define momentum density
		double *mom_density = new double[grid_size];
		double Temperature = 100;
		double mass = 4;
		double arg_scale = -1. / 1000.;

		double norm = 0;
		double tmp = 0;
		for(int i = 0; i < grid_size; i++)
		{
			tmp = i-grid_size/4;
			mom_density[i] = exp(arg_scale * tmp * tmp);
			tmp = i-grid_size*3/4;
			mom_density[i] += exp(arg_scale * tmp * tmp);
			norm += mom_density[i];
		}


		//Normalize momentum density
		for (unsigned int i = 0; i < grid_size; i++)
		{
			mom_density[i] /= norm;
		}

		// TODO Make this use STDEV

		vose * gaussian_sampler = new vose(mom_density, grid_size, mt);
		delete[] mom_density;

		double max_mom = 100;
		unsigned int j = 0;
		double new_mom = 0;
		for (auto &part : particles)
		{
			for (int i = 0; i < dim; i++)
			{
				j = gaussian_sampler->alias_method();         //select momentum index
				
				//convert j to momentum
				new_mom = (j - float(grid_size)/2) * 2*max_mom / float(grid_size);
				part.mom[i] = new_mom;
				
				/*
				//PacMan periodic boundary conditions
				if (fabs(new_mom) > max_mom)
				{
					part.mom[i] = fmod(fmod(part.mom[i], max_mom) +
					                   max_mom, max_mom);
				}
				else
				{
					part.mom[i] = new_mom;
				}
				*/
			}
		}

		delete gaussian_sampler;
	}
	return particles;
}

void update_positions(vector<particle> &particles, const double dt)
{
	double increment = 0;
	for (auto &part: particles)
	{
		increment = dt * part.mom[0] / part.mass;
		part.pos[0] += increment;
	}
}

template <typename T> int sign(T val)
{
	return (T(0) < val) - (val < T(0));
}

void update_state(vector<particle> &particles, const int Ng, const double dt,
                  const double length, double *charge, double *potential, double *EField)
{
	//size_t part_num = particles.size();

	//Calculate Fields
	double dx = length / Ng;

	update_positions(particles, dt);
	//Pac-Man shift particles (Periodic BCs)
	for (auto &part : particles)
	{
		for (int i = 0; i < dim; i++)
		{
			part.pos[i] = fmod(fmod(part.pos[i], length) + length, length);
		}
	}

	calc_charge(Ng, dx, charge, particles);
	calc_potential(Ng, dx, potential, charge);
	calc_EField(Ng, dx, EField, potential);

	for (auto &part: particles)
	{
		for (int i = 0; i < Ng; i++)
		{
			if (fabs(float(i) * dx - part.pos[0]) <= dx * 0.5)
			{
				//printf("EField i: %d\n", i);
				part.mom[0] += dt * 0.5 * EField[i];
				continue;         // To ensure only the closest EField is used.
			}
		}
	}
}

void write_state(vector<particle> &particles, double *charge, double *potential, double *EField, int Ng, int i)
{
	string tstring = to_string(i);
	string filename = "data/positions/pos" + tstring + ".csv";
	
	FILE *fp = fopen(filename.c_str(), "w");

	if (fp == nullptr)             //is it ok?
	{
		printf("Position output file did not open please check it\n");
	}

	for (auto &part : particles)
	{
		fprintf(fp, "%f, %f, %f, %f, %f, %f\n",
		        part.pos[0], part.pos[1], part.pos[2],
		        part.mom[0], part.mom[1], part.mom[2]);
	}
	fclose(fp);

	//--------------------------------------------------//

	filename = "data/field/field" + tstring + ".csv";
	fp = fopen(filename.c_str(), "w");

	if (fp == nullptr)             //is it ok?
	{
		printf("Field output file did not open please check it\n");
	}

	for (int i = 0; i < Ng; i++)
	{
		fprintf(fp, "%f\n", EField[i]);
	}

	fclose(fp);

	//--------------------------------------------------//

	filename = "data/charge/charge" + tstring + ".csv";
	fp = fopen(filename.c_str(), "w");

	if (fp == nullptr)             //is it ok?
	{
		printf("Charge output file did not open please check it\n");
	}

	for (int i = 0; i < Ng; i++)
	{
		fprintf(fp, "%f\n", charge[i]);
	}

	fclose(fp);

	//--------------------------------------------------//

	filename = "data/potential/pot" + tstring + ".csv";
	fp = fopen(filename.c_str(), "w");

	if (fp == nullptr)             //is it ok?
	{
		printf("Potential output file did not open please check it\n");
	}

	for (int i = 0; i < Ng; i++)
	{
		fprintf(fp, "%f\n", potential[i]);
	}

	fclose(fp);
}


double temperature(vector<particle>& particles)
{
	int N = static_cast<int>(particles.size());
	double sum = 0;

	for (auto &part : particles)
	{
		for (int k = 0; k < dim; k++)
		{
			sum += part.mom[k] * part.mom[k] / part.mass;
		}
	}
	return sum / N / dim;
}

void scale_mom(vector<particle> &particles, double scale)
{
	for (auto &part : particles)
	{
		for (int k = 0; k < dim; k++)
		{
			part.mom[k] *= scale;
		}
	}
}

void calc_potential(int Ng, double dx, double *potential, double *charge)
{
	//Allocate matrix
	double** matrix = new double*[Ng];
	for(int i = 0; i<Ng; i++)
	{
		matrix[i] = new double[Ng];
		for(int j = 0; j<Ng; j++)
			matrix[i][j] = 0;
	}

	//Set matrix
	for(int i = 0; i<Ng; i++)
	{
		matrix[i][i] = -2;
		if(i<Ng-1)
			matrix[i][i+1] = 1;
		if(i>0)
			matrix[i][i-1] = 1;
	}
	matrix[0][Ng-1] = 1;
	matrix[Ng-1][0] = 1;

	for (int i = 0; i < Ng; ++i)
		potential[i] = charge[i];

	//Gaussian Elimination
	double div = 1;
	double mult = 1;
	for(int i = 0; i<Ng; i++)
	{
		//Divide by diagonal element
		div = matrix[i][i];
		for(int j = 0; j<Ng; j++)
			matrix[i][j] /= div;
		potential[i] /= div;

		//Subtract from lower rows
		for (int j = i+1; j < Ng; ++j)
		{
			mult = matrix[j][i];
			for (int k = 0; k < Ng; ++k)
				matrix[j][k] -= mult*matrix[i][k];
			potential[j] -= mult*potential[i];
		}
	}

	//Backsubstitution
	for (int i = Ng-1; i>-1; --i)
	{
		for (int j = i+1; j < Ng; ++j)
			potential[i] -= matrix[i][j]*potential[j];
	}

	for (int i = 0; i < Ng; ++i)
	{
		delete[] matrix[i];
	}
	delete[] matrix;

	for (int i = 0; i < Ng; ++i)
		potential[i] *= dx*dx;
}

void calc_EField(int Ng, double dx, double *EField, double *potential)
{
	//Assuming periodic system
	EField[0] = (potential[1] - potential[Ng - 1]) / (2 * dx);
	EField[Ng - 1] = (potential[0] - potential[Ng - 2]) / (2 * dx);
	for (int i = 1; i < Ng - 1; i++)
	{
		EField[i] = (potential[i + 1] - potential[i - 1]) / (2*dx);
	}
}

void calc_charge(int Ng, double dx, double *charge, vector<particle> particles)
{
	double absx = 0;
	float Np = particles.size();
	//printf("Np: %f Ng: %d\n", Np, Ng);

	for (int i = 0; i < Ng; i++)
	{
		charge[i] = 0;
		for (auto &part: particles)
		{
			absx = fabs(part.pos[0] - i * dx);
			if (absx <= dx)
			{
				charge[i] -= float(Ng)/Np*(1 - absx / dx);
			}
		}
		charge[i] += 1;
	}
}
