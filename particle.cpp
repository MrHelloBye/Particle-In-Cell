//Particle.cpp
#include <array>
#include <iostream>
#include <cstdio>
#include <string>

#include "particle.h"

using namespace std;

particle::particle()
{
	particle_ID = 0;
	mass = 1;

    for(int i = 0; i<dim; i++)
    {
        pos[i] = 0;
        mom[i] = 0;
    }
}

particle::particle(unsigned int ID, double m)
{
    particle_ID = ID;
    mass = m;

    for(int i = 0; i<dim; i++)
    {
        pos[i] = 0;
        mom[i] = 0;
    }
}

particle::particle(unsigned int ID, double m, double* init_pos, double* init_mom)
{    
	particle_ID = ID;
	mass = m;
    
    for(int i = 0; i<dim; i++)
    {
    	pos[i] = init_pos[i];
        mom[i] = init_mom[i];
    }
}

void particle::setID(unsigned int ID)
{
	particle_ID = ID;
}

unsigned int particle::getID() const
{
	return particle_ID;
}

void particle::printPos() const
{
	for (int i = 0; i < sizeof(pos)/sizeof(pos[0]); i++)
	{
		std::cout << this->pos[i] << '\n';
	}
}
