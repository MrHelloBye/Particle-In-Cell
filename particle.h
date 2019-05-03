#ifndef particle_H
#define particle_H

#include <array>

constexpr unsigned int dim = 3;

class particle
{
private:
	unsigned int particle_ID;

public:
	particle();
    particle(unsigned int ID, double m);
	particle(unsigned int ID, double m, double* pos, double* mom);
	
	double mass;
    std::array<double, dim> pos;
    std::array<double, dim> mom;
    
	void setID(unsigned int n);
	unsigned int getID() const;

	void printPos() const;
};

#endif
