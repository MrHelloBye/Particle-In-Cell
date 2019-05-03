//Vose Random Number Generation
#include "vose.h"

#include <random>
#include <iostream>
#include <limits>
#include <cmath>
#include <chrono>

using namespace std;

vose::vose(double *array, unsigned int size, mt19937 &mt)
	: local_mt(mt), array_size(size)
{
	//Initialize
	//std::random_device rd; // RNG that uses device for seed
	//auto seed = chrono::high_resolution_clock::now().time_since_epoch().count();
	//mt19937 mt(seed); // Seeds Mersenne Twister with Device RNG

	//mt19937 &local_mt = mt;

	// nextafter(x, y) Returns the next representable
	// value after x in the direction of y.
	uniform_real_distribution<double> uniform(0, 1);

	distribution = new double[array_size];
	distribution = (double *) memcpy(distribution, array, sizeof(double) * array_size);

	//Vose initialization
	unsigned l = 0; unsigned s = 0;
	int *large = new int[size]; int *small = new int[size];
	alias = new unsigned [size];
	prob = new unsigned[size];
	for(int j = 0; j<size; j++)
	{
		large[j] = 0;
		small[j] = 0;
		alias[j] = 0;
		prob[j] = 0;
	}

	for(int j = 0; j<size; j++)
	{
		if (distribution[j]>1/float(size))
		{
			large[l] = j;
			++l;
		}
		else
		{
			small[s] = j;
			++s;
		}
	}

	unsigned j = 0; unsigned k = 0;
	while(s!=0 && l!=0)
	{
		--s;
		j = small[s];
		--l;
		k = large[l];
		
		prob[j] = size*distribution[j];
		alias[j] = k;
		
		distribution[k] += distribution[j] - 1/float(size);
		if (distribution[k] > 1/float(size))
		{
			large[l] = k;
			++l;
		}
		else
		{
			small[s] = k;
			++s;
		}
	}

	/*
	for (int i = 0; i < size; ++i)
	{
		printf("small[%d]: %d\n", i, small[i]);
	}
	for (int i = 0; i < size; ++i)
	{
		printf("large[%d]: %d\n", i, large[i]);
	}
	for (int i = 0; i < size; ++i)
	{
		printf("alias[%d]: %d\n", i, alias[i]);
	}
	cin.get();
	*/

	unsigned temp = 0;
	while(s>0)
	{
		--s;
		temp = small[s];
		prob[temp] = 1;
	}


	while(l>0)
	{
		--l;
		prob[large[l]] = 1;
	}
}

vose::~vose()
{
	cout << "Vose Sampler Destructor." << endl;
	delete[] distribution;
	delete[] alias;
}


void vose::demo()
{
	for (int i = 0; i < 10; ++i)
	{
		cout << uniform(local_mt) << endl;
	}
}

unsigned int vose::alias_method()
{
	double u = array_size*uniform(local_mt);
	unsigned j = floor(u);

	if (j >= array_size)
	{
		cout << "Index out of bounds." << endl;
	}

	if (u-j <= prob[j])
	{
		return j;
	}
	else
	{
		return alias[j];
	}
}
