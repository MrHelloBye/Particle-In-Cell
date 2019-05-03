#ifndef VOSE_H
#define VOSE_H

#include <random>
#include <cstring>

class vose
{
public:
    vose(double*, unsigned int, std::mt19937&);
    ~vose();
    void demo();
    unsigned int alias_method();

private:
    std::mt19937 &local_mt;

    // Don't declare the distribution with any arguments!!
    std::uniform_real_distribution<double> uniform;

    unsigned int array_size = 0;

    double *distribution = NULL;
    unsigned int *alias = NULL;
    unsigned int *prob = NULL;
};

#endif /* VOSE_H */
