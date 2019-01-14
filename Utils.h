#ifndef UTILS_H
#define UTILS_H

#include <gsl/gsl_rng.h>

gsl_rng* initialize_rng(int seed = 0);

#endif