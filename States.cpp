#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_math.h>

#include "States.h"

void printResults(std::vector<double> res) {
  for(int i = 0; i < NUM_OF_STATE_TYPES; i++) {
    fprintf(stderr, "%2s %25.17e\n",
      statestr[i].c_str(),
      res[i]
    );
  }
}

std::vector<double> from_gsl(gsl_vector * gv) {
  return {
    gsl_vector_get(gv, S_STATE),
    gsl_vector_get(gv, I1_STATE),
    gsl_vector_get(gv, R_STATE),
    gsl_vector_get(gv, P_STATE),
    gsl_vector_get(gv, IR_STATE)
  };
}

// MULTIROOT CODE

int multiroot_evolve(const gsl_vector * x, void * p, gsl_vector * f){
    Params * params = (Params *)p;
    const double recovery = (params->recovery);
    const double beta     = (params->beta);
    const double birth    = (params->birth);
    const double death    = (params->death);
    const double kappa    = (params->kappa);
    const double rho      = (params->rho);
    const double Population = 1.0; // (params->Population[0]); // assume const population
    const double S  = gsl_vector_get(x, S_STATE);
    const double I1 = gsl_vector_get(x, I1_STATE);
    const double R  = gsl_vector_get(x, R_STATE);
    const double P  = gsl_vector_get(x, P_STATE);
    const double Ir = gsl_vector_get(x, IR_STATE);
    // assert((S+I1+R+P+Ir) <= 1.0);
    // assert each >= 0

    const double birthdt = birth*Population;
    const double foi = ((beta*(I1+kappa*Ir))/Population);
    const double infection1 = S*foi;
    const double infection2 = P*foi*kappa;
    const double recovery1 = recovery*I1;
    const double recovery2 = (recovery/kappa)*Ir;
    const double waning = rho*R;

    gsl_vector_set(f, S_STATE,  birthdt - infection1 - death*S);
    gsl_vector_set(f, I1_STATE, infection1 - recovery1 - death*I1);
    gsl_vector_set(f, R_STATE,  recovery1 + recovery2 - waning - death*R);
    gsl_vector_set(f, P_STATE,  waning - infection2 - death*P);
    gsl_vector_set(f, IR_STATE, infection2 - recovery2 - death*Ir);

    return GSL_SUCCESS;
}

std::vector<double> multiroot_solver(Params p, int MAXTIMES = 100, bool verbose = false) {

  gsl_multiroot_function F = { &multiroot_evolve, NUM_OF_STATE_TYPES, &p };

  gsl_multiroot_fsolver * solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids, NUM_OF_STATE_TYPES);

  gsl_vector * x = gsl_vector_alloc(NUM_OF_STATE_TYPES); // allocate a zeroes vector
  gsl_vector_set_all(x, 1.0/((double)NUM_OF_STATE_TYPES));

  /* set solver */
  gsl_multiroot_fsolver_set(solver, &F, x);

  int times = 0, status;

  if (verbose) fprintf(stderr, "F solver: %s\n", gsl_multiroot_fsolver_name(solver) );

  do {
    times++;
    status = gsl_multiroot_fsolver_iterate(solver);
    status = gsl_multiroot_test_residual (solver->f, 1e-7);
  } while (status == GSL_CONTINUE && times < MAXTIMES);

  std::vector<double> compartments = from_gsl(solver->x);

  gsl_multiroot_fsolver_free(solver);
  gsl_vector_free(x);

  return compartments;
}

// UNIROOT CODE

std::vector<double> uniroot_solver(Params p) {
  return std::vector<double>();
}

// final definition

std::vector<double> equilibrium_fraction(Params p) { return multiroot_solver(p); }
