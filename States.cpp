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
    const double recovery = (params->recovery),
                 beta     = (params->beta),
                 birth    = (params->birth),
                 death    = (params->death),
                 kappa    = (params->kappa),
                 rho      = (params->rho),
                 Population = 1.0; // (params->Population[0]); // assume const population

    const double S  = gsl_vector_get(x, S_STATE),
                 I1 = gsl_vector_get(x, I1_STATE),
                 R  = gsl_vector_get(x, R_STATE),
                 P  = gsl_vector_get(x, P_STATE),
                 Ir = gsl_vector_get(x, IR_STATE);
    // assert((S+I1+R+P+Ir) <= 1.0);
    // assert each >= 0

    const double foi = ((beta*(I1+kappa*Ir))/Population);
    const double birthdt = birth*Population,
                 recovery1 = recovery*I1,
                 recovery2 = (recovery/kappa)*Ir,
                 waning = rho*R,
                 infection1 = S*foi,
                 infection2 = P*foi*kappa;

    gsl_vector_set(f, S_STATE,  birthdt    - infection1 - death*S);
    gsl_vector_set(f, I1_STATE, infection1 - recovery1  - death*I1);
    gsl_vector_set(f, R_STATE,  recovery1  + recovery2  - death*R - waning);
    gsl_vector_set(f, P_STATE,  waning     - infection2 - death*P);
    gsl_vector_set(f, IR_STATE, infection2 - recovery2  - death*Ir);

    return GSL_SUCCESS;
}

std::vector<double> multiroot_solver(
  Params p, int MAX_ITER = 100, bool verbose = false
) {

  gsl_multiroot_function F = { &multiroot_evolve, NUM_OF_STATE_TYPES, &p };
  gsl_multiroot_fsolver * solver = gsl_multiroot_fsolver_alloc(gsl_multiroot_fsolver_hybrids, NUM_OF_STATE_TYPES);
  gsl_vector * x = gsl_vector_alloc(NUM_OF_STATE_TYPES); // allocate a zeroes vector
  gsl_vector_set_all(x, 1.0/((double)NUM_OF_STATE_TYPES)); // distribute evenly
  gsl_multiroot_fsolver_set(solver, &F, x);

  if (verbose) fprintf(stderr, "F solver: %s\n", gsl_multiroot_fsolver_name(solver) );

  int iter = 0;
  do {
    gsl_multiroot_fsolver_iterate(solver);
  } while (
    (gsl_multiroot_test_residual(solver->f, 1e-7) == GSL_CONTINUE) && (++iter < MAX_ITER)
  );

  auto compartments = from_gsl(solver->x);

  gsl_multiroot_fsolver_free(solver);
  gsl_vector_free(x);

  return compartments;
}

// UNIROOT CODE

double func_s(double I1, Params* p) {
  return 1.0 - ((p->recovery + p->death)/(p->death))*I1;
}

double func_ir(double I1, double S, Params* p) {
  return (((p->recovery+p->death)/(p->beta * S)) - 1.0)*(I1 / p->kappa);
}

double func_r(double I1, double IR, Params* p) {
  return (p->recovery)/(p->rho + p->death)*(I1 + IR/p->kappa);
}

double func_p(double I1, double IR, Params* p) {
  return ((p->recovery)/(p->kappa) + p->death) * IR/(p->kappa * p->beta * (I1 + p->kappa * IR));
}

double uniroot_evolve(double I1, void * p){
    Params * ps = (Params *)p;

    double S  = func_s(I1, ps);
    double IR = func_ir(I1, S, ps);
    double R  = func_r(I1, IR, ps),
           P  = func_p(I1, IR, ps);

    return 1 - S - I1 - IR - R - P;
}

std::vector<double> from_I1(double I1, Params * p) {
  double S = func_s(I1, p);
  double IR = func_ir(I1, S, p);
  return {
    S, I1, func_r(I1, IR, p), func_p(I1, IR, p), IR
  };
}

std::vector<double> uniroot_solver(
  Params params, int MAX_ITER = 100, bool verbose = false
) {

  gsl_function F = { &uniroot_evolve, &params };

  auto solver = gsl_root_fsolver_alloc(gsl_root_fsolver_brent);
  if (verbose) printf("F solver: %s\n", gsl_root_fsolver_name(solver));

  double i1lo = params.death*(1.0/(params.death+params.recovery) - 1.0/params.beta),
         i1hi = (params.death / (params.recovery + params.death)) - 1e-8;

  gsl_root_fsolver_set(solver, &F, i1lo, i1hi);

  int iter = 0;

  /* main loop */
  do {
    gsl_root_fsolver_iterate(solver);
    i1lo = gsl_root_fsolver_x_lower(solver);
    i1hi = gsl_root_fsolver_x_upper(solver);
  } while (
    (gsl_root_test_interval(i1lo, i1hi, 0, 0.001) == GSL_CONTINUE) && (++iter < MAX_ITER)
  );

  auto res = from_I1(gsl_root_fsolver_root(solver), &params);

  gsl_root_fsolver_free(solver);

  return res;

}

// final definition

std::vector<double> equilibrium_fraction(Params p, bool multi) {
  return multi ? multiroot_solver(p) : uniroot_solver(p);
}
