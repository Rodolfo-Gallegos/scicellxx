#ifndef SCICELLXX_LINEARSOLVERS_H
#define SCICELLXX_LINEARSOLVERS_H

// Include the linear solver
#include "ac_linear_solver.h"
#include "cc_lu_solver_numerical_recipes.h"
#ifdef SCICELLXX_USES_ARMADILLO
#include "cc_solver_armadillo.h"
#endif // #ifdef SCICELLXX_USES_ARMADILO

#include "cc_factory_linear_solver.h"

#endif // #ifndef SCICELLXX_LINEARSOLVERS_H
