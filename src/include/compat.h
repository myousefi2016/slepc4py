#ifndef SLEPC4PY_COMPAT_H
#define SLEPC4PY_COMPAT_H

#include <slepc.h>

#if PETSC_VERSION_(3,2,0)
#define PetscShell PetscFwk
#endif

#if SLEPC_VERSION_(3,2,0)
#define IPBILINEAR     "bilinear"
#define IPSESQUILINEAR "sesquilinear"
#define IPINDEFINITE   "indefinite"
#define EPSRQCG "rqcg"
#define EPS_GHIEP ((EPSProblemType)-1)
#endif

#endif/*SLEPC4PY_COMPAT_H*/
