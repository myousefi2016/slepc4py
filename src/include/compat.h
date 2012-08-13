#ifndef SLEPC4PY_COMPAT_H
#define SLEPC4PY_COMPAT_H

#include <slepc.h>

#if PETSC_VERSION_(3,2,0)
#define PetscShell PetscFwk
#endif

#if SLEPC_VERSION_(3,2,0)
#include "compat/slepc-32.h"
#endif


#endif/*SLEPC4PY_COMPAT_H*/
