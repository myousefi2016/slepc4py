#ifndef SLEPC4PY_CUSTOM_H
#define SLEPC4PY_CUSTOM_H

#undef  __FUNCT__
#define __FUNCT__ "SlepcInitializePackage"
static PetscErrorCode SlepcInitializePackage(const char path[])
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = 0; CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#endif/*SLEPC4PY_CUSTOM_H*/
