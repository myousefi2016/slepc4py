#ifndef SLEPC4PY_CUSTOM_H
#define SLEPC4PY_CUSTOM_H

#undef  __FUNCT__
#define __FUNCT__ "SlepcInitializePackage"
static PetscErrorCode SlepcInitializePackage(const char path[])
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
#if SLEPC_VERSION_(3,2,0)
  ierr = PetscClassIdRegister("Direct solver",&DS_CLASSID);CHKERRQ(ierr);
#else
  ierr = 0; CHKERRQ(ierr);
#endif
  PetscFunctionReturn(0);
}

#endif/*SLEPC4PY_CUSTOM_H*/
