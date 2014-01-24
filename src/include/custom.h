#ifndef SLEPC4PY_CUSTOM_H
#define SLEPC4PY_CUSTOM_H

#if SLEPC_VERSION_LE(3,3,0)
#define EPSInitializePackage() EPSInitializePackage(0)
#define SVDInitializePackage() SVDInitializePackage(0)
#define QEPInitializePackage() QEPInitializePackage(0)
#define NEPInitializePackage() NEPInitializePackage(0)
#define MFNInitializePackage() MFNInitializePackage(0)
#define STInitializePackage()  STInitializePackage(0)
#define IPInitializePackage()  IPInitializePackage(0)
#define DSInitializePackage()  DSInitializePackage(0)
#define FNInitializePackage()  FNInitializePackage(0)
#endif

#undef  __FUNCT__
#define __FUNCT__ "SlepcInitializePackageAll"
static PetscErrorCode SlepcInitializePackageAll(void)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = EPSInitializePackage();CHKERRQ(ierr);
  ierr = SVDInitializePackage();CHKERRQ(ierr);
  ierr = QEPInitializePackage();CHKERRQ(ierr);
  ierr = NEPInitializePackage();CHKERRQ(ierr);
  ierr = MFNInitializePackage();CHKERRQ(ierr);
  ierr = STInitializePackage();CHKERRQ(ierr);
  ierr = IPInitializePackage();CHKERRQ(ierr);
  ierr = DSInitializePackage();CHKERRQ(ierr);
  /*ierr = FNInitializePackage();CHKERRQ(ierr);*/
  PetscFunctionReturn(0);
}

#endif/*SLEPC4PY_CUSTOM_H*/
