#ifndef SLEPC4PY_CUSTOM_H
#define SLEPC4PY_CUSTOM_H

#undef  __FUNCT__
#define __FUNCT__ "SlepcInitializePackageAll"
static PetscErrorCode SlepcInitializePackageAll(void)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = EPSInitializePackage();CHKERRQ(ierr);
  ierr = SVDInitializePackage();CHKERRQ(ierr);
  ierr = PEPInitializePackage();CHKERRQ(ierr);
  ierr = NEPInitializePackage();CHKERRQ(ierr);
  ierr = MFNInitializePackage();CHKERRQ(ierr);
  ierr = STInitializePackage();CHKERRQ(ierr);
  ierr = BVInitializePackage();CHKERRQ(ierr);
  ierr = DSInitializePackage();CHKERRQ(ierr);
  ierr = FNInitializePackage();CHKERRQ(ierr);
  ierr = RGInitializePackage();CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#endif/*SLEPC4PY_CUSTOM_H*/
