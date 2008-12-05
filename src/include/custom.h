#if !defined(SLEPC_VERSION_)
#define SLEPC_VERSION_(MAJOR,MINOR,SUBMINOR) \
  ((SLEPC_VERSION_MAJOR == (MAJOR)) &&       \
   (SLEPC_VERSION_MINOR == (MINOR)) &&       \
   (SLEPC_VERSION_SUBMINOR == (SUBMINOR)) && \
   (SLEPC_VERSION_RELEASE  == 1))
#endif


#if SLEPC_VERSION_(2,3,3)
#undef __FUNCT__  
#define __FUNCT__ "IPGetOptionsPrefix"
PETSC_STATIC_INLINE 
PetscErrorCode IPGetOptionsPrefix(IP ip,const char *prefix[])
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ip,IP_COOKIE,1);
  PetscValidPointer(prefix,2);
  ierr = PetscObjectGetOptionsPrefix((PetscObject)ip, prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#endif

#if SLEPC_VERSION_(2,3,3)
#undef __FUNCT__  
#define __FUNCT__ "EPSGetOperators"
PETSC_STATIC_INLINE 
PetscErrorCode EPSGetOperators(EPS eps, Mat *A, Mat *B)
{
  ST st;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_COOKIE,1);
  if (A) PetscValidPointer(A,2);
  if (B) PetscValidPointer(B,3);
  ierr = EPSGetST(eps,&st);CHKERRQ(ierr);
  ierr = STGetOperators(st,A,B);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#endif
