#if !defined(SLEPC_COMPAT_H)
#define SLEPC_COMPAT_H


#if !defined(SLEPC_VERSION_)
#define SLEPC_VERSION_(MAJOR,MINOR,SUBMINOR) \
  ((SLEPC_VERSION_MAJOR == (MAJOR)) &&       \
   (SLEPC_VERSION_MINOR == (MINOR)) &&       \
   (SLEPC_VERSION_SUBMINOR == (SUBMINOR)) && \
   (SLEPC_VERSION_RELEASE  == 1))
#endif

#if SLEPC_VERSION_(3,0,0)
#undef __FUNCT__
#define __FUNCT__ "IPOrthogonalize_300"
PETSC_STATIC_INLINE
PetscErrorCode IPOrthogonalize_300(IP ip,
				   PetscInt nds,Vec *DS,
				   PetscInt n,PetscTruth *which,Vec *V,
				   Vec v,PetscScalar *H,PetscReal *norm,
				   PetscTruth *lindep)
{
  PetscErrorCode ierr;
  Vec *vwork;
  PetscScalar *swork;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ip,IP_COOKIE,1);
  if (nds > 0) {
    SETERRQ(PETSC_ERR_SUP,"operation not supported in this SLEPc version");
  }
  ierr = IPOrthogonalize(ip,n,wich,V,v,H,norm,lindep,
			 PETSC_NULL,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#define IPOrthogonalize IPOrthogonalize_300
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
#define __FUNCT__ "EPSGetOperators_233"
PETSC_STATIC_INLINE
PetscErrorCode EPSGetOperators_233(EPS eps, Mat *A, Mat *B)
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
#define EPSGetOperators EPSGetOperators_233

#undef __FUNCT__
#define __FUNCT__ "EPSSetOperators_233"
PETSC_STATIC_INLINE
PetscErrorCode EPSSetOperators_233(EPS eps, Mat A, Mat B)
{
  ST st;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_COOKIE,1);
  if (A) PetscValidPointer(A,2);
  if (B) PetscValidPointer(B,3);
  ierr = EPSSetOperators(eps,A,B);CHKERRQ(ierr);
  ierr = EPSGetST(eps,&st);CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject)st,"__SLEPc_ST_op_A__",(PetscObject)A);CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject)st,"__SLEPc_ST_op_B__",(PetscObject)B);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#define EPSSetOperators EPSSetOperators_233

#undef __FUNCT__
#define __FUNCT__ "STSetOperators_233"
PETSC_STATIC_INLINE
PetscErrorCode STSetOperators_233(ST st, Mat A, Mat B)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_COOKIE,1);
  if (A) PetscValidPointer(A,2);
  if (B) PetscValidPointer(B,3);
  ierr = STSetOperators(st,A,B);CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject)st,"__SLEPc_ST_op_A__",(PetscObject)A);CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject)st,"__SLEPc_ST_op_B__",(PetscObject)B);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#define STSetOperators STSetOperators_233

#endif


#if SLEPC_VERSION_(2,3,3)

typedef enum {
  EPS_RITZ=1,
  EPS_HARMONIC,
  EPS_REFINED,
  EPS_REFINED_HARMONIC
} EPSExtraction;

#undef __FUNCT__
#define __FUNCT__ "EPSSetExtraction"
PETSC_STATIC_INLINE
PetscErrorCode EPSSetExtraction(EPS eps,EPSExtraction ext)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_COOKIE,1);
  SETERRQ(PETSC_ERR_SUP,"operation not supported in this SLEPc version");
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "EPSGetExtraction"
PETSC_STATIC_INLINE
PetscErrorCode EPSGetExtraction(EPS eps,EPSExtraction *ext)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_COOKIE,1);
  SETERRQ(PETSC_ERR_SUP,"operation not supported in this SLEPc version");
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "EPSSetTarget"
PETSC_STATIC_INLINE
PetscErrorCode EPSSetTarget(EPS eps,PetscScalar target)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_COOKIE,1);
  SETERRQ(PETSC_ERR_SUP,"operation not supported in this SLEPc version");
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "EPSGetTarget"
PETSC_STATIC_INLINE
PetscErrorCode EPSGetTarget(EPS eps,PetscScalar *target)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_COOKIE,1);
  SETERRQ(PETSC_ERR_SUP,"operation not supported in this SLEPc version");
  PetscFunctionReturn(0);
}

#endif


#if SLEPC_VERSION_(2,3,3)

#undef __FUNCT__
#define __FUNCT__ "EPSSetDimensions_233"
PETSC_STATIC_INLINE
PetscErrorCode EPSSetDimensions_233(EPS eps,PetscInt nev,PetscInt ncv,PetscInt mpd)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_COOKIE,1);
  ierr = EPSSetDimensions(eps,nev,ncv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#define EPSSetDimensions EPSSetDimensions_233

#undef __FUNCT__
#define __FUNCT__ "EPSGetDimensions_233"
PETSC_STATIC_INLINE
PetscErrorCode EPSGetDimensions_233(EPS eps, PetscInt *nev,PetscInt *ncv,PetscInt *mpd)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_COOKIE,1);
  ierr = EPSGetDimensions(eps,nev,ncv);CHKERRQ(ierr);
  if (mpd) mpd = 0;
  PetscFunctionReturn(0);
}
#define EPSGetDimensions EPSGetDimensions_233

#undef __FUNCT__
#define __FUNCT__ "SVDSetDimensions_233"
PETSC_STATIC_INLINE
PetscErrorCode SVDSetDimensions_233(SVD svd,PetscInt nev,PetscInt ncv,PetscInt mpd)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_COOKIE,1);
  ierr = SVDSetDimensions(svd,nev,ncv);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#define SVDSetDimensions SVDSetDimensions_233

#undef __FUNCT__
#define __FUNCT__ "SVDGetDimensions_233"
PETSC_STATIC_INLINE
PetscErrorCode SVDGetDimensions_233(SVD svd, PetscInt *nev,PetscInt *ncv,PetscInt *mpd)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_COOKIE,1);
  ierr = SVDGetDimensions(svd,nev,ncv);CHKERRQ(ierr);
  if (mpd) mpd = 0;
  PetscFunctionReturn(0);
}
#define SVDGetDimensions SVDGetDimensions_233

#endif


#if SLEPC_VERSION_(2,3,3)
#undef __FUNCT__
#define __FUNCT__ "IPOrthogonalize_233"
PETSC_STATIC_INLINE
PetscErrorCode IPOrthogonalize_233(IP ip,PetscInt n,PetscTruth *which,Vec *V,Vec v,
				   PetscScalar *H,PetscReal *norm,PetscTruth *lindep,
				   Vec work,PetscScalar* swork)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = IPOrthogonalize(ip,n,which,V,v,H,norm,lindep,work);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#define IPOrthogonalize IPOrthogonalize_233

#endif

#endif /* !SLEPC_COMPAT_H */
