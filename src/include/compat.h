#ifndef SLEPC4PY_COMPAT_H
#define SLEPC4PY_COMPAT_H

#include <slepc.h>

#if !defined(SLEPC_VERSION_)
#define SLEPC_VERSION_(MAJOR,MINOR,SUBMINOR) \
  ((SLEPC_VERSION_MAJOR == (MAJOR))       && \
   (SLEPC_VERSION_MINOR == (MINOR))       && \
   (SLEPC_VERSION_SUBMINOR == (SUBMINOR)) && \
   (SLEPC_VERSION_RELEASE  == 1))
#endif

#if SLEPC_VERSION_(3,1,0) || SLEPC_VERSION_(3,0,0) || SLEPC_VERSION_(2,3,3)

#if SLEPC_VERSION_(3,1,0)
#include <slepcqep.h>
#elif SLEPC_VERSION_(3,0,0) || SLEPC_VERSION_(2,3,3)
#include <petscts.h>
#include <petscda.h>
#include <slepcsvd.h>
#include "compat/slepcqep.h"
#endif

#define PetscBool PetscTruth

#define PetscClassId PetscCookie
#define ST_CLASSID  ST_COOKIE
#define IP_CLASSID  IP_COOKIE
#define EPS_CLASSID EPS_COOKIE
#define SVD_CLASSID SVD_COOKIE
#define QEP_CLASSID QEP_COOKIE

#undef  __FUNCT__
#define __FUNCT__ "SlepcInitializePackage"
static PetscErrorCode SlepcInitializePackage(const char path[])
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  ierr = 0;
#if SLEPC_VERSION_(3,0,0)
  ierr = PetscCookieRegister("Quadratic Eigenproblem Solver",&QEP_COOKIE);CHKERRQ(ierr);
#elif SLEPC_VERSION_(2,3,3)
  ierr = PetscLogClassRegister(&QEP_COOKIE,"Quadratic Eigenproblem Solver");CHKERRQ(ierr);
#endif
  PetscFunctionReturn(0);
}

#if SLEPC_VERSION_(3,0,0) || SLEPC_VERSION_(2,3,3)
/**/
#define EPSGD  "gd"
#define EPSJD  "jd"
/**/
#define EPS_HARMONIC_RELATIVE EPS_HARMONIC
#define EPS_HARMONIC_RIGHT    EPS_HARMONIC
#define EPS_HARMONIC_LARGEST  EPS_HARMONIC
/**/
#endif

#if SLEPC_VERSION_(3,0,0) || SLEPC_VERSION_(2,3,3)
/**/
#define STSINVERT  STSINV
#define STPRECOND  "precond"
/**/
#define ST_MATMODE_COPY    STMATMODE_COPY
#define ST_MATMODE_INPLACE STMATMODE_INPLACE
#define ST_MATMODE_SHELL   STMATMODE_SHELL
/**/
#endif

#if SLEPC_VERSION_(3,0,0) || SLEPC_VERSION_(2,3,3)
/**/
#define IP_ORTH_MGS IP_MGS_ORTH
#define IP_ORTH_CGS IP_CGS_ORTH
/**/
#define IP_INNER_HERMITIAN IPINNER_HERMITIAN
#define IP_INNER_SYMMETRIC IPINNER_SYMMETRIC
/**/
#undef __FUNCT__
#define __FUNCT__ "IPOrthogonalize_300"
static PetscErrorCode IPOrthogonalize_300(IP ip,
                                          PetscInt nds,Vec *DS,
                                          PetscInt n,PetscTruth *which,Vec *V,
                                          Vec v,PetscScalar *H,PetscReal *norm,
                                          PetscTruth *lindep)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ip,IP_COOKIE,1);
  if (nds > 0) {
    SETERRQ(PETSC_ERR_SUP,"operation not supported in this SLEPc version");
  }
  #if SLEPC_VERSION_(3,0,0)
  ierr = IPOrthogonalize(ip,n,which,V,v,H,norm,lindep,
			 PETSC_NULL,0);CHKERRQ(ierr);
  #elif SLEPC_VERSION_(2,3,3)
  ierr = IPOrthogonalize(ip,n,which,V,v,H,norm,lindep,
                         /*work*/0);CHKERRQ(ierr);
  #endif
  PetscFunctionReturn(0);
}
#define IPOrthogonalize IPOrthogonalize_300
/**/
#endif


#if SLEPC_VERSION_(3,0,0) || SLEPC_VERSION_(2,3,3)
/**/
#define EPSDSITRLANCZOS "dsitrlanczos"
/**/
#define EPS_GHIEP ((EPSProblemType)0)
/**/
#define EPS_LANCZOS_REORTHOG_LOCAL     EPSLANCZOS_REORTHOG_LOCAL
#define EPS_LANCZOS_REORTHOG_FULL      EPSLANCZOS_REORTHOG_FULL
#define EPS_LANCZOS_REORTHOG_SELECTIVE EPSLANCZOS_REORTHOG_SELECTIVE
#define EPS_LANCZOS_REORTHOG_PERIODIC  EPSLANCZOS_REORTHOG_PERIODIC
#define EPS_LANCZOS_REORTHOG_PARTIAL   EPSLANCZOS_REORTHOG_PARTIAL
#define EPS_LANCZOS_REORTHOG_DELAYED   EPSLANCZOS_REORTHOG_DELAYED
/**/
#define EPS_POWER_SHIFT_CONSTANT  EPSPOWER_SHIFT_CONSTANT
#define EPS_POWER_SHIFT_RAYLEIGH  EPSPOWER_SHIFT_RAYLEIGH
#define EPS_POWER_SHIFT_WILKINSON EPSPOWER_SHIFT_WILKINSON
/**/
#define EPS_TARGET_MAGNITUDE ((EPSWhich)0)
#define EPS_TARGET_REAL      ((EPSWhich)0)
#define EPS_TARGET_IMAGINARY ((EPSWhich)0)
#define EPS_WHICH_USER       ((EPSWhich)0)
/**/
#undef __FUNCT__  
#define __FUNCT__ "EPSGetLeftVectorsWanted"
static PetscErrorCode EPSGetLeftVectorsWanted(EPS eps,PetscTruth *leftvecs)
{
  EPSClass       epsclass;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_COOKIE,1);
  PetscValidPointer(leftvecs,2);
  ierr = EPSGetClass(eps,&epsclass);CHKERRQ(ierr);
  *leftvecs = (epsclass==EPS_TWO_SIDE)?PETSC_TRUE:PETSC_FALSE;
  PetscFunctionReturn(0);
}
#undef __FUNCT__  
#define __FUNCT__ "EPSSetLeftVectorsWanted"
static PetscErrorCode EPSSetLeftVectorsWanted(EPS eps,PetscTruth leftvecs)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_COOKIE,1);
  ierr = EPSSetClass(eps,leftvecs?EPS_TWO_SIDE:EPS_ONE_SIDE);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
/**/
typedef enum { 
  EPS_BALANCE_NONE=0,
  EPS_BALANCE_ONESIDE,
  EPS_BALANCE_TWOSIDE,
  EPS_BALANCE_USER 
} EPSBalance;
#undef __FUNCT__
#define __FUNCT__ "EPSSetBalance"
static PetscErrorCode EPSSetBalance(EPS eps,EPSBalance b,PetscInt i,PetscReal c)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_COOKIE,1);
  SETERRQ(PETSC_ERR_SUP,"operation not supported in this SLEPc version");
  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "EPSGetBalance"
static PetscErrorCode EPSGetBalance(EPS eps,EPSBalance *b,PetscInt *i,PetscReal *c)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_COOKIE,1);
  PetscValidPointer(b,2);
  PetscValidPointer(i,3);
  PetscValidPointer(c,4);
  SETERRQ(PETSC_ERR_SUP,"operation not supported in this SLEPc version");
  PetscFunctionReturn(0);
}
/**/
#define EPSSetDeflationSpace(eps,n,ds) \
        EPSAttachDeflationSpace(eps,n,ds,PETSC_FALSE)
/**/
#define EPSGetEigenvalue      EPSGetValue
#define EPSGetEigenvector     EPSGetRightVector
#define EPSGetEigenvectorLeft EPSGetLeftVector
/**/
#undef __FUNCT__
#define __FUNCT__ "EPSSetInitialSpace"
static PetscErrorCode EPSSetInitialSpace(EPS eps, PetscInt n, Vec *is)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_COOKIE,1);
  if (n != 1) {
    SETERRQ(PETSC_ERR_SUP,"operation not supported in this SLEPc version");
  }
  PetscValidPointer(is,3);
  PetscValidHeaderSpecific(is[0],VEC_COOKIE,3);
  ierr = EPSSetInitialVector(eps,is[0]);CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject)eps,"__SLEPc_EPS_IniVec__",
                            (PetscObject)is[0]);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "EPSSetInitialSpaceLeft"
static PetscErrorCode EPSSetInitialSpaceLeft(EPS eps, PetscInt n, Vec *is)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_COOKIE,1);
  if (n != 1) {
    SETERRQ(PETSC_ERR_SUP,"operation not supported in this SLEPc version");
  }
  PetscValidPointer(is,3);
  PetscValidHeaderSpecific(is[0],VEC_COOKIE,3);
  ierr = EPSSetLeftInitialVector(eps,is[0]);CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject)eps,"__SLEPc_EPS_LeftIniVec__",
                            (PetscObject)is[0]);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "EPSSetTrackAll"
static PetscErrorCode EPSSetTrackAll(EPS eps,PetscTruth trackall)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_COOKIE,1);
  SETERRQ(PETSC_ERR_SUP,"operation not supported in this SLEPc version");
  PetscFunctionReturn(PETSC_ERR_SUP);
}
#undef __FUNCT__
#define __FUNCT__ "EPSGetTrackAll"
static PetscErrorCode EPSGetTrackAll(EPS eps,PetscTruth *trackall) 
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_COOKIE,1);
  SETERRQ(PETSC_ERR_SUP,"operation not supported in this SLEPc version");
  PetscFunctionReturn(PETSC_ERR_SUP);
}
#endif

#if SLEPC_VERSION_(3,0,0) || SLEPC_VERSION_(2,3,3)
/**/
#undef __FUNCT__
#define __FUNCT__ "SVDSetInitialSpace"
static PetscErrorCode SVDSetInitialSpace(SVD svd, PetscInt n, Vec *is)
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(svd,SVD_COOKIE,1);
  if (n != 1) {
    SETERRQ(PETSC_ERR_SUP,"operation not supported in this SLEPc version");
  }
  PetscValidPointer(is,3);
  PetscValidHeaderSpecific(is[0],VEC_COOKIE,3);
  ierr = SVDSetInitialVector(svd,is[0]);CHKERRQ(ierr);
  ierr = PetscObjectCompose((PetscObject)svd,"__SLEPc_SVD_IniVec__",
                            (PetscObject)is[0]);CHKERRQ(ierr);
  PetscFunctionReturn(0);

}
/**/
#endif

/* ------------------------------------------------------------------------- */

#if SLEPC_VERSION_(2,3,3)

#undef __FUNCT__
#define __FUNCT__ "IPGetOptionsPrefix"
static PetscErrorCode IPGetOptionsPrefix(IP ip,const char *prefix[])
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ip,IP_COOKIE,1);
  PetscValidPointer(prefix,2);
  ierr = PetscObjectGetOptionsPrefix((PetscObject)ip, prefix);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "EPSGetOperators_233"
static PetscErrorCode EPSGetOperators_233(EPS eps, Mat *A, Mat *B)
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
static PetscErrorCode EPSSetOperators_233(EPS eps, Mat A, Mat B)
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
static PetscErrorCode STSetOperators_233(ST st, Mat A, Mat B)
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
static PetscErrorCode EPSSetExtraction(EPS eps,EPSExtraction ext)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_COOKIE,1);
  SETERRQ(PETSC_ERR_SUP,"operation not supported in this SLEPc version");
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "EPSGetExtraction"
static PetscErrorCode EPSGetExtraction(EPS eps,EPSExtraction *ext)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_COOKIE,1);
  SETERRQ(PETSC_ERR_SUP,"operation not supported in this SLEPc version");
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "EPSSetTarget"
static PetscErrorCode EPSSetTarget(EPS eps,PetscScalar target)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_COOKIE,1);
  SETERRQ(PETSC_ERR_SUP,"operation not supported in this SLEPc version");
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "EPSGetTarget"
static PetscErrorCode EPSGetTarget(EPS eps,PetscScalar *target)
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
static PetscErrorCode EPSSetDimensions_233(EPS eps,PetscInt nev,PetscInt ncv,PetscInt mpd)
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
static PetscErrorCode EPSGetDimensions_233(EPS eps, PetscInt *nev,PetscInt *ncv,PetscInt *mpd)
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
static PetscErrorCode SVDSetDimensions_233(SVD svd,PetscInt nev,PetscInt ncv,PetscInt mpd)
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
static PetscErrorCode SVDGetDimensions_233(SVD svd, PetscInt *nev,PetscInt *ncv,PetscInt *mpd)
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

#endif

#endif/*SLEPC4PY_COMPAT_H*/
