#ifndef SLEPC4PY_COMPAT_H
#define SLEPC4PY_COMPAT_H

#include <slepc.h>

#if PETSC_VERSION_(3,2,0)
#define PetscShell PetscFwk
#endif

/* ------------------------------------------------------------------------- */
#if SLEPC_VERSION_(3,1,0)
#include <slepcqep.h>
#elif SLEPC_VERSION_(3,0,0)
#include <slepcsvd.h>
#include "compat/slepcqep.h"
#endif

#undef  __FUNCT__
#define __FUNCT__ "SlepcInitializePackage"
static PetscErrorCode SlepcInitializePackage(const char path[])
{
  PetscErrorCode ierr;
  PetscFunctionBegin;
#if SLEPC_VERSION_(3,0,0)
  ierr = PetscCookieRegister("Quadratic Eigenproblem Solver",&QEP_COOKIE);CHKERRQ(ierr);
#else
  ierr = 0; CHKERRQ(ierr);
#endif
  PetscFunctionReturn(0);
}
/* ------------------------------------------------------------------------- */


/* ------------------------------------------------------------------------- */
#if SLEPC_VERSION_(3,1,0) || SLEPC_VERSION_(3,0,0)

#include "compat/destroy.h"
#include "compat/reset.h"

#define PetscShell void*

#define PetscBool    PetscTruth

#define PetscClassId PetscCookie
#define ST_CLASSID   ST_COOKIE
#define IP_CLASSID   IP_COOKIE
#define EPS_CLASSID  EPS_COOKIE
#define SVD_CLASSID  SVD_COOKIE
#define QEP_CLASSID  QEP_COOKIE

#define EPS_ALL ((EPSWhich)EPS_WHICH_USER)

#undef __FUNCT__
#define __FUNCT__ "EPSSetInterval"
static PetscErrorCode EPSSetInterval(EPS eps,PetscReal inta,PetscReal intb)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  SETERRQ(PETSC_ERR_SUP,"operation not supported in this SLEPc version");
  PetscFunctionReturn(0);
}
#undef __FUNCT__
#define __FUNCT__ "EPSGetInterval"
static PetscErrorCode EPSGetInterval(EPS eps,PetscReal* inta,PetscReal* intb)
{
  PetscFunctionBegin;
  PetscValidHeaderSpecific(eps,EPS_CLASSID,1);
  PetscValidPointer(inta,2);
  PetscValidPointer(intb,3);
  SETERRQ(PETSC_ERR_SUP,"operation not supported in this SLEPc version");
  PetscFunctionReturn(0);
}

#if SLEPC_VERSION_(3,0,0)
/**/
#define EPSGD  "gd"
#define EPSJD  "jd"
/**/
#define EPS_HARMONIC_RELATIVE EPS_HARMONIC
#define EPS_HARMONIC_RIGHT    EPS_HARMONIC
#define EPS_HARMONIC_LARGEST  EPS_HARMONIC
/**/
#endif

#if SLEPC_VERSION_(3,0,0)
/**/
#define STSINVERT  STSINV
#define STPRECOND  "precond"
/**/
#define ST_MATMODE_COPY    STMATMODE_COPY
#define ST_MATMODE_INPLACE STMATMODE_INPLACE
#define ST_MATMODE_SHELL   STMATMODE_SHELL
/**/
#endif


#define IPOrthogType  IPOrthogonalizationType
#define IP_ORTHOG_MGS IP_ORTH_MGS
#define IP_ORTHOG_CGS IP_ORTH_CGS

#define IPOrthogRefineType        IPOrthogonalizationRefinementType
#define IP_ORTHOG_REFINE_NEVER    IP_ORTH_REFINE_NEVER
#define IP_ORTHOG_REFINE_IFNEEDED IP_ORTH_REFINE_IFNEEDED
#define IP_ORTHOG_REFINE_ALWAYS   IP_ORTH_REFINE_ALWAYS

#undef __FUNCT__
#define __FUNCT__ "IPGetMatrix"
static PetscErrorCode IPGetMatrix(IP ip, Mat *mat)
{
  IPBilinearForm form;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ip,IP_COOKIE,1);
  PetscValidPointer(mat,2);
  ierr = IPGetBilinearForm(ip,mat,&form);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#undef __FUNCT__
#define __FUNCT__ "IPSetMatrix"
static PetscErrorCode IPSetMatrix(IP ip, Mat mat)
{
  IPBilinearForm form;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(ip,IP_COOKIE,1);
  ierr = IPGetBilinearForm(ip,0,&form);CHKERRQ(ierr);
  ierr = IPSetBilinearForm(ip,mat,form);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}

#if SLEPC_VERSION_(3,0,0)
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
  if (nds > 0) {SETERRQ(PETSC_ERR_SUP,"operation not supported in this SLEPc version");}
  ierr = IPOrthogonalize(ip,n,which,V,v,H,norm,lindep,
                         PETSC_NULL,0);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#define IPOrthogonalize IPOrthogonalize_300
/**/
#endif


#if SLEPC_VERSION_(3,0,0)
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

#if SLEPC_VERSION_(3,0,0)
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

#endif
/* ------------------------------------------------------------------------- */


#endif/*SLEPC4PY_COMPAT_H*/
