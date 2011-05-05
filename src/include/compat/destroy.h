#ifndef _COMPAT_SLEPC_DESTROY_H
#define _COMPAT_SLEPC_DESTROY_H

#define DESTROY(PetscType, COOKIE)                      \
PETSC_STATIC_INLINE                                     \
PetscErrorCode PetscType##Destroy_new(PetscType *obj)   \
{                                                       \
  PetscType      tmp = 0;                               \
  PetscErrorCode ierr;                                  \
  PetscFunctionBegin;                                   \
  PetscValidPointer((obj),1);                           \
  if (!(*(obj))) PetscFunctionReturn(0);                \
  tmp = *(obj); *(obj) = 0;                             \
  if (COOKIE == -1)                                     \
    {PetscValidPointer(tmp,1); }                        \
  else if (COOKIE == PETSC_OBJECT_COOKIE)               \
    {PetscValidHeader(tmp,1);}                          \
  else                                                  \
    {PetscValidHeaderSpecific(tmp,COOKIE,1);}           \
  ierr = PetscType##Destroy(tmp);CHKERRQ(ierr);         \
  PetscFunctionReturn(0);                               \
}                                                       \
/**/
#undef  __FUNCT__
#define __FUNCT__ "User provided function\0:Destroy"
DESTROY(PetscObject , PETSC_OBJECT_COOKIE )
DESTROY(IP          , IP_COOKIE           )
DESTROY(ST          , ST_COOKIE           )
DESTROY(EPS         , EPS_COOKIE          )
DESTROY(SVD         , SVD_COOKIE          )
DESTROY(QEP         , QEP_COOKIE          )

#undef PetscObjetDestroy
#undef IPDestroy
#undef STDestroy
#undef EPSDestroy
#undef SVDDestroy
#undef QEPDestroy

#define PetscObjectDestroy  PetscObjectDestroy_new
#define IPDestroy           IPDestroy_new
#define STDestroy           STDestroy_new
#define EPSDestroy          EPSDestroy_new
#define SVDDestroy          SVDDestroy_new
#define QEPDestroy          QEPDestroy_new

#endif/*_COMPAT_SLEPC_DESTROY_H*/
