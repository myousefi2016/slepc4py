#ifndef _COMPAT_SLEPC_RESET_H
#define _COMPAT_SLEPC_RESET_H

#define RESET(SlepcType, COOKIE)                                \
PETSC_STATIC_INLINE                                             \
PetscErrorCode SlepcType##Reset(SlepcType obj)                  \
{                                                               \
  PetscFunctionBegin;                                           \
  PetscValidHeaderSpecific(obj,COOKIE,1);                       \
  SETERRQ(PETSC_ERR_SUP,__FUNCT__"() "                          \
          "not supported in this SLEPc version");               \
  PetscFunctionReturn(0);                                       \
}                                                               \

#undef  __FUNCT__
#define __FUNCT__ "IPReset"
RESET(IP, IP_COOKIE)

#undef  __FUNCT__
#define __FUNCT__ "STReset"
RESET(ST, ST_COOKIE)

#undef  __FUNCT__
#define __FUNCT__ "EPSReset"
RESET(EPS, EPS_COOKIE)

#undef  __FUNCT__
#define __FUNCT__ "SVDReset"
RESET(SVD, SVD_COOKIE)

#undef  __FUNCT__
#define __FUNCT__ "QEPReset"
RESET(QEP, QEP_COOKIE)

#endif/*_COMPAT_SLEPC_RESET_H*/
