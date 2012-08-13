#define IPBILINEAR     "bilinear"
#define IPSESQUILINEAR "sesquilinear"
#define IPINDEFINITE   "indefinite"
#define EPSRQCG "rqcg"

#define EPS_GHIEP ((EPSProblemType)-1)

typedef enum { EPS_ORTH_I=1,
               EPS_ORTH_B,
               EPS_ORTH_BOPT } EPSOrthType;

#undef  __FUNCT__
#define __FUNCT__ "EPS_NotSupported"
#define EPS_NOTSUPPORTED(FUNCT,ARGS) \
static PetscErrorCode FUNCT ARGS \
{ \
  PetscFunctionBegin; \
  SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP, \
  #FUNCT"() not supported in this SLEPc version");\
  PetscFunctionReturn(PETSC_ERR_SUP); \
}
EPS_NOTSUPPORTED( EPSIsPositive            , (EPS eps, ...) )
EPS_NOTSUPPORTED( EPSGetDS                 , (EPS eps, ...) )
EPS_NOTSUPPORTED( EPSSetDS                 , (EPS eps, ...) )
EPS_NOTSUPPORTED( EPSKrylovSchurSetRestart , (EPS eps, ...) )
EPS_NOTSUPPORTED( EPSKrylovSchurGetRestart , (EPS eps, ...) )
EPS_NOTSUPPORTED( EPSRQCGSetReset          , (EPS eps, ...) )
EPS_NOTSUPPORTED( EPSRQCGGetReset          , (EPS eps, ...) )

