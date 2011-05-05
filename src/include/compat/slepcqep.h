#ifndef _COMPAT_SLEPC_QEP_H
#define _COMPAT_SLEPC_QEP_H
PETSC_EXTERN_CXX_BEGIN

#undef  __FUNCT__
#define __FUNCT__ "QEP???"

#define SlepcQEP_ERR_SUP                                            \
  PetscFunctionBegin;                                               \
  SETERRQ(PETSC_ERR_SUP,"QEP not supported in this PETSc version"); \
  PetscFunctionReturn(PETSC_ERR_SUP);

static PetscCookie QEP_COOKIE = 0;

/*S
     QEP - Abstract SLEPc object that manages all the quadratic eigenvalue 
     problem solvers.

   Level: beginner

.seealso:  QEPCreate()
S*/
typedef struct _p_QEP* QEP;

/*E
    QEPType - String with the name of a quadratic eigensolver

   Level: beginner

.seealso: QEPSetType(), QEP
E*/
#define QEPType      char*
#define QEPLINEAR    "linear"
#define QEPQARNOLDI  "qarnoldi"

/*E
    QEPProblemType - determines the type of the quadratic eigenproblem

    Level: intermediate

.seealso: QEPSetProblemType(), QEPGetProblemType()
E*/
typedef enum { QEP_GENERAL=1,
               QEP_HERMITIAN,   /* M, C, K  Hermitian */
               QEP_GYROSCOPIC   /* M, K  Hermitian, M>0, C skew-Hermitian */
             } QEPProblemType;

/*E
    QEPWhich - determines which part of the spectrum is requested

    Level: intermediate

.seealso: QEPSetWhichEigenpairs(), QEPGetWhichEigenpairs()
E*/
typedef enum { QEP_LARGEST_MAGNITUDE=1,
               QEP_SMALLEST_MAGNITUDE,
               QEP_LARGEST_REAL,
               QEP_SMALLEST_REAL,
               QEP_LARGEST_IMAGINARY,
               QEP_SMALLEST_IMAGINARY } QEPWhich;

static PetscErrorCode QEPCreate(MPI_Comm c,QEP *o){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPDestroy(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPSetType(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPGetType(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPSetProblemType(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPGetProblemType(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPSetOperators(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPGetOperators(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPSetFromOptions(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPSetUp(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPSolve(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPView(QEP o,...){SlepcQEP_ERR_SUP}

static PetscErrorCode QEPSetIP(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPGetIP(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPSetTolerances(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPGetTolerances(QEP o,...){SlepcQEP_ERR_SUP}

#if 0
static PetscErrorCode QEPSetConvergenceTest(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPDefaultConverged(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPAbsoluteConverged(QEP o,...){SlepcQEP_ERR_SUP}
#endif

static PetscErrorCode QEPSetDimensions(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPGetDimensions(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPSetScaleFactor(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPGetScaleFactor(QEP o,...){SlepcQEP_ERR_SUP}

static PetscErrorCode QEPGetConverged(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPGetEigenpair(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPComputeRelativeError(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPComputeResidualNorm(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPGetErrorEstimate(QEP o,...){SlepcQEP_ERR_SUP}

static PetscErrorCode QEPGetIterationNumber(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPGetOperationCounters(QEP o,...){SlepcQEP_ERR_SUP}

static PetscErrorCode QEPSetInitialSpace(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPSetInitialSpaceLeft(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPSetWhichEigenpairs(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPGetWhichEigenpairs(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPSetLeftVectorsWanted(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPGetLeftVectorsWanted(QEP o,...){SlepcQEP_ERR_SUP}

#if 0
static PetscErrorCode QEPMonitorSet(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPGetMonitorContext(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPMonitorAll(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPMonitorFirst(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPMonitorConverged(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPMonitorLG(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPMonitorLGAll(QEP o,...){SlepcQEP_ERR_SUP}
#endif
static PetscErrorCode QEPMonitorCancel(QEP o,...){SlepcQEP_ERR_SUP}

static PetscErrorCode QEPSetTrackAll(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPGetTrackAll(QEP o,...){SlepcQEP_ERR_SUP}

static PetscErrorCode QEPSetOptionsPrefix(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPAppendOptionsPrefix(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPGetOptionsPrefix(QEP o,...){SlepcQEP_ERR_SUP}

/*E
    QEPConvergedReason - reason an eigensolver was said to 
         have converged or diverged

    Level: beginner

.seealso: QEPSolve(), QEPGetConvergedReason(), QEPSetTolerances()
E*/
typedef enum {/* converged */
              QEP_CONVERGED_TOL                =  2,
              /* diverged */
              QEP_DIVERGED_ITS                 = -3,
              QEP_DIVERGED_BREAKDOWN           = -4,
              QEP_CONVERGED_ITERATING          =  0} QEPConvergedReason;

static PetscErrorCode QEPGetConvergedReason(QEP o,...){SlepcQEP_ERR_SUP}

#if 0
static PetscErrorCode QEPSetEigenvalueComparison(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPSortEigenvalues(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPSortEigenvaluesReal(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPCompareEigenvalues(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPSortDenseSchur(QEP o,...){SlepcQEP_ERR_SUP}
#endif

#if 0
static PetscErrorCode QEPRegister(const char* n,...){SlepcQEP_ERR_SUP}
#if defined(PETSC_USE_DYNAMIC_LIBRARIES)
#define QEPRegisterDynamic(a,b,c,d) QEPRegister(a,b,c,0)
#else
#define QEPRegisterDynamic(a,b,c,d) QEPRegister(a,b,c,d)
#endif
static PetscErrorCode QEPRegisterDestroy(QEP o,...){SlepcQEP_ERR_SUP}
#endif

/* --------- options specific to particular eigensolvers -------- */
#if 0
static PetscErrorCode QEPLinearSetCompanionForm(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPLinearGetCompanionForm(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPLinearSetExplicitMatrix(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPLinearGetExplicitMatrix(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPLinearSetEPS(QEP o,...){SlepcQEP_ERR_SUP}
static PetscErrorCode QEPLinearGetEPS(QEP o,...){SlepcQEP_ERR_SUP}
#endif

PETSC_EXTERN_CXX_END
#endif /* _COMPAT_SLEPC_QEP_H */
