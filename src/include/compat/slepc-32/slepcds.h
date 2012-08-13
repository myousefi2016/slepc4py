typedef struct _p_DS* DS;

static PetscClassId DS_CLASSID = 0;

#define DSType            char*
#define DSHEP             "hep"
#define DSNHEP            "nhep"
#define DSGHEP            "ghep"
#define DSGHIEP           "ghiep"
#define DSGNHEP           "gnhep"
#define DSSVD             "svd"
#define DSQEP             "qep"

typedef enum { DS_STATE_RAW,
               DS_STATE_INTERMEDIATE,
               DS_STATE_CONDENSED,
               DS_STATE_TRUNCATED } DSStateType;

typedef enum { DS_MAT_A,
               DS_MAT_B,
               DS_MAT_C,
               DS_MAT_T,
               DS_MAT_D,
               DS_MAT_Q,
               DS_MAT_Z,
               DS_MAT_X,
               DS_MAT_Y,
               DS_MAT_U,
               DS_MAT_VT,
               DS_MAT_W,
               DS_NUM_MAT } DSMatType;

#undef  __FUNCT__
#define __FUNCT__ "DS_NotSupported"
#define DS_NOTSUPPORTED(FUNCT,ARGS) \
static PetscErrorCode FUNCT ARGS \
{ \
  PetscFunctionBegin; \
  SETERRQ(PETSC_COMM_SELF,PETSC_ERR_SUP, \
  #FUNCT"() not supported in this SLEPc version");\
  PetscFunctionReturn(PETSC_ERR_SUP); \
}

DS_NOTSUPPORTED( DSCreate,  (MPI_Comm comm,DS *ds))
DS_NOTSUPPORTED( DSDestroy, (DS *ds)              )
DS_NOTSUPPORTED( DSView                  , (DS ds,...) )
DS_NOTSUPPORTED( DSReset                 , (DS ds,...) )
DS_NOTSUPPORTED( DSSetType               , (DS ds,...) )
DS_NOTSUPPORTED( DSGetType               , (DS ds,...) )
DS_NOTSUPPORTED( DSSetOptionsPrefix      , (DS ds,...) )
DS_NOTSUPPORTED( DSGetOptionsPrefix      , (DS ds,...) )
DS_NOTSUPPORTED( DSSetFromOptions        , (DS ds,...) )
DS_NOTSUPPORTED( DSAllocate              , (DS ds,...) )
DS_NOTSUPPORTED( DSGetLeadingDimension   , (DS ds,...) )
DS_NOTSUPPORTED( DSSetState              , (DS ds,...) )
DS_NOTSUPPORTED( DSGetState              , (DS ds,...) )
DS_NOTSUPPORTED( DSSetDimensions         , (DS ds,...) )
DS_NOTSUPPORTED( DSGetDimensions         , (DS ds,...) )
DS_NOTSUPPORTED( DSSetMethod             , (DS ds,...) )
DS_NOTSUPPORTED( DSGetMethod             , (DS ds,...) )
DS_NOTSUPPORTED( DSSetCompact            , (DS ds,...) )
DS_NOTSUPPORTED( DSGetCompact            , (DS ds,...) )
DS_NOTSUPPORTED( DSSetExtraRow           , (DS ds,...) )
DS_NOTSUPPORTED( DSGetExtraRow           , (DS ds,...) )
DS_NOTSUPPORTED( DSSetRefined            , (DS ds,...) )
DS_NOTSUPPORTED( DSGetRefined            , (DS ds,...) )
DS_NOTSUPPORTED( DSTruncate              , (DS ds,...) )
DS_NOTSUPPORTED( DSUpdateExtraRow        , (DS ds,...) )
