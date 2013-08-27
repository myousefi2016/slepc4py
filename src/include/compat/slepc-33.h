#define DSGetDimensions(ds,n,m,l,k,t) ((*(t))=0,DSGetDimensions(ds,n,m,l,k))

#undef  __FUNCT__
#define __FUNCT__ "STSetOperators_Compat"
static PetscErrorCode STSetOperators_Compat(ST st,PetscInt n, Mat mats[])
{
  Mat            A=0,B=0;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  if (n>=1) A = mats[0];
  if (n>=2) B = mats[1];
  ierr = STSetOperators(st,A,B);CHKERRQ(ierr);
  PetscFunctionReturn(0);
}
#define STSetOperators STSetOperators_Compat

#undef  __FUNCT__
#define __FUNCT__ "STGetNumMatrices"
static PetscErrorCode STGetNumMatrices(ST st, PetscInt *n)
{
  Mat            A=0,B=0;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  PetscValidPointer(n,2);
  ierr = STGetOperators(st,&A,&B);CHKERRQ(ierr);
  *n = 0; if (A) (*n)++; if (B) (*n)++;
  PetscFunctionReturn(0);
}

#undef  __FUNCT__
#define __FUNCT__ "STGetOperators_Compat"
static PetscErrorCode STGetOperators_Compat(ST st,PetscInt k,Mat *mat)
{
  Mat            A=0,B=0;
  PetscErrorCode ierr;
  PetscFunctionBegin;
  PetscValidHeaderSpecific(st,ST_CLASSID,1);
  ierr = STGetOperators(st,&A,&B);CHKERRQ(ierr);
  *mat = 0; if (k==0) *mat = A; if (k==1) *mat = B;
  PetscFunctionReturn(0);
}
#define STGetOperators STGetOperators_Compat
