/* Author:  Lisandro Dalcin   */
/* Contact: dalcinl@gmail.com */

/* ---------------------------------------------------------------- */

%include petsc4py/petsc4py.i

/* ---------------------------------------------------------------- */

%header %{#include "slepc4py/slepc4py.h"%}
%init   %{import_slepc4py();%}

%define SWIG_TYPECHECK_SLEPC_ST  550 %enddef
%define SWIG_TYPECHECK_SLEPC_BV  551 %enddef
%define SWIG_TYPECHECK_SLEPC_DS  552 %enddef
%define SWIG_TYPECHECK_SLEPC_EPS 553 %enddef
%define SWIG_TYPECHECK_SLEPC_SVD 554 %enddef
%define SWIG_TYPECHECK_SLEPC_PEP 555 %enddef
%define SWIG_TYPECHECK_SLEPC_NEP 556 %enddef
%define SWIG_TYPECHECK_SLEPC_MFN 557 %enddef
%define SWIG_TYPECHECK_SLEPC_FN  558 %enddef

%define %slepc4py_objt(Pkg, PyType, Type, CODE, OBJECT_NULL)
%petsc4py_objt(Pkg, PyType, Type, CODE, OBJECT_NULL)
%enddef /* %slepc4py_typemap */

/* ---------------------------------------------------------------- */

%slepc4py_objt( Slepc , ST ,  ST ,  SLEPC_ST ,  PETSC_NULL )
%slepc4py_objt( Slepc , BV ,  BV ,  SLEPC_BV ,  PETSC_NULL )
%slepc4py_objt( Slepc , DS ,  DS ,  SLEPC_DS ,  PETSC_NULL )
%slepc4py_objt( Slepc , EPS , EPS , SLEPC_EPS , PETSC_NULL )
%slepc4py_objt( Slepc , SVD , SVD , SLEPC_SVD , PETSC_NULL )
%slepc4py_objt( Slepc , PEP , PEP , SLEPC_PEP , PETSC_NULL )
%slepc4py_objt( Slepc , NEP , NEP , SLEPC_NEP , PETSC_NULL )
%slepc4py_objt( Slepc , MFN , MFN , SLEPC_MFN , PETSC_NULL )
%slepc4py_objt( Slepc , FN ,  FN ,  SLEPC_FN ,  PETSC_NULL )

/* ---------------------------------------------------------------- */

/*
 * Local Variables:
 * mode: C
 * End:
 */
