/* Author:  Lisandro Dalcin   */
/* Contact: dalcinl@gmail.com */

/* ---------------------------------------------------------------- */

%include petsc4py/petsc4py.i

/* ---------------------------------------------------------------- */

%header %{#include "slepc4py/slepc4py.h"%}
%init   %{import_slepc4py();%}

%define SWIG_TYPECHECK_SLEPC_ST  550 %enddef
%define SWIG_TYPECHECK_SLEPC_IP  551 %enddef
%define SWIG_TYPECHECK_SLEPC_EPS 552 %enddef
%define SWIG_TYPECHECK_SLEPC_ST  553 %enddef

%define %slepc4py_objt(Pkg, PyType, Type, CODE, OBJECT_NULL)
%petsc4py_objt(Pkg, PyType, Type, CODE, OBJECT_NULL)
%enddef /* %slepc4py_typemap */

/* ---------------------------------------------------------------- */

%slepc4py_objt( Slepc , ST ,  ST ,  SLEPC_ST ,  PETSC_NULL )
%slepc4py_objt( Slepc , IP ,  IP ,  SLEPC_IP ,  PETSC_NULL )
%slepc4py_objt( Slepc , EPS , EPS , SLEPC_EPS , PETSC_NULL )
%slepc4py_objt( Slepc , SVD , SVD , SLEPC_SVD , PETSC_NULL )

/* ---------------------------------------------------------------- */

/*
 * Local Variables:
 * mode: C
 * End:
 */
