# --------------------------------------------------------------------

cdef extern from "mpi.h":

    MPI_Comm MPI_COMM_NULL
    MPI_Comm MPI_COMM_SELF
    MPI_Comm MPI_COMM_WORLD

cdef extern from "petsc.h":

    MPI_Comm PETSC_COMM_SELF
    MPI_Comm PETSC_COMM_WORLD

    MPI_Comm PETSC_COMM_DEFAULT "PETSC_COMM_WORLD"

# --------------------------------------------------------------------

cdef inline MPI_Comm def_Comm(object comm,
                              MPI_Comm deft) except *:
    if comm is None: return deft
    return (<Comm?>comm).comm

# --------------------------------------------------------------------
