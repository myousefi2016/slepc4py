# -----------------------------------------------------------------------------

cdef class Sys:

    @classmethod
    def getVersion(cls, patch=False, devel=False,
                   date=False, author=False):
        cdef int cmajor = SLEPC_VERSION_MAJOR
        cdef int cminor = SLEPC_VERSION_MINOR
        cdef int cmicro = SLEPC_VERSION_SUBMINOR
        cdef int cpatch = SLEPC_VERSION_PATCH
        cdef int cdevel = not SLEPC_VERSION_RELEASE
        cdef const_char *cdate       = SLEPC_VERSION_DATE
        cdef const_char *cauthorinfo = SLEPC_AUTHOR_INFO
        version = (cmajor, cminor, cmicro)
        out = version
        if patch or devel or date or author:
            out = [version]
            if patch:
                out.append(cpatch)
            if devel:
                out.append(<bint>cdevel)
            if date:
                out.append(bytes2str(cdate))
            if author:
                author = bytes2str(cauthorinfo).split('\n')
                author = [s.strip() for s in author if s]
                out.append(author)
        return tuple(out)

    @classmethod
    def getVersionInfo(cls):
        cdef int cmajor = SLEPC_VERSION_MAJOR
        cdef int cminor = SLEPC_VERSION_MINOR
        cdef int cmicro = SLEPC_VERSION_SUBMINOR
        cdef int cpatch = SLEPC_VERSION_PATCH
        cdef int crelease = SLEPC_VERSION_RELEASE
        cdef const_char *cdate       = SLEPC_VERSION_DATE
        cdef const_char *cauthorinfo = SLEPC_AUTHOR_INFO
        author = bytes2str(cauthorinfo).split('\n')
        author = [s.strip() for s in author if s]
        return dict(major      = cmajor,
                    minor      = cminor,
                    subminor   = cmicro,
                    patch      = cpatch,
                    release    = <bint>crelease,
                    date       = bytes2str(cdate),
                    authorinfo = author)

# -----------------------------------------------------------------------------
