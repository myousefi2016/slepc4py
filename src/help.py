# Author:  Lisandro Dalcin
# Contact: dalcinl@gmail.com

"""
Command line access to the SLEPc Options Database.

This module intends to provide command line access to the SLEPc
Options Database. It outputs a listing of the many SLEPc options
indicating option names, default values and descriptions. If you have
Python 2.4 and above, then you can issue at the command line::

  $ python -m slepc4py.help [eps|svd|st|ip] [<slepc-option-list>]

"""

def help(args=None):
    import sys
    # program name
    try:
        prog = sys.argv[0]
    except Exception:
        prog = getattr(sys, 'executable', 'python')
    # arguments
    if args is None:
        args = sys.argv[1:]
    elif isinstance(args, str):
        args = args.split()
    else:
        args = [str(a) for a in args]
    # initialization
    import slepc4py
    slepc4py.init([prog, '-help'] + args)
    from slepc4py import SLEPc
    # and finally ...
    COMM = SLEPc.COMM_SELF
    if 'eps' in args:
        eps = SLEPc.EPS().create(comm=COMM)
        eps.setFromOptions()
        eps.destroy()
        del eps
    if 'svd' in args:
        svd = SLEPc.SVD().create(comm=COMM)
        svd.setFromOptions()
        svd.destroy()
        del svd
    if 'pep' in args:
        pep = SLEPc.PEP().create(comm=COMM)
        pep.setFromOptions()
        pep.destroy()
        del pep
    if 'nep' in args:
        nep = SLEPc.NEP().create(comm=COMM)
        nep.setFromOptions()
        nep.destroy()
        del nep
    if 'mfn' in args:
        mfn = SLEPc.MFN().create(comm=COMM)
        mfn.setFromOptions()
        mfn.destroy()
        del mfn
    if 'st' in args:
        st = SLEPc.ST().create(comm=COMM)
        st.setFromOptions()
        st.destroy()
        del st
    if 'bv' in args:
        bv = SLEPc.BV().create(comm=COMM)
        bv.setFromOptions()
        bv.destroy()
        del bv

if __name__ == '__main__':
    help()
