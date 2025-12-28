      subroutine platfd(uplatc,uplatm)
c
c     This subroutine sets the platform designator strings that are
c     written on the output and screen files of various codes.
c
c     This subroutine is called by:
c
c       EQPT/eqpt.f
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       None
c
c     Output:
c
c       uplatc = platform category (e.g., UNIX, PC, MAC)
c       uplatm = platform machine (e.g., SPARC, SGI, Pentium, 486PC)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      character*8 uplatc,uplatm
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
c       None
c
c-----------------------------------------------------------------------
c
com   BEGIN_UNIX_DEPENDENT_CODE
com
c       uplatc = 'UNIX'
com
c       uplatm = 'SPARC'
cxx     uplatm = 'SGI'
cxx     uplatm = 'HP-UX'
cxx     uplatm = 'AIX'
cxx     uplatm = 'Ultrix'
com
com   END_UNIX_DEPENDENT_CODE
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
com   BEGIN_PC_DEPENDENT_CODE
com
        uplatc = 'PC'
com
cxx     uplatm = 'Pentium II'
cxx     uplatm = 'K6'
cxx     uplatm = 'K6-2'
        uplatm = 'PC'
cxx     uplatm = 'Pentium Pro'
cxx     uplatm = 'P5'
cxx     uplatm = 'Pentium'
cxx     uplatm = '486PC'
cxx     uplatm = '386PC'
com
com   END_PC_DEPENDENT_CODE
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
com   BEGIN_MAC_DEPENDENT_CODE
com
c       uplatc = 'MAC'
com
c       uplatm = 'MAC'
com
com   END_MAC_DEPENDENT_CODE
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
com   BEGIN_VAX_DEPENDENT_CODE
com
c       uplatc = 'VAX/VMS'
com
c       uplatm = 'VAX'
com
com   END_VAX_DEPENDENT_CODE
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
