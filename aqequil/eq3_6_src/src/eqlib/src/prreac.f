      subroutine prreac(cdrs,ndrs,ndrsmx,ndrsr,nf,ns,nstmax,uspec)
c
c     This subroutine writes the n-th reaction in a set on the file
c     whose unit number isf nf.
c
c     This subroutine is called by:
c
c       EQLIB/alters.f
c       EQLIB/echolk.f
c       EQLIB/switch.f
c       EQLIB/swtchb.f
c       EQLIB/swtchk.f
c       EQ3NR/arrsim.f
c       EQ3NR/dawfix.f
c       EQ3NR/echox.f
c       EQ3NR/eq3nr.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       cdrs   = the array of reaction coefficients
c       ndrs   = the array of corresponding species indices
c       ndrsmx = dimension of the cdrs and ndrs arrays
c       ndrsr  = aray giving the range in the cdrs and ndrs arrays
c                  corresponding to the reaction for a given species
c       nf     = the unit number of the file to write on
c       ns     = the index of the species whose reaction is to
c                  be written
c       uspec  = array of species names
c
c     Principal output:
c
c       None
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ndrsmx,nstmax
c
      integer ndrs(ndrsmx),ndrsr(2,nstmax)
      integer nf,ns
c
      character*48 uspec(nstmax)
c
      real*8 cdrs(ndrsmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jlen,n,nr1,nr2,nss,nt
c
      logical qdtach,qfirst
c
      character*56 uspn56
c
      real*8 cx
c
c-----------------------------------------------------------------------
c
      write (nf,1000)
 1000 format(1x)
c
      nr1 = ndrsr(1,ns)
      nr2 = ndrsr(2,ns)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      nt = nr2 - nr1 + 1
      if (nt .lt. 2) then
c
c       Calling sequence substitutions:
c         uspec(ns) for unam48
c
        call fmspnx(jlen,uspec(ns),uspn56)
c
        write (nf,1010) uspn56(1:jlen)
 1010   format(3x,a,' is a strict basis species and has no reaction.',/)
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      qdtach = .false.
      do n = nr1 + 1,nr2
        nss = ndrs(n)
        if (nss .eq. 0) then
          qdtach = .true.
          go to 100
        endif
      enddo
  100 continue
c
      if (qdtach) then
c
c       Calling sequence substitutions:
c         uspec(ns) for unam48
c
        call fmspnx(jlen,uspec(ns),uspn56)
c
        write (nf,1020) uspn56(1:jlen)
 1020   format(3x,a,' is a detached auxiliary basis species and',
     $  ' in effect',/3x,'has no reaction.',/)
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the coefficients and names of the reactants.
c
      qfirst = .true.
      do n = nr1,nr2
        cx = -cdrs(n)
        nss = ndrs(n)
        if (cx .gt. 0.) then
c
c         Calling sequence substitutions:
c           uspec(nss) for unam48
c
          call fmspnx(jlen,uspec(nss),uspn56)
c
          if (qfirst) then
            write (nf,1030) cx,uspn56(1:jlen)
 1030       format(4x,f7.3,2x,a)
            qfirst = .false.
          else
            write (nf,1050) cx,uspn56(1:jlen)
 1050       format(2x,'+ ',f7.3,2x,a)
          endif
        endif
      enddo
c
      write (nf,1070)
 1070 format(10x,'==')
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the coefficients and names of the products.
c
      qfirst = .true.
      do n = nr1,nr2
        cx = cdrs(n)
        nss = ndrs(n)
        if (cx .gt. 0.) then
c
c         Calling sequence substitutions:
c           uspec(nss) for unam48
c
          call fmspnx(jlen,uspec(nss),uspn56)
c
          if (qfirst) then
            write (nf,1030) cx,uspn56(1:jlen)
            qfirst = .false.
          else
            write (nf,1050) cx,uspn56(1:jlen)
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
