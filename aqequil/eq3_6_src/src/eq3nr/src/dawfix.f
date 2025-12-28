      subroutine dawfix(aamatr,cdrs,eps100,gmmatr,iindx1,iodb,
     $ irdxc3,jflag,jjndex,kbt,kkndex,kmax,narn1,nbasp,nbtmax,ncosp,
     $ ndrs,ndrsmx,ndrsr,nelect,nhydr,nodbmx,no2gaq,noutpt,nstmax,
     $ qawfix,uspec)
c
c     This subroutine determines if the activity of water is directly
c     or indirectly fixed. For example, the user might directly specify
c     the activity of water. Alternatively, a set of solubility
c     constraints (e.g., gysum + anhydrite) might indirectly fix it.
c     If the activity of water is fixed directly or indirectly, the
c     the flag 'qawfix' is returned with a value of '.true.'.
c
c     This subroutine is called by:
c
c       EQ3NR/arrset.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c
c     Principal output:
c
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer kmax,nbtmax,ndrsmx,nodbmx,nstmax
c
      integer noutpt
c
      integer iindx1(kmax),iodb(nodbmx),jjndex(nbtmax),jflag(nstmax),
     $ kkndex(nbtmax),ncosp(nbtmax),nbasp(nbtmax),ndrs(ndrsmx),
     $ ndrsr(2,nstmax)
c
      integer irdxc3,kbt,narn1,nelect,nhydr,no2gaq
c
      logical qawfix
c
      character*48 uspec(nstmax)
c
      real(8) aamatr(kmax,kmax),cdrs(ndrsmx),gmmatr(kmax,kmax)
c
      real(8) eps100
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ibt,icol,ihydr,io2gaq,irow,irow1,irow2,iwater,j2,jcol,
     $ jcol1,jcol2,krow,nb,nb1,ns,ns1,ns2
c
      integer ilnobl
c
      logical qldep,qx
c
      real(8) coefdr
c
c-----------------------------------------------------------------------
c
      qawfix = .false.
c
c     Check for directly specified water activity. Note that the
c     following assumes that water is the first aqueous species
c     (e.g., it has the species index narn1).
c
      if (jflag(narn1) .eq. 16) then
        qawfix = .true.
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Build a matrix to test for water activity fixed by simultaneous
c     equilibria.
c
      ibt = 0
      iwater = 0
      ihydr = 0
      io2gaq = 0
c
      do krow = 1,kbt
        nb = iindx1(krow)
        ns = nbasp(nb)
        kkndex(nb) = 0
        qx = .false.
        if (ns .eq. narn1) then
          qx = .true.
        elseif (ns.eq.nelect .or. ns.eq.no2gaq) then
          qx = irdxc3 .ne. 0
        else
          qx = jflag(ns).eq.25 .or. jflag(ns).eq.27
        endif
c
        if (qx) then
          ibt = ibt + 1
          jjndex(ibt) = nb
          kkndex(nb) = 1
          if (ns .eq. narn1) iwater = ibt
          if (ns .eq. nhydr) ihydr = ibt
          if (ns .eq. no2gaq) io2gaq = ibt
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Test the matrix dimension. Unless it is greater than or equal to
c     two, the activity of water is not fixed indirectly.
c
      if (ibt .lt. 2) go to 999
c
c     Build the matrix.
c
      do icol = 1,ibt
        do irow = 1,ibt
          aamatr(irow,icol) = 0.
        enddo
      enddo
c
      do irow = 1,ibt
        nb = jjndex(irow)
        ns = nbasp(nb)
        if (ns .eq. narn1) then
c
c         Water.
c
          aamatr(irow,irow) = 1.0
        elseif (jflag(ns) .eq. 27) then
c
c         Aqueous homogenous equilibrium.
c
          do icol = 1,ibt
            nb1 = jjndex(icol)
            ns1 = nbasp(nb1)
c
c           Calling sequence substitutions:
c             ns1 for nse
c
            aamatr(irow,icol)
     $       = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns1,ns,nstmax)
          enddo
        elseif (jflag(ns) .eq. 25) then
c
c         Heterogeneous equilibrium.
c
          ns2 = ncosp(nb)
          do icol = 1,ibt
            ns1 = jjndex(icol)
c
c           Calling sequence substitutions:
c             ns1 for nse
c             ns2 for ns
c
            aamatr(irow,icol)
     $       = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns1,ns2,nstmax)
          enddo
        elseif (irdxc3 .lt. 0) then
c
c         Eh (an input pe- has been converted to Eh).
c
          aamatr(irow,iwater) = -2.
          aamatr(irow,io2gaq) = 1.
          if (ihydr .gt. 0) aamatr(irow,ihydr) = 4.
        else
c
c         Aqueous redox couple.
c
          ns2 = irdxc3
          do icol = 1,ibt
            ns1 = jjndex(icol)
            if (kkndex(ns1) .ge. 1) then
c
c             Calling sequence substitutions:
c               ns1 for nse
c               ns2 for ns
c
              aamatr(irow,icol)
     $        = coefdr(cdrs,ndrs,ndrsmx,ndrsr,ns1,ns2,nstmax)
            endif
          enddo
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(3) .ge. 3) then
        write (noutpt,1000)
 1000   format(/10x,'--- Matrix to Test for Fixed Activity of Water',
     $  ' ---',/)
        do irow = 1,ibt
          write (noutpt,1010) (aamatr(irow,icol), icol = 1,ibt)
 1010     format(2x,10(f7.2,2x))
        enddo
        write (noutpt,1020)
 1020   format(/1x)
c
        do irow = 1,ibt
          nb = jjndex(irow)
          ns = nbasp(nb)
          j2 = ilnobl(uspec(ns)(1:24))
          write (noutpt,1030) irow,uspec(ns)(1:j2)
 1030     format(2x,i3,2x,a)
c
          if (ns .eq. narn1) then
            write (noutpt,1040)
 1040       format(/10x,'Mole fraction definition',/)
          elseif (jflag(ns) .eq. 27) then
c
c           Calling sequence substitutions:
c             noutpt for nf
c
            call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns,nstmax,uspec)
          elseif (jflag(ns) .eq. 25) then
            ns2 = ncosp(ns)
c
c           Calling sequence substitutions:
c             noutpt for nf
c             ns2 for ns
c
            call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns2,nstmax,uspec)
          elseif (irdxc3 .lt. 0) then
            write (noutpt,1050)
 1050       format(/10x,'2 H2O(l) = 4 H+ + 4 e- + O2(g)',/)
          else
            ns2 = irdxc3
c
c           Calling sequence substitutions:
c             noutpt for nf
c             ns2 for ns
c
            call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns2,nstmax,uspec)
          endif
          write (noutpt,1020)
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make a copy of the matrix.
c
      do irow = 1,ibt
        do jcol = 1,ibt
          gmmatr(irow,jcol) = aamatr(irow,jcol)
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Test for linear dependence omitting the water row and column.
c     The original matrix is overwritten in the process.
c
      irow1 = 2
      irow2 = ibt
      jcol1 = 2
      jcol2 = ibt
      call lindep(aamatr,eps100,irow1,irow2,jcol1,jcol2,kmax,qldep)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qldep) then
c
c       Have linear dependence omitting the water column. Test for
c       linear dependence with the water column (still omitting the
c       water row).
c
c       Restore the matrix from the copy.
c
        do irow = 1,ibt
          do jcol = 1,ibt
            aamatr(irow,jcol) = gmmatr(irow,jcol)
          enddo
        enddo
c
        jcol1 = 1
        call lindep(aamatr,eps100,irow1,irow2,jcol1,jcol2,kmax,qldep)
        qawfix = .not.qldep
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
