      subroutine ptztab(iopr,narn1,narn2,natmax,nmutmx,nmux,nmxi,
     $ nmxmax,nmxx,noprmx,noutpt,nsltmx,nslx,nstmax,nsxi,nsxmax,
     $ nsxx,uspec)
c
c     This subroutine tabulates the species combinations corresponding
c     to coefficients for Pitzer's equations. The tabulation is
c     controlled by iopr(10):
c
c        0 = Don't print
c        1 = Print a summary of the names of the species present and
c              the number of Pitzer interaction coefficients
c        2 = Print a summary of the names of the species present and
c              the number of Pitzer interaction coefficients
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       iopr   = array of print option switches
c       narn1  = start of aqueous species range
c       narn2  = end of aqueous species range
c       nmux   = array of aqueous species indices defining triplets
c                  of species for which mu data are present
c       nmutmx = the maximum number of triplets of aqueous species
c                  for which Pitzer interaction parameters are
c                  defined; the second dimension of the nmux array
c       nmxi   = range pointer array into the nmxx array:
c                  nmxi(1,na) and nmxi(2,na) are the first and last
c                  values of the second subscript (nmx) of the nmxx
c                  array for entries pertaining to the na-th
c                  aqueous species (ns-th species)
c       nmxx   = pointer array:
c                  nmxx(1,nmx) = the species index of the second
c                  species in the nmu-th triplet, nmxx(2,nmx) is the
c                  species index of the third species, and
c                  nmxx(3,nmx) = nmu
c       nslx   = array of aqueous species indices defining pairs
c                  of species for which S-lambda data are present
c       nsltmx = the maximum number of pairs of aqueous species
c                  for which Pitzer interaction parameters are
c                  defined; the second dimension of the nslx array
c       nsxi   = range pointer array into the nsxx array:
c                  nsxi(1,na) and nsxi(2,na) are the first and last
c                  values of the second subscript (nsx) of the nsxx
c                  array for entries pertaining to the na-th
c                  aqueous species (ns-th species)
c       nsxx   = pointer array:
c                  nsxx(1,nsx) = the species index of the second
c                  species in the nsl-th pair, and nsxx(2,nsx) = nsl
c       uspec  = names of species
c
c     Principal output:
c
c        None
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer natmax,nmutmx,nmxmax,noprmx,nsltmx,nstmax,nsxmax
c
      integer noutpt
c
      integer iopr(noprmx),nmux(3,nmutmx),nmxi(2,natmax),nmxx(3,nmxmax),
     $ nslx(2,nsltmx),nsxi(2,natmax),nsxx(2,nsxmax)
      integer narn1,narn2
c
      character*48 uspec(nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,j3,j4,islt,imut,ixf,ixl,na,nmu,nmx,ns,nsl,nsx,ns1,
     $ ns2,ns3
c
      integer ilnobl
c
c-----------------------------------------------------------------------
c
      if (iopr(10) .le. 0) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The following loop assumes that water is the narn1-th species.
c
      if (iopr(10) .eq. 1) then
c
c       Write a table giving the S-lambda and mu tallies for each
c       aqueous species in the current model.
c
        write (noutpt,1000)
 1000   format(/11x,'--- Pitzer Interaction Coefficient Summary ---',
     $  /4x,'Species',15x,'S-lambda Sets',3x,'Mu Sets',/)
c
c       The following loop assumes that water is the narn1-th species.
c
        do ns = narn1 + 1, narn2
          if (uspec(ns)(1:6).ne.'O2(g) ' .and.
     $      uspec(ns)(1:3).ne.'e- ') then
            na = ns - narn1 + 1
            islt = nsxi(2,na) - nsxi(1,na) + 1
            imut = nmxi(2,na) - nmxi(1,na) + 1
            write (noutpt,1010) uspec(ns),islt,imut
 1010       format(2x,a24,2x,i4,7x,i4)
          endif
        enddo
        write (noutpt,1020)
 1020   format(/1x)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopr(10) .ge. 2) then
c
c       Write a table giving the S-lambda and mu tallies for eachi
c       aqueous species in the current model. Also list the other
c       species involved in the S-lambda and mu combinations.
c
        write (noutpt,1100)
 1100   format(/15x,'--- Pitzer Interaction Coefficient Sets ---')
c
c       The following loop assumes that water is the narn1-th species.
c
        do ns = narn1 + 1, narn2
          if (uspec(ns)(1:6).ne.'O2(g) ' .and.
     $      uspec(ns)(1:3).ne.'e- ') then
            na = ns - narn1 + 1
            j2 = ilnobl(uspec(ns)(1:24))
            write (noutpt,1110) uspec(ns)(1:j2)
 1110       format(//' Coefficients for ',a,':',/)
c
            ixf = nsxi(1,na)
            ixl = nsxi(2,na)
            islt = ixl - ixf + 1
            write (noutpt, 1120) islt
 1120       format(3x,'No. of S-lambda sets= ',i4,':',/)
            if (islt .gt. 0) then
              do nsx = ixf,ixl
                nsl = nsxx(2,nsx)
                ns1 = nslx(1,nsl)
                ns2 = nslx(2,nsl)
                j2 = ilnobl(uspec(ns1)(1:24))
                j3 = ilnobl(uspec(ns2)(1:24))
                write (noutpt,1130) uspec(ns1)(1:j2),uspec(ns2)(1:j3)
 1130           format(5x,a,', ',a)
              enddo
            endif
c
            ixf = nmxi(1,na)
            ixl = nmxi(2,na)
            imut = ixl - ixf + 1
            write (noutpt, 2530) imut
 2530       format(/3x,'No. of mu sets= ',i4,':',/)
            if (imut .gt. 0) then
              do nmx = ixf,ixl
                nmu = nmxx(3,nmx)
                ns1 = nmux(1,nmu)
                ns2 = nmux(2,nmu)
                ns3 = nmux(3,nmu)
                j2 = ilnobl(uspec(ns1)(1:24))
                j3 = ilnobl(uspec(ns2)(1:24))
                j4 = ilnobl(uspec(ns3)(1:24))
                write (noutpt,2540) uspec(ns1)(1:j2),uspec(ns2)(1:j3),
     $          uspec(ns3)(1:j4)
 2540           format(5x,a,', ',a,', ',a)
              enddo
            endif
c
          endif
        enddo
        write (noutpt,1020)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
