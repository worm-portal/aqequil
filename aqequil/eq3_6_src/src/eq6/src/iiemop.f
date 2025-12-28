      subroutine iiemop(iemop,iemos,iindx1,ipndx1,jsflag,kdim,kmax,
     $ ncmpe,ncmpr,noutpt,npet,npetmx,npt,nptmax,nset,nsetmx,nstmax,
     $ nttyo,uaqsln,uspec,uphase)
c
c     This subroutine initializes the indexing used for tracking the
c     numbers of moles of phases present in the Equilibrium System (ES).
c     This indexing must be re-initialized whenever a change occurs
c     in the ES phase assemblage. The iemop array contains the indices
c     of the phases in the ES. The iemos array contains the indicies
c     of the active species in these phases, except that only the
c     species H2O(l) is included in the case of the aqueous solution
c     phase. solution. The ncmpe array is a species range pointer
c     array for the phases analogous to ncmpr.
c
c     The iemop and emop arrays respectively contain the indices and
c     numbers of moles of the phases currently in the ES. The iemos
c     and emos array are the analogs for the species of these phases.
c     However, only the species H2O(l) is represented for the aqueous
c     solution.  The ncmpe array is a species range pointer array for
c     the phases analogous to ncmpr.
c
c     The fdpe0 and fdse0 arrays respectively contain the finite
c     differences for the data in the emop and emos arrays. The demop
c     and demos arrays contain the corresponding derivatives with
c     respect to the reaction progress variable, xi.
c
c     This subroutine is called by:
c
c       EQ6/dumpdp.f
c       EQ6/path.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
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
      integer kmax,npetmx,nptmax,nsetmx,nstmax
c
      integer noutpt,nttyo
c
      integer iemop(npetmx),iemos(nsetmx),iindx1(kmax),ipndx1(kmax),
     $ jsflag(nstmax),ncmpe(2,npetmx),ncmpr(2,nptmax)
c
      integer kdim,npet,npt,nset
c
      character*48 uspec(nstmax)
      character*24 uphase(nptmax)
c
      character*24 uaqsln
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jlen,j2,kcol,np,npe,nplast,nr1,nr2,ns,nse
c
      integer ilnobl
c
      character*56 uspn56
      character*8 ux8
c
c-----------------------------------------------------------------------
c
c     Loop over all phases present in the ES.
c
      npe = 0
      nse = 0
      nplast = 0
      do kcol = 1,kdim
        np = ipndx1(kcol)
        if (np .ne. nplast) then
          npe = npe + 1
          iemop(npe) = np
          nplast = np
          nr1 = ncmpr(1,np)
          nr2 = ncmpr(2,np)
          if (uphase(np)(1:24) .eq. uaqsln(1:24)) nr2 = nr1
          ncmpe(1,npe) = nse + 1
          do ns = nr1,nr2
            if (jsflag(ns) .le. 0) then
              nse = nse + 1
c
              if (nse .gt. nsetmx) then
c
c               Calling sequence substitutions:
c                 uspec(ns) for unam48
c
                call fmspnm(jlen,uspec(ns),uspn56)
                write (ux8,'(i5)') nsetmx
                call lejust(ux8)
                j2 = ilnobl(ux8)
                write (noutpt,1000) ux8(1:j2),uspn56(1:jlen)
                write (nttyo,1000) ux8(1:j2),uspn56(1:jlen)
 1000           format(/' * Error - (EQ6/iiemop) Have too many species',
     $          ' to track by finite differences.',/7x,'Exceeded the',
     $          ' dimensioned limit of ',a,' upon trying to set up',
     $          /7x,'tracking for ',a,'. Increase the dimensioning',
     $          /7x,'parameter nsetpa.')
              endif
c
              iemos(nse) = ns
            endif
          enddo
          ncmpe(2,npe) = nse
        endif
      enddo
      npet = npe
      nset = nse
c
      end
