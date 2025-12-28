      subroutine chzrsr(cbsr,elecsr,eps100,jcode,nbasp,nbt,nbtmax,
     $ nbt1mx,noutpt,nrct,nrctmx,nrndex,nsrtmx,nstmax,nttyo,ureac,
     $ uspec,zchar)
c
c     This subroutine checks the reactions for special reactants to
c     ensure that they satisfy charge balance. Reactions made by
c     EQ6/makrsr.f should satisfy this condition, but reactions
c     composed by the user may not. All reactions are checked by this
c     subroutine. All special reactants are presumed to be uncharged.
c
c     This subroutine is called by:
c
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       jcode  = flag array denoting the types of reactants
c       nbasp  = array of indices of species in the basis set
c       nbt    = the number of basis species
c       nbtmax = the maximum number of basis species
c       nbt1mx = the maximum number of basis species plus 1
c       nrctmx = the maximum number of reactants
c       nsrtmx = the maximum number of special reactants
c       ureac  = array of names of reactants
c       uspec  = array of names of species
c       zchar  = array of electrical charges of species
c
c     Principal output:
c
c       cbsr   = array of coefficients for reactions for special
c                  reactants
c       elecsr = array of electrical imbalances of reactions for
c                  special reactants
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nbtmax,nbt1mx,nrctmx,nsrtmx,nstmax
c
      integer noutpt,nttyo
c
      integer jcode(nrctmx),nbasp(nbtmax),nrndex(nrctmx)
      integer nbt,nrct
c
      character*48 uspec(nstmax)
      character*24 ureac(nrctmx)
c
      real*8 cbsr(nbt1mx,nsrtmx),elecsr(nsrtmx),zchar(nstmax)
      real*8 eps100
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer nf
c
      integer j2,nb,nrc,ns,nsr
c
      integer ilnobl
c
      logical qprrsr
c
      character*24 ux24
c
      real*8 acx,aztx,tlersr,ztx
c
c-----------------------------------------------------------------------
c
      tlersr = nbt*eps100
c
      do nrc = 1,nrct
        if (jcode(nrc) .eq. 2) then
          nsr = nrndex(nrc)
          ztx = 0.
          aztx = 0.
          do nb = 1,nbt
            ns = nbasp(nb)
            ztx = ztx + cbsr(nb,nsr)*zchar(ns)
            aztx = aztx + cbsr(nb,nsr)*abs(zchar(ns))
          enddo
          elecsr(nsr) = ztx
c
          qprrsr = .false.
          if (abs(ztx) .gt. tlersr) then
            acx = abs(ztx/aztx)
            if (acx .gt. 1.e-6) then
              j2 = ilnobl(ureac(nrc))
              write (noutpt,1000) ureac(nrc)(1:j2)
              write (nttyo,1000) ureac(nrc)(1:j2)
 1000         format(/' * Warning - (EQ6/chzrsr) The reaction for the',
     $        ' special reactant',/7x,a," isn't well charge balanced.",
     $        /7x,'The reaction is:',/)
              qprrsr = .true.
            elseif (acx .gt. 1.e-10) then
              j2 = ilnobl(ureac(nrc))
              write (noutpt,1010) ureac(nrc)(1:j2)
              write (nttyo,1010) ureac(nrc)(1:j2)
 1010         format(/' * Warning - (EQ6/chzrsr) The reaction for the',
     $        ' special reactant',/7x,a," isn't charge balanced to",
     $        ' high precision.',/7x,'The reaction is:',/)
              qprrsr = .true.
            endif
          endif
c
          if (qprrsr) then
            nf = noutpt
c
c           Calling sequence substitutions:
c             nbasp for nbaspd
c             nbt for nbtd
c
            call prrsr(cbsr,jcode,nbasp,nbt,nbtmax,nbt1mx,
     $      nf,noutpt,nrc,nrctmx,nrndex,nsrtmx,nstmax,nttyo,
     $      ureac,uspec)
c
c           Calling sequence substitutions:
c             nbasp for nbaspd
c             nbt for nbtd
c
            nf = nttyo
            call prrsr(cbsr,jcode,nbasp,nbt,nbtmax,nbt1mx,
     $      nf,noutpt,nrc,nrctmx,nrndex,nsrtmx,nstmax,nttyo,
     $      ureac,uspec)
c
            ux24 = ' '
            write (ux24,'(1pe12.5)') ztx
            j2 = ilnobl(ux24)
            write (noutpt,1050) ux24(1:j2)
            write (nttyo,1050) ux24(1:j2)
 1050       format(/7x,'The imbalance is ',a,'. This will be factored',
     $      ' into the computed',/7x,'charge balance a a compensating',
     $      ' offset. The offset at any point',/7x,'will be',
     $      ' proportional to the extent of reaction. The pH and',
     $      /7x,'redox will be properly calculated, as they are',
     $      ' computed entirely',/7x,'from mass balance relations for',
     $      ' the data file basis species.',/7x,'This version of EQ6',
     $      ' does not compute any of these properties using',
     $      /7x,'the charge balance equation.')
          endif
        endif
      enddo
c
      end
