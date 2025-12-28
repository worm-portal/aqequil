      subroutine bsplkp(axlks,narxmx,nbasp,nbt,nbtmax,ndrs,ndrsmx,
     $ ndrsr,nstmax,ntprmx)
c
c     This subroutine looks at each active auxiliary basis species.
c     It changes the corresponding log K polynomial coefficients so
c     that log K is fixed at a value of -9999999. if any other species
c     in the corresponding dissociation reaction is not in the model.
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
c       nbasp  = array of indices of species in the active basis set
c
c     Principal output:
c
c       axlks  = array of coefficients for computing equilbrium
c                  constants as a function of temperature
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer narxmx,nbtmax,ndrsmx,nstmax,ntprmx
c
      integer nbasp(nbtmax),ndrs(ndrsmx),ndrsr(2,nstmax)
      integer nbt
c
      real*8 axlks(narxmx,ntprmx,nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,j,n,nb,ns,nr1,nr2,nt1,nse
c
c-----------------------------------------------------------------------
c
      do nb = 1,nbt
        ns = nbasp(nb)
        nr1 = ndrsr(1,ns)
        nr2 = ndrsr(2,ns)
        nt1 = nr2 - nr1 + 1
        if (nt1 .ge. 2) then
          do n = nr1,nr2
            nse = ndrs(n)
            if (nse .eq. 0) then
              do j = 1,ntprmx
                do i = 2,narxmx
                  axlks(i,j,ns) = 0.
                enddo
                axlks(1,j,ns) = -9999999.
              enddo
              go to 130
            endif
          enddo
        endif
  130   continue
      enddo
c
      end
