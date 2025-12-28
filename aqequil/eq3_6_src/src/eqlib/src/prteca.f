      subroutine prteca(cteaq,mrmlra,nct,nctmax,noutpt,ppmwe,
     $ qrho,rho,uelem)
c
c     This subroutine prints a table of the elemental composition of the
c     aqueous solution.
c
c     This subroutine is called by:
c
c       EQ3NR/scripx.f
c       EQ6/scripz.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       cteaq  = array of molalities of the chemical elements
c       ppmwe  = array of ppms by weight of the chemical elements
c       rho    = the density of the aqueous solution
c       uelem  = array of names of chemical elements
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
      integer nctmax
c
      integer noutpt
c
      integer nct
c
      logical qrho
c
      character(len=8) uelem(nctmax)
c
      real(8) cteaq(nctmax),ppmwe(nctmax)
      real(8) mrmlra,rho
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer nc
c
      real(8) molrex,ppmve
c
c-----------------------------------------------------------------------
c
      write (noutpt,1000)
 1000 format(/11x,'--- Elemental Composition of the Aqueous Solution',
     $ ' ---')
c
      if (qrho) then
c
c       Have density data, include mg/L and Molarity.
c
        write (noutpt,1010)
 1010   format(/3x,'Element',8x,'mg/L',7x,'mg/kg.sol',4x,
     $  'Molarity',5x,'Molality',/)
c
        do nc = 1,nct
          ppmve = ppmwe(nc)*rho
          molrex = cteaq(nc)*mrmlra
          write (noutpt,1020) uelem(nc),ppmve,ppmwe(nc),molrex,cteaq(nc)
 1020     format(5x,a8,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5)
        enddo
      else
c
c       Have no density data, so do not include mg/L and Molarity.
c
        write (noutpt,1030)
 1030   format(/3x,'Element',6x,'mg/kg.sol',4x,'Molality',/)
c
        do nc = 1,nct
          write (noutpt,1040) uelem(nc),ppmwe(nc),cteaq(nc)
 1040     format(5x,a8,1x,1pe12.5,1x,1pe12.5)
        enddo
      endif
c
      write (noutpt,1050)
 1050 format(1x)
c
      end
