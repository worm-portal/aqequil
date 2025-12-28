      subroutine echolk(axlks,cdrs,ilevel,jsflag,narxmx,ndrs,ndrsmx,
     $ ndrsr,nf,nst,ntprmx,nstmax,press,tempc,uspec,xlks)
c
c     This subroutine prints the species and reactions that are active
c     in the current problem, along with the log K values that
c     correspond to the reactions. Optionally, the coefficients
c     of the interpolating polynomials may also be printed.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/echoz.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       ilevel = print level switch
c                  =  0  no print
c                  =  1  print species and reactions only
c                  =  2  also print equilibrium constants
c                  =  3  also print the coefficients of the
c                        interpolating polynomials
c       nf     = unit number of the file to write on
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
      integer narxmx,ndrsmx,ntprmx,nstmax
c
      integer jsflag(nstmax),ndrs(ndrsmx),ndrsr(2,nstmax)
      integer ilevel,nf,nst
c
      character*48 uspec(nstmax)
c
      real*8 axlks(narxmx,ntprmx,nstmax),cdrs(ndrsmx),xlks(nstmax)
      real*8 press,tempc
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,j,ns
c
c-----------------------------------------------------------------------
c
      if (ilevel .le. 0) go to 999
c
      write (nf,2005)
 2005 format(//' --- Listing of Species and Reactions ---',/)
      if (ilevel .ge. 2) write (nf,2040) tempc,press
 2040 format(7x,'temperature= ',f10.3,'C',
     $ /7x,'pressure= ',g13.6,' bars',/)
c
      do ns = 1,nst
      if (jsflag(ns) .le. 0) then
        write (nf,2050)
 2050   format(' --------------------------------------------------')
c
        call prreac(cdrs,ndrs,ndrsmx,ndrsr,nf,ns,nstmax,uspec)
        if (ilevel .ge. 2) write (nf,2060) xlks(ns)
 2060   format(/10x,'log K= ',f12.4,/)
        if (ilevel .ge. 3) then
          do j = 1,ntprmx
          write (nf,2070) j,(axlks(i,j,ns), i = 1,narxmx)
 2070     format(/3x,'Coefficients for range ',i2,':',/3x,3(2x,g13.6),
     $    /3x,3(2x,g13.6),/)
          enddo
        endif
      endif
      enddo
      write (nf,2050)
c
  999 continue
      end
