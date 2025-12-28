      subroutine prtfug(jgsort,fugac,fugalg,jsflag,ngrn1,ngt,ngtmax,
     $ noutpt,nstmax,uspec)
c
c     This subroutine prints a table giving the equilibrium fugacities
c     of gas species.
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
c       fugac  = gas fugacity array
c       fugalg = gas log fugacity array
c       jsflag = species status flag array
c       ngrn1  = start of gas species range
c       ngt    = number of gas species
c       uspec  = species name array
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
      integer ngtmax,nstmax
c
      integer noutpt
c
      integer jgsort(ngtmax),jsflag(nstmax)
      integer ngrn1,ngt
c
      character*48 uspec(nstmax)
c
      real*8 fugac(ngtmax),fugalg(ngtmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ncount,ng,ngg,ns
c
c-----------------------------------------------------------------------
c
      write(noutpt,1000)
 1000 format(/21x,'--- Fugacities ---',
     $ //5x,'Gas ',20x,'Log Fugacity',4x,'Fugacity',/)
c
      ncount = 0
      do ngg= ngt,1,-1
        ng = jgsort(ngg)
        ns = ngrn1 + ng - 1
        if (jsflag(ns) .lt. 2) then
          ncount = ncount + 1
          write (noutpt,1020) uspec(ns),fugalg(ng),fugac(ng)
 1020     format(3x,a24,3x,f10.5,3x,1pe12.5)
        endif
      enddo
      if (ncount .le. 0) write (noutpt,1030)
 1030 format(3x,'None')
c
      write (noutpt,1040)
 1040 format(1x)
c
      end
