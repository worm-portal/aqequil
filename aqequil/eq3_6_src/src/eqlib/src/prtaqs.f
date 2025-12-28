      subroutine prtaqs(acflg,actlg,conc,conclg,iopr,jcsort,narn1,
     $ narn2,noprmx,noutpt,nstmax,uspec)
c
c     This subroutine prints a table of the concentrations, activities,
c     and activity coefficients of the aqueous solute species. The
c     species are listed in decreasing order of concentration. The
c     level of printing is controlled by the print control flag
c     iopr(4):
c
c       -3  = Omit species with molalities < 1.e-8
c       -2 =  Omit species with molalities < 1.e-12
c       -1 =  Omit species with molalities < 1.e-20
c        0 =  Omit species with molalities < 1.e-100
c        1 =  Include all species
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
c       acflg  = array of species activity coefficients
c       act    = array of species activities
c       actlg  = array of species log activities
c       conc   = array of species concentrations
c       conclg = array of species log concentrations
c       iopr   = array of print control options
c       jcsort = array of species indices, in order of increasing
c                  concentration, but with sorting restricted to within
c                  phase ranges
c       narn1  = start of the range of aqueous species
c       narn2  = end of the range of aqueous species
c       uspec  = array of names of species
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
      integer noprmx,nstmax
c
      integer noutpt
c
      integer iopr(noprmx),jcsort(nstmax)
      integer narn1,narn2
c
      character*48 uspec(nstmax)
c
      real*8 acflg(nstmax),actlg(nstmax),conc(nstmax),conclg(nstmax)
c
      real*8 texp
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ipr,n,ns,nss
c
      real*8 ltrunc(-3:0)
      real*8 mx
c
c-----------------------------------------------------------------------
c
      data (ltrunc(n), n = -3,0)/-8.,-12.,-20.,-100./
c
c-----------------------------------------------------------------------
c
      ipr = iopr(4)
c
      write (noutpt,1000)
 1000 format(//16x,'--- Distribution of Aqueous Solute Species ---',/)
c
      write (noutpt,1010)
 1010 format(4x,'Species',18x,'Molality',4x,'Log Molality',
     $ 3x,'Log Gamma',2x,'Log Activity',/)
c
      if (iopr(4) .ge. 1) then
        do nss = narn1 + 1,narn2
          n = narn2 - nss + narn1
          ns = jcsort(n)
          write (noutpt,1020) uspec(ns),conc(ns),conclg(ns),
     $    acflg(ns),actlg(ns)
 1020     format(1x,a24,2x,1pe11.4,2x,0pf11.4,2x,f11.4,2x,f11.4)
        enddo
      else
        do nss = narn1 + 1,narn2
          n = narn2 - nss + narn1
          ns = jcsort(n)
          if (conclg(ns) .lt. ltrunc(ipr)) go to 110
          write (noutpt,1020) uspec(ns),conc(ns),conclg(ns),
     $    acflg(ns),actlg(ns)
        enddo
        go to 999
  110   mx = texp(ltrunc(ipr))
        write (noutpt,1030) mx
 1030   format(/4x,'Species with molalities less than ',1pe10.3,' are',
     $  ' not listed.')
      endif
      write (noutpt,1040)
 1040 format(1x)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
