      subroutine prtacr(actlg,iopr,jsflag,nbaspd,nbt,nbtmax,nelect,
     $ nhydr,noprmx,no2gaq,noutpt,nstmax,uspec,zchar)
c
c     This subroutine prints a table of cation/H+ activity ratios,
c     anion-H+ activity products, and neutral species activities.
c     Only aqueous basis species are involved. The level of printing
c     is controlled by the print control flag iopr(5):
c
c        0 = Don't print
c        1 = Print cation/H+ activity ratios only
c        2 = Print cation/H+ activity ratios and anion-H+ activity
c              products
c        3 = Print cation/H+ activity ratios, anion-H+ activity
c              products, and neutral species activities
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
c       actlg  = array of log activities of species
c       iopr   = array of print control options
c       jsflag = array of species status flags
c       nbt    = the number of species in the basis set
c       nelect = index of the fictive species aqueous e-
c       nhydr  = index of the species aqueous H+
c       no2gaq = index of the fictive species aqueous O2(g)
c       uspec  = array of species names
c       zchar  = array of species electrical charge numbers
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
      integer nbtmax,noprmx,nstmax
c
      integer noutpt
c
      integer iopr(noprmx),jsflag(nstmax),nbaspd(nbtmax)
      integer nbt,nelect,nhydr,no2gaq
c
      character*48 uspec(nstmax)
c
      real*8 actlg(nstmax),zchar(nstmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer iz,j2,nb,ns
c
      integer ilnobl
c
      real*8 actlh,actrat,zx
c
c-----------------------------------------------------------------------
c
      if (iopr(5) .le. 0) go to 999
      if (iopr(5) .eq. 1) then
        write (noutpt,1000)
 1000   format(/16x,'--- Cation/H+ Activity Ratios ---',/)
      elseif (iopr(5) .eq. 2) then
        write (noutpt,1010)
 1010   format(/16x,'--- Ion-H+ Activity Ratios ---',/)
      elseif (iopr(5) .eq. 3) then
        write (noutpt,1020)
 1020   format(/6x,'--- Ion-H+ Activity Ratios and',
     $  ' Neutral Species Activities ---',/)
      endif
c
      actlh = actlg(nhydr)
c
      do nb = 1,nbt
        ns = nbaspd(nb)
        if (ns .eq. nhydr) go to 100
        if (ns .eq. no2gaq) go to 100
        if (ns .eq. nelect) go to 100
        if (jsflag(ns) .ge. 2) go to 100
        zx = zchar(ns)
        actrat = actlg(ns)
        if (zx .eq. 0.) then
          if (iopr(5) .ge. 3) then
            j2 = ilnobl(uspec(ns)(1:16))
            write (noutpt,1040) uspec(ns)(1:j2),actrat
 1040       format(3x,'Log ( a(',a,') )',t44,'= ',f10.5)
          endif
        else
          actrat = actrat - zx*actlh
          j2 = ilnobl(uspec(ns)(1:16))
          iz = nint(abs(zx))
          if (zx .lt. 0.) then
            if (iopr(5) .ge. 2) then
              write (noutpt,1050) uspec(ns)(1:j2),iz,actrat
 1050         format(3x,'Log ( a(',a,') x a(H+)xx',i2,' )',t44,'= ',
     $        f10.5)
            endif
          else
            write (noutpt,1060) uspec(ns)(1:j2),iz,actrat
 1060       format(3x,'Log ( a(',a,') / a(H+)xx',i2,' )',t44,'= ',
     $      f10.5)
          endif
        endif
  100   continue
      enddo
      write (noutpt,1070)
 1070 format(1x)
c
  999 continue
      end
