      subroutine prtsip(affpd,iopr,ir1,ir2,jpflag,kpsat,kpsst,
     $ noprmx,noutpt,nptmax,sidrph,tolspf,ugroup,uphase)
c
c     This subroutine prints a table of saturation indices and
c     affinities of phases whose indices lie in the range ir1 to ir2.
c     The group corresponding to this range is described by ugroup.
c     The level of printing is controlled by the print control flag
c     iopr(7):
c
c       -1 = Don't print
c        0 = Print for those phases not undersaturated by
c              more than 10 kcal
c        1 = Print for all phases
c
c     This subroutine is called by:
c
c       EQLIB/prtsat.f
c
c-----------------------------------------------------------------------
c
c
c     Principal input:
c
c       affpd  = array of affinities of phases to form (precipitate),
c                  computed from the reactions read from the data file
c       iopr   = array of print control options
c       ir1    = start of the range of pure liquid phases
c       ir2    = end of the range of pure liquid phases
c       jpflag = array of phase status flags
c       kpsat  = number of saturated phases
c       kpsst  = number of supersaturated phases
c       sidrph = array of saturation indices of the various phases,
c                  computed from the reactions read from the data file
c       tolspf = saturation print flag tolerance, used to flag those
c                  phases which are close to saturation
c       ugroup = string describing the group of phases for which
c                  a table is to be printed
c       uphase = array of phase names
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
      integer noprmx,nptmax
c
      integer noutpt
c
      integer iopr(noprmx),jpflag(nptmax)
      integer ir1,ir2,kpsat,kpsst
c
      character*24 uphase(nptmax),ugroup
c
      real*8 affpd(nptmax),sidrph(nptmax)
      real*8 tolspf
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,kount,np
c
      integer ilnobl
c
      real*8 afx
c
c-----------------------------------------------------------------------
c
      j2 = ilnobl(ugroup)
      write (noutpt,1000) ugroup(1:j2)
 1000 format(//11x,'--- Saturation States of ',a,' ---',
     $ //7x,'Phase   ',19x,'Log Q/K',4x,'Affinity, kcal',/)
      kount = 0
c
      do np = ir1,ir2
        if (jpflag(np) .lt. 2) then
          if (iopr(7).gt.0 .or. affpd(np).ge.-10.) then
            kount = kount + 1
            afx = abs(affpd(np))
            if (afx .le. tolspf) then
              kpsat = kpsat + 1
              write (noutpt,1010) uphase(np),sidrph(np),affpd(np)
 1010         format(5x,a24,3x,f10.5,3x,f10.5,5x,'SATD')
            elseif (affpd(np) .gt. tolspf) then
              kpsst = kpsst + 1
              write (noutpt,1020) uphase(np),sidrph(np),affpd(np)
 1020         format(5x,a24,3x,f10.5,3x,f10.5,5x,'SSATD')
            else
              write (noutpt,1030) uphase(np),sidrph(np),affpd(np)
 1030         format(5x,a24,3x,f10.5,3x,f10.5)
            endif
          endif
        endif
      enddo
c
      if (kount .le. 0) write (noutpt,1040)
 1040 format(5x,'None')
c
      if (iopr(7) .eq. 0) write (noutpt,1050)
 1050 format(/8x,'Phases with affinities less than -10 kcal are not',
     $' listed.')
c
      write (noutpt,1060)
 1060 format(1x)
c
      end
