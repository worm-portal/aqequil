      subroutine prtsat(affpd,iern1,iern2,ilrn1,ilrn2,imrn1,imrn2,
     $ iopr,iopt,ixrn1,ixrn2,jpflag,noutpt,noprmx,noptmx,nptmax,
     $ sidrph,tolspf,uphase)
c
c     This subroutine prints tables of saturation indices and affinities
c     for the various non-aqueous phases. The level of printing is
c     controlled by the print control flag iopr(7):
c
c       -1 = Don't print
c        0 = Print for those phases not undersaturated by
c              more than 10 kcal
c        1 = Print for all phases
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
c       affpd  = array of affinities of phases to form (precipitate),
c                  computed from the reactions read from the data file
c       iern1  = start of the range of generic ion exchangers
c       iern2  = end of the range of generic ion exchangers
c       ilrn1  = start of the range of pure liquid phases
c       ilrn2  = end of the range of pure liquid phases
c       imrn1  = start of the range of pure mineral phases
c       imrn2  = end of the range of pure mineral phases
c       iopt   = array of model option switches
c       iopr   = array of print control options
c       ixrn1  = start of the range of solid solution phases
c       ixrn2  = end of the range of solid solution phases
c       jpflag = array of phase status flags
c       sidrph = array of saturation indices of the various phases,
c                  computed from the reactions read from the data file
c       tolspf = saturation print flag tolerance, used to flag those
c                  phases which are close to saturation
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
      integer noprmx,noptmx,nptmax
c
      integer noutpt
c
      integer iopr(noprmx),iopt(noptmx),jpflag(nptmax)
      integer iern1,iern2,ilrn1,ilrn2,imrn1,imrn2,ixrn1,ixrn2
c
      character*24 uphase(nptmax)
c
      real*8 affpd(nptmax),sidrph(nptmax)
      real*8 tolspf
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ir1,ir2,kpsat,kpsst
c
      character*24 ugroup,upusol,upuliq,usosol,usoliq,ugexch,ux24,uy24
c
c-----------------------------------------------------------------------
c
      data upusol /'Pure Solids             '/,
     $     upuliq /'Pure Liquids            '/,
     $     usosol /'Solid Solutions         '/,
     $     usoliq /'Liquid Solutions        '/,
     $     ugexch /'Generic Ion Exchangers  '/
c
c-----------------------------------------------------------------------
c
c     Note: the following statements don't really do anything except
c     cause the compiler not to complain that usoliq is not used.
c
      ux24 = usoliq
      uy24 = ux24
      usoliq = uy24
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopr(7) .le. -1) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Compute and print phase saturation states. these refer to
c     the dissolution reactions as they are written in the
c     'd' set.
c
c       kpsat = the number of saturated phases
c       kpsst = the number of supersaturated phases
c
      kpsat = 0
      kpsst = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print a table for the pure solids.
c
      ugroup = upusol
      ir1 = imrn1
      ir2 = imrn2
      call prtsip(affpd,iopr,ir1,ir2,jpflag,kpsat,kpsst,
     $ noprmx,noutpt,nptmax,sidrph,tolspf,ugroup,uphase)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print a table for the pure liquids.
c
      ugroup = upuliq
      ir1 = ilrn1
      ir2 = ilrn2
      call prtsip(affpd,iopr,ir1,ir2,jpflag,kpsat,kpsst,
     $ noprmx,noutpt,nptmax,sidrph,tolspf,ugroup,uphase)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iopt(4) .ge. 1) then
c
c       Print a table for the solid solutions.
c
        ugroup = usosol
        ir1 = ixrn1
        ir2 = ixrn2
        call prtsip(affpd,iopr,ir1,ir2,jpflag,kpsat,kpsst,
     $  noprmx,noutpt,nptmax,sidrph,tolspf,ugroup,uphase)
c
c       No table is presently printed for non-aqueous liquid
c       solutions.
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iern2 .ge. iern1) then
c
c       Print a table for the generic ion exchangers.
c
        ugroup = ugexch
        ir1 = iern1
        ir2 = iern2
        call prtsip(affpd,iopr,ir1,ir2,jpflag,kpsat,kpsst,
     $  noprmx,noutpt,nptmax,sidrph,tolspf,ugroup,uphase)
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1150)
 1150 format(/11x,'--- Summary of Saturated and Supersaturated',
     $ ' Phases ---',/)
c
      if (kpsat .le. 0) then
        write (noutpt,1160)
 1160   format(16x,'There are no saturated phases.')
      elseif (kpsat .eq. 1) then
        write (noutpt,1170)
 1170   format(16x,'There is 1 saturated phase.')
      else
        write (noutpt,1180) kpsat
 1180   format(16x,'There are ',i4,' saturated phases.')
      endif
c
      if (kpsst .le. 0) then
        write (noutpt,1200)
 1200   format(16x,'There are no supersaturated phases.')
      elseif (kpsst .eq. 1) then
        write (noutpt,1210)
 1210   format(16x,'There is 1 supersaturated phase.')
      else
        write (noutpt,1220) kpsst
 1220   format(16x,'There are ',i4,' supersaturated phases.')
      endif
      write (noutpt,1230)
 1230 format(1x)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
