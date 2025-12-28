      subroutine prtsso(acflg,actlg,affpd,affsd,ixrn1,ixrn2,
     $ jsol,jsomax,ncmpr,noutpt,np,nptmax,nstmax,nxtmax,sidrph,
     $ sidrsp,tolspf,uspec,uphase,uxtype,xbar,xbarlg)
c
c     This subroutine prints tables describing the state and properties
c     of the np-th phase (a solid solution).
c
c     This subroutine is called by:
c
c       EQ3NR/scripx.f
c       EQ6/scripz.f
c
c----------------------------------------------------------------------
c
c     Principal input:
c
c       acflg  = array of log activity coefficients
c       actlg  = array of log activities of species
c       affpd  = array of phase affinities
c       affsd  = array of species affinities
c       sidrph = array of phase saturation indices
c       sidrsp = array of species saturation indices
c       tolspf = saturation print flag tolerance, used to flag those
c                  phases which are close to saturation
c       uphase = array of phase names
c       xbar   = array of mole fractions of species
c       xbarlg = array of log mole fractions of species
c
c     Principal output:
c
c       None
c
c----------------------------------------------------------------------
c
      implicit none
c
c----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer jsomax,nptmax,nstmax,nxtmax
c
      integer noutpt
c
      integer jsol(nxtmax),ncmpr(2,nptmax)
      integer ixrn1,ixrn2,np
c
      character*48 uspec(nstmax)
      character*32 uxtype(jsomax)
      character*24 uphase(nptmax)
c
      real*8 acflg(nstmax),actlg(nstmax),affpd(nptmax),affsd(nstmax),
     $ sidrph(nptmax),sidrsp(nstmax),xbar(nstmax),xbarlg(nstmax)
c
      real*8 tolspf
c
c----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,j3,k,nr1,nr2,ns,nx
c
      integer ilnobl
c
      real*8 aafx,afx
c
c----------------------------------------------------------------------
c
      j2 = ilnobl(uphase(np))
c
      if (np.lt.ixrn1 .or. np.gt. ixrn2) then
        write (noutpt,1000) uphase(np)(1:j2)
        write (noutpt,1000) uphase(np)(1:j2)
 1000   format(/' * Error - (EQLIB/prtsso) Programming error trap:',
     $  /7x,"Can't write an output table for ",a," because",
     $  /7x,"it isn't a solid solution.")
        stop
      endif
c
      nr1 = ncmpr(1,np)
      nr2 = ncmpr(2,np)
c
      nx = np - ixrn1 + 1
      k = jsol(nx)
c
      j3 = ilnobl(uxtype(k))
      write (noutpt,1010) uphase(np)(1:j2),uxtype(k)(1:j3)
 1010 format(/16x,'--- ',a,' ---'//3x,a,/)
c
c     Print mole fractions, activity coefficients, and activities.
c
      write (noutpt,1020)
 1020 format(4x,'Component',20x,'x',11x,'Log x',3x,'Log lambda',2x,
     $ 'Log activity',/)
c
      do ns = nr1,nr2
        if (xbar(ns) .le. 0.) then
          write (noutpt,1030) uspec(ns),xbar(ns)
        else
          write (noutpt,1030) uspec(ns),xbar(ns),xbarlg(ns),
     $    acflg(ns),actlg(ns)
 1030     format(1x,a24,3x,1pe11.4,3(3x,0pf9.4))
        endif
      enddo
      write (noutpt,1040)
 1040 format(1x)
c
c     Print saturation states and affinities.
c
      write (noutpt,1100)
 1100 format(/4x,'Mineral',23x,'Log Q/K',9x,'Aff, kcal',4x,'State',/)
      afx = affpd(np)
      aafx = abs(afx)
      if (aafx .le. tolspf) then
        write (noutpt,1110) uphase(np),sidrph(np),affpd(np)
 1110   format(1x,a24,2x,2(3x,f13.4),3x,'SATD')
      elseif (afx .gt. tolspf) then
        write (noutpt,1120) uphase(np),sidrph(np),affpd(np)
 1120   format(1x,a24,2x,2(3x,f13.4),3x,'SSATD')
      else
        write (noutpt,1130) uphase(np),sidrph(np),affpd(np)
 1130   format(1x,a24,2x,2(3x,f13.4))
      endif
c
      do ns = nr1,nr2
        if (xbar(ns) .gt. 0.) then
          afx = affsd(ns)
          aafx = abs(afx)
          if (aafx .le. tolspf) then
            write (noutpt,1140) uspec(ns),sidrsp(ns),affsd(ns)
 1140       format(3x,a24,2(3x,f13.4),3x,'SATD')
          elseif (afx .gt. tolspf) then
            write (noutpt,1150) uspec(ns),sidrsp(ns),affsd(ns)
 1150       format(3x,a24,2(3x,f13.4),3x,'SSATD')
          else
            write (noutpt,1160) uspec(ns),sidrsp(ns),affsd(ns)
 1160       format(3x,a24,2(3x,f13.4))
          endif
        endif
      enddo
      write (noutpt,1070)
 1070 format(1x)
c
      end
