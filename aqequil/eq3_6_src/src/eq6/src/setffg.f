      subroutine setffg(csts,iindx1,iffg,ipndx1,jpflag,jsflag,kbt,
     $ kdim,kmax,km1,kmt,kx1,kxt,losp,moffg,mtb,mtbaq,nbaspd,nbt,
     $ nbtmax,ncmpr,nffg,nffgmx,noutpt,nphasx,npt,nptmax,nstmax,
     $ nsts,nstsmx,nstsr,nttyo,qloffg,uffg,uspec,uzvec1,zvclg1,zvec1)
c
c     This subroutine puts new fictive fugacity-fixing phases into the
c     matrix if the corresponding masses to add to the equilibrium
c     system (ES) are positive. It keep old such phases in if the
c     corresponding net masses are positive. It deletes any fictive
c     fugacity-fixing phases left over from a previous run for which
c     the corresponding fixed fugacity options no longer apply.
c     To carry out this last function, this subroutine must be called
c     even if nffg is zero. If necessary this subroutine recalculates
c     the mass balance totals.
c
c     This subroutine is called by:
c
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c     Principal output:
c
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer kmax,nbtmax,nffgmx,nptmax,nstmax,nstsmx
c
      integer noutpt,nttyo
c
      integer iindx1(kmax),iffg(nffgmx),ipndx1(kmax),jpflag(nptmax),
     $ jsflag(nstmax),nbaspd(nbtmax),ncmpr(2,nptmax),nphasx(nstmax),
     $ nsts(nstsmx),nstsr(2,nstmax)
c
      integer kbt,kdim,km1,kmt,kx1,kxt,nbt,nffg,npt
c
      logical qloffg
c
      character*48 uspec(nstmax),uzvec1(kmax)
      character*24 uffg(nffgmx)
c
      real*8 csts(nstsmx),losp(nstmax),moffg(nffgmx),mtb(nbtmax),
     $ mtbaq(nbtmax),zvec1(kmax),zvclg1(kmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ier,j2,kcol,n,nb,nerr,np,ns,nse
c
      integer ilnobl
c
      logical qreclc
c
      real*8 cx,mtbo,mtbx,mx,mz
c
      real*8 coefst,texp,tlg
c
c-----------------------------------------------------------------------
c
      qreclc = qloffg
      nerr = 0
c
      do kcol = km1,kxt
        ns = iindx1(kcol)
        losp(ns) = zvclg1(kcol)
      enddo
c
      do n = 1,nffg
        mx = moffg(n)
        ns = iffg(n)
        np = nphasx(ns)
        if (jpflag(np) .ne. -1) then
c
c         The fixed fugacity phase was not already in the matrix.
c
          losp(ns) = -99999.
          if (mx .gt. 0.) then
            qreclc = .true.
            jpflag(np) = -1
            jsflag(ns) = -1
            losp(ns) = tlg(mx)
          elseif (mx .lt. 0.) then
            j2 = ilnobl(uffg(n))
            write (noutpt,1000) uffg(n)(1:j2)
            write (nttyo,1000) uffg(n)(1:j2)
 1000       format(/' * Error - (EQ6/setffg) The new fixed fugacity',
     $      ' phase',
     $      /7x,a," can't have negative moles added.")
            nerr = nerr + 1
          endif
        elseif (mx .ne. 0.) then
c
c         The fixed fugacity phase was already in the matrix.
c
          qreclc = .true.
          mz = texp(losp(ns))
          mz = mz + mx
          if (mz .gt. 0.) then
            losp(ns) = tlg(mz)
          else
            j2 = ilnobl(uffg(n))
            write (noutpt,1010) uffg(n)(1:j2),mz
            write (nttyo,1010) uffg(n)(1:j2),mz
 1010       format(/' * Error - (EQ6/setffg) The adjusted number of',
     $      /7x,'moles of ',a,' is ',g12.5,'. A negative value is not',
     $      /7x,'permitted.')
            nerr = nerr + 1
          endif
        endif
      enddo
c
      if (nerr .gt. 0) stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (.not.qreclc) go to 200
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The following may may combine any of three functions:
c
c       1. Add new fixed fugacity phases to the matrix.
c       2. Delete left over fixed fugacity phases from the matrix.
c       3. Reset the zvec1 and zvclg1 array elements for fixed
c          fugacity phases.
c
      call miidxz(ier,iindx1,ipndx1,jpflag,jsflag,kbt,kdim,
     $ kmax,km1,kmt,kx1,kxt,losp,ncmpr,noutpt,npt,nptmax,nstmax,
     $ nttyo,uspec,uzvec1,zvclg1,zvec1)
c
      if (ier .gt. 0) then
        write (noutpt,1020)
        write (nttyo,1020)
 1020   format(/" * Error - (EQ6/setffg) Can't recover from having",
     $  ' exceeded the maximum ',i4,/7x,'elements of the iindx1',
     $  ' array. Increase the dimensioning parameter kpar.')
        stop
      endif
c
c     Recalculate the mass balance totals. This is done by adding to
c     get the total, not by subtracting from the existing total.
c     This purges the masses of any left over fixed fugacity phases
c     from the mass balance totals and adds any new masses to these
c     totals.
c
      if (nffg .gt. 0) write (noutpt,1030)
 1030 format(//6x,' --- Adding Fictive Fugacity-Fixing Phases to the',
     $ ' ES ---',//13x,'Name',22x,'Moles Added',/)
c
      do n = 1,nffg
        ns = iffg(n)
        write (noutpt,1040) uspec(ns),moffg(n)
 1040   format(11x,a24,3x,1pe12.5)
      enddo
c
      write (noutpt,1050)
 1050 format(//11x,' --- Component Total Numbers of Moles ---',
     $ //35x,'Initial',7x,'Adjusted',/)
c
      do nb = 1,nbt
        nse = nbaspd(nb)
        mtbo = mtb(nb)
        mtbx = mtbaq(nb)
c
        do kcol = km1,kxt
          ns = iindx1(kcol)
          cx = coefst(csts,nsts,nstsmx,nstsr,nb,ns,nstmax)
          mtbx = mtbx + cx*zvec1(kcol)
        enddo
c
        mtb(nb) = mtbx
        write (noutpt,1060) uspec(nse),mtbo,mtbx
 1060   format(5x,a24,3x,1pe12.5,3x,e12.5)
      enddo
      write (noutpt,1070)
 1070 format(1x)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  200 continue
c
      do n = 1,nffg
        moffg(n) = 0.
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
