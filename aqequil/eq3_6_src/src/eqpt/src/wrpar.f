      subroutine wrpar(aamatr,adh,adhh,adhv,aphi,apr,avgrid,bdh,bdhh,
     $ bdhv,bdot,bdoth,bdotv,cco2,cof,dadhh,dadhv,dbdhh,dbdhv,dbdth,
     $ dbdtv,dhfe,dvfe,eps100,gmmatr,ipch,ipchmx,ipcv,ipcvmx,ipivot,
     $ narxmx,narxt,ndata1,ndat1f,noutpt,ntprmx,ntprt,nttyo,presg,
     $ prehw,tdamax,tdamin,tempc,tempcs,tmpcmx,uakey,xhfe,xlke,xvfe,
     $ xvec,yvec)
c
c     This subroutine takes the data read from the DATA0 file by the
c     EQPT/rdpar.f, processes it by converting data on a temperature
c     grid to the equivalent set of coefficients of interpolating
c     polynomials, and writes this data on the DATA1 and DATA1F files.
c
c     This subroutine is called by:
c
c       EQPT/eqpt.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       adh    = array of A(gamma,10) values on the standard temperature
c                  grid
c       aphi   = array of A(phi) values on the standard temperature grid
c       bdh    = array of B(gamma) values on the standard temperature
c                  grid
c       bdot   = array of B-dot values on the standard temperature
c                  grid
c       cco2   = array of coefficients for the Drummond (1981)
c                  equation
c       ndata1 = unit number of the DATA1 file
c       ndat1f = unit number of the DATA1F file
c       presg  = array of pressures on the standard temperature grid
c       tempc  = array of temperatures defining the standard
c                  temperature grid
c       uakey  = string specifying the type of data file ("SEDH" or
c                  "Pitzer") being processed
c       xlke   = array of log K values for the special "Eh" reaction
c                  on the standard temperature grid
c
c     Principal output:
c
c       apr    = array of polynomial coefficients (combined, for all
c                  temperature ranges)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ipchmx,ipcvmx,narxmx,ntprmx
c
      integer ipivot(narxmx),narxt(ntprmx)
c
      integer ndata1,ndat1f,noutpt,nttyo
c
      integer ipch,ipcv,ntprt
c
      character(len=8) uakey
c
      real(8) tdamax,tdamin
      real(8) cco2(5)
c
      real(8) adh(narxmx,ntprmx),adhh(narxmx,ntprmx),
     $ adhv(narxmx,ntprmx),aphi(narxmx,ntprmx),bdh(narxmx,ntprmx),
     $ bdhh(narxmx,ntprmx),bdhv(narxmx,ntprmx),bdot(narxmx,ntprmx),
     $ bdoth(narxmx,ntprmx),bdotv(narxmx,ntprmx),
     $ dadhh(narxmx,ntprmx,ipchmx),dadhv(narxmx,ntprmx,ipcvmx),
     $ dbdhh(narxmx,ntprmx,ipchmx),dbdhv(narxmx,ntprmx,ipcvmx),
     $ dbdth(narxmx,ntprmx,ipchmx),dbdtv(narxmx,ntprmx,ipcvmx),
     $ dhfe(narxmx,ntprmx,ipchmx),dvfe(narxmx,ntprmx,ipcvmx),
     $ prehw(narxmx,ntprmx),presg(narxmx,ntprmx),xhfe(narxmx,ntprmx),
     $ xlke(narxmx,ntprmx),xvfe(narxmx,ntprmx)
c
      real(8) apr(narxmx,ntprmx),avgrid(narxmx,ntprmx)
      real(8) tempc(narxmx,ntprmx),tempcs(narxmx,ntprmx),tmpcmx(ntprmx)
      real(8) aamatr(narxmx,narxmx),gmmatr(narxmx,narxmx)
      real(8) cof(narxmx),xvec(narxmx),yvec(narxmx)
c
      real(8) eps100
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ipc,j2,j3,n,nt,ntpr
c
      integer ilnobl
c
      character(len=10) ux10
      character(len=72) uterm,utermc
      character(len=56) ux56
      character(len=24) ux24,ux24a,ux24b
c
      real(8) adtxmn,adtx11,p01_8,txmn,tx11
c
c-----------------------------------------------------------------------
c
      uterm(1:48) = '+-----------------------------------------------'
      uterm(49:72) = '------------------------'
      utermc = uterm
      utermc(1:1) = '*'
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the temperature limits on the output files.
c
      ux56 = 'Data file maximum and minimum temperatures (C)'
      write (ndata1) ux56
      j2 = ilnobl(ux56)
      write (ndat1f,1100) ux56(1:j2)
 1100 format(a)
      write (ndata1) tdamin,tdamax
      write (ndat1f,1110) tdamin,tdamax
 1110 format(2f10.3)
c
      write (noutpt,1120) tdamin
      write (nttyo,1120) tdamin
 1120 format(//' The minimum temperature is ',f10.3,'C')
      write (noutpt,1130) tdamax
      write (nttyo,1130) tdamax
 1130 format(' The maximum temperature is ',f10.3,'C')
c
      ux56 = 'Maximum temperature (C) by range'
      write (ndata1) ux56
      j2 = ilnobl(ux56)
      write (ndat1f,1140) ux56(1:j2)
 1140 format(a)
      do ntpr = 1,ntprt
        n = narxt(ntpr)
        write (ndata1) tempc(n,ntpr)
        write (ndat1f,1150) tempc(n,ntpr)
 1150   format(f10.3)
      enddo
c
      write (ndat1f,1010) utermc(1:72)
c
      write (noutpt,1160)
      write (nttyo,1160)
 1160 format(/' The maximum temperatures (C) by range are:')
      do ntpr = 1,ntprt
        n = narxt(ntpr)
        write (noutpt,1170) tempc(n,ntpr)
        write (nttyo,1170) tempc(n,ntpr)
 1170   format(3x,f10.3)
      enddo
c
      write (noutpt,1000)
 1000 format(/)
c
c
      tx11 = tempc(1,1)
      txmn = tdamin
c
c     Note: on some data files, 0.01C is used as an approximation
c     to 0C in order to avoid some equation-of-state difficulties.
c     In such cases, consider 0.01C to be the same as 0C for the
c     purpose of the following test.
c
      ux10 = '0.01      '
      read (ux10,'(f10.3)') p01_8
      adtx11 = abs(tx11 - p01_8)
      adtxmn = abs(txmn - p01_8)
      if (adtx11 .le. 1.e-6) tx11 = 0.
      if (adtxmn .le. 1.e-6) txmn = 0.
c
      if ((tx11 - tdamin) .gt. 0.001) then
        write (ux24a,"(f10.3)") tdamin
        call lejust(ux24a)
        j2 = ilnobl(ux24a)
        write (ux24b,"(f10.3)") tempc(1,1)
        call lejust(ux24b)
        j3 = ilnobl(ux24b)
        write (noutpt,1200) ux24a(1:j2),ux24b(1:j3)
        write (nttyo,1200) ux24a(1:j2),ux24b(1:j3)
 1200   format(/' * Warning - (EQPT/rdpar) The minimum temperature for',
     $  ' this data file',/7x,'is less than the lowest temperature in',
     $  ' the first temperature',/7x,'range. If EQ3NR and EQ6 are',
     $  ' required to make any calculations',/7x,'between ',a,'C',
     $  ' and ',a,'C, they will do so by extrapolating',/7x,'the',
     $  ' data from the first range.')
      endif
c
      n = narxt(ntprt)
      if ((tdamax - tempc(n,ntprt)) .gt. 0.001) then
        write (ux24a,"(f10.3)") tempc(n,ntprt)
        call lejust(ux24a)
        j2 = ilnobl(ux24a)
        write (ux24b,"(f10.3)") tdamax
        call lejust(ux24b)
        j3 = ilnobl(ux24b)
        write (noutpt,1210) ux24a(1:j2),ux24b(1:j3)
        write (nttyo,1210) ux24a(1:j2),ux24b(1:j3)
 1210   format(/' * Warning - (EQPT/rdpar) The maximum temperature for',
     $  ' this data file',/7x,'is greater than the highest temperature',
     $  ' in the last temperature',/7x,'range. If EQ3NR and EQ6 are',
     $  ' required to make any calculations',/7x,'between ',a,'C',
     $  ' and ',a,'C, they will do so by extrapolating',/7x,'the',
     $  ' data from the last range.')
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Process and write the data for the standard pressure grid.
c
c     Calling sequence substitutions:
c       presg for avgrid
c
      call intrp(aamatr,apr,presg,cof,eps100,gmmatr,ipivot,
     $ narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $ xvec,yvec)
c
      ux24 = 'presg'
      j2 = ilnobl(ux24)
      write (ndata1) ux24
      write (ndat1f,1010) ux24(1:j2)
 1010 format(a)
      write (noutpt,1020) ux24(1:j2)
 1020 format(7x,a)
c
      do ntpr = 1,ntprt
        nt = narxt(ntpr)
        write (ndata1) (apr(i,ntpr), i = 1,nt)
        write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
 1030   format( 5(1pe16.9) )
      enddo
c
      write (ndat1f,1010) utermc(1:72)
c
      if (ipcv .ge. 0) then
c
c       Process and write the data for the half-width of the standard
c       pressure envelope.
c
c       Calling sequence substitutions:
c         prehw for avgrid
c
        call intrp(aamatr,apr,prehw,cof,eps100,gmmatr,ipivot,
     $  narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $  xvec,yvec)
c
        ux24 = 'prehw'
        j2 = ilnobl(ux24)
        write (ndata1) ux24
        write (ndat1f,1010) ux24(1:j2)
        write (noutpt,1020) ux24(1:j2)
c
        do ntpr = 1,ntprt
          nt = narxt(ntpr)
          write (ndata1) (apr(i,ntpr), i = 1,nt)
          write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
        enddo
c
        write (ndat1f,1010) utermc(1:72)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Process and write the next parameters only if processing a
c     "SEDH" type data file.
c
      if (uakey(1:8) .eq. 'SEDH    ') then
c
c       Process and write the A(gamma,10) data.
c
c       Calling sequence substitutions:
c         adh for avgrid
c
        call intrp(aamatr,apr,adh,cof,eps100,gmmatr,ipivot,
     $  narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $  xvec,yvec)
c
        ux24 = 'adh'
        j2 = ilnobl(ux24)
        write (ndata1) ux24
        write (ndat1f,1010) ux24(1:j2)
        write (noutpt,1020) ux24(1:j2)
c
        do ntpr = 1,ntprt
          nt = narxt(ntpr)
          write (ndata1) (apr(i,ntpr), i = 1,nt)
          write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
        enddo
c
        write (ndat1f,1010) utermc(1:72)
c
        if (ipch .ge. 0) then
c
c         Process and write the A(H) data.
c
c         Calling sequence substitutions:
c           adhh for avgrid
c
          call intrp(aamatr,apr,adhh,cof,eps100,gmmatr,ipivot,
     $    narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $    xvec,yvec)
c
          ux24 = 'adhh'
          j2 = ilnobl(ux24)
          write (ndata1) ux24
          write (ndat1f,1010) ux24(1:j2)
          write (noutpt,1020) ux24(1:j2)
c
          do ntpr = 1,ntprt
            nt = narxt(ntpr)
            write (ndata1) (apr(i,ntpr), i = 1,nt)
            write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
          enddo
c
          write (ndat1f,1010) utermc(1:72)
c
c         Process and write the dnA(H)/dPn data.
c
          do ipc = 1,ipch
c
            do ntpr = 1,ntprt
              do n = 1,narxt(ntpr)
                avgrid(n,ntpr) = dadhh(n,ntpr,ipc)
              enddo
            enddo
c
            call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,
     $      narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $      xvec,yvec)
c
            ux24 = 'dadhh( )'
            write (ux24(7:7),'(i1)') ipc
            j2 = ilnobl(ux24)
            write (ndata1) ux24
            write (ndat1f,1010) ux24(1:j2)
            write (noutpt,1020) ux24(1:j2)
c
            do ntpr = 1,ntprt
              nt = narxt(ntpr)
              write (ndata1) (apr(i,ntpr), i = 1,nt)
              write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
            enddo
c
            write (ndat1f,1010) utermc(1:72)
          enddo
        endif
c
        if (ipcv .ge. 0) then
c
c         Process and write the A(V) data.
c
c         Calling sequence substitutions:
c           adhv for avgrid
c
          call intrp(aamatr,apr,adhv,cof,eps100,gmmatr,ipivot,
     $    narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $    xvec,yvec)
c
          ux24 = 'adhv'
          j2 = ilnobl(ux24)
          write (ndata1) ux24
          write (ndat1f,1010) ux24(1:j2)
          write (noutpt,1020) ux24(1:j2)
c
          do ntpr = 1,ntprt
            nt = narxt(ntpr)
            write (ndata1) (apr(i,ntpr), i = 1,nt)
            write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
          enddo
c
          write (ndat1f,1010) utermc(1:72)
c
c         Process and write the dnA(V)/dPn data.
c
          do ipc = 1,ipcv
c
            do ntpr = 1,ntprt
              do n = 1,narxt(ntpr)
                avgrid(n,ntpr) = dadhv(n,ntpr,ipc)
              enddo
            enddo
c
            call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,
     $      narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $      xvec,yvec)
c
            ux24 = 'dadhv( )'
            write (ux24(7:7),'(i1)') ipc
            j2 = ilnobl(ux24)
            write (ndata1) ux24
            write (ndat1f,1010) ux24(1:j2)
            write (noutpt,1020) ux24(1:j2)
c
            do ntpr = 1,ntprt
              nt = narxt(ntpr)
              write (ndata1) (apr(i,ntpr), i = 1,nt)
              write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
            enddo
c
            write (ndat1f,1010) utermc(1:72)
          enddo
        endif
c
c       Process and write the B(gamma) data.
c
c       Calling sequence substitutions:
c         bdh for avgrid
c
        call intrp(aamatr,apr,bdh,cof,eps100,gmmatr,ipivot,
     $  narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $  xvec,yvec)
c
        ux24 = 'bdh'
        j2 = ilnobl(ux24)
        write (ndata1) ux24
        write (ndat1f,1010) ux24(1:j2)
        write (noutpt,1020) ux24(1:j2)
c
        do ntpr = 1,ntprt
          nt = narxt(ntpr)
          write (ndata1) (apr(i,ntpr), i = 1,nt)
          write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
        enddo
c
        write (ndat1f,1010) utermc(1:72)
c
        if (ipch .ge. 0) then
c
c         Process and write the B(H) data.
c
c         Calling sequence substitutions:
c           bdhh for avgrid
c
          call intrp(aamatr,apr,bdhh,cof,eps100,gmmatr,ipivot,
     $    narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $    xvec,yvec)
c
          ux24 = 'bdhh'
          j2 = ilnobl(ux24)
          write (ndata1) ux24
          write (ndat1f,1010) ux24(1:j2)
          write (noutpt,1020) ux24(1:j2)
c
          do ntpr = 1,ntprt
            nt = narxt(ntpr)
            write (ndata1) (apr(i,ntpr), i = 1,nt)
            write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
          enddo
c
          write (ndat1f,1010) utermc(1:72)
c
c         Process and write the dnB(H)/dPn data.
c
          do ipc = 1,ipch
c
            do ntpr = 1,ntprt
              do n = 1,narxt(ntpr)
                avgrid(n,ntpr) = dbdhh(n,ntpr,ipc)
              enddo
            enddo
c
            call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,
     $      narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $      xvec,yvec)
c
            ux24 = 'dbdhh( )'
            write (ux24(7:7),'(i1)') ipc
            j2 = ilnobl(ux24)
            write (ndata1) ux24
            write (ndat1f,1010) ux24(1:j2)
            write (noutpt,1020) ux24(1:j2)
c
            do ntpr = 1,ntprt
              nt = narxt(ntpr)
              write (ndata1) (apr(i,ntpr), i = 1,nt)
              write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
            enddo
c
            write (ndat1f,1010) utermc(1:72)
          enddo
        endif
c
        if (ipcv .ge. 0) then
c
c         Process and write the B(V) data.
c
c         Calling sequence substitutions:
c           bdhv for avgrid
c
          call intrp(aamatr,apr,bdhv,cof,eps100,gmmatr,ipivot,
     $    narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $    xvec,yvec)
c
          ux24 = 'bdhv'
          j2 = ilnobl(ux24)
          write (ndata1) ux24
          write (ndat1f,1010) ux24(1:j2)
          write (noutpt,1020) ux24(1:j2)
c
          do ntpr = 1,ntprt
            nt = narxt(ntpr)
            write (ndata1) (apr(i,ntpr), i = 1,nt)
            write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
          enddo
c
          write (ndat1f,1010) utermc(1:72)
c
c         Process and write the dnB(V)/dPn data.
c
          do ipc = 1,ipcv
c
            do ntpr = 1,ntprt
              do n = 1,narxt(ntpr)
                avgrid(n,ntpr) = dbdhv(n,ntpr,ipc)
              enddo
            enddo
c
            call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,
     $      narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $      xvec,yvec)
c
            ux24 = 'dbdhv( )'
            write (ux24(7:7),'(i1)') ipc
            j2 = ilnobl(ux24)
            write (ndata1) ux24
            write (ndat1f,1010) ux24(1:j2)
            write (noutpt,1020) ux24(1:j2)
c
            do ntpr = 1,ntprt
              nt = narxt(ntpr)
              write (ndata1) (apr(i,ntpr), i = 1,nt)
              write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
            enddo
c
            write (ndat1f,1010) utermc(1:72)
          enddo
        endif
c
c       Process and write the B-dot data.
c       interpolate and write the polynomial coefficients for bdot
c
c       Calling sequence substitutions:
c         bdot for avgrid
c
        call intrp(aamatr,apr,bdot,cof,eps100,gmmatr,ipivot,
     $  narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $  xvec,yvec)
c
        ux24 = 'bdot'
        j2 = ilnobl(ux24)
        write (ndata1) ux24
        write (ndat1f,1010) ux24(1:j2)
        write (noutpt,1020) ux24(1:j2)
c
        do ntpr = 1,ntprt
          nt = narxt(ntpr)
          write (ndata1) (apr(i,ntpr), i = 1,nt)
          write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
        enddo
c
        write (ndat1f,1010) utermc(1:72)
c
        if (ipch .ge. 0) then
c
c         Process and write the B-dot(H) data.
c
c         Calling sequence substitutions:
c           bdoth for avgrid
c
          call intrp(aamatr,apr,bdoth,cof,eps100,gmmatr,ipivot,
     $    narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $    xvec,yvec)
c
          ux24 = 'bdoth'
          j2 = ilnobl(ux24)
          write (ndata1) ux24
          write (ndat1f,1010) ux24(1:j2)
          write (noutpt,1020) ux24(1:j2)
c
          do ntpr = 1,ntprt
            nt = narxt(ntpr)
            write (ndata1) (apr(i,ntpr), i = 1,nt)
            write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
          enddo
c
          write (ndat1f,1010) utermc(1:72)
c
c         Process and write the dnB-dot(H)/dPn data.
c
          do ipc = 1,ipch
c
            do ntpr = 1,ntprt
              do n = 1,narxt(ntpr)
                avgrid(n,ntpr) = dbdth(n,ntpr,ipc)
              enddo
            enddo
c
            call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,
     $      narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $      xvec,yvec)
c
            ux24 = 'dbdth( )'
            write (ux24(7:7),'(i1)') ipc
            j2 = ilnobl(ux24)
            write (ndata1) ux24
            write (ndat1f,1010) ux24(1:j2)
            write (noutpt,1020) ux24(1:j2)
c
            do ntpr = 1,ntprt
              nt = narxt(ntpr)
              write (ndata1) (apr(i,ntpr), i = 1,nt)
              write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
            enddo
c
            write (ndat1f,1010) utermc(1:72)
          enddo
        endif
c
        if (ipcv .ge. 0) then
c
c         Process and write the B-dot(V) data.
c
c         Calling sequence substitutions:
c           bdotv for avgrid
c
          call intrp(aamatr,apr,bdotv,cof,eps100,gmmatr,ipivot,
     $    narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $    xvec,yvec)
c
          ux24 = 'bdotv'
          j2 = ilnobl(ux24)
          write (ndata1) ux24
          write (ndat1f,1010) ux24(1:j2)
          write (noutpt,1020) ux24(1:j2)
c
          do ntpr = 1,ntprt
            nt = narxt(ntpr)
            write (ndata1) (apr(i,ntpr), i = 1,nt)
            write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
          enddo
c
          write (ndat1f,1010) utermc(1:72)
c
c         Process and write the dnB-dot(V)/dPn data.
c
          do ipc = 1,ipcv
c
            do ntpr = 1,ntprt
              do n = 1,narxt(ntpr)
                avgrid(n,ntpr) = dbdtv(n,ntpr,ipc)
              enddo
            enddo
c
            call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,
     $      narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $      xvec,yvec)
c
            ux24 = 'dbdtv( )'
            write (ux24(7:7),'(i1)') ipc
            j2 = ilnobl(ux24)
            write (ndata1) ux24
            write (ndat1f,1010) ux24(1:j2)
            write (noutpt,1020) ux24(1:j2)
c
            do ntpr = 1,ntprt
              nt = narxt(ntpr)
              write (ndata1) (apr(i,ntpr), i = 1,nt)
              write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
            enddo
c
            write (ndat1f,1010) utermc(1:72)
          enddo
        endif
c
c       Write the CCO2 coefficients.
c
        ux24 = 'cco2'
        j2 = ilnobl(ux24)
        write (ndata1) ux24
        write (ndat1f,1010) ux24(1:j2)
        write (noutpt,1020) ux24(1:j2)
c
        write (ndata1) (cco2(i), i = 1,5)
        write (ndat1f,1030) (cco2(i), i = 1,5)
        write (ndat1f,1010) utermc(1:72)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Process and write the next parameter only if processing a
c     "Pitzer" type data file.
c
      if (uakey(1:8) .eq. 'Pitzer  ') then
c
c       Process and write the A(phi) data.
c
c       Calling sequence substitutions:
c         aphi for avgrid
c
        call intrp(aamatr,apr,aphi,cof,eps100,gmmatr,ipivot,
     $  narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $  xvec,yvec)
c
        ux24 = 'aphi'
        j2 = ilnobl(ux24)
        write (ndata1) ux24
        write (ndat1f,1010) ux24(1:j2)
        write (noutpt,1020) ux24(1:j2)
c
        do ntpr = 1,ntprt
          nt = narxt(ntpr)
          write (ndata1) (apr(i,ntpr), i = 1,nt)
          write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
        enddo
c
        write (ndat1f,1010) utermc(1:72)
c
        if (ipch .ge. 0) then
c
c         Process and write the A(H) data.
c
c         Calling sequence substitutions:
c           adhh for avgrid
c
          call intrp(aamatr,apr,adhh,cof,eps100,gmmatr,ipivot,
     $    narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $    xvec,yvec)
c
          ux24 = 'adhh'
          j2 = ilnobl(ux24)
          write (ndata1) ux24
          write (ndat1f,1010) ux24(1:j2)
          write (noutpt,1020) ux24(1:j2)
c
          do ntpr = 1,ntprt
            nt = narxt(ntpr)
            write (ndata1) (apr(i,ntpr), i = 1,nt)
            write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
          enddo
c
          write (ndat1f,1010) utermc(1:72)
c
c         Process and write the dnA(H)/dPn data.
c
          do ipc = 1,ipch
c
            do ntpr = 1,ntprt
              do n = 1,narxt(ntpr)
                avgrid(n,ntpr) = dadhh(n,ntpr,ipc)
              enddo
            enddo
c
            call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,
     $      narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $      xvec,yvec)
c
            ux24 = 'dadhh( )'
            write (ux24(7:7),'(i1)') ipc
            j2 = ilnobl(ux24)
            write (ndata1) ux24
            write (ndat1f,1010) ux24(1:j2)
            write (noutpt,1020) ux24(1:j2)
c
            do ntpr = 1,ntprt
              nt = narxt(ntpr)
              write (ndata1) (apr(i,ntpr), i = 1,nt)
              write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
            enddo
c
            write (ndat1f,1010) utermc(1:72)
          enddo
        endif
c
        if (ipcv .ge. 0) then
c
c         Process and write the A(V) data.
c
c         Calling sequence substitutions:
c           adhv for avgrid
c
          call intrp(aamatr,apr,adhv,cof,eps100,gmmatr,ipivot,
     $    narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $    xvec,yvec)
c
          ux24 = 'adhv'
          j2 = ilnobl(ux24)
          write (ndata1) ux24
          write (ndat1f,1010) ux24(1:j2)
          write (noutpt,1020) ux24(1:j2)
c
          do ntpr = 1,ntprt
            nt = narxt(ntpr)
            write (ndata1) (apr(i,ntpr), i = 1,nt)
            write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
          enddo
c
          write (ndat1f,1010) utermc(1:72)
c
c         Process and write the dnA(V)/dPn data.
c
          do ipc = 1,ipcv
c
            do ntpr = 1,ntprt
              do n = 1,narxt(ntpr)
                avgrid(n,ntpr) = dadhv(n,ntpr,ipc)
              enddo
            enddo
c
            call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,
     $      narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $      xvec,yvec)
c
            ux24 = 'dadhv( )'
            write (ux24(7:7),'(i1)') ipc
            j2 = ilnobl(ux24)
            write (ndata1) ux24
            write (ndat1f,1010) ux24(1:j2)
            write (noutpt,1020) ux24(1:j2)
c
            do ntpr = 1,ntprt
              nt = narxt(ntpr)
              write (ndata1) (apr(i,ntpr), i = 1,nt)
              write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
            enddo
c
            write (ndat1f,1010) utermc(1:72)
          enddo
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Process and write the thermodynamic data for the "Eh" reaction.
c
c     Begin with the log K data.
c
c     Calling sequence substitutions:
c       xlke for avgrid
c
      call intrp(aamatr,apr,xlke,cof,eps100,gmmatr,ipivot,
     $ narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $ xvec,yvec)
c
      ux24 = 'xlke'
      j2 = ilnobl(ux24)
      write (ndata1) ux24
      write (ndat1f,1010) ux24(1:j2)
      write (noutpt,1020) ux24(1:j2)
c
      do ntpr = 1,ntprt
        nt = narxt(ntpr)
        write (ndata1) (apr(i,ntpr), i = 1,nt)
        write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
      enddo
c
      write (ndat1f,1010) utermc(1:72)
c
      if (ipch .ge. 0) then
c
c       Process and write the enthalpy of reaction data.
c
c       Calling sequence substitutions:
c         xhfe for avgrid
c
        call intrp(aamatr,apr,xhfe,cof,eps100,gmmatr,ipivot,
     $  narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $  xvec,yvec)
c
        ux24 = 'xhfe'
        j2 = ilnobl(ux24)
        write (ndata1) ux24
        write (ndat1f,1010) ux24(1:j2)
        write (noutpt,1020) ux24(1:j2)
c
        do ntpr = 1,ntprt
          nt = narxt(ntpr)
          write (ndata1) (apr(i,ntpr), i = 1,nt)
          write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
        enddo
c
        write (ndat1f,1010) utermc(1:72)
c
c       Process and write the pressure derivatives of the enthalpy
c       of reaction.
c
        do ipc = 1,ipch
c
          do ntpr = 1,ntprt
            do n = 1,narxt(ntpr)
              avgrid(n,ntpr) = dhfe(n,ntpr,ipc)
            enddo
          enddo
c
          call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,
     $    narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $    xvec,yvec)
c
          ux24 = 'dhfe( )'
          write (ux24(7:7),'(i1)') ipc
          j2 = ilnobl(ux24)
          write (ndata1) ux24
          write (ndat1f,1010) ux24(1:j2)
          write (noutpt,1020) ux24(1:j2)
c
          do ntpr = 1,ntprt
            nt = narxt(ntpr)
            write (ndata1) (apr(i,ntpr), i = 1,nt)
            write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
          enddo
c
          write (ndat1f,1010) utermc(1:72)
        enddo
      endif
c
      if (ipcv .ge. 0) then
c
c       Process and write the volume of reaction data.
c
c       Calling sequence substitutions:
c         xvfe for avgrid
c
        call intrp(aamatr,apr,xvfe,cof,eps100,gmmatr,ipivot,
     $  narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $  xvec,yvec)
c
        ux24 = 'xvfe'
        j2 = ilnobl(ux24)
        write (ndata1) ux24
        write (ndat1f,1010) ux24(1:j2)
        write (noutpt,1020) ux24(1:j2)
c
        do ntpr = 1,ntprt
          nt = narxt(ntpr)
          write (ndata1) (apr(i,ntpr), i = 1,nt)
          write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
        enddo
c
        write (ndat1f,1010) utermc(1:72)
c
c       Process and write the pressure derivatives of the volume
c       of reaction.
c
        do ipc = 1,ipcv
c
          do ntpr = 1,ntprt
            do n = 1,narxt(ntpr)
              avgrid(n,ntpr) = dvfe(n,ntpr,ipc)
            enddo
          enddo
c
          call intrp(aamatr,apr,avgrid,cof,eps100,gmmatr,ipivot,
     $    narxmx,narxt,noutpt,ntprmx,ntprt,nttyo,tempc,tempcs,tmpcmx,
     $    xvec,yvec)
c
          ux24 = 'dvfe( )'
          write (ux24(7:7),'(i1)') ipc
          j2 = ilnobl(ux24)
          write (ndata1) ux24
          write (ndat1f,1010) ux24(1:j2)
          write (noutpt,1020) ux24(1:j2)
c
          do ntpr = 1,ntprt
            nt = narxt(ntpr)
            write (ndata1) (apr(i,ntpr), i = 1,nt)
            write (ndat1f,1030) (apr(i,ntpr), i = 1,nt)
          enddo
c
          write (ndat1f,1010) utermc(1:72)
        enddo
      endif
c
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
