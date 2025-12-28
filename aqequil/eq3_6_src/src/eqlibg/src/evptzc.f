      subroutine evptzc(amu,aslm,ipbtmx,jpfcmx,jptffl,nmut,nmutmx,
     $ noutpt,nttyo,nslt,nsltmx,pmu,pslamn,tempc)
c
c     This subroutine computes the S-lambda(n) (pslamn) and mu (pmu)
c     coefficients of Pitzer's equations at the temperature tempc.
c
c     This subroutine is called by:
c
c       EQLIB/evdata.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       amu    = array of coefficients for calculating Pitzer mu
c                coefficients (pmu) as a function of temperature
c       aslm   = array of coefficients for calculating Pitzer
c                  S-lambda(n) coefficients (pslamn) as a function
c                  of temperature
c       jptffl = Pitzer parameter temperature coefficient flag:
c                  -1 = 25C centric Taylor's series truncated at
c                         second order
c                   0 = 25C centric LLNL 5TERM function
c                         (maximal 5th order)
c                   1 = non-25C centric eight-term Greenberg and
c                         Moller (1989) function
c       nmut   = the number of species triples for which there are
c                  mu coefficients
c       nslt   = the number of species pairs for which there are
c                  S-lambda(n) coefficients
c       tempc  = the temperature (C)
c
c     Principal output:
c
c       pmu    = array of Pitzer mu coefficients
c       pslamn = array of Pitzer S-lambda(n) coefficients
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ipbtmx,jpfcmx,nmutmx,nsltmx
c
      integer noutpt,nttyo
c
      integer jptffl,nmut,nslt
c
      real(8) amu(jpfcmx,nmutmx),aslm(jpfcmx,0:ipbtmx,nsltmx),
     $ pmu(nmutmx),pslamn(0:ipbtmx,nsltmx)
c
      real(8) tempc
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j,jpfclm,j2,k,nmu,nsl
c
      integer ilnobl
c
      character(len=8) ux8
c
      real(8) px,t,tr
c
c     Note: the dimension of dtf must match the maximum number of terms
c     in any programmed temperature function.
c
      real(8) dtf(2:8)
c
c-----------------------------------------------------------------------
c
      t = tempc + 273.15
      tr = 298.15
c
      if (jptffl .eq. -1) then
c
c       Classical 25C-centric Taylor's series truncated at second order.
c
        dtf(2) = t - tr
        dtf(3) = 0.5*dtf(2)*dtf(2)
        jpfclm = min(jpfcmx,3)
c
      elseif (jptffl .eq. 0) then
c
c       LLNL maximal 5-term equation.
c
        dtf(2) = (1./t) - (1./tr)
        dtf(3) = log(t/tr)
        dtf(4) = t - tr
        dtf(5) = t**2 - tr**2
        jpfclm = min(jpfcmx,5)
c
      elseif (jptffl .eq. 1) then
c
c       Greenberg and Moller (1989) eight-term equation.
c
        if ((abs(t - 227.) .le. 1.e-3)
     $    .or. (abs(t - 263.) .le. 1.e-3)
     $    .or. (abs(t - 680.) .le. 1.e-3)) then
          write (noutpt,1000) t,tempc
          write (nttyo,1000) t,tempc
 1000     format(/' * Error - (EQLIBG/evptzc) The temperature of ',
     $    f7.4,'K (',f7.4,'C) is',/7x,'too close to a singularity',
     $    ' in the TEQUIL eight-parameter Pitzer',/7x,'parameter',
     $    ' temperature function.')
          stop
        endif
c
        dtf(2) = t
        dtf(3) = 1./t
        dtf(4) = log(t)
        dtf(5) = 1./(t - 263.)
        dtf(6) = t*t
        dtf(7) = 1./(680. - t)
        dtf(8) = 1./(t - 227.)
        jpfclm = min(jpfcmx,8)
c
      else
        write (ux8,'(i5)') jptffl
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1010) ux8(1:j2)
        write (nttyo,1010) ux8(1:j2)
 1010   format(/' * Error - (EQLIBG/evptzc) Programming error trap:',
     $  ' Have an',/7x,'illegal value of ',a,' for the variable',
     $  ' jptffl, which determines',/7x,'the Pitzer parameter',
     $  ' temperature function. Do not have such a',/7x,'function',
     $  ' programmed that matches this value.')
        stop
      endif
c
      do nsl = 1,nslt
        do k = 0,2
          px = aslm(1,k,nsl)
          do j = 2,jpfclm
            px = px + aslm(j,k,nsl)*dtf(j)
          enddo
          pslamn(k,nsl) = px
        enddo
      enddo
c
      do nmu = 1,nmut
        px = amu(1,nmu)
        do j = 2,jpfclm
          px = px + amu(j,nmu)*dtf(j)
        enddo
        pmu(nmu) = px
      enddo
c
      end
