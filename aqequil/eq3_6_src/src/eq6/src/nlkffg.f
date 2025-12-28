      subroutine nlkffg(axlks,cess,cdrs,iffg,ifrn1,ifrn2,jffg,jpflag,
     $ jsflag,mwtsp,narxmx,ncmpr,ndrs,ndrsmx,ndrsr,ness,nessmx,nessr,
     $ nffg,nffgmx,nfrn1,nfrn2,ngrn1,ngrn2,noutpt,nphasx,npt,nptmax,
     $ nst,nstmax,ntpr,ntprmx,nttyo,qcntmp,uffg,ufixf,uphase,uspec,
     $ vosp0,xlkffg)
c
c     This subroutine sets up ficitive minerals, each of which is used
c     to fixing the fugacity of a specified gas. The fugacity is
c     actually fixed only if the corresponding ficitive mineral is in
c     equilibrium with the aqueous system. Otherwise, the fugacity may
c     be less than the specified fugacity value.
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
      integer narxmx,ndrsmx,nessmx,nffgmx,nptmax,nstmax,ntprmx
c
      integer noutpt,nttyo
c
      integer iffg(nffgmx),jffg(nffgmx),jpflag(nptmax),jsflag(nstmax),
     $ ncmpr(2,nptmax),ness(nessmx),nessr(2,nstmax),ndrs(ndrsmx),
     $ ndrsr(2,nstmax),nphasx(nstmax)
c
      integer ifrn1,ifrn2,nffg,nfrn1,nfrn2,ngrn1,ngrn2,npt,nst,ntpr
c
      logical qcntmp
c
      character*48 uspec(nstmax)
      character*24 uffg(nffgmx),uphase(nptmax)
      character*8 ufixf
c
      real*8 axlks(narxmx,ntprmx,nstmax),cess(nessmx),cdrs(ndrsmx),
     $ mwtsp(nstmax),vosp0(nstmax),xlkffg(nffgmx)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,j,j2,k,n,nerr,nerr1,ng,nn,np,nrf1,nrf2,nr1,nr2,ns,
     $ nsl,nss
c
      integer ilnobl
c
      character*24 ux24
c
c-----------------------------------------------------------------------
c
      nerr = 0
      np = npt
      ns = nst
      nfrn1 = nst + 1
      ifrn1 = npt + 1
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do n = 1,nffg
c
        do nss = ngrn1,ngrn2
          if (uffg(n)(1:24) .eq. uspec(nss)(1:24)) go to 105
        enddo
        j2 = ilnobl(uffg(n))
        write (noutpt,1000) uffg(n)(1:j2)
        write (nttyo,1000) uffg(n)(1:j2)
 1000   format (/' * Error - (EQ6/nlkffg) ',a,' was not on the',
     $  /7x,'data file, so its fugacity can not be fixed.')
        nerr = nerr + 1
        go to 200
c
  105   nerr1 = 0
c
        if (npt .ge. nptmax) then
          j2 = ilnobl(uffg(n))
          write (noutpt,1010) nptmax,uffg(n)(1:j2)
          write (nttyo,1010) nptmax,uffg(n)(1:j2)
 1010     format (/' * Error - (EQ6/nlkffg) The maximum ',i3,' phases',
     $    /7x,'would be exceeded by by creating a fugacity fixing',
     $    ' phase',
     $    /7x,'for ',a,'. Increase the dimensioning variable nptpar.')
          nerr1 = nerr1 + 1
          nerr = nerr + 1
        endif
c
        if (nst .ge. nstmax) then
          j2 = ilnobl(uffg(n))
          write (noutpt,1020) nstmax,uffg(n)(1:j2)
          write (nttyo,1020) nstmax,uffg(n)(1:j2)
 1020     format (/' * Error - (EQ6/nlkffg) The maximum ',i3,' species',
     $    /7x,'would be exceeded by creating a fugacity fixing phase',
     $    /7x,'for ',a,'. Increase the dimensioning parameter nstpar.')
          nerr1 = nerr1 + 1
          nerr = nerr + 1
        endif
c
        if (nerr1 .gt. 0) go to 200
c
        np = np + 1
        npt = np
        ns = ns + 1
        nsl = nst
        nst = ns
        iffg(n) = ns
        nphasx(ns) = np
        ncmpr(1,np) = ns
        ncmpr(2,np) = ns
        ng = nss - ngrn1 + 1
        jffg(n) = ng
        jpflag(np) = 0
        jsflag(ns) = 0
c
        ux24 = ufixf(1:5)
        ux24(6:24) = uspec(nss)(1:19)
        uspec(ns)(1:24) = ux24
        uspec(ns)(25:48) = ux24
        uphase(np) = ux24
c
        mwtsp(ns) = mwtsp(nss)
        vosp0(ns) = 0.
c
        nr1 = nessr(1,nss)
        nr2 = nessr(2,nss)
        nrf1 = nessr(2,nsl) + 1
        nrf2 = nrf1 + nr2 - nr1
        if (nrf2 .gt. nessmx) then
          ux24 = uspec(ns)(1:24)
          j2 = ilnobl(ux24)
          write (noutpt,1030) nessmx,uspec(ns)(1:j2)
          write (nttyo,1030) nessmx,uspec(ns)(1:j2)
 1030     format (/' * Error - (EQ6/nlkffg) The maximum ',i5,' entries',
     $    /7x,'in the cess/ness arrays would be exceeded by creating',
     $    /7x,'the species ',a,'. Increase the dimensioning',
     $    /7x,'parameter nesspa.')
          nerr = nerr + 1
          go to 115
        endif
        nessr(1,ns) = nrf1
        nessr(2,ns) = nrf2
        k = nr1 - 1
        do nn = nrf1,nrf2
          k = k + 1
          cess(nn) = cess(k)
          ness(nn) = ness(k)
        enddo
c
  115   nr1 = ndrsr(1,nss)
        nr2 = ndrsr(2,nss)
        nrf1 = ndrsr(2,nsl) + 1
        nrf2 = nrf1 + nr2 - nr1
        if (nrf2 .gt. ndrsmx) then
          ux24 = uspec(ns)(1:24)
          j2 = ilnobl(ux24)
          write (noutpt,1040) ndrsmx,uspec(ns)(1:j2)
          write (nttyo,1040) ndrsmx,uspec(ns)(1:j2)
 1040     format (/' * Error - (EQ6/nlkffg) The maximum ',i5,' entries',
     $    /7x,'in the cdrs/ndrs arrays would be exceeded by creating',
     $    /7x,'the species ',a,'. Increase the dimensioning',
     $    /7x,'parameter ndrspa.')
          nerr = nerr + 1
          go to 125
        endif
        ndrsr(1,ns) = nrf1
        ndrsr(2,ns) = nrf2
        k = nr1 - 1
        do nn = nrf1,nrf2
          k = k + 1
          cdrs(nn) = cdrs(k)
          ndrs(nn) = ndrs(k)
        enddo
        ndrs(nrf1) = ns
c
  125   do i = 1,narxmx
          do j = 1,ntprmx
            axlks(i,j,ns) = axlks(i,j,nss)
          enddo
        enddo
c
c       Alter log K values of fixed fugacity solid.
c
        do j = 1,ntprmx
          if (axlks(1,j,ns) .lt. 9999999.) then
            axlks(1,j,ns) = axlks(1,j,ns) + xlkffg(n)
          else
            if (j.eq.ntpr .or. .not.qcntmp) then
              j2 = ilnobl(uffg(n))
              write (noutpt,1050) uffg(n)(1:j2),j
              write (nttyo,1050) uffg(n)(1:j2),j
 1050         format (/" * Error - (EQ6/nlkffg) Can't fix the fugacity",
     $        /7x,' of ',a,' in temperature range ',i2,' because',
     $        /7x,'actual thermodynamic data are lacking.')
              nerr = nerr + 1
            endif
          endif
        enddo
c
  200   continue
      enddo
c
      nfrn2 = nst
      ifrn2 = npt
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nerr .gt. 0) stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
