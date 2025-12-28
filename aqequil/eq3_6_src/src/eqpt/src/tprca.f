      subroutine tprca(abeta,acphi,alamca,alpha,alphca,amua2c,amuc2a,
     $ icapr,ipbtmx,jpfcmx,natmax,na2ctr,ncapr,ncvca,nc2atr,nerr,
     $ noutpt,npx2mx,npx2t,nttyo,nwarn,pcvca,qpdca,uaqsp,upair,zaqsp)
c
c     Test and process the Pitzer data for ca (cation, anion) pairs
c     read from the DATA0 file. Find and flag errors, such as duplication
c     of data (e.g., two data blocks for the same ca pair). Calculate
c     the conventional primitive Pitzer parameters from the observable
c     compound parameters read from the data file:
c
c       beta(n)(ca) -> lambda(n)(ca)   (n = 0,2)
c       Cphi(ca)    -> mu(cca) and mu(aac)
c
c     Check the coverage of entered Pitzer data against all possible
c     ca pairs that can be composed of the aqueous cations and anions
c     present on the data file.
c
c     This subroutine is called by:
c
c       EQPT/eqpt.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       nat    = the number of aqueous species
c       uaqsp  = array of names of aqueous species
c
c     Principal output:
c
c       nerr   = cumulative error counter
c       nwarn  = cumulative warning counter
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ipbtmx,jpfcmx,natmax,na2ctr,ncapr,nc2atr,npx2mx
c
      integer noutpt,nttyo
c
      integer icapr(2,ncapr)
c
      integer ncvca,nerr,npx2t,nwarn
c
      logical qpdca(ncapr)
c
      character(len=24) uaqsp(natmax),upair(2,npx2mx)
c
      real*8 abeta(jpfcmx,0:ipbtmx,npx2mx),acphi(jpfcmx,npx2mx),
     $ alamca(jpfcmx,0:ipbtmx,ncapr),alphca(ipbtmx,ncapr),
     $ alpha(2,npx2mx),amua2c(jpfcmx,na2ctr),amuc2a(jpfcmx,nc2atr),
     $ zaqsp(natmax)
c
      real*8 pcvca
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,j,jp,jpair,j2,j3,j4,j5,n,ncount,ndupl,nlistl,nn,
     $ nodatc
c
      integer ilnobl
c
      character*56 ustr56
      character*24 unam1,unam2
      character*8 ux8
c
      real*8 fmuc2a,fmua2c,z1,z2
c
c-----------------------------------------------------------------------
c
c     Limit on the list of pairs for which no data were found.
c
      data nlistl / 20 /
c
c-----------------------------------------------------------------------
c
c     Initialize the data arrays.
c
      do n = 1,ncapr
        qpdca(n) = .false.
        do i = 1,ipbtmx
          alphca(i,n) = 0.
        enddo
      enddo
c
      do n = 1,ncapr
        do i = 0,ipbtmx
          do j = 1,jpfcmx
            alamca(j,i,n) = 0.
          enddo
        enddo
      enddo
c
      do n = 1,ncapr
        do j = 1,jpfcmx
          amuc2a(j,n) = 0.
          amua2c(j,n) = 0.
        enddo
      enddo
c
c     Check the entered data for ca pairs.
c
      nodatc = 0
c
      do n = 1,ncapr
        i = icapr(1,n)
        j = icapr(2,n)
        unam1 = uaqsp(i)
        unam2 = uaqsp(j)
        z1 = zaqsp(i)
        z2 = zaqsp(j)
c
c       Search for unam1, unam2 in the upair array.
c       That array corresponds to the species pairs blocks.
c
        call srch22(jpair,unam1,unam2,upair,npx2mx,npx2t)
c
        if (jpair .gt. 0) then
c
c         Have found an entry.
c
          qpdca(n) = .true.
c
c         Store the data.
c
          do i = 1,ipbtmx
            alphca(i,n) = alpha(i,jpair)
          enddo
c
          do i = 0,ipbtmx
            do j = 1,jpfcmx
              alamca(j,i,n) = abeta(j,i,jpair)
            enddo
          enddo
c
          fmuc2a = sqrt(-z1/z2)/6.
          fmua2c = sqrt(-z2/z1)/6.
c
          do j = 1,jpfcmx
            amuc2a(j,n)   = acphi(j,jpair)*fmuc2a
            amua2c(j,n)   = acphi(j,jpair)*fmua2c
          enddo
c
c         Check for duplicate data sets.
c
          ndupl = 0
          do jp = jpair + 1,npx2t
            if (unam1(1:24) .eq. upair(1,jp)(1:24)) then
              if (unam2(1:24) .eq. upair(2,jp)(1:24)) then
                ndupl = ndupl + 1
              endif
            endif
          enddo
c
          if (ndupl .gt. 0) then
            j2 = ilnobl(unam1)
            j3 = ilnobl(unam2)
            if (ndupl .eq. 1) then
              write (noutpt,1010) unam1(1:j2),unam2(1:j3)
              write (nttyo,1010) unam1(1:j2),unam2(1:j3)
 1010         format(/' * Error - (EQPT/tprca) Have found a duplicate',
     $        ' data block on the DATA0 file',/7x,'for the ca pair ',
     $        a,', ',a,'.')
            else
              ux8 = ' '
              write (ux8,'(i5)') ndupl
              call lejust(ux8)
              j5 = ilnobl(ux8)
              write (noutpt,1020) ux8(1:j5),unam1(1:j2),unam2(1:j3)
              write (nttyo,1020) ux8(1:j5),unam1(1:j2),unam2(1:j3)
 1020         format(/' * Error - (EQPT/tprca) Have found ',a,
     $        ' duplicate data blocks on the DATA0 file',
     $        /7x,'for the ca pair ',a,', ',a,'.')
            endif
            nerr = nerr + ndupl
          endif
c
        else
c
c         No data block was found on the DATA0 file.
c         Note that qpdca(n) is left with a value of .false.
c
          nodatc = nodatc + 1
        endif
c
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nodatc .gt. 0) then
        write (noutpt,1040)
        write (nttyo,1040)
 1040   format(/' * Warning - (EQPT/tprca) Did not find a data',
     $  ' block on the DATA0 file',/7x,'for any of the following',
     $  ' ca pairs:',/)
c
        ncount = 0
        do n = 1,ncapr
          if (.not.qpdca(n)) then
            ncount = ncount + 1
            i = icapr(1,n)
            j = icapr(2,n)
            unam1 = uaqsp(i)
            unam2 = uaqsp(j)
            j2 = ilnobl(unam1)
            j3 = ilnobl(unam2)
            ustr56 = unam1(1:j2) // ', ' // unam2(1:j3)
            j4 = ilnobl(ustr56)
            write (noutpt,1050) ustr56(1:j4)
            write (nttyo,1050) ustr56(1:j4)
 1050       format(9x,a)
            if (ncount .eq. nlistl) go to 200
          endif
        enddo
  200   continue
c
        nn = nodatc - ncount
        if (nn .gt. 0) then
          write (ux8,'(i5)') nn
          call lejust(ux8)
          j3 = ilnobl(ux8)
          write (noutpt,1060) ux8(1:j3)
          write (nttyo,1060) ux8(1:j3)
 1060     format(/9x,'plus ',a,' others')
        endif
        write (noutpt,1070)
        write (nttyo,1070)
 1070   format(1x)
        nwarn = nwarn + 1
      endif
c
      ncvca = ncapr - nodatc
      pcvca = (100.*ncvca)/float(ncapr)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
