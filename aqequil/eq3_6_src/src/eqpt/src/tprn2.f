      subroutine tprn2(abeta,acphi,alamn2,amun3,in2pr,ipbtmx,
     $ jpfcmx,natmax,ncvn2,nerr,nn2pr,nn3tr,noutpt,npx2mx,npx2t,
     $ nttyo,nwarn,pcvn2,qpdn2,uaqsp,upair)
c
c     Test and process Pitzer data read from the data file that
c     pertain to a single neutral species (nn and nnn combinations,
c     such as CO2(aq)-CO2(aq) and CO2(aq)-CO2(aq)-CO2(aq). Find and
c     flag errors, such as duplication of data (e.g., two data
c     blocks for the same pair composed of a repeated neutral).
c     The conventional primitive Pitzer parameters are identical
c     to the observable parameters read from the data file:
c
c       lambda(nn) -> lambda(nn)
c       mu(nnn)    -> mu(nnn)
c
c     Check the coverage of entered Pitzer data against all possible
c     pairs that can be composed by repeating a neutral species present
c     on the data file.
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
      integer ipbtmx,jpfcmx,natmax,nn2pr,nn3tr,npx2mx
c
      integer noutpt,nttyo
c
      integer in2pr(nn2pr)
c
      integer ncvn2,nerr,npx2t,nwarn
c
      logical qpdn2(nn2pr)
c
      character(len=24) uaqsp(natmax),upair(2,npx2mx)
c
      real*8 abeta(jpfcmx,0:ipbtmx,npx2mx),acphi(jpfcmx,npx2mx),
     $ alamn2(jpfcmx,0:ipbtmx,nn2pr),amun3(jpfcmx,nn3tr)
c
      real*8 pcvn2
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,j,jp,jpair,j2,j3,j4,j5,n,ncount,ndupl,nlistl,
     $ nn,nodatc
c
      integer ilnobl
c
      character*56 ustr56
      character*24 unam1
      character*8 ux8
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
      do n = 1,nn2pr
        qpdn2(n) = .false.
        do j = 1,jpfcmx
          do i = 0,ipbtmx
            alamn2(j,i,n) = 0.
          enddo
        enddo
      enddo
c
c     Note: nn3tr = nn2pr.
c
      do n = 1,nn3tr
        do j = 1,jpfcmx
          amun3(j,n) = 0.
        enddo
      enddo
c
c     Check the entered data for nn pairs.
c
      nodatc = 0
c
      do n = 1,nn2pr
        i = in2pr(n)
        unam1 = uaqsp(i)
c
c       Search for unam1 in the upair array.
c
        call srch22(jpair,unam1,unam1,upair,npx2mx,npx2t)
c
        if (jpair .gt. 0) then
c
c         Have found an entry.
c
          qpdn2(n) = .true.
c
c         Store the data.
c
          do j = 1,jpfcmx
            do i = 0,ipbtmx
              alamn2(j,i,n) = abeta(j,i,jpair)
            enddo
          enddo
c
          do j = 1,jpfcmx
            amun3(j,n)   = acphi(j,jpair)
          enddo
c
c         Check for duplicate data sets.
c
          ndupl = 0
          do jp = jpair + 1,npx2t
            if (unam1(1:24) .eq. upair(1,jp)(1:24)) then
              if (unam1(1:24) .eq. upair(2,jp)(1:24)) then
                ndupl = ndupl + 1
              endif
            endif
          enddo
c
          if (ndupl .gt. 0) then
            j2 = ilnobl(unam1)
            if (ndupl .eq. 1) then
              write (noutpt,1010) unam1(1:j2),unam1(1:j2)
              write (nttyo,1010) unam1(1:j2),unam1(1:j2)
 1010         format(/' * Error - (EQPT/tprn2) Have found a',
     $        ' duplicate data block on the DATA0 file',
     $        /7x,'for the nn pair ',a,', ',a,'.')
            else
              ux8 = ' '
              write (ux8,'(i5)') ndupl
              call lejust(ux8)
              j5 = ilnobl(ux8)
              write (noutpt,1020) ux8(1:j5),unam1(1:j2),unam1(1:j2)
              write (nttyo,1020) ux8(1:j5),unam1(1:j2),unam1(1:j2)
 1020         format(/' * Error - (EQPT/tprn2) Have found ',a,
     $        ' duplicate data blocks on the DATA0 file',
     $        /7x,'for the nn pair ',a,', ',a,'.')
            endif
            nerr = nerr + ndupl
          endif
c
        else
c
c         No data block was found on the DATA0 file.
c         Note that qpdn2(n) is left with a value of .false.
c
          nodatc = nodatc + 1
        endif
c
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nodatc .gt. 0) then
        write (noutpt,1050)
        write (nttyo,1050)
 1050   format(/' * Warning - (EQPT/tprn2) Did not find a data',
     $  ' block on the DATA0 file',/7x,'for any of the following',
     $  ' nn pairs:',/)
c
        ncount = 0
        do n = 1,nn2pr
          if (.not.qpdn2(n)) then
            ncount = ncount + 1
            i = in2pr(n)
            unam1 = uaqsp(i)
            j2 = ilnobl(unam1)
            ustr56 = unam1(1:j2) // ', ' // unam1(1:j2)
            j4 = ilnobl(ustr56)
            write (noutpt,1060) ustr56(1:j4)
            write (nttyo,1060) ustr56(1:j4)
 1060       format(9x,a)
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
          write (noutpt,1070) ux8(1:j3)
          write (nttyo,1070) ux8(1:j3)
 1070     format(/9x,'plus ',a,' others')
        endif
        write (noutpt,1080)
        write (nttyo,1080)
 1080   format(1x)
        nwarn = nwarn + 1
      endif
c
      ncvn2 = nn2pr - nodatc
      pcvn2 = (100.*ncvn2)/float(nn2pr)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
