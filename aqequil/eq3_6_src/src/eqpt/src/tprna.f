      subroutine tprna(abeta,alamna,inapr,ipbtmx,jpfcmx,natmax,
     $ ncvna,nerr,nnapr,noutpt,npx2mx,npx2t,nttyo,nwarn,pcvna,
     $ qpdna,uaqsp,upair)
c
c     Test and process the na (neutral, anion) pair Pitzer data read
c     from the DATA0 file. Find and flag errors, such as duplication of
c     data (e.g., two data blocks for the same na pair). The
c     conventional primitive Pitzer parameters are identical to the
c     observable parameters read from the data file:
c
c       lambda(na) -> lambda(na)
c
c     Check the coverage of entered Pitzer data against all possible
c     na pairs that can be composed of the aqueous neutral species and
c     anions present on the data file.
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
      integer ipbtmx,jpfcmx,natmax,nnapr,npx2mx
c
      integer noutpt,nttyo
c
      integer inapr(2,nnapr)
c
      integer ncvna,nerr,npx2t,nwarn
c
      logical qpdna(nnapr)
c
      character(len=24) uaqsp(natmax),upair(2,npx2mx)
c
      real*8 abeta(jpfcmx,0:ipbtmx,npx2mx),
     $ alamna(jpfcmx,0:ipbtmx,nnapr)
c
      real*8 pcvna
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
      character*24 unam1,unam2
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
      do n = 1,nnapr
        qpdna(n) = .false.
        do j = 1,jpfcmx
          do i = 0,ipbtmx
            alamna(j,i,n) = 0.
          enddo
        enddo
      enddo
c
c     Check the entered data for na pairs.
c
      nodatc = 0
c
      do n = 1,nnapr
        i = inapr(1,n)
        j = inapr(2,n)
        unam1 = uaqsp(i)
        unam2 = uaqsp(j)
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
          qpdna(n) = .true.
c
c         Store the data.
c
          do j = 1,jpfcmx
            do i = 0,ipbtmx
              alamna(j,i,n) = abeta(j,i,jpair)
            enddo
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
 1010         format(/' * Error - (EQPT/tprna) Have found a',
     $        ' duplicate data block on the DATA0 file',
     $        /7x,'for the na pair ',a,', ',a,'.')
            else
              ux8 = ' '
              write (ux8,'(i5)') ndupl
              call lejust(ux8)
              j5 = ilnobl(ux8)
              write (noutpt,1020) ux8(1:j5),unam1(1:j2),unam2(1:j3)
              write (nttyo,1020) ux8(1:j5),unam1(1:j2),unam2(1:j3)
 1020         format(/' * Error - (EQPT/tprna) Have found ',a,
     $        ' duplicate data blocks on the DATA0 file',
     $        /7x,'for the na pair ',a,', ',a,'.')
            endif
            nerr = nerr + ndupl
          endif
c
        else
c
c         No data block was found on the DATA0 file.
c         Note that qpdna(n) is left with a value of .false.
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
 1040   format(/' * Warning - (EQPT/tprna) Did not find a data',
     $  ' block on the DATA0 file',/7x,'for any of the following',
     $  ' na pairs:',/)
c
        ncount = 0
        do n = 1,nnapr
          if (.not.qpdna(n)) then
            ncount = ncount + 1
            i = inapr(1,n)
            j = inapr(2,n)
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
      ncvna = nnapr - nodatc
      pcvna = (100.*ncvna)/float(nnapr)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
