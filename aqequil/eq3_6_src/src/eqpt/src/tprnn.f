      subroutine tprnn(abeta,acphi,alamnn,innpr,in2pr,ipbtmx,
     $ jpfcmx,natmax,ncvnn,nerr,nnnpr,nn2pr,noutpt,npx2mx,npx2t,
     $ nttyo,nwarn,pcvnn,qpdnn,qpdn2,uaqsp,upair)
c
c     Test and process the Pitzer data for nn' (neutral, different
c     neutral) pairs read from the DATA0 file. Find and flag errors,
c     such as duplication of data (e.g., two data blocks for the same
c     nn' pair). The conventional primitive Pitzer parameters are
c     identical to the observable parameters read from the data file:
c
c       lambda(nn') -> lambda(nn')
c
c     Check the coverage of entered Pitzer data against all possible
c     nn' pairs that can be composed of the neutral species present
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
      integer ipbtmx,jpfcmx,natmax,nnnpr,nn2pr,npx2mx
c
      integer noutpt,nttyo
c
      integer innpr(2,nnnpr),in2pr(nn2pr)
c
      integer ncvnn,nerr,npx2t,nwarn
c
      logical qfound1,qfound2
c
      logical qpdnn(nnnpr),qpdn2(nn2pr)
c
      character(len=24) uaqsp(natmax),upair(2,npx2mx)
c
      real*8 abeta(jpfcmx,0:ipbtmx,npx2mx),acphi(jpfcmx,npx2mx),
     $ alamnn(jpfcmx,0:ipbtmx,nnnpr)
c
      real*8 pcvnn
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
c-----------------------------------------------------------------------
c
      integer i,ii,j,jp,jpair,j2,j3,j4,j5,n,ncount,ndupl,nlistl,
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
      do n = 1,nnnpr
        qpdnn(n) = .false.
        do j = 1,jpfcmx
          do i = 0,ipbtmx
            alamnn(j,i,n) = 0.
          enddo
        enddo
      enddo
c
c     Check the entered data for nn' pairs.
c
      nodatc = 0
c
      do n = 1,nnnpr
        i = innpr(1,n)
        j = innpr(2,n)
        unam1 = uaqsp(i)
        unam2 = uaqsp(j)
        j2 = ilnobl(unam1)
        j3 = ilnobl(unam2)
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
          qpdnn(n) = .true.
c
c         Store the data.
c
          do j = 1,jpfcmx
            do i = 0,ipbtmx
              alamnn(j,i,n) = abeta(j,i,jpair)
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
            if (ndupl .eq. 1) then
              write (noutpt,1010) unam1(1:j2),unam2(1:j3)
              write (nttyo,1010) unam1(1:j2),unam2(1:j3)
 1010         format(/' * Error - (EQPT/tprnn) Have found a',
     $        ' duplicate data block on the DATA0 file',
     $        /7x,"for the nn' pair ",a,', ',a,'.')
            else
              ux8 = ' '
              write (ux8,'(i5)') ndupl
              call lejust(ux8)
              j5 = ilnobl(ux8)
              write (noutpt,1020) ux8(1:j5),unam1(1:j2),unam2(1:j3)
              write (nttyo,1020) ux8(1:j5),unam1(1:j2),unam2(1:j3)
 1020         format(/' * Error - (EQPT/tprnn) Have found ',a,
     $        ' duplicate data blocks on the DATA0 file',
     $        /7x,"for the nn' pair ",a,', ',a,'.')
            endif
            nerr = nerr + ndupl
          endif
c
        else
c
c         No data block was found on the DATA0 file.
c         Note that qpdnn(n) is left with a value of .false.
c
          nodatc = nodatc + 1
        endif
c
c       Check for the presence of the corresponding lamda(nn),
c       and lambda(n'n') data. The presence of these data on
c       the  data file is required.
c
        qfound1 = .false.
c
        do nn = 1,nn2pr
          ii = in2pr(nn)
          if (uaqsp(ii)(1:24) .eq. unam1(1:24)) then
            if (qpdn2(nn)) then
              qfound1 = .true.
              go to 130
            endif
          endif
        enddo
  130   continue
c
        qfound2 = .false.
c
        do nn = 1,nn2pr
          ii = in2pr(nn)
          if (uaqsp(ii)(1:24) .eq. unam2(1:24)) then
            if (qpdn2(nn)) then
              qfound2 = .true.
              go to 140
            endif
          endif
        enddo
  140   continue
c
        if (qpdnn(n) .and. .not.(qfound1 .and. qfound2)) then
c
c         Certain required data were not found on the DATA0 file.
c
          write (noutpt,1030) unam1(1:j2),unam2(1:j3)
          write (nttyo,1030) unam1(1:j2),unam2(1:j3)
 1030     format(/' * Error - (EQPT/tprnn) Have data on the DATA0',
     $    " file for the nn'",/7x,'pair ',a,', ',a,', but',
     $    " don't have the",/7x,'required data for the following',
     $    ' neutral-neutral pair(s):',/)
c
          if (.not.qfound1) then
            write (noutpt,1040) unam1(1:j2),unam1(1:j2)
            write (nttyo,1040) unam1(1:j2),unam1(1:j2)
 1040       format(9x,a,', ',a)
          endif
c
          if (.not.qfound2) then
            write (noutpt,1040) unam2(1:j3),unam2(1:j3)
            write (nttyo,1040) unam2(1:j3),unam2(1:j3)
          endif
c
          nerr = nerr + 1
        endif
c
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nodatc .gt. 0) then
        write (noutpt,1050)
        write (nttyo,1050)
 1050   format(/' * Warning - (EQPT/tprnn) Did not find a data',
     $  ' block on the DATA0 file',/7x,'for any of the following',
     $  " nn' pairs:",/)
c
        ncount = 0
        do n = 1,nnnpr
          if (.not.qpdnn(n)) then
            ncount = ncount + 1
            i = innpr(1,n)
            j = innpr(2,n)
            unam1 = uaqsp(i)
            unam2 = uaqsp(j)
            j2 = ilnobl(unam1)
            j3 = ilnobl(unam2)
            ustr56 = unam1(1:j2) // ', ' // unam2(1:j3)
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
      ncvnn = nnnpr - nodatc
      pcvnn = (100.*ncvnn)/float(nnnpr)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
