      subroutine tprn2n(amun2n,apsi,innpr,in2pr,in2ntr,ipbtmx,
     $ jpfcmx,natmax,ncvn2n,nerr,nnnpr,nn2pr,nn2ntr,noutpt,npx3mx,
     $ npx3t,nttyo,nwarn,pcvn2n,qpdnn,qpdn2,qpdn2n,uaqsp,utripl)
c
c     Test and process the Pitzer data for nnn' (neutral, same neutral,
c     different neutral) triplets read from the DATA0 file. Find and
c     flag errors, such as duplication of data (e.g., two data blocks
c     for the same repeated neutral-distinct neutral triplet). The
c     conventional primitive Pitzer parameters are identical to the
c     observable parameters read from the data file:
c
c       mu(nnn') -> mu(nnn')
c
c     Note that there are two mu parameters mu(nnn') and mu(n'n'n)
c     for each corresponding lambda parameter lambda(nn'). Both are
c     handled in one set.
c
c     Check the coverage of entered Pitzer data against all possible
c     nnn' triplets that can be composed of the neutral species present
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
      integer ipbtmx,jpfcmx,natmax,nnnpr,nn2ntr,nn2pr,npx3mx
c
      integer noutpt,nttyo
c
      integer innpr(2,nnnpr),in2pr(nn2pr),in2ntr(2,nn2ntr)
c
      integer ncvn2n,nerr,npx3t,nwarn
c
      logical qpdnn(nnnpr),qpdn2(nn2pr),qpdn2n(nn2ntr)
c
      character(len=24) uaqsp(natmax),utripl(3,npx3mx)
c
      real*8 amun2n(jpfcmx,nn2ntr),apsi(jpfcmx,npx3mx)
c
      real*8 pcvn2n
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,ii,j,jj,jt,jtripl,j2,j3,j4,j5,j6,k,n,ncount,ndupl,
     $ nlistl,nn,nodatc
c
      integer ilnobl
c
      logical qfound1,qfound2,qfound3
c
      character*80 ustr80
      character*24 unam1,unam2,unam3
      character*8 ux8
c
c-----------------------------------------------------------------------
c
c     Limit on the list of triplets for which no data were found.
c
      data nlistl / 20 /
c
c-----------------------------------------------------------------------
c
c     Initialize the data arrays.
c
      do n = 1,nn2ntr
        qpdn2n(n) = .false.
        do j = 1,jpfcmx
          amun2n(j,n) = 0.
        enddo
      enddo
c
c     Check the entered data for nnn' triplets.
c
      nodatc = 0
c
      do n = 1,nn2ntr
        i = in2ntr(1,n)
        j = in2ntr(1,n)
        k = in2ntr(2,n)
        unam1 = uaqsp(i)
        unam2 = uaqsp(j)
        unam3 = uaqsp(k)
        j2 = ilnobl(unam1)
        j3 = ilnobl(unam2)
        j4 = ilnobl(unam3)
c
c       Search for unam1, unam2, unam3 in the utripl array.
c       That array corresponds to the species triplets blocks.
c
        call srch33(jtripl,unam1,unam2,unam3,utripl,npx3mx,npx3t)
c
        if (jtripl .le. 0) then
c
c
c         No data block was found on the DATA0 file.
c         Note that qpdn2n(n) is left with a value of .false.
c
          nodatc = nodatc + 1
        else
c
c         Have found an entry.
c
          qpdn2n(n) = .true.
c
c         Store the data.
c
          do j = 1,jpfcmx
            amun2n(j,n) = apsi(j,jtripl)
          enddo
c
c         Check for duplicate data sets.
c
          ndupl = 0
          do jt = jtripl + 1,npx3t
            if (unam1(1:24) .eq. utripl(1,jt)(1:24)) then
              if (unam2(1:24) .eq. utripl(2,jt)(1:24)) then
                if (unam3(1:24) .eq. utripl(3,jt)(1:24)) then
                  ndupl = ndupl + 1
                endif
              endif
            endif
          enddo
c
          if (ndupl .gt. 0) then
            if (ndupl .eq. 1) then
              write (noutpt,1010) unam1(1:j2),unam2(1:j3),unam3(1:j4)
              write (nttyo,1010) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1010         format(/' * Error - (EQPT/tprn2n) Have found a',
     $        ' duplicate data block on the DATA0 file',
     $        /7x,"for the nnn' triplet ",a,', ',a,', ',a,'.')
            else
              ux8 = ' '
              write (ux8,'(i5)') ndupl
              call lejust(ux8)
              j6 = ilnobl(ux8)
              write (noutpt,1020) ux8(1:j6),unam1(1:j2),unam2(1:j3),
     $        unam3(1:j4)
              write (nttyo,1020) ux8(1:j6),unam1(1:j2),unam2(1:j3),
     $        unam3(1:j4)
 1020         format(/' * Error - (EQPT/tprn2n) Have found ',a,
     $        ' duplicate data blocks on the DATA0 file',
     $        /7x,"for the nnn' triplet ",a,', ',a,', ',a,'.')
            endif
            nerr = nerr + ndupl
          endif
c
c         Check for the presence of the corresponding lamda(nn),
c         lambda(n'n'), and lambda(nn')  data. The presence of these
c         data on the data file is required.
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
  130     continue
c
          qfound2 = .false.
c
          do nn = 1,nn2pr
            ii = in2pr(nn)
            if (uaqsp(ii)(1:24) .eq. unam3(1:24)) then
              if (qpdn2(nn)) then
                qfound2 = .true.
                go to 140
              endif
            endif
          enddo
  140     continue
c
          qfound3 = .false.
c
          if (unam1 .lt. unam3) then
            do nn = 1,nnnpr
              ii = innpr(1,nn)
              if (uaqsp(ii)(1:24) .eq. unam1(1:24)) then
                jj = innpr(2,nn)
                if (uaqsp(jj)(1:24) .eq. unam3(1:24)) then
                  if (qpdnn(nn)) then
                    qfound3 = .true.
                    go to 150
                  endif
                endif
              endif
            enddo
  150       continue
          else
            do nn = 1,nnnpr
              ii = innpr(1,nn)
              if (uaqsp(ii)(1:24) .eq. unam3(1:24)) then
                jj = innpr(2,nn)
                if (uaqsp(jj)(1:24) .eq. unam1(1:24)) then
                  if (qpdnn(nn)) then
                    qfound3 = .true.
                    go to 160
                  endif
                endif
              endif
            enddo
  160       continue
          endif
c
          if (.not.(qfound1 .and. qfound2 .and. qfound3)) then
c
c           Certain required data were not found on the DATA0 file.
c
            write (noutpt,1030) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1030) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1030       format(/' * Error - (EQPT/tprn2n) Have data on the DATA0',
     $      " file for the nnn'",/7x,'triplet ',a,', ',a,', ',a,', but',
     $      " don't have the",/7x,'required data for the following',
     $      ' neutral-neutral pair(s):',/)
c
            if (.not.qfound1) then
              write (noutpt,1040) unam1(1:j2),unam1(1:j2)
              write (nttyo,1040) unam1(1:j2),unam1(1:j2)
 1040         format(9x,a,', ',a)
            endif
c
            if (.not.qfound2) then
              write (noutpt,1040) unam3(1:j4),unam3(1:j4)
              write (nttyo,1040) unam3(1:j4),unam3(1:j4)
            endif
c
            if (.not.qfound3) then
              write (noutpt,1040) unam1(1:j2),unam3(1:j4)
              write (nttyo,1040) unam1(1:j2),unam3(1:j4)
            endif
c
            nerr = nerr + 1
          endif
        endif
c
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nodatc .gt. 0) then
        write (noutpt,1050)
        write (nttyo,1050)
 1050   format(/' * Warning - (EQPT/tprn2n) Did not find a data',
     $  ' block on the DATA0 file',/7x,'for any of the following',
     $  " nnn' triplets:",/)
c
        ncount = 0
        do n = 1,nn2ntr
          if (.not.qpdn2n(n)) then
            ncount = ncount + 1
            i = in2ntr(1,n)
            j = in2ntr(1,n)
            k = in2ntr(2,n)
            unam1 = uaqsp(i)
            unam2 = uaqsp(j)
            unam3 = uaqsp(k)
            j2 = ilnobl(unam1)
            j3 = ilnobl(unam2)
            j4 = ilnobl(unam3)
            ustr80 = unam1(1:j2) // ', ' // unam2(1:j3) //
     $      ', ' // unam3(1:j4)
            j5 = ilnobl(ustr80)
            write (noutpt,1060) ustr80(1:j5)
            write (nttyo,1060) ustr80(1:j5)
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
      ncvn2n = nn2ntr - nodatc
      pcvn2n = (100.*ncvn2n)/float(nn2ntr)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
