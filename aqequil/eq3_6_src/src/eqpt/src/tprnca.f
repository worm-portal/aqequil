      subroutine tprnca(amunca,apsi,inapr,incatr,incpr,ipbtmx,
     $ jpfcmx,natmax,ncvnca,nerr,nnapr,nncatr,nncpr,noutpt,npx3mx,
     $ npx3t,nttyo,nwarn,pcvnca,qpdna,qpdnca,qpdnc,uaqsp,utripl)
c
c     Test and process the Pitzer data for nca (neutral, cation,
c     anion) triplets read from the DATA0 file. Find and flag errors,
c     such as duplication of data (e.g., two data blocks for the
c     same nca triplet). Calculate the conventional primitive Pitzer
c     parameters from the observable compound parameters read from
c     the data file:
c
c       zeta(nca) -> mu(nca)
c
c     Check the coverage of entered Pitzer data against all possible
c     nca triplets that can be composed of the neutral species, cations,
c     and anions present on the data file.
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
      integer ipbtmx,jpfcmx,natmax,nnapr,nncatr,nncpr,npx3mx
c
      integer noutpt,nttyo
c
      integer inapr(2,nnapr),incatr(3,nncatr),incpr(2,nncpr)
c
      integer ncvnca,nerr,npx3t,nwarn
c
      logical qpdna(nnapr),qpdnca(nncatr),qpdnc(nncpr)
c
      character(len=24) uaqsp(natmax),utripl(3,npx3mx)
c
      real*8 amunca(jpfcmx,nncatr),apsi(jpfcmx,npx3mx)
c
      real*8 pcvnca
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
      logical qfound1,qfound2
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
      do n = 1,nncatr
        qpdnca(n) = .false.
        do j = 1,jpfcmx
          amunca(j,n) = 0.
        enddo
      enddo
c
c     Check the entered data for nca triplets.
c
      nodatc = 0
c
      do n = 1,nncatr
        i = incatr(1,n)
        j = incatr(2,n)
        k = incatr(3,n)
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
c         No data block was found on the DATA0 file.
c         Note that qpdnca(n) is left with a value of .false.
c
          nodatc = nodatc + 1
        else
c
c         Have found an entry.
c
          qpdnca(n) = .true.
c
c         Store the data.
c         Note: here "cphi" is really zeta.
c
          do j = 1,jpfcmx
            amunca(j,n) = apsi(j,jtripl)/6.
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
 1010         format(/' * Error - (EQPT/tprnca) Have found a',
     $        ' duplicate data block on the DATA0 file',
     $        /7x,'for the nca triplet ',a,', ',a,', ',a,'.')
            else
              ux8 = ' '
              write (ux8,'(i5)') ndupl
              call lejust(ux8)
              j6 = ilnobl(ux8)
              write (noutpt,1020) ux8(1:j6),unam1(1:j2),unam2(1:j3),
     $        unam3(1:j4)
              write (nttyo,1020) ux8(1:j6),unam1(1:j2),unam2(1:j3),
     $        unam3(1:j4)
 1020         format(/' * Error - (EQPT/tprnca) Have found ',a,
     $        ' duplicate data blocks on the DATA0 file',
     $        /7x,'for the nca triplet ',a,', ',a,', ',a,'.')
            endif
            nerr = nerr + ndupl
          endif
c
c         Check for the presence of the corresponding lamda(nc) and
c         lambda(na) data. The presence of these data on the data
c         file is required.
c
          qfound1 = .false.
c
          do nn = 1,nncpr
            ii = incpr(1,nn)
            if (uaqsp(ii)(1:24) .eq. unam1(1:24)) then
              jj = incpr(2,nn)
              if (uaqsp(jj)(1:24) .eq. unam2(1:24)) then
                if (qpdnc(nn)) then
                  qfound1 = .true.
                  go to 130
                endif
              endif
            endif
          enddo
  130     continue
c
          qfound2 = .false.
c
          do nn = 1,nnapr
            ii = inapr(1,nn)
            if (uaqsp(ii)(1:24) .eq. unam1(1:24)) then
              jj = inapr(2,nn)
              if (uaqsp(jj)(1:24) .eq. unam3(1:24)) then
                if (qpdna(nn)) then
                  qfound2 = .true.
                  go to 140
                endif
              endif
            endif
          enddo
  140     continue
c
          if (.not.(qfound1 .and. qfound2)) then
c
c           Certain required data were not found on the DATA0 file.
c
            write (noutpt,1030) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1030) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1030       format(/' * Error - (EQPT/tprnca) Have data on the DATA0',
     $      ' file for the nca',/7x,'triplet ',a,', ',a,', ',a,', but',
     $      " don't have the",/7x,'required data for the following',
     $      ' neutral-ion pair(s):',/)
c
            if (.not.qfound1) then
              write (noutpt,1040) unam1(1:j2),unam2(1:j3)
              write (nttyo,1040) unam1(1:j2),unam2(1:j3)
 1040         format(9x,a,', ',a)
            endif
c
            if (.not.qfound2) then
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
 1050   format(/' * Warning - (EQPT/tprnca) Did not find a data',
     $  ' block on the DATA0 file',/7x,'for any of the following',
     $  ' nca triplets:',/)
c
        ncount = 0
        do n = 1,nncatr
          if (.not.qpdnca(n)) then
            ncount = ncount + 1
            i = incatr(1,n)
            j = incatr(2,n)
            k = incatr(3,n)
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
      ncvnca = nncatr - nodatc
      pcvnca = (100.*ncvnca)/float(nncatr)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
