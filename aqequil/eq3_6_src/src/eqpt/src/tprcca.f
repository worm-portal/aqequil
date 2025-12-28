      subroutine tprcca(amucca,amuc2a,apsi,icapr,iccatr,ipbtmx,
     $ jpfcmx,natmax,ncapr,nccatr,ncvcca,nc2atr,nerr,noutpt,
     $ npx3mx,npx3t,nttyo,nwarn,pcvcca,qpdca,qpdcca,uaqsp,
     $ utripl,zaqsp)
c
c     Test and process the Pitzer data for cc'a (cation, different
c     cation, anion) triplets read from the DATA0 file. Find and flag
c     errors, such as duplication of data (e.g., two data blocks for
c     the same cc'a triplet). Calculate the conventional primitive
c     Pitzer parameters from the observable compound parameters read
c     from the data file:
c
c       psi(cc'a) -> mu(cc'a)
c
c     Note: mu(cc'a) depends on cphi(ca) and cphi(c'a) as well as
c     psi(cc'a).
c
c     Check the coverage of entered Pitzer data against all possible
c     cc'a triplets that can be composed of the cations and anions
c     present on the data file.
c
c     This subroutine is complementary to tpraac.f, which does the
c     same thing for aa'c (anion, different anion, cation) triplets.
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
      integer ipbtmx,jpfcmx,natmax,ncapr,nccatr,nc2atr,npx3mx
c
      integer noutpt,nttyo
c
      integer icapr(2,ncapr),iccatr(3,nccatr)
c
      integer ncvcca,nerr,npx3t,nwarn
c
      logical qpdca(ncapr),qpdcca(nccatr)
c
      character(len=24) uaqsp(natmax),utripl(3,npx3mx)
c
      real*8 apsi(jpfcmx,npx3mx),amucca(jpfcmx,nccatr),
     $ amuc2a(jpfcmx,nc2atr),zaqsp(natmax)
c
      real*8 pcvcca
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
      character(len=80) ustr80
      character(len=24) unam1,unam2,unam3
      character(len=8) ux8
c
      real(8), dimension(:), allocatable :: amu1,amu2
c
      real(8) z1,z2
c
c-----------------------------------------------------------------------
c
c     Limit on the list of triplets for which no data were found.
c
      data nlistl / 20 /
c
c-----------------------------------------------------------------------
c
c     Allocate work space arrays.
c
      ALLOCATE(amu1(jpfcmx))
      ALLOCATE(amu2(jpfcmx))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize the data arrays.
c
      do n = 1,nccatr
        qpdcca(n) = .false.
        do j = 1,jpfcmx
          amucca(j,n) = 0.
        enddo
      enddo
c
c     Check the entered data for cc'a triplets.
c
      nodatc = 0
c
      do n = 1,nccatr
        i = iccatr(1,n)
        j = iccatr(2,n)
        k = iccatr(3,n)
        unam1 = uaqsp(i)
        unam2 = uaqsp(j)
        unam3 = uaqsp(k)
        j2 = ilnobl(unam1)
        j3 = ilnobl(unam2)
        j4 = ilnobl(unam3)
        z1 = zaqsp(i)
        z2 = zaqsp(j)
c
c       Search for unam1, unam2, unam3 in the utripl array.
c       That array corresponds to the species triplets blocks.
c
        call srch33(jtripl,unam1,unam2,unam3,utripl,npx3mx,npx3t)
c
        if (jtripl .le. 0) then
c
c         No data block was found on the DATA0 file.
c         Note that qpdcca(n) is left with a value of .false.
c
          nodatc = nodatc + 1
c
        else
c
c         Have found an entry.
c
          qpdcca(n) = .true.
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
 1010         format(/' * Error - (EQPT/tprcca) Have found a',
     $        ' duplicate data block',/7x,'on the DATA0 file',
     $        " for the cc'a triplet",/7x,a,', ',a,', ',a,'.')
            else
              ux8 = ' '
              write (ux8,'(i5)') ndupl
              call lejust(ux8)
              j6 = ilnobl(ux8)
              write (noutpt,1020) ux8(1:j6),unam1(1:j2),unam2(1:j3),
     $        unam3(1:j4)
              write (nttyo,1020) ux8(1:j6),unam1(1:j2),unam2(1:j3),
     $        unam3(1:j4)
 1020         format(/' * Error - (EQPT/tprcca) Have found ',a,
     $        ' duplicate data blocks',/7x,'on the DATA0 file',
     $        " for the cc'a triplet",/7x,a,', ',a,', ',a,'.')
            endif
            nerr = nerr + ndupl
          endif
c
c         Get the data for the two constituent binary systems.
c
          qfound1 = .false.
          do j = 1,jpfcmx
            amu1(j) = 0.
          enddo
c
          do nn = 1,ncapr
            ii = icapr(1,nn)
            if (uaqsp(ii)(1:24) .eq. unam1(1:24)) then
              jj = icapr(2,nn)
              if (uaqsp(jj)(1:24) .eq. unam3(1:24)) then
                if (qpdca(nn)) then
                  qfound1 = .true.
                  do j = 1,jpfcmx
                    amu1(j) = amuc2a(j,nn)
                  enddo
                  go to 130
                endif
              endif
            endif
          enddo
  130     continue
c
          qfound2 = .false.
          do j = 1,jpfcmx
            amu2(j) = 0.
          enddo
c
          do nn = 1,ncapr
            ii = icapr(1,nn)
            if (uaqsp(ii)(1:24) .eq. unam2(1:24)) then
              jj = icapr(2,nn)
              if (uaqsp(jj)(1:24) .eq. unam3(1:24)) then
                if (qpdca(nn)) then
                  qfound2 = .true.
                  do j = 1,jpfcmx
                    amu2(j) = amuc2a(j,nn)
                  enddo
                  go to 140
                endif
              endif
            endif
          enddo
  140     continue
c
          if (qfound1 .and. qfound2) then
            do j = 1,jpfcmx
              amucca(j,n) = (apsi(j,jtripl) +
     $        3.0*((z2/z1)*amu1(j) + (z1/z2)*amu2(j)))/6.
            enddo
          else
c
c           Certain required data were not found on the DATA0 file.
c
            write (noutpt,1030) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1030) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1030       format(/' * Error - (EQPT/tprcca) Have data on the DATA0',
     $      " file for the cc'a",/7x,'triplet ',a,', ',a,', ',a,', but',
     $      " don't have the",/7x,'required data for the following',
     $      ' cation-anion pair(s):',/)
c
            if (.not.qfound1) then
              write (noutpt,1040) unam1(1:j2),unam3(1:j4)
              write (nttyo,1040) unam1(1:j2),unam3(1:j4)
 1040         format(9x,a,', ',a)
            endif
c
            if (.not.qfound2) then
              write (noutpt,1040) unam2(1:j3),unam3(1:j4)
              write (nttyo,1040) unam2(1:j3),unam3(1:j4)
            endif
c
            nerr = nerr + 1
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nodatc .gt. 0) then
        write (noutpt,1060)
        write (nttyo,1060)
 1060   format(/' * Warning - (EQPT/tprcca) Did not find a data',
     $  ' block on the DATA0 file',/7x,'for any of the following',
     $  " cc'a triplets:",/)
c
        ncount = 0
        do n = 1,nccatr
          if (.not.qpdcca(n)) then
            ncount = ncount + 1
            i = iccatr(1,n)
            j = iccatr(2,n)
            k = iccatr(3,n)
            unam1 = uaqsp(i)
            unam2 = uaqsp(j)
            unam3 = uaqsp(k)
            j2 = ilnobl(unam1)
            j3 = ilnobl(unam2)
            j4 = ilnobl(unam3)
            ustr80 = unam1(1:j2) // ', ' // unam2(1:j3) // ', ' //
     $      unam3(1:j4)
            j5 = ilnobl(ustr80)
            write (noutpt,1070) ustr80(1:j5)
            write (nttyo,1070) ustr80(1:j5)
 1070       format(9x,a)
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
          write (noutpt,1080) ux8(1:j3)
          write (nttyo,1080) ux8(1:j3)
 1080     format(/9x,'plus ',a,' others')
        endif
        write (noutpt,1090)
        write (nttyo,1090)
 1090   format(1x)
        nwarn = nwarn + 1
      endif
c
      ncvcca = nccatr - nodatc
      pcvcca = (100.*ncvcca)/float(nccatr)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Deallocate work space arrays.
c
      DEALLOCATE(amu1)
      DEALLOCATE(amu2)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
