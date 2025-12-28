      subroutine wrpz23(alphca,alamaa,alamca,alamcc,alamna,alamnc,
     $ alamnn,alamn2,amuaac,amua2c,amucca,amuc2a,amunca,amun2n,
     $ amun3,ipbtmx,iaapr,icapr,iccpr,inapr,incpr,innpr,in2pr,
     $ iaactr,ia2ctr,iccatr,ic2atr,incatr,in2ntr,in3tr,jpdblo,
     $ jpfcmx,natmax,naapr,ncapr,nccpr,nnapr,nncpr,nnnpr,nn2pr,
     $ naactr,na2ctr,nccatr,nc2atr,nncatr,nn2ntr,nn3tr,ndata1,
     $ ndat1f,noutpt,nttyo,uaqsp,uethfl)
c
c     This subroutine writes the processed Pitzer parameters (in
c     conventional lambda and mu form) on the DATA1 and DATA1F files.
c
c     This subroutine is called by:
c
c       EQPT/eqpt.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       alamaa  = array of coefficients for computing lamda parameters
c                   for aa' pairs
c       alamca  = array of coefficients for computing lamda parameters
c                   for ca pairs
c       alamcc  = array of coefficients for computing lamda parameters
c                   for cc' pairs
c       alamna  = array of coefficients for computing lamda parameters
c                   for na pairs
c       alamnc  = array of coefficients for computing lamda parameters
c                   for nc pairs
c       alamnn  = array of coefficients for computing lamda parameters
c                   for nn' pairs
c       alamn2  = array of coefficients for computing lamda parameters
c                   for nn pairs
c
c       amuaac  = array of coefficients for computing mu parameters
c                   for aa'c triplets
c       amua2c  = array of coefficients for computing mu parameters
c                   for aac triplets
c       amucca  = array of coefficients for computing mu parameters
c                   for cc'a triplets
c       amuc2a  = array of coefficients for computing mu parameters
c                   for cca triplets
c       amunca  = array of coefficients for computing mu parameters
c                   for nca triplets
c       amun2n  = array of coefficients for computing mu parameters
c                   for nnn' triplets
c       amun3   = array of coefficients for computing mu parameters
c                   for nnn triplets
c
c       iaapr   = species index array for aa' pairs
c       icapr   = species index array for ca pairs
c       iccpr   = species index array for cc' pairs
c       inapr   = species index array for na pairs
c       incpr   = species index array for nc pairs
c       innpr   = species index array for nn' pairs
c       in2pr   = species index array for nn pairs
c
c       iaactr  = species index array for aa'c triplets
c       ia2ctr  = species index array for aac triplets
c       iccatr  = species index array for cc'a triplets
c       ic2atr  = species index array for cca triplets
c       incatr  = species index array for nca triplets
c       in2ntr  = species index array for nn'n triplets
c
c       naapr   = number of aa' pairs
c       ncapr   = number of ca pairs
c       nccpr   = number of cc' pairs
c       nnapr   = number of na pairs
c       nncpr   = number of nc pairs
c       nnnpr   = number of nn' pairs
c       nn2pr   = number of nn pairs
c
c       naactr  = number of aa'c triplets
c       na2ctr  = number of aac triplets
c       nccatr  = number of cc'a triplets
c       nc2atr  = number of cca triplets
c       nncatr  = number of nca triplets
c       nn2ntr  = number of nn'n triplets
c
c       ndata1 = unit number of the DATA1 file
c       ndat1f = unit number of the DATA1F file
c
c       uethfl = E-theta flag string
c       upair  = array of names in species pairs
c       utripl = array of names in species triplets
c
c     Principal output:
c
c       none
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ipbtmx,jpfcmx,natmax
c
      integer ndata1,ndat1f,noutpt,nttyo
c
      integer jpdblo
      integer naapr,ncapr,nccpr,nnapr,nncpr,nnnpr,nn2pr
      integer naactr,na2ctr,nccatr,nc2atr,nncatr,nn2ntr,nn3tr
c
      integer iaapr(2,naapr),icapr(2,ncapr),iccpr(2,nccpr),
     $ inapr(2,nnapr),incpr(2,nncpr),innpr(2,nnnpr),in2pr(nn2pr)
      integer iaactr(3,naactr),ia2ctr(2,na2ctr),iccatr(3,nccatr),
     $ ic2atr(2,nc2atr),incatr(3,nncatr),in2ntr(2,nn2ntr),in3tr(nn3tr)
c
      character*24 uaqsp(natmax)
      character*8 uethfl
c
      real*8 alphca(ipbtmx,ncapr)
c
      real*8 alamaa(jpfcmx,0:ipbtmx,naapr),
     $ alamca(jpfcmx,0:ipbtmx,ncapr),alamcc(jpfcmx,0:ipbtmx,nccpr),
     $ alamna(jpfcmx,0:ipbtmx,nnapr),alamnc(jpfcmx,0:ipbtmx,nncpr),
     $ alamnn(jpfcmx,0:ipbtmx,nnnpr),alamn2(jpfcmx,0:ipbtmx,nn2pr)
c
      real*8 amuaac(jpfcmx,naactr),amua2c(jpfcmx,na2ctr),
     $ amucca(jpfcmx,nccatr),amuc2a(jpfcmx,nc2atr),
     $ amunca(jpfcmx,nncatr),amun2n(jpfcmx,nn2ntr),amun3(jpfcmx,nn3tr)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,j,jpfc,j2,j3,j4,k,n
c
      integer ilnobl
c
      character(len=16) ux16
c
      logical qnzdat
c
      character*80 ustr
      character*24 ublk24,ustp24,uend24,unam1,unam2,unam3
      character*72 uterm,utermc
c
      real*8 zero
c
c-----------------------------------------------------------------------
c
      data uend24 / 'endit.                  '  /
      data ublk24 / '                        '  /
      data ustp24 / 'stop.                   '  /
      data zero   /    0.0    /
c
c-----------------------------------------------------------------------
c
      uterm(1:48) = '+-----------------------------------------------'
      uterm(49:72) = '------------------------'
      utermc = uterm
      utermc(1:1) = '*'
c
      j3 = ilnobl(uend24)
      j4 = ilnobl(ustp24)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (jpdblo .eq. -1) then
c
c       Old Pitzer data block organization.
c
c       Write the E-lambda flag.
c
        ustr = 'E-lambda flag = ' // uethfl
        write (ndata1) ustr
        j2 = ilnobl(ustr)
        write (ndat1f,1010) ustr(1:j2)
 1010   format('*',/a,/'*')
      else
c
c       New Pitzer data block organization.
c
c       Do not write the E-lambda flag. It is implicitly "on" in the
c       case of the new Pitzer data block organization.
c
        continue
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the lambda set for ca pairs.
c
      write (ndat1f,1100)
 1100 format("* Coefficients for lamda(ca) parameters")
      write (ndat1f,1110) utermc
 1110 format(a72)
c
      do n = 1,ncapr
c
c       Are there any non-zero data?
c
        qnzdat = .false.
        do i = 0,ipbtmx
          do jpfc = 1,jpfcmx
            if (alamca(jpfc,i,n) .ne. 0.) then
              qnzdat = .true.
              go to 120
            endif
          enddo
        enddo
 120    continue
c
        if (qnzdat) then
c
c         Write a block.
c
          i = icapr(1,n)
          j = icapr(2,n)
          unam1 = uaqsp(i)
          unam2 = uaqsp(j)
c
          write (ndata1) unam1,unam2
          j2 = ilnobl(unam2)
          write (ndat1f,1120) unam1,unam2(1:j2)
 1120     format(a24,2x,a)
c
          if (jpdblo .eq. -1) then
c
c           Old Pitzer data block organization.
c
            write (ndata1) alamca(1,0,n),alamca(1,1,n),alamca(1,2,n)
            write (ndat1f,1130) alamca(1,0,n),alamca(1,1,n),
     $      alamca(1,2,n)
 1130       format(3x,'lambda0 = ',f9.5,3x,'lambda1 = ',f9.5,
     $      3x,'lambda2 = ',f9.5)
c
            write (ndata1) alphca(1,n),alphca(2,n)
            write (ndat1f,1140) alphca(1,n),alphca(2,n)
 1140       format(25x,'alpha1 = ',f5.1,7x,'alpha2 = ',f5.1)
c
            write (ndata1) alamca(2,0,n),alamca(3,0,n)
            write (ndat1f,1150) alamca(2,0,n),alamca(3,0,n)
 1150       format(4x,'dl0/dt = ',1pe10.3,2x,'d2l0/dt2 = ',1pe10.3)
c
            write (ndata1) alamca(2,1,n),alamca(3,1,n)
            write (ndat1f,1160) alamca(2,1,n),alamca(3,1,n)
 1160       format(4x,'dl1/dt = ',1pe10.3,2x,'d2l1/dt2 = ',1pe10.3)
c
            write (ndata1) alamca(2,2,n),alamca(3,2,n)
            write (ndat1f,1170) alamca(2,2,n),alamca(3,2,n)
 1170       format(4x,'dl2/dt = ',1pe10.3,2x,'d2l2/dt2 = ',1pe10.3)
          else
c
c           New Pitzer data block organization.
c
            ux16 = 'alpha( )'
            j2 = 8
            do i = 1,ipbtmx
              write (ux16(7:7),'(i1)') i
              write (ndata1) ux16,alphca(i,n)
              write (ndat1f,'(2x,a," = ",1pe12.5)')
     $        ux16(1:j2),alphca(i,n)
            enddo
c
            ux16 = 'lambda( )'
            j2 = 9
            do i = 0,ipbtmx
              write (ux16(8:8),'(i1)') i
              write (ndata1) ux16
              write (ndat1f,'(2x,a,":")') ux16(1:j2)
              do j = 1,jpfcmx
                write (ndata1) alamca(j,i,n)
                write (ndat1f,'(4x,"a",i1," = ",1pe15.8)')
     $          j,alamca(j,i,n)
              enddo
            enddo
          endif
c
c         Write the block terminator.
c
          write (ndata1) uterm
          write (ndat1f,1180) uterm
 1180     format(a72)
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the lambda set for nn pairs.
c
      write (ndat1f,1200)
 1200 format("* Coefficients for lamda(nn) parameters")
      write (ndat1f,1110) utermc
c
      do n = 1,nn2pr
c
c       Are there any non-zero data?
c
        qnzdat = .false.
        do i = 0,ipbtmx
          do jpfc = 1,jpfcmx
            if (alamn2(jpfc,i,n) .ne. 0.) then
              qnzdat = .true.
              go to 130
            endif
          enddo
        enddo
 130    continue
c
        if (qnzdat) then
c
c         Write a block.
c
          i = in2pr(n)
          unam1 = uaqsp(i)
          unam2 = unam1
c
          write (ndata1) unam1,unam2
          j2 = ilnobl(unam2)
          write (ndat1f,1120) unam1,unam2(1:j2)
c
          if (jpdblo .eq. -1) then
c
c           Old Pitzer data block organization.
c
            write (ndata1) alamn2(1,0,n),alamn2(1,1,n),alamn2(1,2,n)
            write (ndat1f,1130) alamn2(1,0,n),alamn2(1,1,n),
     $      alamn2(1,2,n)
c
c           No alpha parameters for this kind of pair.
c
            write (ndata1) zero,zero
            write (ndat1f,1140) zero, zero
c
            write (ndata1) alamn2(2,0,n),alamn2(3,0,n)
            write (ndat1f,1150) alamn2(2,0,n),alamn2(3,0,n)
c
            write (ndata1) alamn2(2,1,n),alamn2(3,1,n)
            write (ndat1f,1160) alamn2(2,1,n),alamn2(3,1,n)
c
            write (ndata1) alamn2(2,2,n),alamn2(3,2,n)
            write (ndat1f,1170) alamn2(2,2,n),alamn2(3,2,n)
          else
c
c           New Pitzer data block organization.
c
            ux16 = 'lambda(0)'
            j2 = 9
            write (ndata1) ux16
            write (ndat1f,'(2x,a,":")') ux16(1:j2)
            do j = 1,jpfcmx
              write (ndata1) alamn2(j,0,n)
              write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,alamn2(j,0,n)
            enddo
          endif
c
c         Write the block terminator.
c
          write (ndata1) uterm
          write (ndat1f,1180) uterm
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the lambda set for nc pairs.
c
      write (ndat1f,1210)
 1210 format("* Coefficients for lamda(nc) parameters")
      write (ndat1f,1110) utermc
c
      do n = 1,nncpr
c
c       Are there any non-zero data?
c
        qnzdat = .false.
        do i = 0,ipbtmx
          do jpfc = 1,jpfcmx
            if (alamnc(jpfc,i,n) .ne. 0.) then
              qnzdat = .true.
              go to 140
            endif
          enddo
        enddo
 140    continue
c
        if (qnzdat) then
c
c         Write a block.
c
          i = incpr(1,n)
          j = incpr(2,n)
          unam1 = uaqsp(i)
          unam2 = uaqsp(j)
c
          write (ndata1) unam1,unam2
          j2 = ilnobl(unam2)
          write (ndat1f,1120) unam1,unam2(1:j2)
c
          if (jpdblo .eq. -1) then
c
c           Old Pitzer data block organization.
c
            write (ndata1) alamnc(1,0,n),alamnc(1,1,n),alamnc(1,2,n)
            write (ndat1f,1130) alamnc(1,0,n),alamnc(1,1,n),
     $      alamnc(1,2,n)
c
c           No alpha parameters for this kind of pair.
c
            write (ndata1) zero,zero
            write (ndat1f,1140) zero, zero
c
            write (ndata1) alamnc(2,0,n),alamnc(3,0,n)
            write (ndat1f,1150) alamnc(2,0,n),alamnc(3,0,n)
c
            write (ndata1) alamnc(2,1,n),alamnc(3,1,n)
            write (ndat1f,1160) alamnc(2,1,n),alamnc(3,1,n)
c
            write (ndata1) alamnc(2,2,n),alamnc(3,2,n)
            write (ndat1f,1170) alamnc(2,2,n),alamnc(3,2,n)
          else
c
c           New Pitzer data block organization.
c
            ux16 = 'lambda(0)'
            j2 = 9
            write (ndata1) ux16
            write (ndat1f,'(2x,a,":")') ux16(1:j2)
            do j = 1,jpfcmx
              write (ndata1) alamnc(j,0,n)
              write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,alamnc(j,0,n)
            enddo
          endif
c
c         Write the block terminator.
c
          write (ndata1) uterm
          write (ndat1f,1180) uterm
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the lambda set for na pairs.
c
      write (ndat1f,1220)
 1220 format("* Coefficients for lamda(na) parameters")
      write (ndat1f,1110) utermc
c
      do n = 1,nnapr
c
c       Are there any non-zero data?
c
        qnzdat = .false.
        do i = 0,ipbtmx
          do jpfc = 1,jpfcmx
            if (alamna(jpfc,i,n) .ne. 0.) then
              qnzdat = .true.
              go to 150
            endif
          enddo
        enddo
 150    continue
c
        if (qnzdat) then
c
c         Write a block.
c
          i = inapr(1,n)
          j = inapr(2,n)
          unam1 = uaqsp(i)
          unam2 = uaqsp(j)
c
          write (ndata1) unam1,unam2
          j2 = ilnobl(unam2)
          write (ndat1f,1120) unam1,unam2(1:j2)
c
          if (jpdblo .eq. -1) then
c
c           Old Pitzer data block organization.
c
            write (ndata1) alamna(1,0,n),alamna(1,1,n),alamna(1,2,n)
            write (ndat1f,1130) alamna(1,0,n),alamna(1,1,n),
     $      alamna(1,2,n)
c
c           No alpha parameters for this kind of pair.
c
            write (ndata1) zero,zero
            write (ndat1f,1140) zero, zero
c
            write (ndata1) alamna(2,0,n),alamna(3,0,n)
            write (ndat1f,1150) alamna(2,0,n),alamna(3,0,n)
c
            write (ndata1) alamna(2,1,n),alamna(3,1,n)
            write (ndat1f,1160) alamna(2,1,n),alamna(3,1,n)
c
            write (ndata1) alamna(2,2,n),alamna(3,2,n)
            write (ndat1f,1170) alamna(2,2,n),alamna(3,2,n)
          else
c
c           New Pitzer data block organization.
c
            ux16 = 'lambda(0)'
            j2 = 9
            write (ndata1) ux16
            write (ndat1f,'(2x,a,":")') ux16(1:j2)
            do j = 1,jpfcmx
              write (ndata1) alamna(j,0,n)
              write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,alamna(j,0,n)
            enddo
          endif
c
c         Write the block terminator.
c
          write (ndata1) uterm
          write (ndat1f,1180) uterm
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the lambda set for nn' pairs.
c
      write (ndat1f,1230)
 1230 format("* Coefficients for lamda(nn') parameters")
      write (ndat1f,1110) utermc
c
      do n = 1,nnnpr
c
c       Are there any non-zero data?
c
        qnzdat = .false.
        do i = 0,ipbtmx
          do jpfc = 1,jpfcmx
            if (alamnn(jpfc,i,n) .ne. 0.) then
              qnzdat = .true.
              go to 160
            endif
          enddo
        enddo
 160    continue
c
        if (qnzdat) then
c
c         Write a block.
c
          i = innpr(1,n)
          j = innpr(2,n)
          unam1 = uaqsp(i)
          unam2 = uaqsp(j)
c
          write (ndata1) unam1,unam2
          j2 = ilnobl(unam2)
          write (ndat1f,1120) unam1,unam2(1:j2)
c
          if (jpdblo .eq. -1) then
c
c           Old Pitzer data block organization.
c
            write (ndata1) alamnn(1,0,n),alamnn(1,1,n),alamnn(1,2,n)
            write (ndat1f,1130) alamnn(1,0,n),alamnn(1,1,n),
     $      alamnn(1,2,n)
c
c           No alpha parameters for this kind of pair.
c
            write (ndata1) zero,zero
            write (ndat1f,1140) zero, zero
c
            write (ndata1) alamnn(2,0,n),alamnn(3,0,n)
            write (ndat1f,1150) alamnn(2,0,n),alamnn(3,0,n)
c
            write (ndata1) alamnn(2,1,n),alamnn(3,1,n)
            write (ndat1f,1160) alamnn(2,1,n),alamnn(3,1,n)
c
            write (ndata1) alamnn(2,2,n),alamnn(3,2,n)
            write (ndat1f,1170) alamnn(2,2,n),alamnn(3,2,n)
          else
c
c           New Pitzer data block organization.
c
            ux16 = 'lambda(0)'
            j2 = 9
            write (ndata1) ux16
            write (ndat1f,'(2x,a,":")') ux16(1:j2)
            do j = 1,jpfcmx
              write (ndata1) alamnn(j,0,n)
              write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,alamnn(j,0,n)
            enddo
          endif
c
c         Write the block terminator.
c
          write (ndata1) uterm
          write (ndat1f,1180) uterm
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (ndata1) uend24,ublk24,ublk24
      write (ndat1f,1270) uend24(1:j3)
 1270 format(a)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the lambda set for cc' pairs.
c
      write (ndat1f,1320)
 1320 format("* Coefficients for lamda(cc') parameters")
      write (ndat1f,1110) utermc
c
      do n = 1,nccpr
c
c       Are there any non-zero data?
c
        qnzdat = .false.
        do i = 0,ipbtmx
          do jpfc = 1,jpfcmx
            if (alamcc(jpfc,i,n) .ne. 0.) then
              qnzdat = .true.
              go to 170
            endif
          enddo
        enddo
 170    continue
c
        if (qnzdat) then
c
c         Write a block.
c
          i = iccpr(1,n)
          j = iccpr(2,n)
          unam1 = uaqsp(i)
          unam2 = uaqsp(j)
c
          write (ndata1) unam1,unam2
          j2 = ilnobl(unam2)
          write (ndat1f,1120) unam1,unam2(1:j2)
c
          if (jpdblo .eq. -1) then
c
c           Old Pitzer data block organization.
c
            write (ndata1) alamcc(1,0,n),alamcc(1,1,n),alamcc(1,2,n)
            write (ndat1f,1130) alamcc(1,0,n),alamcc(1,1,n),
     $      alamcc(1,2,n)
c
c           No alpha parameters for this kind of pair.
c
            write (ndata1) zero,zero
            write (ndat1f,1140) zero, zero
c
            write (ndata1) alamcc(2,0,n),alamcc(3,0,n)
            write (ndat1f,1150) alamcc(2,0,n),alamcc(3,0,n)
c
            write (ndata1) alamcc(2,1,n),alamcc(3,1,n)
            write (ndat1f,1160) alamcc(2,1,n),alamcc(3,1,n)
c
            write (ndata1) alamcc(2,2,n),alamcc(3,2,n)
            write (ndat1f,1170) alamcc(2,2,n),alamcc(3,2,n)
          else
c
c           New Pitzer data block organization.
c
            ux16 = 'lambda(0)'
            j2 = 9
            write (ndata1) ux16
            write (ndat1f,'(2x,a,":")') ux16(1:j2)
            do j = 1,jpfcmx
              write (ndata1) alamcc(j,0,n)
              write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,alamcc(j,0,n)
            enddo
          endif
c
c         Write the block terminator.
c
          write (ndata1) uterm
          write (ndat1f,1180) uterm
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the lambda set for aa' pairs.
c
      write (ndat1f,1330)
 1330 format("* Coefficients for lamda(aa') parameters")
      write (ndat1f,1110) utermc
c
      do n = 1,naapr
c
c       Are there any non-zero data?
c
        qnzdat = .false.
        do i = 0,ipbtmx
          do jpfc = 1,jpfcmx
            if (alamaa(jpfc,i,n) .ne. 0.) then
              qnzdat = .true.
              go to 180
            endif
          enddo
        enddo
 180    continue
c
        if (qnzdat) then
c
c         Write a block.
c
          i = iaapr(1,n)
          j = iaapr(2,n)
          unam1 = uaqsp(i)
          unam2 = uaqsp(j)
c
          write (ndata1) unam1,unam2
          j2 = ilnobl(unam2)
          write (ndat1f,1120) unam1,unam2(1:j2)
c
          if (jpdblo .eq. -1) then
c
c           Old Pitzer data block organization.
c
            write (ndata1) alamaa(1,0,n),alamaa(1,1,n),alamaa(1,2,n)
            write (ndat1f,1130) alamaa(1,0,n),alamaa(1,1,n),
     $      alamaa(1,2,n)
c
c           No alpha parameters for this kind of pair.
c
            write (ndata1) zero,zero
            write (ndat1f,1140) zero, zero
c
            write (ndata1) alamaa(2,0,n),alamaa(3,0,n)
            write (ndat1f,1150) alamaa(2,0,n),alamaa(3,0,n)
c
            write (ndata1) alamaa(2,1,n),alamaa(3,1,n)
            write (ndat1f,1160) alamaa(2,1,n),alamaa(3,1,n)
c
            write (ndata1) alamaa(2,2,n),alamaa(3,2,n)
            write (ndat1f,1170) alamaa(2,2,n),alamaa(3,2,n)
          else
c
c           New Pitzer data block organization.
c
            ux16 = 'lambda(0)'
            j2 = 9
            write (ndata1) ux16
            write (ndat1f,'(2x,a,":")') ux16(1:j2)
            do j = 1,jpfcmx
              write (ndata1) alamaa(j,0,n)
              write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,alamaa(j,0,n)
            enddo
          endif
c
c         Write the block terminator.
c
          write (ndata1) uterm
          write (ndat1f,1180) uterm
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (ndata1) uend24,ublk24,ublk24
      write (ndat1f,1270) uend24(1:j3)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the mu set for cca triplets.
c
      write (ndat1f,1340)
 1340 format("* Coefficients for mu(cca) parameters")
      write (ndat1f,1110) utermc
c
      do n = 1,nc2atr
c
c       Are there any non-zero data?
c
        qnzdat = .false.
        do jpfc = 1,jpfcmx
          if (amuc2a(jpfc,n) .ne. 0.) then
            qnzdat = .true.
            go to 220
          endif
        enddo
 220    continue
c
        if (qnzdat) then
c
c         Write a block.
c
          i = ic2atr(1,n)
          k = ic2atr(2,n)
          unam1 = uaqsp(i)
          unam2 = unam1
          unam3 = uaqsp(k)
c
          write (ndata1) unam1,unam2,unam3
          j2 = ilnobl(unam3)
          write (ndat1f,1370) unam1,unam2,unam3(1:j2)
 1370     format(a24,2x,a24,2x,a)
c
          if (jpdblo .eq. -1) then
c
c           Old Pitzer data block organization.
c
            write (ndata1) amuc2a(1,n),amuc2a(2,n),amuc2a(3,n)
            write (ndat1f,1380) amuc2a(1,n),amuc2a(2,n),amuc2a(3,n)
 1380       format(5x,'mummx = ',f9.5,2x,'dmummx/dt = ',1pe10.3,
     $      2x,'d2mummx/dt2 = ',1pe10.3)
          else
c
c           New Pitzer data block organization.
c
            ux16 = 'mu'
            j2 = 2
            write (ndata1) ux16
            write (ndat1f,'(2x,a,":")') ux16(1:j2)
            do j = 1,jpfcmx
              write (ndata1) amuc2a(j,n)
              write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,amuc2a(j,n)
            enddo
          endif
c
c         Write the block terminator.
c
          write (ndata1) uterm
          write (ndat1f,1180) uterm
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the mu set for aac triplets.
c
      write (ndat1f,1430)
 1430 format("* Coefficients for mu(aac) parameters")
      write (ndat1f,1110) utermc
c
      do n = 1,na2ctr
c
c       Are there any non-zero data?
c
        qnzdat = .false.
        do jpfc = 1,jpfcmx
          if (amua2c(jpfc,n) .ne. 0.) then
            qnzdat = .true.
            go to 230
          endif
        enddo
 230    continue
c
        if (qnzdat) then
c
c         Write a block.
c
          i = ia2ctr(1,n)
          k = ia2ctr(2,n)
          unam1 = uaqsp(i)
          unam2 = unam1
          unam3 = uaqsp(k)
c
          write (ndata1) unam1,unam2,unam3
          j2 = ilnobl(unam3)
          write (ndat1f,1370) unam1,unam2,unam3(1:j2)
c
          if (jpdblo .eq. -1) then
c
c           Old Pitzer data block organization.
c
            write (ndata1) amua2c(1,n),amua2c(2,n),amua2c(3,n)
            write (ndat1f,1380) amua2c(1,n),amua2c(2,n),amua2c(3,n)
          else
c
c           New Pitzer data block organization.
c
            ux16 = 'mu'
            j2 = 2
            write (ndata1) ux16
            write (ndat1f,'(2x,a,":")') ux16(1:j2)
            do j = 1,jpfcmx
              write (ndata1) amua2c(j,n)
              write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,amua2c(j,n)
            enddo
          endif
c
c         Write the block terminator.
c
          write (ndata1) uterm
          write (ndat1f,1180) uterm
       endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (ndata1) uend24,ublk24,ublk24
      write (ndat1f,1270) uend24(1:j3)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the mu set for cc'a triplets.
c
      write (ndat1f,1470)
 1470 format("* Coefficients for mu(cc'a) parameters")
      write (ndat1f,1110) utermc
c
      do n = 1,nccatr
c
c       Are there any non-zero data?
c
        qnzdat = .false.
        do jpfc = 1,jpfcmx
          if (amucca(jpfc,n) .ne. 0.) then
            qnzdat = .true.
            go to 240
          endif
        enddo
 240    continue
c
        if (qnzdat) then
c
c         Write a block.
c
          i = iccatr(1,n)
          j = iccatr(2,n)
          k = iccatr(3,n)
          unam1 = uaqsp(i)
          unam2 = uaqsp(j)
          unam3 = uaqsp(k)
c
          write (ndata1) unam1,unam2,unam3
          j2 = ilnobl(unam3)
          write (ndat1f,1370) unam1,unam2,unam3(1:j2)
c
          if (jpdblo .eq. -1) then
c
c           Old Pitzer data block organization.
c
            write (ndata1) amucca(1,n),amucca(2,n),amucca(3,n)
            write (ndat1f,1380) amucca(1,n),amucca(2,n),amucca(3,n)
          else
c
c           New Pitzer data block organization.
c
            ux16 = 'mu'
            j2 = 2
            write (ndata1) ux16
            write (ndat1f,'(2x,a,":")') ux16(1:j2)
            do j = 1,jpfcmx
              write (ndata1) amucca(j,n)
              write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,amucca(j,n)
            enddo
          endif
c
c         Write the block terminator.
c
          write (ndata1) uterm
          write (ndat1f,1180) uterm
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the mu set for aa'c triplets.
c
      write (ndat1f,1490)
 1490 format("* Coefficients for mu(aa'c) parameters")
      write (ndat1f,1110) utermc
c
      do n = 1,naactr
c
c       Are there any non-zero data?
c
        qnzdat = .false.
        do jpfc = 1,jpfcmx
          if (amuaac(jpfc,n) .ne. 0.) then
            qnzdat = .true.
            go to 260
          endif
        enddo
 260    continue
c
        if (qnzdat) then
c
c         Write a block.
c
          i = iaactr(1,n)
          j = iaactr(2,n)
          k = iaactr(3,n)
          unam1 = uaqsp(i)
          unam2 = uaqsp(j)
          unam3 = uaqsp(k)
c
          write (ndata1) unam1,unam2,unam3
          j2 = ilnobl(unam3)
          write (ndat1f,1370) unam1,unam2,unam3(1:j2)
c
          if (jpdblo .eq. -1) then
c
c           Old Pitzer data block organization.
c
            write (ndata1) amuaac(1,n),amuaac(2,n),amuaac(3,n)
            write (ndat1f,1380) amuaac(1,n),amuaac(2,n),amuaac(3,n)
          else
c
c           New Pitzer data block organization.
c
            ux16 = 'mu'
            j2 = 2
            write (ndata1) ux16
            write (ndat1f,'(2x,a,":")') ux16(1:j2)
            do j = 1,jpfcmx
              write (ndata1) amuaac(j,n)
              write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,amuaac(j,n)
            enddo
          endif
c
c         Write the block terminator.
c
          write (ndata1) uterm
          write (ndat1f,1180) uterm
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the mu set for nnn triplets.
c
      write (ndat1f,1510)
 1510 format("* Coefficients for mu(nnn) parameters")
      write (ndat1f,1110) utermc
c
      do n = 1,nn3tr
c
c       Are there any non-zero data?
c
        qnzdat = .false.
        do jpfc = 1,jpfcmx
          if (amun3(jpfc,n) .ne. 0.) then
            qnzdat = .true.
            go to 270
          endif
        enddo
 270    continue
c
        if (qnzdat) then
c
c         Write a block.
c
          i = in3tr(n)
          unam1 = uaqsp(i)
          unam2 = unam1
          unam3 = unam1
c
          write (ndata1) unam1,unam2,unam3
          j2 = ilnobl(unam3)
          write (ndat1f,1370) unam1,unam2,unam3(1:j2)
c
          if (jpdblo .eq. -1) then
c
c           Old Pitzer data block organization.
c
            write (ndata1) amun3(1,n),amun3(2,n),amun3(3,n)
            write (ndat1f,1380) amun3(1,n),amun3(2,n),amun3(3,n)
          else
c
c           New Pitzer data block organization.
c
            ux16 = 'mu'
            j2 = 2
            write (ndata1) ux16
            write (ndat1f,'(2x,a,":")') ux16(1:j2)
            do j = 1,jpfcmx
              write (ndata1) amun3(j,n)
              write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,amun3(j,n)
            enddo
          endif
c
c         Write the block terminator.
c
          write (ndata1) uterm
          write (ndat1f,1180) uterm
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the mu set for nnn' triplets.
c
      write (ndat1f,1520)
 1520 format("* Coefficients for mu(nnn') parameters")
      write (ndat1f,1110) utermc
c
      do n = 1,nn2ntr
c
c       Are there any non-zero data?
c
        qnzdat = .false.
        do jpfc = 1,jpfcmx
          if (amun2n(jpfc,n) .ne. 0.) then
            qnzdat = .true.
            go to 280
          endif
        enddo
 280    continue
c
        if (qnzdat) then
c
c         Write a block.
c
          i = in2ntr(1,n)
          k = in2ntr(2,n)
          unam1 = uaqsp(i)
          unam2 = unam1
          unam3 = uaqsp(k)
c
          write (ndata1) unam1,unam2,unam3
          j2 = ilnobl(unam3)
          write (ndat1f,1370) unam1,unam2,unam3(1:j2)
c
          if (jpdblo .eq. -1) then
c
c           Old Pitzer data block organization.
c
            write (ndata1) amun2n(1,n),amun2n(2,n),amun2n(3,n)
            write (ndat1f,1380) amun2n(1,n),amun2n(2,n),amun2n(3,n)
          else
c
c           New Pitzer data block organization.
c
            ux16 = 'mu'
            j2 = 2
            write (ndata1) ux16
            write (ndat1f,'(2x,a,":")') ux16(1:j2)
            do j = 1,jpfcmx
              write (ndata1) amun2n(j,n)
              write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,amun2n(j,n)
            enddo
          endif
c
c         Write the block terminator.
c
          write (ndata1) uterm
          write (ndat1f,1180) uterm
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Write the mu set for nca triplets.
c
      write (ndat1f,1530)
 1530 format("* Coefficients for mu(nca) parameters")
      write (ndat1f,1110) utermc
c
      do n = 1,nncatr
c
c       Are there any non-zero data?
c
        qnzdat = .false.
        do jpfc = 1,jpfcmx
          if (amunca(jpfc,n) .ne. 0.) then
            qnzdat = .true.
            go to 290
          endif
        enddo
 290    continue
c
        if (qnzdat) then
c
c         Write a block.
c
          i = incatr(1,n)
          j = incatr(2,n)
          k = incatr(3,n)
          unam1 = uaqsp(i)
          unam2 = uaqsp(j)
          unam3 = uaqsp(k)
c
          write (ndata1) unam1,unam2,unam3
          j2 = ilnobl(unam3)
          write (ndat1f,1370) unam1,unam2,unam3(1:j2)
c
          if (jpdblo .eq. -1) then
c
c           Old Pitzer data block organization.
c
            write (ndata1) amunca(1,n),amunca(2,n),amunca(3,n)
            write (ndat1f,1380) amunca(1,n),amunca(2,n),amunca(3,n)
          else
c
c           New Pitzer data block organization.
c
            ux16 = 'mu'
            j2 = 2
            write (ndata1) ux16
            write (ndat1f,'(2x,a,":")') ux16(1:j2)
            do j = 1,jpfcmx
              write (ndata1) amunca(j,n)
              write (ndat1f,'(4x,"a",i1," = ",1pe15.8)') j,amunca(j,n)
            enddo
          endif
c
c         Write the block terminator.
c
          write (ndata1) uterm
          write (ndat1f,1180) uterm
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (ndata1) uend24,ublk24,ublk24
      write (ndat1f,1270) uend24(1:j3)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (ndata1) ustp24,ublk24,ublk24
      write (ndat1f,1270) ustp24(1:j4)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
