      subroutine bldspc(iaapr,icapr,iccpr,inapr,incpr,innpr,
     $ in2pr,iaactr,ia2ctr,iccatr,ic2atr,incatr,in2ntr,in3tr,
     $ jassan,jassca,jassne,nat,natmax,naapr,ncapr,nccpr,nnapr,
     $ nncpr,nnnpr,nn2pr,naactr,na2ctr,nncatr,nccatr,nc2atr,
     $ nn2ntr,nn3tr,uaqsp,zaqsp)
c
c     Build index arrays for the pair and triplet combinations of
c     aqueous solute species for use with Pitzer parameters.
c
c     All pairs and triplets are distinct in the sense that the
c     ordering of the members is not relevant; e.g., ca = ac,
c     cc' = c'c, nca = anc = can = cna = acn = nac.
c
c     This subroutine is called by:
c
c       EQPT/eqpt.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       naapr  = number of aa' pairs
c       ncapr  = number of ca pairs
c       nccpr  = number of cc' pairs
c       nnapr  = number of na pairs
c       nncpr  = number of nc pairs
c       nnnpr  = number of nn' pairs
c       nn2pr  = number of nn pairs
c
c       naactr = number of aa'c triplets
c       na2ctr = number of aac triplets
c       nccatr = number of cc'a triplets
c       nc2atr = number of cca triplets
c       nncatr = number of nca triplets
c       nn3tr = number of nnn triplets
c       nn2ntr = number of nnn' triplets
c
c       nat    = the number of aqueous species
c       uaqsp  = array of names of aqueous species
c       zaqsp  = array of electrical charages of aqueous species
c
c     Principal output:
c
c       iaapr  = array of indices of species in aa' pairs
c       icapr  = array of indices of species in ca pairs
c       iccpr  = array of indices of species in cc' pairs
c       inapr  = array of indices of species in na pairs
c       incpr  = array of indices of species in nc pairs
c       innpr  = array of indices of species in nn' pairs
c       in2pr  = array of indices of species in nn pairs
c
c       iaactr = array of indices of species in aa'c triplets
c       ia2ctr = array of indices of species in aac triplets
c       iccatr = array of indices of species in cc'a triplets
c       ic2atr = array of indices of species in cca triplets
c       incatr = array of indices of species in nca triplets
c       in3tr = array of indices of species in nnn triplets
c       in2ntr = array of indices of species in nnn' triplets
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer natmax
c
      integer naapr,ncapr,nccpr,nnapr,nncpr,nnnpr,nn2pr,
     $ naactr,na2ctr,nccatr,nc2atr,nncatr,nn2ntr,nn3tr
c
      integer iaapr(2,naapr),icapr(2,ncapr),iccpr(2,nccpr),
     $ inapr(2,nnapr),incpr(2,nncpr),innpr(2,nnnpr),in2pr(nn2pr)
c
      integer iaactr(3,naactr),ia2ctr(2,na2ctr),iccatr(3,nccatr),
     $ ic2atr(2,nc2atr),incatr(3,nncatr),in2ntr(2,nn2ntr),in3tr(nn3tr)
c
      integer jassan,jassca,jassne,nat
c
      character(len=24) uaqsp(natmax)
c
      real*8 zaqsp(natmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer, dimension(:), allocatable :: iworka,iworkc,iworkn
c
      character(len=24), dimension(:), allocatable :: uworka,uworkc,
     $ uworkn
c
      integer i,j,k,jn,jc,ja,n,npr
c
c-----------------------------------------------------------------------
c
c     Allocate work arrays.
c
      ALLOCATE(iworkn(jassne))
      ALLOCATE(iworkc(jassca))
      ALLOCATE(iworka(jassan))
c
      ALLOCATE(uworkn(jassne))
      ALLOCATE(uworkc(jassca))
      ALLOCATE(uworka(jassan))
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Build scratch lists for cations, anions, and neutrals.
c     Do not include solvent water or O2(g) in the neutrals list.
c     Do not include e- on the anions list.
c
      jn = 0
      jc = 0
      ja = 0
      do n = 2,nat
        if (zaqsp(n) .eq. 0.) then
          if (uaqsp(n)(1:6) .ne. 'O2(g) ') then
            jn = jn + 1
            uworkn(jn) = uaqsp(n)
            iworkn(jn) = n
          endif
        elseif (zaqsp(n) .gt. 0.) then
          jc = jc + 1
          uworkc(jc) = uaqsp(n)
          iworkc(jc) = n
        else
          if (uaqsp(n)(1:3) .ne. 'e- ') then
            ja = ja + 1
            uworka(ja) = uaqsp(n)
            iworka(ja) = n
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Construct the index array for ca pairs.
c
      n = 0
      do i = 1,jc
        do j = 1,ja
          n = n + 1
          icapr(1,n) = iworkc(i)
          icapr(2,n) = iworka(j)
        enddo
      enddo
c
c     Construct the index array for cc' pairs (c' != c).
c
      n = 0
      do i = 1,jc
        do j = i + 1,jc
          n = n + 1
c
c         Store the pair in alphabetical order.
c
          if (uworkc(i) .lt. uworkc(j)) then
            iccpr(1,n) = iworkc(i)
            iccpr(2,n) = iworkc(j)
          else
            iccpr(1,n) = iworkc(j)
            iccpr(2,n) = iworkc(i)
          endif
        enddo
      enddo
c
c     Construct the index array for anion-distinct anion pairs.
c     Ignore pairs where an anion is paired with itself, such as
c     Cl-, Cl-, because the associated Pitzer coefficients are
c     either implicitly undefined or explicitly set to zero.
c     Ignore pairs where an anion is paired with itself.
c
      n = 0
      do i = 1,ja
        do j = i + 1,ja
          n = n + 1
c
c         Store the pair in alphabetical order.
c
          if (uworka(i) .lt. uworka(j)) then
            iaapr(1,n) = iworka(i)
            iaapr(2,n) = iworka(j)
          else
            iaapr(1,n) = iworka(j)
            iaapr(2,n) = iworka(i)
          endif
        enddo
      enddo
c
c     Construct the index array for nn pairs.
c
      n = 0
      do i = 1,jn
        n = n + 1
c
c       Store the "pair."
c
        in2pr(n) = iworkn(i)
      enddo
c
c     Construct the index array for nn' pairs.
c
      n = 0
      do i = 1,jn
        do j = i + 1,jn
          n = n + 1
c
c         Store the pair in alphabetical order.
c
          if (uworkn(i) .le. uworkn(j)) then
            innpr(1,n) = iworkn(i)
            innpr(2,n) = iworkn(j)
          else
            innpr(1,n) = iworkn(j)
            innpr(2,n) = iworkn(i)
          endif
        enddo
      enddo
c
c     Construct the index array for nc pairs.
c
      n = 0
      do i = 1,jn
        do j = 1,jc
          n = n + 1
          incpr(1,n) = iworkn(i)
          incpr(2,n) = iworkc(j)
        enddo
      enddo
c
c     Construct the index array for na pairs.
c
      n = 0
      do i = 1,jn
        do j = 1,ja
          n = n + 1
          inapr(1,n) = iworkn(i)
          inapr(2,n) = iworka(j)
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Construct the index array for cca triplets. Note that we will
c     generally just use the icapr array whenever the ic2atr array
c     is required.
c
      do n = 1,nc2atr
        ic2atr(1,n) = icapr(1,n)
        ic2atr(2,n) = icapr(2,n)
      enddo
c
c     Construct the index array for aac triplets. Note that we will
c     generally just use the icapr array whenever the ia2ctr array
c     is required.
c
      do n = 1,na2ctr
        ia2ctr(1,n) = icapr(2,n)
        ia2ctr(2,n) = icapr(1,n)
      enddo
c
c     Construct the index array for nnn triplets. Note that we will
c     generally just use the in2pr array whenever the ic2atr array
c     is required.
c
      do n = 1,nn3tr
        in3tr(n) = in2pr(n)
      enddo
c
c     Construct the index array for nnn' triplets. Note that we could
c     just use the innpr array (with some manipulation to account for
c     the fact that there are two triplets nnn' and n'n'n for each pair
c     nn') whenever the in2ntr array is required. See comments above
c     for the ic2atr array. Also note that nn2ntr = 2*nn2pr.
c
      n = 0
      do npr = 1,nn2pr
c
c       Store the two triplet combinations corresponding to a given
c       pair together.
c
        n = n + 1
        in2ntr(1,n) = innpr(1,npr)
        in2ntr(2,n) = innpr(2,npr)
        n = n + 1
        in2ntr(1,n) = innpr(2,npr)
        in2ntr(2,n) = innpr(1,npr)
      enddo
c
c     Construct the index array for cc'a triplets.
c
      n = 0
      do i = 1,jc
        do j = i + 1,jc
          do k = 1,ja
            n = n + 1
c
c           Store the two cations in alphabetical order.
c
            if (uworkc(i) .le. uworkc(j)) then
              iccatr(1,n) = iworkc(i)
              iccatr(2,n) = iworkc(j)
            else
              iccatr(1,n) = iworkc(j)
              iccatr(2,n) = iworkc(i)
            endif
            iccatr(3,n) = iworka(k)
          enddo
        enddo
      enddo
c
c     Construct the index array for aac triplets.
c
      n = 0
      do i = 1,ja
        do j = i + 1,ja
          do k = 1,jc
            n = n + 1
c
c           Store the two anions in alphabetical order.
c
            if (uworka(i) .le. uworka(j)) then
              iaactr(1,n) = iworka(i)
              iaactr(2,n) = iworka(j)
            else
              iaactr(1,n) = iworka(j)
              iaactr(2,n) = iworka(i)
            endif
            iaactr(3,n) = iworkc(k)
          enddo
        enddo
      enddo
c
c     Construct the index array for nca triplets.
c
      n = 0
      do i = 1,jn
        do j = 1,jc
          do k = 1,ja
            n = n + 1
            incatr(1,n) = iworkn(i)
            incatr(2,n) = iworkc(j)
            incatr(3,n) = iworka(k)
          enddo
        enddo
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Deallocate work arrays.
c
      DEALLOCATE(iworka)
      DEALLOCATE(iworkc)
      DEALLOCATE(iworkn)
c
      DEALLOCATE(uworka)
      DEALLOCATE(uworkc)
      DEALLOCATE(uworkn)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
