      subroutine tripck(na,nc,nerr,nn,noutpt,nttyo,n1,n2,n3,
     $ qdup12,qdup13,qdup23,unam1,unam2,unam3)
c
c     This suboutine checks the species triplets that were read
c     from the DATA0 file for illegal combinations.
c
c     The only legal combinations here correspond to ternary
c     systems and maximal third order in the expansion of the
c     Gibbs energy. These are:
c
c       nnn'      mu(nnn') maps to itself
c       nca       zeta(nca) maps to mu(nca)
c       cc'a      psi(cc'a) maps to mu(cc'a)
c       aa'c      psi(aa'c) maps to mu(aa'c)
c
c     Here n = neutral, n' a different neutral, c = cation, c' = a
c     different cation, a = anion, and a' = a different anion.
c     Solvent water (w) may not appear in any combination in the
c     normal Pitzer treatment of electrolyte solutions.
c
c     Other possible combinations are not allowed because of one or
c     more of the following reasons:
c
c       1. The combination corresponds to systems of higher order
c          (e.g., cc'c'')
c       2. The combination corresponds to systems of lesser order
c          (e.g., nnn and cca)
c       3. The combination corresponds to parameters that are defined
c          by convention to be zero (i.e., to unused parameters).
c          Examples include ccc and ncc.
c
c     This suboutine is called by:
c
c       EQPT/rdpz3.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       na     = number of anions in a species triplet
c       nc     = number of cations in a species triplet
c       nn     = number of neutral species in a species triplet
c       n1     = index of the first species in a triplet
c       n2     = index of the second species in a triplet
c       n3     = index of the third species in a triplet
c       unam1  = first name in a species triplet
c       unam2  = second name in a species triplet
c       unam3  = third name in a species triplet
c
c     Principal output:
c
c       nerr   = error counter
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer noutpt,nttyo
c
      integer na,nc,nn,nerr,n1,n2,n3
c
      logical qdup12,qdup13,qdup23
c
      character*24 unam1,unam2,unam3
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,j3,j4
c
      integer ilnobl
c
c-----------------------------------------------------------------------
c
      j2 = ilnobl(unam1)
      j3 = ilnobl(unam2)
      j4 = ilnobl(unam3)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for illegal appearance of solvent water.
c
      if (unam1(1:4) .eq.'H2O ' .or. unam2(1:4) .eq.'H2O ' .or.
     $  unam3(1:4).eq.'H2O ' .or. n1.eq.1 .or. n2.eq.1 .or.
     $  n3.eq.1) then
        write (noutpt,1430) unam1(1:j2),unam2(1:j3),unam3(1:j4)
        write (nttyo,1430) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1430   format(/' * Error - (EQPT/tripck) Have found an illegal data',
     $  ' block for the',/7x,'species triplet ',a,', ',a,', ',a,
     $  ' among the blocks for',/7x,'"mixture parameters".',
     $  ' Solvent water may not appear in such a',/7x,'triplet in',
     $  ' the normal Pitzer treatment of electrolyte',
     $  /7x,'solutions.')
        nerr = nerr + 1
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qdup12 .and. qdup13) then
c
c       The same species appears three times.
c
        if (nc.gt.0 .or. na.gt.0) then
c
c         The same cation or anion appears three times.(ccc or aaa).
c
          write (noutpt,1110) unam1(1:j2),unam2(1:j3),unam3(1:j4)
          write (nttyo,1110) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1110     format(/' * Error - (EQPT/tripck) Have found an illegal data',
     $    ' block for the',/7x,'species triplet ',a,', ',a,', ',a,
     $    ' among the blocks',/7x,'for "mixture parameters". The same',
     $    ' ion appears three times.',/7x,'This is not a valid',
     $    ' combination, as the corrresponding lambda and mu',
     $    /7x,'data are defined to be zero by convention in the normal',
     $    ' Pitzer',/7x,'treatment of electrolyte solutions.')
        endif
c
        if (nn .gt. 0) then
c
c         The same neutral appears three times (nnn).
c
          write (noutpt,1112) unam1(1:j2),unam2(1:j3),unam3(1:j4)
          write (nttyo,1112) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1112     format(/' * Error - (EQPT/tripck) Have found an illegal data',
     $    ' block for the',/7x,'species triplet ',a,', ',a,', ',a,
     $    ' among the blocks',/7x,'for "mixture parameters". The',
     $    ' same neutral appears three times.',/7x,'This is not',
     $    ' a valid combination for the present kind of block.',
     $    /7x,'Enter these data in a block for "single-salt',
     $    ' parameters".')
        endif
c
        nerr = nerr + 1
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for other illegal combinations involving two or more
c     neutrals.
c
      if (nn .eq. 3) then
        if (.not.(qdup12 .or. qdup13 .or. qdup23)) then
c
c         Have three distinct neutrals (nn'n'').
c
          write (noutpt,1130) unam1(1:j2),unam2(1:j3),unam3(1:j4)
          write (nttyo,1130) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1130     format(/' * Error - (EQPT/tripck) Have found an illegal data',
     $    ' block for the',/7x,'species triplet ',a,', ',a,', ',a,
     $    ' among the blocks',/7x,'for "mixture parameters". Three',
     $    ' distinct neutrals appear. This is not',/7x,'a valid',
     $    ' combination in the context of the normal Pitzer',
     $    /7x,'treatment of electrolyte solutions as it corresponds',
     $    ' to a',/7x,'quaternary system.')
          nerr = nerr + 1
          go to 999
        endif
      endif
c
      if (nn .eq. 2) then
        if (nc.eq.1 .or. na.eq.1) then
c
          if (qdup12 .or. qdup23 .or. qdup13) then
c
c           Have one of the following combinations: nnc or nna.
c
            write (noutpt,1132) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1132) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1132       format(/' * Error - (EQPT/tripck) Have found an illegal',
     $      ' data block for the',/7x,'species triplet ',a,', ',a,', ',
     $      a,' among the blocks',/7x,'for "mixture parameters". One',
     $      ' neutral appears twice with an ion.',/7x,'This is not a',
     $      ' valid combination in the context of the normal Pitzer',
     $      /7x,'treatment of electrolyte solutions. The corresponding',
     $      ' parameters',/7x,'are not used.')
          else
c
c           Have one of the following combinations: nn'c or nn'a.
c
            write (noutpt,1136) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1136) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1136       format(/' * Error - (EQPT/tripck) Have found an illegal',
     $      ' data block for the',/7x,'species triplet ',a,', ',a,', ',
     $      a,' among the blocks',/7x,'for "mixture parameters". Two',
     $      ' neutrals appear with an ion.',/7x,'This is not a',
     $      ' valid combination in the context of the normal Pitzer',
     $      /7x,'treatment of electrolyte solutions as it corresponds',
     $      ' to a',/7x,'quaternary system.')
          endif
          nerr = nerr + 1
          go to 999
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for other illegal combinations involving two or more
c     cations.
c
      if (nc .eq. 3) then
        if (.not.(qdup12 .or. qdup13 .or. qdup23)) then
c
c         Have three distinct cations (cc'c'').
c
          write (noutpt,1140) unam1(1:j2),unam2(1:j3),unam3(1:j4)
          write (nttyo,1140) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1140     format(/' * Error - (EQPT/tripck) Have found an illegal data',
     $    ' block for the',/7x,'species triplet ',a,', ',a,', ',a,
     $    ' among the blocks',/7x,'for "mixture parameters". Three',
     $    ' distinct cations appear. This is',/7x,'not a valid',
     $    ' combination in the context of the normal Pitzer',
     $    /7x,'treatment of electrolyte solutions as it corresponds',
     $    ' to a',/7x,'quaternary system.')
          nerr = nerr + 1
          go to 999
        endif
c
        if (qdup12 .and. qdup13) then
c
c         Have one cation appearing three times (ccc).
c
          write (noutpt,1142) unam1(1:j2),unam2(1:j3),unam3(1:j4)
          write (nttyo,1142) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1142     format(/' * Error - (EQPT/tripck) Have found an illegal data',
     $    ' block for the',/7x,'species triplet ',a,', ',a,', ',a,
     $    ' among the blocks',/7x,'for "mixture parameters". A single',
     $    ' cation appears three times. This',/7x,'is not a valid',
     $    ' combination in the context of the normal Pitzer',
     $    /7x,'treatment of electrolyte solutions. The corresponding',
     $    ' parameters are',/7x,'zero by convention.')
          nerr = nerr + 1
          go to 999
        endif
c
        if (qdup12 .or. qdup13 .or. qdup23) then
c
c         Have one cation appearing two times (ccc').
c
          write (noutpt,1144) unam1(1:j2),unam2(1:j3),unam3(1:j4)
          write (nttyo,1144) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1144     format(/' * Error - (EQPT/tripck) Have found an illegal data',
     $    ' block for the',/7x,'species triplet ',a,', ',a,', ',a,
     $    ' among the blocks',/7x,'for "mixture parameters". A single',
     $    ' cation appears twice with another',/7x,'cation. This is',
     $    ' not a valid combination in the context of the normal',
     $    /7x,'Pitzer treatment of electrolyte solutions. The',
     $    ' corresponding',/7x,'parameters are not used.')
          nerr = nerr + 1
          go to 999
        endif
      endif
c
      if (nc .eq. 2) then
c
        if (qdup12 .or. qdup23 .or. qdup13) then
c
          if (nn .eq. 1) then
c
c           Have a repeated cation and a neutral (ccn).
c
            write (noutpt,1150) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1150) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1150       format(/' * Error - (EQPT/tripck) Have found an illegal',
     $      ' data block for the',/7x,'species triplet ',a,', ',a,', ',
     $      a,' among the blocks',/7x,'for "mixture parameters". One',
     $      ' cation appears twice with a neutral.',/7x,'This is not a',
     $      ' valid combination in the context of the normal Pitzer',
     $      /7x,'treatment of electrolyte solutions. The',
     $      ' corresponding parameters',/7x,'are not used.')
          elseif (na .eq. 1) then
c
c           Have a repeated cation and an anion (cca).
c
            write (noutpt,1152) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1152) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1152       format(/' * Error - (EQPT/tripck) Have found an illegal',
     $      ' data block for the',/7x,'species triplet ',a,', ',a,', ',
     $      a,' among the blocks',/7x,'for "mixture parameters". One',
     $      ' cation appears twice with an anion.',/7x,'This is not a',
     $      ' valid combination for the present type of block. The',
     $      /7x,'corresponding parameters properly appear in a block',
     $      ' under',/7x,'"single-salt parameters".')
          endif
          nerr = nerr + 1
          go to 999
        else
c
          if (nn .eq. 1) then
c
c           Have two cations and a neutral (c'cn).
c
            write (noutpt,1154) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1154) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1154       format(/' * Error - (EQPT/tripck) Have found an illegal',
     $      ' data block for the',/7x,'species triplet ',a,', ',a,', ',
     $      a,' among the blocks',/7x,'for "mixture parameters". Two',
     $      ' cations appear with a neutral.',/7x,'This is not a',
     $      ' valid combination in the context of the normal Pitzer',
     $      /7x,'treatment of electrolyte solutions as it corresponds',
     $      ' to a',/7x,'quaternary system.')
            nerr = nerr + 1
            go to 999
          endif
c
c         Have two cations and an anion (cc'a). That is a legal
c         combination.
c
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for other illegal combinations involving two or more
c     anions.
c
      if (na .eq. 3) then
        if (.not.(qdup12 .or. qdup13 .or. qdup23)) then
c
c         Have three distinct anions (aa'a'').
c
          write (noutpt,1160) unam1(1:j2),unam2(1:j3),unam3(1:j4)
          write (nttyo,1160) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1160     format(/' * Error - (EQPT/tripck) Have found an illegal data',
     $    ' block for the',/7x,'species triplet ',a,', ',a,', ',a,
     $    ' among the blocks',/7x,'for "mixture parameters". Three',
     $    ' distinct anions appear. This is',/7x,'not a valid',
     $    ' combination in the context of the normal Pitzer',
     $    /7x,'treatment of electrolyte solutions as it corresponds',
     $    ' to a',/7x,'quaternary system.')
          nerr = nerr + 1
          go to 999
        endif
c
        if (qdup12 .and. qdup13) then
c
c         Have one anion appearing three times (aaa).
c
          write (noutpt,1162) unam1(1:j2),unam2(1:j3),unam3(1:j4)
          write (nttyo,1162) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1162     format(/' * Error - (EQPT/tripck) Have found an illegal data',
     $    ' block for the',/7x,'species triplet ',a,', ',a,', ',a,
     $    ' among the blocks',/7x,'for "mixture parameters". A single',
     $    ' anion appears three times. This',/7x,'is not a valid',
     $    ' combination in the context of the normal Pitzer',
     $    /7x,'treatment of electrolyte solutions. The corresponding',
     $    ' parameters are',/7x,'zero by convention.')
          nerr = nerr + 1
          go to 999
        endif
c
        if (qdup12 .or. qdup13 .or. qdup23) then
c
c         Have one anion appearing two times (aaa').
c
          write (noutpt,1164) unam1(1:j2),unam2(1:j3),unam3(1:j4)
          write (nttyo,1164) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1164     format(/' * Error - (EQPT/tripck) Have found an illegal data',
     $    ' block for the',/7x,'species triplet ',a,', ',a,', ',a,
     $    ' among the blocks',/7x,'for "mixture parameters". A single',
     $    ' anion appears twice with another',/7x,'anion. This is',
     $    ' not a valid combination in the context of the normal',
     $    /7x,'Pitzer treatment of electrolyte solutions. The',
     $    ' corresponding',/7x,'parameters are not used.')
          nerr = nerr + 1
          go to 999
        endif
      endif
c
      if (na .eq. 2) then
c
        if (qdup12 .or. qdup23 .or. qdup13) then
c
          if (nn .eq. 1) then
c
c           Have a repeated anion and a neutral (aan).
c
            write (noutpt,1170) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1170) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1170       format(/' * Error - (EQPT/tripck) Have found an illegal',
     $      ' data block for the',/7x,'species triplet ',a,', ',a,', ',
     $      a,' among the blocks',/7x,'for "mixture parameters". One',
     $      ' anion appears twice with a neutral.',/7x,'This is not a',
     $      ' valid combination in the context of the normal Pitzer',
     $      /7x,'treatment of electrolyte solutions. The',
     $      ' corresponding parameters',/7x,'are not used.')
          elseif (nc .eq. 1) then
c
c           Have a repeated anion and a cation (aac).
c
            write (noutpt,1173) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1173) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1173       format(/' * Error - (EQPT/tripck) Have found an illegal',
     $      ' data block for the',/7x,'species triplet ',a,', ',a,', ',
     $      a,' among the blocks',/7x,'for "mixture parameters". One',
     $      ' anion appears twice with a cation.',/7x,'This is not a',
     $      ' valid combination for the present type of block. The',
     $      /7x,'corresponding parameters properly appear in a block',
     $      ' under',/7x,'"single-salt parameters".')
          endif
          nerr = nerr + 1
          go to 999
        else
c
          if (nn .eq. 1) then
c
c           Have two anions and a neutral (a'an).
c
            write (noutpt,1174) unam1(1:j2),unam2(1:j3),unam3(1:j4)
            write (nttyo,1174) unam1(1:j2),unam2(1:j3),unam3(1:j4)
 1174       format(/' * Error - (EQPT/tripck) Have found an illegal',
     $      ' data block for the',/7x,'species triplet ',a,', ',a,', ',
     $      a,' among the blocks',/7x,'for "mixture parameters". Two',
     $      ' anions appear with a neutral.',/7x,'This is not a',
     $      ' valid combination in the context of the normal Pitzer',
     $      /7x,'treatment of electrolyte solutions as it corresponds',
     $      ' to a',/7x,'quaternary system.')
            nerr = nerr + 1
            go to 999
          endif
c
c         Have two anions and a cation (aa'c). That is a legal
c         combination.
c
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
c
      end
