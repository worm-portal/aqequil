      subroutine intnsp(coval,covali,ier,jflag,narn1,narn2,nbasp,nbt,
     $ nbti,nbtmax,nchlor,ncosp,ndecsp,nhydr,noutpt,nst,nstmax,nttyo,
     $ ucospi,uspec)
c
c     This subroutine finds the indices of species required to evaluate
c     certain types of constraints, such as a mean activity or a
c     heterogeneous equilibrium.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c     Principal output:
c
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer nbtmax,nstmax
c
      integer noutpt,nttyo
c
      integer jflag(nstmax),nbasp(nbtmax),ncosp(nbtmax),ndecsp(nbtmax)
c
      integer ier,narn1,narn2,nbt,nbti,nchlor,nhydr,nst
c
      character(len=48) ucospi(nbtmax),uspec(nstmax)
c
      real(8) coval(nbtmax),covali(nbtmax)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jfl,jlen,j2,n,nbi,nb1,nb2,ns1,ns2
c
      character(len=56) uspn56
      character(len=24) ublk24
c
c-----------------------------------------------------------------------
c
      data ublk24 /'                        '/
c
c-----------------------------------------------------------------------
c
      ier = 0
c
      do nbi = 1,nbti
        nb1 = ndecsp(nbi)
        if (nb1 .le. 0) go to 300
        ns1 = nbasp(nb1)
        jfl = jflag(ns1)
        coval(nb1) = covali(nbi)
        if (jfl.eq.17 .or. jfl.eq.18) then
c
c         Find the species index of the counterion for a mean log
c         activity or neutral log activity combination.
c
          do nb2 = 1,nbt
            ns2 = nbasp(nb2)
            if (ucospi(nbi)(1:24) .eq. uspec(ns2)(1:24)) then
              ncosp(nb1) = ns2
              go to 300
            endif
          enddo
c
          ier = ier + 1
c
c         Calling sequence substitutions:
c           ucospi(nbi) for unam48
c
          call fmspnx(jlen,ucospi(nbi),uspn56)
          if (jfl .eq. 17) then
            write (noutpt,1000) uspn56(1:jlen)
            write (nttyo,1000) uspn56(1:jlen)
 1000       format(/' * Error - (EQ3NR/intnsp) The species ',a,
     $      /7x,'is required for an abs(zj)*log ai + abs(zi)*log aj',
     $      /7x,"constraint (jflag= 17), but it isn't otherwise",
     $      ' present in the',/7x,'current system.')
          elseif (jfl .eq. 18) then
            write (noutpt,1010) uspn56(1:jlen)
            write (nttyo,1010) uspn56(1:jlen)
 1010       format(/' * Error - (EQ3NR/intnsp) The species ',a,
     $      /7x,'is required for a log a(+/-,ij) constraint',
     $      ' (jflag= 18),',/7x,"but it isn't otherwise present",
     $      ' in the current system.')
          endif
c
        elseif (jfl .eq. 21) then
c
c         Find the counterion for the pHCl constraint. Here the
c         constraint may be applied to either H+ or Cl-.
c
          if (ns1 .eq. nhydr) then
            ucospi(nbi) = 'Cl-'
            j2 = 3
            if (nchlor .gt. 0) then
              ncosp(nb1) = nchlor
            else
              ier = ier + 1
              write (noutpt,1020) ucospi(nbi)(1:j2)
              write (nttyo,1020) ucospi(nbi)(1:j2)
 1020         format(/' * Error - (EQ3NR/intnsp) The species ',a,
     $        /7x,'is required for a pHCl constraint (jflag= 21), but',
     $        /7x,"it isn't otherwise present in the current system.")
            endif
          elseif (ns1 .eq. nchlor) then
            ucospi(nbi) = 'H+'
            j2 = 2
            if (nhydr .gt. 0) then
              ncosp(nb1) = nhydr
            else
              ier = ier + 1
              write (noutpt,1020) ucospi(nbi)(1:j2)
              write (nttyo,1020) ucospi(nbi)(1:j2)
            endif
          else
            ier = ier + 1
c
c           Calling sequence substitutions:
c             uspec(ns1) for unam48
c
            call fmspnx(jlen,uspec(ns1),uspn56)
            write (noutpt,1040) uspn56(1:jlen)
            write (nttyo,1040) uspn56(1:jlen)
 1040       format(/' * Error - (EQ3NR/intnsp) The pHCl constraint',
     $      ' (jflag= 21)',/7x,"can't be applied to the species ",a,
     $      '.',/7x,'It can only be applied to H+ or Cl-.')
          endif
c
        elseif (jfl .eq. 25) then
c
c         Find the species index of the species whose heterogeneous
c         reaction is used to define an equilibrium constraint.
c
          n = narn1 - 1
          if (n .ge. 1) then
            do ns2 = 1,n
              if ( ucospi(nbi)(1:24) .eq. uspec(ns2)(1:24) .and.
     $           ( ucospi(nbi)(25:48) .eq. uspec(ns2)(25:48) .or.
     $             ucospi(nbi)(25:48) .eq. ublk24(1:24) ) ) then
                ncosp(nb1) = ns2
                go to 300
              endif
            enddo
          endif
c
          n = narn2 + 1
          if (nst .ge. n) then
            do ns2 = n,nst
              if ( ucospi(nbi)(1:24) .eq. uspec(ns2)(1:24) .and.
     $           ( ucospi(nbi)(25:48) .eq. uspec(ns2)(25:48) .or.
     $             ucospi(nbi)(25:48) .eq. ublk24(1:24) ) ) then
                ncosp(nb1) = ns2
                go to 300
              endif
            enddo
          endif
c
          ier = ier + 1
c
c         Calling sequence substitutions:
c           ucospi(nbi) for unam48
c
          call fmspnx(jlen,ucospi(nbi),uspn56)
          write (noutpt,1050) uspn56(1:jlen)
          write (nttyo,1050) uspn56(1:jlen)
 1050     format(/' * Error - (EQ3NR/intnsp) The species ',a,
     $    /7x,'is required for a heterogeneous equilibrium',
     $    ' constraint',/7x,"(jflag= 25), but it isn't otherwise",
     $    ' present in the current system.')
        endif
  300   continue
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
