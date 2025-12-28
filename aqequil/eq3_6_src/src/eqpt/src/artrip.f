      subroutine artrip(iz1,iz2,iz3,na,nc,nn,n1,n2,n3,
     $ u1,u2,u3,z1,z2,z3)
c
c     This suboutine arranges the members of a species triplet
c     according to the following rules:
c
c       If exactly one neutral is present:
c         neutral, cation, anion
c
c       If exactly two neutrals are present:
c         neutral 1, neutral 1, neutral2
c
c       If two cations are present:
c         cation 1, cation 2, anion
c
c       If two anions are present:
c         anion 1, anion 2, cation
c
c     This suboutine is called by:
c
c       EQPT/rdpz3.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       iz1    = (integer) charge of the first species in a triplet
c       iz2    = (integer) charge of the second species in a triplet
c       iz3    = (integer) charge of the third species in a triplet
c       na     = number of anions in a species triplet
c       nc     = number of cations in a species triplet
c       nn     = number of neutral species in a species triplet
c       u1     = first name in a species triplet
c       u2     = second name in a species triplet
c       u3     = third name in a species triplet
c       z1     = charge of the first species in a triplet
c       z2     = charge of the second species in a triplet
c       z3     = charge of the third species in a triplet
c
c     Principal output:
c
c       none     (u1, u2, u3 are rearranged)
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer iz1,iz2,iz3,na,nc,nn,n1,n2,n3
c
      character(len=24) u1,u2,u3
c
      real(8) z1,z2,z3
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer nu,nx,izx
c
      character(len=24) ux24
c
      real(8) zx
c
c-----------------------------------------------------------------------
c
      if (nn .eq. 3) then
        if (n1.ne.n2 .or. n1.ne.n3) then
c
c         Have the combination nnn'.
c
c         Make sure that the unique neutral is in third place.
c
c         Find the unique neutral.
c
          if (n2 .eq. n3) nu = 1
          if (n1 .eq. n3) nu = 2
          if (n1 .eq. n2) nu = 3
c
          if (nu .eq. 1) then
c
c           The unique neutral is the first species. Switch it
c           with the third one.
c
            ux24 = u1
            u1 = u3
            u3 = ux24
c
            nx = n1
            n1 = n3
            n3 = nx
c
            zx = z1
            z1 = z3
            z3 = zx
c
            izx = iz1
            iz1 = iz3
            iz3 = izx
          elseif (nu .eq. 2) then
c
c           The unique neutral is the second species. Switch it
c           with the third one.
c
            ux24 = u2
            u2 = u3
            u3 = ux24
c
            nx = n2
            n2 = n3
            n3 = nx
c
            zx = z2
            z2 = z3
            z3 = zx
c
            izx = iz2
            iz2 = iz3
            iz3 = izx
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nn.eq.1 .and. nc.eq.1 .and. na.eq.1) then
c
c       Have the combination nca.
c
c       Make sure that the neutral is in first place.
c
        if (iz2 .eq. 0) then
c
          ux24 = u1
          u1 = u2
          u2 = ux24
c
          nx = n1
          n1 = n2
          n2 = nx
c
          zx = z1
          z1 = z2
          z2 = zx
c
          izx = iz1
          iz1 = iz2
          iz2 = izx
        elseif (iz3 .eq. 0) then
c
          ux24 = u1
          u1 = u3
          u3 = ux24
c
          nx = n1
          n1 = n3
          n3 = nx
c
          zx = z1
          z1 = z3
          z3 = zx
c
          izx = iz1
          iz1 = iz3
          iz3 = izx
        endif
c
c       Make sure that the cation is in second place. It
c       must presently be in second or third place.
c
        if (iz3 .gt. 0) then
c
          ux24 = u2
          u2 = u3
          u3 = ux24
c
          nx = n2
          n2 = n3
          n3 = nx
c
          zx = z2
          z2 = z3
          z3 = zx
c
          izx = iz2
          iz2 = iz3
          iz3 = izx
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nc.eq.2 .and. na.eq.1) then
c
c       Have the combination cc'a.
c
c       Make sure that the anion is in third place.
c
        if (iz1 .lt. 0) then
c
          ux24 = u1
          u1 = u3
          u3 = ux24
c
          nx = n1
          n1 = n3
          n3 = nx
c
          zx = z1
          z1 = z3
          z3 = zx
c
          izx = iz1
          iz1 = iz3
          iz3 = izx
        elseif (iz2 .lt. 0) then
c
          ux24 = u2
          u2 = u3
          u3 = ux24
c
          nx = n2
          n2 = n3
          n3 = nx
c
          zx = z2
          z2 = z3
          z3 = zx
c
          izx = iz2
          iz2 = iz3
          iz3 = izx
        endif
c
c       Put the two cations in alphabetical order.
c
        if (u1(1:24) .gt. u2(1:24)) then
c
          ux24 = u1
          u1 = u2
          u2 = ux24
c
          nx = n1
          n1 = n2
          n2 = nx
c
          zx = z1
          z1 = z2
          z2 = zx
c
          izx = iz1
          iz1 = iz2
          iz2 = izx
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (na.eq.2 .and. nc.eq.1) then
c
c       Have the combination aa'c.
c
C       Make sure that the cation is in third place.
c
        if (iz1 .gt. 0) then
c
          ux24 = u1
          u1 = u3
          u3 = ux24
c
          nx = n1
          n1 = n3
          n3 = nx
c
          zx = z1
          z1 = z3
          z3 = zx
c
          izx = iz1
          iz1 = iz3
          iz3 = izx
        elseif (iz2 .gt. 0) then
c
          ux24 = u2
          u2 = u3
          u3 = ux24
c
          nx = n2
          n2 = n3
          n3 = nx
c
          zx = z2
          z2 = z3
          z3 = zx
c
          izx = iz2
          iz2 = iz3
          iz3 = izx
        endif
c
c       Put the two anions in alphabetical order.
c
        if (u1(1:24) .gt. u2(1:24)) then
c
          ux24 = u1
          u1 = u2
          u2 = ux24
c
          nx = n1
          n1 = n2
          n2 = nx
c
          zx = z1
          z1 = z2
          z2 = zx
c
          izx = iz1
          iz1 = iz2
          iz2 = izx
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
