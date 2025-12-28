      subroutine ldlxrc(al10,delxi,dlxmin,dzvc0,iodb,iindx1,kbt,kdim,
     $ kelect,khydr,khydx,km1,kmax,ko2gaq,krdxsp,kwater,kxt,nbasp,
     $ nbtmax,nodbmx,nord,noutpt,nrd1mx,nstmax,nttyo,qrapch,uspec,
     $ zklogu,zvclg0,zvclg1,zvec0,zvec1)
c
c     This subroutine limits delxi when a variable corresponding to a
c     basis species is rapidly changing. Special forms of this
c     constraint apply to the following species:
c
c       H2O
c       H+ (or OH- if that is in the basis set instead of H+)
c       O2(g,aq) (or any other auxiliary basis species such as e-,
c                 HS-, or Fe+++, which might be used as the redox
c                 defining basis species)
c
c     In tracing a reaction path, it is undesirable, for example,
c     to suddenly run out of solvent water, or to skip over a redox
c     jump, without obtaining some level of resolution. This subroutine
c     forces smaller steps to be taken to avoid such problems. Also,
c     it insures that some level of descriptive detail is obtained
c     when a signficant, abrupt change is occurring in the system.
c
c     It is not the purpose of this subroutine to produce step size
c     values that correspond exactly to the specified constraints.
c     The methodology here is rather approximate. Changes in activity
c     coefficients along the reaction path are not accounted for.
c
c     In reducing delxi, the basic approach is based only on the
c     value at the base point and the first derivative. The full power
c     of the finite difference representations is therefore not
c     brought to bear.
c
c     This subroutine is called by:
c
c       EQ6/path.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c
c     Principal output:
c
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
      integer kmax,nbtmax,nodbmx,nrd1mx,nstmax
c
      integer noutpt,nttyo
c
      integer iodb(nodbmx),iindx1(kmax),nbasp(nbtmax)
c
      integer kbt,kdim,kelect,khydr,khydx,km1,ko2gaq,krdxsp,kwater,kxt,
     $ nord
c
      logical qrapch
c
      character*48 uspec(nstmax)
c
      real*8 dzvc0(nrd1mx,kmax),zvclg0(kmax),zvclg1(kmax),zvec0(kmax),
     $ zvec1(kmax)
c
      real*8 al10,delxi,dlxmin,zklogu
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,kcol,nb,ns,nstep,nstepl
c
      integer ilnobl
c
      logical qztayl
c
      character*24 unamsp
c
      real*8 lobasp,lelect,lhydr,lo2gaq,lwater,lrdxsp
c
      real*8 dlxic,dlxsv,zdif,zdzid,zdzspi,zdzwa,zdzwai,zx
c
c-----------------------------------------------------------------------
c
c     Maximum number of cycles cutting the step size in this subroutine.
c
      data nstepl /20/
c
c     Change limit values.
c
      data lelect /8.0/,lhydr /0.5/,lo2gaq /8.0/,lobasp /8.0/,
     $ lrdxsp /8.0/,lwater /0.50/
c
c-----------------------------------------------------------------------
c
      qrapch = .false.
      dlxsv = delxi
      unamsp = 'Error'
      nstep = -1
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Here is a return point if the value of delxi has been reduced by
c     this subroutine. This provides an additional check, as the
c     step-size reduction algorithms are only approximate.
c
  100 nstep = nstep + 1
c
      if (nstep .ge. nstepl) then
        if (iodb(1).gt.0 .or. iodb(5).gt.0) then
          write (noutpt,1000) nstepl
          write (nttyo,1000) nstepl
 1000     format(/' * Note - (EQ6/ldlxrc) The step size has been cut',
     $    ' the maximum ',i2,' times',/7x,'to limit rapid change in',
     $    ' one or more of the variables associated with the',
     $    /7x,'aqueous basis species. Cutting the step size to the',
     $    ' minimum value.')
        endif
        delxi = dlxmin
        go to 990
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Make a Taylor's series expansion of the z vector, without applying
c     change limits.
c
      qztayl = .false.
      call ztaylr(delxi,dzvc0,kdim,kmax,km1,kxt,nord,nrd1mx,
     $ qztayl,zklogu,zvclg0,zvclg1,zvec0,zvec1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (kwater .gt. 0) then
c
c       Number of moles of H2O. The limit is +/- lwater %.
c
        nb = iindx1(kwater)
        ns = nbasp(nb)
        dlxic = delxi
        zdif = zvec0(kwater) - zvec1(kwater)
        zx = zdif/zvec0(kwater)
        if (abs(zx) .gt. lwater) then
c
c         The predicted change is too big. Reduce the step size.
c
          if (dzvc0(1,kwater) .ne. 0.) then
c
c           Estimate the new step size from the finite-difference
c           data.
c
            zdzwa = zvec0(kwater)/dzvc0(1,kwater)
            dlxic = lwater*abs(zdzwa)
            if (dlxic .ge. delxi) dlxic = 0.5*delxi
          else
c
c           Just cut the step size.
c
            dlxic = 0.5*delxi
          endif
          dlxic = max(dlxic,dlxmin)
          if (dlxic .lt. delxi) then
            delxi = dlxic
            unamsp = uspec(ns)
            go to 100
          endif
        endif
      endif
c
      zdzwai = dzvc0(1,kwater)/zvec0(kwater)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (kwater .gt. 0) then
        if (khydr .gt. 0) then
c
c         The pH, if H+ is in the basis set. The limit is +/- lhydr
c         pH unit.
c
          nb = iindx1(khydr)
          ns = nbasp(nb)
          dlxic = delxi
          zdif = zvclg0(khydr) - zvclg1(khydr)
          if (abs(zdif) .gt. lhydr) then
c
c           The predicted change is too big. Reduce the step size.
c
            zdzspi = dzvc0(1,khydr)/zvec0(khydr)
            zdzid  = zdzspi - zdzwai
            if (zdzid .ne. 0.) then
c
c             Estimate the new step size from the finite-difference
c             data.
c
              dlxic = al10*lhydr/abs(zdzid)
              if (dlxic .ge. delxi) dlxic = 0.5*delxi
            else
c
c             Just cut the step size.
c
              dlxic = 0.5*delxi
            endif
            dlxic = max(dlxic,dlxmin)
            if (dlxic .lt. delxi) then
              delxi = dlxic
              unamsp = uspec(ns)
              go to 100
            endif
          endif
        elseif (khydx .gt. 0) then
c
c         The pH, if OH- is in the basis set instead of H+.
c         The limit is +/- lhydr pH unit.
c
          nb = iindx1(khydx)
          ns = nbasp(nb)
          dlxic = delxi
          zdif = zvclg0(khydx) - zvclg1(khydx)
          if (abs(zdif) .gt. lhydr) then
c
c           The predicted change is too big. Reduce the step size.
c
            zdzspi = dzvc0(1,khydx)/zvec0(khydx)
            zdzid  = zdzspi - zdzwai
            if (zdzid .ne. 0.) then
c
c             Estimate the new step size from the finite-difference
c             data.
c
              dlxic = al10*lhydr/abs(zdzid)
              if (dlxic .ge. delxi) dlxic = 0.5*delxi
            else
c
c             Just cut the step size.
c
              dlxic = 0.5*delxi
            endif
            dlxic = max(dlxic,dlxmin)
            if (dlxic .lt. delxi) then
              delxi = dlxic
              unamsp = uspec(ns)
              go to 100
            endif
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (ko2gaq .gt. 0) then
c
c       The log fO2. The limit is +/- lo2gaq log units.
c
        nb = iindx1(ko2gaq)
        ns = nbasp(nb)
        dlxic = delxi
        zdif = zvclg0(ko2gaq) - zvclg1(ko2gaq)
        if (abs(zdif) .gt. lo2gaq) then
          if (dzvc0(1,ko2gaq) .ne. 0.) then
c
c           Estimate the new step size from the finite-difference
c           data.
c
            dlxic = al10*lo2gaq*abs(zvec0(ko2gaq)/dzvc0(1,ko2gaq))
            if (dlxic .ge. delxi) dlxic = 0.5*delxi
          else
c
c           Just cut the step size.
c
            dlxic = 0.5*delxi
          endif
          dlxic = max(dlxic,dlxmin)
          if (dlxic .lt. delxi) then
            delxi = dlxic
            unamsp = uspec(ns)
            go to 100
          endif
        endif
      endif
c
      if (kelect .gt. 0) then
c
c       The pe. The limit is +/- lelect log units.
c
        nb = iindx1(kelect)
        ns = nbasp(nb)
        dlxic = delxi
        zdif = zvclg0(kelect) - zvclg1(kelect)
        if (abs(zdif) .gt. lelect) then
          if (dzvc0(1,kelect) .ne. 0.) then
c
c           Estimate the new step size from the finite-difference
c           data.
c
            dlxic = al10*lelect*abs(zvec0(kelect)/dzvc0(1,kelect))
            if (dlxic .ge. delxi) dlxic = 0.5*delxi
          else
c
c           Just cut the step size.
c
            dlxic = 0.5*delxi
          endif
          dlxic = max(dlxic,dlxmin)
          if (dlxic .lt. delxi) then
            delxi = dlxic
            unamsp = uspec(ns)
            go to 100
          endif
        endif
      endif
c
      if (krdxsp.gt.0 .and. kwater.gt.0) then
        if (krdxsp.ne.ko2gaq .and. krdxsp.ne.kelect) then
c
c         The log molality of an auxiliary basis species. The limit is
c         +/- lrdxsp log units.
c
          nb = iindx1(krdxsp)
          ns = nbasp(nb)
          dlxic = delxi
          zdif = zvclg0(krdxsp) - zvclg1(krdxsp)
          if (abs(zdif) .gt. lrdxsp) then
            zdzspi = dzvc0(1,krdxsp)/zvec0(krdxsp)
            zdzid  = zdzspi - zdzwai
            if (zdzid .ne. 0.) then
c
c             Estimate the new step size from the finite-difference
c             data.
c
              dlxic = al10*lrdxsp/abs(zdzid)
              if (dlxic .ge. delxi) dlxic = 0.5*delxi
            else
c
c             Just cut the step size.
c
              dlxic = 0.5*delxi
            endif
            dlxic = max(dlxic,dlxmin)
            if (dlxic .lt. delxi) then
              if (dlxic .gt. dlxmin) then
                delxi = dlxic
                unamsp = uspec(ns)
                go to 100
              endif
            endif
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      do kcol = 1,kbt
        if (kcol.ne.kwater .and. kcol.ne.khydr .and. kcol.ne.khydx
     $    .and. kcol.ne.krdxsp) then
c
c         The log molality of any other basis species. The limit is
c         +/- lobasp log units.
c
          nb = iindx1(kcol)
          ns = nbasp(nb)
          dlxic = delxi
          zdif = zvclg0(kcol) - zvclg1(kcol)
          if (abs(zdif) .gt. lobasp) then
            zdif = zvclg0(kcol) - zvclg1(kcol)
            if (abs(zdif) .gt. lobasp) then
              zdzspi = dzvc0(1,kcol)/zvec0(kcol)
              zdzid  = zdzspi - zdzwai
              if (zdzid .ne. 0.) then
c
c               Estimate the new step size from the finite-difference
c               data.
c
                dlxic = al10*lobasp/abs(zdzid)
                if (dlxic .ge. delxi) dlxic = 0.5*delxi
              else
c
c               Just cut the step size.
c
                dlxic = 0.5*delxi
              endif
              dlxic = max(dlxic,dlxmin)
              if (dlxic .lt. delxi) then
                delxi = dlxic
                unamsp = uspec(ns)
                go to 100
              endif
            endif
          endif
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  990 qrapch = delxi .lt. dlxsv
      if (qrapch) then
        if (iodb(5) .gt. 0) then
          j2 = ilnobl(unamsp)
          write (noutpt,1010) dlxsv,delxi,unamsp(1:j2)
          write (nttyo,1010) dlxsv,delxi,unamsp(1:j2)
 1010     format(/' * Note - (EQ6/ldlxrc) The step size has been cut',
     $    ' from ',1pe12.5,/7x,'to ',e12.5,' to limit rapid change',
     $    ' associated with ',a,'.')
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
