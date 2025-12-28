      subroutine setpk6(actwlg,awmax,awmaxi,awmin,awmini,eh,ehmax,
     $ ehmaxi,ehmin,ehmini,fo2lg,iindx1,jflag,jflgi,kbt,kdim,kmax,
     $ kprs,mprph,mprphi,mprsp,mprspi,mtb,mtbi,mtbaq,mtbaqi,nbasp,
     $ nbaspd,nbaspi,nbti,nbtmax,ncmpr,nobswt,noutpt,nprpmx,nprpti,
     $ nprsmx,nprsti,npt,nptmax,nttyo,nstmax,o2max,o2maxi,o2min,
     $ o2mini,ph,phmax,phmaxi,phmin,phmini,prcinf,press,pressi,
     $ tempc,tempci,time1,timemx,timmxi,tistti,ubmtbi,uobsw,uphase,
     $ uprphi,uprspi,uspec,uzveci,uzvec1,xi1,ximax,ximaxi,xistti,
     $ zvclgi,zvclg1)
c
c     This subroutine sets up certain variables and arrays for writing
c     on the pickup file. This subroutine must be called prior to
c     writing a pickup file.
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
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer kmax,nbtmax,nprpmx,nprsmx,nptmax,nstmax
c
      integer noutpt,nttyo
c
      integer iindx1(kmax),jflag(nstmax),jflgi(nbtmax),nbasp(nbtmax),
     $ nbaspd(nbtmax),nbaspi(nbtmax),ncmpr(2,nptmax)
c
      integer kbt,kdim,kprs,nbti,nobswt,nprpti,nprsti,npt
c
      character*48 ubmtbi(nbtmax),uobsw(2,nbtmax),uprspi(nprsmx),
     $ uspec(nstmax),uzveci(kmax),uzvec1(kmax)
      character*24 uphase(nptmax),uprphi(nprpmx)
c
      real*8 mprph(nptmax),mprphi(nprpmx),mprsp(nstmax),mprspi(nprsmx),
     $ mtb(nbtmax),mtbi(nbtmax),mtbaq(nbtmax),mtbaqi(nbtmax),
     $ zvclgi(kmax),zvclg1(kmax)
c
      real*8 actwlg,awmax,awmaxi,awmin,awmini,eh,ehmax,ehmaxi,ehmin,
     $ ehmini,fo2lg,o2max,o2maxi,o2min,o2mini,ph,phmax,phmaxi,phmin,
     $ phmini,prcinf,press,pressi,tempc,tempci,time1,timemx,timmxi,
     $ tistti,xi1,ximax,ximaxi,xistti
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer j2,kcol,krow,n,nb,nbi,nerr,nmax,np,nr1,nr2,ns,nsi,ns1,ns2
c
      integer ilnobl
c
      character*8 ux8
      character*24 uaqsln
c
c-----------------------------------------------------------------------
c
      data uaqsln /'Aqueous solution        '/
c
c-----------------------------------------------------------------------
c
      nerr = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Starting values of reaction progress and time.
c
      xistti = xi1
      tistti = time1
c
c     Maximum values of reaction progress and time.
c
      if (xi1 .ge. ximax) ximaxi = prcinf
      if (time1 .ge. timemx) timmxi = prcinf
c
c     Minimum and maximum values of other quantities.
c
      if (pH .le. phmin) phmini = -prcinf
      if (pH .ge. phmax) phmaxi = prcinf
      if (eh .le. ehmin) ehmini = -prcinf
      if (eh .ge. ehmax) ehmaxi = prcinf
      if (fo2lg .le. o2min) o2mini = -prcinf
      if (fo2lg .ge. o2max) o2maxi = prcinf
      if (actwlg .le. awmin) awmini = -prcinf
      if (actwlg .ge. awmax) awmaxi = prcinf
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Original temperature.
c
      tempci = tempc
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Original pressure.
c
      pressi = press
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Mass balance totals and jflag values.
c
      call initaz(mtbi,nbtmax)
      call initaz(mtbaqi,nbtmax)
      call initiz(jflgi,nbtmax)
c
      do krow = 1,kbt
        nb = iindx1(krow)
        ns = nbaspd(nb)
        nsi = nbaspi(nb)
        nbi = krow
        ubmtbi(nbi) = uspec(nsi)
        mtbi(nbi) = mtb(nb)
        mtbaqi(nbi) = mtbaq(nb)
        jflgi(nbi) = jflag(ns)
      enddo
      nbti = kbt
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Ordinary basis switching directives.
c
      nmax = 2*nbtmax
      call initcb(uobsw,nmax)
c
      n = 0
      do kcol = 1,kbt
        nb = iindx1(kcol)
        ns1 = nbaspd(nb)
        ns2 = nbasp(nb)
        if (ns1 .ne. ns2) then
          n = n + 1
          uobsw(1,n) = uspec(ns1)
          uobsw(2,n) = uspec(ns2)
        endif
      enddo
      nobswt = n
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Matrix column variables and corresponding values.
c
      call initcb(uzveci,kmax)
      call initaz(zvclgi,kmax)
c
      do kcol = 1,kdim
        uzveci(kcol) = uzvec1(kcol)
        zvclgi(kcol) = zvclg1(kcol)
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Phases and species in the PRS.
c
      kprs = 0
      nprpti = 0
      nprsti = 0
c
      call initcb(uprphi,nprpmx)
      call initcb(uprspi,nprsmx)
      call initaz(mprphi,nprpmx)
      call initaz(mprspi,nprsmx)
c
      do np = 1,npt
        if (uphase(np)(1:24) .ne. uaqsln(1:24)) then
          if (mprph(np) .gt. 0.) then
            nprpti = nprpti + 1
c
            if (nprpti .gt. nprpmx) then
              write (ux8,'(i5)') nprpmx
              call lejust(ux8)
              j2 = ilnobl(ux8)
              write (noutpt,1000) ux8(1:j2)
              write (nttyo,1000) ux8(1:j2)
 1000         format(/' * Error - (EQ6/setpk6) Have too many phases',
     $        ' in the physically removed system',/7x,'(PRS) to write',
     $        ' them all on the pickup file. The code is only',
     $        ' dimensioned',/7x,'to allow ',a,' such phases on the',
     $        ' input and pickup files. Increase',/7x,'the',
     $        ' dimensioning parameter nprppa.')
              nerr = nerr + 1
            endif
c
            uprphi(nprpti) = uphase(np)
            mprphi(nprpti) = mprph(np)
            nr1 = ncmpr(1,np)
            nr2 = ncmpr(2,np)
            do ns = nr1,nr2
              if (mprsp(ns) .gt. 0.) then
                nprsti = nprsti + 1
c
                if (nprsti .gt. nprsmx) then
                  write (ux8,'(i5)') nprsmx
                  call lejust(ux8)
                  j2 = ilnobl(ux8)
                  write (noutpt,1010) ux8
                  write (nttyo,1010) ux8
 1010             format(/' * Error - (EQ6/setpk6) Have too many',
     $            ' species in the physically removed system',
     $            /7x,'(PRS) to write them all on the pickup file.',
     $            ' The code is only dimensioned',/7x,'to allow ',a,
     $            ' such species on the input and pickup files.',
     $            ' Increase',/7x,'the dimensioning parameter nprspa.')
                  nerr = nerr + 1
                endif
c
                uprspi(nprsti) = uspec(ns)
                mprspi(nprsti) = mprsp(ns)
              endif
            enddo
          endif
        endif
      enddo
      if (nprpti .gt. 0) kprs = 1
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nerr .gt. 0) stop
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
