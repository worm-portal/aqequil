      subroutine eqphas(aamatr,abar,acflg,acflgo,act,actlg,
     $ adh,adhh,adhv,afcnst,affp,affs,alpha,al10,amtb,aphi,
     $ apx,avcnst,azero,a3bar,a3bars,bacfmx,bbig,bdh,bdhh,bdhv,
     $ bdot,bdoth,bdotv,beta,betamx,betao,bgamx,bneg,bpx,cco2,
     $ cegexs,cess,cdrs,cdrsd,cdrsx,cdrtw,cdrw,cjbasp,cnufac,
     $ conc,conclg,cpgexs,cscale,csts,delvco,delvec,d1zvc1,dlogxw,
     $ egexjc,egexjf,egexs,eh,ehfac,eps100,farad,fje,fjeo,fo2,
     $ fo2lg,fsort,fugac,fugalg,fxi,fxio,gmmatr,iapxt,ibpxt,ibswx,
     $ ielam,ier,iern1,iern2,ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,
     $ ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,iindx0,iindx1,ilcphi1,
     $ ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,ilzeta,imrn1,
     $ imrn2,insgf,iodb,iopg,iopt,ipch,ipivot,ipndx1,ipcv,istack,
     $ iter,itermx,ixbasp,ixrn1,ixrn2,izmax,jcsort,jflag,jgsort,
     $ jgstak,jjsort,jpflag,jsflag,jsitex,jsol,jssort,jstack,kbt,
     $ kction,kdim,kelect,khydr,khydx,km1,km10,kmt,kmt0,ko2gaq,
     $ kpsat,kpsst,krdxsp,kwater,kx1,kx10,kxt0,kxt,loph,losp,lsort,
     $ moph,mosp,mrgexs,mtb,mtbaq,narn1,narn2,narxt,nat,nbasp,nbaspd,
     $ nbaspx,nbt,nbtd,nbw,nchlor,ncmpr,nct,ndrs,ndrsd,ndrsx,ndrsr,
     $ ndrsrd,ndrsrx,nelect,nern1,nern2,ness,nessr,net,nfrn1,nfrn2,
     $ ngrn1,ngrn2,ngt,nhydr,nhydx,nlrn1,nlrn2,nlt,nmrn1,nmrn2,
     $ nmt,nord,no2gaq,noutpt,npchk,nphasx,npt,nrdxsp,nst,nsts,
     $ nstsr,ntpr,ntrymx,nttyo,nxrn1,nxrn2,nxt,omega,omeglg,prcinf,
     $ press,qbassw,qbseqc,qbye,qcnpre,qcntmp,qhawep,qmod,qoptmz,
     $ qpit75,qredox,qsspgb,qstart,qxknph,q6mode,rcnstv,rconst,
     $ rhsvec,rtcnst,screwd,sidrph,sidrsp,sigmam,sigmmo,smp100,tempc,
     $ tempk,tolbt,toldl,tolsat,tolsst,ubacmx,ubgamx,ulbeta,uldel,
     $ uphase,uspec,uzvec1,weight,wfac,xbar,xbarlg,xbarw,xbarwc,
     $ xbrwlc,xbrwlg,zchar,zchcu6,zchsq2,zvclg1,zvec1)
c
c     This subroutine attempts to calculate the equilibrium state of
c     the equilibrium system. As necessary, it makes repeated calls
c     to EQ6/eqcalc.f, which attempts to make this calculation for a
c     given specified phase assemblage. The present subroutine adjusts
c     this assemblage as needed.
c
c     This subroutine is called by:
c
c       EQ6/eqshel.f
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
c       qbye   = flag which is .true. if EQ6/path.f just changed the
c                  phase assemblage by deleting one or more phases.
c                  This condition instructs the present subroutine to
c                  print the index structure of the system of equations
c                  as it does on the starting call and whenever this
c                  subroutine itself adds or deletes a phase.
c
c       iter   = the number of hydrbid Newton-Raphson iterations
c                  done by EQLIB/newton.f
c
c       ier    = error flag:
c
c                  Values returned from EQ6/eqcalc.f:
c                    =    0  Okay
c                    =    1  Encountered a zero matrix
c                    =    2  Encountered a non-zero, computationally
c                              singular matrix
c                    =    3  Iteration was diverging
c                    =    4  Hit the maximum number of iterations
c                              (itermx)
c
c                  Values returned from EQ6/miidxz.f:
c                    =    0  Okay
c                    =    8  Exceeded the dimension of the iindx1
c                              array
c
c                  Values returned by the present subroutine:
c                    =    0  Okay
c                    =    8  Exceeded the dimension of the iindx1
c                              array; go back and cut the step size
c                              if possible
c                    =   10  Go back and take a smaller step size to
c                              avoid exceeding the supersaturation
c                              tolerance (tolsst). This is done only
c                              when an appropriate set of conditions
c                              is satisified. It isn't necessary for
c                              EQ6/path.f to analyze the situation
c                              when this value is returned to it by
c                              EQ6/ eqshel.f.
c                    =   20  Hit the maximum number of tries to find
c                              the correct phase assemblage
c                    =   30  Caught in a region of computational
c                              instability about the phase boundary for
c                              an appearing phase. Iteration fails when
c                              the phase is added to the phase
c                              assemblage.
c                    =   40  Caught in a region of computational
c                              instability about the phase boundary for
c                              a disappearing phase. The system is
c                              too supersaturated with respect to a
c                              phase if that phase is deleted from the
c                              phase assemblage.
c                    =   50  Caught in a region of computational
c                              instability associated with the
c                              master redox variable. This is probably
c                              associated with a redox jump.
c                    =   60  Caught in a region of computational
c                              instability associated with solvent
c                              water. The amount of water in the
c                              system is probably very low.
c                    =   70  Caught in a region of computational
c                              instability associated with the
c                              aqueous activity coefficient model.
c                    =  100  Detected out-of-range values for variables
c                              associated with the basis species before
c                              starting iteration.
c                    =  110  Detected out-of-range values for variables
c                              associated with the basis species after
c                              an iteration crash.
c                    =  150  Calculation failed, no diagnostics were
c                              generated.
c
c-----------------------------------------------------------------------
c
c     Modules.
c
c     The module mod6pt contains data required to evaluate Pitzer's
c     equations.
c
      use mod6pt
c
c     The module mod6xf contains most of the standard-state
c     thermodynamic data.
c
      use mod6xf
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
      include 'eqlib/eqlpar.h'
      include 'eqlib/eqldv.h'
      include 'eqlib/eqlge.h'
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer noutpt,nttyo
c
      integer iapxt(nxtmax),ibpxt(nxtmax),ibswx(nbtmax),
     $ igstak(ngtmax),iindx0(kmax),iindx1(kmax),ipivot(kmax),
     $ ipndx1(kmax),insgf(natmax),iodb(nodbmx),iopg(nopgmx),
     $ iopt(noptmx),istack(nstmax),ixbasp(nbtmax),jcsort(nstmax),
     $ jflag(nstmax),jgsort(ngtmax),jgstak(ngtmax),jjsort(nstmax),
     $ jpflag(nptmax),jsflag(nstmax),jsitex(nstmax),jsol(nxtmax),
     $ jssort(nstmax),jstack(nstmax),kction(nbtmax)
c
      integer narxt(ntprmx),nbasp(nbtmax),nbaspd(nbtmax),nbaspx(nbtmax),
     $ ncmpr(2,nptmax),ness(nessmx),nessr(2,nstmax),ndrs(ndrsmx),
     $ ndrsd(ndrsmx),ndrsx(ndrsmx),ndrsr(2,nstmax),ndrsrd(2,nstmax),
     $ ndrsrx(2,nstmax),nfac(nbtmax),npchk(nptmax),nphasx(nstmax),
     $ nsts(nstsmx),nstsr(2,nstmax)
c
      integer iern1,iern2,ifrn1,ifrn2,ilrn1,ilrn2,imrn1,imrn2,
     $ ixrn1,ixrn2
c
      integer ifcphi1,ifcphi2,ifnnn,ifn2n,ifpsi1,ifpsi2,ifzeta,
     $ ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilzeta
c
      integer nat,nbt,nct,net,ngt,nlt,nmt,npt,nst,nxt
c
      integer narn1,narn2,nern1,nern2,nfrn1,nfrn2,ngrn1,ngrn2,
     $ nlrn1,nlrn2,nmrn1,nmrn2,nxrn1,nxrn2
c
      integer ielam,ier,igas,ipch,ipcv,iter,itermx,izmax,kbt,kdim,
     $ kelect,khydr,khydx,km1,km10,kmt,kmt0,ko2gaq,kpsat,kpsst,krdxsp,
     $ kwater,kx1,kx10,kxt,kxt0,nbtd,nbw,nchlor,nelect,nhydr,nhydx,
     $ nord,no2gaq,nrdxsp,ntpr,ntrymx
c
      logical qxknph(nptmax)
c
      logical qbassw,qbseqc,qbye,qcnpre,qcntmp,qhawep,qmod,qoptmz,
     $ qpit75,qredox,qsspgb,qstart,q6mode
c
      character*48 uspec(nstmax),uzvec1(kmax)
      character*48 ubacmx,ubgamx
      character*24 uphase(nptmax)
      character*8 ulbeta(kmax),uldel(kmax)
c
      real*8 aamatr(kmax,kmax),acflg(nstmax),acflgo(nstmax),
     $ act(nstmax),actlg(nstmax),affp(nptmax),affs(nstmax),
     $ alpha(kmax),amtb(nbtmax),apx(iapxmx,nxtmax),azero(natmax),
     $ a3bars(natmax),beta(kmax),betao(kmax),bpx(ibpxmx,nxtmax),
     $ cco2(5),cegexs(ietmax,jetmax,netmax),cess(nessmx),
     $ cdrs(ndrsmx),cdrsd(ndrsmx),cdrsx(ndrsmx),cdrtw(nstmax),
     $ cdrw(nstmax),cjbasp(nbtmax),cnufac(nstmax),conc(nstmax),
     $ conclg(nstmax),cpgexs(ietmax,jetmax,netmax),cscale(nstmax),
     $ csts(nstsmx),delvco(kmax),delvec(kmax),dlogxw(nbtmax),
     $ d1zvc1(kmax),egexjc(jetmax,netmax),egexjf(jetmax,netmax),
     $ egexs(ietmax,jetmax,netmax)
c
      real*8 fsort(ngtmax),fugac(ngtmax),fugalg(ngtmax),
     $ gmmatr(kmax,kmax),loph(nptmax),
     $ losp(nstmax),lsort(nstmax),moph(nptmax),mosp(nstmax),
     $ mrgexs(ietmax,jetmax,netmax),mtb(nbtmax),mtbaq(nbtmax),
     $ rhsvec(kmax),sidrph(nptmax),sidrsp(nstmax),weight(nstmax),
     $ wfac(iktmax,nxtmax),xbar(nstmax),xbarlg(nstmax),
     $ zchar(nstmax),zchcu6(nstmax),zchsq2(nstmax),
     $ zvclg1(kmax),zvec1(kmax)
c
      real*8 adh,adhh,adhv,aphi,bdh,bdhh,bdhv,bdot,bdoth,bdotv
c
      real*8 abar,afcnst,al10,avcnst,a3bar,bacfmx,bbig,betamx,bgamx,
     $ bneg,eh,ehfac,eps100,farad,fje,fjeo,fo2,fo2lg,fxi,fxio,omega,
     $ omeglg,prcinf,press,rcnstv,rconst,rtcnst,screwd,sigmam,sigmmo,
     $ smp100,tempc,tempk,tolbt,toldl,tolsat,tolsst,xbarw,xbarwc,xbrwlc,
     $ xbrwlg
c
c-----------------------------------------------------------------------
c
c     Local variable declarations with global dimensioning.
c
c     Variables for restoring the entering configuration, including
c     the entering phase assemblage.
c
      integer isv_kmax,isv_nbtmax,isv_nstmax
c
      SAVE isv_kmax,isv_nbtmax,isv_nstmax
c
      integer, dimension(:), allocatable :: iindxs,ipndxs,nbasps
c
      SAVE iindxs,ipndxs,nbasps
c
      real(8), dimension(:), allocatable :: acflgs,zvclgs
c
      SAVE acflgs,zvclgs
c
c     The following do not need to be SAVEd.
c
      integer kdims,km1s,kmts,kx1s,kxts
c
      real(8) xbarws,xbrwls
c
c-----------------------------------------------------------------------
c
c     Local variable declarations with special dimensioning.
c
c     Data for the eight most supersaturated minerals.
c
      integer nssppa
      parameter (nssppa = 8)
c
      integer nsspmx
c
      integer nssp(nssppa)
c
      character*24 ussp(nssppa)
c
      real*8 afssp(nssppa),afssps(nssppa),msspmx(nssppa)
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer ibetmx,idelmx,ifail,inext,irow1,irow2,isave,isspt,j,
     $ jcol1,jcol2,jj,jkl,jlen,j2,j3,k,kcol,n,nb,nerr,nords,np,npadd,
     $ npaddi,npdel,npdeli,nrn1,nrn2,nr1,nr2,ns,nss,nsspt,nswtch,ns2,
     $ nt,ntry,ntryai,ntrydi,nxx
c
      integer ilnobl
c
      logical qadd,qbswok,qphruv
c
      character*56 uspn56
      character*8 ux8
c
      real*8 afscal,bfje,bfxi,bsigmm,cx,cxx,delmax,mophmx,mophn,mx,
     $ mxmadd,mpmxaq,mpmxes,msmxaq,msmxes
c
      real*8 tlg
c
c-----------------------------------------------------------------------
c
c     Local dimensioning variables.
c
      nsspmx = nssppa
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Allocate or reallocate local work arrays as needed.
c
      if (.not.ALLOCATED(iindxs)) then
c
c       Local work arrays are not allocated. Zero the saved
c       array size variables. Note that only one array is tested
c       to see if it is allocated. It is assumed that all local
c       work arrays are either allocated or not.
c
        isv_kmax = 0
        isv_nbtmax = 0
        isv_nstmax = 0
      else
c
c       Local work arrays are allocated. Check to see if any of the
c       array size variables have changed. If so, deallocate
c       the corresponding local work arrays and zero the corresponding
c       saved size variables.
c
        if (kmax .ne. isv_kmax) then
          DEALLOCATE(iindxs,ipndxs)
          DEALLOCATE(zvclgs)
          isv_kmax = 0
        endif
c
        if (nbtmax .ne. isv_nbtmax) then
          DEALLOCATE(nbasps)
          isv_nbtmax = 0
        endif
c
        if (nstmax .ne. isv_nstmax) then
          DEALLOCATE(acflgs)
          isv_nstmax = 0
        endif
      endif
c
c     At this point, the saved array size values are zero if the
c     corresponding arrays need to be allocated.
c
      if (isv_kmax .eq. 0) then
        ALLOCATE(iindxs(kmax),ipndxs(kmax))
        ALLOCATE(zvclgs(kmax))
        isv_kmax = kmax
      endif
c
      if (isv_nbtmax .eq. 0) then
        ALLOCATE(nbasps(nbtmax))
        isv_nbtmax = nbtmax
      endif
c
      if (isv_nstmax .eq. 0) then
        ALLOCATE(acflgs(nstmax))
        isv_nstmax = nstmax
      endif
c
c     Zero the contents of the local work arrays.
c
      do k = 1,kmax
        iindxs(k) = 0
        ipndxs(k) = 0
      enddo
c
      do k = 1,kmax
        zvclgs(k) = -99999.
      enddo
c
      do n = 1,nbtmax
        nbasps(n) = 0
      enddo
c
      do n = 1,nstmax
        acflgs(n) = 0.
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Initialize some parameters.
c
      ier = 0
c
      qadd = .false.
      qmod = .false.
      qbseqc = .false.
c
      ifail = 0
c
      npadd = 0
      npaddi = 0
c
      npdel = 0
      npdeli = 0
c
      ntry = 1
      ntryai = 0
      ntrydi = 0
c
      kpsst = 0
      nsspt = 0
      call initiz(nssp,nsspmx)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Here is a return point if the phase assemblage has been modified.
c
  100 bgamx = 0.
      bsigmm = 0.
      bfxi = 0.
      bfje = 0.
c
      if (ntry.gt.1 .or. qstart .or. qbye .or. iodb(1).ge.3) then
c
c       Print the current phase assemblage.
c
        write (noutpt,1000) ntry
 1000   format(/' Attempted phase assemblage number ',i3,/)
        do kcol = 1,kbt
          nb = iindx1(kcol)
          ns = nbaspd(nb)
c
c         Calling sequence substitutions:
c           uspec(ns) for unam48
c
          call fmspnx(jlen,uspec(ns),uspn56)
          write (noutpt,1010) kcol,uspn56(1:jlen)
 1010     format(2x,i3,2x,a)
        enddo
        do kcol = km1,kxt
          ns = iindx1(kcol)
c
c         Calling sequence substitutions:
c           uspec(ns) for unam48
c
          call fmspnm(jlen,uspec(ns),uspn56)
          write (noutpt,1010) kcol,uspn56(1:jlen)
        enddo
        write (noutpt,1020)
 1020   format(1x)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Save the "kernel" description for the current equilibrium system.
c     This will be used to recover if an equilibrium calculation made
c     by EQ6/eqcalc.f fails for any reason. This kernel contains the
c     minimum of information required to calculate the fully expanded
c     description using EQLIB/ncmpex.f. Note that EQ6/eqshel.f, which
c     calls the present subroutine, has its own backup kernel.
c     EQ6/path.f, which calls EQ6/eqshel.f, also has its own backup
c     kernel, which is used to restore the equilibrium system at the
c     previous point of reaction progress.
c
c     Note that the present backup kernel is not a compact description
c     of the equilibrium system from the last call to EQ6/eqcalc.f.
c     Upon entry to the present subroutine, the backup kernel
c     corresponds to a description of the system either as read from
c     the input file (at the starting value of reaction progress) or
c     as obtained by finite-difference prediction in stepping to the
c     current point of reaction progress. Otherwise, it corresponds to
c     what resulted from the last successful call to EQ6/eqcalc.f from
c     the present subroutine, modified by the addition, deletion, or
c     replacement of one phase.
c
      km1s = km1
      kmts = kmt
      kx1s = kx1
      kxts = kxt
      kdims = kdim
c
      call copyia(nbasp,nbasps,nbt)
      call copyia(iindx1,iindxs,kdim)
      call copyia(ipndx1,ipndxs,kdim)
c
      call copyaa(zvclg1,zvclgs,kdim)
      call copyaa(acflg,acflgs,nst)
c
      xbarws = xbarwc
      xbrwls = xbrwlc
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Save the order of the finite differences.
c
      nords = nord
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check for aqueous basis species (or any other species switched
c     with such) for cases in which the log number of moles variable
c     is -99999. This implies a number of moles values which is zero,
c     which in turn implies that the species is not actually present in
c     the equilibrium system. Note that in this code -99999. is the
c     conventional value for log10(0). Conversely, 10**(-99999.) returns
c     a value of zero.
c
      nerr = 0
      do kcol = 1,kbt
        if (zvclg1(kcol) .le. -99999.) then
c
c         Calling sequence substitutions:
c           uzvec1 for unam48
c
          call fmspnx(jlen,uzvec1(kcol),uspn56)
          if (kcol .ne. krdxsp) then
            write (noutpt,1050) uspn56(1:jlen),zvclg1(kcol)
            write (nttyo,1050) uspn56(1:jlen),zvclg1(kcol)
 1050       format(/' * Warning - (EQ6/eqphas) The log number of moles',
     $      ' variable for basis species',/7x,a,' has the out-of-range',
     $      ' value ',1pe12.5,'.')
          elseif (qredox) then
            if (krdxsp .eq. ko2gaq) then
              write (noutpt,1060) zvclg1(kcol)
              write (nttyo,1060) zvclg1(kcol)
 1060         format(/' * Warning - (EQ6/eqphas) The log oxygen',
     $        ' fugacity variable has the out-of-range',/7x,'value ',
     $        1pe12.5,'.')
            elseif (krdxsp .eq. kelect) then
              write (noutpt,1070) zvclg1(kcol)
              write (nttyo,1070) zvclg1(kcol)
 1070         format(/' * Warning - (EQ6/eqphas) The log electron',
     $        ' activity variable has the',/7x,'out-of-range value',
     $        ' ',1pe12.5,'.')
            else
              write (noutpt,1050) uspn56(1:jlen),zvclg1(kcol)
              write (nttyo,1050) uspn56(1:jlen),zvclg1(kcol)
            endif
          endif
          nerr = nerr + 1
        endif
      enddo
      if (nerr .gt. 0) then
        ier = 100
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check the current phase assemblage for mineral species for which
c     the log number of moles variable is -99999. This implies that the
c     corresponding number of moles variable is zero, and in turn that
c     the species is not actually present in the equilibrium system.
c     In such a case, the corresponding phase must be deleted from the
c     current phase assemblage.
c
      do kcol = km1,kxt
        if (zvclg1(kcol) .le. -99999.) then
          npdel = ipndx1(kcol)
          go to 700
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      call eqcalc(aamatr,abar,acflg,acflgo,act,actlg,adh,
     $ adhh,adhv,afcnst,alpha,al10,amtb,aphi,apx,avcnst,azero,
     $ a3bar,a3bars,bacfmx,bbig,bdh,bdhh,bdhv,bdot,bdoth,bdotv,
     $ beta,betamx,betao,bfje,bfxi,bgamx,bneg,bpx,bsigmm,cco2,
     $ cegexs,cess,cdrs,cdrsd,cdrsx,cdrtw,cdrw,cjbasp,cnufac,conc,
     $ conclg,cpgexs,cscale,csts,delmax,delvco,delvec,dlogxw,
     $ egexjc,egexjf,egexs,eh,ehfac,eps100,farad,fje,fjeo,fo2,
     $ fo2lg,fsort,fugac,fugalg,fxi,fxio,gmmatr,iapxt,ibetmx,ibpxt,
     $ ibswx,idelmx,ielam,ier,iern1,iern2,ifcphi1,ifcphi2,ifnnn,
     $ ifn2n,ifpsi1,ifpsi2,ifrn1,ifrn2,ifzeta,igas,igstak,iindx1,
     $ ilcphi1,ilcphi2,ilnnn,iln2n,ilpsi1,ilpsi2,ilrn1,ilrn2,
     $ ilzeta,imrn1,imrn2,insgf,iodb,iopg,iopt,ipch,ipcv,ipivot,
     $ ipndx1,istack,iter,itermx,ixbasp,ixrn1,ixrn2,izmax,jcsort,
     $ jflag,jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jsol,
     $ jssort,jstack,kbt,kction,kdim,kelect,khydr,khydx,km1,kmt,
     $ ko2gaq,krdxsp,kwater,kx1,kxt,loph,losp,lsort,moph,mosp,
     $ mrgexs,mtb,narn1,narn2,narxt,nat,nbasp,nbaspd,nbaspx,nbt,
     $ nbtd,nbw,nchlor,ncmpr,nct,ndrs,ndrsd,ndrsx,ndrsr,ndrsrd,
     $ ndrsrx,nelect,nern1,nern2,ness,nessr,net,nfac,nfrn1,nfrn2,
     $ ngrn1,ngrn2,ngt,nhydr,nhydx,nlrn1,nlrn2,nlt,nmrn1,nmrn2,
     $ nmt,noutpt,no2gaq,nphasx,npt,nrdxsp,nst,nsts,nstsr,ntpr,
     $ nttyo,nxrn1,nxrn2,nxt,omega,omeglg,press,qbassw,qhawep,
     $ qoptmz,qpit75,qredox,q6mode,rhsvec,screwd,sigmam,sigmmo,
     $ smp100,tempc,tempk,tolbt,toldl,ubacmx,ubgamx,ulbeta,uldel,
     $ uphase,uspec,uzvec1,weight,wfac,xbar,xbarlg,xbarw,xbarwc,
     $ xbrwlc,xbrwlg,zchar,zchcu6,zchsq2,zvclg1,zvec1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Were any basis switches made by EQ6/optmzr.f?
c
      nswtch = 0
      do kcol = 1,kbt
        nb = iindx1(kcol)
        ns = nbasp(nb)
        ns2 = nbasps(nb)
        if (ns2 .ne. ns) then
          nswtch = nswtch + 1
        endif
      enddo
      qbseqc = nswtch .gt. 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qbseqc) then
c
c       Drop the order to zero. The finite difference data are tied
c       to the old basis set. The original order and the use of these
c       data can be restored later in the present subroutine if the
c       switches just made are undone.
c
        nord = nords
      endif
c
      if (ier .gt. 0) go to 500
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Hybrid Newton-Raphson iteration has converged. Check for super-
c     saturations (other than those which are permitted to exist
c     metastably).
c
      ifail = 0
      npadd = 0
      nsspt = 0
      call initiz(nssp,nsspmx)
c
c     Calculate the total number of moles of each basis species
c     present in the aqueous phase.
c
      call initaz(mtbaq,nbt)
c
      do nss = narn1,narn2
        ns = jcsort(nss)
        if (mosp(ns) .ne. 0.) then
          nr1 = nstsr(1,ns)
          nr2 = nstsr(2,ns)
          do n = nr1,nr2
            nb = nsts(n)
            mtbaq(nb) = mtbaq(nb) + csts(n)*mosp(ns)
          enddo
        endif
      enddo
c
c     Check for supersaturations.
c
      call satchk(acflg,act,actlg,afcnst,affp,affs,apx,bpx,
     $ cdrs,eps100,iindx1,iodb,iopt,iapxmx,ibpxmx,iktmax,ixrn1,jflag,
     $ jpflag,jsflag,jsol,kmax,km1,kpsat,kpsst,kxt,nbasp,nbt,nbtmax,
     $ ncmpr,ndrs,ndrsmx,ndrsr,nodbmx,noptmx,noutpt,npchk,npt,nptmax,
     $ nstmax,nttyo,nxrn1,nxrn2,nxtmax,qxknph,sidrph,sidrsp,tolsat,
     $ uphase,uspec,wfac,xbar,xbarlg,xlks)
c
      if (kpsst .gt. 0) then
        if (iodb(4) .le. 0) then
          ux8 = ' '
          write (ux8,'(i7)') iter
          call lejust(ux8)
          j2 = ilnobl(ux8)
          write (noutpt,1090) ux8
 1090     format(' iter= ',a)
        endif
c
        ux8 = ' '
        write (ux8,'(i7)') kpsst
        call lejust(ux8)
        j2 = ilnobl(ux8)
        write (noutpt,1100) ux8(1:j2)
 1100   format(/'   Have ',a,' supersaturated phases.')
      endif
c
      if (kpsst .le. 0) then
c
c       The last equilibrium calculation converged and there are no
c       supersaturations to eliminate. The calculation overseen by
c       the present subroutine has completed successfully.
c
        go to 990
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Have one or more supersaturations. Find a phase to add to the
c     current phase assemblage. To do this, scale the affinities of the
c     supersaturated phases according to the number of ions produced
c     and destroyed in the respective dissolution reactions. The
c     scaled affinity is arbitrarily reduced if the maximum possible
c     number of moles of the new phase as estimated from the aqueous
c     solution chemistry is very small. The supersaturated phases are
c     arranged in decreasing order of scaled affinity. The phase with
c     the greatest scaled affinity is added to the phase assemblage.
c     In general, this algorithm picks the correct phase about four
c     out of five times.
c
      call initcb(ussp,nsspmx)
      call initaz(afssp,nsspmx)
      call initaz(afssps,nsspmx)
      call initaz(msspmx,nsspmx)
c
      isspt = 0
      do np = 1,npt
        if (jpflag(np) .eq. -2) then
          jpflag(np) = 0
          isspt = isspt + 1
c
          cxx = 0.
          nxx = 0
          nrn1 = ncmpr(1,np)
          nrn2 = ncmpr(2,np)
          nt = nrn2 - nrn1 + 1
          mpmxaq = 0.
          mpmxes = 0.
c
          do ns = nrn1,nrn2
            if (jsflag(ns) .le. 0) then
              cxx = cxx + cscale(ns)
              nxx = nxx + 1
              msmxaq = 1.e+38
              msmxes = 1.e+38
              nr1 = nstsr(1,ns)
              nr2 = nstsr(2,ns)
              do n = nr1,nr2
                cx = csts(n)
                if (cx .ne. 0.) then
                  nb = nsts(n)
                  nss = nbaspd(nb)
                  if (nss.ne.nhydr .and. nss.ne.nhydx .and.
     $              nss.ne.no2gaq .and. nss.ne.nelect) then
c
c                   Note: the total number of moles of a basis species
c                   can not be used to bound the number of moles of
c                   the species currently being examined (and hence
c                   also the number of moles of the corresponding phase)
c                   if the current basis species is H+, OH-, O2(g,aq),
c                   or e-. The latter two species are fictive. More to
c                   the point, the total number of moles of such
c                   basis species may be zero or even a negative number.
c
                    if (cx .gt. 0.) then
                      if (mtbaq(nb) .gt. 0.) then
                        mx = mtbaq(nb)/cx
                        msmxaq = min(msmxaq,mx)
                      endif
                      if (mtb(nb) .gt. 0.) then
                        mx = mtb(nb)/cx
                        msmxes = min(msmxes,mx)
                      endif
                    endif
                  endif
                endif
              enddo
              if (msmxaq .ge. 1.e+38) msmxaq = 1.0
              if (msmxes .ge. 1.e+38) msmxes = 1.0
              mpmxaq = max(mpmxaq,msmxaq)
              mpmxes = max(mpmxes,msmxes)
            endif
          enddo
c
c         Pick a maximum number of moles for the current phase based
c         on the moles of components in the aqueous solution and in
c         the equilibrium system.
c
          mophmx = max(mpmxaq,(0.05*mpmxes))
c
c         If the phase consists of more than one species, compute and
c         use the average affinity scaling factor.
c
          cxx = cxx/nxx
          afscal = affp(np)/cxx
c
c         Reduce the scaled affinity if the maximum precipitable
c         number of moles is small.
c
          if (mpmxaq .le. 5.e-6) afscal = 0.01*afscal
c
          do jj = 1,nsspmx
            if (afscal .ge. afssps(jj)) then
              jkl = nsspmx - jj
c
              do j = 1,jkl
                k = nsspmx - j
                nssp(k + 1) = nssp(k)
                ussp(k + 1) = ussp(k)
                afssps(k + 1) = afssps(k)
                afssp(k + 1) = afssp(k)
                msspmx(k + 1) = msspmx(k)
              enddo
c
              nssp(jj) = np
              ussp(jj) = uphase(np)
              afssps(jj) = afscal
              afssp(jj) = affp(np)
              msspmx(jj) = mophmx
              go to 200
            endif
          enddo
  200     continue
        endif
      enddo
      nsspt = min(isspt,nsspmx)
c
      if (iodb(1) .ge. 2) then
        write (noutpt,1110)
 1110   format(//5x,'The most supersaturated phases:',
     $  //31x,'Scaled',6x,'Affinity,',5x,'Maximum',
     $  /3x,'Name',23x,'Affinity',6x,'kcal/mol',4x,'Precipitable',
     $  /59x,'Moles',/)
        do n = 1,nsspt
          write (noutpt,1120) ussp(n),afssps(n),afssp(n),
     $    msspmx(n)
 1120     format(1x,a24,2x,f12.7,2x,f12.7,2x,1pe12.5)
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (afssp(1) .gt. tolsst) then
        if (qsspgb .and. ntry.le.1) then
c
c         Return a special error flag so that subroutine path will cut
c         the step size in order to be able to home in on the value of
c         reaction progress which corresponds to the target
c         supersaturation.
c
          write (noutpt,1150)
 1150     format(/' --- The extent of supersaturation exceeds the',
     $    ' normal tolerance ---')
          ier = 10
          go to 999
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      npadd = nssp(1)
      mxmadd = msspmx(1)
c
      if (npadd.eq.npaddi .and. ntry.eq.(ntryai + 2)) then
c
c       The phase chosen to add was added two steps previously and
c       removed on the previous step.
c
        if (iodb(1) .gt. 0) then
          write (noutpt,1160)
          write (nttyo,1160)
 1160     format(/' * Note - (EQ6/eqphas) Caught in a region of',
     $    /7x,'critical instability in the ES phase assemblage. A',
     $    /7x,"supersaturated phase can't be precipitated due to",
     $    /7x,'lack of convergence when it is included in the ES',
     $    /7x,'phase assemblage.')
        endif
c
        ier = 30
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     The following is a return point to add one new member to the
c     current phase assemblage and try again.
c
  300 j2 = ilnobl(uphase(npadd))
      write (ux8,'(i7)') npadd
      call lejust(ux8)
      j3 = ilnobl(ux8)
      write (noutpt,1180) uphase(npadd)(1:j2),ux8(1:j3)
 1180 format(/'   The phase to be added is ',a,1x,'(',a,')')
c
      qadd = .true.
      npaddi = npadd
      ntryai = ntry
c
      do kcol = 1,kmax
        delvec(kcol) = 0.
      enddo
c
c     The starting number of moles is set to some fraction of the
c     probable maximum value. Here it is set to 100% of that.
c
      mophn = mxmadd
      np = npadd
      jpflag(np) = -1
      loph(np) = tlg(mophn)
      nr1 = ncmpr(1,np)
      nr2 = ncmpr(2,np)
      nt = nr2 - nr1 + 1
      if (nt .eq. 1) then
        jsflag(nr1) = -1
        losp(nr1) = loph(np)
      else
        do ns = nr1,nr2
          if (jsflag(ns) .le. 0) then
            jsflag(ns) = -1
            losp(ns) = loph(np) + xbarlg(ns)
          endif
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Modify the current phase assemblage.
c
  400 ntry = ntry + 1
c
      if (ntry .gt. ntrymx) then
        write (noutpt,1200) ntrymx
        write (nttyo,1200) ntrymx
 1200   format(/' * Note - (EQ6/eqphas) Have done the maximum ',i3,
     $  ' tries to find the',/7x,'correct phase assemblage without',
     $  ' succeeding. Suggest increasing the',/7x,'value of the ntrymx',
     $  ' variable on the input file.')
        ier = 20
        go to 999
      endif
c
      qmod = .true.
      kpsst = 0
c
c     Modify the indexing of the Jacobian.
c
      call miidxz(ier,iindx1,ipndx1,jpflag,jsflag,kbt,kdim,
     $ kmax,km1,kmt,kx1,kxt,losp,ncmpr,noutpt,npt,nptmax,nstmax,
     $ nttyo,uspec,uzvec1,zvclg1,zvec1)
      if (ier .gt. 0) then
        ier = 8
        go to 999
      endif
c
      write (noutpt,1020)
c
c     Go back and try to make the equilibrium calculation with the
c     new phase assemblage.
c
      go to 100
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  500 continue
c
c     Iteration has terminated without achieving convergence. Attempt to
c     diagnose the cause. It may be that a phase must be removed from
c     the current phase assemblage.
c
      ifail = ifail + 1
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     If there are no mineral phases in the equilibrium system phase
c     assemblage to consider deleting, look for another cause for
c     the failure to converge.
c
      if (kxt .le. kbt) go to 980
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Try to find a phase do delete from the phase assemblage.
c
      npdel = 0
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Trap out-of-range values for the variables associated with the
c     mineral species. An out-of-range value ("-99999.") is likely to
c     be encountered only as a result of having read it from the input
c     file.
c
      do kcol = km1,kxt
c
c       If a value of -99999. is found for the log number of moles
c       variable of any mineral species, go remove the associated
c       phase from the current phase assemblage and try again.
c
        if (zvclg1(kcol) .le. -99999.) then
          npdel = ipndx1(kcol)
          go to 700
        endif
      enddo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if ((qadd .or. qstart) .and. (kxt - km1 + 1).ge.2) then
c
c       Check for the possibility of a mineralogic phase rule violation.
c       The mineralogic phase rule is 'f = c - p', whereas the phase
c       rule of thermodynamics is 'f = c - p + 2'. The mineralogic
c       phase rule applies here because the present computations
c       are made for fixed temperature and pressure. A violation of
c       the mineralogic phase rule is equivalent to a condition of
c       linear dependence among the mineral mass action rows
c       (krow = km1,kxt).
c
        irow1 = km1
        irow2 = kmt
cXX     irow2 = kxt
        jcol1 = 1
        jcol2 = kmt
cXX     jcol2 = kxt
c
c       Calling sequence substitutions:
c         qphruv for qldep
c
        call lindep(aamatr,eps100,irow1,irow2,jcol1,jcol2,kmax,
     $  qphruv)
        if (qphruv) then
          write (noutpt,1400)
          write (nttyo,1400)
 1400     format(/' * Note - (EQ6/eqphas) The current phase assemblage',
     $    ' violates the',/7x,'mineralogic phase rule.')
c
          if (iopt(1).eq.2 .and. .not.qstart .and.
     $      .not.(qcntmp.and.qcnpre)) then
c
c           A mineralogic phase rule violation occurred while running
c           a non-isothermal fluid-centered flow-through open system
c           simulation. The violation may be due to crossing a
c           univariant curve.
c
            write (noutpt,1500)
 1500       format('   --- Possibly attempting to cross a univariant ',
     $      'curve ---',/)
          endif
c
c         Identify a phase which must be deleted in order to avoid the
c         violation. Ordinarily, there would only be one.
c
          call jgibbs(aamatr,afcnst,affp,cdrs,cscale,csts,delvec,
     $    eps100,gmmatr,iindx1,iodb,ipivot,ipndx1,jpflag,kbt,kdim,kmax,
     $    km1,kmt,kx1,kxt,mtb,nbasp,nbtmax,ndrs,ndrsmx,ndrsr,nodbmx,
     $    noutpt,npadd,npdel,nptmax,nstmax,nsts,nstsmx,nstsr,nttyo,
     $    rhsvec,uphase,uspec,xlks)
c
          if (npdel .gt. 0) then
c
c           Go delete the offending phase and try again.
c
            go to 700
          endif
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Try to find a phase to delete using the following two empirical
c     criteria:
c
c        1. The log number of moles of a phase was becoming very
c           negative (or was already very negative if iter .le. 1).
c
c        2. Its derivative with respect to Xi is very negative
c           (must have nord.ge.1 in order to apply this).
c
      call phsdrp(d1zvc1,iindx0,iindx1,iodb,ipndx1,iter,kmax,
     $ km1,km1s,kxt,kxts,nodbmx,nord,noutpt,npadd,npdel,nptmax,ntry,
     $ uphase,zvclgs,zvclg1)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (npdel .le. 0) then
c
c       Couldn't find a phase to delete.
c
        go to 980
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  700 continue
c
c     Have found a phase to delete from the phase assemblage.
c
      if (npdel.eq.npdeli .and. ntry.eq.(ntrydi + 2)) then
c
c       The phase chosen to delete was deleted two steps previously and
c       then added on the previous step.
c
        if (iodb(1) .gt. 0) then
          write (noutpt,1510)
          write (nttyo,1510)
 1510     format(/' * Note - (EQ6/eqphas) Caught in a region of',
     $    /7x,'critical instability in the ES phase assemblage. A',
     $    /7x,"phase can't be deleted because when it is not included",
     $    /7x,'in the ES phase assemblage the system is supersaturated',
     $    /7x,'beyond the allowed tolerance.')
        endif
c
        ier = 40
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Delete a phase from the phase assemblage and try again.
c
      ier = 0
      npdeli = npdel
      ntrydi = ntry
c
      j2 = ilnobl(uphase(npdel))
      write (ux8,'(i7)') npdel
      call lejust(ux8)
      j3 = ilnobl(ux8)
      write (noutpt,1520) uphase(npdel)(1:j2),ux8(1:j3)
 1520 format(/'   The phase to be dropped is ',a,1x,'(',a,')')
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (qbseqc) then
c
c       If any basis switches were done by EQ6/optmzr.f, undo them.
c       Restore the kernel for the equilibrium sytem from the backup.
c
        do kcol = 1,kbt
          nb = iindx1(kcol)
          ns = nbasp(nb)
          ns2 = nbasps(nb)
          if (ns2 .ne. ns) then
            call switch(adhfs,adhfsx,advfs,advfsx,axhfs,axhfsx,
     $      axlks,axlksx,axvfs,axvfsx,cdrs,cdrsx,eps100,ipch,ipchmx,
     $      ipcv,ipcvmx,jflag,jsflag,narn1,narxmx,nbasp,nbaspd,nbaspx,
     $      nb,nbt,nbtmax,nbw,ndrs,ndrsmx,ndrsx,ndrsr,ndrsrx,noutpt,
     $      ns2,nst,nstmax,ntprmx,nttyo,qbassw,qbswok,uspec)
          endif
        enddo
c
        qbseqc = .false.
        nord = nords
c
        km1 = km1s
        kmt = kmts
        kx1 = kx1s
        kxt = kxts
        kdim = kdims
c
        call copyia(nbasps,nbasp,nbt)
        call copyia(iindxs,iindx1,kdim)
        call copyia(ipndxs,ipndx1,kdim)
c
        call copyaa(zvclgs,zvclg1,kdim)
        call copyaa(acflgs,acflg,nst)
c
        xbarwc = xbarws
        xbrwlc = xbrwls
c
c       Reset the ixbasp and cjbasp arrays. The former is a flag
c       array, each member of which denotes whether the
c       thermodynamic activity of the corresponding basis species
c       is defined in terms of molality (= 0) or mole fraction (= 1).
c       The cjbasp array contains any site stoichiometric factors
c       associated with the operational basis species.
c
        call gibasp(cgexj,cjbasp,iern1,ixbasp,jern1,jern2,
     $  jetmax,jgext,narn1,narn2,nbasp,nbt,nbtmax,nern1,nern2,
     $  netmax,nphasx,nstmax)
c
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Reset the log number of moles variables (losp) for all mineral
c     species.
c
      do kcol = km1,kxt
        ns = iindx1(kcol)
        losp(ns) = zvclg1(kcol)
      enddo
c
c     Delete the desired phase.
c
      np = npdel
      jpflag(np) = 0
      loph(np) = -99999.
      moph(np) = 0.
      nr1 = ncmpr(1,np)
      nr2 = ncmpr(2,np)
      do ns = nr1,nr2
        if (jsflag(ns) .le. 0) then
          jsflag(ns) = 0
          losp(ns) = -99999.
          mosp(ns) = 0.
        endif
      enddo
c
c     If the phase to be deleted was the last one added, it is replaced
c     by the phase which had the next highest scaled affinity (as long
c     as such can be found in the limited list maintained for this
c     purpose). Otherwise, the phase to be deleted is just removed.
c
      qadd = .false.
      if (npdel .eq. npaddi) then
        inext = ifail + 1
        if (inext .le. nsspt) then
          npadd = nssp(inext)
c
c         Replace the phase to be deleted with another.
c
          go to 300
        else
c
c         Just delete the phase.
c
c         Note: from here the calculation needs to slide forward
c         in delxi in order to get past a region of computational
c         instability associated with a phase appearance boundary.
c
          go to 400
        endif
      else
c
c       Just delete the phase.
c
        go to 400
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  980 continue
c
c     Test for critical redox instability.
c
      if (qredox) then
        if (idelmx.eq.krdxsp .or. abs(delvec(krdxsp)).ge.2.0) then
          write (noutpt,2000)
          write (nttyo,2000)
 2000     format(/' * Note - (EQ6/eqphas) Caught in a region of',
     $    ' critical redox instability.',/7x,'The value of the redox',
     $    ' master variable (log fO2, Eh, pe-, etc.)',/7x,'is',
     $    ' essentially indifferent to the constraints defining',
     $    ' the state',/7x,'of the ES.')
          ier = 50
          go to 999
        endif
      endif
c
c     Couldn't find a phase to delete from the current phase assemblage.
c     Look for other explanations for the failure to converge.
c
c     Test for critical instability associated with solvent water.
c
      if (kwater .gt. 0) then
        if (idelmx.eq.kwater .or. abs(delvec(kwater)).ge.2.0) then
          write (noutpt,2010)
          write (nttyo,2010)
 2010     format(/' * Note - (EQ6/eqphas) Caught in a region of',
     $    ' critical instability associated',/7x,'with solvent water.')
          if (zvclg1(kwater) .le. -14.0) then
            write (noutpt,2020)
            write (nttyo,2020)
 2020       format(/7x,'There is almost no solvent water present.')
          endif
          ier = 60
          go to 999
        endif
      endif
c
c     Test for instability that appears to be associated with the
c     aqueous activity coefficient model.
c
      if (bgamx.ge.1.0 .or. bfxi.ge.1.0 .or. bsigmm.ge.1.0) then
        write (noutpt,2050) bgamx,fxi,fxio,sigmam,sigmmo
        write (nttyo,2050) bgamx,fxi,fxio,sigmam,sigmmo
 2050   format(/' * Note - (EQ6/eqphas) Caught in a region of',
     $  ' instability that may be',/7x,'associated with the',
     $  ' aqueous activity coefficient model. The max',/7x,'norm of',
     $  ' the activity coefficients for the last change was ',1pe12.5,
     $  /7x,'The last value of the ionic strength in the iteration',
     $  ' was',/7x,e12.5,' molal, the immediately preceding value',
     $  ' was',/7x,e12.5,'. The last value of the sum of the solute',
     $  ' molalities in',/7x,'the iteration was ',e12.5,', the',
     $  ' immediately preceding value was',/7x,e12.5,'.')
        if (zvclg1(kwater) .le. -14.0) then
          write (noutpt,2060)
          write (nttyo,2060)
 2060     format(/7x,'There is almost no solvent water present.')
        endif
        ier = 80
        go to 999
      endif
c
c     Trap out-of-range values for the variables associated with the
c     basis species.
c
      nerr = 0
      do kcol = 1,kbt
        if (zvclg1(kcol) .le. -99999.) then
c
c         Calling sequence substitutions:
c           uzvec1 for unam48
c
          call fmspnx(jlen,uzvec1(kcol),uspn56)
          if (kcol .ne. krdxsp) then
            write (noutpt,2100) uspn56(1:jlen),zvclg1(kcol)
            write (nttyo,2100) uspn56(1:jlen),zvclg1(kcol)
 2100       format(/' * Warning - (EQ6/eqphas) The log number of moles',
     $      ' variable for basis species',/7x,a,' has the out-of-range',
     $      ' value ',1pe12.5,/7x,'after an iteration crash.')
          elseif (qredox) then
            if (krdxsp .eq. ko2gaq) then
              write (noutpt,2110) zvclg1(kcol)
              write (nttyo,2110) zvclg1(kcol)
 2110         format(/' * Warning - (EQ6/eqphas) The log oxygen',
     $        ' fugacity variable has the out-of-range',/7x,'value ',
     $        1pe12.5,' after an iteration crash.')
            elseif (krdxsp .eq. kelect) then
              write (noutpt,2120) zvclg1(kcol)
              write (nttyo,2120) zvclg1(kcol)
 2120         format(/' * Warning - (EQ6/eqphas) The log electron',
     $        ' activity variable has the',/7x,'out-of-range value',
     $        ' ',1pe12.5,' after an iteration crash.')
            else
              write (noutpt,2100) uspn56(1:jlen),zvclg1(kcol)
              write (nttyo,2100) uspn56(1:jlen),zvclg1(kcol)
            endif
          endif
          nerr = nerr + 1
        endif
      enddo
      if (nerr .gt. 0) then
        ier = 110
        go to 999
      endif
c
c     Test simply for almost no solvent water present. This can
c     sometimes be the problem even if some of the other tests
c     applied above fail to detect it.
c
      if (kwater .gt. 0) then
        if (zvclg1(kwater) .le. -14.0) then
          write (noutpt,2030)
          write (nttyo,2030)
 2030     format(/' * Note - (EQ6/eqphas) There appears to be almost',
     $    ' no solvent water present.')
          ier = 70
          go to 999
        endif
      endif
c
c     Okay, can't diagnose what caused the iteration to crash.
c
      write (noutpt,2150)
      write (nttyo,2150)
 2150 format(/' * Note - (EQ6/eqphas) An equilibrium calculation',
     $  ' failed for reasons',/7x,"that couldn't be diagnosed. Will",
     $  ' try to recover.')
      ier = 150
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  990 continue
c
c     Check to see if the new phase assemblage is the same as the old
c     one. If qmod is true at this point, it only means that the
c     putative phase assemblage has changed in the course of the
c     equilibrium calculations. It may, however, have ended up being
c     the original assemblage.
c
      if (qmod) then
        if (km1.eq.km10 .and. kmt.eq.kmt0) then
          if (kx1.eq.kx10 .and. kxt.eq.kxt0) then
            do kcol = km1,kxt
              if (iindx1(kcol) .ne. iindx0(kcol)) go to 995
            enddo
            qmod = .false.
          endif
        endif
      endif
  995 continue
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
