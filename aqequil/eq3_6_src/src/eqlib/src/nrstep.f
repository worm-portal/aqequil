      subroutine nrstep(aamatr,acflg,act,actlg,afcnst,alpha,al10,amtb,
     $ bbig,beta,betamx,betao,betfnc,betmxo,bneg,bpx,btfcnr,cdrs,cdrtw,
     $ cdrw,cegexs,cgexj,cjbasp,cnufac,conc,conclg,cpgexs,csts,coval,
     $ delfnc,delmax,delvco,delvec,dlogxw,egexjc,egexjf,egexs,eh,ehfac,
     $ eps100,fo2,fo2lg,fsort,fugac,fugalg,gmmatr,ibetmx,ibpxmx,idelmx,
     $ iebal,ier,iern1,iern2,ietmax,ifrn1,ifrn2,igas,igstak,iindx1,
     $ ilrn1,ilrn2,imrn1,imrn2,iodb,ipivot,ipndx1,irdxc3,istack,iter,
     $ itermx,ixbasp,ixrn1,ixrn2,jcsort,jern1,jern2,jetmax,jflag,jgext,
     $ jgsort,jgstak,jjsort,jpflag,jsflag,jsitex,jssort,jstack,kbt,
     $ kction,kdim,kelect,khydr,kmax,kmt,km1,ko2gaq,kwater,kxt,kx1,
     $ loph,losp,lsort,mgext,moph,mosp,mrgexs,mtb,narn1,narn2,nbasp,
     $ nbt,nbtmax,nbw,ncmpr,ncosp,ndrs,ndrsmx,ndrsr,negbfc,negdfc,
     $ negxfc,nelect,nern1,nern2,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,
     $ ngtmax,nhydr,nodbmx,noibfc,noutpt,no2gaq,npconv,nphasx,npobfc,
     $ npodfc,npt,nptmax,nredox,nst,nstmax,nsts,nstsmx,nstsr,ntfx,
     $ ntfxmx,ntfxt,nttyo,nxtmax,omega,omeglg,press,qcacf,qcbeta,
     $ qredox,qxbarw,q6mode,rhsvec,rlxfac,screwd,screwn,sigmam,sigmmo,
     $ tfx,ubetmx,ubbig,ubneg,ugexj,ugexmo,uphase,uspec,uzvec1,weight,
     $ xbar,xbarlg,xbarw,xbarwc,xbrwlc,xbrwlg,xlke,xlks,zchar,zgexj,
     $ zvclg1,zvec1)
c
c     This subroutine performs one Newton-Raphson step.
c
c     This subroutine is called by:
c
c       EQLIB/newton.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       iodb   = array of debugging print options
c       ipivot = the pivot vector
c       screwd = under-relaxation control parameter. It is used to
c                  reduce the magnitude of the delvec vector, if
c                  necessary, so that the magnitude of the largest
c                  element of that vector does not exceed screwd
c       screwn = under-relaxation control parameter
c       itermx = maximum number of iterations
c       kdim   = dimension of aamatr
c       kmax   = maximum dimension of aamatr
c       qcbeta = true if the convergence criterion on betamx was met
c       qcacf  = true if the convergence criterion on activity
c                  coefficients was met
c
c     Principal output:
c
c       aamatr = Jacobian matrix
c       gmmatr = copy of aamatr
c       rhsvec = right hand side vector
c       alpha  = residual function array
c       idelmx = kcol index corresponding to delmax
c       delfnc = convergence function, defined by reference to delmax
c       betfnc = convergence function, defined by reference to betamx
c       uzvec1 = name array corresponding to zvclg1
c       ier    = error flag:
c                  =  0  Okay
c                  =  1  Encountered a zero matrix
c                  =  2  Encountered a non-zero, computationally
c                          singular matrix
c                  =  3  Iteration was diverging
c                  =  4  Hit the maximum number of iterations (itermx)
c
c     Principal input/output:
c
c       iter   = the number of Newton-Raphson iterations
c       negdfc = the number of successive iterations that delmax
c                  has been greater than zero and the corresponding
c                  convergence function delfnc has been less than
c                  or equal to zero
c       negbfc = the number of successive iterations that betamx
c                  has been greater than zero and the corresponding
c                  convergence function betfnc has been less than
c                  or equal to zero
c       negxfc = the number of successive iterations that either
c                  delmax has been greater than zero while delfnc
c                  has been less than or equal to zero, or betamx
c                  has been greater than zero while betfnc has
c                  been less than or equal to zero
c       noibfc = the number of successive iterations that betamx
c                  has been greater than zero and the corresponding
c                  convergence function betfnc has been less than
c                  1.e-6
c       npconv = number of successive steps in which the convergence
c                  criteria on the residual norm has been satisfied
c       npodfc = the number of successive iterations that the
c                  convergence function delfnc has been greater
c                  than zero
c       npobfc = the number of successive iterations that the
c                  convergence function betfnc has been greater
c                  than zero
c       zvclg1 = the 'log z' array, the array corrected by
c                  Newton-Raphson iteration
c       beta   = normalized residual function array
c       delvec = correction array
c       betao  = old beta array
c       delvco = old delvec array
c       betamx = max norm of the beta array
c       delmax = max norm of the delvec array
c       rlxfac = under-relaxation factor
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer ibpxmx,ietmax,jetmax,kmax,nbtmax,ndrsmx,netmax,ngtmax,
     $ nodbmx,nptmax,nstmax,nstsmx,ntfxmx,nxtmax
c
      integer noutpt,nttyo
c
      integer igstak(ngtmax),iindx1(kmax),iodb(nodbmx),ipivot(kmax),
     $ ipndx1(kmax),istack(nstmax),ixbasp(nbtmax),jcsort(nstmax),
     $ jern1(jetmax,netmax),jern2(jetmax,netmax),jflag(nstmax),
     $ jgext(netmax),jgsort(ngtmax),jgstak(ngtmax),jjsort(nstmax),
     $ jpflag(nptmax),jsflag(nstmax),jsitex(nstmax),jssort(nstmax),
     $ jstack(nstmax),kction(nbtmax),nbasp(nbtmax),ncmpr(2,nptmax),
     $ ncosp(nbtmax),ndrs(ndrsmx),ndrsr(2,nstmax),
     $ ngexsa(ietmax,jetmax,netmax),ngext(jetmax,netmax),
     $ nphasx(nstmax),nsts(nstsmx),nstsr(2,nstmax),ntfx(ntfxmx)
c
      integer ibetmx,idelmx,iebal,ier,iern1,iern2,ifrn1,ifrn2,igas,
     $ ilrn1,ilrn2,imrn1,imrn2,irdxc3,iter,itermx,ixrn1,ixrn2,kbt,kdim,
     $ kelect,khydr,kmt,km1,ko2gaq,kwater,kxt,kx1,narn1,narn2,nbt,nbw,
     $ negbfc,negdfc,negxfc,nelect,nern1,nern2,ngrn1,ngrn2,ngt,nhydr,
     $ noibfc,no2gaq,npconv,npobfc,npodfc,npt,nredox,nst,ntfxt
c
      logical qcacf,qcbeta,qredox,qxbarw,q6mode
c
      character*48 uspec(nstmax),uzvec1(kmax)
      character*48 ubbig,ubneg,ubetmx
      character*24 ugexmo(netmax),uphase(nptmax)
      character*8 ugexj(jetmax,netmax)
c
      real*8 aamatr(kmax,kmax),acflg(nstmax),act(nstmax),actlg(nstmax),
     $ alpha(kmax),amtb(nbtmax),beta(kmax),betao(kmax),
     $ bpx(ibpxmx,nxtmax),cdrs(ndrsmx),cdrtw(nstmax),cdrw(nstmax),
     $ cegexs(ietmax,jetmax,netmax),cgexj(jetmax,netmax),cjbasp(nbtmax),
     $ conc(nstmax),conclg(nstmax),cnufac(nstmax),coval(nbtmax),
     $ cpgexs(ietmax,jetmax,netmax),csts(nstsmx),delvco(kmax),
     $ delvec(kmax),dlogxw(nbtmax)
c
      real*8 egexjc(jetmax,netmax),egexjf(jetmax,netmax),
     $ egexs(ietmax,jetmax,netmax),fsort(ngtmax),fugac(ngtmax),
     $ fugalg(ngtmax),gmmatr(kmax,kmax),loph(nptmax),losp(nstmax),
     $ lsort(nstmax),mgext(jetmax,netmax),moph(nptmax),mosp(nstmax),
     $ mrgexs(ietmax,jetmax,netmax),mtb(nbtmax),rhsvec(kmax),
     $ tfx(ntfxmx),weight(nstmax),xbar(nstmax),xbarlg(nstmax),
     $ xlks(nstmax),zchar(nstmax),zgexj(jetmax,netmax),zvclg1(kmax),
     $ zvec1(kmax)
c
      real*8 afcnst,al10,bbig,betamx,betfnc,betmxo,bneg,btfcnr,
     $ delfnc,delmax,eh,ehfac,eps100,fo2,fo2lg,omega,omeglg,press,
     $ rlxfac,screwd,screwn,sigmam,sigmmo,xbarw,xbarwc,xbrwlc,
     $ xbrwlg,xlke
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer i,jlen,kcol,krow,ncut
c
      integer iarmxn
c
      logical qconc1,qconc2,qconc3,qpr
c
      character*56 uspn56
c
      real*8 delmxo,divfmx,divfnc,fx,rdx,sxu
c
c-----------------------------------------------------------------------
c
      data qpr    /.false./
c
c-----------------------------------------------------------------------
c
c     Save the current values of the beta and delvec arrays for use in
c     under-relaxation schemes.
c
      call copyaa(beta,betao,kdim)
      call copyaa(delvec,delvco,kdim)
      delmxo = delmax
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Check to see if iteration should be stopped because there is
c     little probability of further improvement. If any of the
c     following conditions cause iteration to stop, the result is
c     considered acceptable if it satisfies the convergence criteria
c     on residual functions, but not on the delvec vector.
c
c     Note: betfnc is the convergence function based on the beta
c     residual vector, and delfnc is that based on the delvec
c     correction vector. If convergence is occurring, each function
c     should approach unity (the maximum possible value).
c
c     Update counters. These all count the number of times in a row
c     that a condition has been satisfied.
c
      if (qcbeta .and. qcacf) then
        npconv = npconv + 1
      else
        npconv = 0
      endif
c
      if (delmax.gt.0. .and. delfnc.le.0.) then
        negdfc = negdfc + 1
      else
        negdfc = 0
      endif
c
      if (betamx.gt.0. .and. betfnc.le.0.) then
        negbfc = negbfc + 1
      else
        negbfc = 0
      endif
c
      if ((betamx.gt.0. .and. betfnc.le.0.) .or.
     $  (delmax.gt.0. .and. delfnc.le.0.)) then
        negxfc = negxfc + 1
      else
        negxfc = 0
      endif
c
      if (delfnc .gt. 0.) then
        npodfc = npodfc + 1
      else
        npodfc = 0
      endif
c
      if (betfnc .gt. 0.) then
        npobfc = npobfc + 1
      else
        npobfc = 0
      endif
c
      if (betamx.gt.0. .and. betfnc.le.1.e-6) then
        noibfc = noibfc + 1
      else
        noibfc = 0
      endif
c
c     Stop iteration if the convergence criteria on residual functions
c     (derived from the beta vector) have been satisfied for four
c     consecutive iterations.
c
      if (npconv .ge. 4) go to 300
c
c     Stop iteration if the behavior of delfnc indicates that the
c     iteration is diverging, and this indication is not contradicted
c     by the behavior of betfnc. Technically, delfnc is negative or
c     zero for twelve consecutive iterations and betfnc has not been
c     positive for three or more consecutive iterations.
c
      if (negdfc.ge.12 .and. npobfc.lt.3) go to 300
c
c     Stop iteration if the behavior of betfnc indicates that the
c     iteration is diverging, and this indication is not contradicted
c     by the behavior of delfnc. Technically, betfnc is negative or
c     zero for twelve consecutive iterations and delfnc has not been
c     positive for three or more consecutive iterations.
c
      if (negbfc.ge.12 .and. npodfc.lt.3) go to 300
c
c     Stop iteration if the behavior of delfnc and betfnc together
c     suggests that the iteration is diverging. Technically, either
c     delfnc or betfnc is negative or zero for twelve consecutive
c     iterations, and iter is greater than or equal to 40.
c
      if (negxfc.ge.12 .and. iter.ge.40) go to 300
c
c     Stop iteration if the behavior of betfnc indicates that no
c     significant improvement is occurring. Technically, betfnc is
c     less than or equal to 1.e-6 for twelve consecutive iterations,
c     and iter is greater than or equal to 20.
c
      if (noibfc.ge.12 .and. iter.ge.20) go to 300
c
c     Stop iteration if the maximum number of iterations have
c     been done.
c
      if (iter .ge. itermx) then
        write (noutpt,1000) itermx
        write (nttyo,1000) itermx
 1000   format(/' * Note - (EQLIB/nrstep) Have completed ',i4,' hybrid',
     $  /7x,'Newton-Raphson iterations. This is the maximum number',
     $  ' permitted.')
        ier = 4
        go to 999
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Do another Newton-Raphson iteration.
c
      iter = iter + 1
c
c     Calculate the Jacobian matrix (aamatr).
c
      call matrix(aamatr,al10,bpx,cdrs,cdrtw,cdrw,cjbasp,
     $ cnufac,conc,csts,dlogxw,eps100,ibpxmx,iebal,iern1,ietmax,
     $ iindx1,ipndx1,irdxc3,ixbasp,ixrn1,jcsort,jern1,jern2,jetmax,
     $ jflag,jjsort,jsitex,kbt,kction,kdim,kelect,khydr,kmax,kmt,
     $ km1,ko2gaq,kwater,kxt,kx1,mosp,narn1,narn2,nbasp,nbt,nbtmax,
     $ nbw,ncosp,ndrs,ndrsmx,ndrsr,nelect,nern1,nern2,netmax,
     $ noutpt,no2gaq,nphasx,nredox,nst,nstmax,nsts,nstsmx,nstsr,
     $ ntfx,ntfxmx,ntfxt,nttyo,nxtmax,omega,qredox,q6mode,tfx,
     $ ugexmo,uspec,weight,xbar,xbarw,xbarwc,zchar)
c
c      Get the right-hand-side vector (rhsvec = -alpha).
c
       do krow = 1,kdim
         rhsvec(krow) = -alpha(krow)
       enddo
c
      if (iodb(4) .ge. 3) then
        write (noutpt,1010)
 1010   format(/16x,'--- The aamatr array, with rhsvec ---',/)
        do krow = 1,kdim
c
c         Calling sequence substitutions:
c           uzvec1(krow) for unam48
c
          call fmspnx(jlen,uzvec1(krow),uspn56)
          write (noutpt,1020) krow,uspn56(1:jlen)
 1020     format(1x,i3,2x,a)
          write (noutpt,1030) (aamatr(krow,i),i = 1,kdim),rhsvec(krow)
 1030     format(( 6(2x,g10.3)) )
          write(noutpt,1040)
 1040     format(' ')
        enddo
      endif
c
c     Solve the Jacobian system. If EQLIBU/msolvr.f can't solve the
c     matrix, it is because the matrix is either zero (ier = 1) or
c     non-zero, but computationally singular (ier = 2).
c
      call msolvr(aamatr,delvec,gmmatr,ier,ipivot,kdim,kmax,
     $ noutpt,nttyo,qpr,rhsvec)
c
      if (ier .gt. 0) go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Find the max norm of the correction vector (delvec).
c     In the limit of the solution, the (unrelaxed) correction
c     term bounds the error in the variables being calculated.
c
      idelmx = iarmxn(delvec,kdim)
      delmax = 0.
      if (idelmx .gt. 0) delmax = abs(delvec(idelmx))
c
c     Calculate the delvec improvement function (delfnc).
c
      delfnc = 0.
      if (delmxo .gt. 0.) delfnc = (delmxo -delmax)/delmxo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Apply under-relaxation according to some simple criteria.
c     The under-relaxation at this point is done before recomputing
c     the residual functions. A subsequent under-relaxation
c     will be done to insure that the residual functions do not
c     increase too much.
c
      rlxfac = 1.0
c
c     Limit the magnitude of change per iteration.
c
      if (delmax .gt. screwd) rlxfac = screwd/delmax
c
c     If within the first eight iterations and the residual norm
c     (betamx) is large, go slow.
c
      if (iter.le.8 .and. betamx.gt.screwn) then
        rlxfac = 0.5*rlxfac
      endif
c
      qconc1 = sigmam.ge.2.0 .or. sigmmo.ge.2.0
      qconc2 = sigmam.ge.8.0 .or. sigmmo.ge.8.0
      qconc3 = sigmam.ge.12.0 .or. sigmmo.ge.12.0
c
      if (qconc1) then
c
c       Force some minimum under-relaxation for modestly concentrated
c       solutions.
c
        if (iter .le. 15) then
          sxu = 0.50
          rlxfac = min(sxu,rlxfac)
        endif
      endif
c
      if (qconc2) then
c
c       Force some minimum under-relaxation for more concentrated
c       solutions.
c
        if (iter .le. 15) then
          sxu = 0.35
          rlxfac = min(sxu,rlxfac)
        endif
      endif
c
      if (qconc3) then
c
c       Force some minimum under-relaxation for highly concentrated
c       solutions.
c
        if (iter .le. 15) then
          sxu = 0.25
          rlxfac = min(sxu,rlxfac)
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Apply correction terms. The label below is a return point
c     under-relaxation which restricts the growth of the residual
c     functions.
c
      ncut = 0
  140 continue
c
      if (iodb(4) .ge. 2) write (noutpt,1050) ncut
 1050 format(/3x,'ncut= ',i3,/)
c
      do kcol = 1,kdim
        zvclg1(kcol) = zvclg1(kcol) + rlxfac*delvec(kcol)
      enddo
c
      if (iodb(4) .ge. 2) then
        write (noutpt,1060)
 1060   format(/8x,'Name',38x,'del',9x,'zvclg1'/)
        do kcol = 1,kdim
c
c         Calling sequence substitutions:
c           uzvec1(kcol) for unam48
c
          call fmspnx(jlen,uzvec1(kcol),uspn56)
          jlen = min(jlen,38)
          rdx = rlxfac*delvec(kcol)
          write (noutpt,1062) kcol,uspn56(1:jlen),rdx,zvclg1(kcol)
 1062     format(1x,i3,2x,a,t46,1pe12.5,2x,e12.5)
        enddo
        write (noutpt,1064)
 1064   format(1x)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Expand the system and recalculate the Newton-Raphson residual
c     functions. It is important here to keep the activity coefficients
c     constant,that is, not to update them until all under-relaxation
c     has been completed.
c
c     Recalculate the concentrations, etc., of dependent species.
c
      call ncmpex(acflg,act,actlg,cdrs,cegexs,cgexj,conc,
     $ conclg,cpgexs,egexjc,egexjf,egexs,eps100,fo2,fo2lg,fsort,
     $ fugac,fugalg,iern1,iern2,ietmax,ifrn1,ifrn2,igas,igstak,
     $ iindx1,ilrn1,ilrn2,imrn1,imrn2,istack,ixrn1,ixrn2,jcsort,
     $ jern1,jern2,jetmax,jflag,jgext,jgsort,jgstak,jjsort,jpflag,
     $ jsflag,jsitex,jssort,jstack,kbt,kdim,kelect,kmax,km1,ko2gaq,
     $ kwater,kxt,loph,losp,lsort,mgext,mrgexs,mtb,moph,mosp,narn1,
     $ narn2,nbasp,nbt,nbtmax,ncmpr,ndrs,ndrsmx,ndrsr,nelect,nern1,
     $ nern2,netmax,ngexsa,ngext,ngrn1,ngrn2,ngt,ngtmax,noutpt,
     $ no2gaq,nphasx,npt,nptmax,nst,nstmax,nttyo,omega,omeglg,
     $ press,qxbarw,q6mode,ugexj,ugexmo,uphase,uspec,xbar,xbarlg,
     $ xbarw,xbarwc,xbrwlc,xbrwlg,xlks,zchar,zgexj,zvclg1,zvec1)
c
      call betas(acflg,actlg,afcnst,alpha,amtb,bbig,beta,
     $ betamx,bneg,cdrs,conc,conclg,coval,csts,eh,ehfac,fo2lg,
     $ ibetmx,iebal,iindx1,irdxc3,jcsort,jflag,jsflag,jssort,kbt,
     $ kdim,kelect,khydr,kmax,km1,ko2gaq,kwater,kxt,mtb,mosp,
     $ narn1,narn2,nbasp,nbtmax,ncosp,ndrs,ndrsmx,ndrsr,nelect,
     $ nern1,nern2,nhydr,noutpt,no2gaq,nredox,nst,nstmax,nsts,
     $ nstsmx,nstsr,ntfx,ntfxmx,ntfxt,nttyo,omega,qredox,q6mode,
     $ tfx,ubbig,ubneg,ubetmx,uspec,uzvec1,weight,xbrwlg,xlke,
     $ xlks,zchar)
c
      if (iodb(4) .ge. 2) then
        write (noutpt,1080)
 1080   format(/8x,'Name',37x,'Beta',/)
        do kcol = 1,kdim
c
c         Calling sequence substitutions:
c           uzvec1(kcol) for unam48
c
          call fmspnx(jlen,uzvec1(kcol),uspn56)
          jlen = min(jlen,38)
          write (noutpt,1085) kcol,uspn56(1:jlen),beta(kcol)
 1085     format(1x,i3,2x,a,t46,1pe12.5)
        enddo
        write (noutpt,1090)
 1090   format(/1x)
      endif
c
      btfcnr = 0.
      if (betmxo .gt. 0.) btfcnr = (betmxo - betamx)/betmxo
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Use under-relaxation to control residual growth, if any.
c     Terminate iteration if residual tracking indicates that the
c     iteration is diverging.
c
      if (iter .le. 12) then
        divfmx = 0.
        do kcol = 1,kdim
          divfnc = abs(beta(kcol)) - betao(kcol)
          if (divfmx .gt. divfnc) divfmx = divfnc
        enddo
        if (divfmx .gt. 0.10) then
c
          if (ncut .le. 7) then
            ncut = ncut + 1
            if (iodb(4) .ge. 2) write (noutpt,1100)
 1100       format('  Residual tracking requires under-relaxation.')
            fx = 1./rlxfac
            do kcol = 1,kdim
              zvclg1(kcol) = zvclg1(kcol) - fx*delvec(kcol)
            enddo
            rlxfac = 0.25*rlxfac
            go to 140
          elseif (iter .ge. 8) then
            write (noutpt,1110)
 1110       format('  Residual tracking indicates that iteration',
     $      ' is diverging.')
            go to 300
          endif
c
        endif
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (iodb(4) .ge. 2) write (noutpt,1120)
 1120 format(' A Newton-Raphson correction has been completed.',/)
      go to 999
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  300 continue
c
c     The code detected that iteration was diverging.
c
      ier = 3
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
  999 continue
      end
