      subroutine echox(azero,cdrs,covali,eh,fo2lg,iebal3,iodb,iopg,
     $ iopr,iopt,irdxc3,itdsf3,itermx,ixrn1,ixrn2,jflag,jflgi,jpres3,
     $ narn1,narn2,nat,nata,natmax,nbt,nbta,nbti,nbtmax,nct,ncta,
     $ nctmax,ncmpr,ncosp,ndecsp,ndrs,ndrsmx,ndrsr,nelect,ngt,ngta,
     $ ngtmax,nhydr,nhydx,njfmax,nlt,nlta,nltmax,nmt,nmta,nmtmax,
     $ nodbmx,nopgmx,noprmx,noptmx,noutpt,no2gaq,npt,npta,nptmax,
     $ nredox,nst,nsta,nstmax,nxt,nxta,nxti,nxtmax,pe,presg,press,
     $ rho,scamas,tdspkg,tdspl,tempc,tolbt,toldl,tolspf,uactop,
     $ ucospi,uebal,ujflls,uphase,uredox,uspec,uspeci,xbar)
c
c     This subroutine echoes the problem input defining the
c     compositional constraints on the aqueous solution being modeled
c     and the options and tolerances selected.
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
      integer natmax,nbtmax,nctmax,ndrsmx,ngtmax,njfmax,nltmax,
     $ nmtmax,nodbmx,nopgmx,noprmx,noptmx,nptmax,nstmax,nxtmax
c
      integer noutpt
c
      integer iodb(nodbmx),iopg(nopgmx),iopr(noprmx),iopt(noptmx),
     $ jflag(nstmax),jflgi(nbtmax),ncmpr(2,nptmax),ncosp(nbtmax),
     $ ndecsp(nbtmax),ndrs(ndrsmx),ndrsr(2,nstmax)
c
      integer iebal3,irdxc3,itdsf3,itermx,ixrn1,ixrn2,jpres3,narn1,
     $ narn2,nat,nata,nbt,nbta,nbti,nct,ncta,nelect,ngt,ngta,nhydr,
     $ nhydx,nlt,nlta,nmt,nmta,no2gaq,npt,npta,nredox,nst,nsta,nxt,
     $ nxta,nxti
c
      character*48 ucospi(nbtmax),uspec(nstmax),uspeci(nbtmax)
      character*32 ujflls(0:njfmax)
      character*32 uactop
      character*24 uphase(nptmax)
      character*24 uebal,uredox
c
      real(8) azero(natmax),cdrs(ndrsmx),covali(nbtmax),xbar(nstmax)
c
      real(8) dp,eh,fo2lg,pe,presg,press,rho,scamas,tdspkg,tdspl,tempc,
     $ tolbt,toldl,tolspf
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jfl,j2,j3,j4,n,na,nbi,nb1,np,nr1,nr2,ns,nss
c
      integer ilnobl
c
      character*24 ux24
c
c-----------------------------------------------------------------------
c
      j2 = ilnobl(uactop)
      write (noutpt,1000) uactop(1:j2)
 1000 format(/' The activity coefficients of aqueous species will be',
     $ /' calculated using ',a,'.'/)
c
      if (iopr(3) .gt. 1) then
        write (noutpt,1010)
 1010   format(/7x,'Species',17x,'HC Diameter',/)
        do ns = narn1 + 1,narn2
          na = ns - narn1 + 1
          write (noutpt,1020) uspec(ns),azero(na)
 1020     format(5x,a24,3x,f7.3)
        enddo
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1030) tempc
 1030 format(/' Temperature= ',f6.2,' C',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1032) jpres3
 1032 format(/' jpres3=  ',i3,' (Pressure option switch)')
      if (jpres3 .eq. 0) then
        write (noutpt,1034) press
 1034   format(/'   Pressure= ',1pg12.5,' bars (data file reference',
     $  ' curve value)',/)
      elseif (jpres3 .eq. 1) then
        dp = press - presg
        write (noutpt,1036) press,presg,dp
 1036   format(/'   Pressure= ',1pg12.5,' bars (1.013-bar/steam-',
     $  'saturation curve value)',
     $  /'   Data file reference curve pressure= ',g12.5,' bars',
     $  /'   Pressure difference= ',g12.5,' bars',/)
      elseif (jpres3 .eq. 2) then
        dp = press - presg
        write (noutpt,1038) press,presg,dp
 1038   format(/'   Pressure= ',1pg12.5,' bars (specified value)',
     $  /'   Data file reference curve pressure= ',g12.5,' bars',
     $  /'   Pressure difference= ',g12.5,' bars',/)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Print a table showing for phases, species, and groups thereof,
c     the number of each of entity on the data base, the number the
c     software is dimensioned for, and the number appearing in the
c     current problem.
c
      call prtntt(nat,nata,natmax,nbt,nbta,nbtmax,nct,ncta,
     $ nctmax,ngt,ngta,ngtmax,nlt,nlta,nltmax,nmt,nmta,nmtmax,noutpt,
     $ npt,npta,nptmax,nst,nsta,nstmax,nxt,nxta,nxtmax)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1040) (iopt(n), n = 1,10)
 1040 format(/' iopt(1)=  ',i2,' (Used only by EQ6)',
     $ /' iopt(2)=  ',i2,' (Used only by EQ6)',
     $ /' iopt(3)=  ',i2,' (Used only by EQ6)',
     $ /' iopt(4)=  ',i2,' (Solid solutions)',
     $ /' iopt(5)=  ',i2,' (Used only by EQ6)',
     $ /' iopt(6)=  ',i2,' (Used only by EQ6)',
     $ /' iopt(7)=  ',i2,' (Not used)',
     $ /' iopt(8)=  ',i2,' (Not used)',
     $ /' iopt(9)=  ',i2,' (Not used)',
     $ /' iopt(10)= ',i2,' (Not used)')
      write (noutpt,1045) (iopt(n), n = 11,18)
 1045 format(' iopt(11)= ',i2,' (Auto basis switching, in',
     $ ' pre-Newton-Raphson optimization)',
     $ /' iopt(12)= ',i2,' (Used only by EQ6)',
     $ /' iopt(13)= ',i2,' (Not used)',
     $ /' iopt(14)= ',i2,' (Not used)',
     $ /' iopt(15)= ',i2,' (Used only by EQ6)',
     $ /' iopt(16)= ',i2,' (Not used)',
     $ /' iopt(17)= ',i2,' (pickup file options)',
     $ /' iopt(18)= ',i2,' (Used only by EQ6)',
     $ /' iopt(19)= ',i2,' (Advanced EQ3NR pickup file options)',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1050) (iopg(n), n = 1,2)
 1050 format(/'   iopg(1)=  ',i2,' (Aqueous species activity',
     $ ' coefficient model)',
     $ /'   iopg(2)=  ',i2,' (pH scale)',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1060) (iopr(n), n = 1,17)
 1060 format(/' iopr(1)=  ',i2,' (List all species)',
     $ /' iopr(2)=  ',i2,' (List all reactions)',
     $ /' iopr(3)=  ',i2,' (List HC diameters)',
     $ /' iopr(4)=  ',i2,' (Aqueous species concentration print',
     $ ' cut-off)',
     $ /' iopr(5)=  ',i2,' (Ion/H+ activity ratios)',
     $ /' iopr(6)=  ',i2,' (Mass balance percentages)',
     $ /' iopr(7)=  ',i2,' (Affinity print cut-off)',
     $ /' iopr(8)=  ',i2,' (Fugacities)',
     $ /' iopr(9)=  ',i2,' (Mean molal activity coefficients)',
     $ /' iopr(10)= ',i2,' (Pitzer coefficients tabulation)',
     $ /' iopr(11)= ',i2,' (Not used)',
     $ /' iopr(12)= ',i2,' (Not used)',
     $ /' iopr(13)= ',i2,' (Not used)',
     $ /' iopr(14)= ',i2,' (Not used)',
     $ /' iopr(15)= ',i2,' (Not used)',
     $ /' iopr(16)= ',i2,' (Not used)',
     $ /' iopr(17)= ',i2,' (pickup file format)',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1070) (iodb(n), n = 1,7)
 1070 format(/'   iodb(1)=  ',i2,' (General diagnostics)',
     $ /'   iodb(2)=  ',i2,' (Used only by EQ6)',
     $ /'   iodb(3)=  ',i2,' (pre-Newton-Raphson optimization',
     $ ' iterations)',
     $ /'   iodb(4)=  ',i2,' (Newton-Raphson iterations)',
     $ /'   iodb(5)=  ',i2,' (Used only by EQ6)',
     $ /'   iodb(6)=  ',i2,' (Hypothetical affinity iterations)',
     $ /'   iodb(7)=  ',i2,' (Used only by EQ6)',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      jfl = jflag(no2gaq)
      write (noutpt,1080) irdxc3
 1080 format(/' irdxc3=  ',i3,' (Default redox constraint switch)')
      if (irdxc3 .eq. 1) then
        if (uredox(1:7).ne.'O2(aq) ' .and.
     $    uredox(1:7).ne.'H2(aq) ') then
          nr1 = ndrsr(1,nredox)
          nr2 = ndrsr(2,nredox)
          do n = nr1,nr2
            nss = ndrs(n)
            if (nss.ne.nredox .and. nss.ne.narn1 .and. nss.ne.nhydr
     $      .and. nss.ne.nhydx .and. nss.ne.no2gaq .and. nss.ne.nelect)
     $      go to 100
          enddo
  100     continue
        else
          nss = narn1
        endif
        j2 = ilnobl(uredox)
        j3 = ilnobl(uspec(nss)(1:24))
        write (noutpt,1100) uredox(1:j2),uspec(nss)(1:j3)
 1100   format(/'   The default redox state is constrained by the',
     $   /'   ',a,'/',a,' couple.',/)
      elseif (irdxc3 .eq. 0) then
        write (noutpt,1110) fo2lg
 1110   format(/'   The default redox state is constrained by',
     $   ' Log fO2 = ',g12.5,' (log bars).',/)
      elseif (irdxc3 .eq. -1) then
        write (noutpt,1120) eh
 1120   format(/'   The default redox state is constrained by',
     $  ' Eh = ',f8.5,' volts.',/)
      elseif (irdxc3 .eq. -2) then
        write (noutpt,1130) pe
 1130   format(/'   The default redox state is constrained by',
     $  ' pe- = ',1pe12.5,'.',/)
      elseif (irdxc3 .eq. -3) then
        write (noutpt,1140)
 1140   format(/'   The default redox state is controlled by a',
     $  /' heterogeneous reaction (see below).',/)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1142) iebal3
 1142 format(/' iebal3=  ',i3,' (Electrical balancing option switch)')
      if (iebal3 .le. 0) then
        write (noutpt,1144)
 1144   format(/'   No electrical balancing adjustment will be made.',
     $  /'   The imbalance will be calculated.',/)
      elseif (iebal3 .eq. 1) then
        j2 = ilnobl(uebal)
        write (noutpt,1146) uebal(1:j2)
 1146   format(/'   The species ',a,' will be adjusted to',
     $  /'   achieve electrical balance.',/)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1200) rho
 1200 format(/' Solution density = ',f8.5,' g/ml',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1202) itdsf3
 1202 format(/' itdsf3=  ',i3,' (Total dissolved solutes option',
     $ ' switch)')
      if (itdsf3 .le. 0) then
        write (noutpt,1204) tdspkg
 1204   format( /'   Total dissolved salts = ',f10.2,' mg/kg.sol',/)
      else
        write (noutpt,1206) tdspl
 1206  format(/'   Total dissolved salts = ',f10.2,' mg/L',/)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1210) tolbt,toldl,tolspf
 1210 format(/' tolbt  = ',1pe12.5,' (convergence tolerance on',
     $ ' residual functions)',/' toldl  = ',e12.5,' (convergence',
     $ ' tolerance on correction terms)',/' tolspf = ',e12.5,
     $ ' (saturation print flag tolerance, does not affect',/25x,
     $ 'convergence)',/)
c
      write (noutpt,1220) itermx
 1220 format(/' itermx = ',i3,' (maximum number of iterations)',/)
c
      write (noutpt,1230) scamas
 1230 format(/' scamas = ',1pe12.5,' (scale factor for aqueous',
     $ ' solution',/25x,'mass written on the pickup file)',/)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1300)
 1300 format(/21x,'--- Original Input Constraints ---',
     $ //5x,'Species',20x,'coval   jflag   Type of Input',/)
c
      do nbi = 1,nbti
        jfl = jflgi(nbi)
        nb1 = ndecsp(nbi)
        if (nb1 .gt. 0) then
          j2 = ilnobl(ujflls(jfl))
c
          if (jfl.eq.17 .or. jfl.eq.18) then
            write (noutpt,1310) uspeci(nbi),covali(nbi),jfl,
     $      ujflls(jfl)(1:j2)
 1310       format(2x,a24,2x,1pe12.5,2x,i2,2x,a)
            j3 = ilnobl(ucospi(nbi)(1:24))
            write (noutpt,1315) ucospi(nbi)(1:j3)
 1315       format(39x,'Counterion= ',a)
          elseif (jfl .eq. 21) then
            write (noutpt,1310) uspeci(nbi),covali(nbi),jfl,
     $      ujflls(jfl)(1:j2)
          elseif (jfl .eq. 25) then
            write (noutpt,1310) uspeci(nbi),covali(nbi),jfl,
     $      ujflls(jfl)(1:j2)
            ns = ncosp(nb1)
            j3 = ilnobl(ucospi(nbi)(1:24))
            ux24 = ucospi(nbi)(25:48)
            j4 = ilnobl(ux24)
            if (j4 .le. 0) then
              ux24 = ucospi(nbi)(1:24)
              j4 = j3
            endif
            write (noutpt,1320) ucospi(nbi)(1:j3),ux24(1:j4)
 1320       format(46x,'Species= ',a,/48x,'Phase= ',a)
c
c           Calling sequence substitutions:
c             noutpt for nf
c
            call prreac(cdrs,ndrs,ndrsmx,ndrsr,noutpt,ns,nstmax,uspec)
            write (noutpt,1322)
 1322       format(1x)
          elseif (jfl.eq.27 .or. jfl.eq.30) then
            write (noutpt,1323) uspeci(nbi),jfl,ujflls(jfl)(1:j2)
 1323       format(2x,a24,16x,i2,2x,a)
          elseif (uspeci(nbi)(1:6) .eq. 'O2(g) ') then
            write (noutpt,1325) uspeci(nbi),covali(nbi),jfl
 1325       format(2x,a24,2x,1pe12.5,2x,i2,2x,'Log fO2')
          else
            if (jfl .ge. 0) then
              write (noutpt,1310) uspeci(nbi),covali(nbi),jfl,
     $        ujflls(jfl)(1:j2)
            else
              write (noutpt,1327) uspeci(nbi)
 1327         format(2x,a24,16x,'Not present in the model')
            endif
          endif
c
        endif
      enddo
c
      write (noutpt,1340)
 1340 format(1x)
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      if (nxti .gt. 0) then
        write (noutpt,1600)
 1600   format(/21x,'--- Input Solid Solution Compositions ---',/)
        do np = ixrn1,ixrn2
          nr1 = ncmpr(1,np)
          nr2 = ncmpr(2,np)
          do ns = nr1,nr2
            if (xbar(ns) .ne. 0.) go to 200
          enddo
          go to 190
c
  200     j2 = ilnobl(uphase(np))
          write (noutpt,1610) uphase(np)(1:j2)
 1610     format(/6x,a,/)
          write (noutpt,1620)
 1620     format(30x,'Mole Fraction')
          do ns = nr1,nr2
            write (noutpt,1630) uspec(ns),xbar(ns)
 1630       format(12x,a24,3x,f6.4)
          enddo
  190     continue
        enddo
        write (noutpt,1340)
      endif
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      write (noutpt,1900)
 1900 format(/' - - - - - - - - - - - - - - - - - - - - - - - - ',
     $ '- - - - - - ',/)
      end
