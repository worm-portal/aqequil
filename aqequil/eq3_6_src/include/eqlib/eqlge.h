c eqlge.h
c
c     Data specific to the generic ion-exchange model.
c
      character*56 ugexr
      character*24 ugexmo,ugexp,ugexs,ugexsr
      character*8 ugexj,uhfgex,uvfgex,uxkgex
c
      common /eqlgec/ ugexj(jet_par,net_par),ugexmo(net_par),
     $ ugexp(net_par),ugexr(iet_par,jet_par,net_par),
     $ ugexs(iet_par,jet_par,net_par),ugexsr(2,iet_par,jet_par,net_par),
     $ uhfgex(iet_par,jet_par,net_par),uvfgex(iet_par,jet_par,net_par),
     $ uxkgex(iet_par,jet_par,net_par)
c
      integer jern1,jern2,jgext,kern1,kern2,kgexsa,ngexro,ngexrt,
     $ ngexsa,ngexso,ngext
c
      common /eqlgei/ jern1(jet_par,net_par),jern2(jet_par,net_par),
     $ jgext(net_par),kern1(net_par),kern2(net_par),
     $ kgexsa(ket_par,net_par),ngexro(iet_par,jet_par,net_par),
     $ ngexrt(jet_par,net_par),ngexsa(iet_par,jet_par,net_par),
     $ ngexso(iet_par,jet_par,net_par),ngext(jet_par,net_par)
c
      real*8 cgexj,mgext,mwtges,tgexp,xhfgex,xlkgex,xvfgex,zgexj
c
      common /eqlge/ cgexj(jet_par,net_par),mgext(jet_par,net_par),
     $ mwtges(net_par),tgexp(net_par),xhfgex(iet_par,jet_par,net_par),
     $ xlkgex(iet_par,jet_par,net_par),xvfgex(iet_par,jet_par,net_par),
     $ zgexj(jet_par,net_par)
c
c     End of INCLUDE file eqlge.h
c-----------------------------------------------------------------------
