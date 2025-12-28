      subroutine cpcomb(jassan,jassca,jassne,naapr,nat,ncapr,nccpr,
     $ nnapr,nncpr,nnnpr,nn2pr,naactr,na2ctr,nncatr,nccatr,nc2atr,
     $ nn2ntr,nn3tr)
c
c     Compute all the pair and triplet combinations of aqueous
c     solute species for use with Pitzer's equations. parameters.
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
c       jassca = number of aqueous cation species
c       jassan = number of aqueous anion species (excluding aqueous e-)
c       jassne = number of aqueous neutral species (excluding solvent
c                  water and aqeuous O2(g))
c       nat    = the number of aqueous species
c
c     Principal output:
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
c       nn3tr  = number of nnn triplets
c       nn2ntr = number of nnn' triplets
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer jassan,jassca,jassne,nat
c
      integer naapr,ncapr,nccpr,nnapr,nncpr,nnnpr,nn2pr,
     $ naactr,na2ctr,nccatr,nc2atr,nncatr,nn2ntr,nn3tr
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      integer jn,jc,ja
c
c-----------------------------------------------------------------------
c
      jn = jassne
      jc = jassca
      ja = jassan
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Pair combinations.
c
c     Get the number of ca pairs.
c
      ncapr = jc*ja
c
c     Get the number of cc' pairs.
c
      nccpr = jc*(jc -1)/2
c
c     Get the number of aa' pairs.
c
      naapr = ja*(ja -1)/2
c
c     Get the number of nn pairs.
c
      nn2pr = jn
c
c     Get the number of nn' pairs.
c
      nnnpr = jn*(jn -1)/2
c
c     Get the number of nc pairs.
c
      nncpr = jn*jc
c
c     Get the number of na pairs.
c
      nnapr = jn*ja
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c     Triplet combinations.
c
c     Get the number of cca triplets.
c
      nc2atr = ncapr
c
c     Get the number of aac triplets.
c
      na2ctr = ncapr
c
c     Get the number of nnn triplets.
c
      nn3tr = nn2pr
c
c     Get the number of nnn' triplets.
c
      nn2ntr = 2*nn2pr
c
c     Get the number of cc'a triplets.
c
      nccatr = nccpr*ja
c
c     Get the number of aa'c triplets.
c
      naactr = naapr*jc
c
c     Get the number of nca triplets.
c
      nncatr = jn*jc*ja
c
c* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      end
