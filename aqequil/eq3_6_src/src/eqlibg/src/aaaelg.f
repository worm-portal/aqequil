      subroutine aaaelg(ustelg,uveelg)
c
c     EQLIBG: EQ3/6 Library - Activity Coefficients
c     EQ3/6 version 8.0a R43a (Patched 10/01/2009)
c
c     Last revised 04/24/02 by TJW
c
c
c     Copyright (c) 1987, 1990-1993, 1995, 1997, 2002 The Regents of
c     the University of California, Lawrence Livermore National
c     Laboratory. All rights reserved.
c
c     This work was produced at the University of California,
c     Lawrence Livermore National Laboratory (UC LLNL) under
c     contract no. W-7405-ENG-48 between the U.S. Department of
c     Energy (DOE) and The Regents of the University of California
c     (University) for the operation of UC LLNL. Copyright is
c     reserved to the University for purposes of controlled
c     dissemination, commercialization through formal licensing,
c     or other disposition under terms of Contract 48; DOE
c     policies, regulations, and orders; and U.S. statutes.
c
c
c                            DISCLAIMER
c
c     This computer code was prepared as an account of work
c     sponsored by an agency of the United States Government.
c     Neither the United States Government nor the University of
c     California nor any of their employees, makes any warranty,
c     express or implied, or assumes any liability or responsibility
c     for the accuracy, completeness, or usefulness of any
c     information, apparatus, product, or process disclosed, or
c     represents that its use would not infringe privately-owned
c     rights. Reference herein to any specific commercial products,
c     process, or service by trade name, trademark, manufacturer,
c     or otherwise, does not necessarily constitute or imply its
c     endorsement, recommendation, or favoring by the United States
c     Government or the University of California. The views and
c     opinions of authors expressed herein do not necessarily state
c     or reflect those of the United States government or the
c     University of California, and shall not be used for
c     advertising or product endorsement purposes.
c
c
c-----------------------------------------------------------------------
c
c     See the readme.txt file that came with this software for further
c     information, including contact information and references.
c
c
c-----------------------------------------------------------------------
c
c     This library supports the EQ3NR and EQ6 codes. It contains
c     subroutines that deal with calculations specific to activity
c     coefficient models. The EQLIB library contains all other
c     subroutines that address elements of geochemical modeling
c     calculations that are common to EQ3NR and EQ6. The EQLIBU
c     library contains subroutines with generic utility and
c     mathematical functions.
c
c-----------------------------------------------------------------------
c
c     This subroutine is designed to ensure that the copyright statement
c     and legal and other disclaimers appear at the beginning of the
c     concatenated source code for EQLIBG. Concatenation is normally
c     alphabetical, hence the form of the name. This subroutine is
c     generally called by the main program of any code that utilizes
c     subroutines from this library. It returns the stage and version
c     numbers of this library.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       None
c
c     Output:
c
c       ustelg = EQLIBG stage number
c       uveelg = EQLIBG version number
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      character*8 ustelg,uveelg
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      character*8 ust,uve
c
c-----------------------------------------------------------------------
c
      data ust /'R43a    '/
      data uve /'8.0a    '/
c
c-----------------------------------------------------------------------
c
      ustelg = ust
      uveelg = uve
c
      end
