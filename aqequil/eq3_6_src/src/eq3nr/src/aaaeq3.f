      subroutine aaaeq3(usteq3,uveeq3)
c
c     EQ3NR: EQ3NR Speciation-Solubility Code
c     EQ3/6 version 8.0a R43a (Patched 10/01/2009)
c
c     Last revised 04/24/02 by TJW
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
c     This code computes chemical models of aqueous solutions. Problems
c     are composed using mostly analytical data. Thermodynamic data
c     are used in the calculations. This code is similar in function
c     to the codes in the WATEQ series. It can be used by itself to
c     determine aqueous speciation and solution-mineral saturation
c     relations. It can also be used to initialize EQ6 calculations.
c
c     The user defines the run parameters on the file called input.
c     The file called data1 provides the basic supporting thermodynamic,
c     etc., data. EQ3NR writes a pickup file and an output file called
c     output. The pickup file is used to pass data to EQ6. This file
c     is the second half of an EQ6 input file.
c
c     This code requires subroutines from the following EQ3/6 libraries:
c
c       EQLIB
c       EQLIBG
c       EQLIBU
c
c
c-----------------------------------------------------------------------
c
c     This subroutine is designed to ensure that the copyright statement
c     and legal and other disclaimers appear at the beginning of the
c     concatenated source code for EQ3NR. Concatenation is normally
c     alphabetical, hence the form of the name. This subroutine returns
c     the stage and version numbers of this computer code.
c
c     This subroutine is called by:
c
c       EQ3NR/eq3nr.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       None
c
c     Principal output:
c
c       usteq3 = EQ3NR stage number
c       uveeq3 = EQ3NR version number
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      character(len=8) usteq3,uveeq3
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
      character(len=8) ust,uve
c
c-----------------------------------------------------------------------
c
      data ust /'R43a    '/
      data uve /'8.0a    '/
c
c-----------------------------------------------------------------------
c
      usteq3 = ust
      uveeq3 = uve
c
      end
