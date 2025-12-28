      subroutine aaaeq6(usteq6,uveeq6)
c
c     EQ6: EQ6 Reaction-Path Code
c     EQ3/6 version 8.0a R43a (Patched 10/01/2009)
c
c     Last revised 04/24/02 by TJW
c
c
c     Copyright (c) 1987, 1990-1993, 1995, 1997, 2002 The Regents
c     of the University of California, Lawrence Livermore National
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
c     This code computes reaction-path models of water/rock and water/
c     rock/waste systems. A model of an aqueous solution must first
c     be created using the EQ3NR code. Problems are composed by choosing
c     substances to react with this water. Such a substance can include
c     another aqueous solution. Thermodynamic data are used in the
c     calculations. Kinetic data of some sort is also required. The
c     code predicts rates of reaction, changes in fluid composition
c     (including pH and Eh), and the precipitation and possible
c     re-dissolution of secondary phases. The temperature and pressure
c     may vary along a reaction path. This code primarily models the
c     chemistry in a small box. It can also model the reaction of a
c     packet of water (the first packet) traversing a reactive medium.
c     This system differs in that secondary minerals are left behind
c     and can't re-dissolve. The fugacities of selected gases may be
c     fixed to simulate reaction in the presence of a large gas volume,
c     such as the atmosphere. This code is similar in function to such
c     codes as CHILLER, REACT, and SOLMINEQ.88.
c
c     The user defines the run parameters on the file called input.
c     The file called DATA1 provides the basic supporting thermodynamic,
c     etc., data. EQ6 writes a pickup file and an output file called
c     output. The pickup file is an EQ6 input file which can be used
c     to restart a calculation. This code also writes data to the
c     tab file. This is an ASCII file containing tables which can be
c     used in graphics post-processing. This code also creates a tabx
c     (scrambled tab) file. Normally this is just used at the end of
c     a run to write the tab file. However, an option exists to use it
c     to create a tab file that spans more than one run.
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
c     concatenated source code for EQ6. Concatenation is normally
c     alphabetical, hence the form of the name. This subroutine returns
c     the stage and version numbers of this computer code.
c
c     This subroutine is called by:
c
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
c       usteq6 = EQ6 stage number
c       uveeq6 = EQ6 version number
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      character*8 usteq6,uveeq6
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
      usteq6 = ust
      uveeq6 = uve
c
      end
