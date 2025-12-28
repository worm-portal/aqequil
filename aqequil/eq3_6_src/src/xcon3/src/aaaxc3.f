      subroutine aaaxc3(ustxc3,uvexc3)
c
c     XCON3: EQ3NR Input File Converter
c     EQ3/6 version 8.0 R43a (Patched 10/01/2009)
c
c     Last revised 04/24/02 by TJW
c
c
c     Copyright (c) 1993, 1995, 1997, 2002 The Regents of the University
c     of California, Lawrence Livermore National Laboratory. All rights
c     reserved.
c
c     This work was produced at the University of California,
c     Lawrence Livermore National Laboratory (UC LLNL) under
c     contract no. W-7405-ENG-48 between the U.S. Department of
c     Energy (DOE) and The Regents of the University of California
c     (University) for the operation of UC LLNL.  Copyright is
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
c     rights.  Reference herein to any specific commercial products,
c     process, or service by trade name, trademark, manufacturer,
c     or otherwise, does not necessarily constitute or imply its
c     endorsement, recommendation, or favoring by the United States
c     Government or the University of California.  The views and
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
c     This is a utility program to convert EQ3NR input files from
c     one version level to another and from one format to another.
c     This program is most conveniently run using the XCIF3 utility
c     discussed below. It may also be run directly.
c
c     The old input file is known as 'input'. It is preserved in the
c     conversion. The new input file is known as 'newin'. The version
c     levels are:
c
c       '6.0' = versions 6.0-6.1 (3245.0288, 3245.0888)
c       '7.0' = versions 7.0-7.1
c       '7.2' = version 7.2
c       '8.0' = version 8.0
c
c     The formats are:
c
c       'W' = compact
c       'D' = menu-style
c
c     'D' format does not exist for version level '6.0'.
c
c     To specify the choice of version level and format of the new
c     input file, use the control file IXCON. The same control file
c     also governs the EQ6 input file converter XCON6. A copy of the
c     contents of this file can be found in the comments in
c     EQLIBU/rddixc.f.
c
c     This code automatically determines the format of the old input
c     file. It may have trouble determining the version level, however.
c     It will first look for a string of the form 'Version level= X.X'
c     or 'Version number= X.X' in the first title of the first problem
c     on this file. Here 'X.X' denotes the version level. If none is
c     found, the code will use a default value specified on the IXCON
c     file. If no IXCON file is present, the code will use a default
c     specified in the main program.
c
c     If a version level marker as described above exists on the old
c     input file, it will be modified on the new input file to give
c     its correct version level. If no version level marker is on the
c     old input file, one will be added to the new input file.
c
c     The version level and format of the new input file may match
c     those of the old input file. The new input file may be identical
c     to the old one, or made somewhat more pretty.
c
c     As noted above, it is more convenient to run XCON3 using the
c     interface software XCIF3 than to run XCON3 directly. On a UNIX
c     system, XCIF3 is a symbolic link to the shell script xcif36RXX,
c     where RXX is a revision number (R01, R02, R03, etc.). On a
c     PC system, XCIF3 is a Fortran program (XCIF3.EXE; the source
c     code file is XCIF3.FOR) containing DOS system calls. XCIF3
c     replaces the specified input files. XCIF3 also creates the
c     IXCON control file used by XCON3, and deletes it prior to
c     completing execution. To convert a set of EQ3NR input files from
c     version level 7.2 ('W' or 'D' format) to version level 8.0,
c     'W' format, enter:
c
c       xcif3 7.2 8.0 W *.3i
c
c     To convert the same files to 'D' format instead, enter:
c
c       xcif3 7.2 8.0 D *.3i
c
c     XCON6 and XCIF6 are analogues of XCON3 and XCIF3, respectively.
c     They perform the analogous conversions of EQ6 input files.
c
c     This code requires subroutines from the following EQ3/6 libraries:
c
c       EQLIBU
c
c
c-----------------------------------------------------------------------
c
c     This subroutine is designed to ensure that the copyright statement
c     and legal and other disclaimers appear at the beginning of the
c     concatenated source code for XCON3. Concatenation is normally
c     alphabetical, hence the form of the name. This subroutine returns
c     the stage and version numbers of this computer code.
c
c     This subroutine is called by:
c
c       XCON3/xcon3.f
c
c-----------------------------------------------------------------------
c
c     Principal input:
c
c       None
c
c     Principal output:
c
c       ustxc3 = XCON3 stage number
c       uvexc3 = XCON3 version number
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      character*8 ustxc3,uvexc3
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
      data uve /'8.0     '/
c
c-----------------------------------------------------------------------
c
      ustxc3 = ust
      uvexc3 = uve
c
      end
