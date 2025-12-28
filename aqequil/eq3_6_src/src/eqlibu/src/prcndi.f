      subroutine prcndi(noutpt,nttyo)
c
c     This subroutine writes legal statements and disclaimers to the
c     output and screen files.
c
c     This subroutine is called by:
c
c       EQPT/eqpt.f
c       EQ3NR/eq3nr.f
c       EQ6/eq6.f
c
c-----------------------------------------------------------------------
c
c     Input:
c
c       noutpt = the unit number of the output file
c       nttyo  = the unit number of the screen file
c
c     Output:
c
c       None
c
c-----------------------------------------------------------------------
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     Calling sequence variable declarations.
c
      integer noutpt,nttyo
c
c-----------------------------------------------------------------------
c
c     Local variable declarations.
c
c       None
c
c-----------------------------------------------------------------------
c
c     Write notice of applicable statements and disclaimers.
c
      write (noutpt,1000)
      write (nttyo,1000)
 1000 format(' This work is subject to additional statements and',
     $ /' disclaimers which may be found in the README.txt file',
     $ /' included in the EQ3/6 software transmittal package.',//)
c
c     Write additional copyright notice paragraph.
c
c     write (noutpt,1010)
c     write (nttyo,1010)
c
 1010 format(' This work was produced at the University of California,',
     $ /' Lawrence Livermore National Laboratory (UC LLNL) under',
     $ /' contract no. W-7405-ENG-48 between the U.S. Department of',
     $ /' Energy (DOE) and The Regents of the University of California',
     $ /' (University) for the operation of UC LLNL. Copyright is',
     $ /' reserved to the University for purposes of controlled',
     $ /' dissemination, commercialization through formal licensing,',
     $ /' or other disposition under terms of Contract 48; DOE',
     $ /' policies, regulations, and orders; and U.S. statutes.',//)
c
c     Write standard LLNL disclaimer.
c
c     write (noutpt,1020)
c     write (nttyo,1020)
c
 1020 format(24x,'DISCLAIMER',//
     $ ' This computer code was prepared as an account of work',
     $ /' sponsored by an agency of the United States Government.',
     $ /' Neither the United States Government nor the University of',
     $ /' California nor any of their employees, makes any warranty,',
     $ /' express or implied, or assumes any liability or responsi-',
     $ /' bility for the accuracy, completeness, or usefulness of any',
     $ /' information, apparatus, product, or process disclosed, or',
     $ /' represents that its use would not infringe privately-owned',
     $ /' rights. Reference herein to any specific commercial,',
     $ /' product, process, or service by trade name, trademark,',
     $ /' manufacturer, or otherwise, does not necessarily constitute',
     $ /' or imply its endorsement, recommendation, or favoring by the',
     $ /' United States Government or the University of California.',
     $ /' The views and opinions of authors expressed herein do not',
     $ /' necessarily state or reflect those of the United States',
     $ /' government or the University of California, and shall not',
     $ /' be used for advertising or product endorsement purposes.'//)
c
      end
