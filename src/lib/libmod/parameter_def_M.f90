!====================================================================
      module parameter_def
!    THis module defines some global parameters for the
!    application
!...Translated by Pacific-Sierra Research 77to90 4.3E 11:01:42 1/15/07 
!...Modified by Charlotte Froese Fischer 
!                     Gediminas Gaigalas  10/05/17
!====================================================================
!    KEYORB: Is a packing parameter defined by a packing algorithm.  
!            This is the largest integer for which KEYORB**4 < 2^31 
!            (the largest positive  INTEGER*4)
!    NNNP:   Number  of points in the radial grid 
!    NNN1:   = NNNP+10 
!    NNNW    Maximum number or orbitals (previously 120). An n=10 
!            calculation has 100 orbitals.  This parameter is used
!            to assign array dimensions such as NNNW*NCFG where
!            NCFG is the total number of CSFs.
!    NNNWM1: = NNNW-1
!    NNNWM2: = NNNW-2
!------------------------------------------------------------------  

         integer, parameter :: KEYORB = 215
         integer, parameter :: NNNP   = 590
         integer, parameter :: NNN1   = 600
         integer, parameter :: NNNW   = 127
         integer, parameter :: NNNWM1 = 126
         integer, parameter :: NNNWM2 = 125
      end module parameter_def
!======================================================================
