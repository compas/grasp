!***********************************************************************
!                                                                      *
      SUBROUTINE SETCON
!                                                                      *
!   This  subprogram  sets the values of the fundamental and derived   *
!   physical constants, and other useful constants.                    *
!                                                                      *
!   Written by Farid A Parpia, at Oxford   Last updated: 06 Oct 1992   *
!                                                                      *
!***********************************************************************
!...Translated by Pacific-Sierra Research 77to90  4.3E  10:50:28   2/14/04
!...Modified by Charlotte Froese Fischer
!                     Gediminas Gaigalas  10/05/17
!-----------------------------------------------
!   M o d u l e s
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE
      USE DEBUG_C
      USE DEF_C
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(DOUBLE), PARAMETER :: AINFCM = 0.52917721067D-08
      REAL(DOUBLE), PARAMETER :: ALFAI  = 137.035999139D00
      REAL(DOUBLE), PARAMETER :: CCMPS  = 2.99792458D10
      REAL(DOUBLE), PARAMETER :: EESU   = 4.803204673D-10
      REAL(DOUBLE), PARAMETER :: EMEG   = 9.10938356D-28
      REAL(DOUBLE), PARAMETER :: EMEAMU = 5.48579909070D-04
      REAL(DOUBLE), PARAMETER :: EMPAMU = 1.007276466879D00
      REAL(DOUBLE), PARAMETER :: HBARES = 1.054571800D-27
      REAL(DOUBLE), PARAMETER :: RINFEV = 13.605693009D00
      REAL(DOUBLE), PARAMETER :: RINFK  = 109737.31568508D00
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
!
!
!   Physical constants: see description below.
!
!
!   The values of the physical constants used here are taken from
!    (see setcon.f-????)
!   E R Cohen and B N Taylor, The 1986 Adjstment of the Fundamental
!   Physical Constants, Report of the CODATA Task Group on Funda-
!   mental COnstants, CODATA Bulletin 63, Pergamon, Elmsford, NY
!   (1986)
!
!      AINFCM: Bohr radius in cm
!      ALFAI : Inverse of the fine-structure constant
!      CCMPS : Speed of light in cm/s
!      EESU  : Electron charge in esu
!              conversion : EESU = EEC * 2997924580
!      EMEG  : Electron mass in g
!      EMEAMU: Electron mass in amu
!      EMPAMU: Proton mass in amu
!      HBARES: Rationalized Planck's constant in erg s
!      RINFEV: Rydberg in eV
!      RINFK : Rydberg in Kaysers
!
!   In the context of new results, Taylor warns that ... since the
!   output values of a least-squares adjustment are related in a
!   complex way and a change in the measured value of one constant
!   usually leads to corresponding changes in the adjusted values
!   of others, one must be cautious in carrying out calculations
!   using both the [above values] and the results of more recent
!   measurements.
!
!   Calculate constants for /DEF3/:
!
      EMPAM = EMPAMU
      RBCM = AINFCM
!
!   Calculate constants for  /DEF11/:
!
!      CVAC is the speed of light in the vacuum in atomic units
!      B1 is the conversion factor between the amu and the atomic
!         unit of mass
!
      CVAC = ALFAI
      AUMAMU = EMEAMU
!
!   Calculate constants for  /DEF10/:
!
!      AUCM converts au to kaysers assuming an infinitely-heavy
!           nucleus
!      AUEV Converst au to eV assuming an infinitely-heavy
!           nucleus
!      CCMS is the speed of light in the vacuum in centimetres
!            per second
!      FASI converts the Einstein A coefficients from atomic to SI
!            units
!      FBSI converts the Einstein B coefficients from atomic to SI
!            units
!
      AUCM = 2.0D00*RINFK
      AUEV = 2.0D00*RINFEV
      CCMS = CCMPS
      FASI = (EMEG/HBARES)*(EESU*EESU/HBARES)**2
      FBSI = (10.0D00*AINFCM**3/HBARES)*FASI
!
!   Calculate conversion factor for Fermis to Bohr radii
!
      FMTOAU = 1.0D-13/AINFCM
!
!   Calculate \pi from FORTRAN function
!
      PI = 4.0D00*ATAN(1.0D00)
!
!   Printouts
!
      IF (LDBPG(2)) WRITE (99, 300) AINFCM, ALFAI, CCMPS, EESU, EMEG, EMEAMU, &
         EMPAMU, HBARES, RINFEV, RINFK
!
      RETURN
!
  300 FORMAT(/,'From SUBROUTINE SETCON:'/,' AINFCM (Bohr radius in cm): ',0P,1D&
         16.9,','/,' ALFAI (Inverse of the fine-structure constant): ',3P,1D&
         16.9,','/,' CCMPS (Speed of light in cm/s): ',1P,1D15.8,','/,&
         ' EESU (Electron charge in esu): ',1P,1D15.8,','/,&
         ' EMEG (Electron mass in g): ',1P,1D14.7,','/,&
         ' EMEAMU (Electron mass in u): ',1P,1D16.9,','/,&
         ' EMPAMU (Proton mass in u): ',1P,1D16.9,','/,&
         ' HBARES (Rationalized Planck constant in erg s): ',1P,1D15.8,','/,&
         ' RINFEV (Rydberg in eV): ',2P,1D15.8,','/,&
         ' RINFK (Rydberg in Kaysers): ',6P,1D17.10,'.')
      RETURN
!
      END SUBROUTINE SETCON
