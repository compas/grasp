!     A GENERAL HARTREE-FOCK PROGRAM
!
!                    C O P Y R I G H T 2008
!
!     by C. Froese Fischer
!        Vanderbilt University
!        Nashville, TN 37235 USA
!
!     August, 2008
!
!     This is a FORTRAN 90/95 version of the HF Fortran 77  program
!     published in Computer Physics Communications, Vol. 43, 355-365 (1987)
!
!     The Pacific Sierra f77tof90 translater was used to produce an initial
!     version in which
!        i) Named COMMON were replacd by a module
!       ii) INTERFACE modules are used to check calling sequences
!      iii) Many obsolete features were removed.
!     
!     Incomplete translations and the removal of EQUIVALENCE  and DATA
!     statements were dealt with by the author.  
!     ------------------------------------------------------------------
!
!     All comments in the program listing assume the radial function P
!     is the solution of an equation of the form
!
!      P" + ( 2Z/R - Y - L(L+1)/R**2 - E)P = X + T
!
!     where Y is called a potential function
!           X is called an exchange function, and
!           T includes contributions from off-diagonal energy parameter
!
!     The program uses LOG(Z*R) as independent variable and
!                      P/SQRT(R) as dependent variable.
!     As a result all equations must be transformed as described in
!     Sec. 6-2 and 6-4 of the book - ``The Hartree-Fock Method for
!     Atoms'',Wiley Interscience, 1977, by Charlotte FROESE FISCHER.
!     (All references to equations and sections refer to this book.)
!
!     Examples on the use of this version may be found in the book
!     "Computational Atomic Structure: An MCHF approach"
!     by C. Froese Fischer, T. Brage, and P. J\"onssn
!     Institute of Physics Publishing, 1997
!
!     Numerical procedures are the same as those for MCHF and are
!     described in Computer Physics Reports, Vol. 3, 273--326 (1986).
!
!     ------------------------------------------------------------------
!     M O D U L E S
!     ------------------------------------------------------------------
!     Defining data types
      module vast_kind_param                                        
         integer, parameter :: byte_log = selected_int_kind(2)      
         integer, parameter :: short_log = selected_int_kind(4)     
         integer, parameter :: long_log = selected_int_kind(18)     
         integer, parameter :: byte = selected_int_kind(2)          
         integer, parameter :: short = selected_int_kind(4)         
         integer, parameter :: long = selected_int_kind(18)         
         integer, parameter :: double = selected_real_kind(13)      
         integer, parameter :: extended = selected_real_kind(30)    
         integer, parameter :: double_ext = selected_real_kind(50)  
         integer, parameter :: dble_complex = selected_real_kind(13)
         integer, parameter :: ext_complex = selected_real_kind(30) 
      end module vast_kind_param                                    
!     ------------------------------------------------------------------
!     C O M M O N   M O D U L E S  
!     ------------------------------------------------------------------
      MODULE blume_C 
        USE vast_kind_param, ONLY:  DOUBLE 
        REAL(DOUBLE), DIMENSION(4) :: COEFN2, COEFNK, COEFVK 
      END MODULE blume_C 
      MODULE coeff_C 
        USE vast_kind_param, ONLY:  DOUBLE 
        INTEGER, DIMENSION(5,5) :: IJPTR 
        REAL(DOUBLE), DIMENSION(200) :: COEF 
      END MODULE coeff_C 
      MODULE DE_C 
        USE vast_kind_param, ONLY:  DOUBLE 
        INTEGER, PARAMETER :: NWFD = 20, NOD = 220 
        INTEGER, DIMENSION(NWFD) :: IND 
        INTEGER :: M, NODE, MK, KK, NJ
        REAL(DOUBLE), DIMENSION(NOD) :: P2, HQ, XX
        REAL(DOUBLE):: V, B4, CN, C, XY, XP, AZZ, PP, FN, EM, FM, EU, FU,DELTAE
      END MODULE de_C 
      MODULE eav_C 
        USE vast_kind_param, ONLY:  DOUBLE 
        REAL(DOUBLE), DIMENSION(10) :: CCA 
        REAL(DOUBLE), DIMENSION(35) :: CCB 
      END MODULE eav_C 
      MODULE estP_C 
        USE vast_kind_param, ONLY:  DOUBLE 
        INTEGER, PARAMETER :: NWFD = 20 
        INTEGER, DIMENSION(NWFD) :: IND
        REAL(DOUBLE), DIMENSION(NWFD) :: ZZ
      END MODULE estP_C
      MODULE fact_C 
        USE vast_kind_param, ONLY:  DOUBLE 
        REAL(DOUBLE), DIMENSION(100) :: GAM 
      END MODULE fact_C 
      MODULE inout_C 
        USE vast_kind_param, ONLY:  DOUBLE 
        INTEGER :: IUF, OUF 
      END MODULE inout_C 
      MODULE label_C 
        USE vast_kind_param, ONLY:  DOUBLE 
        INTEGER, PARAMETER :: NWFD = 20 
        CHARACTER, DIMENSION(NWFD) :: EL*3 
        CHARACTER :: CONFIG*50, ATOM*6, TERM*6 
      END MODULE label_C 
      MODULE param_C 
        USE vast_kind_param, ONLY:  DOUBLE 
        INTEGER :: NO, ND, NWF, NP, NCFG, IB, IC, ID, NSCF, NCLOSD 
        REAL(DOUBLE) :: H, H1, H3, CH, EH, RHO, Z, TOL, D0, D1, D2, D3, D4, D5, &
           D6, D8, D10, D12, D16, D30, FINE 
      END MODULE param_C 
      MODULE radial_C 
        USE vast_kind_param, ONLY:  DOUBLE 
        INTEGER, PARAMETER :: NOD = 220 
        INTEGER, PARAMETER :: NWFD = 20 
        INTEGER, DIMENSION(NWFD) :: L, MAX, N 
        REAL(DOUBLE), DIMENSION(NOD) :: R, RR, R2 
        REAL(DOUBLE), DIMENSION(NOD,NWFD) :: P 
        REAL(DOUBLE), DIMENSION(NOD) :: YK, YR, X
        REAL(DOUBLE), DIMENSION(NWFD) :: AZ 
      END MODULE radial_C 
      MODULE test_C 
        USE vast_kind_param, ONLY:  DOUBLE 
        LOGICAL :: FAIL, OMIT, REL, ALL, TRACE 
      END MODULE test_C 
      MODULE wave_C 
        USE vast_kind_param, ONLY:  DOUBLE 
        INTEGER, PARAMETER :: NOD = 220 
        INTEGER, PARAMETER :: NWFD = 20 
        INTEGER, DIMENSION(NWFD) :: METH, IORD 
        INTEGER :: IPR 
        REAL(DOUBLE), DIMENSION(NOD) :: PDE
        REAL(DOUBLE), DIMENSION(NWFD) :: EK 
        REAL(DOUBLE), DIMENSION(NWFD,NWFD) :: E 
        REAL(DOUBLE), DIMENSION(NWFD) :: SUM, S, DPM, ACC 
        REAL(DOUBLE) :: ED, AZD 
      END MODULE wave_C 
!     ------------------------------------------------------------------
!     I N T E R F A C E  M O D U L E S
!     ------------------------------------------------------------------
      MODULE a_I   
        INTERFACE
        REAL(KIND(0.0D0)) FUNCTION a (I, J, K) 
        INTEGER NOD, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        INTEGER, INTENT(IN) :: I, J, K 
        END FUNCTION  
        END INTERFACE 
      END MODULE 
      MODULE add_I   
        INTERFACE
        SUBROUTINE add (C, K, I, J, FIRST) 
        USE vast_kind_param,ONLY: DOUBLE 
        INTEGER NOD, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        REAL(DOUBLE), INTENT(IN) :: C 
        INTEGER, INTENT(IN) :: K, I, J 
        LOGICAL, INTENT(IN) :: FIRST 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE array_I   
        INTERFACE
        SUBROUTINE array 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE b_I   
        INTERFACE
        REAL(KIND(0.0D0)) FUNCTION b (I, J, K) 
        INTEGER NOD, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        INTEGER, INTENT(IN) :: I, J, K 
        END FUNCTION  
        END INTERFACE 
      END MODULE 
      MODULE bwint_I   
        INTERFACE
        SUBROUTINE bwint (LC, LO) 
        INTEGER, INTENT(IN) :: LC, LO 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE bwzeta_I   
        INTERFACE
        REAL(KIND(0.0D0)) FUNCTION bwzeta (I1) 
        INTEGER NOD, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        INTEGER, INTENT(IN) :: I1 
        END FUNCTION  
        END INTERFACE 
      END MODULE 
      MODULE ca_I   
        INTERFACE
        REAL(KIND(0.0D0)) FUNCTION ca (L, K) 
        INTEGER, INTENT(IN) :: L, K 
        END FUNCTION  
        END INTERFACE 
      END MODULE 
      MODULE cb_I   
        INTERFACE
        REAL(KIND(0.0D0)) FUNCTION cb (L, LP, K) 
        INTEGER, INTENT(IN) :: L, LP, K 
        END FUNCTION  
        END INTERFACE 
      END MODULE 
      MODULE data_I   
        INTERFACE
        SUBROUTINE data 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE de_I   
        INTERFACE
        SUBROUTINE de (I1) 
        INTEGER NOD, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        INTEGER, INTENT(IN) :: I1 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE dev_I   
        INTERFACE
        SUBROUTINE dev (IEL, L, Q, I, DONE) 
        USE vast_kind_param,ONLY: DOUBLE 
        INTEGER :: IEL 
        INTEGER, INTENT(IN) :: L 
        REAL(DOUBLE), INTENT(IN) :: Q 
        INTEGER, INTENT(INOUT) :: I 
        LOGICAL, INTENT(OUT) :: DONE 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE dyk_I   
        INTERFACE
        SUBROUTINE dyk (I, J, K) 
        INTEGER NODi, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        INTEGER, INTENT(IN) :: I, J, K 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE ekin_I   
        INTERFACE
        REAL(KIND(0.0D0)) FUNCTION ekin (I, II, REL) 
        INTEGER NOD, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        INTEGER :: I 
        INTEGER, INTENT(IN) :: II 
        LOGICAL :: REL 
        END FUNCTION  
        END INTERFACE 
      END MODULE 
      MODULE energy_I   
        INTERFACE
        SUBROUTINE energy (ETOTAL) 
        USE vast_kind_param,ONLY: DOUBLE 
        INTEGER NOD, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        REAL(DOUBLE), INTENT(OUT) :: ETOTAL 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE enexpr_I   
        INTERFACE
        SUBROUTINE enexpr (TERM, DONE) 
        INTEGER NOD, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        CHARACTER (LEN = 6), INTENT(IN) :: TERM 
        LOGICAL, INTENT(OUT) :: DONE 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE eptr_I   
        INTERFACE
        SUBROUTINE eptr (EL, ELSYMB, IEL, J1) 
        CHARACTER (LEN = 3), DIMENSION(*), INTENT(IN) :: EL 
        CHARACTER (LEN = 3), INTENT(IN) :: ELSYMB 
        INTEGER, INTENT(OUT) :: IEL, J1 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE factrl_I   
        INTERFACE
        SUBROUTINE factrl (NFACT) 
        INTEGER, INTENT(IN) :: NFACT 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE fk_I   
        INTERFACE
        REAL(KIND(0.0D0)) FUNCTION fk (I, J, K, REL) 
        INTEGER :: I, J, K 
        LOGICAL :: REL 
        END FUNCTION  
        END INTERFACE 
      END MODULE 
      MODULE gk_I   
        INTERFACE
        REAL(KIND(0.0D0)) FUNCTION gk (I, J, K, REL) 
        INTEGER :: I, J, K 
        LOGICAL :: REL 
        END FUNCTION  
        END INTERFACE 
      END MODULE 
      MODULE grange_I   
        INTERFACE
        SUBROUTINE grange 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE help_I   
        INTERFACE
        SUBROUTINE help (CASE) 
        INTEGER, INTENT(IN) :: CASE 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE hl_I   
        INTERFACE
        REAL(KIND(0.0D0)) FUNCTION hl (EL, I, J, REL) 
        INTEGER NOD, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        CHARACTER (LEN = 3), DIMENSION(*), INTENT(IN) :: EL 
        INTEGER, INTENT(IN) :: I, J 
        LOGICAL, INTENT(IN) :: REL 
        END FUNCTION  
        END INTERFACE 
      END MODULE 
      MODULE hnorm_I   
        INTERFACE
        REAL(KIND(0.0D0)) FUNCTION hnorm (N, L, ZZ) 
        USE vast_kind_param,ONLY: DOUBLE 
        INTEGER, INTENT(IN) :: N, L 
        REAL(DOUBLE), INTENT(IN) :: ZZ 
        END FUNCTION  
        END INTERFACE 
      END MODULE 
      MODULE hwf_I   
        INTERFACE
        REAL(KIND(0.0D0)) FUNCTION hwf (N, L, ZZ, R) 
        USE vast_kind_param,ONLY: DOUBLE 
        INTEGER, INTENT(IN) :: N, L 
        REAL(DOUBLE), INTENT(IN) :: ZZ, R 
        END FUNCTION  
        END INTERFACE 
      END MODULE 
      MODULE init_I   
        INTERFACE
        SUBROUTINE init 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE la_I   
        INTERFACE
        INTEGER FUNCTION la (A) 
        character (LEN = 1), INTENT(IN) :: A 
        END FUNCTION  
        END INTERFACE 
      END MODULE 
      MODULE looktm_I   
        INTERFACE
        SUBROUTINE looktm (L, SL, SEN, Q, IP, NSL) 
        USE vast_kind_param,ONLY: DOUBLE 
        INTEGER, INTENT(IN) :: L 
        CHARACTER (LEN = 2), INTENT(IN) :: SL 
        CHARACTER (LEN = 1), INTENT(IN) :: SEN 
        REAL(DOUBLE), INTENT(IN) :: Q 
        INTEGER, INTENT(OUT) :: IP 
        INTEGER, INTENT(OUT) :: NSL 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE lookup_I   
        INTERFACE
        SUBROUTINE lookup (TAB, P1, P2, IND, NO, KEY) 
        INTEGER, DIMENSION(*), INTENT(IN) :: TAB 
        INTEGER, INTENT(IN) :: P1, P2 
        INTEGER, INTENT(OUT) :: IND 
        INTEGER, INTENT(INOUT) :: NO 
        INTEGER, INTENT(IN) :: KEY 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE lval_I   
        INTERFACE
        INTEGER FUNCTION lval (SYMBOL) 
        CHARACTER (LEN = 1), INTENT(IN) :: SYMBOL 
        END FUNCTION  
        END INTERFACE 
      END MODULE 
      MODULE menu_I   
        INTERFACE
        SUBROUTINE menu 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE methd1_I   
        INTERFACE
        SUBROUTINE methd1 (I) 
        INTEGER NOD, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        INTEGER, INTENT(IN) :: I 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE nmrvs_I   
        INTERFACE
        SUBROUTINE nmrvs (NJ, DELTA, MM, PDE, F) 
        USE vast_kind_param,ONLY: DOUBLE 
        INTEGER NOD, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        INTEGER, INTENT(IN) :: NJ 
        REAL(DOUBLE), INTENT(OUT) :: DELTA 
        INTEGER, INTENT(OUT) :: MM 
        REAL(DOUBLE), DIMENSION(NOD), INTENT(INOUT) :: PDE 
        REAL(DOUBLE), DIMENSION(NOD), INTENT(IN) :: F 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE nodec_I   
        INTERFACE
        INTEGER FUNCTION nodec (M) 
        INTEGER NOD, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        INTEGER, INTENT(INOUT) :: M 
        END FUNCTION  
        END INTERFACE 
      END MODULE 
      MODULE orthog_I   
        INTERFACE
        SUBROUTINE orthog 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE output_I   
        INTERFACE
        SUBROUTINE output (PRINT) 
        INTEGER NODi, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        LOGICAL, INTENT(IN) :: PRINT 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE potl_I   
        INTERFACE
        SUBROUTINE potl (I, REL) 
        INTEGER NODi, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        INTEGER, INTENT(IN) :: I 
        LOGICAL :: REL 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE quad_I   
        INTERFACE
        REAL(KIND(0.0D0)) FUNCTION quad (I, M, F, G) 
        USE vast_kind_param,ONLY: DOUBLE 
        INTEGER NOD, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        INTEGER, INTENT(IN) :: I, M 
        REAL(DOUBLE), DIMENSION(NOD), INTENT(IN) :: F, G 
        END FUNCTION  
        END INTERFACE 
      END MODULE 
      MODULE quadr_I   
        INTERFACE
        REAL(KIND(0.0D0)) FUNCTION quadr (I, J, KK) 
        INTEGER NOD, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        INTEGER, INTENT(IN) :: I, J, KK 
        END FUNCTION  
        END INTERFACE 
      END MODULE 
      MODULE quads_I   
        INTERFACE
        REAL(KIND(0.0D0)) FUNCTION quads (I, J, KK) 
        INTEGER NOD, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        INTEGER, INTENT(IN) :: I, J, KK 
        END FUNCTION  
        END INTERFACE 
      END MODULE 
      MODULE reform_I   
        INTERFACE
        SUBROUTINE reform (STR1, STR2) 
        CHARACTER (LEN = 50), INTENT(IN) :: STR1 
        CHARACTER (LEN = 50), INTENT(OUT) :: STR2 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE reord_I   
        INTERFACE
        SUBROUTINE reord (OF, ELC, NWF, IERR) 
        CHARACTER (LEN = 3), DIMENSION(:), INTENT(INOUT) :: OF 
        CHARACTER (LEN = 3), INTENT(IN) :: ELC 
        INTEGER, INTENT(IN) :: NWF 
        INTEGER, INTENT(OUT) :: IERR 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE rk_I   
        INTERFACE
        REAL(KIND(0.0D0)) FUNCTION rk (I, J, II, JJ, K, REL) 
        INTEGER :: I, J, II, JJ, K 
        LOGICAL :: REL 
        END FUNCTION  
        END INTERFACE 
      END MODULE 
      MODULE rlshft_I   
        INTERFACE
        REAL(KIND(0.0D0)) FUNCTION rlshft (I1, I2) 
        INTEGER NODi, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        INTEGER, INTENT(IN) :: I1, I2 
        END FUNCTION  
        END INTERFACE 
      END MODULE 
      MODULE rme_I   
        INTERFACE
        REAL(KIND(0.0D0)) FUNCTION rme (L, LP, K) 
        INTEGER, INTENT(IN) :: L, LP, K
        END FUNCTION  
        END INTERFACE 
      END MODULE 
      MODULE rotate_I   
        INTERFACE
        SUBROUTINE rotate (I, J) 
        INTEGER NOD, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        INTEGER, INTENT(IN) :: I, J 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE scale_I   
        INTERFACE
        SUBROUTINE scale (ZZ) 
        USE vast_kind_param,ONLY: DOUBLE 
        INTEGER NOD 
        PARAMETER (NOD = 220) 
        INTEGER NWFD 
        PARAMETER (NWFD = 20) 
        REAL(DOUBLE), INTENT(IN) :: ZZ 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE scf_I   
        INTERFACE
        SUBROUTINE scf (ETOTAL, SCFTOL, EREL) 
        USE vast_kind_param,ONLY: DOUBLE 
        INTEGER NOD, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        REAL(DOUBLE), INTENT(OUT) :: ETOTAL 
        REAL(DOUBLE), INTENT(IN) :: SCFTOL 
        REAL(DOUBLE), INTENT(OUT) :: EREL 
        END SUBROUTINE  
      END INTERFACE 
      END MODULE 
        MODULE search_I   
        INTERFACE
        SUBROUTINE search (NJ, I) 
        INTEGER NODi, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        INTEGER, INTENT(OUT) :: NJ 
        INTEGER, INTENT(IN) :: I 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE sn_I   
        INTERFACE
        REAL(KIND(0.0D0)) FUNCTION sn (I, J, II, JJ, K) 
        INTEGER :: I, J, II, JJ, K 
        END FUNCTION  
        END INTERFACE 
      END MODULE 
      MODULE solve_I   
        INTERFACE
        SUBROUTINE solve (I, FIRST, REL) 
        INTEGER NOD, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        INTEGER, INTENT(IN) :: I 
        LOGICAL, INTENT(IN) :: FIRST 
        LOGICAL :: REL 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE summry_I   
        INTERFACE
        SUBROUTINE summry (ET, EREL) 
        USE vast_kind_param,ONLY: DOUBLE 
        INTEGER NOD, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        REAL(DOUBLE), INTENT(IN) :: ET 
        REAL(DOUBLE), INTENT(IN) :: EREL 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE vk_I   
        INTERFACE
        REAL(KIND(0.0D0)) FUNCTION vk (I, J, II, JJ, K) 
        INTEGER :: I, J, II, JJ, K 
        END FUNCTION  
        END INTERFACE 
      END MODULE 
      MODULE wavefn_I   
        INTERFACE
        SUBROUTINE wavefn 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE xch_I   
        INTERFACE
        SUBROUTINE xch (I, IOPT) 
        INTEGER NOD, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        INTEGER, INTENT(IN) :: I 
        INTEGER, INTENT(IN) :: IOPT 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE ykf_I   
        INTERFACE
        SUBROUTINE ykf (I, J, K, REL) 
        INTEGER NOD, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        INTEGER, INTENT(IN) :: I, J, K
        LOGICAL, INTENT(IN) :: REL
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
      MODULE zk_I   
        INTERFACE
        SUBROUTINE zk (I, J, K) 
        INTEGER NOD, NWFD 
        PARAMETER (NOD = 220) 
        PARAMETER (NWFD = 20) 
        INTEGER, INTENT(IN) :: I, J, K 
        END SUBROUTINE  
        END INTERFACE 
      END MODULE 
!     ------------------------------------------------------------------
!       M A I N    P R O G R A M
!     ------------------------------------------------------------------
!
!       The MAIN program controls the overall calculation  that
!   may  consist  of  a  series of atoms or ions in an iso-electronic
!   sequence.  Initial  estimates  for  the first are obtained either
!   from a file WFN.INP, if it exists, and scaled for the appropriate
!   Z, if necessary, or from a screened hydrogenic approximation.  All
!   others are  obtained  by  scaling  the previous results using the
!   scaling of Sec.  (7-2).
!
      PROGRAM MAIN 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE TEST_C 
      USE LABEL_C 
      USE WAVE_C, ONLY: E, DPM, ACC, NOD 
      USE PARAM_C 
      USE INOUT_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE init_I 
      USE data_I 
      USE help_I 
      USE scf_I 
      USE output_I 
      USE summry_I 
      USE menu_I 
      USE scale_I 
      USE orthog_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J 
      REAL(DOUBLE) :: SCFTOL, ETOTAL, EREL, ZZ 
      LOGICAL :: PRINT, STRONG, OLD 
      CHARACTER(LEN=1) :: ANS, ASTER = '*' 
!-----------------------------------------------
!
!  *****  WRITE OUT HEADER
!
      WRITE (6, 9) 
    9 FORMAT(/,/,/,/,/,/,/,22X,'=============================',/,22X,&
         ' H A R T R E E - F O C K . 86',/,22X,'=============================') 
!
!  *****  WRITE OUT DIMENSION INFORMATION
!
      WRITE (6, 99) 'NWF', NWFD, 'NO', NOD 
   99 FORMAT(/,/,15X,'THE DIMENSIONS FOR THE CURRENT VERSION ARE:'/,13X,2(10X,2&
         (A6,'=',I3,4X),/),/) 
!
!  *****  INITIALIZE
!
      CALL INIT 
!
!  ***** SET UNIT NUMBERS AND OPEN FILES
!
      WRITE (6, '(2/A/A,2/)') ' START OF CASE', ' =============' 
      INQUIRE(FILE='wfn.inp', EXIST=OLD) 
      IF (OLD) THEN 
         IUF = 21 
         OPEN(UNIT=IUF, FILE='wfn.inp', STATUS='OLD', FORM='UNFORMATTED', &
            POSITION='asis') 
      ELSE 
         IUF = 0 
      ENDIF 
      OUF = 31 
      OPEN(UNIT=OUF, FILE='wfn.out', STATUS='UNKNOWN', FORM='UNFORMATTED', &
         POSITION='asis') 
      OPEN(UNIT=3, FILE='hf.log', STATUS='UNKNOWN', POSITION='asis') 
!
      FAIL = .FALSE. 
      DPM(:NWFD) = D10 
      E(:NWFD,:NWFD) = D0 
!
!  *****  DETERMINE DATA ABOUT THE PROBLEM
!
      CALL DATA 
!
!  *****  SET PARAMETERS TO THEIR DEFAULT VALUE
!
   13 CONTINUE 
      PRINT = .FALSE.  ! plot.dat is generated when PRINT = .TRUE.
      PRINT = .TRUE.  ! plot.dat is generated when PRINT = .TRUE.
      SCFTOL = 1.D-8 
      NSCF = 12 
      IC = 2 + (NWF + 1 - IB)/4 
      TRACE = .FALSE.       !! Set to .TRUE. if tracing is desired.
      IF (IB <= NWF) THEN 
         WRITE (0, '(/A)') ' Default values for remaining parameters? (Y/N/H) ' 
         READ(5,'(A)') ANS
         IF (ANS=='H' .OR. ANS=='h') THEN 
            CALL HELP (4) 
            GO TO 13 
         ENDIF 
         IF (ANS/='Y' .AND. ANS/='y') THEN 
!
!  *****  ADDITIONAL PARAMETERS
!
   50       CONTINUE 
            WRITE (0, '(/A)') ' Default values (NO,STRONG) ? (Y/N/H) ' 
            READ (5, '(A)') ANS 
            IF (ANS=='H' .OR. ANS=='h') THEN 
               CALL HELP (3) 
               GO TO 50 
            ENDIF 
            IF (ANS/='Y' .AND. ANS/='y') THEN 
               WRITE (0, *) ' Enter values in FORMAT(I3,1X,L1) ' 
               READ (5, '(I3,1X,L1)') NO, STRONG 
               IF (NO > NOD) THEN 
                  WRITE (0, '(A,A,I4)') ' TOO MANY POINTS: the allowed', &
                     ' MAXIMUM is ', NOD 
                  GO TO 50 
               ENDIF 
               ND = NO - 2 
               OMIT = .NOT.STRONG 
            ENDIF 
   16       CONTINUE 
            WRITE (0, '(A)') ' Default values for PRINT, SCFTOL ? (Y/N/H)' 
            READ (5, '(A)') ANS 
            IF (ANS=='H' .OR. ANS=='h') THEN 
               CALL HELP (5) 
               GO TO 16 
            ENDIF 
            IF (ANS/='Y' .AND. ANS/='y') THEN 
!CFF      W>= D+7  (was E6.1)
               WRITE (0, '(A)') ' Input FORMAT(L1, 1X, E8.1) ' 
               READ (5, '(L1,1X,E8.1)') PRINT, SCFTOL 
            ENDIF 
   17       CONTINUE 
            WRITE (0, '(A)') ' Default values for NSCF, IC ? (Y/N/H) ' 
            READ (5, '(A)') ANS 
            IF (ANS=='H' .OR. ANS=='h') THEN 
               CALL HELP (6) 
               GO TO 17 
            ENDIF 
            IF (ANS/='Y' .AND. ANS/='y') THEN 
               WRITE (0, '(A)') ' Input FORMAT(I2, 1X, I1) ' 
               READ (5, '(I2,1X,I1)') NSCF, IC 
            ENDIF 
   18       CONTINUE 
            WRITE (0, '(A)') ' Default values for TRACE ? (Y/N/H) ' 
            READ (5, '(A)') ANS 
            IF (ANS=='H' .OR. ANS=='h') THEN 
               CALL HELP (7) 
               GO TO 18 
            ENDIF 
            IF (ANS=='N' .OR. ANS=='n') TRACE = .TRUE. 
         ENDIF 
      ENDIF 
!
!
!  *****  PERFORM THE MCHF ITERATION
!
      CALL SCF (ETOTAL, SCFTOL, EREL) 
!
!  *****  OUTPUT RESULTS IF PRINT = .TRUE.
!
      CALL OUTPUT (PRINT) 
      IF (.NOT.FAIL) THEN 
         CALL SUMMRY (ETOTAL, EREL) 
   19    CONTINUE 
         WRITE (0, '(/A)') ' Additional parameters ? (Y/N/H) ' 
         READ(5,'(A1)') ANS
         IF (ANS=='H' .OR. ANS=='h') THEN 
            CALL HELP (8) 
            GO TO 19 
         ENDIF 
         IF (ANS=='Y' .OR. ANS=='y') CALL MENU 
!
!  *****  CHECK FOR ISOELECTRONIC SEQUENCE OR END OF CASE.
!
   20    CONTINUE 
         WRITE (0, '(/A)') ' Do you wish to continue along the sequence ? ' 
         READ(5,'(A1)') ANS
         IF (ANS=='H' .OR. ANS=='h') THEN 
            CALL HELP (9) 
            GO TO 20 
         ENDIF 
         IF (ANS=='Y' .OR. ANS=='y') THEN 
            WRITE (0, *) '   Enter: ATOM, ZZ, (ACC(I),I=1,NWF) in ', &
               ' format(A6,F6.0,(20F3.1))' 
            READ (5, '(A6,F6.0,(20F3.1))') ATOM, ZZ, (ACC(I),I=1,NWF) 
!
!  *****  SCALE RESULTS FOR ANOTHER MEMBER OF THE ISOELECTRONIC SEQUENCE
!
            CALL SCALE (ZZ) 
            WRITE (3, 14) ATOM, TERM 
   14       FORMAT('1',9X,2A6) 
            CALL ORTHOG 
            GO TO 13 
         ENDIF 
      ENDIF 
!
!  *****  DETERMINE END OF CASE
!
      WRITE (6, '(2/A/A,2/)') ' END OF CASE', ' ===========' 
      STOP  
      END PROGRAM MAIN 


!
!     ----------------------------------------------------------------
!               A
!     ----------------------------------------------------------------
!
!       Determine the coefficient in the potential for electron i of
!       Y^k(j,j)
 
      REAL(KIND(0.0D0)) FUNCTION A (I, J, K) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE COEFF_C 
      USE PARAM_C 
      USE RADIAL_C 
      USE WAVE_C, ONLY: SUM 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ca_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I 
      INTEGER , INTENT(IN) :: J 
      INTEGER  :: K 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ISTART 
      REAL(DOUBLE) :: C 
!-----------------------------------------------
!
      IF (I>NCLOSD .AND. J>NCLOSD) THEN 
         ISTART = IJPTR(I-NCLOSD,J-NCLOSD) + 1 
         A = COEF(ISTART+K/2) 
      ELSE IF (I == J) THEN 
         C = SUM(I) - D1 
         IF (K == 0) THEN 
            A = C 
         ELSE 
            A = -C*CA(L(I),K) 
         ENDIF 
      ELSE IF (K == 0) THEN 
         A = SUM(J) 
      ELSE 
         A = D0 
      ENDIF 
      RETURN  
      END FUNCTION A 
!
!     ----------------------------------------------------------------
!               A D D
!     ----------------------------------------------------------------
!
!     Add a Slater integral to the data structure associated with the
!     energy expression
!
      SUBROUTINE ADD(C, K, I, J, FIRST) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE COEFF_C 
      USE PARAM_C 
      USE RADIAL_C 
      USE WAVE_C, ONLY: SUM 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: K 
      INTEGER , INTENT(IN) :: I 
      INTEGER , INTENT(IN) :: J 
      REAL(DOUBLE) , INTENT(IN) :: C 
      LOGICAL , INTENT(IN) :: FIRST 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IP 
!-----------------------------------------------
!
      IP = IJPTR(I-NCLOSD,J-NCLOSD) 
 
      IF (FIRST) THEN 
         COEF(IP+K/2+1) = C/SUM(I) + COEF(IP+K/2+1) 
      ELSE 
         IP = IP + MIN(L(I),L(J)) + 1 + (K - ABS(L(I)-L(J)))/2 + 1 
         COEF(IP) = COEF(IP) + C/SUM(I) 
      ENDIF 
      RETURN  
      END SUBROUTINE ADD 
!
!     ----------------------------------------------------------------
!               A R R A Y
!     ----------------------------------------------------------------
!
!     Set up the data structure associated with the average energy
!
      SUBROUTINE ARRAY 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE COEFF_C 
      USE PARAM_C 
      USE RADIAL_C 
      USE WAVE_C, ONLY: SUM 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ca_I 
      USE cb_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IP, I, ISUMI, J, ISUMJ, K 
      REAL(DOUBLE) :: DSUMI, DSUMJ, C 
!-----------------------------------------------
!
      IP = 0 
      DO I = NCLOSD + 1, NWF 
         ISUMI = SUM(I) 
         DSUMI = SUM(I) - ISUMI 
         DO J = NCLOSD + 1, NWF 
            ISUMJ = SUM(J) 
            DSUMJ = SUM(J) - ISUMJ 
            IF (I /= J) THEN 
               C = SUM(J) 
               IF (DSUMI/=D0 .AND. DSUMJ/=D0) C = (DSUMI*(ISUMI + 1)*ISUMJ + &
                  DSUMJ*(ISUMJ + 1)*ISUMI)/SUM(I) 
            ELSE 
               C = SUM(I) - D1 
               IF (DSUMI /= D0) C = (ISUMI*(SUM(I)+DSUMI-1))/SUM(I) 
            ENDIF 
!
            IJPTR(I-NCLOSD,J-NCLOSD) = IP 
!
!  *****        Direct contribution
!
            DO K = 0, 2*MIN0(L(I),L(J)), 2 
               IP = IP + 1 
               IF (IP > 200) STOP ' COEF array too small: MAX = (200)' 
               COEF(IP) = D0 
               IF (K == 0) THEN 
                  COEF(IP) = C 
               ELSE IF (I == J) THEN 
                  COEF(IP) = -C*CA(L(I),K) 
               ENDIF 
            END DO 
!
!  *****        Exchange contribution
!
            IF (I == J) CYCLE  
            DO K = ABS(L(I)-L(J)), L(I) + L(J), 2 
               IP = IP + 1 
               IF (IP > 200) STOP ' COEF array too small: MAX = (200)' 
               COEF(IP) = -C*CB(L(I),L(J),K) 
            END DO 
         END DO 
      END DO 
      RETURN  
      END SUBROUTINE ARRAY 
!
!     ----------------------------------------------------------------
!               B
!     ----------------------------------------------------------------
!
!     Determine the coefficient of the Y^k(i,j)P(j) term in the exchange
!     expression of electron i
!
      REAL(KIND(0.0D0)) FUNCTION B (I, J, K) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE COEFF_C 
      USE PARAM_C 
      USE RADIAL_C 
      USE WAVE_C, ONLY: SUM 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE cb_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I 
      INTEGER , INTENT(IN) :: J 
      INTEGER  :: K 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LL, ISTART, KK 
!-----------------------------------------------
!
      IF (I == J) THEN 
         B = D0 
      ELSE IF (I>NCLOSD .AND. J>NCLOSD) THEN 
!
!   ..... LL is the number of direct terms
!         ISTART the beginning of the exchange terms
!
         LL = MIN(L(I),L(J)) + 1 
         ISTART = IJPTR(I-NCLOSD,J-NCLOSD) + 1 + LL 
         KK = (K - ABS(L(I)-L(J)))/2 
         B = COEF(ISTART+KK) 
      ELSE 
         B = -SUM(J)*CB(L(I),L(J),K) 
      ENDIF 
      RETURN  
      END FUNCTION B 

!
!     ------------------------------------------------------------------
!               B W I N T
!     ------------------------------------------------------------------
!
      SUBROUTINE BWINT(LC, LO) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE BLUME_C 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: LC, LO 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LC1 
!-----------------------------------------------
!
! ... LC IS THE L-VALUE OF THE FILLED SUBSHELL, LO IS THE L-VALUE
!     OF THE PARTIALLY-FILLED SUBSHELL.
!
      IF (LC>3 .OR. LO>4) THEN 
         WRITE (0, 100) LC, LO 
  100    FORMAT(' INCORRECT CALLING OF BWINT WITH LC =',I2,', LO =',I2) 
      ENDIF 
      IF (LC == 0) then
!
! ... S-P
!
         IF (LO == 1) THEN
            COEFNK(1) = 1.D0 
            COEFN2(1) = -2.D0 
            COEFVK(1) = 1.D0 
            RETURN  
!
! ... S-D
!
         ELSE IF (LO == 2) THEN
            COEFNK(1) = 6.D0/5.D0 
            COEFN2(1) = -9.D0/5.D0 
            COEFVK(1) = 3.D0/5.D0 
            RETURN  
!
! ... S-F
!
         ELSE IF (LO == 3) THEN
            COEFNK(1) = 9.D0/7.D0 
            COEFN2(1) = -12.D0/7.D0 
            COEFVK(1) = 3.D0/7.D0 
            RETURN  
!
! ... S-G
!
         ELSE IF (LO == 4) THEN
            COEFNK(1) = 4.D0/3.D0 
            COEFN2(1) = -5.D0/3.D0 
            COEFVK(1) = 1.D0/3.D0 
            RETURN  
         END IF
      ELSE IF (LC == 1) then
!
! ... P-P
!
         IF (LO ==1 ) THEN
            COEFNK(1) = 0.D0 
            COEFN2(1) = 3.D0 
            COEFVK(1) = 9.D0/5.D0 
            RETURN  
!
! ... P-D
!        
         ELSE IF (LO == 2) THEN
            COEFNK(1) = 3.D0/7.D0 
            COEFNK(2) = 36.D0/35.D0 
            COEFN2(1) = -12.D0/5.D0 
            COEFN2(2) = 0.D0 
            COEFVK(1) = 3.D0/5.D0 
            COEFVK(2) = 36.D0/35.D0 
            RETURN  
!
! ... P-F
!
         ELSE IF (LO == 3) THEN
            COEFNK(1) = 1.D0/7.D0 
            COEFNK(2) = 10.D0/7.D0 
            COEFN2(1) = -18.D0/7.D0 
            COEFN2(2) = 0.D0 
            COEFVK(1) = 18.D0/35.D0 
            COEFVK(2) = 5.D0/7.D0 
            RETURN  
!
! ... P-G
!
!         
         ELSE IF (LO == 4) THEN
            COEFNK(1) = 5.D0/77.D0 
            COEFNK(2) = 18.D0/11.D0 
            COEFN2(1) = -18.D0/7.D0 
            COEFN2(2) = 0.D0 
            COEFVK(1) = 3.D0/7.D0 
            COEFVK(2) = 6.D0/11.D0 
            RETURN  
         END IF
      ELSE  IF (LC == 2) then
!
! ... D-P
!
         IF (LO == 1) then
            COEFNK(1) = 59.D0/7.D0 
            COEFNK(2) = -18.D0/7.D0 
            COEFN2(1) = -4.D0 
            COEFN2(2) = 0.D0 
            COEFVK(1) = -1.D0 
            COEFVK(2) = 18.D0/7.D0 
            RETURN  
!
! ... D-D
!
         ELSE IF (LO == 2) THEN
            COEFNK(1) = 6.D0/7.D0 
            COEFNK(2) = 0.D0 
            COEFN2(1) = 3.D0 
            COEFN2(2) = 0.D0 
            COEFVK(1) = 3.D0/7.D0 
            COEFVK(2) = 10.D0/7.D0 
            RETURN  
!
! ... D-F
!
         ELSE IF (LO == 3) THEN
            COEFNK(1) = 9.D0/7.D0 
            COEFNK(2) = -13.D0/77.D0 
            COEFNK(3) = 75.D0/77.D0 
            COEFN2(1) = -18.D0/7.D0 
            COEFN2(2) = 0.D0 
            COEFN2(3) = 0.D0 
            COEFVK(1) = 3.D0/7.D0 
            COEFVK(2) = 3.D0/7.D0 
            COEFVK(3) = 75.D0/77.D0 
            RETURN  
!
! ... D-G
!
         ELSE IF (LO == 4) THEN
            COEFNK(1) = 741.D0/693.D0 
            COEFNK(2) = -215.D0/429.D0 
            COEFNK(3) = 210.D0/143.D0 
            COEFN2(1) = -3.D0 
            COEFN2(2) = 0.D0 
            COEFN2(3) = 0.D0 
            COEFVK(1) = 3.D0/7.D0 
            COEFVK(2) = 255.D0/693.D0 
            COEFVK(3) = 105.D0/143.D0 
            RETURN  
         END IF
      ELSE IF (LC == 3) THEN
!
! ... F-P
         IF (LO == 1) THEN
            COEFNK(1) = 52.D0/3.D0 
            COEFNK(2) = -20.D0/3.D0 
            COEFN2(1) = -9.D0 
            COEFN2(2) = 0.D0 
            COEFVK(1) = -9.D0/5.D0 
            COEFVK(2) = 10.D0/3.D0 
            RETURN  
   !
   ! ... F-D
   !
         ELSE IF (LO == 2) THEN
            COEFNK(1) = 5.D0 
            COEFNK(2) = 142.D0/55.D0 
            COEFNK(3) = -20.D0/11.D0 
            COEFN2(1) = -18.D0/5.D0 
            COEFN2(2) = 0.D0 
            COEFN2(3) = 0.D0 
            COEFVK(1) = -3.D0/5.D0 
            COEFVK(2) = 2.D0/5.D0 
            COEFVK(3) = 20.D0/11.D0 
            RETURN  
!
! ... F-F
!
         ELSE IF (LO == 3) THEN
            COEFNK(1) = 1.D0 
            COEFNK(2) = 5.D0/11.D0 
            COEFNK(3) = 0.D0 
            COEFN2(1) = 3.D0 
            COEFN2(2) = 0.D0 
            COEFN2(3) = 0.D0 
            COEFVK(1) = 1.D0/5.D0 
            COEFVK(2) = 5.D0/11.D0 
            COEFVK(3) = 175.D0/143.D0 
            RETURN  
!
! ... F-G
!
         ELSE IF (LO == 4) THEN
            COEFNK(1) = 53.D0/33.D0 
            COEFNK(2) = 57.D0/143.D0 
            COEFNK(3) = -115.D0/429.D0 
            COEFNK(4) = 392.D0/429.D0 
            COEFN2(1) = -8.D0/3.D0 
            COEFN2(2) = 0.D0 
            COEFN2(3) = 0.D0 
            COEFN2(4) = 0.D0 
            COEFVK(1) = 1.D0/3.D0 
            COEFVK(2) = 3.D0/11.D0 
            COEFVK(3) = 57.D0/143.D0 
            COEFVK(4) = 392.D0/429.D0 
            RETURN  
         END IF
      ELSE
         WRITE (0, 100) LC, LO 
      END IF
      END SUBROUTINE BWINT 
!
!     ----------------------------------------------------------------
!               B W Z E T A
!     ----------------------------------------------------------------
!  ***** COMPUTES THE NUCLEAR SPIN-ORBIT PARAMETER AND THE
!        CORRECTIONS FOR THE OTHER ELECTRONS
!        USING THE FORMULA DERIVED BY Blume and Watson.
!
      REAL(KIND(0.0D0)) FUNCTION BWZETA (I1) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE BLUME_C 
      USE PARAM_C 
      USE RADIAL_C 
      USE WAVE_C, ONLY: SUM 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE quadr_I 
      USE sn_I 
      USE bwint_I 
      USE vk_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I1 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LB, I, LA, KE1, IP, K 
      REAL(DOUBLE), DIMENSION(3) :: SS 
      REAL(DOUBLE) :: ZETA, C 
!-----------------------------------------------
!
      ZETA = FINE*Z*QUADR(I1,I1,-3) 
 
      LB = L(I1) 
      DO I = 1, NWF 
         IF (I == I1) CYCLE  
         LA = L(I) 
         ZETA = ZETA - SUM(I)*SN(I1,I,I1,I,0) 
         IF (SUM(I) /= 4*L(I) + 2) CYCLE  
         CALL BWINT (LA, LB) 
         KE1 = 2 
         IF (LA /= LB) KE1 = IABS(LA - LB) 
         IP = 0 
         DO K = KE1, LA + LB, 2 
            IP = IP + 1 
            ZETA = ZETA + COEFN2(IP)*SN(I1,I,I,I1,K - 2) + COEFNK(IP)*SN(I,I1,&
               I1,I,K) + COEFVK(IP)*(VK(I1,I,I,I1,K - 1) - VK(I,I1,I1,I,K - 1)) 
         END DO 
         WRITE (*, *) 'zeta,i', ZETA, I 
 
      END DO 
      ZETA = D2*ZETA 
      WRITE (*, *) 'zeta', ZETA 
      C = SUM(I1) 
      IF (C /= D1) THEN 
         SS(1) = SN(I1,I1,I1,I1,0) 
         C = C + C - D3 
         ZETA = ZETA - C*SS(1) 
         IF (LB == 2) THEN 
            SS(2) = SN(I1,I1,I1,I1,2) 
            ZETA = ZETA + SS(2)*6.D0/7.D0 
         ELSE IF (LB == 3) THEN 
            SS(2) = SN(I1,I1,I1,I1,2) 
            SS(3) = SN(I1,I1,I1,I1,4) 
            ZETA = ZETA + SS(2) + SS(3)/2.2D0 
         ENDIF 
      ENDIF 
      WRITE (*, *) 'zeta', ZETA 
 
      BWZETA = ZETA 
      RETURN  
      END FUNCTION BWZETA 
!
!     ------------------------------------------------------------------
!               C A
!     ------------------------------------------------------------------
!
      REAL(KIND(0.0D0)) FUNCTION CA (L, K) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE EAV_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rme_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: L 
      INTEGER  :: K 
!-----------------------------------------------
!
      IF (L <= 4) THEN 
         CA = CCA((L*(L-1)+K)/2) 
      ELSE 
!        Corrected according to Prof. P. Bogdanovich 1996.03.18.
!        CA = RME(L,L,K)**2
         CA = RME(L,L,K)**2/((2*L + 1)*(4*L + 1)) 
      ENDIF 
      RETURN  
      END FUNCTION CA 
!
!     -----------------------------------------------------------------
!                 C B
!     -----------------------------------------------------------------
      REAL(KIND(0.0D0)) FUNCTION CB (L, LP, K) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE EAV_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rme_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: L, LP, K 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(0:4) :: ICBPTR = (/1, 6, 14, 23, 31/)
      INTEGER :: L1, L2 
!-----------------------------------------------
!
      IF (L <= LP) THEN 
         L1 = L 
         L2 = LP 
      ELSE 
         L1 = LP 
         L2 = L 
      ENDIF 
      IF (L2 <= 4) THEN 
         CB = CCB(ICBPTR(L1)+(K+L1-L2)/2+(L1+1)*(L2-L1)) 
      ELSE 
         CB = RME(L,LP,K)**2/(2*(2*L + 1)*(2*LP + 1)) 
      ENDIF 
      RETURN  
      END FUNCTION CB 
!
!     ------------------------------------------------------------------
!               D A T A
!     ------------------------------------------------------------------
!       Data concerning the number of configurations (NCFG), the number
!   and type of electrons in each  configuration,  as  well  as  data
!   associated with the energy expression are read and stored.
!
      SUBROUTINE DATA 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE TEST_C 
      USE LABEL_C 
      USE PARAM_C 
      USE WAVE_C, ONLY: EK, E, SUM, S, ACC, METH, IORD, NOD 
      USE RADIAL_C, ONLY: AZ, L, N 
      USE ESTP_C, ONLY: IND 
      USE INOUT_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE reform_I 
      USE lval_I 
      USE help_I 
      USE reord_I 
      USE array_I 
      USE enexpr_I 
      USE eptr_I 
      USE add_I 
      USE wavefn_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J1, J2, I, J, IFULL, MAXORB, K, NIT, JJ, NEXT, IERR, KKK, &
         ISELEC, KFG, IFG, JFG, JP, N1, L1, N2, L2, ITEMP 
      REAL(DOUBLE) :: SS, CFG 
      LOGICAL :: FIRST, STRONG, DONE, ORDERD 
      CHARACTER :: ANS, STRING*50, EL1*3, EL2*3 
      CHARACTER , DIMENSION(18) :: ELCSD*3 
      CHARACTER(LEN=1) :: ASTER = '*', W 
!-----------------------------------------------
!
    1 FORMAT(18(1X,A3)) 
    7 FORMAT(A3,F6.0,I3,I3,F3.1) 
    5 CONTINUE 
      WRITE (0, '(A/A)') ' Enter ATOM,TERM,Z', &
         ' Examples: O,3P,8. or Oxygen,AV,8.' 
      READ (5, '(A50)') STRING 
      I = INDEX(STRING,',') 
      IF (I == 0) THEN 
         WRITE (0, *) ' ATOM, TERM, and Z must be separated by commas ' 
         GO TO 5 
      ENDIF 
      ATOM = STRING(1:I-1) 
      J = INDEX(STRING(I+1:),',') 
      IF (J == 0) THEN 
         WRITE (0, *) ' ATOM, TERM, and Z must be separated by commas ' 
         GO TO 5 
      ENDIF 
      TERM = STRING(I+1:I+J-1) 
      READ (STRING(I+J+1:), '(F3.0)') Z 
!
!  *****  INPUT COMMON CLOSED SHELLS
!
    2 CONTINUE 
      WRITE (0, *) 
      WRITE (0, '(A,A)') ' List the CLOSED shells in the fields indicated', &
         ' (blank line if none)' 
      WRITE (0, '(A)') ' ... ... ... ... ... ... ... ... etc.' 
      READ (5, 1) (ELCSD(I),I=1,18) 
!
!  *****  INPUT THE CONFIGURATION
!
      WRITE (0, '(/A,A/A)') ' Enter electrons outside CLOSED shells ', &
         '(blank line if none)', ' Example: 2s(1)2p(3)' 
      READ (5, '(A)') STRING 
      CALL REFORM (STRING, CONFIG) 
!
!      Determine the number of closed shells
!
      I = 0 
      SS = D0 
   12 CONTINUE 
      IF (ELCSD(I+1) /= '   ') THEN 
         I = I + 1 
         EL(I) = ADJUSTR(ELCSD(I)) 
         J = 3 
         IF (EL(I)(1:1) /= ' ') J = 2 
         L(I) = LVAL(EL(I)(J:J)) 
         N(I) = ICHAR(EL(I)(J-1:J-1)) - ICHAR('1') + 1 
         IFULL = 2*(2*L(I)+1) 
         SUM(I) = IFULL 
         S(I) = SS + IFULL/2 
         SS = SS + IFULL 
         METH(I) = 1 
         ACC(I) = D0 
         IND(I) = 0 
         IF (IUF /= 0) IND(I) = -1 
         IF (I < 18) GO TO 12 
         STOP ' TOO MANY CLOSED SHELLS: MAX = 18' 
      ENDIF 
      NCLOSD = I 
!
!  *****  DETERMINE THE OTHER ELECTRONS
!
      MAXORB = NCLOSD 
      STRING = CONFIG 
      J = 2 
      I = 0 
   16 CONTINUE 
      IF (STRING(J:J+2) /= '   ') THEN 
!
!  --------- An electron has been found; is it a new one?
!
         I = I + 1 
         IF (I > 5) STOP ' TOO MANY SHELLS: MAX= (5)' 
         EL1 = STRING(J:J+2) 
         K = NCLOSD + 1 
   17    CONTINUE 
         IF (K <= MAXORB) THEN 
            IF (EL(K) /= EL1) THEN 
               K = K + 1 
               IF (K > NWFD) THEN 
                  WRITE (0, '(A,I4)') ' TOO MANY ELECTRONS: MAX =', NWFD 
                  GO TO 2 
               ELSE 
                  GO TO 17 
               ENDIF 
            ENDIF 
         ELSE 
!
!  ------------  A new electron has been found; add it to the list
!
            MAXORB = K 
            EL(MAXORB) = EL1 
            READ (STRING(J+4:J+7), '(F4.0)') SUM(K) 
         ENDIF 
         J = J + 10 
         IF (J < 50) GO TO 16 
      ENDIF 
!
!  -----  The list of electrons has been determined
!
      WRITE (0, 19) MAXORB, (EL(J),J=1,MAXORB) 
   19 FORMAT(/,' There are ',I3,' orbitals as follows:'/(1X,18(1X,A3))) 
      NWF = MAXORB 
      IF (NIT < 0) NIT = NWF 
   21 CONTINUE 
      WRITE (0, '(/A,A)') ' Orbitals to be varied: ', &
         'ALL/NONE/=i (last i)/comma delimited list/H' 
      READ (5, '(A)') STRING 
      IF (STRING(1:1)=='h' .OR. STRING(1:1)=='H') THEN 
         CALL HELP (1) 
         GO TO 21 
      ELSE IF (STRING(1:3)=='ALL' .OR. STRING(1:3)=='all') THEN 
         NIT = NWF 
      ELSE IF (STRING(1:4)=='NONE' .OR. STRING(1:4)=='none') THEN 
         NIT = 0 
      ELSE IF (INDEX(STRING,'=') /= 0) THEN 
         J = INDEX(STRING,'=') 
         JJ = INDEX(STRING,' ') 
         READ (STRING(J+1:), '(I2)') NIT 
         IF (JJ == J + 2) NIT = MOD(NIT,10) 
      ELSE 
         NIT = 0 
         J = 1 
   22    CONTINUE 
         NEXT = INDEX(STRING(J:),',') 
!
!        ***  Search for last electron label which need not be followed
!             by a comma
!
         IF (NEXT==0 .AND. STRING(J:J+2)/='   ') NEXT = INDEX(STRING(J+1:),' ')&
             + 1 
         IF (NEXT >= 1) THEN 
            IF (NEXT == 4) THEN 
               EL1 = STRING(J:J+2) 
            ELSE IF (NEXT == 3) THEN 
               EL1 = ' '//STRING(J:J+1) 
            ELSE 
               WRITE (0, *) 'Electron labels must be separated by commas;' 
               WRITE (0, *) ' each label must contain 2 or 3 characters' 
               GO TO 21 
            ENDIF 
            CALL REORD (EL, EL1, NWF, IERR) 
            IF (IERR == 0) THEN 
               NIT = NIT + 1 
               J = J + NEXT 
               IF (J < 72) GO TO 22 
            ELSE 
               WRITE (0, *) ' Case must match as well as position of', &
                  ' imbedded blanks' 
               WRITE (0, *) ' For 3rd character of label to be blank', &
                  ' follow blank with comma' 
 
               GO TO 21 
            ENDIF 
         ENDIF 
      ENDIF 
!
      IB = NWF - NIT + 1 
      IF (NIT /= 0) THEN 
   23    CONTINUE 
         WRITE (0, '(/A)') ' Default electron parameters ? (Y/N/H) ' 
         READ(5,'(A)') ANS
         IF (ANS=='H' .OR. ANS=='h') THEN 
            CALL HELP (2) 
            GO TO 23 
         ENDIF 
      ELSE 
         ANS = 'Y' 
      ENDIF 
      IF (ANS/='Y' .AND. ANS/='y') WRITE (0, '(A,A)') &
         ' S, IND, METH, ACC for non-closed Shell electrons: ' 
      DO I = NCLOSD + 1, NWF 
         IF (ANS=='Y' .OR. ANS=='y') THEN 
            S(I) = SS + (SUM(I)-D1)/D2 
            SS = SS + SUM(I) 
            METH(I) = 1 
            ACC(I) = D0 
            IND(I) = 0 
            IF (IUF /= 0) IND(I) = -1 
         ELSE 
            WRITE (0, '(A,A)') EL(I), ':  ' 
            READ (5, *) S(I), IND(I), METH(I), ACC(I) 
         ENDIF 
         J = 2 
         IF (EL(I)(1:1) == ' ') J = 3 
         L(I) = LVAL(EL(I)(J:J))
         N(I) = ICHAR(EL(I)(J-1:J-1)) - ICHAR('1') + 1
         IF (IND(I) == 1) CYCLE  
         EK(I) = D0 
         AZ(I) = D0 
      END DO 
!
!  *****  DEFINE ALL ORBITALS IN THE CONFIGURATION TO BE ORTHOGONAL
!
      DO I = 1, NWF 
         E(I,I) = D0 
         E(I,:I-1) = D0 
         WHERE (L(I) == L(:I-1))  
            E(I,:I-1) = 1.D-5 
         END WHERE 
         E(:I-1,I) = E(I,:I-1) 
      END DO 
      IB = NWF - NIT + 1 
      NO = NOD 
      ND = NO - 2 
      STRONG = .FALSE. 
      WRITE (3, 62) ATOM, TERM, Z, (EL(I),INT(SUM(I)),I=1,NCLOSD) 
   62 FORMAT('1'/,/,/,9X,'HARTREE-FOCK WAVE FUNCTIONS FOR  ',2A6,' Z =',F5.1,/,&
         /,14X,'Core =',5(1X,A3,'(',I4,')'),/(20X,5(1X,A3,'(',I4,')'))) 
      WRITE (3, '(5X,A15,A50)') 'Configuration =', CONFIG 
      WRITE (3, 71) 
   71 FORMAT(/,/,9X,'INPUT DATA'/,9X,'----- ----'/,/,13X,'WAVE FUNCTION',&
         '  PROCEDURE'/,17X,'NL  SIGMA METH ACC OPT'/,/,/) 
      DO I = 1, NWF 
         WRITE (3, 78) I, EL(I), N(I), L(I), S(I), METH(I), ACC(I), IND(I) 
   78    FORMAT(I8,2X,A3,2I3,F7.1,I4,F4.1,I4) 
      END DO 
      OMIT = .NOT.STRONG 
!
      CALL ARRAY 
      CALL ENEXPR (TERM, DONE) 
      IF (.NOT.DONE) THEN 
!
!  ---  Case needs additional data
!
         WRITE (0, 85) 
   85    FORMAT(/,' The program could not derive the energy expression'/,&
            ' Select one of the following options and enter:'/,&
            '    1  Re-enter the term and configuration'/,&
            '    2  Enter the deviations from Eav as input'/,'    3  STOP'/) 
         READ (5, *) ISELEC 
         GO TO (5,86,99) ISELEC 
   86    CONTINUE 
         WRITE (0, 83) 
   83    FORMAT(/,' Input data for deviations from the average energy'/,&
            ' First FK integrals, then GK integrals in indicated format'/,&
            '  cc.ccccccccccFkk(el1,el2)  - terminate each list with an *',&
            ' in the F column') 
         FIRST = .TRUE. 
!
!  *****  READ 'FK' AND 'GK' CARDS, OMITTING THE HEADER IF A FILE
!
   82    CONTINUE 
         READ (5, 84) CFG, W, KFG, EL1, EL2 
   84    FORMAT(F14.8,A1,I2,1X,A3,1X,A3) 
         IF (W /= ASTER) THEN 
            CALL EPTR (EL, EL1, IFG, J1) 
            IF (J1 == 1) GO TO 99 
            CALL EPTR (EL, EL2, JFG, J2) 
            IF (J2 == 1) GO TO 99 
            CALL ADD (CFG, KFG, IFG, JFG, FIRST) 
            CALL ADD (CFG, KFG, JFG, IFG, FIRST) 
            GO TO 82 
         ELSE IF (FIRST) THEN 
            FIRST = .FALSE. 
            GO TO 82 
         ENDIF 
      ENDIF 
!
!  *****  COMPUTE THE INITIAL ARRAY AND INITIAL RADIAL FUNCTIONS
!
      CALL WAVEFN 
!
!      ... Define an order for the functions to be iterated
!
      DO JP = 1, NWF 
         IORD(JP) = JP 
      END DO 
   91 CONTINUE 
      ORDERD = .TRUE. 
      DO JP = IB, NWF - 1 
         N1 = N(IORD(JP)) 
         L1 = L(IORD(JP)) 
         N2 = N(IORD(JP+1)) 
         L2 = L(IORD(JP+1)) 
         IF (.NOT.(N1>N2 .OR. N1==N2 .AND. L1>L2)) CYCLE  
         ITEMP = IORD(JP) 
         IORD(JP) = IORD(JP+1) 
         IORD(JP+1) = ITEMP 
         ORDERD = .FALSE. 
      END DO 
      IF (.NOT.ORDERD) GO TO 91 
      RETURN  
   99 CONTINUE 
      STOP  
      END SUBROUTINE DATA 
!
!     ------------------------------------------------------------------
!                       D E
!     ------------------------------------------------------------------
!
!       This routine controls the solution of the differenttial equation
!   for the radial function P  .  One of three methods is selected -
!                            I1
!   M1, M2, or M3 -  for solving the equations,  the  initial  choice
!   being determined by an input paramter, METH(I), except when no
!   exchange is present, in which case M2 is selected. (For further
!   information see Sec. 7-4)
!
!        Value of METH(I)     Method
!        ---------------      ------
!        < or =1            M1 with search for an acceptable solution
!             =2            M2 with search for an acceptable solution
!             =3            M3 without any checking
!
!   If M1 fails to find an acceptable solution, the radial  functions
!   are  orthogonalized,  off-diagonal  energy parameters recomputed,
!   and the method tried again.   Should it continue to fail, METH(I)
!   is set to 2.
!
!
      SUBROUTINE DE(I1) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE de_C
      USE TEST_C 
      USE COEFF_C 
      USE PARAM_C 
      USE RADIAL_C 
      USE WAVE_C, ONLY: PDE, E, SUM, DPM, ACC, METH, IPR 
      USE LABEL_C, ONLY: EL 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE methd1_I 
      USE quad_I 
      USE quadr_I 
      USE orthog_I 
      USE grange_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I1 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J, NN, IJ, JJ 
      REAL(DOUBLE) :: ED2, PN, ED1, CD, DP, DIFF, DPW, PNN 
      LOGICAL :: CHANGE 
      CHARACTER(LEN=3), DIMENSION(3) :: ASTER = (/ '  ', '* ', '**'/)
!-----------------------------------------------
!
      I = I1 
      ED2 = E(I,I) 
      KK = MAX0(1,METH(I)) 
      IF (NWF == 1) KK = 2 
      NODE = N(I) - L(I) - 1 

!
!  *****  CALL METHD1 TO SOLVE THE DIFFERENTIAL EQUATION
!
      CALL METHD1 (I) 
      IF (FAIL) GO TO 25 
!
   12 CONTINUE 
      PN = DSQRT(QUAD(I,M,PDE,PDE)) 
      PDE(:M) = PDE(:M)/PN 
      AZZ = AZZ/PN 
!
!  *****  CHECK IF DIFFERENT METHOD SHOULD BE USED
!
      IF (KK == 1) THEN 
         IF (DABS(D1 - ED2/E(I,I))<0.005D0 .AND. DMAX1(DABS(D1-PN),DABS(D1/PN-&
            D1))>0.20D0) THEN 
            METH(I) = 2 
            KK = 2 
            GO TO 25 
         ENDIF 
      ELSE 
         IF (DABS(D1 - ED2/E(I,I))<0.0001D0 .AND. IC>1) IC = IC - 1 
      ENDIF 
!
!  *****  SET THE ACCELERATING PARAMETER
!
      IF (IPR /= I) THEN 
         ACC(I) = 0.75*ACC(I) 
      ELSE 
         ED2 = ED2 - E(I,I) 
         IF (ED1*ED2 > D0) THEN 
            ACC(I) = 0.75*ACC(I) 
         ELSE 
            ACC(I) = (D1 + D3*ACC(I))/D4 
         ENDIF 
      ENDIF 
      C = ACC(I) 
      CD = D1 - C 
!
!   *****  IMPROVE THE ESTIMATES
!
      MAX(I) = M 
      DP = D0 
      DO J = 1, M 
         DIFF = P(J,I) - PDE(J) 
         DP = DMAX1(DP,DABS(DIFF)*R2(J)) 
         P(J,I) = PDE(J) + C*DIFF 
      END DO 
      IF (M /= NO) THEN 
         M = M + 1 
         P(M:NO,I) = D0 
         AZ(I) = CD*AZZ + C*AZ(I) 
         AZZ = AZ(I) 
      ENDIF 
!
!  *****  CHECK THE ORTHOGONALIZATION
!
      NN = NWF 
      IF (OMIT) NN = IB - 1 
      IJ = 0 
      DPW = DP/DSQRT(SUM(I)) 
      M = MAX(I) 
      CHANGE = .FALSE. 
      DO J = 1, NN 
         IF (E(I,J)==D0 .OR. I==J) CYCLE  
         IF (DPM(J)>=DSQRT(SUM(J))*DPW .AND. J>=IB) CYCLE  
!
!        ORTHOGONALITY CONDITION APPLIES
!
         C = QUADR(I,J,0) 
         WRITE (6, 63) EL(J), EL(I), C 
   63    FORMAT(6X,'<',A3,'|',A3,'>=',1P,D8.1) 
         M = MAX0(M,MAX(J)) 
         P(:M,I) = P(:M,I) - C*P(:M,J) 
         AZZ = AZZ - C*AZ(J) 
         CHANGE = .TRUE. 
      END DO 
      IF (CHANGE .OR. C/=D0) THEN 
         PNN = DSQRT(QUADR(I,I,0)) 
         P(:M,I) = P(:M,I)/PNN 
         AZZ = AZZ/PNN 
      ENDIF 
      M = NO 
   67 CONTINUE 
      IF (DABS(P(M,I)) < 1.D-15) THEN 
         P(M,I) = D0 
         M = M - 1 
         GO TO 67 
      ENDIF 
      MAX(I) = M 
      IF (AZZ > D0) AZ(I) = DMAX1(AZZ,D5*AZ(I)) 
      WRITE (6, 17) EL(I), E(I,I), AZ(I), PN, ASTER(KK), DP 
   17 FORMAT(20X,A3,2F15.7,F12.7,A2,1P,D10.2) 
      DPM(I) = DP 
      IF (IPR == I1) THEN 
         ED1 = ED2 
      ELSE 
         ED1 = ED2 - E(I1,I1) 
      ENDIF 
      IPR = I1 
      RETURN  
!
!  *****  IF METHD1 FAILED TO FIND AN ACCEPTABLE SOLUTION, ORTHOGONALIZE
!  *****  THE ESTIMATES AND TRY AGAIN
!
   25 CONTINUE 
      IF (I /= IB) THEN 
         CALL ORTHOG 
         CALL GRANGE 
      ENDIF 
   27 CONTINUE 
      CALL METHD1 (I) 
      IF (FAIL) THEN 
!
!  *****  ERROR RETURN FROM SECOND TRY.  IF M1 WAS USED,SWITCH TO
!         M2 AND TRY ONCE MORE.
!
         IF (KK == 2) RETURN  
         KK = 2 
         GO TO 27 
      ELSE 
         GO TO 12 
      ENDIF 
      RETURN  
      END SUBROUTINE DE 
!
!     ------------------------------------------------------------------
!                       D E V
!     ------------------------------------------------------------------
!
!     Add the deviations to the average energy for a partially filled
!       p- or d- shell
!
      SUBROUTINE DEV(IEL, L, Q, I, DONE) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!...Translated by Pacific-Sierra Research 77to90  4.3E  17:25:54  12/28/06  
!...Switches:                     
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE add_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: IEL 
      INTEGER , INTENT(IN) :: L 
      INTEGER , INTENT(INOUT) :: I 
      REAL(DOUBLE) , INTENT(IN) :: Q 
      LOGICAL , INTENT(OUT) :: DONE 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(6) :: F2PP = (/ -3, 3, 12, -9, 0, 6/)  
      INTEGER , DIMENSION(45) :: F2DD = &
       (/ -58, 77, 50, -13, 140, -93, 42, -12, -57, 123, 105, 69, -12, &
          -105, -69, -24, 66, 12, 39, 21, 57, -51, 30, 48, 84, 219,    &
           111, 210, 138, -175, -85, 23, -22, -112, -76, -58, 167, 23, &
           -85,  59, 140, 104, 86, 320, 113/)    
      INTEGER , DIMENSION(45) :: F4DD = &
       (/ 5, -70, 15, 50, 140, -30, -105, 30, 55, -45, 105, -15, 30, &
         -105, 15, -10, 45, -30, -45, 70, -55, 75, 135, 20, 0, 30, -15,&
          210, -30, -175, -50, -40, -85, 35, 50, 110, -15, -5, 125,    &
          -25, 140, 20, -40, -100, -55/)
      INTEGER :: N 
!-----------------------------------------------
!
      DONE = .TRUE. 
      N = Q 
      IF (N > 2*L + 1) N = 4*L + 2 - N 
      IF (N > 1) THEN 
         IF (L == 1) THEN 
            CALL ADD (2*F2PP(I)/25.D0, 2, IEL, IEL, .TRUE.) 
         ELSE IF (L == 2) THEN 
            I = I - 6 
            CALL ADD (2*F2DD(I)/441.D0, 2, IEL, IEL, .TRUE.) 
            CALL ADD (2*F4DD(I)/441.D0, 4, IEL, IEL, .TRUE.) 
         ELSE 
            DONE = .FALSE. 
         ENDIF 
      ENDIF 
      RETURN  
      END SUBROUTINE DEV 
!
!     ------------------------------------------------------------------
!                              D Y K
!     ------------------------------------------------------------------
!
!       Stores in YK the values of the integral of
!              k
!       P (s/r) (dP /ds - P /s) integrated over the interval (0,r)
!        i         j       j
!
!   which enter into the spin-orbit calculation.
!
!
      SUBROUTINE DYK(I, J, K) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE PARAM_C 
      USE RADIAL_C 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I, J, K 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: JJ, M, MM 
      REAL(DOUBLE) :: DEN, FACT, A, AA, F1, F2, F3, C, HH 
!-----------------------------------------------
      DEN = L(I) + L(J) + 2 + K 
      FACT = L(J) 
      YK(:2) = FACT*P(:2,I)*P(:2,J)*R(:2)/DEN 
      A = EH**K 
      AA = A*A 
      A = D4*A 
      F1 = FACT*P(1,I)*P(1,J)*R(1) 
      F2 = FACT*P(2,I)*P(2,J)*R(2) 
      DO M = 3, ND 
         F3 = ((-P(M+2,J))+D8*(P(M+1,J)-P(M-1,J))+P(M-2,J))/(D6*H) 
         F3 = D5*P(M,I)*(F3 - P(M,J))*R(M) 
         YK(M) = YK(M-2)*AA + H3*(F3 + A*F2 + AA*F1) 
         F1 = F2 
         F2 = F3 
      END DO 
      A = A*EH**3 
      AA = A*A/D16 
      C = 2*K + 3 
      HH = C*H3 
      YK(NO) = YK(ND) 
      F1 = YK(NO) 
      F2 = F1 
      DO MM = 3, NO 
         M = NO - MM + 1 
         F3 = YK(M) 
         YK(M) = YK(M+2)*AA + HH*(F3 + A*F2 + AA*F1) 
         F1 = F2 
         F2 = F3 
      END DO 
      RETURN  
      END SUBROUTINE DYK 
!
!     ------------------------------------------------------------------
!                       E K I N
!     ------------------------------------------------------------------
!
!       Returns the value of the integral of
!
!         (2/r)P (Y P  + X )
!               j  i i    i
!
!   integrated with respect to r.
!
!
      REAL(KIND(0.0D0)) FUNCTION EKIN (I, II, REL) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE RADIAL_C 
      USE PARAM_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE xch_I 
      USE potl_I 
      USE quads_I 
      USE quad_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I 
      INTEGER  :: II 
      LOGICAL  :: REL 
!-----------------------------------------------
      CALL XCH (I, 2) 
      CALL POTL (I, REL) 
      YK(:NO) = YR(:NO) 
      YR(:NO) = P(:NO,II) 
      EKIN = D2*QUADS(I,II,1) + QUAD(II,NO,YR,X) 
      RETURN  
      END FUNCTION EKIN 
!
!     ------------------------------------------------------------------
!                       E N E R G Y
!     ------------------------------------------------------------------
!
!       Determines the position of the electron in the electron list
!
      SUBROUTINE ENERGY(ETOTAL) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE TEST_C 
      USE LABEL_C 
      USE RADIAL_C, ONLY: L, NOD 
      USE WAVE_C, ONLY: EK, SUM 
      USE PARAM_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE hl_I 
      USE a_I 
      USE fk_I 
      USE b_I 
      USE gk_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE) , INTENT(OUT) :: ETOTAL 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J, K 
      REAL(DOUBLE) :: C 
!-----------------------------------------------
!
!  *****  COMPUTE KINETIC ENERGY IF NECESSARY
!
      DO I = 1, NWF 
         EK(I) = -D5*HL(EL,I,I,REL) 
      END DO 
!
      ETOTAL = D0 
      DO I = 1, NWF 
         ETOTAL = ETOTAL + SUM(I)*EK(I) 
         DO J = 1, I 
            DO K = 0, 2*MIN0(L(I),L(J)), 2 
               C = A(I,J,K)*SUM(I) 
               IF (I == J) C = C/D2 
               IF (ABS(C) == D0) CYCLE  
               ETOTAL = ETOTAL + C*FK(I,J,K,REL) 
            END DO 
         END DO 
         DO J = 1, I - 1 
            DO K = ABS(L(I)-L(J)), L(I) + L(J), 2 
               C = B(I,J,K)*SUM(I) 
               IF (ABS(C) == D0) CYCLE  
               ETOTAL = ETOTAL + C*GK(I,J,K,REL) 
            END DO 
         END DO 
      END DO 
      RETURN  
      END SUBROUTINE ENERGY 
!
!     ------------------------------------------------------------------
!                       E N E X P R
!     ------------------------------------------------------------------
!
!     Determine the deviations to the average energy for the following:
!        i) an open p- or d-shell
!       ii) a single electron or hole, any l
!      iii) an s-electron and a single electron, any l
!       iv) an s-electron and an open p- or d-shell
!        v) an open p-shell and a single electron, any l
!
      SUBROUTINE ENEXPR(TERM, DONE) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE PARAM_C 
      USE WAVE_C 
      USE RADIAL_C, ONLY: L 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE looktm_I 
      USE dev_I 
      USE add_I 
      USE lval_I 
      USE lookup_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      LOGICAL  :: DONE 
      CHARACTER , INTENT(IN) :: TERM*6 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER , DIMENSION(5) :: SUMTAB = (/ 1, 4, 7, 10, 11/)  
      INTEGER , DIMENSION(11) :: PARTAB = &
       (/ 2, 3, 1, 1, 4, 2, 2, 3, 1, 1, 2/)
      INTEGER , DIMENSION(11) :: PTRTAB = &
       (/ 6, 12, 17, 18, 20, 30, 36, 42, 47, 48, 54/)  
      CHARACTER, DIMENSION(11) :: PARCH = &
       (/ 'P', 'P', 'D', 'S', 'S', 'D', 'P', 'P', 'D', 'S', 'P'/)  
!
!     ... Encoded term value -- S = LTAB/10
!                         Lterm = L + (LTAB mod 10 - 5)
!         Example:  LTAB = 36 with L = 2  is 3F
      INTEGER , DIMENSION(54) :: LTAB = & 
       (/ 36, 35, 34, 16, 15, 14, 46, 45, 44, 26, 25, 24, 27, 26, 25,&
          24, 23, 25, 55, 35, 37, 36, 35, 34, 33, 17, 16, 15, 14, 13,&
          36, 35, 34, 16, 15, 14, 46, 45, 44, 26, 25, 24, 27, 26, 25,&
          24, 23, 25, 36, 35, 34, 16, 15, 14/)  
      INTEGER , DIMENSION(11) :: PLVAL = &
       (/ 1, 1, 2, 0, 0, 2, 1, 1, 2, 0, 1/)  
      INTEGER :: PACVAL, SP, PS1, PS2 
!
!     ... FINT, GINT1, and GINT2 are coefficients of Slater integrals
!         in l, tabulated by Slater,
!
!
!     ... coefficients of F2 integrals for p(n)l(1) configurations
!
      INTEGER , DIMENSION(3,54) :: FINT = RESHAPE( &
       (/ 2, -1, 0, -4, -4, 3, 2, 5, 3, 2, -1, 0, -4, -4, 3, 2, 5, 3, &
         -2, 1, 0,  4, 4, -3, -2, -5, -3, -2, 1, 0, 4, 4, -3, -2, -5, &
         -3, 4, -2, 0, -2, -11, 6, -4, -4, 15, -2, 7, 15, 4, 10, 6, 0,&
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, &
          0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, -1, 0,   &
         -4, -4, 3, 2, 5, 3, 2, -1, 0, -4, -4, 3, 2, 5, 3, -4, 2, 0,  &
          2, 11, -6, 4, 4, -15, 2, -7, -15, -4, -10, -6, 0, 0, 0, -2, &
          1, 0, 4, 4, -3, -2, -5, -3, -2, 1, 0, 4, 4, -3, -2, -5, -3/)&
          , (/ 3, 54/) )  
!
!     ... coefficients of G(l-1) integrals
!
      INTEGER , DIMENSION(3,54) :: GINT1 = RESHAPE( &
       (/-10, 5, 0, 2, 11, -6, 2, -1, -6, 14, -7, 0, 2, -13, 6, 2, -1,&
          6, -8, 4, 0, -8, 4, 0, 4, 10, 0, 10, -5, 0, 10, -5, 0, 4,   &
         -8, 0, -8, 4, 0, -2, 13, -6, 2, 11, -12, 4, 4, -12, 4, -2,   &
         -12, 0, 0, 0, -6, 3, 0, 10, -5, 0, -6, 3, 0, -6, 3, 0, -4, 8,&
          3, 0, 12, 3, 6, 9, 0, -6, 3, 0, 6, 21, 12, 12, 12, -27, 12, &
         -6, 27, 6, -15, -12, -6, 3, 0, 0, 6, -3, 0, 0, -3, 6, -3, 0, &
          0, -6, 3, 12, 6, 3, -4, 2, 0, -4, 2, 0, -4, 2, 0, -4, 2, 0, &
         14, 11, -9, 14, -7, -9, -4, -2, 0, -4, 2, 0, -2, 7, 3, 2, 11,&
          3, 8, 8, 0, 0, 0, 0, -2, 1, 0, -2, 1, 0, -2, 1, 0, -2, 1, 0,&
         -2, 1, 0, 22, 13, 0/), (/ 3, 54 /) )  
!
!     ... coefficients of G(l+1) integrals
!
      INTEGER , DIMENSION(3,54) :: GINT2 = RESHAPE( &
       (/ 2, 5, -3, 2, -7, -15, -10, -25, -15, 2, 5, 9, 2, 17, 21, 14,&
         35, 21, 4, -2, -6, -8, -20, -12, -8, -20, -12, 4, 16, 12, 10,&
         25, 15, 10, 25, 15, 4, 10, 0, 4, 4, -12, 2, -7, -21, -2, -17,&
        -21, -8, -20, -12, 0, 0, 0, -6, -15, -9, 10, 25, 15, 6, 3, -3,&
          0, -12, -9, -4, -16, -9, -6, -15, -9, -6, -15, -9, 6, 27, 9,&
         12, 30, -9, 12, 12, -27, 6, -9, -27, -6, -15, -9, 0, 0, -3,  &
          0, -6, -9, -6, -15, -9, 12, 18, 9, 0, 6, 9, 6, 15, 9, -4,   &
        -10, -6, -4, -10, -6, -4, -10, -6, 14, 35, 12, 14, 17, -6, -4,&
        -10, -6, 8, 8, 0, 2, -7, -6, -2, -11, -6, -4, -10, -6, -4,    &
        -10, -6, 0, 0, 0, -2, -5, -3, -2, -5, -3, -2, -5, -3, 22, 31, &
          9, -2, -5, -3, -2, -5, -3/), (/ 3, 54/) )  
      INTEGER , DIMENSION(2) :: NOS 
      INTEGER :: IP, IL, IS, J, I, IIS, IIL, NSL, IPM, NSLM, IPP, NSLP,&
                 ISUMP,  NL, LP, IPTR1, IPTR2, NOMACH, IND, LV 
      REAL(DOUBLE) :: C, CSP, VAL1, VAL2, VAL3 
      CHARACTER :: SL*2, SENOR, PSL*2, SLM*2, SLP*2 
!-----------------------------------------------
!
      IP = 1 
    1 CONTINUE 
      IF (TERM(IP:IP) == ' ') THEN 
         IP = IP + 1 
         GO TO 1 
      ENDIF 
      SL = TERM(IP:IP+1) 
      SENOR = ' ' 
      IF (IP <= 4) SENOR = TERM(IP+2:IP+2) 
!
!   ---  convert lowercase L symbol to uppercase
!
      IF (SL(2:2)>'a' .AND. SL(2:2)<'z') SL(2:2) = CHAR(ICHAR(SL(2:2))+ICHAR(&
         'A')-ICHAR('a')) 
!
!  ---  determine if FK or GK data needs to be input
!
      IL = 0 
      IS = 0 
      J = 1 
      IF (SL/='AV' .AND. SL/='aV') THEN 
         DO I = NCLOSD + 1, NWF 
            IF (SUM(I)==4*L(I) + 2 .OR. SUM(I)==0.D0) CYCLE  
            IF (J > 2) THEN 
               DONE = .FALSE. 
               RETURN  
            ENDIF 
            NOS(J) = I 
            J = J + 1 
            IF (L(I)==0 .AND. IS==0) THEN 
               IS = IS + 1 
               IIS = I 
            ELSE 
               IL = IL + 1 
               IIL = I 
            ENDIF 
         END DO 
      ELSE 
         DO I = NCLOSD + 1, NWF 
            IF (SUM(I)==4*L(I) + 2 .OR. SUM(I)==0.D0) CYCLE  
            IF (J > 2) THEN 
               DONE = .TRUE. 
               RETURN  
            ENDIF 
            NOS(J) = I 
            J = J + 1 
            IF (L(I)==0 .AND. IS==0) THEN 
               IS = IS + 1 
               IIS = I 
            ELSE 
               IL = IL + 1 
               IIL = I 
            ENDIF 
         END DO 
      ENDIF 
      IF (SL/='AV' .AND. SL/='aV' .AND. IS+IL/=0) THEN 
         DONE = .FALSE. 
         C = 0.D0 
         IF (IS + IL<=2 .AND. IL<=1) THEN 
            IF (IS==0 .AND. IL==1) THEN 
    3          CONTINUE 
               CALL LOOKTM (L(IIL), SL, SENOR, SUM(IIL), IP, NSL) 
               IF (NSL > 1) THEN 
                  WRITE (0, *) ' Ambiguous term: enter seniority' 
                  READ (5, '(A1)') SENOR 
                  GO TO 3 
               ENDIF 
               CALL DEV (IIL, L(IIL), SUM(IIL), IP, DONE) 
            ELSE IF (IS==1 .AND. IL==1) THEN 
               SLM = SL 
               SLP = SL 
               SLM(1:1) = CHAR(ICHAR(SLM(1:1))-1) 
               SLP(1:1) = CHAR(ICHAR(SLP(1:1))+1) 
               CALL LOOKTM (L(IIL), SLM, SENOR, SUM(IIL), IPM, NSLM) 
               CALL LOOKTM (L(IIL), SLP, SENOR, SUM(IIL), IPP, NSLP) 
               IF (NSLM + NSLP == 0) THEN 
                  DONE = .FALSE. 
                  RETURN  
               ELSE IF (NSLM==1 .AND. NSLP==0) THEN 
                  SL = SLM 
                  IP = IPM 
               ELSE IF (NSLM==0 .AND. NSLP==1) THEN 
                  SL = SLP 
                  IP = IPP 
               ELSE IF (NSLM==1 .AND. NSLP==1) THEN 
    4             CONTINUE 
                  WRITE (0, '(A,A3,A,A3)') ' Ambiguous l**n term: enter', SLM, &
                     ' or ', SLP 
                  READ (5, '(A2)') SL 
                  IF (SL == SLM) THEN 
                     IP = IPM 
                  ELSE IF (SL == SLP) THEN 
                     IP = IPP 
                  ELSE 
                     WRITE (0, *) ' Term not allowed: re-enter' 
                     GO TO 4 
                  ENDIF 
               ELSE 
    5             CONTINUE 
                  WRITE (0, '(A,A)') ' Ambiguous l**n parent term:', &
                     'Enter term and seniority' 
                  READ (5, '(A2,A1)') SL, SENOR 
                  IF (SENOR == ' ') THEN 
                     WRITE (0, *) 'Seniority is needed' 
                     GO TO 5 
                  ENDIF 
                  CALL LOOKTM (L(IIL), SL, SENOR, SUM(IIL), IP, NSL) 
                  IF (NSL /= 1) THEN 
                     WRITE (0, '(A,A3,A,A3,A)') ' Allowed terms are ', SLM, &
                        ' or ', SLP, ' plus seniority' 
                     GO TO 5 
                  ENDIF 
               ENDIF 
               CALL DEV (IIL, L(IIL), SUM(IIL), IP, DONE) 
               IF (DONE) THEN 
                  SP = ICHAR(SL(1:1)) - ICHAR('0') 
                  CSP = (SP - 1)/2. 
                  IF (SL == SLM) THEN 
                     C = -CSP/(2*L(IIL)+1) 
                  ELSE 
                     C = (CSP + 1)/(2*L(IIL)+1) 
                  ENDIF 
                  CALL ADD (C, L(IIL), IIS, IIL, .FALSE.) 
                  CALL ADD (C, L(IIL), IIL, IIS, .FALSE.) 
               ENDIF 
            ELSE IF (IS==1 .AND. IL==0) THEN 
               DONE = .TRUE. 
            ENDIF 
         ELSE 
            IF (L(NOS(1))==1 .AND. SUM(NOS(2))==1.D0 .OR. L(NOS(2))==1 .AND. &
               SUM(NOS(1))==1.D0) THEN 
               IF (L(NOS(1))==1 .AND. SUM(NOS(2))==1.D0) THEN 
                  ISUMP = SUM(NOS(1)) 
                  NP = NOS(1) 
                  NL = NOS(2) 
               ELSE 
                  ISUMP = SUM(NOS(2)) 
                  NP = NOS(2) 
                  NL = NOS(1) 
               ENDIF 
               SP = ICHAR(SL(1:1)) - ICHAR('0') 
               LP = LVAL(SL(2:2)) 
               PS1 = SP + 1 
               PS2 = SP - 1 
               IF (ISUMP == 1) THEN 
                  IPTR1 = 1 
               ELSE 
                  IPTR1 = SUMTAB(ISUMP-1) + 1 
               ENDIF 
               IPTR2 = SUMTAB(ISUMP) 
               NOMACH = 0 
               CALL LOOKUP (PARTAB, IPTR1, IPTR2, IND, NOMACH, PS1) 
               CALL LOOKUP (PARTAB, IPTR1, IPTR2, IND, NOMACH, PS2) 
               PSL(1:1) = CHAR(PARTAB(IND)+ICHAR('0')) 
               PSL(2:2) = PARCH(IND) 
               IF (NOMACH > 1) THEN 
                  WRITE (0, *) ' AMBIGUOUS PARENT CASE' 
   10             CONTINUE 
                  WRITE (0, *) ' ENTER THE SL TERM FOR p(n) SUBSHELL' 
                  READ (5, '(A)') PSL 
                  IF (PSL(2:2)>'a' .AND. PSL(2:2)<'z') PSL(2:2) = CHAR(ICHAR(&
                     PSL(2:2))+ICHAR('A')-ICHAR('a')) 
                  PS1 = ICHAR(PSL(1:1)) - ICHAR('0') 
                  PS2 = LVAL(PSL(2:2)) 
                  CALL LOOKUP (PLVAL, IPTR1, IPTR2, IND, NOMACH, PS2) 
                  IF (NOMACH/=1 .AND. PARTAB(IND)/=PS1) GO TO 10 
               ENDIF 
               IF (ISUMP == 1) THEN 
                  IPTR1 = 1 
               ELSE 
                  IPTR1 = PTRTAB(IND-1) + 1 
               ENDIF 
               IPTR2 = PTRTAB(IND) 
               LV = L(NL) 
               PACVAL = SP*10 + LP - LV + 5 
               NOMACH = 0 
               CALL LOOKUP (LTAB, IPTR1, IPTR2, IND, NOMACH, PACVAL) 
               IF (NOMACH /= 1) THEN 
                  DONE = .FALSE. 
                  RETURN  
               ENDIF 
               VAL1 = ((FINT(1,IND)*LV+FINT(2,IND))*LV+FINT(3,IND))/(5.D0*(2*LV&
                   - 1)*(2*LV + 3)) 
               VAL2 = ((GINT1(1,IND)*LV+GINT1(2,IND))*LV+GINT1(3,IND))/(2.D0*(2&
                  *LV + 1)*(2*LV - 1)**2) 
               VAL3 = ((GINT2(1,IND)*LV+GINT2(2,IND))*LV+GINT2(3,IND))/(2.D0*(2&
                  *LV + 1)*(2*LV + 3)**2) 
!
!     ...  Add contributions from between p-subshell and l-electron
!
               CALL ADD (VAL1, 2, NP, NL, .TRUE.) 
               CALL ADD (VAL1, 2, NL, NP, .TRUE.) 
               CALL ADD (VAL2, LV - 1, NP, NL, .FALSE.) 
               CALL ADD (VAL2, LV - 1, NL, NP, .FALSE.) 
               CALL ADD (VAL3, LV + 1, NP, NL, .FALSE.) 
               CALL ADD (VAL3, LV + 1, NL, NP, .FALSE.) 
!
!     ... Add deviations for p-subshell
!
               CALL LOOKTM (1, PSL, ' ', SUM(NP), IP, NSL) 
               CALL DEV (NP, 1, SUM(NP), IP, DONE) 
            ELSE 
               DONE = .FALSE. 
            ENDIF 
         ENDIF 
      ELSE 
         DONE = .TRUE. 
      ENDIF 
      RETURN  
      END SUBROUTINE ENEXPR 
!
!     ------------------------------------------------------------------
!                       E P T R
!     ------------------------------------------------------------------
!
!       Determines the position of the electron in the electron list
!
      SUBROUTINE EPTR(EL, ELSYMB, IEL, J2) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE PARAM_C 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(OUT) :: IEL 
      CHARACTER , INTENT(IN) :: ELSYMB*3 
      CHARACTER , INTENT(IN) :: EL(*)*3 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J2, I, J1 
      CHARACTER(LEN=3) :: BL = '   '
!-----------------------------------------------
      J2 = 0 
!
! ***** SEARCH ELECTRON LIST FOR LSYMB
!
      IF (ELSYMB == BL) THEN 
         IEL = 0 
         RETURN  
      ENDIF 
      DO I = 1, NWF 
         IF (EL(I) /= ELSYMB) CYCLE  
         IEL = I 
         RETURN  
      END DO 
      IEL = -1 
      WRITE (0, 20) ELSYMB 
   20 FORMAT(/,10X,A3,' NOT FOUND IN ELECTRON LIST') 
      J2 = 1 
      RETURN  
      END SUBROUTINE EPTR 
!
!     -----------------------------------------------------------------
!           F A C T R L
!     -----------------------------------------------------------------
!
!
      SUBROUTINE FACTRL(NFACT) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE FACT_C 
!
!      GAM(I) = LOG( GAMMA(I-1) ), WHERE GAMMA(I) = FACTORIAL I-1
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NFACT 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
      REAL(DOUBLE) :: ZERO = 0.D0, ONE = 1.D0 , TWO = 2.d0 , GAMMA, X 
!-----------------------------------------------
!
      GAMMA = ONE 
      GAM(1) = ZERO 
      DO I = 1, NFACT - 1 
         GAMMA = I*GAMMA 
         GAM(I+1) = DLOG(GAMMA) 
      END DO 
      DO I = NFACT + 1, 100 
         X = I - 1 
         GAM(I) = GAM(I-1) + DLOG(X) 
      END DO 
      RETURN  
      END SUBROUTINE FACTRL 
!
!     ------------------------------------------------------------------
!                 F K
!     ------------------------------------------------------------------
!                             k
!       Returns the value of F (i,j)
!
!
      REAL(KIND(0.0D0)) FUNCTION FK (I, J, K, REL) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ykf_I 
      USE quads_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I, J, K 
      LOGICAL  :: REL 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CALL YKF (I, I, K, REL) 
      FK = QUADS(J,J,1) 
      RETURN  
      END FUNCTION FK 
!
!     ------------------------------------------------------------------
!                 G K
!     ------------------------------------------------------------------
!                             k
!       Returns the value of G (i,j).
!
!
      REAL(KIND(0.0D0)) FUNCTION GK (I, J, K, REL) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ykf_I 
      USE quads_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I, J, K 
      LOGICAL  :: REL 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!-----------------------------------------------
      CALL YKF (I, J, K, REL) 
      GK = QUADS(I,J,1) 
      RETURN  
      END FUNCTION GK 
!
!     ------------------------------------------------------------------
!               G R A N G E
!     ------------------------------------------------------------------
!
!       Controls the calculation of off-diagonal energy parameters.
!   It searches for all pairs (i,j) which are constrained through  an
!   orthogonality requirement.   Eq. (7-10) is used to calculate the
!   parameter.
!
      SUBROUTINE GRANGE 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE TEST_C 
      USE PARAM_C 
      USE RADIAL_C 
      USE WAVE_C, ONLY: E, SUM 
      USE LABEL_C, ONLY: EL 
      USE COEFF_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rotate_I 
      USE hl_I 
      USE ekin_I 
      USE a_I 
      USE b_I 
      USE rk_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, J, II, K, KK 
      REAL(DOUBLE) :: C, RES 
!-----------------------------------------------
!
!  *****  ROTATE PAIRS CONNECTED BY ORTHOGONALITY BUT NOT WHEN ONE OF
!         THE ORBITALS IS SIMULTANEOUSLY ORTHOGONAL TO A NON-ORTHOGONAL
!         PAIR
!
      DO I = IB, NWF - 1 
         DO J = I + 1, NWF 
            IF (DABS(E(I,J)) <= 1.D-10) CYCLE  
            CALL ROTATE (I, J) 
         END DO 
      END DO 
!
!  *****   COMPUTE OFF-DIAGONAL ENERGY PARAMETERS
!
      DO I = MAX0(2,IB), NWF 
         DO J = 1, I - 1 
            IF (DABS(E(I,J)) > 1.D-10) THEN 
               IF (J < IB) THEN 
                  E(I,J) = HL(EL,I,J,REL) - EKIN(I,J,REL) 
                  E(J,I) = D0 
               ELSE IF (SUM(I) == SUM(J)) THEN 
                  C = HL(EL,I,J,REL) - (EKIN(I,J,REL) + EKIN(J,I,REL))/D2 
                  E(I,J) = C 
                  E(J,I) = C 
               ELSE 
                  RES = D0 
                  DO II = 1, NWF 
                     IF (II==I .OR. II==J) THEN 
                        DO K = 0, 2*L(I), 2 
                           IF (II == I) THEN 
                              C = A(I,I,K) - A(J,I,K) - B(J,I,K) 
                              IF (DABS(C) > 1.D-10) RES = RES + C*RK(I,I,I,J,K,&
                                 REL) 
                           ELSE IF (II == J) THEN 
                              C = A(J,J,K) - A(I,J,K) - B(I,J,K) 
                              IF (DABS(C) > 1.D-10) RES = RES - C*RK(J,J,J,I,K,&
                                 REL) 
                           ENDIF 
                        END DO 
                     ELSE 
                        DO K = 0, 2*MIN0(L(I),L(II)), 2 
                           C = A(I,II,K) - A(J,II,K) 
                           IF (DABS(C) > 1.D-10) RES = RES + C*RK(I,II,J,II,K,&
                              REL) 
                           KK = ABS(L(I)-L(II)) + K 
                           C = B(I,II,KK) - B(J,II,KK) 
                           IF (DABS(C) <= 1.D-10) CYCLE  
                           RES = RES + C*RK(I,II,II,J,KK,REL) 
                        END DO 
                     ENDIF 
                  END DO 
                  E(I,J) = D2*SUM(J)*RES/(SUM(I)-SUM(J)) 
                  E(J,I) = SUM(I)*E(I,J)/SUM(J) 
               ENDIF 
            ENDIF 
            IF (DABS(E(I,J)) <= 1.D-10) CYCLE  
            WRITE (6, 35) EL(I), EL(J), E(I,J), EL(J), EL(I), E(J,I) 
   35       FORMAT(7X,2(3X,'E(',2A3,') =',F12.5)) 
         END DO 
      END DO 
      RETURN  
      END SUBROUTINE GRANGE 
!
!     ------------------------------------------------------------------
!               H E L P
!     ------------------------------------------------------------------
!
!       Provide HELP information about the data requested
!
      SUBROUTINE HELP(CASE) 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: CASE 
!-----------------------------------------------
!
      SELECT CASE (CASE)  
!
      CASE DEFAULT 
         WRITE (0, 11) 
!  ***** Which orbitals varied?
   11    FORMAT(/,/,1X,'Response  ALL will vary all orbitals'/,11X,&
            'NONE will not vary any orbitals'/,11X,&
            '=n (integer n), will vary last n'/,11X,&
            'with comma delimited list will vary',1X,&
            'only the orbitals in the list'/,/) 
         RETURN  
 
      CASE (2)  
         WRITE (0, 21) 
!  ***** Default electron parameters ?
   21    FORMAT(/,/,1X,'Response  N  will prompt the user for:'/,/,7X,&
            'S     : Screening parameter  (Real number) '/,7X,&
            'IND   : Indicator specifying the type of initial estimate'/,15X,&
            '0 - Screened hydrogenic functions'/,14X,&
            '-1 - Search for functions in wavefunction file; if'/,20X,&
            'not present use screened hydrogenic.'/,7X,&
            'METH  : Method for solving the differential equation'/,15X,&
            '1 - Method 1 solves the boundary value problem for an'/,20X,&
            'acceptable solution which need not be normalized.'/,15X,&
            '2 - Method 2 solves the boundary value problem for '/,20X,&
            'an acceptable solution which is normalized to first'/,20X,&
            'order. If the exchange function is identically zero the '/,20X,&
            'program will automatically select Method 2.'/,15X,&
            '3 - Method 3 is similar to Method 1 but omits all checks'/,20X,&
            'for acceptability.'/,7X,'ACC   : Inital accelerating factor '/,18X&
            ,'( Real number such that 0 .LE. ACC .LT. 1 )'/,/) 
         RETURN  
 
      CASE (3)  
         WRITE (0, 31) 
!  ***** Default values (NO,STRONG) ?
   31    FORMAT(/,/,1X,'Response  Y will set default values as follows'/,/,10X,&
            'NO - Maximum number of points in the range of the '/,15X,&
            'function is set to 220.'/,10X,&
            'Strong - is logically set to .FALSE.'/,/,1X,&
            'Response N  will prompt the user for'/,/,10X,&
            'NO - Maximum number of points in the range of '/,15X,&
            'the function which should be a positive integer'/,15X,&
            'from 160 for a small atom to 220 for a large atom.'/,10X,&
            'Strong - may be set to .TRUE. or .FALSE by user.'/,/) 
         RETURN  
 
      CASE (4)  
         WRITE (0, 41) 
!  ***** Default values for remaining parameters ?
   41    FORMAT(/,/,1X,1X,'Response  Y --Sets the following default values'/,/,&
            10X,'PRINT=.FALSE.'/,10X,'SCFTOL=1.D-8'/,10X,'NSCF=12'/,10X,&
            'IC=2 + (NWF + 1 - IB)/4'/,10X,'TRACE=.FALSE.'/,/,1X,&
            'Response  N --prompts user for new parameter values.'/,/) 
         RETURN  
 
      CASE (5)  
         WRITE (0, 51) 
!  ***** Default values for PRINT, SCFTOL ?
   51    FORMAT(/,/,1X,1X,&
            'Response  Y --Default value for PRINT is .FALSE. thus'/,16X,&
            'radial functions are NOT printed.'/,/,16X,&
            'SCFTOL -the initial value of the parameter defining the'/,16X,&
            'self-consistency tolerance for radial functions is set'/,16X,&
            'to a default value of  1.D-8 .'/,/,1X,&
            'Response  N --prompts user for new values of PRINT and ',&
            'SCFTOL .'/,/) 
         RETURN  
 
      CASE (6)  
         WRITE (0, 61) 
!  ***** Default values for NSCF,  IC ?
   61    FORMAT(/,/,1X,1X,&
            'Response  Y --NSCF, the maximum number of cycles for'/,16X,&
            'the SCF process is set to default value of 12 and'/,&
            'IC is set to 2 + (NWF + 1 -IB)/4 .'/,/,1X,&
            'Response  N --user prompted for new NSCF and IC values'/,/) 
         RETURN  
 
      CASE (7)  
         WRITE (0, 71) 
!  ***** Default values for TRACE ?
   71    FORMAT(/,/,1X,&
            'Response  Y --sets trace to default value of .FALSE. thus a'/,16X,&
            'trace of energy adjustment will  NOT  be printed.'/,/,1X,&
            'Response  N --prompts the user for the new value of trace.'/,16X,&
            'If .TRUE. a trace will be printed showing the energy ',&
            'adjustment'/,16X,'process used by METHD1 for finding an',&
            ' acceptable solution'/,16X,'with the correct number of',' nodes.'/&
            ,/) 
         RETURN  
 
      CASE (8)  
         WRITE (0, 81) 
!  ***** Additional parameters ?
   81    FORMAT(/,/,1X,'Response  Y --additional values may be computed for '/,&
            16X,'SLATER OR MAGNETIC INTEGRALS'/,16X,&
            'EXPECTATION VALUES OF R**K'/,/,16X,&
            'ELECTRON DENSITY AT THE NUCLEUS'/,16X,'SPIN-ORBIT ','PARAMETER'/,&
            16X,'TRANSITION INTEGRALS'/,/,1X,&
            'Response  N --additional parameter computation is ','skipped'/,/) 
         RETURN  
 
      CASE (9)  
         WRITE (0, 91) 
!  ***** Do you wish to continue along the sequence ?
   91    FORMAT(/,/,1X,'Response  Y --sequence is continued'/,/,1X,&
            'Response  N --current case is ended, but new case may be'/,16X,&
            'started'/,/) 
         RETURN  
 
      CASE (10)  
         WRITE (0, 101) 
!  ***** Do you wish to continue ?
  101    FORMAT(/,/,1X,&
            'Response  Y --prompts user for additional iterations and'/,16X,&
            'new IC.  Then performs additional NSCF self-consistent'/,16X,&
            'field iterations.'/,/,1X,&
            'Response  N --terminates the calculation.'/,/) 
         RETURN  
      END SELECT 
      END SUBROUTINE HELP 
!
!     ------------------------------------------------------------------
!                 H L
!     ------------------------------------------------------------------
!
!       Returns the value of <i^L^j>, using a special formula to
!  preserve symmetry.
!
      REAL(KIND(0.0D0)) FUNCTION HL (EL, I, J, REL) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE PARAM_C 
      USE RADIAL_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE rlshft_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I 
      INTEGER  :: J 
      LOGICAL , INTENT(IN) :: REL 
      CHARACTER , INTENT(IN) :: EL(*)*3 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LI, MM, K 
      REAL(DOUBLE) :: C, A1, A2, A3, ZR, DI1, DI2, DJ1, DJ2, DI4, DI3, DI5, DI6&
         , DJ4, DJ3, DJ5, DJ6, TZ, HL2 
!-----------------------------------------------
      IF (IABS(L(I)-L(J)) /= 0) THEN 
         WRITE (0, 4) EL(I), L(I), EL(J), L(J) 
    4    FORMAT(10X,'UNALLOWED L VALUES OCCURRED IN HL SUBROUTINE'/,2(10X,A3,&
            ' HAS L = ',I3)) 
         STOP  
      ENDIF 
      LI = L(I) 
      C = 2*LI + 1 
      A1 = -D2/(C*(LI + 1)) 
      A2 = A1/((C + D2)*(LI + 1)) 
      A3 = A2/((LI + 2)*(LI + 1)) 
      ZR = Z*R(1) 
      HL = H*C*P(1,I)*P(1,J)*(D1 + ZR*(A1 + ZR*(A2 + ZR*A3))) 
      MM = MIN0(MAX(I)+3,MAX(J)+3,ND-1) 
      K = 2 
      C = D4/D3 
      DI1 = P(K+1,I) - P(K-1,I) 
      DI2 = P(K+1,I) - D2*P(K,I) + P(K-1,I) 
      DJ1 = P(K+1,J) - P(K-1,J) 
      DJ2 = P(K+1,J) - D2*P(K,J) + P(K-1,J) 
      HL = HL + DI1*DJ1 + C*DI2*DJ2 
      DO K = 4, MM, 2 
         DI1 = P(K+1,I) - P(K-1,I) 
         DI2 = P(K+1,I) - D2*P(K,I) + P(K-1,I) 
         DI4 = P(K+2,I) - D4*(P(K+1,I)+P(K-1,I)) + D6*P(K,I) + P(K-2,I) 
         DI3 = P(K+2,I) - P(K-2,I) - D2*DI1 
         DI5 = P(K+3,I) - P(K-3,I) - D4*(P(K+2,I)-P(K-2,I)) + 5.D0*(P(K+1,I)-P(&
            K-1,I)) 
         DI6 = P(K+3,I) + P(K-3,I) - D6*(P(K+2,I)+P(K-2,I)) + 15.D0*(P(K+1,I)+P&
            (K-1,I)) - 20.D0*P(K,I) 
         DJ1 = P(K+1,J) - P(K-1,J) 
         DJ2 = P(K+1,J) - D2*P(K,J) + P(K-1,J) 
         DJ4 = P(K+2,J) - D4*(P(K+1,J)+P(K-1,J)) + D6*P(K,J) + P(K-2,J) 
         DJ3 = P(K+2,J) - P(K-2,J) - D2*DJ1 
         DJ5 = P(K+3,J) - P(K-3,J) - D4*(P(K+2,J)-P(K-2,J)) + 5.D0*(P(K+1,J)-P(&
            K-1,J)) 
         DJ6 = P(K+3,J) + P(K-3,J) - D6*(P(K+2,J)+P(K-2,J)) + 15.D0*(P(K+1,J)+P&
            (K-1,J)) - 20.D0*P(K,J) 
         HL = HL + DI1*DJ1 + C*DI2*DJ2 + (DI3*DJ3 + DI2*DJ4 + DI4*DJ2)/45.D0 - &
            (DI3*DJ5 + DI5*DJ3)/252.D0 - (DI2*DJ6 + DI6*DJ2 - 1.1*DI4*DJ4)/&
            378.D0 
      END DO 
      TZ = Z + Z 
      C = (LI + D5)**2 
      HL2 = D5*(TZ*R(1)-C)*P(1,I)*P(1,J) 
      HL2 = HL2 + SUM(D2*(TZ*R(2:MM:2)-C)*P(2:MM:2,I)*P(2:MM:2,J)+(TZ*R(3:MM+1:&
         2)-C)*P(3:MM+1:2,I)*P(3:MM+1:2,J)) 
      HL = (-HL/(D2*H)) + HL2*H1 
      IF (REL) HL = HL - D2*RLSHFT(I,J) 
 
      RETURN  
      END FUNCTION HL 
!
!     ------------------------------------------------------------------
!               H N O R M
!     ------------------------------------------------------------------
!
!       Returns the value of the normalization constant for an (nl)
!   hydrogenic function with nuclear charge ZZ.
!
!
      REAL(KIND(0.0D0)) FUNCTION HNORM (N, L, ZZ) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE PARAM_C 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: L 
      REAL(DOUBLE) , INTENT(IN) :: ZZ 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: M, I 
      REAL(DOUBLE) :: A, B, T, D 
!-----------------------------------------------
      M = L + L + 1 
      A = N + L 
      B = M 
      T = A 
      D = B 
      M = M - 1 
      IF (M /= 0) THEN 
         DO I = 1, M 
            A = A - D1 
            B = B - D1 
            T = T*A 
            D = D*B 
         END DO 
      ENDIF 
      HNORM = DSQRT(ZZ*T)/(N*D) 
      RETURN  
      END FUNCTION HNORM 
!    MCHF_HF (Part2 of 2)
!     ------------------------------------------------------------------
!               H W F
!     ------------------------------------------------------------------
!
!       Returns the value of an unnormalized (nl) hydrogenic function
!   with nuclear charge ZZ and radius r.
!
!
      REAL(KIND(0.0D0)) FUNCTION HWF (N, L, ZZ, R) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE PARAM_C 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: N 
      INTEGER , INTENT(IN) :: L 
      REAL(DOUBLE) , INTENT(IN) :: ZZ 
      REAL(DOUBLE) , INTENT(IN) :: R 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, I 
      REAL(DOUBLE) :: P, A, B, C, X 
!-----------------------------------------------
      K = N - L - 1 
      P = D1 
      A = D1 
      B = K 
      C = N + L 
      X = -D2*ZZ*R/N 
!
!  *****  TEST IF UNDERFLOW MAY OCCUR, IF SO SET HWF = 0
!
      IF (X >= (-150.D0)) THEN 
         IF (K >= 0) THEN 
            IF (K /= 0) THEN 
               DO I = 1, K 
                  P = D1 + A/B*P/C*X 
                  A = A + D1 
                  B = B - D1 
                  C = C - D1 
               END DO 
            ENDIF 
            HWF = P*DEXP(X/D2)*(-X)**(L + 1) 
            RETURN  
         ENDIF 
         WRITE (0, 7) N, L, ZZ, R 
    7    FORMAT(' FORBIDDEN COMBINATION OF N AND L IN HWF SUBPROGRAM'/,' N =',&
            I4,'   L =',I4,'   Z =',F6.1,'   R =',F8.4) 
         STOP  
      ENDIF 
      HWF = D0 
      RETURN  
      END FUNCTION HWF 
!
!     ------------------------------------------------------------------
!               I N I T
!     ------------------------------------------------------------------
!
!       Initializes basic constants of the program including those
!   which define the average energy of a configuration.
!
!
      SUBROUTINE INIT 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE EAV_C 
      USE FACT_C 
      USE PARAM_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE factrl_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
!
!  *****  SET THE COMMONLY USED DOUBLE PRECISION CONSTANTS
!
      D0 = 0.D0 
      D1 = 1.D0 
      D2 = 2.D0 
      D3 = 3.D0 
      D4 = 4.D0 
      D5 = 1.D0/2.D0 
      D6 = 6.D0 
      D8 = 8.D0 
      D10 = 10.D0 
      D12 = 12.D0 
      D16 = 16.D0 
      D30 = 30.D0 
!
!  ***** Set the factorial needed by RME
!
      CALL FACTRL (32) 
!
!  *****  SET FINE STRUCTURE CONSTANT
!
      FINE = 0.25D0/137.036**2 
!
!  *****  SET THE STARTING POINT, STEP SIZE, AND RELATED PARAMETERS
!
      RHO = -4.D0 
      H = 1./16.D0 
      H1 = H/1.5 
      H3 = H/3. 
      CH = H*H/12. 
      EH = DEXP((-H)) 
      NO = 220 
      ND = NO - 2 
!
!  *****  AVERAGE INTERACTIONS FOR EQUIVALENT ELECTRONS
!
!  *****  P - P
!
      CCA(1) = 2.D0/25.D0 
!
!  *****  D - D
!
      CCA(2) = 2.D0/63.D0 
      CCA(3) = 2.D0/63.D0 
!
!  *****  F - F
!
      CCA(4) = 4.D0/195.D0 
      CCA(5) = 2.D0/143.D0 
      CCA(6) = 100.D0/5577.D0 
!
!  *****  G - G
!
      CCA(7) = 20.D0/1309.D0 
      CCA(8) = 162.D0/17017.D0 
      CCA(9) = 20.D0/2431.D0 
      CCA(10) = 4410.D0/371943.D0 
!
!
!  ***** AVERAGE INTERACTIONS FOR NON-EQUIVALENT ELECTRONS
!
!  *****  S - ( S, P, D, F, G )
!
      CCB(1) = 1.D0/2.D0 
      CCB(2) = 1.D0/6.D0 
      CCB(3) = 1.D0/10.D0 
      CCB(4) = 1.D0/14.D0 
      CCB(5) = 1.D0/18.D0 
!
!  *****  P - ( P, D, F, G )
!
      CCB(6) = 1.D0/6.D0 
      CCB(7) = 1.D0/15.D0 
      CCB(8) = 1.D0/15.D0 
      CCB(9) = 3.D0/70.D0 
      CCB(10) = 3.D0/70.D0 
      CCB(11) = 2.D0/63.D0 
      CCB(12) = 2.D0/63.D0 
      CCB(13) = 5.D0/198.D0 
!
!  *****  D - ( D, F, G )
!
      CCB(14) = 1.D0/10.D0 
      CCB(15) = 1.D0/35.D0 
      CCB(16) = 1.D0/35.D0 
      CCB(17) = 3.D0/70.D0 
      CCB(18) = 2.D0/105.D0 
      CCB(19) = 5.D0/231.D0 
      CCB(20) = 1.D0/35.D0 
      CCB(21) = 10.D0/693.D0 
      CCB(22) = 5.D0/286.D0 
!
!  *****  F - ( F, G )
!
      CCB(23) = 1.D0/14.D0 
      CCB(24) = 2.D0/105.D0 
      CCB(25) = 1.D0/77.D0 
      CCB(26) = 50.D0/3003.D0 
      CCB(27) = 2.D0/63.D0 
      CCB(28) = 1.D0/77.D0 
      CCB(29) = 10.D0/1001.D0 
      CCB(30) = 35.D0/2574.D0 
!
!  *****  G - ( G )
!
      CCB(31) = 1.D0/18.D0 
      CCB(32) = 10.D0/693.D0 
      CCB(33) = 9.D0/1001.D0 
      CCB(34) = 10.D0/1287.D0 
      CCB(35) = 245.D0/21879.D0 
      RETURN  
      END SUBROUTINE INIT 
! 
!     ------------------------------------------------------------------
!                       L O O K - T M
!     ------------------------------------------------------------------
!
!     Add the deviations to the average energy for a partially filled
!       p- or d- shell
!
      SUBROUTINE LOOKTM(L, SL, SEN, Q, IP, NSL) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: L 
      INTEGER , INTENT(OUT) :: IP 
      INTEGER , INTENT(OUT) :: NSL 
      REAL(DOUBLE) , INTENT(IN) :: Q 
      CHARACTER , INTENT(IN) :: SL*2 
      CHARACTER , INTENT(IN) :: SEN 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I1, J2 
      INTEGER , DIMENSION(5) :: IPTR = (/ 6, 11, 19, 35, 51 /)
      INTEGER :: N, IBEGIN, IEND, I 
      CHARACTER(LEN=3), DIMENSION(51) :: TERMS = (/ &
       '3P2', '1D2', '1S0', '4S3', '2D3', '2P1', '3F2', '3P2', '1G2',&
       '1D2', '1S0', '4F3', '4P3', '2H3', '2G3', '2F3', '2D1', '2D3',&
       '2P3', '5D4', '3H4', '3G4', '3F2', '3F4', '3D4', '3P2', '3P4',&
       '1I4', '1G2', '1G4', '1F4', '1D2', '1D4', '1S0', '1S4', '6S5',& 
       '4G5', '4F3', '4D5', '4P3', '2I5', '2H3', '2G3', '2G5', '2F3',&
       '2F5', '2D1', '2D3', '2D5', '2P3', '2S5'/)  
 
!
!  --- search for a partially unfilled p- or d-shell
!
      N = Q 
      IF (N > 2*L + 1) N = 4*L + 2 - N 
      IP = 0 
      NSL = 0 
      IF (N>1 .AND. L<=2) THEN 
         IF (L == 1) THEN 
            IBEGIN = 1 
            IEND = 6 
         ELSE 
            IBEGIN = IPTR(N-1) + 1 
            IEND = IPTR(N) 
         ENDIF 
         I = IBEGIN 
         I1 = I 
         J2 = MAX(IEND,I1) 
         DO I = I1, J2 
            IF (SL /= TERMS(I)(1:2)) CYCLE  
            IF (SEN/=' ' .AND. SEN/=TERMS(I)(3:3)) CYCLE  
            NSL = NSL + 1 
            IP = I 
         END DO 
      ELSE IF (N==1 .AND. SL(1:1)=='2') THEN 
         NSL = 1 
      ENDIF 
      RETURN  
      END SUBROUTINE LOOKTM 
!
!     -----------------------------------------------------------------
!                L O O K - U P
!     -----------------------------------------------------------------
!
      SUBROUTINE LOOKUP(TAB, P1, P2, IND, NO, KEY) 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: P1, P2 
      INTEGER , INTENT(OUT) :: IND 
      INTEGER , INTENT(INOUT) :: NO 
      INTEGER , INTENT(IN) :: KEY 
      INTEGER , INTENT(IN) :: TAB(*) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
!-----------------------------------------------
      DO I = P1, P2 
         IF (TAB(I) /= KEY) CYCLE  
         NO = NO + 1 
         IND = I 
      END DO 
      RETURN  
      END SUBROUTINE LOOKUP 
!
!     -----------------------------------------------------------------
!                 L V A L
!     -----------------------------------------------------------------
!
!
      INTEGER FUNCTION LVAL (SYMBOL) 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER , INTENT(IN) :: SYMBOL 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LOCATE 
      CHARACTER(LEN=22) :: SET = 'spdfghiklmnSPDFGHIKLMN'  

      LOCATE = INDEX(SET,SYMBOL) 
      IF (LOCATE <= 11) THEN 
         LVAL = LOCATE - 1 
      ELSE 
         LVAL = LOCATE - 12 
      ENDIF 
      IF (LVAL <0) then
         Write (0,*) 'Symbol ',SYMBOL,' was not found'
         STOP
      END IF
      RETURN  
      END FUNCTION LVAL 
!
!     ------------------------------------------------------------------
!               M E N U
!     ------------------------------------------------------------------
!
!
!     This routine evaluates a variety of atomic parameters as
!     requested by the user.
!
!
      SUBROUTINE MENU 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE LABEL_C 
      USE RADIAL_C, ONLY: AZ, L, NOD 
      USE TEST_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE eptr_I 
      USE quadr_I 
      USE fk_I 
      USE gk_I 
      USE rk_I 
      USE sn_I 
      USE vk_I 
      USE bwzeta_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J1, J2, J3, J4, J5, J6, J7, J8, J9, IFUNC, K, I, I1, I2, I3, &
         I4, LL 
      REAL(DOUBLE) :: RKEV, SI, D, ZETA, ZETACM, TI 
      CHARACTER :: EL1*3, EL2*3, EL3*3, EL4*3, FUNC 
!-----------------------------------------------
!
 
    4 CONTINUE 
      WRITE (0, 5) 
    5 FORMAT(/,/,5X,' These various functions are available:',/,/,10X,&
         '1 - EXPECTATION VALUES OF R**K'/,10X,&
         '2 - SLATER OR MAGNETIC INTEGRALS'/,10X,&
         '3 - ELECTRON DENSITY AT THE NUCLEUS'/,10X,'4 - SPIN-ORBIT PARAMETER'/&
         ,10X,'5 - TRANSITION INTEGRALS'/,10X,'6 - EXIT TO MAIN PROGRAM'/) 
      WRITE (0, '(5X,A)') 'Input number corresponding to your selection:' 
      READ (5, '(I1)') IFUNC 
      GO TO (10,20,30,40,50,60) IFUNC 
 
!  ****  COMPUTE EXPECTATION VALUES
 
   10 CONTINUE 
      WRITE (0, '(/5X,A,/A,T22,A)') &
         'INPUT LABEL FOR ELECTRON FOLLOWED BY k: Example', '  2p  3', &
         'FORMAT(1X,A3,I3) ' 
      READ (5, '(1X,A3,I3)') EL1, K 
      CALL EPTR (EL, EL1, I, J1) 
      IF (J1 == 1) GO TO 10 
      RKEV = QUADR(I,I,K) 
      WRITE (3, 12) EL1, K, EL1, RKEV 
      WRITE (0, 12) EL1, K, EL1, RKEV 
   12 FORMAT(/,15X,' VALUE OF <',A3,'|R**',I2,'|',A3,'> = ',1P,D14.7,' a.u.'/) 
      GO TO 4 
 
!  ****  DETERMINE SLATER INTEGRALS  FK, GK, RK, NK, MK, VK
 
 
   20 CONTINUE 
      WRITE (0, '(/5X,A/A,T22,A)') &
         'INPUT PARAMETERS FOR  Fk,Gk,Rk,Nk,Mk or Vk INTEGRAL: Example', &
         'F 0( 1s, 2s)', 'FORMAT:  (A1,I2,1X,4(A3,1X)) ' 
      READ (5, '(A1,I2,1X,4(A3,1X))') FUNC, K, EL1, EL2, EL3, EL4 
      IF (FUNC>='a' .AND. FUNC<='z') FUNC = CHAR(ICHAR(FUNC) + ICHAR('A') - &
         ICHAR('a')) 
      CALL EPTR (EL, EL1, I1, J2) 
      IF (J2 == 1) GO TO 20 
      CALL EPTR (EL, EL2, I2, J3) 
      IF (J3 == 1) GO TO 20 
      IF (EL3 /= ' ') THEN 
         CALL EPTR (EL, EL3, I3, J4) 
         IF (J4 == 1) GO TO 20 
      ENDIF 
      IF (EL4 /= ' ') THEN 
         CALL EPTR (EL, EL4, I4, J5) 
         IF (J5 == 1) GO TO 20 
      ENDIF 
      SELECT CASE (FUNC)  
      CASE ('F')  
         SI = FK(I1,I2,K,REL) 
      CASE ('G')  
         SI = GK(I1,I2,K,REL) 
      CASE ('R')  
         SI = RK(I1,I2,I3,I4,K,REL) 
      CASE ('N')  
         SI = SN(I1,I2,I2,I1,K) 
      CASE ('M')  
         SI = SN(I1,I2,I1,I2,K) 
      CASE ('V')  
         SI = VK(I1,I2,I2,I1,K) - VK(I2,I1,I1,I2,K) 
      CASE DEFAULT 
         WRITE (0, 41) 
   41    FORMAT(15X,'INTEGRAL UNKNOWN: RE-ENTER') 
         GO TO 20 
      END SELECT 
      IF (FUNC /= 'R') THEN 
         WRITE (3, 25) FUNC, K, EL1, EL2, SI, 219474.D0*SI 
         WRITE (0, 25) FUNC, K, EL1, EL2, SI, 219474.D0*SI 
   25    FORMAT(/,15X,'INTEGRAL  ',A1,I2,'(',A3,',',A3,') = ',1P,D14.7,' a.u.'/&
            ,40X,0P,F14.3,' cm-1'/) 
      ELSE 
         WRITE (3, 26) FUNC, K, EL1, EL2, EL3, EL4, SI, 219474.D0*SI 
         WRITE (0, 26) FUNC, K, EL1, EL2, EL3, EL4, SI, 219474.D0*SI 
   26    FORMAT(/,15X,'INTEGRAL  ',A1,I2,'(',2A3,',',2A3,') = ',1P,D14.7,&
            ' a.u.'/,46X,0P,F14.3,' cm-1'/) 
      ENDIF 
      GO TO 4 
 
!  ****  COMPUTE ELECTRON DENSITY AT THE NUCLEUS
 
   30 CONTINUE 
      WRITE (0, '(/5X,A/A,T22,A)') &
         'INPUT IDENTIFYING LABEL FOR ELECTRON: Example', '  1s', &
         'FORMAT(1X,A3) ' 
      READ (5, '(1X,A3)') EL1 
      CALL EPTR (EL, EL1, I, J6) 
      IF (J6 == 1) GO TO 30 
      LL = L(I) 
      IF (LL == 0) THEN 
         D = AZ(I)**2 
      ELSE 
         D = 0 
      ENDIF 
      WRITE (3, 32) EL1, D 
      WRITE (0, 32) EL1, D 
   32 FORMAT(/,15X,'DENSITY AT THE NUCLEUS FOR ',A3,' = ',1P,D14.7,' a.u.'/) 
      GO TO 4 
 
!  ****  COMPUTE SPIN-ORBIT PARAMETER
 
   40 CONTINUE 
      WRITE (0, '(/,5X,A/A,T22,A)') &
         'INPUT IDENTIFYING LABEL FOR ELECTRON: Example', '  2p', &
         'FORMAT(1X,A3) ' 
      READ (5, '(1X,A3)') EL1 
      CALL EPTR (EL, EL1, I, J7) 
      IF (J7 == 1) GO TO 40 
      ZETA = 0.D0 
      IF (L(I) /= 0) ZETA = BWZETA(I) 
      ZETACM = 219474*ZETA 
      WRITE (3, 43) EL1, ZETA, ZETACM 
      WRITE (0, 43) EL1, ZETA, ZETACM 
   43 FORMAT(/,15X,'SPIN-ORBIT PARAMETER FOR ',A3,' = ',1P,D14.7,' a.u.'/,46X,&
         0P,F14.3,' cm-1'/) 
      GO TO 4 
 
!  ****  COMPUTE TRANSITION INTEGRALS
 
   50 CONTINUE 
      WRITE (0, '(/5X,A/A,T22,A)') &
         'INPUT IDENTIFYING LABELS AND POWER OF R: Example', 'T 1( 2s, 2p)', &
         'FORMAT:  (A1,I2,2(1X,A3)) ' 
      READ (5, '(A1,I2,1X,A3,1X,A3)') FUNC, K, EL1, EL2 
      CALL EPTR (EL, EL1, I1, J8) 
      IF (J8 == 1) GO TO 50 
      CALL EPTR (EL, EL2, I2, J9) 
      IF (J9 == 1) GO TO 50 
      TI = QUADR(I1,I2,K) 
      WRITE (3, 52) FUNC, K, EL1, EL2, TI 
      WRITE (0, 52) FUNC, K, EL1, EL2, TI 
   52 FORMAT(/,15X,'INTEGRAL  ',A1,I2,'(',A3,',',A3,') = ',1P,D14.7,' a.u.'/) 
      GO TO 4 
   60 CONTINUE 
      RETURN  
      END SUBROUTINE MENU 
!
!     ------------------------------------------------------------------
!               M E T H O D
!     ------------------------------------------------------------------
!
!       Uses M1, M2, or M3 to solve the radial equation. If the input
!   data indicated METH(I) = 3, then this  solution  is  returned  to
!   DE.  Otherwise,  the routine searches for an acceptable  solution
!   which  is  both  positive  near  the  origin and has the required
!   number  of nodes.
!
      SUBROUTINE METHD1(I) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE TEST_C 
      USE COEFF_C 
      USE WAVE_C 
      USE RADIAL_C, ONLY: L, N 
      USE PARAM_C 
      USE LABEL_C, ONLY: EL 
      USE DE_C
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE solve_I 
      USE nodec_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: MN, NC, J 
      REAL(DOUBLE) :: DEL, EDP 
      LOGICAL :: V2, FIRST 
!-----------------------------------------------
      FIRST = .TRUE. 
      FAIL = .FALSE. 
      EM = D0 
      EU = ((Z - DMIN1(D5*S(I),D2*S(I)))/N(I))**2 
      FU = EU 
      MK = 0 
   17 CONTINUE 
      CALL SOLVE (I, FIRST, REL) 
!
!  *****  IF KK EQUALS 3, OMIT THE NODE CHECKING
!
      IF (KK /= 3) THEN 
!
!  *****  COUNT THE NUMBER OF NODES
!
         MN = M 
         NC = NODEC(MN) 
         IF (TRACE) WRITE (6, 99) EL(I), NC, MN, NJ, PDE(MN), ED, EU, EM, &
            DELTAE 
   99    FORMAT(2X,A3,' NC =',I3,' MN =',I3,' NJ =',I3,' PDE(MN) =',D10.2,&
            ' ED =',D10.2,' EU =',D10.2,' EM =',D10.2,' DELTAE =',D10.2) 
!
!  *****  IF NODE COUNT IS OFF BY NO MORE THAN 1 AND DELTAE IS STILL
!  *****  QUITE LARGE, APPLY THE DELTAE CORRECTION
!
         IF (IABS(NC - NODE)==1 .AND. DABS(DELTAE/ED)>0.02D0) GO TO 46 
!
!  *****  BRANCH ACCORDING TO WHETHER THE NODE COUNT IS TOO SMALL,
!  *****  JUST RIGHT, OR TOO LARGE
!
         IF (NC - NODE < 0) GO TO 8 
         IF (NC - NODE > 0) GO TO 10 
         V2 = DABS(DELTAE)<1.D-3 .OR. DABS(DELTAE)/ED<1.D-5 
         IF (PDE(MN)<D0 .AND. .NOT.V2) GO TO 46 
         IF (PDE(MN) <= D0) THEN 
            PDE(:NO) = -PDE(:NO) 
            PP = (-D2) - PP 
         ENDIF 
      ENDIF 
      AZZ = AZD*(D1 + PP) 
      E(I,I) = ED 
      RETURN  
!
!  *****  THE SOLUTION HAS TOO FEW NODES
!
    8 CONTINUE 
      IF (PDE(MN) <= D0) GO TO 11 
      DEL = D1 - ED/EU 
      EU = ED 
      IF (DEL < .05D0) THEN 
         FU = FU*((L(I)+1+NC)/FN)**2.5 
      ELSE 
         FU = ED*((L(I)+1+NC)/FN)**2.5 
      ENDIF 
      IF (FU < EM) FU = D5*(EU + EM) 
      IF (DABS(FU - ED) < 0.001D0) GO TO 27 
      ED = FU 
      GO TO 33 
!
!  *****  TRY A NEW VALUE OF ED WHICH MUST LIE WITHIN THE UPPER AND
!  *****  LOWER BOUND
!
   11 CONTINUE 
      EDP = ED 
      ED = ED*((L(I)+1+NC)/FN)**2.5 
      IF (ED >= EU) ED = D5*(EU + EDP) 
      IF (ED <= EM) ED = D5*(EM + EDP) 
   33 CONTINUE 
      MK = MK + 1 
      IF (EU <= EM) WRITE (6, 30) EM, EU, ED 
   30 FORMAT(6X,'WARNING: DIFFICULTY WITH NODE COUNTING PROCEDURE'/,6X,&
         'LOWER BOUND ON ED GREATER THAN UPPER BOUND'/,6X,'EL = ',F10.6,&
         '  EU = ',F10.6,'  ED = ',F10.6) 
      FIRST = .FALSE. 
      IF (MK>3*N(I) .OR. EU-EM<FN**(-3)) GO TO 27 
      GO TO 17 
!
!  *****  THE SOLUTION HAS TOO MANY NODES
!
   10 CONTINUE 
      IF (PDE(MN) < D0) GO TO 11 
      DEL = D1 - EM/ED 
      EM = ED 
      IF (DEL < 0.05D0) THEN 
         FM = FM*((L(I)+1+NC)/FN)**2.5 
      ELSE 
         FM = ED*((L(I)+1+NC)/FN)**2.5 
      ENDIF 
      IF (FM > EU) FM = D5*(EU + EM) 
      IF (DABS(FM - ED) < 0.001D0) GO TO 27 
      ED = FM 
      GO TO 33 
!
!  *****  ADJUST ENERGY TO LIE BETWEEN UPPER AND LOWER BOUND
!
   46 CONTINUE 
      ED = ED - DELTAE 
      IF (ED>=EM .AND. ED<=EU) GO TO 33 
      EDP = ED 
      IF (NC - NODE /= 0) ED = (ED + DELTAE)*((L(I)+1+NC)/FN)**2.5 
      IF (ED>=EM .AND. ED<=EU) GO TO 33 
      ED = EDP + DELTAE + DELTAE 
      IF (ED>=EM .AND. ED<=EU) GO TO 33 
      ED = ED - DELTAE 
      DELTAE = D5*DELTAE 
      GO TO 46 
!
!  *****  METHOD WAS UNABLE TO FIND AN ACCEPTABLE SOLUTION
!
   27 CONTINUE 
      WRITE (6, 28) KK, EL(I), NC, NJ, ED, EM, EU 
   28 FORMAT(10X,'METHOD',I2,' UNABLE TO SOLVE EQUATION FOR ELECTRON',A3,/,10X,&
         'NC = ',I3,3X,'NJ = ',I3,3X,'ED = ',F10.6,3X,'EL = ',F10.6,3X,'EU = ',&
         F10.6) 
      FAIL = .TRUE. 
      RETURN  
      END SUBROUTINE METHD1 
!
!     ------------------------------------------------------------------
!               N M R V S
!     ------------------------------------------------------------------
!
!       Given two starting values, PDE(1) and PDE(2), values of PDE(j),
!   j=3,4,...,NJ+1 are obtained by outward integration of
!               Y" = YR y + F
!   using the discretization  of  Eq.  (6-27 )  with  the  difference
!   correction.  With PDE(NJ) given, the tail procedure is applied to
!   PDE(j),j=NJ+1,  NJ+2,...,MM, where MM is determined automatically
!   and DELTA is the difference between  PDE(NJ+1)  for  outward  and
!   inward integration. (See Eq 6-32, 6-33, and 6-37 for further
!   details.)
!
!
      SUBROUTINE NMRVS(NJ, DELTA, MM, PDE, F) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE RADIAL_C 
      USE PARAM_C 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NJ 
      INTEGER , INTENT(OUT) :: MM 
      REAL(DOUBLE) , INTENT(OUT) :: DELTA 
      REAL(DOUBLE) , INTENT(INOUT) :: PDE(NOD)
      REAL(DOUBLE) , INTENT(IN) :: F(NOD)
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: M, I, K, J, II 
      REAL(DOUBLE), DIMENSION(150) :: A, D 
      REAL(DOUBLE) :: G3, Y1, Y2, G1, G2, Y3, RATIO, CON 
!-----------------------------------------------
      Y1 = PDE(1) 
      Y2 = PDE(2) 
      G1 = YR(1) 
      G2 = YR(2) 
      M = NJ + 1 
      DO I = 3, M 
         G3 = YR(I) 
         Y3 = (Y2 + Y2 - Y1 + (D10*G2*Y2 + G1*Y1) + F(I-1))/(D1 - G3) 
         PDE(I) = Y3 
         Y1 = Y2 
         Y2 = Y3 
         G1 = G2 
         G2 = G3 
      END DO 
      DELTA = Y3 
!
!  *****  APPLY THE TAIL PROCEDURE
!
      K = 1 
      PDE(M) = (-(D1 - G1)*Y1) + F(M) 
      A(1) = D1 - G3 
      D(1) = -(D2 + D10*G3) 
      RATIO = A(K)/D(K) 
!
!  *****  THE INTEGER 149 IN THE NEXT STATEMENT IS THE DIMENSION OF A
!  *****  MINUS 1
!
      IF (K>=150 - 1 .OR. M==ND) GO TO 23 
      K = K + 1 
      M = M + 1 
      G3 = YR(M) 
      A(K) = D1 - G3 
      D(K) = (-(D2 + D10*G3)) - A(K)*RATIO 
      PDE(M) = (-PDE(M-1)*RATIO) + F(M) 
      DO WHILE(DABS(PDE(M)) + DABS(PDE(M-1))>TOL .OR. K<9) 
         RATIO = A(K)/D(K) 
!
!  *****  THE INTEGER 149 IN THE NEXT STATEMENT IS THE DIMENSION OF A
!  *****  MINUS 1
!
         IF (K>=150 - 1 .OR. M==ND) GO TO 23 
         K = K + 1 
         M = M + 1 
         G3 = YR(M) 
         A(K) = D1 - G3 
         D(K) = (-(D2 + D10*G3)) - A(K)*RATIO 
         PDE(M) = (-PDE(M-1)*RATIO) + F(M) 
      END DO 
   20 CONTINUE 
      CON = DSQRT(EH)*DEXP((-DSQRT(DABS(G3/CH - 0.25)/RR(M))*(R(M+1)-R(M)))) 
      PDE(M) = PDE(M)/(D(K)+CON*(D1-YR(M+1))) 
      J = M + 1 
      I = J 
      IF (NO - J + 1 > 0) THEN 
         PDE(J:NO) = D0 
         I = NO + 1 
      ENDIF 
      DO J = 2, K 
         I = M - J + 1 
         II = K - J + 1 
         PDE(I) = (PDE(I)-A(II+1)*PDE(I+1))/D(II) 
      END DO 
!
!  *****  SET DELTA = DIFFERENCE OF THE TWO SOLUTIONS AT NJ+1
!  *****         MM = NUMBER OF POINTS IN THE RANGE OF THE SOLUTION
!
      DELTA = DELTA - PDE(I) 
      MM = M 
      RETURN  
   23 CONTINUE 
      WRITE (6, 24) 
   24 FORMAT(6X,'WARNING: FUNCTIONS TRUNCATED BY NMRVS IN TAIL REGION') 
      GO TO 20 
      END SUBROUTINE NMRVS 
!
!     ------------------------------------------------------------------
!               N O D E C
!     ------------------------------------------------------------------
!
!      Counts the number of nodes of the function PDE(j) in the range
!   j = 40,...,M-10.   The node counting procedure counts the local max
!   and min values.   Only nodes between sufficiently large max and
!   min values are counted.
!
!
      INTEGER FUNCTION NODEC (M) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE WAVE_C 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(INOUT) :: M 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: MM, J, NCC 
      REAL(DOUBLE) :: DM, SIGN, DIFF1, DIFF2 
!-----------------------------------------------
!
!   ***** FIND MAX|PDE(J)|
!
      MM = M - 10 
      DM = 0.D0 
      DO J = 40, MM 
         DM = DMAX1(DM,DABS(PDE(J))) 
      END DO 
!
!   *****  COUNT THE NUMBER OF LOCAL MAX OR MIN'S
!
      NCC = 0 
      SIGN = 0.D0 
      DIFF1 = PDE(40) - PDE(39) 
      DO J = 40, MM 
         DIFF2 = PDE(J+1) - PDE(J) 
         IF (DIFF2*DIFF1<=0.D0 .AND. DIFF1/=0.D0) THEN 
!
!   *****  A MAX OR MIN HAS BEEN FOUND.   TEST IF IT IS
!          SUFFICIENTLY LARGE
!
            IF (DABS(PDE(J))/DM >= 0.05D0) THEN 
!
!   ***** CHECK IF THIS IS THE FIRST SIGNIFICANT MAXIMUM
!
               IF (SIGN == 0.D0) THEN 
                  M = J 
               ELSE 
!
!   ***** IF NOT THE FIRST, TEST WHETHER A SIGN CHANGE HAS
!         OCCURRED SINCE THE LAST SIGNIFICANT MAX OR MIN
!
                  IF (PDE(J)*SIGN > 0.D0) GO TO 1002 
                  NCC = NCC + 1 
               ENDIF 
!
!   ***** RESET FOR THE NEXT NODE
!
               SIGN = PDE(J) 
            ENDIF 
         ENDIF 
 1002    CONTINUE 
         DIFF1 = DIFF2 
      END DO 
      NODEC = NCC 
      RETURN  
      END FUNCTION NODEC 


!
!     ------------------------------------------------------------------
!               O R T H O G
!     ------------------------------------------------------------------
!
!       This routine orthogonalizes the set of radial functions when an
!   orthogonality constraint applies.  A Gram-Schmidt type of  process
!   is used.
!
      SUBROUTINE ORTHOG 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE de_C
      USE TEST_C 
      USE PARAM_C 
      USE RADIAL_C 
      USE WAVE_C, ONLY: E 
      USE LABEL_C, ONLY: EL 
      USE COEFF_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE quadr_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: II, I, J, JJ 
      REAL(DOUBLE) :: PNN 
      LOGICAL :: CHANGE 
!-----------------------------------------------
!
      IF (NWF==1 .OR. IB>NWF) RETURN  
      WRITE (6, 26) 
   26 FORMAT(/) 
      II = MAX0(2,IB) 
      DO I = II, NWF 
         CHANGE = .FALSE. 
         AZZ = AZ(I) 
         DO J = 1, I - 1 
            IF (E(I,J) == D0) CYCLE  
!
!        ORTHOGONALITY CONDITION APPLIES
!
            C = QUADR(I,J,0) 
            IF (DABS(C) <= 1.D-10) CYCLE  
            WRITE (6, 63) EL(J), EL(I), C 
   63       FORMAT(6X,'<',A3,'|',A3,'>=',1P,D8.1) 
            M = MAX0(M,MAX(J)) 
            P(:M,I) = P(:M,I) - C*P(:M,J) 
            AZZ = AZZ - C*AZ(J) 
            CHANGE = .TRUE. 
         END DO 
         IF (.NOT.CHANGE) CYCLE  
         PNN = DSQRT(QUADR(I,I,0)) 
         IF (P(1,I) < D0) PNN = -PNN 
         P(:M,I) = P(:M,I)/PNN 
         AZZ = AZZ/PNN 
         M = NO 
   67    CONTINUE 
         IF (DABS(P(M,I)) < 1.D-15) THEN 
            P(M,I) = D0 
            M = M - 1 
            GO TO 67 
         ENDIF 
         MAX(I) = M 
         AZ(I) = AZZ 
      END DO 
      RETURN  
      END SUBROUTINE ORTHOG 
!
!     ------------------------------------------------------------------
!               O U T P U T
!     ------------------------------------------------------------------
!
!       The radial functions and orthogonality integrals are printed,
!   if PRINT is .TRUE.   The  functions  will  also  be  punched  (or
!   stored) on unit OUF, if OUF .NE. 0.
!
      SUBROUTINE OUTPUT(PRINT) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE INOUT_C 
      USE WAVE_C 
      USE PARAM_C 
      USE RADIAL_C, ONLY: R, R2, P, AZ, MAX 
      USE LABEL_C, ONLY: EL, ATOM, TERM 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      LOGICAL , INTENT(IN) :: PRINT 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ML, MU, I, MX, J, K, KK, JJ, IJ 
      REAL(DOUBLE), DIMENSION(8) :: OUT 
!-----------------------------------------------
      IF (PRINT) THEN 
         OPEN(Unit=50,FILE='plot.dat',STATUS='UNKNOWN',FORM='FORMATTED')
!
!  *****  PRINT RADIAL FUNCTIONS, 8 PER PAGE
!
         ML = IB 
    2    CONTINUE 
         MU = MIN0(ML + 7,NWF) 
         I = MU - ML + 1 
         MX = 0 
         DO J = ML, MU 
            MX = MAX0(MX,MAX(J)) 
         END DO 
         WRITE (3, 5) ATOM, TERM, (EL(J),J=ML,MU) 
         WRITE (50, *)
         WRITE (50, '(5X,"R",4X,8(3X,A3,4X))') (EL(J),J=ML,MU)
         WRITE (50, 11)  D0,  (D0, J=ML,MU) 
    5    FORMAT('1',9X,'WAVE FUNCTIONS FOR ',2A6,/,/,10X,'R',8(10X,A3)) 
         K = 0 
         KK = 0 
         DO J = 1, MX 
            OUT(:MU-ML+1) = P(J,ML:MU)*R2(J) 
            K = K + 1 
            IF (K > 10) THEN 
               K = 1 
               KK = KK + 1 
               IF (KK >= 5) THEN 
                  KK = 0 
                  WRITE (3, 23) 
   23             FORMAT('1'/,/) 
               ELSE 
                  WRITE (3, 8) 
    8             FORMAT(1X) 
               ENDIF 
            ENDIF 
            WRITE (3, 10) R(J), (OUT(JJ),JJ=1,I) 
            WRITE (50,11) R(J), (OUT(JJ),JJ=1,I) 
         END DO 
   10    FORMAT(F12.5,F12.6,7F11.6) 
   11    FORMAT(1P,9E10.2) 
         OUT(:MU-ML+1) = DPM(ML:MU) 
         WRITE (3, 16) (OUT(J),J=1,I) 
   16    FORMAT(3X,'MAX. DIFF.',F12.7,7F11.7) 
         ML = ML + 8 
         IF (ML <= NWF) GO TO 2 
         CLOSE (UNIT=50)
      ENDIF 
      
      IF (OUF /= 0) THEN 
!
!  *****  OUTPUT FUNCTIONS ON UNIT OUF FOR FUTURE INPUT
!
         DO I = 1, NWF 
            MX = MAX(I) 
            WRITE (OUF) ATOM, TERM, EL(I), MX, Z, E(I,I), EK(I), AZ(I), (P(J,I)&
               ,J=1,MX) 
         END DO 
      ENDIF 
!
      RETURN  
      END SUBROUTINE OUTPUT 
!
!     ------------------------------------------------------------------
!               P O T L
!     ------------------------------------------------------------------
!
!       Computes and stores the potential function
!                              2(k-1)
!              YR = SUM  a    Y      (j,j;r)
!                   j,k   ijk
!
      SUBROUTINE POTL(I, REL) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE RADIAL_C 
      USE PARAM_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE a_I 
      USE ykf_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I 
      LOGICAL  :: REL 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, K, JJ 
      REAL(DOUBLE) :: C 
!-----------------------------------------------
      YR(:NO) = D0 
      DO J = 1, NWF 
         DO K = 0, 2*MIN0(L(I),L(J)), 2 
            C = A(I,J,K) 
            IF (DABS(C) <= 1.D-8) CYCLE  
            CALL YKF (J, J, K, REL) 
            YR(:NO) = YR(:NO) + C*YK(:NO) 
         END DO 
      END DO 
      RETURN  
      END SUBROUTINE POTL 
!
!     ------------------------------------------------------------------
!               Q U A D
!     ------------------------------------------------------------------
!
!       Evaluates the integral of F(r)G(r) with respect to r , where
!   F(r) and G(r) have the same asymptotic properties as P (r).   The
!                                                         i
!   composite Simpson's rule is used.   The integrand is zero for r >
!   r  .
!    M
!
      REAL(KIND(0.0D0)) FUNCTION QUAD (I, M, F, G) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE PARAM_C 
      USE RADIAL_C 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I, M 
      REAL(DOUBLE) , INTENT(IN) :: F(NOD), G(NOD) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J 
      REAL(DOUBLE) :: D, QUAD2 
!-----------------------------------------------
      D = (D1 + D5*Z*R(1))/(H1*(2*L(I)+3)) 
      QUAD = RR(1)*F(1)*G(1)*(D - D5) 
      QUAD2 = D0 
      QUAD = QUAD + DOT_PRODUCT(RR(:M-1:2)*F(:M-1:2),G(:M-1:2)) 
      QUAD2 = QUAD2 + DOT_PRODUCT(RR(2:M:2)*F(2:M:2),G(2:M:2)) 
      QUAD = H1*(QUAD + D2*QUAD2) 
      RETURN  
      END FUNCTION QUAD 
!
!     ------------------------------------------------------------------
!                 Q U A D R
!     ------------------------------------------------------------------
!
!                                   kk
!       Evaluates the integral of  r   P (r) P (r) with respect to r
!                                       i     j
!
      REAL(KIND(0.0D0)) FUNCTION QUADR (I, J, KK) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE PARAM_C 
      USE RADIAL_C 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I, J, KK 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, LI, LJ, M, JJ, JP 
      REAL(DOUBLE) :: DEN, ZR, BI, BJ, ALPHA, BETA, D, DD 
!-----------------------------------------------
      K = KK + 2 
      LI = L(I) 
      LJ = L(J) 
      DEN = LI + LJ + 1 + K 
      ZR = Z*R(4) 
      BI = (P(4,I)/(AZ(I)*R2(4)*R(4)**LI)-D1+ZR/(LI+1))/ZR**2 
      BJ = (P(4,J)/(AZ(J)*R2(4)*R(4)**LJ)-D1+ZR/(LJ+1))/ZR**2 
      ALPHA = (D1/(LI + 1) + D1/(LJ + 1))/(DEN + D1) 
      ZR = Z*R(1) 
      BETA = (DEN + D1)*ALPHA**2 - D2*(BI + BJ + D1/((LI + 1)*(LJ + 1)))/(DEN&
          + D2) 
      D = P(1,I)*P(1,J)*R(1)**K*(((BETA*ZR + ALPHA)*ZR + D1)/(DEN*H1) + D5) 
      DD = D0 
      M = MIN0(MAX(I),MAX(J)) - 1 
      D = D + DOT_PRODUCT(P(3:M+1:2,I)*P(3:M+1:2,J),R(3:M+1:2)**K) 
      DD = DD + DOT_PRODUCT(P(2:M:2,I)*P(2:M:2,J),R(2:M:2)**K) 
      QUADR = H1*(D + D2*DD) 
      RETURN  
      END FUNCTION QUADR 
!
!     ------------------------------------------------------------------
!                 Q U A D S
!     ------------------------------------------------------------------
!
!                                       kk
!       Evaluates the integral of  (1/r)   YK(r) P (r) P (r)  with
!                                                 i     j
!   respect to r.
!
      REAL(KIND(0.0D0)) FUNCTION QUADS (I, J, KK) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE PARAM_C 
      USE RADIAL_C 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I, J, KK 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, MX, M 
      REAL(DOUBLE) :: DEN, CD, D, DD 
!-----------------------------------------------
      DEN = L(I) + L(J) + 3 
      K = 2 - KK 
      CD = D1 + Z*R(1)*(DEN - D1)/((DEN + D1)*((L(I)+1)*(L(J)+1))) 
      D = YK(1)*P(1,I)*P(1,J)*R(1)**K*(CD/(DEN*H1) + D5) 
      DD = D0 
      MX = MIN0(MAX(I),MAX(J)) - 1 
      DD = DD + DOT_PRODUCT(YK(2:MX:2)*P(2:MX:2,I)*P(2:MX:2,J),R(2:MX:2)**K) 
      D = D + DOT_PRODUCT(YK(3:MX+1:2)*P(3:MX+1:2,I)*P(3:MX+1:2,J),R(3:MX+1:2)&
         **K) 
      QUADS = H1*(D + D2*DD) 
      RETURN  
      END FUNCTION QUADS 
!
!     ------------------------------------------------------------------
!               R E F O R M
!     ------------------------------------------------------------------
!
!     Convert the free-format STR1 to the fixed 5(1X,A3,1X,I4,1X) format
!     for STR2
!
      SUBROUTINE REFORM(STR1, STR2) 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      CHARACTER  :: STR1*50 
      CHARACTER , INTENT(OUT) :: STR2*50 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I, IS, JS 
      CHARACTER(LEN=50) :: BLANK = '   '
!-----------------------------------------------
!
    1 CONTINUE 
      I = 0 
      STR2 = BLANK 
      IS = 0 
    2 CONTINUE 
      JS = INDEX(STR1(IS+1:),'(') 
      IF (JS /= 0) THEN 
         IF (JS > 5) GO TO 10 
         I = I + 5 
         STR2(I-JS+1:I) = STR1(IS+1:IS+JS) 
         IS = IS + JS 
         JS = INDEX(STR1(IS+1:),')') 
         IF (JS==0 .OR. JS>5) GO TO 10 
         I = I + 5 
         STR2(I-JS+1:I) = STR1(IS+1:IS+JS) 
         IS = IS + JS 
         GO TO 2 
      ENDIF 
      RETURN  
   10 CONTINUE 
      WRITE (0, *) ' Error in ', STR1, ': Re-enter' 
      READ (5, '(A)') STR1 
      GO TO 1 
      END SUBROUTINE REFORM 
!
!   --------------------------------------------------------------------
!               R E O R D
!   --------------------------------------------------------------------
!
!       Reorder the list of first appearance so that the functions to be
!   iterated appear last in the list.
!
      SUBROUTINE REORD(OF, ELC, NWF, IERR) 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE eptr_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: NWF 
      INTEGER , INTENT(OUT) :: IERR 
      CHARACTER(LEN=3)  :: ELC 
      CHARACTER(LEN=3), DIMENSION(:) :: OF
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J1, I, J 
!-----------------------------------------------
!
      IERR = 1 
      CALL EPTR (OF, ELC, I, J1) 
      IF (J1 == 1) GO TO 99 
      OF(I:NWF-1) = OF(I+1:NWF) 
      OF(NWF) = ELC 
      IERR = 0 
   99 CONTINUE 
      RETURN  
      END SUBROUTINE REORD 
!
!
!     ------------------------------------------------------------------
!                 R K
!     ------------------------------------------------------------------
!
!                   k
!       Evaluates  R (i, j; ii, jj)
!
      REAL(KIND(0.0D0)) FUNCTION RK (I, J, II, JJ, K, REL) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE ykf_I 
      USE quads_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I, J, II, JJ, K 
      LOGICAL  :: REL 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CALL YKF (I, II, K, REL) 
      RK = QUADS(J,JJ,1) 
      RETURN  
      END FUNCTION RK 
!
!     ------------------------------------------------------------------
!                R M E
!     ------------------------------------------------------------------
!
!
      REAL(KIND(0.0D0)) FUNCTION RME (L, LP, K) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE FACT_C 
!
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: L, LP, K 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I2G, IG, I1, I2, I3 
      REAL(DOUBLE) :: QUSQRT 
!-----------------------------------------------
!
!--- EVALUATES THE REDUCED MATRIX ELEMENT (L//C(K)//LP)  -  SEE FANO
!    AND RACAH, IRREDUCIBLE TENSORIAL SETS, CHAP. 14, P. 81
!
      IF (MIN0(L,LP) == 0) THEN 
         RME = 1.D0 
      ELSE IF (K == 0) THEN 
         RME = 2*L + 1 
         RME = DSQRT(RME) 
      ELSE 
         I2G = L + LP + K 
         IG = I2G/2 
         IF (I2G - 2*IG /= 0) THEN 
            RME = 0.D0 
         ELSE 
            I1 = IG - L 
            I2 = IG - LP 
            I3 = IG - K 
            QUSQRT = (2*L + 1)*(2*LP + 1) 
            RME = DSQRT(QUSQRT)*DEXP((GAM(2*I1+1)+GAM(2*I2+1)+GAM(2*I3+1)-GAM(&
               I2G+2))/2.D0+GAM(IG+1)-GAM(I1+1)-GAM(I2+1)-GAM(I3+1)) 
         ENDIF 
      ENDIF 
      RETURN  
      END FUNCTION RME 
!
!     ------------------------------------------------------------------
!               R O T A T E
!     ------------------------------------------------------------------
!
!        This routine analyses the energy expression to determine the
!   stationary condition with respect to rotation of orbials i and j.
!   If the condition is zero, the off-diagonal energy parameters may
!   be set to zero;  otherwise the orbials are rotated so as to satisfy
!   the stationay condition to first order in the perturbation.
!
!
      SUBROUTINE ROTATE(I, J) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE TEST_C 
      USE PARAM_C 
      USE RADIAL_C 
      USE WAVE_C, ONLY: E, SUM 
      USE LABEL_C, ONLY: EL 
      USE COEFF_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE hl_I 
      USE a_I 
      USE b_I 
      USE rk_I 
      USE fk_I 
      USE gk_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I 
      INTEGER  :: J 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: K, M, KK, JJ 
      REAL(DOUBLE) :: DG, C, FKII, FKIJ, GKIJ, CJ, FKJJ, EPS, DD, PI 
      REAL(DOUBLE) :: QI, QJ, G, DI, DJ, DII, DJJ, DIJ, DJI 
!-----------------------------------------------
      ALL = .TRUE. 
      G = D0 
      DG = D0 
      QI = SUM(I) 
      QJ = SUM(J) 
      IF (QI/=D2*(2*L(I)+1) .OR. QJ/=D2*(2*L(J)+1)) THEN 
         IF (DABS(QI - QJ) >= 1.D-14) THEN 
            C = D5*(QI - QJ) 
            G = G - C*HL(EL,I,J,REL) 
            DG = DG - C*(HL(EL,I,I,REL) - HL(EL,J,J,REL)) 
         ENDIF 
!
         DO K = 0, 2*L(I), 2 
            C = QI*(A(I,I,K) - A(I,J,K) - B(I,J,K)) 
            IF (DABS(C) >= 1.D-8) THEN 
               G = G + C*RK(I,I,I,J,K,REL) 
               FKII = FK(I,I,K,REL) 
               FKIJ = FK(I,J,K,REL) 
               GKIJ = GK(I,J,K,REL) 
               DG = DG + C*(FKII - FKIJ - D2*GKIJ) 
            ENDIF 
            CJ = QJ*(A(J,J,K) - A(J,I,K) - B(J,I,K)) 
            IF (DABS(CJ) < 1.D-8) CYCLE  
            FKJJ = FK(J,J,K,REL) 
            IF (DABS(C) < 1.D-8) THEN 
               FKIJ = FK(I,J,K,REL) 
               GKIJ = GK(I,J,K,REL) 
            ENDIF 
            G = G - CJ*RK(J,J,J,I,K,REL) 
            DG = DG + CJ*(FKJJ - FKIJ - D2*GKIJ) 
         END DO 
         DO M = 1, NWF 
            IF (M==I .OR. M==J) CYCLE  
            DO K = 0, 2*MIN0(L(I),L(M)), 2 
               C = A(I,M,K)*QI - A(J,M,K)*QJ 
               IF (DABS(C) >= 1.D-8) THEN 
                  G = G + C*RK(I,M,J,M,K,REL) 
                  DG = DG + C*(FK(I,M,K,REL) - FK(J,M,K,REL)) 
               ENDIF 
               KK = IABS(L(I)-L(M)) + K 
               C = B(I,M,KK)*QI - B(J,M,KK)*QJ 
               IF (DABS(C) < 1.D-8) CYCLE  
               G = G + C*RK(I,J,M,M,KK,REL) 
               DG = DG + C*(GK(I,M,KK,REL) - GK(J,M,KK,REL)) 
            END DO 
         END DO 
         IF (DABS(QI - QJ) + DABS(G) + DABS(DG) > 2.D-8) THEN 
            IF (DABS(G) + DABS(DG)>1.D-8 .OR. DABS(E(I,J))>2.D-5) THEN 
               EPS = G/DG 
               EPS = DSIGN(DMIN1(DABS(EPS),0.2D0),EPS) 
               DD = DSQRT(D1 + EPS*EPS) 
               DO JJ = 1, NO 
                  PI = (P(JJ,I)+EPS*P(JJ,J))/DD 
                  P(JJ,J) = (P(JJ,J)-EPS*P(JJ,I))/DD 
                  P(JJ,I) = PI 
               END DO 
            ELSE 
               EPS = D0 
            ENDIF 
            WRITE (6, 100) EL(I), EL(J), G, EL(I), EL(J), DG, EPS 
  100       FORMAT(10X,'C(',2A3,') =',F12.5,3X,'V(',2A3,') =',F12.5,3X,'EPS =',&
               F9.6) 
            RETURN  
!
!  *****  THE ENERGY IS STATIONARY WITH RESPECT TO ROTATIONS
!
         ENDIF 
      ENDIF 
      E(I,J) = 1.D-10 
      E(J,I) = 1.D-10 
      RETURN  
      END SUBROUTINE ROTATE 
!
!     ------------------------------------------------------------------
!                 S H I F T
!     ------------------------------------------------------------------
!
!       Computes the mass velocity  and one-body   Darwin term
!   corrections for the relativistic shift in the energy of the electron
!   including non-diagonal corrections
!
!
      REAL(KIND(0.0D0)) FUNCTION RLSHFT (I1, I2) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE PARAM_C 
      USE RADIAL_C 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I1 
      INTEGER , INTENT(IN) :: I2 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: LL, L2, L3, MX, J, I, KK, MM, K 
      REAL(DOUBLE) :: FL, C, ZZ, HH, A1, B1, B2, YY, A2, A, B, RELSH, RELSH2, &
         RLM, RLD 
!-----------------------------------------------
!
!  *****  FORM DD  -L(L+1)/RR|P(I)>
!
      FL = L(I1) 
      C = (FL + D5)**2 
      LL = L(I1) + 1 
      L2 = 2*L(I1) + 1 
      L3 = 2*L(I1) + 3 
      ZZ = Z*Z 
      HH = 180.D0*H*H 
      MX = MAX0(MAX(I1),MAX(I2)) 
      YK(2:MX) = -D1/RR(2:MX) 
!
!  *****  FORM THE INTEGRAND
!
      I = I1 
      A1 = D0 
      B1 = D0 
      DO KK = 1, 2 
         B2 = B1 
         YY = (P(3,I)+P(1,I)-D2*P(2,I))/(H*H) - C*P(2,I) 
         YK(2) = YY*YK(2) 
         YK(3) = YK(3)*(((-(P(5,I)+P(1,I)))+D16*(P(4,I)+P(2,I))-D30*P(3,I))/(&
            D12*H*H)-C*P(3,I)) 
         MM = MAX(I) - 3 
         DO K = 4, MM 
            YY = D2*(P(K+3,I)+P(K-3,I)) 
            YY = YY - 27.D0*(P(K+2,I)+P(K-2,I)) 
            YY = YY + 270.D0*(P(K+1,I)+P(K-1,I)) - 490.D0*P(K,I) 
            YY = YY/HH - C*P(K,I) 
            YK(K) = YY*YK(K) 
            IF (K /= 4) CYCLE  
            B1 = (YY/(D2*Z*P(4,I)*R(4))+D1)/R(4) 
         END DO 
         MM = MM + 1 
         YK(MM) = YK(MM)*(((-(P(MM+2,I)+P(MM-2,I)))+D16*(P(MM+1,I)+P(MM-1,I))-&
            D30*P(MM,I))/(D12*H*H)-C*P(MM,I)) 
         MM = MM + 1 
         YK(MM) = YK(MM)*((P(MM+1,I)+P(MM-1,I)-D2*P(MM,I))/(H*H)-C*P(MM,I)) 
         A2 = A1 
         A1 = (P(1,I)/(AZ(I)*R(1)**L(I)*R2(1))-D1+Z*R(1)/LL)/RR(1) 
         I = I2 
      END DO 
!
!  ***** DETERMINE CONTRIBUTION FROM NEAR THE NUCLEUS
!
      A = (Z/LL - L2*(B1 + B2)/D2)/LL 
      B = (L2*B1*B2 - D2*(A1 + A2) + (Z/LL**2)*(D2*Z*(D1 + D1/LL) - L2*(B1 + B2&
         )))/L3 
      RELSH = -P(4,I1)*P(4,I2)*(D1 + A*R(4)+B*RR(4))*D4*ZZ/L2 
      RELSH = RELSH/H1 - D5*YK(4) 
      RELSH2 = D0 
!
!  *****  INTEGRATE
!
      RELSH2 = RELSH2 + SUM(YK(5:MX:2)) 
      RELSH = RELSH + SUM(YK(4:MX-1:2)) 
      RELSH = (RELSH + D2*RELSH2)*H1 
 
!      IF ( L(I1) .EQ. 0 ) RELSH = RELSH + Z*AZ(I1)*AZ(I2)
!      RLSHFT = RELSH*D5*FINE
 
      RLM = RELSH*D5*FINE 
      RLD = 0.D0 
      IF (L(I1) == 0) RLD = Z*AZ(I1)*AZ(I2)*D5*FINE 
      RLSHFT = RLM + RLD 
 
!      write(3,'(a,4i4,a,3F10.5)') 'NL = ', N(i1),L(i1),N(i2),L(i2),
!     :               '  RLSHFT =', RLM,RLD,RLSHFT
 
      RETURN  
      END FUNCTION RLSHFT 
!
!     ------------------------------------------------------------------
!               S C A L E
!     ------------------------------------------------------------------
!
!       The current radial functions are scaled according to the
!   procedures of Sec. 7-2  .   Values of AZ and E(I,I), the starting
!   values and the diagonal energy parameters are also scaled.
!
!
      SUBROUTINE SCALE(ZZ) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE PARAM_C 
      USE RADIAL_C 
      USE WAVE_C, ONLY: E, S 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE quadr_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE) , INTENT(IN) :: ZZ 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, I ,K
      REAL(DOUBLE) :: RATIO, SR, SC, SS, F0, F1, F2, F3, PNORM, THETA 
      REAL(DOUBLE), DIMENSION(NOD) :: RS, PS 
!-----------------------------------------------
      RATIO = Z/ZZ 
      SR = DSQRT(RATIO) 
      R(:NO) = R(:NO)*RATIO 
      RR(:NO) = R(:NO)*R(:NO) 
      R2(:NO) = R2(:NO)*SR 
      DO I = 1, NWF 
         SC = (ZZ - S(I))/(Z - S(I)) 
         SS = SC*RATIO 
         E(I,I) = E(I,I)*SC**2 
         RS(:NO) = R(:NO)/SS 
         PS(:NO) = P(:NO,I)*SC 
         SC = (ZZ - D5*S(I))/(Z - D5*S(I)) 
         AZ(I) = AZ(I)*SC**(L(I)+1)*DSQRT(SC) 
         K = 3 
!
!  *****  INTERPOLATE THE (RS,PS) FUNCTIONS FOR VALUES OF P AT THE SET
!  *****  OF POINTS R
!
         DO J = 1, NO 
!
!  *****  SEARCH FOR THE NEAREST ENTRIES IN THE (RS,PS) TABLE
!
    5       CONTINUE 
            IF (K /= ND) THEN 
               IF (RS(K) > R(J)) GO TO 6 
               K = K + 1 
               GO TO 5 
!
!  *****  INTERPOLATE
!
    6          CONTINUE 
               THETA = DLOG(R(J)/RS(K-1))/H 
               F0 = PS(K-2) 
               F1 = PS(K-1) 
               F2 = PS(K) 
               F3 = PS(K+1) 
               P(J,I) = D5*(F1 + F2) + (THETA - D5)*(F2 - F1) + THETA*(THETA - &
                  D1)*(F0 - F1 - F2 + F3)/D4 
            ELSE 
               P(J,I) = D0 
            ENDIF 
         END DO 
         MAX(I) = NO 
!
!  ***** NORMALIZE THE INTERPOLATED FUNCTION
!
         PNORM = DSQRT(QUADR(I,I,0)) 
         P(:NO,I) = P(:NO,I)/PNORM 
      END DO 
      Z = ZZ 
      RETURN  
      END SUBROUTINE SCALE 


!
!     ------------------------------------------------------------------
!               S C F
!     -----------------------------------------------------------------
!
!       This routine controls the SCF procedure described in Chapter
!   7.  If certain input parameters are zero (or blank) they will  be
!   set to their default value.
!
!          Parameter       Default Value
!          --------        -------------
!          SCFTOL          1.D-7
!          I*              (NWF + 1 - IB)/4 + 3
!          NSCF            12
!
!   The self-consistency convergence criterion is
!
!          Z2 = SQRT( SCFTOL*(Z*NWF/2) )
!
!   It is increased by a factor two at the end of each iteration.
!
!
      SUBROUTINE SCF(ETOTAL, SCFTOL, EREL) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE TEST_C 
      USE LABEL_C 
      USE PARAM_C 
      USE WAVE_C, ONLY: SUM, DPM, IORD, IPR, NOD 
      USE INOUT_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE grange_I 
      USE de_I 
      USE orthog_I 
      USE help_I 
      USE energy_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE)  :: ETOTAL 
      REAL(DOUBLE) , INTENT(IN) :: SCFTOL 
      REAL(DOUBLE) , INTENT(OUT) :: EREL 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: ICYCLE, I, JP, J, JJ, II, NIT 
      REAL(DOUBLE) :: Z2, DP1, DP, CFGTOL, ENONR 
      LOGICAL :: LAST 
      CHARACTER :: ANS 
!-----------------------------------------------
!
!  *****  SET THE SCF CONVERGENCE PARAMETER TO AN OPTIMISTIC VALUE
!
      REL = .FALSE. 
      TOL = DSQRT(Z)*1.D-10 
      Z2 = SCFTOL*DSQRT(Z*NWF) 
      WRITE (6, 15) 
   15 FORMAT(/,/) 
      WRITE (6, 16) OMIT, SCFTOL, NO 
   16 FORMAT(10X,'WEAK ORTHOGONALIZATION DURING THE SCF CYCLE=',L4,/,10X,&
         'SCF CONVERGENCE TOLERANCE (FUNCTIONS)      =',1P,D9.2,/,10X,&
         'NUMBER OF POINTS IN THE MAXIMUM RANGE      =',I4) 
!
!  *****  SET ITERATION PARAMETERS
!
      IPR = 0 
      DP1 = D0 
      ETOTAL = D0 
      ICYCLE = 0 
      IF (IB <= NWF) THEN 
!
!  *****  PERFORM NSCF SELF-CONSISTENT FIELD ITERATIONS
!
         LAST = .FALSE. 
    9    CONTINUE 
         DO I = 1, NSCF 
            ICYCLE = ICYCLE + 1 
            WRITE (6, 7) ICYCLE, Z2 
    7       FORMAT(/,/,10X,'ITERATION NUMBER ',I2,/,10X,'----------------'/,/,&
               10X,'SCF CONVERGENCE CRITERIA (SCFTOL*SQRT(Z*NWF)) = ',1P,D9.1,/&
               ) 
            DP1 = D0 
            CALL GRANGE 
!
!  *****  SOLVE EACH DIFFERENTIAL EQUATION IN TURN
!
            WRITE (6, 14) 
   14       FORMAT(/,20X,' EL',9X,'ED',13X,'AZ',11X,'NORM',7X,'DPM') 
            DO JP = IB, NWF 
               J = IORD(JP) 
               CALL DE (J) 
               IF (FAIL) RETURN  
               DP = DPM(J)*DSQRT(SUM(J)) 
               IF (DP1 >= DP) CYCLE  
               DP1 = DP 
               JJ = J 
            END DO 
            IF (DP1 >= Z2) THEN 
!
!  *****  SOLVE IC DIFFERENTIAL EQUATIONS EACH TIME SELECTING THE
!  *****  ONE WITH THE LARGEST DPM
!
               DO II = 1, IC 
                  CALL DE (JJ) 
                  IF (FAIL) RETURN  
                  DP1 = D0 
                  DO JP = IB, NWF 
                     J = IORD(JP) 
                     DP = DSQRT(SUM(J))*DPM(J) 
                     IF (DP1 > DP) CYCLE  
                     JJ = J 
                     DP1 = DP 
                  END DO 
                  IF (DP1 >= Z2) CYCLE  
                  EXIT  
               END DO 
            ENDIF 
            CALL ORTHOG 
            IF (DP1<Z2 .AND. LAST) GO TO 17 
            IF (I /= NSCF) THEN 
!
!  *****  IF FUNCTIONS APPEAR TO HAVE CONVERGED,SOLVE EACH AGAIN, IN
!  *****  TURN, AND TEST AGAIN
!
               IF (DP1 <= Z2) LAST = .TRUE. 
            ENDIF 
!
!  *****  INCREASE THE CONVERGENCE CRITERION FOR SELF-CONSISTENCY
!
            Z2 = D2*Z2 
            WRITE (3, 8) EL(JJ), DP1 
            WRITE (0, 8) EL(JJ), DP1 
    8       FORMAT(/,6X,'LEAST SELF-CONSISTENT FUNCTION IS ',A3,&
               ' :WEIGHTED MAXIMUM CHANGE =',1P,D10.2) 
            CFGTOL = 1.4D0*CFGTOL 
         END DO 
         WRITE (0, 13) 
   13    FORMAT(10X,/,' SCF ITERATIONS HAVE CONVERGED TO THE ABOVE ACCURACY') 
         WRITE (3, 13) 
   20    CONTINUE 
         WRITE (0, '(A)') ' Do you wish to continue ? (Y/N/H) ' 
         READ (5, '(A)') ANS 
         IF (ANS=='H' .OR. ANS=='h') THEN 
            CALL HELP (10) 
            GO TO 20 
         ENDIF 
         IF (ANS=='Y' .OR. ANS=='y') THEN 
            WRITE (0, '(A,A)') ' Enter the additional iterations and new IC ', &
               'in FORMAT(I2, 1X, I2)' 
            READ (5, '(I2, 1X, I2)') NSCF, IC 
            GO TO 9 
         ENDIF 
         FAIL = .TRUE. 
!
!  *****  PERFORM RELATIVISTIC AND NON-RELATIVISTIC CALCULATIONS
!
      ENDIF 
   17 CONTINUE 
      CALL ENERGY (ENONR) 
      REL = .TRUE. 
      CALL ENERGY (ETOTAL) 
      REL = .FALSE. 
      EREL = ETOTAL - ENONR 
      NIT = NWF - IB + 1 
      WRITE (3, 105) NIT, DP1 
  105 FORMAT(/,/,10X,'NUMBER OF FUNCTIONS ITERATED          =',I6,/,10X,&
         'MAXIMUM WEIGHTED CHANGE IN FUNCTIONS  =',D10.2,/) 
      RETURN  
      END SUBROUTINE SCF 
!
!     ------------------------------------------------------------------
!               S E A R C H
!     ------------------------------------------------------------------
!
!       This routine searches for the NJ>70 such that YR(j) > 0 for all
!   j > NJ.
!
!
      SUBROUTINE SEARCH(NJ, I) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE RADIAL_C 
      USE PARAM_C 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(OUT) :: NJ 
      INTEGER , INTENT(IN) :: I 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: IA, IL, NK 
!-----------------------------------------------
      IA = 70 
      IL = NO 
    4 CONTINUE 
      IF (YR(IA) >= D0) THEN 
         IA = IA + 2 
         IF (IA < IL) GO TO 4 
         NJ = MAX0(70,MAX(I)-100) 
         RETURN  
      ENDIF 
      NK = (IA + IL)/2 
      IF (YR(NK) >= D0) THEN 
         IL = NK 
      ELSE 
         IA = NK 
      ENDIF 
      DO WHILE(IL - IA > 1) 
         NK = (IA + IL)/2 
         IF (YR(NK) >= D0) THEN 
            IL = NK 
         ELSE 
            IA = NK 
         ENDIF 
      END DO 
      NJ = IL - 7 
      RETURN  
      END SUBROUTINE SEARCH 

!
!     ------------------------------------------------------------------
!                 S N
!     ------------------------------------------------------------------
!
!                                      3              k
!       Evaluates the integral of (1/r)  P (r) P (r) Z (i, j; r)  with
!                                         i     j
!   respect to r.
!
      REAL(KIND(0.0D0)) FUNCTION SN (I, J, II, JJ, K) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE PARAM_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE zk_I 
      USE quads_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I, J, II, JJ, K 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CALL ZK (J, JJ, K) 
      SN = QUADS(I,II,3)*FINE 
      RETURN  
      END FUNCTION SN 
!
!     ------------------------------------------------------------------
!               S O L V E
!     ------------------------------------------------------------------
!
!       When FIRST is .TRUE., SOLVE computes the potential and exchange
!   function and initializes variables for the i'th radial  equation.
!   The vector P1 is the solution of the radial equation and P2 the
!   variation of the solution with respect to the energy parameter
!   E(I,I).
!
!
      SUBROUTINE SOLVE(I, FIRST, REL) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE de_C
      USE WAVE_C, ONLY: PDE, ED, E, AZD 
      USE RADIAL_C, ONLY: R, RR, R2, P, YK, YR, X, AZ, L, MAX, N 
      USE PARAM_C 
      USE LABEL_C, ONLY: EL 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE potl_I 
      USE xch_I 
      USE hnorm_I 
      USE hwf_I 
      USE search_I 
      USE nmrvs_I 
      USE quad_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I 
      LOGICAL , INTENT(IN) :: FIRST 
      LOGICAL  :: REL 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, JJ, MH, M1, M2 
      REAL(DOUBLE), DIMENSION(NOD) :: ZERO 
      REAL(DOUBLE) :: ZINF, FL, CD, X1, X2, X3, X4, X5, RL, F1, C11, FNUM, DEL1&
         , PN, B3, HW, DELH, PNORM, Y1, Y2, Y3, DELTA, DEL2, AA, A11, B11, DISC&
         , DE1, DE2 
!-----------------------------------------------
!
!  *****  IF FIRST IS 'TRUE', CALL POTL AND XCH AND SET UP ARRAYS
!
      IF (FIRST) THEN 
         CALL POTL (I, REL) 
         CALL XCH (I, 3) 
         ZINF = DMAX1(0.05D0,Z - YR(ND)) 
         FN = N(I) 
         FL = L(I) 
         V = YR(1)/R(1) 
         B4 = Z*(FL + D4/D3)/((FL + D1)*(FL + D2)) 
         CN = (D2*Z/FN)**(L(I)+1) 
         C = D4*FL + D6 
         CD = (FL + D5)**2 
         XY = X(1) 
         XP = X(2) 
         ED = E(I,I) 
         X1 = X(1) 
         X2 = X(2) 
         X3 = X(3) 
         X4 = X(4) 
         DO J = 3, ND 
            X5 = X(J+2) 
            X(J) = CH*((-X5) + 24.D0*(X4 + X2) + 194.D0*X3 - X1)/20.D0 
            X1 = X2 
            X2 = X3 
            X3 = X4 
            X4 = X5 
         END DO 
         X(NO-1) = CH*(X4 + D10*X3 + X2) 
         YK(:NO) = (-D2*(Z - YR(:NO))*R(:NO)) + CD 
         X1 = CH*P(1,I)*(YK(1)+ED*RR(1)) 
         X2 = CH*P(2,I)*(YK(2)+ED*RR(2)) 
         X3 = CH*P(3,I)*(YK(3)+ED*RR(3)) 
         X4 = CH*P(4,I)*(YK(4)+ED*RR(4)) 
         DO J = 3, ND 
            X5 = CH*P(J+2,I)*(YK(J+2)+ED*RR(J+2)) 
            X(J) = X(J) - (X5 - D4*(X2 + X4) + D6*X3 + X1)/20.D0 
            X1 = X2 
            X2 = X3 
            X3 = X4 
            X4 = X5 
         END DO 
         RL = L(I) + 2.5 
         X(2) = R(2)**RL*(X(5)/R(5)**RL-D3*(X(4)/R(4)**RL-X(3)/R(3)**RL)) 
!
!  *****  DETERMINE LOWER BOUND ON THE ENERGY PARAMETER
!
         IF (KK /= 3) GO TO 80 
         DO JJ = 15, ND 
            J = NO - JJ 
            IF (YK(J) < D0) GO TO 63 
         END DO 
         WRITE (6, 12) 
   12    FORMAT(10X,'POTENTIAL FUNCTION TOO SMALL - 2R*(Z-Y)<(L+.5)**2') 
!     STOP
         GO TO 80 
   63    CONTINUE 
         EM = -YK(J)/RR(J) 
         GO TO 81 
   80    CONTINUE 
         EM = (ZINF/(FN + D5))**2 
   81    CONTINUE 
         FM = EM 
!
!  *****  DETERMINE DIAGONAL ENERGY PARAMETER
!
         F1 = D0 
         C11 = D0 
         M = MIN0(MAX(I),NO-1) 
         DO J = 2, M 
            FNUM = P(J+1,I) - P(J,I) - P(J,I) + P(J-1,I) 
            FNUM = FNUM - CH*(YK(J+1)*P(J+1,I)+D10*YK(J)*P(J,I)+YK(J-1)*P(J-1,I&
               )) - X(J) 
            DEL1 = RR(J+1)*P(J+1,I) + D10*RR(J)*P(J,I) + RR(J-1)*P(J-1,I) 
            F1 = F1 + P(J,I)*FNUM 
            C11 = C11 + P(J,I)*DEL1 
         END DO 
         ED = F1/(C11*CH) 
         IF (ED <= EM) THEN 
!
!  *****  ERROR MESSAGE AND ENERGY ADJUSTMENT FOR AN ENERGY PARAMETER
!  *****  TOO SMALL FOR THE RANGE OF THE FUNCTION
!
            WRITE (6, 24) ED 
   24       FORMAT(10X,'ED = ',F10.6,'; ADJUSTED TO ALLOWED MINIMUM ENERGY') 
            ED = EM 
            IF (DABS(FM - E(I,I))<=1.D-6 .AND. KK/=3) THEN 
!
!  ***** RETURN HYDROGENIC FUNCTION
!
               PN = HNORM(N(I),L(I),ZINF) 
               DO J = 1, NO 
                  PDE(J) = PN*HWF(N(I),L(I),ZINF,R(J))/R2(J) 
               END DO 
               AZD = PN*(D2*ZINF/N(I))**(L(I)+1) 
               PP = D0 
               DELTAE = D0 
               ED = (ZINF/N(I))**2 
               WRITE (6, 66) EL(I), ZINF 
   66          FORMAT(/,/,10X,'RETURN HYDROGENIC FUNCTION FOR ',A3,&
                  ' WITH EFFECTIVE CHARGE ',F10.3) 
               RETURN  
!
!  *****  CHECK IF UPPER BOUND IS CORRECT
!
            ENDIF 
         ENDIF 
         IF (D10*ED >= EU) THEN 
            EU = D10*ED 
            FU = EU 
         ENDIF 
         AZD = AZ(I) 
      ENDIF 
      YR(:NO) = (YK(:NO)+ED*RR(:NO))*CH 
      ZERO(:NO) = D0 
!
!  *****  SEARCH FOR THE POINT AT WHICH YR BECOMES POSITIVE
!
      CALL SEARCH (NJ, I) 
!
!  *****  COMPUTE STARTING VALUES FROM SERIES EXPANSION
!
      B3 = (V + V + ED - (Z/FN)**2)/C 
      DO J = 1, 2 
         HW = HWF(N(I),L(I),Z,R(J))/CN 
         HQ(J) = AZD*(HW + R(J)**(L(I)+3)*B3*(D1-R(J)*B4))/R2(J) 
      END DO 
!
!  *****  OBTAIN HOMOGENEOUS SOLUTION
! 
      CALL NMRVS (NJ, DELH, MH, HQ, ZERO) 
      PDE(1) = HQ(1) + XY/C 
      PDE(2) = HQ(2) + XP/C 
!
!  *****  OBTAIN PARTICULAR SOLUTION
!
      CALL NMRVS (NJ, DEL1, M1, PDE, X) 
!
!  *****  DETERMINE THE ENERGY ADJUSTMENT REQUIRED FOR A SOLUTION WITH
!  *****  GIVEN A0
!
      M = MAX0(M1,MH) 
      PNORM = D0 
      PNORM = PNORM + DOT_PRODUCT(RR(:M)*HQ(:M),PDE(:M)) 
      Y1 = PDE(NJ-1) 
      Y2 = PDE(NJ) 
      Y3 = PDE(NJ+1) 
      DELTA = Y2 - Y1 + Y2 - Y3 + YR(NJ-1)*Y1 + D10*YR(NJ)*Y2 + YR(NJ+1)*Y3 + X&
         (NJ) 
      DELTAE = HQ(NJ)*DELTA/(H*H*PNORM) 
      PP = -DEL1/DELH 
!
!  *****  MATCH AT THE JOIN FOR A SOLUTION OF THE DIFFERENTIAL EQUATION
!
      PDE(:NO) = PDE(:NO) + PP*HQ(:NO) 
!
!  *****  IF  THE EQUATIONS APPEAR TO BE NEARLY
!  ****  SINGULAR, SOLVE THE VARIATIONAL EQUATIONS
!
      IF (KK /= 2) RETURN  
      X1 = P(1,I)*RR(1) 
      X2 = P(2,I)*RR(2) 
      P2(1) = X1/C 
      P2(2) = X2/C 
      DO J = 3, NO 
         X3 = P(J,I)*RR(J) 
         XX(J-1) = (D10*X2 + X1 + X3)*CH 
         X1 = X2 
         X2 = X3 
      END DO 
      CALL NMRVS (NJ, DEL2, M2, P2, XX) 
      AA = -DEL2/DELH 
      M = MAX0(M,M2) 
      P2(:NO) = P2(:NO) + AA*HQ(:NO) 
      A11 = QUAD(I,M,P2,P2) 
      B11 = QUAD(I,M,PDE,P2) 
      C11 = QUAD(I,M,PDE,PDE) - D1 
      DISC = B11*B11 - A11*C11 
      IF (DISC >= D0) THEN 
         DE1 = -(B11 + DSQRT(DISC))/A11 
         DE2 = C11/A11/DE1 
         IF (PDE(3) + DE1*P2(3) < D0) DE1 = DE2 
      ELSE 
         DE1 = C11/A11 
      ENDIF 
      PDE(:NO) = PDE(:NO) + DE1*P2(:NO) 
      PP = PP + DE1*AA 
      RETURN  
      END SUBROUTINE SOLVE 
!
!     ------------------------------------------------------------------
!               S U M M R Y
!     ------------------------------------------------------------------
!
!       The results of a calculation are summarized.   These include
!   the following for each electron:
!
!          E(NL)   - diagonal energy parameter
!          I(NL)   - -(1/2)<nl|L|nl>
!          KE      - I(NL) + Z <r>
!          REL     - Relativistic shift (mass-velocity, Darwin term,
!                    spin-spin contact term)
!          SIGMA   - screening parameter as defined by Eq. (6-  ).
!          AZ(NL)  - starting parameter, P(r)/r**(l+1) as r -> 0.
!          1/R**3  - expected value of <1/r**3>
!          1/R     - expected value of <1/r>
!          R       - expected mean radius
!          R**2    - expected value of <r**2>
!
!   These results are followed by:
!
!          KINETIC ENERGY (EK)
!          POTENTIAL ENERGY (EP) = ET - EN
!          RATIO                 = EP/EN
!          NON- RELATIVISTIC ENERGY (ET - EREL)
!          RELATIVISTIC SHIFT (EREL) FOR THE STATE
!          TOTAL ENERGY (ET)
!
      SUBROUTINE SUMMRY(ET, EREL) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE TEST_C 
      USE INOUT_C 
      USE LABEL_C 
      USE RADIAL_C, ONLY: AZ, L, N, NOD 
      USE WAVE_C, ONLY: EK, E, SUM, S 
      USE PARAM_C 
      USE COEFF_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE quadr_I 
      USE hl_I 
      USE rlshft_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(DOUBLE) , INTENT(IN) :: ET 
      REAL(DOUBLE) , INTENT(IN) :: EREL 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: I 
      REAL(DOUBLE) :: PI, EN, EKINP, RH, SC, RELS, RM3, RP2, RZ, ETN, EPOTL, &
         RATIO 
      REAL(DOUBLE), DIMENSION(NWFD) :: R1, RM1
!-----------------------------------------------
!
      PI = ACOS((-D1)) 
      WRITE (3, 9) ATOM, TERM 
    9 FORMAT(/,/,/,24X,'ATOM ',A6,3X,'TERM ',A6,/,/,2X,'nl',7X,'E(nl)',8X,&
         'I(nl)',7X,'KE(nl)',8X,'Rel(nl)',3X,'S(nl)',7X,'Az(nl)') 
      EN = D0 
      REL = .FALSE. 
!
!  *****  COMPUTE AND PRINT ONE-ELECTRON PARAMETERS
!
      DO I = 1, NWF 
         R1(I) = QUADR(I,I,1) 
         EK(I) = -D5*HL(EL,I,I,REL) 
         RM1(I) = QUADR(I,I,-1) 
         EKINP = EK(I) + Z*RM1(I) 
         EN = EN + SUM(I)*EKINP 
         RH = 3*N(I)*N(I) - L(I)*(L(I)+1) 
         SC = Z - D5*RH/R1(I) 
         S(I) = SC 
         RELS = RLSHFT(I,I) 
         WRITE (3, 15) EL(I), E(I,I), EK(I), EKINP, RELS, S(I), AZ(I) 
   15    FORMAT(1X,A3,F14.7,3F13.6,F8.3,F14.6) 
      END DO 
!
!  *****  Compute Moments
!
      WRITE (3, 8) 'Delta(R)' 
    8 FORMAT(/,/,2X,'nl',6X,A8,5X,'1/R**3',7X,'1/R',9X,'R',8X,'R**2') 
      DO I = 1, NWF 
         RM3 = 0 
         IF (L(I) /= 0) RM3 = QUADR(I,I,-3) 
         RP2 = QUADR(I,I,2) 
         RZ = 0. 
         IF (L(I) == 0) RZ = AZ(I)**2/(4.*PI) 
         WRITE (3, 16) EL(I), RZ, RM3, RM1(I), R1(I), RP2 
   16    FORMAT(1X,A3,F14.3,F13.4,F11.5,F10.5,F11.5) 
      END DO 
      ETN = ET - EREL 
      EPOTL = ETN - EN 
      RATIO = EPOTL/EN 
      WRITE (0, 26) ETN, EN, EREL, EPOTL, ET, RATIO 
      WRITE (3, 26) ETN, EN, EREL, EPOTL, ET, RATIO 
   26 FORMAT(/,/,5X,'TOTAL ENERGY (a.u.)'/,5X,'----- ------'/,10X,&
         ' Non-Relativistic   ',F15.8,T50,'Kinetic   ',F15.8,/,10X,&
         ' Relativistic Shift ',F15.8,T50,'Potential ',F15.8,/,10X,&
         ' Relativistic       ',F15.8,T50,'Ratio     ',F15.9) 
      RETURN  
      END SUBROUTINE SUMMRY 
!
!     ------------------------------------------------------------------
!                 V
!     ------------------------------------------------------------------
!
!                  k
!       Evaluates V (i,j) as defined by Blume and Watson (1962).
!
      REAL(KIND(0.0D0)) FUNCTION VK (I, J, II, JJ, K) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE PARAM_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE dyk_I 
      USE quads_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I 
      INTEGER  :: J 
      INTEGER  :: II 
      INTEGER  :: JJ 
      INTEGER  :: K 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      CALL DYK (I, II, K) 
      VK = QUADS(J,JJ,2)*FINE 
      RETURN  
      END FUNCTION VK 
!
!     ------------------------------------------------------------------
!               W A V E F N
!     ------------------------------------------------------------------
!
!       This routine initializes radial functions by the procedure
!   indicated by IND(I).
!
!         Value of IND(I)     Method
!         ---------------     ------
!             -1           Functions read from unit IU2
!              0           Screened hydrogenic functions with ZZ=Z-S(I)
!              1           Functions in memory left unchanged
!
!   The set of functions are then orthogonalized.
!
!
      SUBROUTINE WAVEFN 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE INOUT_C 
      USE PARAM_C 
      USE TEST_C 
      USE RADIAL_C 
      USE WAVE_C, ONLY: EK, E, S 
      USE LABEL_C, ONLY: EL 
      USE COEFF_C 
      USE ESTP_C, ONLY: ZZ, IND
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE eptr_I 
      USE hnorm_I 
      USE hwf_I 
      USE quadr_I 
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J1, I, J, IM, II, M, K, MT 
      REAL(DOUBLE) :: C, PNN, PN, Z2, FN, ZT, ETI, EKI, AZI  
      REAL(DOUBLE) :: PT(NOD)
      CHARACTER :: EL1*3, AT*6, TT*6 
      CHARACTER , DIMENSION(NWFD) :: ATM*6, TRM*6 
      CHARACTER :: TITLE*24 
!-----------------------------------------------
!
!  *****  GENERATE ARRAYS FOR R,R*R AND SQRT(R) WITH A CONSTANT MESH
!  *****  SIZE IN THE LOG(Z*R) VARIABLE
!
      DO I = 1, NO 
         R(I) = DEXP(RHO)/Z 
         RR(I) = R(I)*R(I) 
         R2(I) = DSQRT(R(I)) 
         RHO = RHO + H 
      END DO 
      RHO = RHO - NO*H 
!
!  ***** READ THE WAVEFUNCTIONS
!
      IF (IUF /= 0) THEN 
!   2    CONTINUE 
    2    READ (IUF, END=5) AT, TT, EL1, M, ZT, ETI, EKI, AZI, (PT(J),J=1,M) 
         CALL EPTR (EL, EL1, I, J1) 
         IF (J1 == 1) GO TO 2 
         IF (I>0 .AND. IND(I)==(-1)) THEN 
            ATM(I) = AT 
            TRM(I) = TT 
            MAX(I) = M 
            ZZ(I) = ZT 
            C = D1 
            IF (Z /= ZT) C = Z/ZT 
!
!  *****  SCALE RESULTS IF CARDS ARE FOR AN ATOM WITH A DIFFERENT Z
!
            E(I,I) = C*C*ETI 
            EK(I) = C*C*EKI 
            AZ(I) = AZI*C**(L(I)+1)*DSQRT(C) 
            P(:M,I) = C*PT(:M) 
!
!  *****  SET REMAINING VALUES IN THE RANGE = 0.
!
            IF (M /= NO) THEN 
               M = M + 1 
               P(M:NO,I) = D0 
            ENDIF 
            IND(I) = -2 
         ENDIF 
         GO TO 2 
!
!  *****  SET PARAMTERS FOR ELECTRONS AND INITIALIZE FUNCTIONS
!
      ENDIF 
    5 CONTINUE 
      DO I = 1, NWF 
         IF (IND(I) > 0) CYCLE  
         IF (IND(I) /= 0) THEN 
            IF (IND(I) == (-2)) GO TO 4 
            IND(I) = 0 
            WRITE (6, 27) EL(I) 
   27       FORMAT(8X,'WAVE FUNCTIONS NOT FOUND FOR ',A3) 
         ENDIF 
!
!  *****  DETERMINE ESTIMATES OF THE WAVE FUNCTIONS BY THE SCREENED
!  *****  HYDROGENIC APPROXIMATION
!
         IF (Z - S(I) > D0) THEN 
            PN = HNORM(N(I),L(I),Z-S(I)) 
         ELSE 
            WRITE (0, '(A,A)') ' Effective nuclear charge zero for ', EL(I) 
            STOP  
         ENDIF 
         DO J = 1, NO 
            P(J,I) = PN*HWF(N(I),L(I),Z-S(I),R(J))/R2(J) 
         END DO 
 
         M = NO 
   30    CONTINUE 
         IF (DABS(P(M,I)) > 1.D-15) GO TO 31 
         P(M,I) = D0 
         M = M - 1 
         GO TO 30 
   31    CONTINUE 
         MAX(I) = M 
!
!  ***** SET THE AZ(I) VALUE
!
         AZ(I) = PN*(D2*(Z - D5*S(I))/N(I))**(L(I)+1) 
!
!  *****  ORTHOGONALIZE TO INNER FUNCTIONS
!
    4    CONTINUE 
         IF (I == 1) CYCLE  
         IM = I - 1 
         DO II = 1, IM 
            IF (E(I,II) == D0) CYCLE  
            PN = QUADR(I,II,0) 
            IF (DABS(PN) <= 1.D-8) CYCLE  
            PNN = DSQRT(D1 - PN*PN) 
            IF (P(50,I) - PN*P(50,II) < D0) PNN = -PNN 
            M = MAX0(MAX(I),MAX(II)) 
            P(:M,I) = (P(:M,I)-PN*P(:M,II))/PNN 
         END DO 
      END DO 
      WRITE (3, 14) 
   14 FORMAT(/,/,/,8X,'INITIAL ESTIMATES '/,/,10X,'NL',4X,'SIGMA',6X,'E(NL)',4X&
         ,'AZ(NL)',4X,'FUNCTIONS'/,/) 
!
!  *****  COMPUTE ONE-ELECTRON ENERGY PARAMETERS IF THEY WERE NOT
!  *****  SPECIFIED ON INPUT.
!
      DO I = 1, NWF 
         K = IND(I) + 2 
         IF (IND(I) == (-2)) THEN 
            TITLE = ' SCALED '//ATM(I)//TRM(I) 
         ELSE IF (IND(I) == 0) THEN 
            TITLE = ' SCREENED HYDROGENIC' 
         ELSE 
            TITLE = ' UNCHANGED' 
         ENDIF 
         WRITE (3, 19) EL(I), S(I), E(I,I), AZ(I), TITLE 
         WRITE (6, 19) EL(I), S(I), E(I,I), AZ(I), TITLE 
   19    FORMAT(9X,A3,F9.2,F11.3,F10.3,3X,A24) 
      END DO 
      RETURN  
      END SUBROUTINE WAVEFN 
!
!     ------------------------------------------------------------------
!               X C H
!     ------------------------------------------------------------------
!
!       This routine computes functions associated with the exchange
!   function for the i'th radial equation,  including   contributions
!   from   the  interactions.    The   exact   form  of  the function
!   depends on the value of IOPT.
!
!          Value of IOPT      Function
!          -------------      --------
!              1             SQRT(r) X(r)
!              2             X(r)/SQRT(r)
!              3             r SQRT(r) ( X(r) + SUM e   P )
!                                                    ij  j
!
      SUBROUTINE XCH(I, IOPT) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE TEST_C 
      USE PARAM_C 
      USE WAVE_C 
      USE RADIAL_C, ONLY: R, RR, P, YK, X, L 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE b_I 
      USE ykf_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I 
      INTEGER , INTENT(IN) :: IOPT 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: J, K, JJ 
      REAL(DOUBLE) :: C 
!-----------------------------------------------
      X(:NO) = D0 
      DO J = 1, NWF 
         IF (J == I) CYCLE  
         DO K = IABS(L(I)-L(J)), L(I) + L(J), 2 
            C = B(I,J,K)*D2 
            IF (DABS(C) <= 1.D-10) CYCLE  
            CALL YKF (I, J, K, REL) 
            X(:NO) = X(:NO) + C*YK(:NO)*P(:NO,J) 
         END DO 
      END DO 
      GO TO (75,76,77) IOPT 
   76 CONTINUE 
      X(:NO) = X(:NO)/R(:NO) 
      GO TO 75 
   77 CONTINUE 
      X(:NO) = R(:NO)*X(:NO) 
      DO J = 1, NWF 
         C = E(I,J) 
         IF (DABS(C)<=1.D-10 .OR. J==I) CYCLE  
         X(:NO) = X(:NO) + C*P(:NO,J)*RR(:NO) 
      END DO 
!
!  *****  CHECK IF EXCHANGE IS ZERO: IF SO, METHOD 2 SHOULD BE USED.
!
   75 CONTINUE 
      IF (METH(I) == 2) RETURN  
      IF (DABS(X(1)) + DABS(X(2)) + DABS(X(3)) == D0) METH(I) = 2 
      RETURN  
      END SUBROUTINE XCH 
!
!     ------------------------------------------------------------------
!                 Y K F
!     ------------------------------------------------------------------
!
!               k
!       Stores Y (i, j; r) in the array YK
!
!
      SUBROUTINE YKF(I, J, K, REL) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE PARAM_C 
      USE RADIAL_C 
!-----------------------------------------------
!   I n t e r f a c e   B l o c k s
!-----------------------------------------------
      USE zk_I 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER  :: I, J, K 
      LOGICAL , INTENT(IN) :: REL 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: MX, M 
      REAL(DOUBLE) :: A, C, A2, H90, A3, AI, AN, A34, F1, F2, F3, F4, F5 
!-----------------------------------------------
      CALL ZK (I, J, K) 
      A = EH**(K + 1) 
      C = 2*K + 1 
      A2 = A*A 
      H90 = C*H3/D30 
      A3 = A2*A*H90 
      AI = H90/A 
      AN = 114.D0*A*H90 
      A34 = 34.D0*H90 
      MX = (MIN0(MAX(I),MAX(J))/2)*2 
      F1 = YK(MX)*EH**K 
      F2 = YK(MX) 
      F3 = YK(MX-1) 
      F4 = YK(MX-2) 
      DO M = MX - 2, 2, -1 
         F5 = YK(M-1) 
         YK(M) = YK(M+2)*A2 + (AN*F3 + A34*(F4 + A2*F2) - F5*AI - F1*A3) 
         F1 = F2 
         F2 = F3 
         F3 = F4 
         F4 = F5 
      END DO 
      YK(1) = YK(3)*A2 + C*H3*(F4 + D4*A*F3 + A2*F2) 
      IF (.NOT.REL) RETURN  
      C = C*FINE 
      YK(:MX) = YK(:MX) + C*P(:MX,I)*P(:MX,J) 
      RETURN  
      END SUBROUTINE YKF 
!
!     ------------------------------------------------------------------
!                 Z K
!     ------------------------------------------------------------------
!
!               k
!       Stores Z (i, j; r) in the array YK.
!
!
      SUBROUTINE ZK(I, J, K) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
      USE vast_kind_param, ONLY:  DOUBLE 
      USE PARAM_C 
      USE RADIAL_C 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER , INTENT(IN) :: I, J, K 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: MX, M, M1, M2 
      REAL(DOUBLE) :: DEN, FACT, A, A2, H90, A3, AI, AN, A34, F1, F2, F3, F4, &
         F5, C1, C2 
!-----------------------------------------------
      DEN = L(I) + L(J) + 3 + K 
      FACT = (D1/(L(I)+1)+D1/(L(J)+1))/(DEN + D1) 
      A = EH**K 
      A2 = A*A 
      H90 = H/90.D0 
      A3 = A2*A*H90 
      AI = H90/A 
      AN = 114.D0*A*H90 
      A34 = 34.D0*H90 
      F1 = RR(1)*P(1,I)*P(1,J) 
      F2 = RR(2)*P(2,I)*P(2,J) 
      F3 = RR(3)*P(3,I)*P(3,J) 
      F4 = RR(4)*P(4,I)*P(4,J) 
      YK(1) = F1*(D1 + Z*R(1)*FACT)/DEN 
      YK(2) = F2*(D1 + Z*R(2)*FACT)/DEN 
      YK(3) = YK(1)*A2 + H3*(F3 + D4*A*F2 + A2*F1) 
      MX = (MIN0(MAX(I),MAX(J))/2)*2 
      DO M = 5, MX 
         F5 = (RR(M)*P(M,I))*P(M,J) 
         YK(M-1) = YK(M-3)*A2 + (AN*F3 + A34*(F4 + A2*F2) - F5*AI - F1*A3) 
         F1 = F2 
         F2 = F3 
         F3 = F4 
         F4 = F5 
      END DO 
      M1 = MX - 1 
      IF (IABS(I - J) + IABS(K) == 0) THEN 
!
!  *****  FOR Y0(I,I) SET THE LIMIT TO 1 AND REMOVE OSCILLATIONS
!  *****  INTRODUCED BY THE USE OF SIMPSON'S RULE
!
         M2 = M1 - 1 
         C1 = D1 - YK(M1) 
         C2 = D1 - YK(M2) 
         YK(:M1:2) = YK(:M1:2) + C1 
         YK(2:M1+1:2) = YK(2:M1+1:2) + C2 
      ENDIF 
      DO M = M1 + 1, NO 
         YK(M) = A*YK(M-1) 
      END DO 
      RETURN  
      END SUBROUTINE ZK 
