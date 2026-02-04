C=======================================================================
C  READ_XYZIN
C
C  Read Merlino enriched XYZ file "xyzin".
C  Accept both atomic symbols and atomic numbers in the XYZ block.
C  Normalize symbols (case-insensitive) into correct form (Cl, Si, ...).
C=======================================================================

      SUBROUTINE READ_XYZIN(
     &   NAT, SYM, X, Y, Z,
     &   CHARGE, MULT, PG,
     &   A, B, C, KAPPA,
     &   SIGMA, T,
     &   QTRANS, QROT_ISO, QROT_AVG,
     &   IERR
     & )
      IMPLICIT NONE

C-----------------------------------------------------------------------
C Arguments
C-----------------------------------------------------------------------
      INTEGER NAT, CHARGE, MULT, SIGMA, IERR
      CHARACTER*2 SYM(*)
      CHARACTER*8 PG
      DOUBLE PRECISION X(*), Y(*), Z(*)
      DOUBLE PRECISION A, B, C, KAPPA
      DOUBLE PRECISION T, QTRANS, QROT_ISO, QROT_AVG

C-----------------------------------------------------------------------
C Locals
C-----------------------------------------------------------------------
      INTEGER I, ZNUM
      CHARACTER*160 LINE
      CHARACTER*32  KEY
      CHARACTER*16  ATOK
      DOUBLE PRECISION XX, YY, ZZ
      LOGICAL IN_BASIC, IN_ROT

C-----------------------------------------------------------------------
C Defaults
C-----------------------------------------------------------------------
      IERR     = 0
      CHARGE   = 0
      MULT     = 1
      PG       = 'C1'
      SIGMA    = 1
      T        = 298.15D0
      QTRANS   = 0.D0
      QROT_ISO = 0.D0
      QROT_AVG = 0.D0
      A        = 0.D0
      B        = 0.D0
      C        = 0.D0
      KAPPA    = 0.D0

      IN_BASIC = .FALSE.
      IN_ROT   = .FALSE.

C-----------------------------------------------------------------------
C Open input file (hard-coded)
C-----------------------------------------------------------------------
      OPEN(10, FILE='xyzin', STATUS='OLD', ERR=900)

C=======================================================================
C Read standard XYZ header
C   1st line: NAT
C   2nd line: comment (can be empty!)
C=======================================================================
      READ(10,*,ERR=910) NAT
      READ(10,'(A)',ERR=910) LINE

C=======================================================================
C Read coordinates
C   Accept:
C     C   x y z
C     cl  x y z
C     6   x y z
C     17  x y z
C=======================================================================
      DO 10 I = 1, NAT

C --- Read full line (robust also with blank/comment oddities)
   20    CONTINUE
         READ(10,'(A)',ERR=920) LINE
         IF (LINE .EQ. ' ') GO TO 20

C --- Parse first token + xyz
         ATOK = ' '
         XX   = 0.D0
         YY   = 0.D0
         ZZ   = 0.D0
         READ(LINE,*,ERR=920) ATOK, XX, YY, ZZ

C --- Try atomic number first
         ZNUM = 0
         CALL TRY_READ_INT(ATOK, ZNUM)

         IF (ZNUM .GT. 0) THEN
            CALL Z2SYMBOL(ZNUM, SYM(I), IERR)
            IF (IERR .NE. 0) THEN
               IERR = 3
               CLOSE(10)
               RETURN
            ENDIF
         ELSE
            CALL NORM_SYMBOL(ATOK, SYM(I), IERR)
            IF (IERR .NE. 0) THEN
               IERR = 3
               CLOSE(10)
               RETURN
            ENDIF
         ENDIF

         X(I) = XX
         Y(I) = YY
         Z(I) = ZZ

   10 CONTINUE

C=======================================================================
C Read remaining file line by line (optional sections)
C=======================================================================
 100  CONTINUE
      READ(10,'(A)',END=800) LINE
      CALL UPCASE(LINE)

C --- Skip empty lines
      IF (LINE .EQ. ' ') GO TO 100

C --- Section headers
      IF (LINE(1:1) .EQ. '#') THEN
         IN_BASIC = .FALSE.
         IN_ROT   = .FALSE.

         IF (LINE(1:6)  .EQ. '#BASIC')      IN_BASIC = .TRUE.
         IF (LINE(1:11) .EQ. '#ROTATIONAL') IN_ROT   = .TRUE.

         GO TO 100
      ENDIF

C --- BASIC section
      IF (IN_BASIC) THEN
         READ(LINE,*) KEY

         IF (KEY .EQ. 'CHARGE') THEN
            READ(LINE,*) KEY, CHARGE
         ELSE IF (KEY .EQ. 'SPIN_MULTIPLICITY') THEN
            READ(LINE,*) KEY, MULT
         ELSE IF (KEY .EQ. 'POINT_GROUP') THEN
            READ(LINE,*) KEY, PG
         ENDIF

         GO TO 100
      ENDIF

C --- ROTATIONAL section
      IF (IN_ROT) THEN
         READ(LINE,*) KEY

         IF (KEY .EQ. 'A_MHZ') THEN
            READ(LINE,*) KEY, A
         ELSE IF (KEY .EQ. 'B_MHZ') THEN
            READ(LINE,*) KEY, B
         ELSE IF (KEY .EQ. 'C_MHZ') THEN
            READ(LINE,*) KEY, C
         ELSE IF (KEY .EQ. 'KAPPA') THEN
            READ(LINE,*) KEY, KAPPA
         ELSE IF (KEY .EQ. 'SIGMA') THEN
            READ(LINE,*) KEY, SIGMA
         ELSE IF (KEY .EQ. 'T_K') THEN
            READ(LINE,*) KEY, T
         ELSE IF (KEY .EQ. 'Q_TRANS') THEN
            READ(LINE,*) KEY, QTRANS
         ELSE IF (KEY .EQ. 'Q_ROT_ISO') THEN
            READ(LINE,*) KEY, QROT_ISO
         ELSE IF (KEY .EQ. 'Q_ROT_AVG') THEN
            READ(LINE,*) KEY, QROT_AVG
         ENDIF

         GO TO 100
      ENDIF

      GO TO 100

C=======================================================================
C Normal end
C=======================================================================
 800  CONTINUE
      CLOSE(10)
      RETURN

C=======================================================================
C Errors
C=======================================================================
 900  CONTINUE
      IERR = 1
      RETURN

 910  CONTINUE
      IERR = 2
      CLOSE(10)
      RETURN

 920  CONTINUE
      IERR = 3
      CLOSE(10)
      RETURN

      END


C=======================================================================
C  TRY_READ_INT
C  Try to read integer from a token (ATOK).
C  If fails, returns 0.
C=======================================================================
      SUBROUTINE TRY_READ_INT(ATOK, IVAL)
      IMPLICIT NONE
      CHARACTER*(*) ATOK
      INTEGER IVAL

      IVAL = 0
      READ(ATOK,*,ERR=10) IVAL
      RETURN
   10 CONTINUE
      IVAL = 0
      RETURN
      END


C=======================================================================
C  NORM_SYMBOL
C  Normalize element symbol capitalization:
C    "cl" "CL" "cL" -> "Cl"
C    "c"  "C"       -> "C "
C  Output is CHARACTER*2.
C=======================================================================
      SUBROUTINE NORM_SYMBOL(ATOK, SYM, IERR)
      IMPLICIT NONE
      CHARACTER*(*) ATOK
      CHARACTER*2 SYM
      INTEGER IERR
      CHARACTER*2 TMP

      IERR = 0
      TMP  = '  '
      SYM  = '  '

      CALL TRIM2(ATOK, TMP)

      IF (TMP(1:1) .EQ. ' ') THEN
         IERR = 1
         RETURN
      ENDIF

C --- Normalize case:
C     first letter -> uppercase
C     second letter -> lowercase
      CALL TOUP1(TMP(1:1))
      CALL TOLO1(TMP(2:2))

      CALL CHECK_SYMBOL(TMP, IERR)
      IF (IERR .NE. 0) THEN
         SYM = '  '
         RETURN
      ENDIF

      SYM = TMP
      RETURN
      END


C=======================================================================
C  CHECK_SYMBOL
C  Validate that a CHARACTER*2 symbol is a real element (1..118).
C=======================================================================
      SUBROUTINE CHECK_SYMBOL(SYM, IERR)
      IMPLICIT NONE
      CHARACTER*2 SYM
      INTEGER IERR, ZNUM
      CHARACTER*2 TMP

      IERR = 0
      CALL SYMBOL2Z(SYM, ZNUM, IERR)
      IF (IERR .NE. 0) THEN
         IERR = 1
      ELSE
         CALL Z2SYMBOL(ZNUM, TMP, IERR)
         IF (IERR .NE. 0) IERR = 1
      ENDIF
      RETURN
      END


C=======================================================================
C  Z2SYMBOL
C  Convert atomic number (1..118) to element symbol (CHAR*2).
C=======================================================================
      SUBROUTINE Z2SYMBOL(ZNUM, SYM, IERR)
      IMPLICIT NONE
      INTEGER ZNUM, IERR
      CHARACTER*2 SYM

      IERR = 0
      SYM  = '  '

      IF (ZNUM .EQ.   1) SYM = 'H '
      IF (ZNUM .EQ.   2) SYM = 'He'
      IF (ZNUM .EQ.   3) SYM = 'Li'
      IF (ZNUM .EQ.   4) SYM = 'Be'
      IF (ZNUM .EQ.   5) SYM = 'B '
      IF (ZNUM .EQ.   6) SYM = 'C '
      IF (ZNUM .EQ.   7) SYM = 'N '
      IF (ZNUM .EQ.   8) SYM = 'O '
      IF (ZNUM .EQ.   9) SYM = 'F '
      IF (ZNUM .EQ.  10) SYM = 'Ne'
      IF (ZNUM .EQ.  11) SYM = 'Na'
      IF (ZNUM .EQ.  12) SYM = 'Mg'
      IF (ZNUM .EQ.  13) SYM = 'Al'
      IF (ZNUM .EQ.  14) SYM = 'Si'
      IF (ZNUM .EQ.  15) SYM = 'P '
      IF (ZNUM .EQ.  16) SYM = 'S '
      IF (ZNUM .EQ.  17) SYM = 'Cl'
      IF (ZNUM .EQ.  18) SYM = 'Ar'
      IF (ZNUM .EQ.  19) SYM = 'K '
      IF (ZNUM .EQ.  20) SYM = 'Ca'
      IF (ZNUM .EQ.  21) SYM = 'Sc'
      IF (ZNUM .EQ.  22) SYM = 'Ti'
      IF (ZNUM .EQ.  23) SYM = 'V '
      IF (ZNUM .EQ.  24) SYM = 'Cr'
      IF (ZNUM .EQ.  25) SYM = 'Mn'
      IF (ZNUM .EQ.  26) SYM = 'Fe'
      IF (ZNUM .EQ.  27) SYM = 'Co'
      IF (ZNUM .EQ.  28) SYM = 'Ni'
      IF (ZNUM .EQ.  29) SYM = 'Cu'
      IF (ZNUM .EQ.  30) SYM = 'Zn'
      IF (ZNUM .EQ.  31) SYM = 'Ga'
      IF (ZNUM .EQ.  32) SYM = 'Ge'
      IF (ZNUM .EQ.  33) SYM = 'As'
      IF (ZNUM .EQ.  34) SYM = 'Se'
      IF (ZNUM .EQ.  35) SYM = 'Br'
      IF (ZNUM .EQ.  36) SYM = 'Kr'
      IF (ZNUM .EQ.  37) SYM = 'Rb'
      IF (ZNUM .EQ.  38) SYM = 'Sr'
      IF (ZNUM .EQ.  39) SYM = 'Y '
      IF (ZNUM .EQ.  40) SYM = 'Zr'
      IF (ZNUM .EQ.  41) SYM = 'Nb'
      IF (ZNUM .EQ.  42) SYM = 'Mo'
      IF (ZNUM .EQ.  43) SYM = 'Tc'
      IF (ZNUM .EQ.  44) SYM = 'Ru'
      IF (ZNUM .EQ.  45) SYM = 'Rh'
      IF (ZNUM .EQ.  46) SYM = 'Pd'
      IF (ZNUM .EQ.  47) SYM = 'Ag'
      IF (ZNUM .EQ.  48) SYM = 'Cd'
      IF (ZNUM .EQ.  49) SYM = 'In'
      IF (ZNUM .EQ.  50) SYM = 'Sn'
      IF (ZNUM .EQ.  51) SYM = 'Sb'
      IF (ZNUM .EQ.  52) SYM = 'Te'
      IF (ZNUM .EQ.  53) SYM = 'I '
      IF (ZNUM .EQ.  54) SYM = 'Xe'
      IF (ZNUM .EQ.  55) SYM = 'Cs'
      IF (ZNUM .EQ.  56) SYM = 'Ba'
      IF (ZNUM .EQ.  57) SYM = 'La'
      IF (ZNUM .EQ.  58) SYM = 'Ce'
      IF (ZNUM .EQ.  59) SYM = 'Pr'
      IF (ZNUM .EQ.  60) SYM = 'Nd'
      IF (ZNUM .EQ.  61) SYM = 'Pm'
      IF (ZNUM .EQ.  62) SYM = 'Sm'
      IF (ZNUM .EQ.  63) SYM = 'Eu'
      IF (ZNUM .EQ.  64) SYM = 'Gd'
      IF (ZNUM .EQ.  65) SYM = 'Tb'
      IF (ZNUM .EQ.  66) SYM = 'Dy'
      IF (ZNUM .EQ.  67) SYM = 'Ho'
      IF (ZNUM .EQ.  68) SYM = 'Er'
      IF (ZNUM .EQ.  69) SYM = 'Tm'
      IF (ZNUM .EQ.  70) SYM = 'Yb'
      IF (ZNUM .EQ.  71) SYM = 'Lu'
      IF (ZNUM .EQ.  72) SYM = 'Hf'
      IF (ZNUM .EQ.  73) SYM = 'Ta'
      IF (ZNUM .EQ.  74) SYM = 'W '
      IF (ZNUM .EQ.  75) SYM = 'Re'
      IF (ZNUM .EQ.  76) SYM = 'Os'
      IF (ZNUM .EQ.  77) SYM = 'Ir'
      IF (ZNUM .EQ.  78) SYM = 'Pt'
      IF (ZNUM .EQ.  79) SYM = 'Au'
      IF (ZNUM .EQ.  80) SYM = 'Hg'
      IF (ZNUM .EQ.  81) SYM = 'Tl'
      IF (ZNUM .EQ.  82) SYM = 'Pb'
      IF (ZNUM .EQ.  83) SYM = 'Bi'
      IF (ZNUM .EQ.  84) SYM = 'Po'
      IF (ZNUM .EQ.  85) SYM = 'At'
      IF (ZNUM .EQ.  86) SYM = 'Rn'
      IF (ZNUM .EQ.  87) SYM = 'Fr'
      IF (ZNUM .EQ.  88) SYM = 'Ra'
      IF (ZNUM .EQ.  89) SYM = 'Ac'
      IF (ZNUM .EQ.  90) SYM = 'Th'
      IF (ZNUM .EQ.  91) SYM = 'Pa'
      IF (ZNUM .EQ.  92) SYM = 'U '
      IF (ZNUM .EQ.  93) SYM = 'Np'
      IF (ZNUM .EQ.  94) SYM = 'Pu'
      IF (ZNUM .EQ.  95) SYM = 'Am'
      IF (ZNUM .EQ.  96) SYM = 'Cm'
      IF (ZNUM .EQ.  97) SYM = 'Bk'
      IF (ZNUM .EQ.  98) SYM = 'Cf'
      IF (ZNUM .EQ.  99) SYM = 'Es'
      IF (ZNUM .EQ. 100) SYM = 'Fm'
      IF (ZNUM .EQ. 101) SYM = 'Md'
      IF (ZNUM .EQ. 102) SYM = 'No'
      IF (ZNUM .EQ. 103) SYM = 'Lr'
      IF (ZNUM .EQ. 104) SYM = 'Rf'
      IF (ZNUM .EQ. 105) SYM = 'Db'
      IF (ZNUM .EQ. 106) SYM = 'Sg'
      IF (ZNUM .EQ. 107) SYM = 'Bh'
      IF (ZNUM .EQ. 108) SYM = 'Hs'
      IF (ZNUM .EQ. 109) SYM = 'Mt'
      IF (ZNUM .EQ. 110) SYM = 'Ds'
      IF (ZNUM .EQ. 111) SYM = 'Rg'
      IF (ZNUM .EQ. 112) SYM = 'Cn'
      IF (ZNUM .EQ. 113) SYM = 'Nh'
      IF (ZNUM .EQ. 114) SYM = 'Fl'
      IF (ZNUM .EQ. 115) SYM = 'Mc'
      IF (ZNUM .EQ. 116) SYM = 'Lv'
      IF (ZNUM .EQ. 117) SYM = 'Ts'
      IF (ZNUM .EQ. 118) SYM = 'Og'

      IF (SYM .EQ. '  ') IERR = 1
      RETURN
      END


C=======================================================================
C  SYMBOL2Z
C  Convert element symbol (CHAR*2) to atomic number (1..118)
C=======================================================================
      SUBROUTINE SYMBOL2Z(SYM, ZNUM, IERR)
      IMPLICIT NONE
      CHARACTER*2 SYM
      INTEGER ZNUM, IERR

      IERR = 0
      ZNUM = 0

      IF (SYM .EQ. 'H ')  ZNUM =   1
      IF (SYM .EQ. 'He')  ZNUM =   2
      IF (SYM .EQ. 'Li')  ZNUM =   3
      IF (SYM .EQ. 'Be')  ZNUM =   4
      IF (SYM .EQ. 'B ')  ZNUM =   5
      IF (SYM .EQ. 'C ')  ZNUM =   6
      IF (SYM .EQ. 'N ')  ZNUM =   7
      IF (SYM .EQ. 'O ')  ZNUM =   8
      IF (SYM .EQ. 'F ')  ZNUM =   9
      IF (SYM .EQ. 'Ne')  ZNUM =  10
      IF (SYM .EQ. 'Na')  ZNUM =  11
      IF (SYM .EQ. 'Mg')  ZNUM =  12
      IF (SYM .EQ. 'Al')  ZNUM =  13
      IF (SYM .EQ. 'Si')  ZNUM =  14
      IF (SYM .EQ. 'P ')  ZNUM =  15
      IF (SYM .EQ. 'S ')  ZNUM =  16
      IF (SYM .EQ. 'Cl')  ZNUM =  17
      IF (SYM .EQ. 'Ar')  ZNUM =  18
      IF (SYM .EQ. 'K ')  ZNUM =  19
      IF (SYM .EQ. 'Ca')  ZNUM =  20
      IF (SYM .EQ. 'Sc')  ZNUM =  21
      IF (SYM .EQ. 'Ti')  ZNUM =  22
      IF (SYM .EQ. 'V ')  ZNUM =  23
      IF (SYM .EQ. 'Cr')  ZNUM =  24
      IF (SYM .EQ. 'Mn')  ZNUM =  25
      IF (SYM .EQ. 'Fe')  ZNUM =  26
      IF (SYM .EQ. 'Co')  ZNUM =  27
      IF (SYM .EQ. 'Ni')  ZNUM =  28
      IF (SYM .EQ. 'Cu')  ZNUM =  29
      IF (SYM .EQ. 'Zn')  ZNUM =  30
      IF (SYM .EQ. 'Ga')  ZNUM =  31
      IF (SYM .EQ. 'Ge')  ZNUM =  32
      IF (SYM .EQ. 'As')  ZNUM =  33
      IF (SYM .EQ. 'Se')  ZNUM =  34
      IF (SYM .EQ. 'Br')  ZNUM =  35
      IF (SYM .EQ. 'Kr')  ZNUM =  36
      IF (SYM .EQ. 'Rb')  ZNUM =  37
      IF (SYM .EQ. 'Sr')  ZNUM =  38
      IF (SYM .EQ. 'Y ')  ZNUM =  39
      IF (SYM .EQ. 'Zr')  ZNUM =  40
      IF (SYM .EQ. 'Nb')  ZNUM =  41
      IF (SYM .EQ. 'Mo')  ZNUM =  42
      IF (SYM .EQ. 'Tc')  ZNUM =  43
      IF (SYM .EQ. 'Ru')  ZNUM =  44
      IF (SYM .EQ. 'Rh')  ZNUM =  45
      IF (SYM .EQ. 'Pd')  ZNUM =  46
      IF (SYM .EQ. 'Ag')  ZNUM =  47
      IF (SYM .EQ. 'Cd')  ZNUM =  48
      IF (SYM .EQ. 'In')  ZNUM =  49
      IF (SYM .EQ. 'Sn')  ZNUM =  50
      IF (SYM .EQ. 'Sb')  ZNUM =  51
      IF (SYM .EQ. 'Te')  ZNUM =  52
      IF (SYM .EQ. 'I ')  ZNUM =  53
      IF (SYM .EQ. 'Xe')  ZNUM =  54
      IF (SYM .EQ. 'Cs')  ZNUM =  55
      IF (SYM .EQ. 'Ba')  ZNUM =  56
      IF (SYM .EQ. 'La')  ZNUM =  57
      IF (SYM .EQ. 'Ce')  ZNUM =  58
      IF (SYM .EQ. 'Pr')  ZNUM =  59
      IF (SYM .EQ. 'Nd')  ZNUM =  60
      IF (SYM .EQ. 'Pm')  ZNUM =  61
      IF (SYM .EQ. 'Sm')  ZNUM =  62
      IF (SYM .EQ. 'Eu')  ZNUM =  63
      IF (SYM .EQ. 'Gd')  ZNUM =  64
      IF (SYM .EQ. 'Tb')  ZNUM =  65
      IF (SYM .EQ. 'Dy')  ZNUM =  66
      IF (SYM .EQ. 'Ho')  ZNUM =  67
      IF (SYM .EQ. 'Er')  ZNUM =  68
      IF (SYM .EQ. 'Tm')  ZNUM =  69
      IF (SYM .EQ. 'Yb')  ZNUM =  70
      IF (SYM .EQ. 'Lu')  ZNUM =  71
      IF (SYM .EQ. 'Hf')  ZNUM =  72
      IF (SYM .EQ. 'Ta')  ZNUM =  73
      IF (SYM .EQ. 'W ')  ZNUM =  74
      IF (SYM .EQ. 'Re')  ZNUM =  75
      IF (SYM .EQ. 'Os')  ZNUM =  76
      IF (SYM .EQ. 'Ir')  ZNUM =  77
      IF (SYM .EQ. 'Pt')  ZNUM =  78
      IF (SYM .EQ. 'Au')  ZNUM =  79
      IF (SYM .EQ. 'Hg')  ZNUM =  80
      IF (SYM .EQ. 'Tl')  ZNUM =  81
      IF (SYM .EQ. 'Pb')  ZNUM =  82
      IF (SYM .EQ. 'Bi')  ZNUM =  83
      IF (SYM .EQ. 'Po')  ZNUM =  84
      IF (SYM .EQ. 'At')  ZNUM =  85
      IF (SYM .EQ. 'Rn')  ZNUM =  86
      IF (SYM .EQ. 'Fr')  ZNUM =  87
      IF (SYM .EQ. 'Ra')  ZNUM =  88
      IF (SYM .EQ. 'Ac')  ZNUM =  89
      IF (SYM .EQ. 'Th')  ZNUM =  90
      IF (SYM .EQ. 'Pa')  ZNUM =  91
      IF (SYM .EQ. 'U ')  ZNUM =  92
      IF (SYM .EQ. 'Np')  ZNUM =  93
      IF (SYM .EQ. 'Pu')  ZNUM =  94
      IF (SYM .EQ. 'Am')  ZNUM =  95
      IF (SYM .EQ. 'Cm')  ZNUM =  96
      IF (SYM .EQ. 'Bk')  ZNUM =  97
      IF (SYM .EQ. 'Cf')  ZNUM =  98
      IF (SYM .EQ. 'Es')  ZNUM =  99
      IF (SYM .EQ. 'Fm')  ZNUM = 100
      IF (SYM .EQ. 'Md')  ZNUM = 101
      IF (SYM .EQ. 'No')  ZNUM = 102
      IF (SYM .EQ. 'Lr')  ZNUM = 103
      IF (SYM .EQ. 'Rf')  ZNUM = 104
      IF (SYM .EQ. 'Db')  ZNUM = 105
      IF (SYM .EQ. 'Sg')  ZNUM = 106
      IF (SYM .EQ. 'Bh')  ZNUM = 107
      IF (SYM .EQ. 'Hs')  ZNUM = 108
      IF (SYM .EQ. 'Mt')  ZNUM = 109
      IF (SYM .EQ. 'Ds')  ZNUM = 110
      IF (SYM .EQ. 'Rg')  ZNUM = 111
      IF (SYM .EQ. 'Cn')  ZNUM = 112
      IF (SYM .EQ. 'Nh')  ZNUM = 113
      IF (SYM .EQ. 'Fl')  ZNUM = 114
      IF (SYM .EQ. 'Mc')  ZNUM = 115
      IF (SYM .EQ. 'Lv')  ZNUM = 116
      IF (SYM .EQ. 'Ts')  ZNUM = 117
      IF (SYM .EQ. 'Og')  ZNUM = 118

      IF (ZNUM .LE. 0) IERR = 1
      RETURN
      END


C=======================================================================
C  TRIM2
C  Take first 2 non-blank chars of ATOK into OUT (CHAR*2).
C=======================================================================
      SUBROUTINE TRIM2(ATOK, OUT)
      IMPLICIT NONE
      CHARACTER*(*) ATOK
      CHARACTER*2 OUT
      INTEGER I, L, K

      OUT = '  '
      L = LEN(ATOK)
      K = 0

      DO 10 I = 1, L
         IF (ATOK(I:I) .NE. ' ') THEN
            K = K + 1
            IF (K .LE. 2) OUT(K:K) = ATOK(I:I)
            IF (K .EQ. 2) RETURN
         ENDIF
   10 CONTINUE

      RETURN
      END


C=======================================================================
C  TOUP1 / TOLO1
C  Convert single character using ASCII code.
C=======================================================================
      SUBROUTINE TOUP1(CH)
      IMPLICIT NONE
      CHARACTER*1 CH
      INTEGER IC

      IC = ICHAR(CH)
      IF (IC .GE. ICHAR('a') .AND. IC .LE. ICHAR('z')) THEN
         CH = CHAR(IC - 32)
      ENDIF
      RETURN
      END

      SUBROUTINE TOLO1(CH)
      IMPLICIT NONE
      CHARACTER*1 CH
      INTEGER IC

      IC = ICHAR(CH)
      IF (IC .GE. ICHAR('A') .AND. IC .LE. ICHAR('Z')) THEN
         CH = CHAR(IC + 32)
      ENDIF
      RETURN
      END


