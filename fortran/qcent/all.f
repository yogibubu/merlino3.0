C======================================================================
C  QCENT_ALL.F
C
C  Harmonic quartic centrifugal distortion constants from Gaussian LOG/FCHK
C  + Watson representations (Ir, IIr, IIIr, Il, IIl, IIIl)
C  + Watson A and S reductions for chosen representation
C  + rotational constants A,B,C from principal moments
C
C  Requires:
C    - Gaussian output file (.log or .fchk) with harmonic modes
C    - mass.dat (NAtoms lines, masses in amu)
C
C  Compile:
C    gfortran -std=legacy -ffixed-line-length-none qcent_all.f -o qcent.x
C======================================================================

      PROGRAM TESTQCENT_HARM
      IMPLICIT REAL*8(A-H,O-Z)

      INTEGER MaxAt, MaxVib
      PARAMETER (MaxAt=300, MaxVib=900)

      INTEGER IOut, IPrint, IForm
      INTEGER NAtoms, NVib
      INTEGER IUnit
      INTEGER REP

      CHARACTER*200 FName

      REAL*8 Mass(MaxAt)
      REAL*8 Xeq(3,MaxAt)
      REAL*8 Freq(MaxVib)
      REAL*8 Lmw(3,MaxAt,MaxVib)

      REAL*8 PMom(3)
      LOGICAL AxOK(3)

      REAL*8 dIdQ(6,MaxVib)
      REAL*8 Tau(3,3,3,3)
      REAL*8 AlpXI(6,MaxVib,MaxVib)

      REAL*8 QA_sel_cm1(5), QA_sel_MHz(5)
      REAL*8 QS_sel_cm1(5), QS_sel_MHz(5)

      REAL*8 ABC_MHz(3), ABC_cm1(3)

      CALL INIT_PHYCON()

      IOut   = 6
      IPrint = 1

      WRITE(IOut,'(/,A)') ' TESTQCENT_HARM - quartic CD constants'
      WRITE(IOut,'(A)')   ' Enter Gaussian file name:'
      READ(*,'(A)') FName

      WRITE(IOut,'(A)') ' File type? 1=LOG, 2=FCHK'
      READ(*,*) IForm

      WRITE(IOut,'(A)') ' Enter NAtoms:'
      READ(*,*) NAtoms

      NVib = 3*NAtoms - 6
      WRITE(IOut,'(A,I6)') ' NVib = 3N-6 = ', NVib

      IF(NAtoms.GT.MaxAt) THEN
        WRITE(IOut,'(A)') ' ERROR: NAtoms > MaxAt'
        STOP
      ENDIF
      IF(NVib.GT.MaxVib) THEN
        WRITE(IOut,'(A)') ' ERROR: NVib > MaxVib'
        STOP
      ENDIF

      WRITE(IOut,'(/,A)') ' Choose representation:'
      WRITE(IOut,'(A)')   '   1 = Ir   (z x y)'
      WRITE(IOut,'(A)')   '   2 = IIr  (x y z)'
      WRITE(IOut,'(A)')   '   3 = IIIr (y z x)'
      WRITE(IOut,'(A)')   '   4 = Il   (z y x)'
      WRITE(IOut,'(A)')   '   5 = IIl  (y x z)'
      WRITE(IOut,'(A)')   '   6 = IIIl (x z y)'
      READ(*,*) REP

      IF(REP.LT.1 .OR. REP.GT.6) THEN
        WRITE(IOut,'(A)') ' Invalid REP, using IIr (=2)'
        REP=2
      ENDIF

      CALL READ_MASSDAT(NAtoms, Mass)

      IUnit = 10
      OPEN(IUnit,FILE=FName,STATUS='OLD',ERR=901)

      IF(IForm.EQ.1) THEN
        CALL READ_GAUSSIAN_LOG_HARM(IUnit, IOut, IPrint,
     $       NAtoms, NVib, Xeq, Freq, Lmw)
      ELSEIF(IForm.EQ.2) THEN
        CALL READ_GAUSSIAN_FCHK_HARM(IUnit, IOut, IPrint,
     $       NAtoms, NVib, Xeq, Freq, Lmw)
      ELSE
        WRITE(IOut,'(A)') ' ERROR: unknown IForm'
        STOP
      ENDIF

      CLOSE(IUnit)

C     LOG coords are Angstrom -> bohr
      IF(IForm.EQ.1) THEN
        CALL SCALE_GEOM(NAtoms, Xeq, Ang2Bohr)
      ENDIF

C     Build Tau (harmonic) + PMom etc.
      CALL QCENT_BUILD_TAU_ONLY(IOut, IPrint,
     $     NAtoms, NVib,
     $     Mass, Xeq, Freq, Lmw,
     $     AxOK, PMom,
     $     dIdQ, Tau, AlpXI)

C     rotational constants
      CALL PMOM_TO_ABC(PMom, AxOK, ABC_MHz, ABC_cm1)

C     chosen representation A/S
      CALL QCENT_ONE_REP(Tau, AxOK, REP,
     $     QA_sel_cm1, QA_sel_MHz,
     $     QS_sel_cm1, QS_sel_MHz)

      CALL PRINT_QCENT_ONE_REP(IOut, REP,
     $     QA_sel_cm1, QA_sel_MHz,
     $     QS_sel_cm1, QS_sel_MHz)

      WRITE(IOut,'(/,A)') ' Principal moments (amu*bohr^2):'
      WRITE(IOut,'(3(1X,1PE16.8))') PMom(1),PMom(2),PMom(3)

      CALL PRINT_ABC(IOut, ABC_MHz, ABC_cm1)

      CALL PRINT_SPFIT_BLOCK(IOut, ABC_MHz, QS_sel_MHz, QA_sel_MHz)

      STOP

  901 CONTINUE
      WRITE(IOut,'(A)') ' ERROR: cannot open Gaussian file'
      STOP
      END

C======================================================================
C  ROUTINES (alphabetical order)
C======================================================================

      SUBROUTINE BUILD_ALPXI_HARM_AU(NVib, Freq_au,
     $     PMom, AxOK, dIdQ, AlpXI)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER NVib
      INTEGER C,K,L

      REAL*8 Freq_au(NVib), PMom(3)
      LOGICAL AxOK(3)
      REAL*8 dIdQ(6,NVib)
      REAL*8 AlpXI(6,NVib,NVib)

      REAL*8 DEN

      DO C=1,6
        DO K=1,NVib
          DO L=1,NVib
            AlpXI(C,K,L)=0.0D0
          END DO
        END DO
      END DO

C     Harmonic-only simple contraction (kept for workflow completeness)
      DO C=1,6
        DO K=1,NVib
          DO L=1,NVib
            DEN = Freq_au(K)*Freq_au(L)
            IF(DEN.NE.0.0D0) THEN
              AlpXI(C,K,L)=dIdQ(C,K)*dIdQ(C,L)/DEN
            ENDIF
          END DO
        END DO
      END DO

      RETURN
      END


      SUBROUTINE BUILD_TAU_WILSON_AU(NVib, Freq_au,
     $     PMom, AxOK, dIdQ, Tau)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER NVib
      INTEGER IX,JX,KX,LX,IJX,KLX,IMODE

      REAL*8 Freq_au(NVib), PMom(3)
      LOGICAL AxOK(3)
      REAL*8 dIdQ(6,NVib)
      REAL*8 Tau(3,3,3,3)

      REAL*8 Pt5, Const1, Zero
      REAL*8 A, X
      PARAMETER(Zero=0.0D0, Pt5=0.5D0)

      Const1 = 1.0D0

C     init
      DO IX=1,3
        DO JX=1,3
          DO KX=1,3
            DO LX=1,3
              Tau(IX,JX,KX,LX)=0.0D0
            END DO
          END DO
        END DO
      END DO

C     Wilson Tau
      DO 100 IX=1,3
        IF(AxOK(IX)) THEN
          DO 110 JX=1,3
            IF(AxOK(JX)) THEN
              DO 120 KX=1,3
                IF(AxOK(KX)) THEN
                  DO 130 LX=1,3
                    IF(AxOK(LX)) THEN

                      IJX = IJ2SYM6(IX,JX)
                      KLX = IJ2SYM6(KX,LX)

                      A = Zero
                      DO 140 IMODE=1,NVib
                        X = Freq_au(IMODE)**2 *
     $                      PMom(IX)*PMom(JX)*PMom(KX)*PMom(LX)
                        A = A + dIdQ(IJX,IMODE)*dIdQ(KLX,IMODE)/X
  140                 CONTINUE

                      Tau(IX,JX,KX,LX) = -Pt5*Const1*A

                    ENDIF
  130             CONTINUE
                ENDIF
  120         CONTINUE
              ENDIF
  110     CONTINUE
            ENDIF
  100 CONTINUE

      RETURN
      END


      SUBROUTINE DIDQ_SYM6(NAtoms, NVib, Mass, Xeq, Lmw, dIdQ)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER NAtoms, NVib
      INTEGER IAT, K, C

      REAL*8 Mass(NAtoms), Xeq(3,NAtoms)
      REAL*8 Lmw(3,NAtoms,NVib)
      REAL*8 dIdQ(6,NVib)

      REAL*8 M, SQRM
      REAL*8 R1,R2,R3, DR1,DR2,DR3
      REAL*8 RDOTDR
      REAL*8 dI11,dI22,dI33,dI12,dI13,dI23

      DO K=1,NVib
        DO C=1,6
          dIdQ(C,K)=0.0D0
        END DO
      END DO

      DO K=1,NVib
        DO IAT=1,NAtoms

          M    = Mass(IAT)
          SQRM = DSQRT(M)

          R1 = Xeq(1,IAT)
          R2 = Xeq(2,IAT)
          R3 = Xeq(3,IAT)

C         dr/dQ from mass-weighted normal modes
          DR1 = Lmw(1,IAT,K)/SQRM
          DR2 = Lmw(2,IAT,K)/SQRM
          DR3 = Lmw(3,IAT,K)/SQRM

          RDOTDR = R1*DR1 + R2*DR2 + R3*DR3

          dI11 = M*( 2.0D0*RDOTDR - 2.0D0*R1*DR1 )
          dI22 = M*( 2.0D0*RDOTDR - 2.0D0*R2*DR2 )
          dI33 = M*( 2.0D0*RDOTDR - 2.0D0*R3*DR3 )

          dI12 = M*( -R1*DR2 - R2*DR1 )
          dI13 = M*( -R1*DR3 - R3*DR1 )
          dI23 = M*( -R2*DR3 - R3*DR2 )

          dIdQ(1,K)=dIdQ(1,K)+dI11
          dIdQ(2,K)=dIdQ(2,K)+dI22
          dIdQ(3,K)=dIdQ(3,K)+dI33
          dIdQ(4,K)=dIdQ(4,K)+dI12
          dIdQ(5,K)=dIdQ(5,K)+dI13
          dIdQ(6,K)=dIdQ(6,K)+dI23

        END DO
      END DO

      RETURN
      END


      SUBROUTINE FREQ_CM1_TO_AU(NVib, Fcm1, Fau)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER NVib, K
      REAL*8 Fcm1(NVib), Fau(NVib)

      COMMON /PHYCON/
     $    Pi, TwoPi, Avog, Clight, Planck, Hartree2J,
     $    Bohr2Ang, Ang2Bohr, AMU2AU, AU2AMU,
     $    Cm12AU, AU2Cm1, Cm12MHz

      REAL*8 Pi, TwoPi, Avog, Clight, Planck, Hartree2J
      REAL*8 Bohr2Ang, Ang2Bohr, AMU2AU, AU2AMU
      REAL*8 Cm12AU, AU2Cm1, Cm12MHz

      DO K=1,NVib
        Fau(K)=Fcm1(K)*Cm12AU
      END DO

      RETURN
      END


      SUBROUTINE GET_REP_PERM(REP, Perm)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER REP, Perm(3)

      IF(REP.EQ.1) THEN
C       Ir = (z,x,y) -> (3,1,2)
        Perm(1)=3
        Perm(2)=1
        Perm(3)=2
      ELSEIF(REP.EQ.2) THEN
C       IIr = (x,y,z) -> (1,2,3)
        Perm(1)=1
        Perm(2)=2
        Perm(3)=3
      ELSEIF(REP.EQ.3) THEN
C       IIIr = (y,z,x) -> (2,3,1)
        Perm(1)=2
        Perm(2)=3
        Perm(3)=1
      ELSEIF(REP.EQ.4) THEN
C       Il = (z,y,x) -> (3,2,1)
        Perm(1)=3
        Perm(2)=2
        Perm(3)=1
      ELSEIF(REP.EQ.5) THEN
C       IIl = (y,x,z) -> (2,1,3)
        Perm(1)=2
        Perm(2)=1
        Perm(3)=3
      ELSEIF(REP.EQ.6) THEN
C       IIIl = (x,z,y) -> (1,3,2)
        Perm(1)=1
        Perm(2)=3
        Perm(3)=2
      ELSE
        Perm(1)=1
        Perm(2)=2
        Perm(3)=3
      ENDIF

      RETURN
      END


      SUBROUTINE INERTIA_TENSOR(NAtoms, Mass, X, I)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER NAtoms
      INTEGER IAT, IA, IB

      REAL*8 Mass(NAtoms), X(3,NAtoms), I(3,3)
      REAL*8 M, X1, Y1, Z1

      DO IA=1,3
        DO IB=1,3
          I(IA,IB)=0.0D0
        END DO
      END DO

      DO IAT=1,NAtoms
        M  = Mass(IAT)
        X1 = X(1,IAT)
        Y1 = X(2,IAT)
        Z1 = X(3,IAT)

        I(1,1) = I(1,1) + M*(Y1*Y1 + Z1*Z1)
        I(2,2) = I(2,2) + M*(X1*X1 + Z1*Z1)
        I(3,3) = I(3,3) + M*(X1*X1 + Y1*Y1)

        I(1,2) = I(1,2) - M*(X1*Y1)
        I(1,3) = I(1,3) - M*(X1*Z1)
        I(2,3) = I(2,3) - M*(Y1*Z1)
      END DO

      I(2,1)=I(1,2)
      I(3,1)=I(1,3)
      I(3,2)=I(2,3)

      RETURN
      END


      SUBROUTINE INIT_PHYCON()
      IMPLICIT REAL*8(A-H,O-Z)

      COMMON /PHYCON/
     $    Pi, TwoPi, Avog, Clight, Planck, Hartree2J,
     $    Bohr2Ang, Ang2Bohr, AMU2AU, AU2AMU,
     $    Cm12AU, AU2Cm1, Cm12MHz

      REAL*8 Pi, TwoPi, Avog, Clight, Planck, Hartree2J
      REAL*8 Bohr2Ang, Ang2Bohr, AMU2AU, AU2AMU
      REAL*8 Cm12AU, AU2Cm1, Cm12MHz

      Pi     = 4.0D0*DATAN(1.0D0)
      TwoPi  = 2.0D0*Pi

      Clight = 2.99792458D10
      Planck = 6.62607015D-34
      Avog   = 6.02214076D23

      Hartree2J = 4.3597447222071D-18

      Bohr2Ang = 0.529177210903D0
      Ang2Bohr = 1.0D0/Bohr2Ang

      AMU2AU = 1822.888486209D0
      AU2AMU = 1.0D0/AMU2AU

      AU2Cm1 = 219474.6313705D0
      Cm12AU = 1.0D0/AU2Cm1

      Cm12MHz = 29979.2458D0

      RETURN
      END


      INTEGER FUNCTION IJ2SYM6(I,J)
      INTEGER I,J
      IF(I.EQ.1 .AND. J.EQ.1) THEN
        IJ2SYM6=1
      ELSEIF(I.EQ.2 .AND. J.EQ.2) THEN
        IJ2SYM6=2
      ELSEIF(I.EQ.3 .AND. J.EQ.3) THEN
        IJ2SYM6=3
      ELSEIF((I.EQ.1 .AND. J.EQ.2) .OR. (I.EQ.2 .AND. J.EQ.1)) THEN
        IJ2SYM6=4
      ELSEIF((I.EQ.1 .AND. J.EQ.3) .OR. (I.EQ.3 .AND. J.EQ.1)) THEN
        IJ2SYM6=5
      ELSEIF((I.EQ.2 .AND. J.EQ.3) .OR. (I.EQ.3 .AND. J.EQ.2)) THEN
        IJ2SYM6=6
      ELSE
        IJ2SYM6=0
      ENDIF
      RETURN
      END


      SUBROUTINE JACOBI3(A, D, V)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER IP, IQ, J, N

      REAL*8 A(3,3), D(3), V(3,3)
      REAL*8 TRESH, THETA, TAU, T, SM, S, H, G, C

      DO IP=1,3
        DO IQ=1,3
          V(IP,IQ)=0.0D0
        END DO
        V(IP,IP)=1.0D0
      END DO

      DO IP=1,3
        D(IP)=A(IP,IP)
      END DO

      DO N=1,50

        SM = DABS(A(1,2)) + DABS(A(1,3)) + DABS(A(2,3))
        IF(SM.LT.1.0D-15) RETURN

        TRESH = 0.0D0
        IF(N.LT.4) TRESH = 0.2D0*SM/9.0D0

        DO IP=1,2
          DO IQ=IP+1,3

            G = 100.0D0*DABS(A(IP,IQ))

            IF(N.GT.4 .AND.
     $       (DABS(D(IP))+G).EQ.DABS(D(IP)) .AND.
     $       (DABS(D(IQ))+G).EQ.DABS(D(IQ))) THEN

              A(IP,IQ)=0.0D0

            ELSE IF(DABS(A(IP,IQ)).GT.TRESH) THEN

              H = D(IQ) - D(IP)

              IF((DABS(H)+G).EQ.DABS(H)) THEN
                T = A(IP,IQ)/H
              ELSE
                THETA = 0.5D0*H/A(IP,IQ)
                T = 1.0D0/(DABS(THETA)+DSQRT(1.0D0+THETA*THETA))
                IF(THETA.LT.0.0D0) T = -T
              ENDIF

              C = 1.0D0/DSQRT(1.0D0+T*T)
              S = T*C
              TAU = S/(1.0D0+C)

              H = T*A(IP,IQ)
              A(IP,IQ)=0.0D0

              D(IP)=D(IP)-H
              D(IQ)=D(IQ)+H

              DO J=1,IP-1
                CALL ROT(A, J,IP, J,IQ, S, TAU)
              END DO
              DO J=IP+1,IQ-1
                CALL ROT(A, IP,J, J,IQ, S, TAU)
              END DO
              DO J=IQ+1,3
                CALL ROT(A, IP,J, IQ,J, S, TAU)
              END DO
              DO J=1,3
                CALL ROT(V, J,IP, J,IQ, S, TAU)
              END DO

            ENDIF

          END DO
        END DO

      END DO

      RETURN
      END


      SUBROUTINE PMOM_TO_ABC(PMom, AxOK, ABC_MHz, ABC_cm1)
      IMPLICIT REAL*8(A-H,O-Z)

      LOGICAL AxOK(3)
      REAL*8 PMom(3), ABC_MHz(3), ABC_cm1(3)

      COMMON /PHYCON/
     $    Pi, TwoPi, Avog, Clight, Planck, Hartree2J,
     $    Bohr2Ang, Ang2Bohr, AMU2AU, AU2AMU,
     $    Cm12AU, AU2Cm1, Cm12MHz

      REAL*8 Pi, TwoPi, Avog, Clight, Planck, Hartree2J
      REAL*8 Bohr2Ang, Ang2Bohr, AMU2AU, AU2AMU
      REAL*8 Cm12AU, AU2Cm1, Cm12MHz

      REAL*8 AMU_kg, BOHR_m
      REAL*8 I_SI, B_Hz
      REAL*8 EIGHTPI2
      INTEGER IAX

      AMU_kg = 1.66053906660D-27
      BOHR_m = 0.529177210903D-10
      EIGHTPI2 = 8.0D0*Pi*Pi

      DO IAX=1,3
        ABC_MHz(IAX)=0.0D0
        ABC_cm1(IAX)=0.0D0

        IF(AxOK(IAX)) THEN
          I_SI = PMom(IAX) * AMU_kg * (BOHR_m*BOHR_m)
          B_Hz = Planck / (EIGHTPI2 * I_SI)
          ABC_MHz(IAX) = B_Hz * 1.0D-6
          ABC_cm1(IAX) = B_Hz / Clight
        ENDIF
      END DO

      RETURN
      END


      SUBROUTINE PERMUTE_TAU(TauIn, Perm, TauOut)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER Perm(3)
      INTEGER I,J,K,L
      REAL*8 TauIn(3,3,3,3), TauOut(3,3,3,3)

      DO I=1,3
        DO J=1,3
          DO K=1,3
            DO L=1,3
              TauOut(I,J,K,L) =
     $          TauIn(Perm(I),Perm(J),Perm(K),Perm(L))
            END DO
          END DO
        END DO
      END DO

      RETURN
      END


      SUBROUTINE PRINCIPAL_MOMENTS(I0, PMom, AxOK)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER IA, IB
      REAL*8 I0(3,3), PMom(3)
      LOGICAL AxOK(3)

      REAL*8 A(3,3), Eig(3), V(3,3)

      DO IA=1,3
        DO IB=1,3
          A(IA,IB)=I0(IA,IB)
        END DO
      END DO

      CALL JACOBI3(A, Eig, V)
      CALL SORT3(Eig)

      PMom(1)=Eig(1)
      PMom(2)=Eig(2)
      PMom(3)=Eig(3)

      DO IA=1,3
        AxOK(IA)=.TRUE.
        IF(PMom(IA).LT.1.0D-14) AxOK(IA)=.FALSE.
      END DO

      RETURN
      END


      SUBROUTINE PRINT_ABC(IOut, ABC_MHz, ABC_cm1)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER IOut
      REAL*8 ABC_MHz(3), ABC_cm1(3)

      WRITE(IOut,'(/,A)') ' Rotational constants (from PMom):'
      WRITE(IOut,'(A)')   '   Axis      MHz                 cm^-1'
      WRITE(IOut,'(A,2X,1PE16.8,3X,1PE16.8)') '   a   ',
     $     ABC_MHz(1), ABC_cm1(1)
      WRITE(IOut,'(A,2X,1PE16.8,3X,1PE16.8)') '   b   ',
     $     ABC_MHz(2), ABC_cm1(2)
      WRITE(IOut,'(A,2X,1PE16.8,3X,1PE16.8)') '   c   ',
     $     ABC_MHz(3), ABC_cm1(3)

      RETURN
      END


      SUBROUTINE PRINT_QCENT_ONE_REP(IOut, REP,
     $     QA_cm1, QA_MHz, QS_cm1, QS_MHz)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER IOut, REP

      REAL*8 QA_cm1(5), QA_MHz(5)
      REAL*8 QS_cm1(5), QS_MHz(5)

      CHARACTER*5 RNAME(6)
      DATA RNAME /'Ir   ','IIr  ','IIIr ','Il   ','IIl  ','IIIl '/

      WRITE(IOut,'(/,A,A)') ' == Selected representation: ', RNAME(REP)

      WRITE(IOut,'(A)') ' A-reduction '
      WRITE(IOut,'(A,5(1X,1PE16.8))') '   cm^-1:',
     $   QA_cm1(1),QA_cm1(2),QA_cm1(3),QA_cm1(4),QA_cm1(5)
      WRITE(IOut,'(A,5(1X,1PE16.8))') '   MHz :',
     $   QA_MHz(1),QA_MHz(2),QA_MHz(3),QA_MHz(4),QA_MHz(5)

      WRITE(IOut,'(A)') ' S-reduction (DJ DJK DK d1 d2)'
      WRITE(IOut,'(A,5(1X,1PE16.8))') '   cm^-1:',
     $   QS_cm1(1),QS_cm1(2),QS_cm1(3),QS_cm1(4),QS_cm1(5)
      WRITE(IOut,'(A,5(1X,1PE16.8))') '   MHz :',
     $   QS_MHz(1),QS_MHz(2),QS_MHz(3),QS_MHz(4),QS_MHz(5)

      RETURN
      END


      SUBROUTINE PRINT_SPFIT_BLOCK(IOut, ABC_MHz,
     $     QS_MHz, QA_MHz)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER IOut
      REAL*8 ABC_MHz(3)
      REAL*8 QS_MHz(5)
      REAL*8 QA_MHz(5)

      WRITE(IOut,'(/,A)') ' --- SPFIT-like block (MHz) ---'
      WRITE(IOut,'(A,1X,1PE16.8)') ' A  =', ABC_MHz(1)
      WRITE(IOut,'(A,1X,1PE16.8)') ' B  =', ABC_MHz(2)
      WRITE(IOut,'(A,1X,1PE16.8)') ' C  =', ABC_MHz(3)

      WRITE(IOut,'(/,A)') ' S-reduction [MHz]'
      WRITE(IOut,'(A,1X,1PE16.8)') ' DJ  =', QS_MHz(1)
      WRITE(IOut,'(A,1X,1PE16.8)') ' DJK =', QS_MHz(2)
      WRITE(IOut,'(A,1X,1PE16.8)') ' DK  =', QS_MHz(3)
      WRITE(IOut,'(A,1X,1PE16.8)') ' d1  =', QS_MHz(4)
      WRITE(IOut,'(A,1X,1PE16.8)') ' d2  =', QS_MHz(5)

      WRITE(IOut,'(/,A)') ' A-reduction [MHz]'
      WRITE(IOut,'(A,1X,1PE16.8)') ' DelJ  =', QA_MHz(1)
      WRITE(IOut,'(A,1X,1PE16.8)') ' DelJK =', QA_MHz(2)
      WRITE(IOut,'(A,1X,1PE16.8)') ' DelK  =', QA_MHz(3)
      WRITE(IOut,'(A,1X,1PE16.8)') ' delJ  =', QA_MHz(4)
      WRITE(IOut,'(A,1X,1PE16.8)') ' delK  =', QA_MHz(5)

      RETURN
      END


      SUBROUTINE QCENT_ARED_TO_SRED_IR(QA_AU, QS_AU)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 QA_AU(5), QS_AU(5)

      REAL*8 DELTAJ, DELTAJK, DELTAK, SMALLJ, SMALLK

      DELTAJ  = QA_AU(1)
      DELTAJK = QA_AU(2)
      DELTAK  = QA_AU(3)
      SMALLJ  = QA_AU(4)
      SMALLK  = QA_AU(5)

      QS_AU(1) = DELTAJ
      QS_AU(2) = DELTAJK + DELTAJ
      QS_AU(3) = DELTAK  + DELTAJK + DELTAJ
      QS_AU(4) = SMALLJ
      QS_AU(5) = SMALLK

      RETURN
      END


      SUBROUTINE QCENT_AU_TO_CM1(Qcent)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 Qcent(5)
      INTEGER I

      COMMON /PHYCON/
     $    Pi, TwoPi, Avog, Clight, Planck, Hartree2J,
     $    Bohr2Ang, Ang2Bohr, AMU2AU, AU2AMU,
     $    Cm12AU, AU2Cm1, Cm12MHz

      REAL*8 Pi, TwoPi, Avog, Clight, Planck, Hartree2J
      REAL*8 Bohr2Ang, Ang2Bohr, AMU2AU, AU2AMU
      REAL*8 Cm12AU, AU2Cm1, Cm12MHz

      DO I=1,5
        Qcent(I)=Qcent(I)*AU2Cm1
      END DO

      RETURN
      END


      SUBROUTINE QCENT_BUILD_TAU_ONLY(IOut, IPrint,
     $     NAtoms, NVib,
     $     Mass, Xeq, Freq_cm1, Lmw,
     $     AxOK, PMom,
     $     dIdQ, Tau, AlpXI)
      IMPLICIT REAL*8(A-H,O-Z)

      INTEGER IOut, IPrint
      INTEGER NAtoms, NVib

      REAL*8 Mass(NAtoms)
      REAL*8 Xeq(3,NAtoms)
      REAL*8 Freq_cm1(NVib)
      REAL*8 Lmw(3,NAtoms,NVib)

      LOGICAL AxOK(3)
      REAL*8 PMom(3)

      REAL*8 dIdQ(6,NVib)
      REAL*8 Tau(3,3,3,3)
      REAL*8 AlpXI(6,NVib,NVib)

      REAL*8 I0(3,3)
      REAL*8 Freq_au(900)

      CALL INERTIA_TENSOR(NAtoms, Mass, Xeq, I0)
      CALL PRINCIPAL_MOMENTS(I0, PMom, AxOK)

      CALL FREQ_CM1_TO_AU(NVib, Freq_cm1, Freq_au)

      CALL DIDQ_SYM6(NAtoms, NVib, Mass, Xeq, Lmw, dIdQ)

      CALL BUILD_TAU_WILSON_AU(NVib, Freq_au, PMom, AxOK, dIdQ, Tau)

      CALL BUILD_ALPXI_HARM_AU(NVib, Freq_au, PMom, AxOK, dIdQ, AlpXI)

      IF(IPrint.GT.0) THEN
        WRITE(IOut,'(A)') ' QCENT_BUILD_TAU_ONLY: done'
      ENDIF

      RETURN
      END


      SUBROUTINE QCENT_CM1_TO_MHZ(Qcm1, QMHz)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 Qcm1(5), QMHz(5)
      INTEGER I

      COMMON /PHYCON/
     $    Pi, TwoPi, Avog, Clight, Planck, Hartree2J,
     $    Bohr2Ang, Ang2Bohr, AMU2AU, AU2AMU,
     $    Cm12AU, AU2Cm1, Cm12MHz

      REAL*8 Pi, TwoPi, Avog, Clight, Planck, Hartree2J
      REAL*8 Bohr2Ang, Ang2Bohr, AMU2AU, AU2AMU
      REAL*8 Cm12AU, AU2Cm1, Cm12MHz

      DO I=1,5
        QMHz(I)=Qcm1(I)*Cm12MHz
      END DO

      RETURN
      END


      SUBROUTINE QCENT_ONE_REP(Tau, AxOK, REP,
     $     QA_cm1, QA_MHz, QS_cm1, QS_MHz)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL AxOK(3)
      INTEGER REP, Perm(3)
      INTEGER I

      REAL*8 Tau(3,3,3,3)
      REAL*8 TauP(3,3,3,3)

      REAL*8 QA_cm1(5), QA_MHz(5)
      REAL*8 QS_cm1(5), QS_MHz(5)

      REAL*8 QA_au(5), QS_au(5)

      CALL GET_REP_PERM(REP, Perm)
      CALL PERMUTE_TAU(Tau, Perm, TauP)

      CALL TAU_TO_QCENT_ARED_IR(TauP, AxOK, QA_au)

      CALL QCENT_ARED_TO_SRED_IR(QA_au, QS_au)

      DO I=1,5
        QA_cm1(I)=QA_au(I)
        QS_cm1(I)=QS_au(I)
      END DO

      CALL QCENT_AU_TO_CM1(QA_cm1)
      CALL QCENT_AU_TO_CM1(QS_cm1)

      CALL QCENT_CM1_TO_MHZ(QA_cm1, QA_MHz)
      CALL QCENT_CM1_TO_MHZ(QS_cm1, QS_MHz)

      RETURN
      END


      SUBROUTINE READ_DISP3(Line, KMODE, NVib, IAT, NAtoms, Lmw)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER KMODE, NVib, IAT, NAtoms

      CHARACTER*(*) Line
      REAL*8 Lmw(3,NAtoms,NVib)

      INTEGER IA, AN
      REAL*8 X1,Y1,Z1, X2,Y2,Z2, X3,Y3,Z3

      X1=0.0D0
      Y1=0.0D0
      Z1=0.0D0
      X2=0.0D0
      Y2=0.0D0
      Z2=0.0D0
      X3=0.0D0
      Y3=0.0D0
      Z3=0.0D0

      READ(Line,*,ERR=99) IA, AN,
     $   X1,Y1,Z1, X2,Y2,Z2, X3,Y3,Z3

      IF(KMODE+1.LE.NVib) THEN
        Lmw(1,IAT,KMODE+1)=X1
        Lmw(2,IAT,KMODE+1)=Y1
        Lmw(3,IAT,KMODE+1)=Z1
      ENDIF

      IF(KMODE+2.LE.NVib) THEN
        Lmw(1,IAT,KMODE+2)=X2
        Lmw(2,IAT,KMODE+2)=Y2
        Lmw(3,IAT,KMODE+2)=Z2
      ENDIF

      IF(KMODE+3.LE.NVib) THEN
        Lmw(1,IAT,KMODE+3)=X3
        Lmw(2,IAT,KMODE+3)=Y3
        Lmw(3,IAT,KMODE+3)=Z3
      ENDIF

      RETURN
  99  CONTINUE
      RETURN
      END


      SUBROUTINE READ_FREQ3(Line, KMODE, NVib, Freq)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER KMODE, NVib
      CHARACTER*(*) Line
      REAL*8 Freq(NVib)
      REAL*8 F1,F2,F3

      F1=0.0D0
      F2=0.0D0
      F3=0.0D0

      READ(Line(16:),*,ERR=99) F1,F2,F3

      IF(KMODE+1.LE.NVib) Freq(KMODE+1)=F1
      IF(KMODE+2.LE.NVib) Freq(KMODE+2)=F2
      IF(KMODE+3.LE.NVib) Freq(KMODE+3)=F3

      RETURN
  99  CONTINUE
      RETURN
      END


      SUBROUTINE READ_GAUSSIAN_FCHK_HARM(IUnit, IOut, IPrint,
     $     NAtoms, NVib, Xeq, Freq, Lmw)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER IUnit, IOut, IPrint
      INTEGER NAtoms, NVib
      INTEGER IAT, K, N, NTOT

      REAL*8 Xeq(3,NAtoms), Freq(NVib), Lmw(3,NAtoms,NVib)

      CHARACTER*200 Line
      INTEGER FoundCoord, FoundFreq, FoundModes

      FoundCoord=0
      FoundFreq=0
      FoundModes=0

      REWIND(IUnit)

  10  CONTINUE
      READ(IUnit,'(A)',END=200) Line

      IF(INDEX(Line,'Current cartesian coordinates').GT.0) THEN
        FoundCoord=1
        NTOT = 3*NAtoms
        DO N=1,NTOT
          READ(IUnit,*) Xeq( MOD(N-1,3)+1, (N-1)/3 + 1 )
        END DO
      ENDIF

      IF(INDEX(Line,'Vibrational Frequencies').GT.0) THEN
        FoundFreq=1
        DO K=1,NVib
          READ(IUnit,*) Freq(K)
        END DO
      ENDIF

      IF(INDEX(Line,'Vibrational Normal Modes').GT.0) THEN
        FoundModes=1
        DO K=1,NVib
          DO IAT=1,NAtoms
            READ(IUnit,*) Lmw(1,IAT,K)
            READ(IUnit,*) Lmw(2,IAT,K)
            READ(IUnit,*) Lmw(3,IAT,K)
          END DO
        END DO
      ENDIF

      GOTO 10

  200 CONTINUE

      IF(FoundCoord.EQ.0) THEN
        WRITE(IOut,'(A)') 'ERROR: FCHK coords not found'
        STOP
      ENDIF
      IF(FoundFreq.EQ.0) THEN
        WRITE(IOut,'(A)') 'ERROR: FCHK frequencies not found'
        STOP
      ENDIF
      IF(FoundModes.EQ.0) THEN
        WRITE(IOut,'(A)') 'ERROR: FCHK normal modes not found'
        STOP
      ENDIF

      IF(IPrint.GT.0) THEN
        WRITE(IOut,'(A)') ' FCHK harmonic reader: OK'
      ENDIF

      RETURN
      END


      SUBROUTINE READ_GAUSSIAN_LOG_HARM(IUnit, IOut, IPrint,
     $     NAtoms, NVib, Xeq, Freq, Lmw)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER IUnit, IOut, IPrint
      INTEGER NAtoms, NVib
      INTEGER IAT, KMODE, IA, ATN, ATTYPE

      REAL*8 Xeq(3,NAtoms), Freq(NVib), Lmw(3,NAtoms,NVib)
      CHARACTER*200 Line

      DO K=1,NVib
        Freq(K)=0.0D0
        DO IAT=1,NAtoms
          Lmw(1,IAT,K)=0.0D0
          Lmw(2,IAT,K)=0.0D0
          Lmw(3,IAT,K)=0.0D0
        END DO
      END DO

C     last Standard orientation
      REWIND(IUnit)

  10  CONTINUE
      READ(IUnit,'(A)',END=100) Line
      IF(INDEX(Line,'Standard orientation').GT.0) THEN

        READ(IUnit,'(A)') Line
        READ(IUnit,'(A)') Line
        READ(IUnit,'(A)') Line
        READ(IUnit,'(A)') Line
        READ(IUnit,'(A)') Line

        DO IAT=1,NAtoms
          READ(IUnit,*) IA, ATN, ATTYPE,
     $      Xeq(1,IAT), Xeq(2,IAT), Xeq(3,IAT)
        END DO

      ENDIF
      GOTO 10

  100 CONTINUE

C     modes
      REWIND(IUnit)
      KMODE=0

  20  CONTINUE
      READ(IUnit,'(A)',END=200) Line

      IF(INDEX(Line,'Frequencies --').GT.0) THEN

        CALL READ_FREQ3(Line, KMODE, NVib, Freq)

  30    CONTINUE
        READ(IUnit,'(A)',END=200) Line
        IF(INDEX(Line,'Atom').GT.0 .AND. INDEX(Line,'AN').GT.0) THEN

          READ(IUnit,'(A)') Line

          DO IAT=1,NAtoms
            READ(IUnit,'(A)') Line
            CALL READ_DISP3(Line, KMODE, NVib, IAT, NAtoms, Lmw)
          END DO

          KMODE = KMODE + 3

        ELSE
          GOTO 30
        ENDIF

      ENDIF

      GOTO 20

  200 CONTINUE

      IF(IPrint.GT.0) THEN
        WRITE(IOut,'(A)') ' LOG harmonic reader: OK'
      ENDIF

      RETURN
      END


      SUBROUTINE READ_MASSDAT(NAtoms, Mass)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER NAtoms, I
      REAL*8 Mass(NAtoms)

      OPEN(99,FILE='mass.dat',STATUS='OLD',ERR=901)
      DO I=1,NAtoms
        READ(99,*,END=902) Mass(I)
      END DO
      CLOSE(99)
      RETURN

  901 CONTINUE
      WRITE(6,'(A)') ' ERROR: cannot open mass.dat'
      STOP
  902 CONTINUE
      WRITE(6,'(A)') ' ERROR: not enough masses in mass.dat'
      STOP
      END


      SUBROUTINE ROT(A, I, J, K, L, S, TAU)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER I, J, K, L
      REAL*8 A(3,3), S, TAU
      REAL*8 G, H

      G = A(I,J)
      H = A(K,L)

      A(I,J) = G - S*(H + G*TAU)
      A(K,L) = H + S*(G - H*TAU)

      RETURN
      END


      SUBROUTINE SCALE_GEOM(NAtoms, X, Fact)
      IMPLICIT REAL*8(A-H,O-Z)
      INTEGER NAtoms, I
      REAL*8 X(3,NAtoms), Fact

      DO I=1,NAtoms
        X(1,I)=X(1,I)*Fact
        X(2,I)=X(2,I)*Fact
        X(3,I)=X(3,I)*Fact
      END DO

      RETURN
      END


      SUBROUTINE SORT3(X)
      IMPLICIT REAL*8(A-H,O-Z)
      REAL*8 X(3), TMP

      IF(X(1).GT.X(2)) THEN
        TMP=X(1)
        X(1)=X(2)
        X(2)=TMP
      ENDIF
      IF(X(2).GT.X(3)) THEN
        TMP=X(2)
        X(2)=X(3)
        X(3)=TMP
      ENDIF
      IF(X(1).GT.X(2)) THEN
        TMP=X(1)
        X(1)=X(2)
        X(2)=TMP
      ENDIF

      RETURN
      END


      SUBROUTINE TAU_TO_QCENT_ARED_IR(Tau, AxOK, Q_AU)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL AxOK(3)
      REAL*8 Tau(3,3,3,3)
      REAL*8 Q_AU(5)

      INTEGER I
      REAL*8 DELTAJ, DELTAJK, DELTAK, SMALLJ, SMALLK
      REAL*8 Pt25, Pt125
      PARAMETER(Pt25=0.25D0, Pt125=0.125D0)

      IF(.NOT.(AxOK(1).AND.AxOK(2).AND.AxOK(3))) THEN
        DO I=1,5
          Q_AU(I)=0.0D0
        END DO
        RETURN
      ENDIF

      DELTAJ  = Pt125*( Tau(2,2,2,2) + Tau(3,3,3,3)
     $                + 2.0D0*Tau(2,2,3,3) )

      DELTAK  = Pt125*Tau(1,1,1,1)

      DELTAJK = -Pt25*( Tau(1,1,2,2) + Tau(1,1,3,3) )

      SMALLJ  = Pt125*( Tau(2,2,2,2) + Tau(3,3,3,3)
     $                - 2.0D0*Tau(2,2,3,3) )

      SMALLK  = -Pt25*( Tau(1,2,1,2) + Tau(1,3,1,3) )

      Q_AU(1)=DELTAJ
      Q_AU(2)=DELTAJK
      Q_AU(3)=DELTAK
      Q_AU(4)=SMALLJ
      Q_AU(5)=SMALLK

      RETURN
      END



