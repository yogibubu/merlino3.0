c=======================================================================
c BDPCS3 updated backend, aligned with survibfit.modify_geom.py
c
c The numerical model is intentionally the same as the Python function
c bdpcs3_delta_and_order_updated:
c   - r_cov: H/C/O/S explicit BDPCS3 values, otherwise Mantina-Truhlar
c   - f_coord = 0.5*(1-erf((r - 1.3*r_cov)/0.057))
c   - DeltaCV = -0.0025*(max(n_i,n_j)-1), with n capped at 3
c   - C-C and C-S delocalization compensation from Pyykko radii
c   - bond order is reported as max(override, exp((r_cov-r)/0.30))
c   - H...Y hydrogen-bond targets use the same angle gate and distance
c     damping as the Python backend.
c=======================================================================
      PROGRAM BDPCS3_XYZ_DRIVER
      INTEGER MAXAT
      PARAMETER (MAXAT=2000)
      INTEGER N,I,J,IZ
      INTEGER Z(MAXAT)
      DOUBLE PRECISION XYZ(3,MAXAT)
      DOUBLE PRECISION DX,DY,DZ,R,DELTA,RBD,BO
      DOUBLE PRECISION RX,RY,RZ
      CHARACTER*256 FNAME,TITLE
      CHARACTER*2 SYM,ZSYM
      INTEGER ATOMZ

      WRITE(*,*) 'Enter XYZ file name:'
      READ(*,'(A)') FNAME
      OPEN(10,FILE=FNAME,STATUS='OLD',ERR=900)

      READ(10,*,ERR=901,END=901) N
      IF (N .GT. MAXAT) THEN
         WRITE(*,*) 'Error: N > MAXAT.'
         STOP
      ENDIF
      READ(10,'(A)',ERR=902,END=902) TITLE

      DO 10 I=1,N
         READ(10,*,ERR=902,END=902) SYM,RX,RY,RZ
         IZ=ATOMZ(SYM)
         IF (IZ .LE. 0) THEN
            WRITE(*,*) 'Error: unknown element at atom ',I
            STOP
         ENDIF
         Z(I)=IZ
         XYZ(1,I)=RX
         XYZ(2,I)=RY
         XYZ(3,I)=RZ
   10 CONTINUE
      CLOSE(10)

      WRITE(*,'(/,A)') 'BDPCS3 updated backend (Python-aligned)'
      WRITE(*,'(A)') 'Units: Angstrom'
      WRITE(*,'(/,A)')
     & '  i   j  Zi  Zj  Si Sj       r_DPCS3       Delta      r_BDPCS3          BO'
      WRITE(*,'(A)')
     & '----------------------------------------------------------------------------'

      DO 30 I=1,N-1
         DO 20 J=I+1,N
            DX=XYZ(1,J)-XYZ(1,I)
            DY=XYZ(2,J)-XYZ(2,I)
            DZ=XYZ(3,J)-XYZ(3,I)
            R=DSQRT(DX*DX+DY*DY+DZ*DZ)
            CALL BDPCS3_UPDATED_BOND(Z(I),Z(J),R,DELTA,BO,RBD)
            WRITE(*,'(1X,I3,1X,I3,2(1X,I3),2(1X,A2),4(1X,F12.6))')
     &        I,J,Z(I),Z(J),ZSYM(Z(I)),ZSYM(Z(J)),R,DELTA,RBD,BO
   20    CONTINUE
   30 CONTINUE
      CALL BDPCS3_REPORT_HBONDS(N,Z,XYZ)
      STOP

  900 WRITE(*,*) 'Error: cannot open file.'
      STOP
  901 WRITE(*,*) 'Error: cannot read atom count.'
      STOP
  902 WRITE(*,*) 'Error: malformed XYZ file.'
      STOP
      END

c=======================================================================
c Public API without external topology bond-order override.
c=======================================================================
      SUBROUTINE BDPCS3_UPDATED_BOND(ZI,ZJ,R,DELTA,BONDORD,RBD)
      INTEGER ZI,ZJ
      DOUBLE PRECISION R,DELTA,BONDORD,RBD
      CALL BDPCS3_UPDATED_BOND_OVR(ZI,ZJ,R,0.0D0,0,
     &                             DELTA,BONDORD,RBD)
      RETURN
      END

c=======================================================================
c Public API with optional topology bond-order override.
c IOVR=0: use geometric bond order only.
c IOVR<>0: BONDORD=max(BOOVR, BO_GEOM), exactly as Python.
c=======================================================================
      SUBROUTINE BDPCS3_UPDATED_BOND_OVR(ZI,ZJ,R,BOOVR,IOVR,
     &                                   DELTA,BONDORD,RBD)
      INTEGER ZI,ZJ,IOVR,NI,NJ,NMAX
      DOUBLE PRECISION R,BOOVR,DELTA,BONDORD,RBD
      DOUBLE PRECISION VAL0,BOGEOM,FCOORD,DCV,DDELOC
      DOUBLE PRECISION SIGMA,RCOVBD,DERF_AS
      INTEGER NPRINC

      SIGMA=0.057D0
      VAL0=RCOVBD(ZI)+RCOVBD(ZJ)
      DELTA=0.0D0
      RBD=R
      BONDORD=0.0D0

      IF (VAL0 .LE. 0.0D0) RETURN
      BOGEOM=DEXP((VAL0-R)/0.30D0)
      BONDORD=BOGEOM
      IF (IOVR .NE. 0) THEN
         IF (BOOVR .GT. BONDORD) BONDORD=BOOVR
      ENDIF

      NI=NPRINC(ZI)
      NJ=NPRINC(ZJ)
      IF (NI .GT. 3) NI=3
      IF (NJ .GT. 3) NJ=3
      NMAX=NI
      IF (NJ .GT. NMAX) NMAX=NJ
      IF (NMAX .LE. 1) RETURN

      FCOORD=0.5D0*(1.0D0-DERF_AS((R-(1.3D0*VAL0))/SIGMA))
      DCV=-0.0025D0*DBLE(NMAX-1)
      CALL BDPCS3_DELOC(ZI,ZJ,R,DCV,DDELOC)
      DELTA=(DCV+DDELOC)*FCOORD
      RBD=R+DELTA
      RETURN
      END

c=======================================================================
c Hydrogen-bond BDPCS3 target for the H...Y distance.
c Donor is X in X-H...Y.  The angular gate is hard; the distance is
c smoothly damped by an error function so long contacts vanish.
c=======================================================================
      SUBROUTINE BDPCS3_HBOND_TARGET(ZD,ZA,RHY,ANGXHY,DELTA,RBD)
      INTEGER ZD,ZA
      DOUBLE PRECISION RHY,ANGXHY,DELTA,RBD
      DOUBLE PRECISION BASE,ARG,DAMP,DERF_AS,BDPCS3_HBOND_BASE
      INCLUDE 'bdpcs3_hbond_params.inc'
      DELTA=0.0D0
      RBD=RHY
      IF (ANGXHY .LT. BDPCS3_HB_ANGLE_MIN) RETURN
      BASE=BDPCS3_HBOND_BASE(ZD,ZA)
      IF (BASE .EQ. 0.0D0) RETURN
      ARG=(RHY-BDPCS3_HB_DIST_CUTOFF)/BDPCS3_HB_DIST_WIDTH
      DAMP=0.5D0*(1.0D0-DERF_AS(ARG))
      IF (DAMP .LT. 0.0D0) DAMP=0.0D0
      IF (DAMP .GT. 1.0D0) DAMP=1.0D0
      DELTA=BASE*DAMP
      RBD=RHY+DELTA
      RETURN
      END

      DOUBLE PRECISION FUNCTION BDPCS3_HBOND_BASE(ZD,ZA)
      INTEGER ZD,ZA
      INCLUDE 'bdpcs3_hbond_params.inc'
      BDPCS3_HBOND_BASE=0.0D0
      IF (ZD .EQ. 8 .AND. ZA .EQ. 8) THEN
         BDPCS3_HBOND_BASE=BDPCS3_HB_DELTA_OO
      ELSE IF (ZD .EQ. 7 .AND. ZA .EQ. 7) THEN
         BDPCS3_HBOND_BASE=BDPCS3_HB_DELTA_NN
      ELSE IF ((ZD .EQ. 7 .OR. ZD .EQ. 8) .AND.
     &         (ZA .EQ. 7 .OR. ZA .EQ. 8)) THEN
         BDPCS3_HBOND_BASE=0.5D0*
     &      (BDPCS3_HB_DELTA_OO+BDPCS3_HB_DELTA_NN)
      ENDIF
      RETURN
      END

      SUBROUTINE BDPCS3_REPORT_HBONDS(N,Z,XYZ)
      INTEGER N,Z(*)
      INTEGER MAXHB
      PARAMETER (MAXHB=2000)
      DOUBLE PRECISION XYZ(3,*)
      DOUBLE PRECISION RHY,ANG,DELTA,RBD,BESTR
      DOUBLE PRECISION DHD,DDA,DIST3,ANGLE_XHY
      INTEGER H,I,J,DON,ACC,IH,FOUND,NHB,ID(MAXHB),IA(MAXHB)
      LOGICAL COV_BOND,SEEN
      CHARACTER*2 ZSYM
      INCLUDE 'bdpcs3_hbond_params.inc'

      WRITE(*,'(/,A)') 'BDPCS3 hydrogen-bond targets'
      WRITE(*,'(A,F6.1,A,F5.2,A)')
     & 'Criteria: angle X-H-Y >= ',BDPCS3_HB_ANGLE_MIN,
     & ' deg; distance damping centered at ',BDPCS3_HB_DIST_CUTOFF,
     & ' A'
      WRITE(*,'(A)')
     & '  X   H   Y  ZX ZY       H...Y       Angle       Delta      Target'
      WRITE(*,'(A)')
     & '--------------------------------------------------------------------'
      NHB=0
      DO 40 H=1,N
         IF (Z(H) .NE. 1) GOTO 40
         DON=0
         BESTR=1.0D20
         DO 10 I=1,N
            IF (I .EQ. H) GOTO 10
            IF (Z(I) .NE. 7 .AND. Z(I) .NE. 8 .AND.
     &          Z(I) .NE. 16) GOTO 10
            DHD=DIST3(XYZ(1,H),XYZ(1,I))
            IF (.NOT. COV_BOND(Z(H),Z(I),DHD)) GOTO 10
            IF (DHD .LT. BESTR) THEN
               BESTR=DHD
               DON=I
            ENDIF
   10    CONTINUE
         IF (DON .EQ. 0) GOTO 40

         ACC=0
         BESTR=1.0D20
         DO 20 J=1,N
            IF (J .EQ. H .OR. J .EQ. DON) GOTO 20
            IF (Z(J) .NE. 7 .AND. Z(J) .NE. 8) GOTO 20
            DDA=DIST3(XYZ(1,DON),XYZ(1,J))
            IF (COV_BOND(Z(DON),Z(J),DDA)) GOTO 20
            RHY=DIST3(XYZ(1,H),XYZ(1,J))
            IF (RHY .GT. BDPCS3_HB_SEARCH_CUTOFF) GOTO 20
            ANG=ANGLE_XHY(XYZ(1,DON),XYZ(1,H),XYZ(1,J))
            IF (ANG .LT. BDPCS3_HB_ANGLE_MIN) GOTO 20
            IF (RHY .LT. BESTR) THEN
               BESTR=RHY
               ACC=J
            ENDIF
   20    CONTINUE
         IF (ACC .EQ. 0) GOTO 40
         SEEN=.FALSE.
         DO 30 IH=1,NHB
            IF (ID(IH) .EQ. DON .AND. IA(IH) .EQ. ACC) SEEN=.TRUE.
   30    CONTINUE
         IF (SEEN) GOTO 40
         IF (NHB .LT. MAXHB) THEN
            NHB=NHB+1
            ID(NHB)=DON
            IA(NHB)=ACC
         ENDIF
         RHY=DIST3(XYZ(1,H),XYZ(1,ACC))
         ANG=ANGLE_XHY(XYZ(1,DON),XYZ(1,H),XYZ(1,ACC))
         CALL BDPCS3_HBOND_TARGET(Z(DON),Z(ACC),RHY,ANG,DELTA,RBD)
         WRITE(*,'(3(1X,I3),2(1X,A2),4(1X,F11.6))')
     &      DON,H,ACC,ZSYM(Z(DON)),ZSYM(Z(ACC)),RHY,ANG,DELTA,RBD
   40 CONTINUE
      IF (NHB .EQ. 0) WRITE(*,'(A)') '  none'
      RETURN
      END

      DOUBLE PRECISION FUNCTION DIST3(A,B)
      DOUBLE PRECISION A(3),B(3),DX,DY,DZ
      DX=A(1)-B(1)
      DY=A(2)-B(2)
      DZ=A(3)-B(3)
      DIST3=DSQRT(DX*DX+DY*DY+DZ*DZ)
      RETURN
      END

      DOUBLE PRECISION FUNCTION ANGLE_XHY(X,H,Y)
      DOUBLE PRECISION X(3),H(3),Y(3),V1(3),V2(3),N1,N2,COSV
      INTEGER I
      N1=0.0D0
      N2=0.0D0
      COSV=0.0D0
      DO 10 I=1,3
         V1(I)=X(I)-H(I)
         V2(I)=Y(I)-H(I)
         N1=N1+V1(I)*V1(I)
         N2=N2+V2(I)*V2(I)
         COSV=COSV+V1(I)*V2(I)
   10 CONTINUE
      IF (N1 .LE. 1.0D-24 .OR. N2 .LE. 1.0D-24) THEN
         ANGLE_XHY=0.0D0
         RETURN
      ENDIF
      COSV=COSV/DSQRT(N1*N2)
      IF (COSV .GT. 1.0D0) COSV=1.0D0
      IF (COSV .LT. -1.0D0) COSV=-1.0D0
      ANGLE_XHY=DACOS(COSV)*180.0D0/(4.0D0*DATAN(1.0D0))
      RETURN
      END

      LOGICAL FUNCTION COV_BOND(ZI,ZJ,R)
      INTEGER ZI,ZJ
      DOUBLE PRECISION R,RCOVBD,RLIM
      RLIM=RCOVBD(ZI)+RCOVBD(ZJ)+0.45D0
      COV_BOND=(R .GT. 0.1D0 .AND. R .LE. RLIM)
      RETURN
      END

c=======================================================================
c C-C / C-S delocalization compensation.
c=======================================================================
      SUBROUTINE BDPCS3_DELOC(ZI,ZJ,R,DCV,DDELOC)
      INTEGER ZI,ZJ
      DOUBLE PRECISION R,DCV,DDELOC
      DOUBLE PRECISION R1S,R1D,R1T,R2S,R2D,R2T
      DOUBLE PRECISION RS,RD,RT,D1,D2,SIGMA,DMIN,XARG

      DDELOC=0.0D0
      IF (.NOT. ((ZI .EQ. 6 .AND. ZJ .EQ. 6) .OR.
     &    (ZI .EQ. 6 .AND. ZJ .EQ. 16) .OR.
     &    (ZI .EQ. 16 .AND. ZJ .EQ. 6))) RETURN

      CALL BDPCS3_PYY_REF(ZI,R1S,R1D,R1T)
      CALL BDPCS3_PYY_REF(ZJ,R2S,R2D,R2T)
      RS=R1S+R2S
      RD=R1D+R2D
      RT=R1T+R2T
      D1=DABS(RD-RS)
      D2=DABS(RT-RD)
      DMIN=D1
      IF (D2 .LT. DMIN) DMIN=D2
      IF (DMIN .GT. 0.0D0) THEN
         SIGMA=DMIN/1.5D0
         IF (SIGMA .LT. 1.0D-6) SIGMA=1.0D-6
      ELSE
         SIGMA=0.05D0
      ENDIF

      XARG=((R-RD)/SIGMA)**2
      IF (XARG .GT. 700.0D0) THEN
         DDELOC=0.0D0
      ELSE
         DDELOC=-DCV*DEXP(-XARG)
      ENDIF
      RETURN
      END

c=======================================================================
c Pyykko radii references exactly needed by the Python delocalization path.
c Carbon: single=coord3, double=coord2, triple=coord1.
c Sulfur: single=coord4, double=coord2, triple=coord6.
c=======================================================================
      SUBROUTINE BDPCS3_PYY_REF(Z,RSG,RDB,RTR)
      INTEGER Z
      DOUBLE PRECISION RSG,RDB,RTR,RCOVBD

      IF (Z .EQ. 6) THEN
         RSG=0.76D0
         RDB=0.67D0
         RTR=0.60D0
      ELSE IF (Z .EQ. 16) THEN
         RSG=1.05D0
         RDB=1.02D0
         RTR=1.09D0
      ELSE
         RSG=RCOVBD(Z)
         RDB=RSG
         RTR=RSG
      ENDIF
      RETURN
      END

c=======================================================================
c Principal quantum number as in Python _principal_quantum().
c=======================================================================
      INTEGER FUNCTION NPRINC(Z)
      INTEGER Z
      NPRINC=3
      IF (Z .LE. 10) NPRINC=2
      IF (Z .LE. 2) NPRINC=1
      RETURN
      END

c=======================================================================
c BDPCS3 updated covalent radius:
c explicit H/C/O/S values, then Mantina-Truhlar fallback table.
c=======================================================================
      DOUBLE PRECISION FUNCTION RCOVBD(Z)
      INTEGER Z
      DOUBLE PRECISION R
      R=0.0D0

c     Mantina-Truhlar fallback values from merlino_fit/topology/covalent_radii.py
      IF (Z .EQ.   1) R=0.32D0
      IF (Z .EQ.   2) R=0.37D0
      IF (Z .EQ.   3) R=1.30D0
      IF (Z .EQ.   4) R=0.99D0
      IF (Z .EQ.   5) R=0.84D0
      IF (Z .EQ.   6) R=0.75D0
      IF (Z .EQ.   7) R=0.71D0
      IF (Z .EQ.   8) R=0.64D0
      IF (Z .EQ.   9) R=0.60D0
      IF (Z .EQ.  10) R=0.62D0
      IF (Z .EQ.  11) R=1.60D0
      IF (Z .EQ.  12) R=1.40D0
      IF (Z .EQ.  13) R=1.24D0
      IF (Z .EQ.  14) R=1.14D0
      IF (Z .EQ.  15) R=1.09D0
      IF (Z .EQ.  16) R=1.04D0
      IF (Z .EQ.  17) R=1.00D0
      IF (Z .EQ.  18) R=1.01D0
      IF (Z .EQ.  19) R=2.00D0
      IF (Z .EQ.  20) R=1.74D0
      IF (Z .EQ.  21) R=1.59D0
      IF (Z .EQ.  22) R=1.48D0
      IF (Z .EQ.  23) R=1.44D0
      IF (Z .EQ.  24) R=1.30D0
      IF (Z .EQ.  25) R=1.29D0
      IF (Z .EQ.  26) R=1.24D0
      IF (Z .EQ.  27) R=1.18D0
      IF (Z .EQ.  28) R=1.17D0
      IF (Z .EQ.  29) R=1.22D0
      IF (Z .EQ.  30) R=1.20D0
      IF (Z .EQ.  31) R=1.23D0
      IF (Z .EQ.  32) R=1.20D0
      IF (Z .EQ.  33) R=1.20D0
      IF (Z .EQ.  34) R=1.18D0
      IF (Z .EQ.  35) R=1.17D0
      IF (Z .EQ.  36) R=1.16D0
      IF (Z .EQ.  37) R=2.15D0
      IF (Z .EQ.  38) R=1.90D0
      IF (Z .EQ.  39) R=1.76D0
      IF (Z .EQ.  40) R=1.64D0
      IF (Z .EQ.  41) R=1.56D0
      IF (Z .EQ.  42) R=1.46D0
      IF (Z .EQ.  43) R=1.38D0
      IF (Z .EQ.  44) R=1.36D0
      IF (Z .EQ.  45) R=1.34D0
      IF (Z .EQ.  46) R=1.30D0
      IF (Z .EQ.  47) R=1.36D0
      IF (Z .EQ.  48) R=1.40D0
      IF (Z .EQ.  49) R=1.42D0
      IF (Z .EQ.  50) R=1.40D0
      IF (Z .EQ.  51) R=1.40D0
      IF (Z .EQ.  52) R=1.37D0
      IF (Z .EQ.  53) R=1.36D0
      IF (Z .EQ.  54) R=1.36D0
      IF (Z .EQ.  55) R=2.38D0
      IF (Z .EQ.  56) R=2.06D0
      IF (Z .EQ.  57) R=1.94D0
      IF (Z .EQ.  58) R=1.84D0
      IF (Z .EQ.  59) R=1.90D0
      IF (Z .EQ.  60) R=1.88D0
      IF (Z .EQ.  61) R=1.86D0
      IF (Z .EQ.  62) R=1.85D0
      IF (Z .EQ.  63) R=1.83D0
      IF (Z .EQ.  64) R=1.82D0
      IF (Z .EQ.  65) R=1.81D0
      IF (Z .EQ.  66) R=1.80D0
      IF (Z .EQ.  67) R=1.79D0
      IF (Z .EQ.  68) R=1.77D0
      IF (Z .EQ.  69) R=1.77D0
      IF (Z .EQ.  70) R=1.78D0
      IF (Z .EQ.  71) R=1.74D0
      IF (Z .EQ.  72) R=1.64D0
      IF (Z .EQ.  73) R=1.58D0
      IF (Z .EQ.  74) R=1.50D0
      IF (Z .EQ.  75) R=1.41D0
      IF (Z .EQ.  76) R=1.36D0
      IF (Z .EQ.  77) R=1.32D0
      IF (Z .EQ.  78) R=1.30D0
      IF (Z .EQ.  79) R=1.30D0
      IF (Z .EQ.  80) R=1.32D0
      IF (Z .EQ.  81) R=1.44D0
      IF (Z .EQ.  82) R=1.45D0
      IF (Z .EQ.  83) R=1.50D0
      IF (Z .EQ.  84) R=1.42D0
      IF (Z .EQ.  85) R=1.48D0
      IF (Z .EQ.  86) R=1.46D0
      IF (Z .EQ.  87) R=2.42D0
      IF (Z .EQ.  88) R=2.11D0
      IF (Z .EQ.  89) R=2.01D0
      IF (Z .EQ.  90) R=1.90D0
      IF (Z .EQ.  91) R=1.84D0
      IF (Z .EQ.  92) R=1.83D0
      IF (Z .EQ.  93) R=1.80D0
      IF (Z .EQ.  94) R=1.80D0
      IF (Z .EQ.  95) R=1.73D0
      IF (Z .EQ.  96) R=1.68D0
      IF (Z .EQ.  97) R=1.68D0
      IF (Z .EQ.  98) R=1.68D0
      IF (Z .EQ.  99) R=1.65D0
      IF (Z .EQ. 100) R=1.67D0
      IF (Z .EQ. 101) R=1.73D0
      IF (Z .EQ. 102) R=1.76D0
      IF (Z .EQ. 103) R=1.61D0
      IF (Z .EQ. 104) R=1.57D0
      IF (Z .EQ. 105) R=1.49D0
      IF (Z .EQ. 106) R=1.43D0
      IF (Z .EQ. 107) R=1.41D0
      IF (Z .EQ. 108) R=1.34D0
      IF (Z .EQ. 109) R=1.29D0
      IF (Z .EQ. 110) R=1.28D0
      IF (Z .EQ. 111) R=1.21D0
      IF (Z .EQ. 112) R=1.22D0
      IF (Z .EQ. 113) R=1.36D0
      IF (Z .EQ. 114) R=1.43D0
      IF (Z .EQ. 115) R=1.62D0
      IF (Z .EQ. 116) R=1.75D0
      IF (Z .EQ. 117) R=1.65D0
      IF (Z .EQ. 118) R=1.57D0

c     Explicit BDPCS3 updated overrides from Python.
      IF (Z .EQ. 1) R=0.31D0
      IF (Z .EQ. 6) R=0.76D0
      IF (Z .EQ. 8) R=0.66D0
      IF (Z .EQ. 16) R=1.05D0

      RCOVBD=R
      RETURN
      END

c=======================================================================
c erf(x) approximation, Abramowitz-Stegun 7.1.26.
c=======================================================================
      DOUBLE PRECISION FUNCTION DERF_AS(X)
      DOUBLE PRECISION X,AX,T,POLY,EX2
      DOUBLE PRECISION P,A1,A2,A3,A4,A5,ONE
      ONE=1.0D0
      P=0.3275911D0
      A1=0.254829592D0
      A2=-0.284496736D0
      A3=1.421413741D0
      A4=-1.453152027D0
      A5=1.061405429D0
      AX=DABS(X)
      T=ONE/(ONE+P*AX)
      POLY=(((A5*T+A4)*T+A3)*T+A2)*T+A1
      IF (AX .GT. 26.0D0) THEN
         EX2=0.0D0
      ELSE
         EX2=DEXP(-AX*AX)
      ENDIF
      DERF_AS=DSIGN(ONE-POLY*T*EX2,X)
      RETURN
      END

c=======================================================================
c Element parsing utilities.
c=======================================================================
      INTEGER FUNCTION ATOMZ(S)
      CHARACTER*(*) S
      CHARACTER*2 U,UCASE2
      INTEGER IZ
      IZ=-1
      READ(S,*,ERR=10,END=10) IZ
      IF (IZ .GT. 0) THEN
         ATOMZ=IZ
         RETURN
      ENDIF
   10 CONTINUE
      U=UCASE2(S)
      ATOMZ=0
      IF (U .EQ. 'H ') ATOMZ=1
      IF (U .EQ. 'C ') ATOMZ=6
      IF (U .EQ. 'N ') ATOMZ=7
      IF (U .EQ. 'O ') ATOMZ=8
      IF (U .EQ. 'F ') ATOMZ=9
      IF (U .EQ. 'SI') ATOMZ=14
      IF (U .EQ. 'P ') ATOMZ=15
      IF (U .EQ. 'S ') ATOMZ=16
      IF (U .EQ. 'CL') ATOMZ=17
      IF (U .EQ. 'BR') ATOMZ=35
      IF (U .EQ. 'I ') ATOMZ=53
      RETURN
      END

      CHARACTER*2 FUNCTION ZSYM(Z)
      INTEGER Z
      ZSYM='??'
      IF (Z .EQ. 1) ZSYM='H '
      IF (Z .EQ. 6) ZSYM='C '
      IF (Z .EQ. 7) ZSYM='N '
      IF (Z .EQ. 8) ZSYM='O '
      IF (Z .EQ. 9) ZSYM='F '
      IF (Z .EQ. 14) ZSYM='Si'
      IF (Z .EQ. 15) ZSYM='P '
      IF (Z .EQ. 16) ZSYM='S '
      IF (Z .EQ. 17) ZSYM='Cl'
      IF (Z .EQ. 35) ZSYM='Br'
      IF (Z .EQ. 53) ZSYM='I '
      RETURN
      END

      CHARACTER*2 FUNCTION UCASE2(S)
      CHARACTER*(*) S
      CHARACTER*2 T
      INTEGER I
      T='  '
      DO 10 I=1,2
         IF (I .LE. LEN(S)) T(I:I)=S(I:I)
   10 CONTINUE
      DO 20 I=1,2
         IF (T(I:I) .GE. 'a' .AND. T(I:I) .LE. 'z') THEN
            T(I:I)=CHAR(ICHAR(T(I:I))-32)
         ENDIF
   20 CONTINUE
      UCASE2=T
      RETURN
      END
