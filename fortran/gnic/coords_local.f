
C=======================================================================
C  NOMENCLATURA DELLE COORDINATE INTERNE E COLLETTIVE (CONGELATA)
C=======================================================================
C  Primitive:
C    QTYPE=1  Bond
C    QTYPE=2  Angle (non lineare)
C    QTYPE=3  Dihedral (torsione vera)
C    QTYPE=5  Linear bend (2 componenti)
C    QTYPE=6  Out-of-plane (Wilson)
C
C  Collettive (QTYPE=4): breathing, cyclic valence bending,
C  cyclic torsional, puckering, butterfly/hinge.
C=======================================================================

C=======================================================================
C  BUILD_ATOM_CLASSES
C=======================================================================
      SUBROUTINE BUILD_ATOM_CLASSES(
     &     MAXAT, MAXBND, MAXDEG,
     &     NAT, ZAT, NBOND, BND,
     &     ACLASS, NCLASS, IERR)

      IMPLICIT NONE
      INTEGER MAXAT, MAXBND, MAXDEG
      INTEGER NAT, NBOND, NCLASS, IERR
      INTEGER ZAT(MAXAT)
      INTEGER BND(2,MAXBND)
      INTEGER ACLASS(MAXAT)

      INTEGER I, J, K, IT, MAXITER
      INTEGER DEG(MAXAT), NEI(MAXDEG,MAXAT)
      INTEGER COL(MAXAT), COLNEW(MAXAT)
      INTEGER SIG(MAXDEG+1,MAXAT)
      INTEGER IDX(MAXAT)
      INTEGER CHANGED, NCUR

      IERR = 0
      MAXITER = 50

      DO I = 1, NAT
         DEG(I) = 0
      END DO

      DO K = 1, NBOND
         I = BND(1,K)
         J = BND(2,K)
         DEG(I) = DEG(I) + 1
         DEG(J) = DEG(J) + 1
         NEI(DEG(I),I) = J
         NEI(DEG(J),J) = I
      END DO

      DO I = 1, NAT
         COL(I) = 1000*ZAT(I) + DEG(I)
      END DO

      DO IT = 1, MAXITER
         DO I = 1, NAT
            SIG(1,I) = COL(I)
            DO K = 1, DEG(I)
               SIG(1+K,I) = COL(NEI(K,I))
            END DO
            CALL ISORT_INT(SIG(2,I), DEG(I))
         END DO

         DO I = 1, NAT
            IDX(I) = I
         END DO
         CALL SORT_IDX_BY_SIG(NAT, IDX, SIG, MAXDEG)

         NCUR = 0
         DO I = 1, NAT
            J = IDX(I)
            IF (I .EQ. 1) THEN
               NCUR = 1
            ELSE
               IF (SIG(1,J) .NE. SIG(1,IDX(I-1))) NCUR = NCUR + 1
            END IF
            COLNEW(J) = NCUR
         END DO

         CHANGED = 0
         DO I = 1, NAT
            IF (COLNEW(I) .NE. COL(I)) CHANGED = 1
            COL(I) = COLNEW(I)
         END DO
         IF (CHANGED .EQ. 0) EXIT
      END DO

      NCLASS = NCUR
      DO I = 1, NAT
         ACLASS(I) = COL(I)
      END DO
      RETURN
      END

C=======================================================================
C  BUILD_PRIMITIVE_COORDS
C=======================================================================
      SUBROUTINE BUILD_PRIMITIVE_COORDS(
     &     MAXAT, MAXBND, MAXQ, MAXDEG,
     &     NAT, ZAT, NBOND, BND,
     &     X, Y, ZC,
     &     NQ, QTYPE, QSUB, QATOM, QCENTER,
     &     IERR)

      IMPLICIT NONE
      INTEGER MAXAT, MAXBND, MAXQ, MAXDEG
      INTEGER NAT, NBOND, NQ, IERR
      INTEGER ZAT(MAXAT)
      INTEGER BND(2,MAXBND)
      DOUBLE PRECISION X(MAXAT), Y(MAXAT), ZC(MAXAT)

      INTEGER QTYPE(MAXQ), QSUB(MAXQ), QCENTER(MAXQ)
      INTEGER QATOM(4,MAXQ)

      INTEGER I, J, K, L, M
      INTEGER DEG(MAXAT), NEI(MAXDEG,MAXAT)
      DOUBLE PRECISION ANG, PI
      LOGICAL ISLINEAR

      PI = 3.141592653589793D0
      IERR = 0
      NQ = 0

C ---- adjacency
      DO I = 1, NAT
         DEG(I) = 0
      END DO

      DO K = 1, NBOND
         I = BND(1,K)
         J = BND(2,K)
         DEG(I) = DEG(I) + 1
         DEG(J) = DEG(J) + 1
         NEI(DEG(I),I) = J
         NEI(DEG(J),J) = I
      END DO

C ---- bonds
      DO K = 1, NBOND
         NQ = NQ + 1
         QTYPE(NQ) = 1
         QATOM(1,NQ) = BND(1,K)
         QATOM(2,NQ) = BND(2,K)
         QATOM(3,NQ) = 0
         QATOM(4,NQ) = 0
         QCENTER(NQ) = 0
         IF (ZAT(BND(1,K)) .EQ. 1 .OR. ZAT(BND(2,K)) .EQ. 1) THEN
            QSUB(NQ) = 1
         ELSE
            QSUB(NQ) = 2
         END IF
      END DO

C ---- angles / linear bends
      DO J = 1, NAT
         DO I = 1, DEG(J)-1
            DO K = I+1, DEG(J)
               CALL COMPUTE_ANGLE(
     &              NEI(I,J), J, NEI(K,J),
     &              X, Y, ZC, ANG)
               ISLINEAR = (ABS(ANG-PI) .LT. 1.0D-3)
               IF (.NOT. ISLINEAR) THEN
                  NQ = NQ + 1
                  QTYPE(NQ) = 2
                  QSUB(NQ) = 0
                  QATOM(1,NQ) = NEI(I,J)
                  QATOM(2,NQ) = J
                  QATOM(3,NQ) = NEI(K,J)
                  QATOM(4,NQ) = 0
                  QCENTER(NQ) = J
               ELSE
                  DO M = 1, 2
                     NQ = NQ + 1
                     QTYPE(NQ) = 5
                     QSUB(NQ) = M
                     QATOM(1,NQ) = NEI(I,J)
                     QATOM(2,NQ) = J
                     QATOM(3,NQ) = NEI(K,J)
                     QATOM(4,NQ) = 0
                     QCENTER(NQ) = J
                  END DO
               END IF
            END DO
         END DO
      END DO

C ---- out-of-plane (Wilson)
      DO J = 1, NAT
         IF (DEG(J) .GE. 3) THEN
            DO I = 1, DEG(J)
               DO K = 1, DEG(J)-1
                  DO L = K+1, DEG(J)
                     IF (I.NE.K .AND. I.NE.L) THEN
                        NQ = NQ + 1
                        QTYPE(NQ) = 6
                        QSUB(NQ) = 0
                        QATOM(1,NQ) = NEI(I,J)
                        QATOM(2,NQ) = J
                        QATOM(3,NQ) = NEI(K,J)
                        QATOM(4,NQ) = NEI(L,J)
                        QCENTER(NQ) = J
                     END IF
                  END DO
               END DO
            END DO
         END IF
      END DO

C ---- dihedrals
      DO K = 1, NBOND
         J = BND(1,K)
         L = BND(2,K)
         DO I = 1, DEG(J)
            IF (NEI(I,J) .NE. L) THEN
               DO M = 1, DEG(L)
                  IF (NEI(M,L) .NE. J) THEN
                     NQ = NQ + 1
                     QTYPE(NQ) = 3
                     QSUB(NQ) = 0
                     QATOM(1,NQ) = NEI(I,J)
                     QATOM(2,NQ) = J
                     QATOM(3,NQ) = L
                     QATOM(4,NQ) = NEI(M,L)
                     QCENTER(NQ) = 0
                  END IF
               END DO
            END IF
         END DO
      END DO

      RETURN
      END

C=======================================================================
C  UTILITIES
C=======================================================================
      SUBROUTINE ISORT_INT(A, N)
      IMPLICIT NONE
      INTEGER A(*), N
      INTEGER I, J, KEY
      DO I = 2, N
         KEY = A(I)
         J = I - 1
 10      IF (J .GE. 1 .AND. A(J) .GT. KEY) THEN
            A(J+1) = A(J)
            J = J - 1
            GO TO 10
         END IF
         A(J+1) = KEY
      END DO
      RETURN
      END

      SUBROUTINE SORT_IDX_BY_SIG(N, IDX, SIG, MAXDEG)
      IMPLICIT NONE
      INTEGER N, MAXDEG
      INTEGER IDX(*)
      INTEGER SIG(MAXDEG+1,*)
      INTEGER I, J, TMP
      DO I = 1, N-1
         DO J = 1, N-I
            IF (SIG(1,IDX(J)) .GT. SIG(1,IDX(J+1))) THEN
               TMP = IDX(J)
               IDX(J) = IDX(J+1)
               IDX(J+1) = TMP
            END IF
         END DO
      END DO
      RETURN
      END

C=======================================================================
C  COMPUTE_ANGLE
C=======================================================================
      SUBROUTINE COMPUTE_ANGLE(I, J, K, X, Y, ZC, ANG)
      IMPLICIT NONE
      INTEGER I, J, K
      DOUBLE PRECISION X(*), Y(*), ZC(*)
      DOUBLE PRECISION V1X, V1Y, V1Z, V2X, V2Y, V2Z
      DOUBLE PRECISION DOT, N1, N2
      DOUBLE PRECISION ANG

      V1X = X(I) - X(J)
      V1Y = Y(I) - Y(J)
      V1Z = ZC(I) - ZC(J)

      V2X = X(K) - X(J)
      V2Y = Y(K) - Y(J)
      V2Z = ZC(K) - ZC(J)

      DOT = V1X*V2X + V1Y*V2Y + V1Z*V2Z
      N1 = DSQRT(V1X*V1X + V1Y*V1Y + V1Z*V1Z)
      N2 = DSQRT(V2X*V2X + V2Y*V2Y + V2Z*V2Z)

      ANG = DACOS(DOT/(N1*N2))
      RETURN
      END

C=======================================================================
C  BUILD_SIGNATURES
C
C  Costruisce le firme strutturali QSIG per le coordinate primitive.
C=======================================================================
      SUBROUTINE BUILD_SIGNATURES(
     &     MAXQ,
     &     NQ, QTYPE, QSUB, QATOM, QCENTER,
     &     ACLASS,
     &     QSIG)

      IMPLICIT NONE
      INTEGER MAXQ, NQ
      INTEGER QTYPE(MAXQ), QSUB(MAXQ), QCENTER(MAXQ)
      INTEGER QATOM(4,MAXQ)
      INTEGER ACLASS(*)
      INTEGER QSIG(6,MAXQ)

      INTEGER IQ, I1, I2, I3, I4

      DO IQ = 1, NQ
         I1 = QATOM(1,IQ)
         I2 = QATOM(2,IQ)
         I3 = QATOM(3,IQ)
         I4 = QATOM(4,IQ)

         QSIG(1,IQ) = QTYPE(IQ)
         QSIG(2,IQ) = QSUB(IQ)

         IF (QCENTER(IQ) .GT. 0) THEN
            QSIG(3,IQ) = ACLASS(QCENTER(IQ))
         ELSE
            QSIG(3,IQ) = 0
         END IF

         IF (I1 .GT. 0) THEN
            QSIG(4,IQ) = ACLASS(I1)
         ELSE
            QSIG(4,IQ) = 0
         END IF

         IF (I3 .GT. 0) THEN
            QSIG(5,IQ) = ACLASS(I3)
         ELSE
            QSIG(5,IQ) = 0
         END IF

         IF (I4 .GT. 0) THEN
            QSIG(6,IQ) = ACLASS(I4)
         ELSE
            QSIG(6,IQ) = 0
         END IF
      END DO

      RETURN
      END

C=======================================================================
C  GROUP_PRIMITIVES
C
C  Raggruppa le coordinate primitive combinabili in base a QSIG.
C=======================================================================
      SUBROUTINE GROUP_PRIMITIVES(
     &     MAXQ, MAXGRP, MAXORD,
     &     NQ, QSIG, QTYPE, QSUB, QATOM,
     &     NGRP, GRP_SIZE, GRP_Q)

      IMPLICIT NONE
      INTEGER MAXQ, MAXGRP, MAXORD
      INTEGER NQ, NGRP
      INTEGER QSIG(6,MAXQ)
      INTEGER QTYPE(MAXQ), QSUB(MAXQ)
      INTEGER QATOM(4,MAXQ)
      INTEGER GRP_SIZE(MAXGRP)
      INTEGER GRP_Q(MAXORD,MAXGRP)

      INTEGER IQ, JQ, G, K
      LOGICAL USED(MAXQ)
      LOGICAL SAME

      DO IQ = 1, NQ
         USED(IQ) = .FALSE.
      END DO

      NGRP = 0

      DO IQ = 1, NQ
         IF (USED(IQ)) CYCLE

         NGRP = NGRP + 1
         GRP_SIZE(NGRP) = 1
         GRP_Q(1,NGRP) = IQ
         USED(IQ) = .TRUE.

         DO JQ = IQ+1, NQ
            IF (USED(JQ)) CYCLE

            SAME = .TRUE.
            DO K = 1, 6
               IF (QSIG(K,JQ) .NE. QSIG(K,IQ)) SAME = .FALSE.
            END DO

            IF (SAME) THEN
               GRP_SIZE(NGRP) = GRP_SIZE(NGRP) + 1
               GRP_Q(GRP_SIZE(NGRP),NGRP) = JQ
               USED(JQ) = .TRUE.
            END IF
         END DO
      END DO

      RETURN
      END
