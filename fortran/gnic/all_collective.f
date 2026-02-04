C=======================================================================
C  BUILD_ALL_COLLECTIVE_COORDS
C
C  Costruttore generale delle coordinate collettive:
C   - gruppi locali di primitive
C   - coordinate cicliche di anello (breathing, CVB, CT)
C   - riduzione N-1 sui cluster di anelli
C
C  Fortran77 / old-style F90 compliant
C=======================================================================

      SUBROUTINE BUILD_ALL_COLLECTIVE_COORDS(
     &     MAXQ, MAXORD, MAXCOL, MAXRSZ, MAXRING,
     &     NQ, QTYPE, QSUB, QATOM, QCENTER,
     &     NGRP, GRP_SIZE, GRP_Q,
     &     NRING, RING_SIZE, RING_ATOM,
     &     NCLUST, CLUST_SIZE, CLUST_RING,
     &     NCOL, COL_NTERM, COL_Q, COL_COEFF, COL_NORM)

      IMPLICIT NONE

C-----------------------------------------------------------------------
C Arguments
C-----------------------------------------------------------------------
      INTEGER MAXQ, MAXORD, MAXCOL, MAXRSZ, MAXRING
      INTEGER NQ, NGRP, NRING, NCLUST, NCOL

      INTEGER QTYPE(MAXQ), QSUB(MAXQ), QCENTER(MAXQ)
      INTEGER QATOM(4,MAXQ)

      INTEGER GRP_SIZE(*)
      INTEGER GRP_Q(MAXORD,*)

      INTEGER RING_SIZE(*)
      INTEGER RING_ATOM(MAXRSZ,*)

      INTEGER CLUST_SIZE(*)
      INTEGER CLUST_RING(MAXRING,*)

      INTEGER COL_NTERM(MAXCOL)
      INTEGER COL_Q(MAXORD,MAXCOL)
      DOUBLE PRECISION COL_COEFF(MAXORD,MAXCOL)
      DOUBLE PRECISION COL_NORM(MAXCOL)

C-----------------------------------------------------------------------
C Locals
C-----------------------------------------------------------------------
      INTEGER IG, IR, I, K
      INTEGER NTERM
      LOGICAL USE_RING(MAXRING)

      INTEGER TMP_Q(MAXORD)
      DOUBLE PRECISION TMP_COEFF(MAXORD)
      DOUBLE PRECISION TMP_NORM

C=======================================================================
C Initialization
C=======================================================================
      NCOL = 0

C=======================================================================
C 1) COLLECTIVE COORDINATES FROM LOCAL GROUPS
C=======================================================================
      DO IG = 1, NGRP

         CALL BUILD_COLLECTIVE_FROM_GROUP(
     &      MAXORD,
     &      IG, GRP_Q(1,IG), GRP_SIZE(IG),
     &      NTERM, TMP_COEFF, TMP_NORM)

         IF (NTERM .GT. 0) THEN
            NCOL = NCOL + 1
            COL_NTERM(NCOL) = NTERM

            DO I = 1, NTERM
               COL_Q(I,NCOL)     = GRP_Q(I,IG)
               COL_COEFF(I,NCOL)= TMP_COEFF(I)
            END DO

            COL_NORM(NCOL) = TMP_NORM
         END IF

      END DO

C=======================================================================
C 2) REDUCTION N-1 ON RING CLUSTERS
C=======================================================================
      CALL REDUCE_RING_CLUSTER_N1(
     &     MAXRING,
     &     NCLUST, CLUST_SIZE, CLUST_RING,
     &     USE_RING)

C=======================================================================
C 3) RING COLLECTIVE COORDINATES
C=======================================================================
      DO IR = 1, NRING

         IF (.NOT. USE_RING(IR)) GOTO 300

C-----------------------------------------------------------------------
C 3a) Ring breathing
C-----------------------------------------------------------------------
         NCOL = NCOL + 1
         COL_NTERM(NCOL) = RING_SIZE(IR)

         DO I = 1, RING_SIZE(IR)
            COL_Q(I,NCOL)      = RING_ATOM(I,IR)
            COL_COEFF(I,NCOL) = 1.0D0
         END DO

         COL_NORM(NCOL) = DBLE(RING_SIZE(IR))

C-----------------------------------------------------------------------
C 3b) Cyclic Valence Bending (CVB)
C-----------------------------------------------------------------------
         CALL BUILD_CVB_FOR_RING(
     &      MAXQ,
     &      RING_SIZE(IR), RING_ATOM(1,IR),
     &      NQ, QTYPE, QATOM, QCENTER,
     &      NTERM, TMP_Q, TMP_COEFF, TMP_NORM)

         IF (NTERM .GT. 0) THEN
            NCOL = NCOL + 1
            COL_NTERM(NCOL) = NTERM

            DO I = 1, NTERM
               COL_Q(I,NCOL)      = TMP_Q(I)
               COL_COEFF(I,NCOL) = TMP_COEFF(I)
            END DO

            COL_NORM(NCOL) = TMP_NORM
         END IF

C-----------------------------------------------------------------------
C 3c) Cyclic torsion
C-----------------------------------------------------------------------
         CALL BUILD_CT_FOR_RING(
     &      MAXQ,
     &      RING_SIZE(IR), RING_ATOM(1,IR),
     &      NQ, QTYPE, QATOM,
     &      NTERM, TMP_Q, TMP_COEFF, TMP_NORM)

         IF (NTERM .GT. 0) THEN
            NCOL = NCOL + 1
            COL_NTERM(NCOL) = NTERM

            DO I = 1, NTERM
               COL_Q(I,NCOL)      = TMP_Q(I)
               COL_COEFF(I,NCOL) = TMP_COEFF(I)
            END DO

            COL_NORM(NCOL) = TMP_NORM
         END IF

 300     CONTINUE
      END DO

      RETURN
      END

