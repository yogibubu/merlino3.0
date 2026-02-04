
C=======================================================================
C  STEP 3.6 - CANONICALIZZAZIONE E CLUSTERING DEGLI ANELLI
C
C  Contiene:
C   - CANONICALIZE_RING      : rende canonica la numerazione di un anello
C   - NORMALIZE_RINGS        : applica la canonizzazione a tutti gli anelli
C   - REMOVE_DUPLICATE_RINGS : elimina duplicati topologici
C   - BUILD_RING_ATOM_MAP    : costruisce mappa atomo-anello
C   - BUILD_RING_CLUSTERS    : cluster di anelli con atomi condivisi
C
C  Fortran77 rigido, solo topologia.
C=======================================================================

C=======================================================================
C  CANONICALIZE_RING
C=======================================================================
      SUBROUTINE CANONICALIZE_RING(LENG, RIN, ROUT)
      IMPLICIT NONE
      INTEGER LENG
      INTEGER RIN(*), ROUT(*)
      INTEGER TMP(100)
      INTEGER I, J, MAXRSZ, K, BEST

C --- inizializza con prima rotazione
      DO I = 1, LENG
         ROUT(I) = RIN(I)
      END DO

C --- rotazioni dirette
      DO K = 1, LENG
         DO I = 1, LENG
            J = I + K - 1
            IF (J .GT. LENG) J = J - LENG
            TMP(I) = RIN(J)
         END DO
         CALL CMP_AND_SET(LENG, TMP, ROUT)
      END DO

C --- rotazioni inverse
      DO K = 1, LENG
         DO I = 1, LENG
            J = K - I + 1
            IF (J .LE. 0) J = J + LENG
            TMP(I) = RIN(J)
         END DO
         CALL CMP_AND_SET(LENG, TMP, ROUT)
      END DO

      RETURN
      END

C=======================================================================
C  CMP_AND_SET
C  Sostituisce ROUT se TMP è lessicograficamente minore
C=======================================================================
      SUBROUTINE CMP_AND_SET(LENG, TMP, ROUT)
      IMPLICIT NONE
      INTEGER LENG
      INTEGER TMP(*), ROUT(*)
      INTEGER I, J

      DO I = 1, LENG
         IF (TMP(I) .LT. ROUT(I)) THEN
            DO J = 1, LENG
               ROUT(J) = TMP(J)
            END DO
            RETURN
         ELSE IF (TMP(I) .GT. ROUT(I)) THEN
            RETURN
         END IF
      END DO

      RETURN
      END

C=======================================================================
C  NORMALIZE_RINGS
C=======================================================================
      SUBROUTINE NORMALIZE_RINGS(
     &     MAXRING, MAXRSZ,
     &     NRING, RING_SIZE, RING_ATOM)

      IMPLICIT NONE
      INTEGER MAXRING, MAXRSZ, NRING
      INTEGER RING_SIZE(MAXRING)
      INTEGER RING_ATOM(MAXRSZ,MAXRING)

      INTEGER IR, TMP(100)

      DO IR = 1, NRING
         CALL CANONICALIZE_RING(
     &        RING_SIZE(IR),
     &        RING_ATOM(1,IR),
     &        TMP)
         CALL COPY_RING(RING_SIZE(IR), TMP, RING_ATOM(1,IR))
      END DO

      RETURN
      END

C=======================================================================
C  COPY_RING
C=======================================================================
      SUBROUTINE COPY_RING(LENG, SRC, DST)
      IMPLICIT NONE
      INTEGER LENG
      INTEGER SRC(*), DST(*)
      INTEGER I
      DO I = 1, LENG
         DST(I) = SRC(I)
      END DO
      RETURN
      END

C=======================================================================
C  REMOVE_DUPLICATE_RINGS
C=======================================================================
      SUBROUTINE REMOVE_DUPLICATE_RINGS(
     &     MAXRING, MAXRSZ,
     &     NRING, RING_SIZE, RING_ATOM)

      IMPLICIT NONE
      INTEGER MAXRING, MAXRSZ, NRING
      INTEGER RING_SIZE(MAXRING)
      INTEGER RING_ATOM(MAXRSZ,MAXRING)

      INTEGER IAT, I, J, K, NRNEW
      LOGICAL DUP

      NRNEW = 0

      DO I = 1, NRING
         DUP = .FALSE.
         DO J = 1, NRNEW
            IF (RING_SIZE(I) .EQ. RING_SIZE(J)) THEN
               DUP = .TRUE.
               DO K = 1, RING_SIZE(I)
                  IF (RING_ATOM(K,I) .NE. RING_ATOM(K,J)) DUP = .FALSE.
               END DO
               IF (DUP) EXIT
            END IF
         END DO

         IF (.NOT. DUP) THEN
            NRNEW = NRNEW + 1
            RING_SIZE(NRNEW) = RING_SIZE(I)
            DO K = 1, RING_SIZE(I)
               RING_ATOM(K,NRNEW) = RING_ATOM(K,I)
            END DO
         END IF
      END DO

      NRING = NRNEW
      RETURN
      END

C=======================================================================
C  BUILD_RING_ATOM_MAP
C=======================================================================
      SUBROUTINE BUILD_RING_ATOM_MAP(
     &     MAXAT, MAXRING, MAXRSZ,
     &     NAT, NRING, RING_SIZE, RING_ATOM,
     &     A2R, NA2R)

      IMPLICIT NONE
      INTEGER MAXAT, MAXRING, MAXRSZ
      INTEGER NAT, NRING
      INTEGER RING_SIZE(MAXRING)
      INTEGER RING_ATOM(MAXRSZ,MAXRING)
      INTEGER A2R(MAXAT,MAXRING)
      INTEGER NA2R(MAXAT)

      INTEGER IAT, J, K

      DO IAT = 1, NAT
         NA2R(IAT) = 0
      END DO

      DO J = 1, NRING
         DO K = 1, RING_SIZE(J)
            IAT = RING_ATOM(K,J)
            NA2R(IAT) = NA2R(IAT) + 1
            A2R(IAT,NA2R(IAT)) = J
         END DO
      END DO

      RETURN
      END

C=======================================================================
C  BUILD_RING_CLUSTERS
C=======================================================================
      SUBROUTINE BUILD_RING_CLUSTERS(
     &     MAXRING, MAXRSZ,
     &     NRING, RING_SIZE, RING_ATOM,
     &     NCLUST, CLUST_SIZE, CLUST_RING)

      IMPLICIT NONE
      INTEGER MAXRING, MAXRSZ
      INTEGER NRING, NCLUST
      INTEGER RING_SIZE(MAXRING)
      INTEGER RING_ATOM(MAXRSZ,MAXRING)
      INTEGER CLUST_SIZE(MAXRING)
      INTEGER CLUST_RING(MAXRING,MAXRING)

      LOGICAL USED(MAXRING), SHARE_ATOM
      INTEGER I, J

      DO I = 1, NRING
         USED(I) = .FALSE.
      END DO

      NCLUST = 0

      DO I = 1, NRING
         IF (USED(I)) CYCLE
         NCLUST = NCLUST + 1
         CLUST_SIZE(NCLUST) = 1
         CLUST_RING(1,NCLUST) = I
         USED(I) = .TRUE.

         DO J = I+1, NRING
            IF (.NOT. USED(J)) THEN
               IF (SHARE_ATOM(I, J, MAXRSZ, RING_SIZE, RING_ATOM)) THEN
                  CLUST_SIZE(NCLUST) = CLUST_SIZE(NCLUST) + 1
                  CLUST_RING(CLUST_SIZE(NCLUST),NCLUST) = J
                  USED(J) = .TRUE.
               END IF
            END IF
         END DO
      END DO

      RETURN
      END

C=======================================================================
C  SHARE_ATOM
C=======================================================================
      LOGICAL FUNCTION SHARE_ATOM(I, J, MAXRSZ, RING_SIZE, RING_ATOM)
      IMPLICIT NONE
      INTEGER I, J, MAXRSZ
      INTEGER RING_SIZE(*)
      INTEGER RING_ATOM(MAXRSZ,*)
      INTEGER A, B, K, L

      SHARE_ATOM = .FALSE.
      DO K = 1, RING_SIZE(I)
         A = RING_ATOM(K,I)
         DO L = 1, RING_SIZE(J)
            B = RING_ATOM(L,J)
            IF (A .EQ. B) THEN
               SHARE_ATOM = .TRUE.
               RETURN
            END IF
         END DO
      END DO

      RETURN
      END
