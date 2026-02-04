
C=======================================================================
C  STEP 4 - COSTRUZIONE DELLE COORDINATE COLLETTIVE
C
C  Convenzione:
C   - Le combinazioni sono ORTOGONALI ma NON normalizzate
C   - I fattori di normalizzazione sono restituiti separatamente
C
C  Dipendenze:
C   - GROUP_PRIMITIVES (STEP 3)
C   - RING_ATOM, CLUST_RING (STEP 3.5-3.6)
C=======================================================================

C=======================================================================
C  BUILD_COLLECTIVE_FROM_GROUP
C
C  Costruisce combinazioni ortogonali:
C   - 1 coordinata media
C   - (N-1) deviazioni indipendenti
C=======================================================================
      SUBROUTINE BUILD_COLLECTIVE_FROM_GROUP(
     &     MAXORD,
     &     NG, GRP_Q, GRP_SIZE,
     &     NLC, LC_COEFF, LC_NORM)

      IMPLICIT NONE

      INTEGER MAXORD
      INTEGER NG, GRP_SIZE
      INTEGER GRP_Q(MAXORD)

      INTEGER NLC
      DOUBLE PRECISION LC_COEFF(MAXORD,MAXORD)
      DOUBLE PRECISION LC_NORM(MAXORD)

      INTEGER I, J

C --- numero di coordinate collettive
      NLC = GRP_SIZE

C --- media
      DO I = 1, GRP_SIZE
         LC_COEFF(1,I) = 1.0D0
      END DO
      LC_NORM(1) = DBLE(GRP_SIZE)

C --- deviazioni ortogonali
      DO I = 2, GRP_SIZE
         DO J = 1, GRP_SIZE
            LC_COEFF(I,J) = 0.0D0
         END DO
         LC_COEFF(I,1) = 1.0D0
         LC_COEFF(I,I) = -1.0D0
         LC_NORM(I) = 2.0D0
      END DO

      RETURN
      END

C=======================================================================
C  BUILD_RING_COLLECTIVES
C
C  Costruisce:
C   - ring breathing
C   - cyclic valence bending (CVB)
C   - cyclic torsional
C=======================================================================
      SUBROUTINE BUILD_RING_COLLECTIVES(
     &     MAXRSZ,
     &     RING_SIZE, RING_ATOM,
     &     NLC, LC_COEFF, LC_NORM)

      IMPLICIT NONE

      INTEGER MAXRSZ
      INTEGER RING_SIZE
      INTEGER RING_ATOM(MAXRSZ)

      INTEGER NLC
      DOUBLE PRECISION LC_COEFF(MAXRSZ)
      DOUBLE PRECISION LC_NORM

      INTEGER I

C --- ring breathing (somma)
      NLC = 1
      DO I = 1, RING_SIZE
         LC_COEFF(I) = 1.0D0
      END DO
      LC_NORM = DBLE(RING_SIZE)

      RETURN
      END

C=======================================================================
C  REDUCE_CLUSTER_DOF
C
C  Elimina una combinazione ridondante per cluster (regola N-1)
C=======================================================================
      SUBROUTINE REDUCE_CLUSTER_DOF(
     &     MAXRING,
     &     CLUST_SIZE, CLUST_RING,
     &     USE_RING)

      IMPLICIT NONE

      INTEGER MAXRING
      INTEGER CLUST_SIZE
      INTEGER CLUST_RING(MAXRING)
      LOGICAL USE_RING(MAXRING)

      INTEGER I

C --- elimina la prima per convenzione
      USE_RING(CLUST_RING(1)) = .FALSE.
      DO I = 2, CLUST_SIZE
         USE_RING(CLUST_RING(I)) = .TRUE.
      END DO

      RETURN
      END
