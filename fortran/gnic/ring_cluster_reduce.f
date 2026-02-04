
C=======================================================================
C  STEP 4c - RIDUZIONE N-1 PER COORDINATE CICLICHE DI ANELLO
C
C  Applica la regola:
C     cluster con M anelli -> M-1 coordinate indipendenti
C
C  La routine NON ricostruisce le coordinate:
C   - seleziona quali coordinate usare/scartare
C   - lascia invariati coefficienti e norme
C
C  Convenzione:
C   - per ogni cluster, si scarta l'anello con indice più piccolo
C   - deterministica, riproducibile
C=======================================================================

C=======================================================================
C  REDUCE_RING_CLUSTER_N1
C
C  Apply N-1 reduction to ring clusters:
C     cluster with M rings -> M-1 independent coordinates
C
C  Strategy (deterministic):
C   - for each cluster, discard the ring with the smallest index
C
C  Fortran77 / old-style F90 compliant
C=======================================================================

      SUBROUTINE REDUCE_RING_CLUSTER_N1(
     &     MAXRING,
     &     NCLUST, CLUST_SIZE, CLUST_RING,
     &     USE_RING)

      IMPLICIT NONE

C-----------------------------------------------------------------------
C Arguments
C-----------------------------------------------------------------------
      INTEGER MAXRING
      INTEGER NCLUST
      INTEGER CLUST_SIZE(*)
      INTEGER CLUST_RING(MAXRING,*)
      LOGICAL USE_RING(MAXRING)

C-----------------------------------------------------------------------
C Locals
C-----------------------------------------------------------------------
      INTEGER IR, IC, K
      INTEGER IRMIN

C=======================================================================
C Initialization: all rings usable
C=======================================================================
      DO IR = 1, MAXRING
         USE_RING(IR) = .TRUE.
      END DO

C=======================================================================
C Apply N-1 reduction cluster by cluster
C=======================================================================
      DO IC = 1, NCLUST

         IF (CLUST_SIZE(IC) .LE. 1) GOTO 100

C ---- find smallest ring index in this cluster
         IRMIN = CLUST_RING(1,IC)
         DO K = 2, CLUST_SIZE(IC)
            IF (CLUST_RING(K,IC) .LT. IRMIN) THEN
               IRMIN = CLUST_RING(K,IC)
            END IF
         END DO

C ---- discard that ring
         USE_RING(IRMIN) = .FALSE.

 100     CONTINUE
      END DO

      RETURN
      END

