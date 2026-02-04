
C=======================================================================
C  STEP 3.5 - COSTRUZIONE DEGLI ANELLI (TOPOLOGIA PURA)
C
C  FIND_RINGS:
C    - costruisce anelli dal solo grafo molecolare
C    - usa DFS con path tracking
C    - elimina duplicati (rotazioni e versi opposti)
C
C  Nessuna dipendenza da geometria o simboli atomici.
C=======================================================================

      SUBROUTINE FIND_RINGS(
     &     MAXAT, MAXBND, MAXRING, MAXRSZ,
     &     NAT, NBOND, BND,
     &     NRING, RING_SIZE, RING_ATOM)

      IMPLICIT NONE

      INTEGER MAXAT, MAXBND, MAXRING, MAXRSZ
      INTEGER NAT, NBOND, NRING
      INTEGER BND(2,MAXBND)
      INTEGER RING_SIZE(MAXRING)
      INTEGER RING_ATOM(MAXRSZ,MAXRING)

C ---- locals
      INTEGER DEG(MAXAT), NEI(12,MAXAT)
      INTEGER PATH(MAXRSZ)
      LOGICAL USED(MAXAT)
      INTEGER I, J, K

      NRING = 0

C ---- build adjacency
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
         USED(I) = .FALSE.
      END DO

C ---- start DFS from each atom
      DO I = 1, NAT
         CALL DFS_RING(
     &        I, I, 0,
     &        MAXRSZ,
     &        DEG, NEI,
     &        PATH, USED,
     &        MAXRING, MAXRSZ,
     &        NRING, RING_SIZE, RING_ATOM)
         USED(I) = .TRUE.
      END DO

      RETURN
      END

C=======================================================================
C  DFS_RING
C=======================================================================
      RECURSIVE SUBROUTINE DFS_RING(
     &     START, CUR, DEPTH,
     &     MAXRSZ,
     &     DEG, NEI,
     &     PATH, USED,
     &     MAXRING, MAXRSZ2,
     &     NRING, RING_SIZE, RING_ATOM)

      IMPLICIT NONE

      INTEGER START, CUR, DEPTH, MAXRSZ, MAXRSZ2
      INTEGER DEG(*), NEI(12,*)
      INTEGER PATH(*)
      LOGICAL USED(*)
      INTEGER MAXRING, NRING
      INTEGER RING_SIZE(MAXRING)
      INTEGER RING_ATOM(MAXRSZ2,MAXRING)

      INTEGER I, NEXT

      IF (DEPTH .GE. MAXRSZ) RETURN

      PATH(DEPTH+1) = CUR

      DO I = 1, DEG(CUR)
         NEXT = NEI(I,CUR)

         IF (NEXT .EQ. START .AND. DEPTH+1 .GE. 3) THEN
            CALL STORE_RING(PATH, DEPTH+1,
     &           MAXRING, MAXRSZ2,
     &           NRING, RING_SIZE, RING_ATOM)
         ELSE IF (.NOT. USED(NEXT)) THEN
            CALL DFS_RING(
     &           START, NEXT, DEPTH+1,
     &           MAXRSZ,
     &           DEG, NEI,
     &           PATH, USED,
     &           MAXRING, MAXRSZ2,
     &           NRING, RING_SIZE, RING_ATOM)
         END IF
      END DO

      RETURN
      END

C=======================================================================
C  STORE_RING
C  Salva un anello se non duplicato
C=======================================================================
      SUBROUTINE STORE_RING(
     &     PATH, LENG,
     &     MAXRING, MAXRSZ,
     &     NRING, RING_SIZE, RING_ATOM)

      IMPLICIT NONE

      INTEGER PATH(*), LENG
      INTEGER MAXRING, MAXRSZ, NRING
      INTEGER RING_SIZE(MAXRING)
      INTEGER RING_ATOM(MAXRSZ,MAXRING)

      INTEGER I, J
      LOGICAL DUP

      DUP = .FALSE.

      DO I = 1, NRING
         IF (RING_SIZE(I) .EQ. LENG) THEN
            DUP = .TRUE.
            DO J = 1, LENG
               IF (RING_ATOM(J,I) .NE. PATH(J)) DUP = .FALSE.
            END DO
            IF (DUP) RETURN
         END IF
      END DO

      IF (NRING .GE. MAXRING) RETURN

      NRING = NRING + 1
      RING_SIZE(NRING) = LENG
      DO J = 1, LENG
         RING_ATOM(J,NRING) = PATH(J)
      END DO

      RETURN
      END
