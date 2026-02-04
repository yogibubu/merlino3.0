
C=======================================================================
C  STEP 4a - CYCLIC VALENCE BENDING (CVB)
C
C  Definizione:
C   - CVB = combinazione ciclica orientata degli angoli di valenza
C     lungo un anello
C   - Analogo alle torsioni cicliche
C
C  Proprietà:
C   - 1 coordinata per anello
C   - ortogonale alle combinazioni locali
C   - non normalizzata (norma restituita separatamente)
C
C  Input:
C   - RING_ATOM, RING_SIZE
C   - lista globale degli angoli di valenza (primitive)
C
C  Output:
C   - coefficienti CVB sugli angoli primitivi
C=======================================================================

C=======================================================================
C  BUILD_CVB_FOR_RING
C
C  Cyclic Valence Bending (CVB) for a single ring
C
C  Fortran77 / old-style F90 compliant
C=======================================================================

      SUBROUTINE BUILD_CVB_FOR_RING(
     &     MAXQ,
     &     RING_SIZE, RING_ATOM,
     &     NQ, QTYPE, QATOM, QCENTER,
     &     NTERM, CVB_Q, CVB_COEFF, CVB_NORM)

      IMPLICIT NONE

C-----------------------------------------------------------------------
C Arguments
C-----------------------------------------------------------------------
      INTEGER MAXQ
      INTEGER RING_SIZE
      INTEGER RING_ATOM(*)

      INTEGER NQ
      INTEGER QTYPE(MAXQ), QCENTER(MAXQ)
      INTEGER QATOM(4,MAXQ)

      INTEGER NTERM
      INTEGER CVB_Q(MAXQ)
      DOUBLE PRECISION CVB_COEFF(MAXQ)
      DOUBLE PRECISION CVB_NORM

C-----------------------------------------------------------------------
C Locals
C-----------------------------------------------------------------------
      INTEGER I, IQ
      INTEGER A, B, C
      INTEGER IA, IC

C=======================================================================
C Initialization
C=======================================================================
      NTERM    = 0
      CVB_NORM = 0.0D0

C=======================================================================
C Loop over ring atoms as central atoms
C=======================================================================
      DO I = 1, RING_SIZE

         B = RING_ATOM(I)

         IF (I .LT. RING_SIZE) THEN
            A = RING_ATOM(I)
            C = RING_ATOM(I+1)
         ELSE
            A = RING_ATOM(I)
            C = RING_ATOM(1)
         END IF

C-----------------------------------------------------------------------
C Find matching valence angle A-B-C
C-----------------------------------------------------------------------
         DO IQ = 1, NQ

            IF (QTYPE(IQ) .EQ. 2) THEN
               IF (QCENTER(IQ) .EQ. B) THEN

                  IA = QATOM(1,IQ)
                  IC = QATOM(3,IQ)

                  IF ((IA .EQ. A .AND. IC .EQ. C) .OR.
     &                (IA .EQ. C .AND. IC .EQ. A)) THEN

                     NTERM = NTERM + 1
                     CVB_Q(NTERM)      = IQ
                     CVB_COEFF(NTERM)  = 1.0D0
                     CVB_NORM          = CVB_NORM + 1.0D0
                     GOTO 100

                  END IF
               END IF
            END IF

         END DO

 100     CONTINUE
      END DO

      RETURN
      END

