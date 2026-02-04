
C=======================================================================
C  STEP 4b - CYCLIC TORSIONAL COORDINATE (CT)
C
C  Definizione:
C   - CT = combinazione ciclica delle torsioni lungo un anello
C   - Analogo diretto del CVB per le torsioni
C
C  Proprietà:
C   - 1 coordinata per anello
C   - ortogonale alle torsioni locali
C   - non normalizzata (norma restituita separatamente)
C
C  Input:
C   - RING_ATOM, RING_SIZE
C   - lista globale delle torsioni primitive (QTYPE = 3)
C
C  Output:
C   - coefficienti CT sulle torsioni primitive
C=======================================================================

C=======================================================================
C  BUILD_CT_FOR_RING
C
C  Cyclic torsional coordinate for a single ring
C
C  Fortran77 / old-style F90 compliant
C=======================================================================

      SUBROUTINE BUILD_CT_FOR_RING(
     &     MAXQ,
     &     RING_SIZE, RING_ATOM,
     &     NQ, QTYPE, QATOM,
     &     NTERM, CT_Q, CT_COEFF, CT_NORM)

      IMPLICIT NONE

C-----------------------------------------------------------------------
C Arguments
C-----------------------------------------------------------------------
      INTEGER MAXQ
      INTEGER RING_SIZE
      INTEGER RING_ATOM(*)

      INTEGER NQ
      INTEGER QTYPE(MAXQ)
      INTEGER QATOM(4,MAXQ)

      INTEGER NTERM
      INTEGER CT_Q(MAXQ)
      DOUBLE PRECISION CT_COEFF(MAXQ)
      DOUBLE PRECISION CT_NORM

C-----------------------------------------------------------------------
C Locals
C-----------------------------------------------------------------------
      INTEGER I, IQ
      INTEGER A, B, C, D
      INTEGER I1, I2, I3, I4

C=======================================================================
C Initialization
C=======================================================================
      NTERM   = 0
      CT_NORM = 0.0D0

C=======================================================================
C Loop over cyclic sequences of four atoms in the ring
C=======================================================================
      DO I = 1, RING_SIZE

C ---- determine A-B-C-D with cyclic indexing
         IF (I .LE. RING_SIZE-3) THEN
            A = RING_ATOM(I)
            B = RING_ATOM(I+1)
            C = RING_ATOM(I+2)
            D = RING_ATOM(I+3)

         ELSE IF (I .EQ. RING_SIZE-2) THEN
            A = RING_ATOM(I)
            B = RING_ATOM(I+1)
            C = RING_ATOM(I+2)
            D = RING_ATOM(1)

         ELSE IF (I .EQ. RING_SIZE-1) THEN
            A = RING_ATOM(I)
            B = RING_ATOM(I+1)
            C = RING_ATOM(1)
            D = RING_ATOM(2)

         ELSE
            A = RING_ATOM(I)
            B = RING_ATOM(1)
            C = RING_ATOM(2)
            D = RING_ATOM(3)
         END IF

C-----------------------------------------------------------------------
C Find matching torsion primitive A-B-C-D
C-----------------------------------------------------------------------
         DO IQ = 1, NQ

            IF (QTYPE(IQ) .EQ. 3) THEN

               I1 = QATOM(1,IQ)
               I2 = QATOM(2,IQ)
               I3 = QATOM(3,IQ)
               I4 = QATOM(4,IQ)

               IF ((I1 .EQ. A .AND. I2 .EQ. B .AND.
     &              I3 .EQ. C .AND. I4 .EQ. D) .OR.
     &             (I1 .EQ. D .AND. I2 .EQ. C .AND.
     &              I3 .EQ. B .AND. I4 .EQ. A)) THEN

                  NTERM = NTERM + 1
                  CT_Q(NTERM)     = IQ
                  CT_COEFF(NTERM) = 1.0D0
                  CT_NORM         = CT_NORM + 1.0D0
                  GOTO 100

               END IF
            END IF

         END DO

 100     CONTINUE
      END DO

      RETURN
      END

