
C=======================================================================
C  STEP 3 - FIRME STRUTTURALI E GROUPING DELLE COORDINATE PRIMITIVE
C
C  Questo file contiene SOLO le routine di STEP 3:
C   - BUILD_SIGNATURES
C   - GROUP_PRIMITIVES
C
C  Dipendenze:
C   - ACLASS (classi atomiche)
C   - QTYPE, QSUB, QATOM, QCENTER (coordinate primitive)
C
C  Nessuna dipendenza da simboli atomici o gruppo puntuale.
C=======================================================================

C=======================================================================
C  BUILD_SIGNATURES
C
C  Costruisce le firme strutturali QSIG per le coordinate primitive.
C
C  QSIG(1) = QTYPE
C  QSIG(2) = QSUB
C  QSIG(3) = classe dell'atomo centrale (se esiste)
C  QSIG(4) = classe del primo atomo coinvolto
C  QSIG(5) = classe del terzo atomo coinvolto
C  QSIG(6) = classe del quarto atomo coinvolto
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

      INTEGER IQ, I1, I3, I4

      DO IQ = 1, NQ

         QSIG(1,IQ) = QTYPE(IQ)
         QSIG(2,IQ) = QSUB(IQ)

         IF (QCENTER(IQ) .GT. 0) THEN
            QSIG(3,IQ) = ACLASS(QCENTER(IQ))
         ELSE
            QSIG(3,IQ) = 0
         END IF

         I1 = QATOM(1,IQ)
         IF (I1 .GT. 0) THEN
            QSIG(4,IQ) = ACLASS(I1)
         ELSE
            QSIG(4,IQ) = 0
         END IF

         I3 = QATOM(3,IQ)
         IF (I3 .GT. 0) THEN
            QSIG(5,IQ) = ACLASS(I3)
         ELSE
            QSIG(5,IQ) = 0
         END IF

         I4 = QATOM(4,IQ)
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
C  Raggruppa le coordinate primitive combinabili in base alla firma QSIG.
C  Ogni gruppo contiene coordinate con firme IDENTICHE.
C
C  Le regole fisiche aggiuntive (X-H geminali, nodi policiclici, ecc.)
C  NON sono applicate qui e verranno gestite nello STEP successivo.
C=======================================================================
      SUBROUTINE GROUP_PRIMITIVES(
     &     MAXQ, MAXGRP, MAXORD,
     &     NQ, QSIG,
     &     NGRP, GRP_SIZE, GRP_Q)

      IMPLICIT NONE

      INTEGER MAXQ, MAXGRP, MAXORD
      INTEGER NQ, NGRP
      INTEGER QSIG(6,MAXQ)
      INTEGER GRP_SIZE(MAXGRP)
      INTEGER GRP_Q(MAXORD,MAXGRP)

      INTEGER IQ, JQ, K
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
               IF (QSIG(K,JQ) .NE. QSIG(K,IQ)) THEN
                  SAME = .FALSE.
               END IF
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
