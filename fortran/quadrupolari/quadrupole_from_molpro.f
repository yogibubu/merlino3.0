      PROGRAM PARSE_EFG_ONE
      IMPLICIT NONE

      DOUBLE PRECISION EXTRACT_LAST_NUMBER
      DOUBLE PRECISION QNUC

      CHARACTER*256 LINE
      CHARACTER*80 FNAME,OUTFILE
      INTEGER UNITIN,UNITOUT
      INTEGER FOUND
      DOUBLE PRECISION VXXi,VYYi,VZZi
      DOUBLE PRECISION ETA
      DOUBLE PRECISION CHIXX,CHIYY,CHIZZ
      DOUBLE PRECISION Q
      INTEGER Z

      WRITE(*,*) 'Nome file MOLPRO (.out):'
      READ(*,'(A)') FNAME

      WRITE(*,*) 'Numero atomico:'
      READ(*,*) Z

      OUTFILE='risultati_quadrupolo.dat'

      UNITIN=10
      UNITOUT=20

      OPEN(UNITIN,FILE=FNAME,STATUS='OLD')
      OPEN(UNITOUT,FILE=OUTFILE,STATUS='UNKNOWN')

      VXXi=0.D0
      VYYi=0.D0
      VZZi=0.D0
      FOUND=0

10    CONTINUE
      READ(UNITIN,'(A)',END=100) LINE
      CALL CLEAN(LINE)

      IF (INDEX(LINE,'(ED,relax)').GT.0) THEN

         IF (INDEX(LINE,'FGXX').GT.0) THEN
            VXXi = EXTRACT_LAST_NUMBER(LINE)
            FOUND = 1
         ENDIF

         IF (INDEX(LINE,'FGYY').GT.0) THEN
            VYYi = EXTRACT_LAST_NUMBER(LINE)
            FOUND = 1
         ENDIF

         IF (INDEX(LINE,'FGZZ').GT.0) THEN
            VZZi = EXTRACT_LAST_NUMBER(LINE)
            FOUND = 1
         ENDIF

      ENDIF

      GOTO 10

100   CONTINUE

      IF (FOUND.EQ.0) THEN
         WRITE(*,*) 'Nessun EFG trovato!'
         STOP
      ENDIF

C ---- calcolo eta
      CALL ORDER_EFG(VXXi, VYYi, VZZi, ETA)

C ---- conversione quadrupolare diretta
      Q = QNUC(Z)
      CHIXX = 234.9647D0 * Q * VXXi
      CHIYY = 234.9647D0 * Q * VYYi
      CHIZZ = 234.9647D0 * Q * VZZi

C ---- output
      WRITE(UNITOUT,*) 'Risultati quadrupolari (MHz)'
      WRITE(UNITOUT,*) '---------------------------'
      WRITE(UNITOUT,*) 'VXX = ',VXXi
      WRITE(UNITOUT,*) 'VYY = ',VYYi
      WRITE(UNITOUT,*) 'VZZ = ',VZZi
      WRITE(UNITOUT,*) ' '
      WRITE(UNITOUT,*) 'Eta = ',ETA
      WRITE(UNITOUT,*) 'ChiXX (MHz) = ',CHIXX
      WRITE(UNITOUT,*) 'ChiYY (MHz) = ',CHIYY
      WRITE(UNITOUT,*) 'ChiZZ (MHz) = ',CHIZZ

      CLOSE(UNITOUT)
      WRITE(*,*) 'File scritto:',OUTFILE

      STOP
      END


C ==============================================================

      SUBROUTINE CLEAN(L)
      IMPLICIT NONE
      CHARACTER*(*) L
      INTEGER I,K
      DO I=1,LEN(L)
         K=ICHAR(L(I:I))
         IF (K.LT.32 .OR. K.GT.126) L(I:I)=' '
      ENDDO
      RETURN
      END

C ==============================================================

      DOUBLE PRECISION FUNCTION EXTRACT_LAST_NUMBER(L)
      IMPLICIT NONE
      CHARACTER*(*) L
      CHARACTER*64 TMP
      INTEGER I,J
      DOUBLE PRECISION X

      EXTRACT_LAST_NUMBER=0.D0
      TMP=' '

      DO I=LEN(L),1,-1
         IF (L(I:I).NE.' ') GOTO 10
      ENDDO
      RETURN

10    CONTINUE
      J=I
      DO WHILE (I.GT.1 .AND. L(I-1:I-1).NE.' ')
         I=I-1
      ENDDO

      TMP=L(I:J)

      READ(TMP,*,ERR=20) X
      EXTRACT_LAST_NUMBER=X
      RETURN

20    CONTINUE
      EXTRACT_LAST_NUMBER=0.D0
      RETURN
      END

C ==============================================================

      DOUBLE PRECISION FUNCTION QNUC(Z)
      IMPLICIT NONE
      INTEGER Z

      QNUC=0.D0

      IF (Z.EQ.3)  QNUC=-4.00D0
      IF (Z.EQ.4)  QNUC=5.29D0
      IF (Z.EQ.7)  QNUC=0.02044D0
      IF (Z.EQ.8)  QNUC=-0.02578D0
      IF (Z.EQ.11) QNUC=10.4D0
      IF (Z.EQ.13) QNUC=14.66D0
      IF (Z.EQ.14) QNUC=-0.10D0
      IF (Z.EQ.15) QNUC=0.05D0
      IF (Z.EQ.17) QNUC=-0.08165D0
      IF (Z.EQ.19) QNUC=5.85D0
      IF (Z.EQ.20) QNUC=7.1D0
      IF (Z.EQ.23) QNUC=-5.2D0
      IF (Z.EQ.25) QNUC=0.33D0
      IF (Z.EQ.35) QNUC=0.33D0
      IF (Z.EQ.37) QNUC=13.2D0
      IF (Z.EQ.53) QNUC=-0.30D0
      IF (Z.EQ.55) QNUC=-0.0035D0

      RETURN
      END

      SUBROUTINE ORDER_EFG(VX, VY, VZ, ETA)
      IMPLICIT NONE
      DOUBLE PRECISION VX, VY, VZ
      DOUBLE PRECISION VMAX, VMED, VMIN, ETA
      DOUBLE PRECISION A(3), B(3)
      INTEGER I, J
      DOUBLE PRECISION TEMP

C --- Copia i valori in un array di appoggio
      A(1) = VX
      A(2) = VY
      A(3) = VZ

C --- Copia i valori assoluti in B
      B(1) = DABS(VX)
      B(2) = DABS(VY)
      B(3) = DABS(VZ)

C --- Ordina A e B per valore assoluto decrescente (bubble sort 3x3)
      DO I = 1, 2
        DO J = I+1, 3
          IF (B(J) .GT. B(I)) THEN
            TEMP = B(I)
            B(I) = B(J)
            B(J) = TEMP

            TEMP = A(I)
            A(I) = A(J)
            A(J) = TEMP
          ENDIF
        ENDDO
      ENDDO

C --- Dopo l’ordinamento:
C     A(1) = VMAX (|V| più grande)
C     A(2) = VMED (intermedio)
C     A(3) = VMIN (più piccolo)

      VMAX = A(1)
      VMED = A(2)
      VMIN = A(3)

C --- Calcola eta
      IF (VMAX .NE. 0.D0) THEN
         ETA = (VMIN - VMED) / VMAX
      ELSE
         ETA = 0.D0
      ENDIF

      RETURN
      END

