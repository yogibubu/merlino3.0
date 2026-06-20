C=======================================================================
C  TEST_READER
C
C  Test program for READ_XYZIN.
C  Reads Merlino enriched XYZ file "xyzin" and prints parsed data.
C=======================================================================

      PROGRAM TEST_READER
      IMPLICIT NONE

C-----------------------------------------------------------------------
C Parameters
C-----------------------------------------------------------------------
      INTEGER MAXAT
      PARAMETER (MAXAT = 200)

C-----------------------------------------------------------------------
C Variables
C-----------------------------------------------------------------------
      INTEGER NAT, I, IERR
      INTEGER CHARGE, MULT, SIGMA

      CHARACTER*2 SYM(MAXAT)
      CHARACTER*8 PG

      DOUBLE PRECISION XYZ(3,MAXAT)
      DOUBLE PRECISION A, B, C, KAPPA
      DOUBLE PRECISION T, QTRANS, QROT_ISO, QROT_AVG

C-----------------------------------------------------------------------
C Call reader
C-----------------------------------------------------------------------
      CALL READ_XYZIN(
     &   NAT, SYM, XYZ,
     &   CHARGE, MULT, PG,
     &   A, B, C, KAPPA,
     &   SIGMA, T,
     &   QTRANS, QROT_ISO, QROT_AVG,
     &   IERR
     & )

C-----------------------------------------------------------------------
C Error handling
C-----------------------------------------------------------------------
      IF (IERR .NE. 0) THEN
         WRITE(*,*) 'ERROR: READ_XYZIN failed. Code = ', IERR
         IF (IERR .EQ. 1) THEN
            WRITE(*,*) 'Reason: file "xyzin" not found.'
         ELSE IF (IERR .EQ. 2) THEN
            WRITE(*,*) 'Reason: error reading XYZ header.'
         ELSE IF (IERR .EQ. 3) THEN
            WRITE(*,*) 'Reason: error reading XYZ coordinates.'
         ENDIF
         STOP
      ENDIF

C-----------------------------------------------------------------------
C Print results
C-----------------------------------------------------------------------
      WRITE(*,*) ' '
      WRITE(*,*) 'READ_XYZIN test successful.'
      WRITE(*,*) ' '

      WRITE(*,*) 'Number of atoms: ', NAT
      WRITE(*,*) 'Atoms and coordinates:'

      DO I = 1, NAT
         WRITE(*,'(I3,2X,A2,3F12.6)') I, SYM(I),
     &        XYZ(1,I), XYZ(2,I), XYZ(3,I)
      END DO

      WRITE(*,*) ' '
      WRITE(*,*) 'BASIC section:'
      WRITE(*,*) '  CHARGE            = ', CHARGE
      WRITE(*,*) '  SPIN MULTIPLICITY = ', MULT
      WRITE(*,*) '  POINT GROUP       = ', PG

      WRITE(*,*) ' '
      WRITE(*,*) 'ROTATIONAL section:'
      WRITE(*,*) '  A (MHz)           = ', A
      WRITE(*,*) '  B (MHz)           = ', B
      WRITE(*,*) '  C (MHz)           = ', C
      WRITE(*,*) '  KAPPA             = ', KAPPA
      WRITE(*,*) '  SIGMA             = ', SIGMA
      WRITE(*,*) '  TEMPERATURE (K)   = ', T
      WRITE(*,*) '  Q_TRANS           = ', QTRANS
      WRITE(*,*) '  Q_ROT_ISO         = ', QROT_ISO
      WRITE(*,*) '  Q_ROT_AVG         = ', QROT_AVG

      WRITE(*,*) ' '
      WRITE(*,*) 'End of test.'

      END

