      SUBROUTINE READ_ROTATIONAL_SECTION(
     &     A_MHZ, B_MHZ, C_MHZ,
     &     KAPPA,
     &     ROTOR,
     &     FOUND
     & )
C=========================================================
C  Read #ROTATIONAL section from Merlino enriched XYZ file
C  (file name is fixed: "xyzin")
C
C  Output:
C     A_MHZ, B_MHZ, C_MHZ  - rotational constants (MHz)
C     KAPPA               - Ray's asymmetry parameter
C     ROTOR               - rotor type (string)
C     FOUND               - .TRUE. if section was found
C
C  Notes:
C     - Section order is irrelevant
C     - Unknown keywords are ignored
C     - Reading stops at next '#' or EOF
C=========================================================

      IMPLICIT NONE

      DOUBLE PRECISION A_MHZ, B_MHZ, C_MHZ, KAPPA
      CHARACTER*64 ROTOR
      LOGICAL FOUND

      INTEGER IOS
      CHARACTER*256 STR, KEY
      DOUBLE PRECISION VAL

C --- Defaults
      FOUND = .FALSE.
      A_MHZ = 0.0D0
      B_MHZ = 0.0D0
      C_MHZ = 0.0D0
      KAPPA = 0.0D0
      ROTOR = ' '

C --- Open fixed input file
      OPEN(UNIT=20, FILE='xyzin', STATUS='OLD', IOSTAT=IOS)
      IF (IOS .NE. 0) THEN
         WRITE(*,*) 'ERROR: cannot open xyzin'
         RETURN
      ENDIF

C --- Scan file for #ROTATIONAL
10    CONTINUE
      READ(20,'(A)',IOSTAT=IOS) STR
      IF (IOS .NE. 0) GO TO 99

      IF (STR(1:11) .EQ. '#ROTATIONAL') THEN
         FOUND = .TRUE.
         GO TO 20
      ENDIF

      GO TO 10

C --- Read section content
20    CONTINUE
      READ(20,'(A)',IOSTAT=IOS) STR
      IF (IOS .NE. 0) GO TO 99

C --- Stop if new section starts
      IF (STR(1:1) .EQ. '#') GO TO 99

C --- Parse known keywords
      READ(STR,*,IOSTAT=IOS) KEY

      IF (KEY .EQ. 'A_MHz') THEN
         READ(STR,*,IOSTAT=IOS) KEY, A_MHZ
      ELSE IF (KEY .EQ. 'B_MHz') THEN
         READ(STR,*,IOSTAT=IOS) KEY, B_MHZ
      ELSE IF (KEY .EQ. 'C_MHz') THEN
         READ(STR,*,IOSTAT=IOS) KEY, C_MHZ
      ELSE IF (KEY .EQ. 'KAPPA') THEN
         READ(STR,*,IOSTAT=IOS) KEY, KAPPA
      ELSE IF (KEY .EQ. 'ROTOR') THEN
         READ(STR,*,IOSTAT=IOS) KEY, ROTOR
      ENDIF

      GO TO 20

99    CONTINUE
      CLOSE(20)
      RETURN
      END

