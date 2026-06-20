      PROGRAM SMILES_TO_ZMATRIX
C     Trasforma una stringa SMILES in una Z-Matrix standard
C     (elenco di atomi, legami, angoli e diedri)
      IMPLICIT NONE

      INTEGER MAX_ATOMS, MAX_STACK, pos, N, top, i
      PARAMETER (MAX_ATOMS = 100, MAX_STACK = 100)
      CHARACTER*100 SMILES
      CHARACTER*2 atom_symbol_arr(MAX_ATOMS)
      INTEGER ref1(MAX_ATOMS), ref2(MAX_ATOMS), ref3(MAX_ATOMS)
      REAL bond_length(MAX_ATOMS), bond_angle(MAX_ATOMS), dihedral_angle(MAX_ATOMS)
      INTEGER stack(MAX_STACK)
      INTEGER len_smiles
      CHARACTER*1 ch, nextch
      CHARACTER*2 atom_symbol

C     Definisci la stringa SMILES (modifica qui secondo le tue esigenze)
C     Esempio: "CC(O)C" rappresenta una catena lineare con una ramificazione ossidrilica
      DATA SMILES /'CC(O)C'/
      len_smiles = LEN_TRIM(SMILES)

C     Inizializza gli array
      DO i = 1, MAX_ATOMS
         atom_symbol_arr(i) = '  '
         ref1(i) = 0
         ref2(i) = 0
         ref3(i) = 0
         bond_length(i) = 1.0
         bond_angle(i) = 109.5
         dihedral_angle(i) = 180.0
      END DO
      N = 0
      top = 0

      pos = 1
10    IF (pos .GT. len_smiles) GOTO 100

      ch = SMILES(pos:pos)
      IF (ch .EQ. '(') THEN
C        Apertura di una ramificazione: salva l'indice dell'atomo corrente nello stack
         top = top + 1
         stack(top) = N
      ELSE IF (ch .EQ. ')') THEN
C        Chiusura di una ramificazione: ripristina l'indice dell'atomo salvato
         N = stack(top)
         top = top - 1
      ELSE IF (ch .EQ. '-') THEN
C        Legame singolo esplicito: ignorato (usa valore di default per bond_length)
      ELSE IF (ch .EQ. '=') THEN
C        Doppio legame: ignorato per la Z-Matrix
      ELSE IF (ch .EQ. '#') THEN
C        Triplo legame: ignorato per la Z-Matrix
      ELSE IF (ch .GE. '0' .AND. ch .LE. '9') THEN
C        Chiusura di anello: non gestita in questo esempio
      ELSE IF (ch .GE. 'A' .AND. ch .LE. 'Z') THEN
C        Riconosce l'inizio di un simbolo atomico
         atom_symbol(1:1) = ch
         IF (pos .LT. len_smiles) THEN
            nextch = SMILES(pos+1:pos+1)
            IF (nextch .GE. 'a' .AND. nextch .LE. 'z') THEN
               atom_symbol(2:2) = nextch
C              Incrementa pos per consumare il carattere minuscolo
               pos = pos + 1
            ELSE
               atom_symbol(2:2) = ' '
            END IF
         ELSE
            atom_symbol(2:2) = ' '
         END IF
C        Aggiunge il nuovo atomo
         N = N + 1
         atom_symbol_arr(N) = atom_symbol
C        Determina i riferimenti per la Z-Matrix:
C         - Per il primo atomo non ci sono riferimenti.
C         - Il secondo atomo usa il primo.
C         - Il terzo usa il secondo e il primo.
C         - Dal quarto in poi si usano in maniera elementare gli atomi precedenti.
         IF (N .EQ. 1) THEN
            ref1(N) = 0
            ref2(N) = 0
            ref3(N) = 0
         ELSE IF (N .EQ. 2) THEN
            ref1(N) = 1
            ref2(N) = 0
            ref3(N) = 0
         ELSE IF (N .EQ. 3) THEN
            ref1(N) = 2
            ref2(N) = 1
            ref3(N) = 0
         ELSE
            ref1(N) = N - 1
            ref2(N) = ref1(N-1)
            IF (ref2(N) .EQ. 0) THEN
               ref2(N) = 1
            END IF
            ref3(N) = ref1(ref2(N))
            IF (ref3(N) .EQ. 0) THEN
               ref3(N) = 1
            END IF
         END IF
      END IF
      pos = pos + 1
      GOTO 10

100   PRINT *, 'Z-Matrix:'
C   Stampa della Z-Matrix in formato standard:
C   - Righe: numero, simbolo atomo, [ref1, distanza], [ref2, angolo], [ref3, diedro]
      DO i = 1, N
         IF (i .EQ. 1) THEN
            WRITE(*,*) i, atom_symbol_arr(i)
         ELSE IF (i .EQ. 2) THEN
            WRITE(*,*) i, atom_symbol_arr(i), ref1(i), bond_length(i)
         ELSE IF (i .EQ. 3) THEN
            WRITE(*,*) i, atom_symbol_arr(i), ref1(i), bond_length(i), ref2(i), bond_angle(i)
         ELSE
            WRITE(*,*) i, atom_symbol_arr(i), ref1(i), bond_length(i), ref2(i), bond_angle(i), ref3(i), dihedral_angle(i)
         END IF
      END DO

      END


