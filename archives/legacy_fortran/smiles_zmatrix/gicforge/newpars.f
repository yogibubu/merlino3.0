      Implicit None 
      Integer MxAt,MxBnd,MxRng,MxAtRn,MxBrn,In,IOut,IPunch,NAtoms,NRing
      Integer NBrnch,NSmi,I,II,J,JJ,IAt,JAt,IXYZ,ICyc,Num,NZ,NBTot
      Parameter(MxAt=100,MxBnd=4,MxRng=10,MxBrn=10,MxAtRn=10)
      Integer IAn(MxAt),IArom(MxAt),IRing1(MxAt),IBrnch(MxAt)
      Integer NPi(MxAt),NBond(MxAt),IBond(MxBnd,MxAt)
      Integer NH(MxAt),ISter1(MxAt),ISter2(3,MxAt)
      Integer Isotop(MxAt),IZ(4,MxAt),NAtRn1(MxRng)
      Integer LBl(MxAt),LAlpha(MxAt),LBeta(MxAt)
      Integer IRING(3,MxAt),NatRng(MxRng),IAtRng(MxAtRn,MxRng)
      Integer IAnZ(MxAt),IAtomB(2,MxAt)
      Real*8 Charge(MxAt),BndOrd(MxBnd,MxAt)
      Real*8 BL(MxAt),Alpha(MxAt),Beta(MxAt)
      Real*8 Pi,ToAng,ToDeg,ToRad,ScalI
      Real*8 C(3,MxAt),CZ(3,MxAt)
      Real*8 A(MxAt),B(MxAt),D(MxAt),Alpha1(MxAt),Beta1(MxAt)
      Character*100 TITLE
      Character*1 SMILES(100)
      Character*2 El,IAnEl2
      Character*80 LinScr
      Logical ActivR(MxRng)
      Logical FndSp,Error,TstAng,Ttest
C Define I/O files
      In=5
      Iout=6
      IPunch=7
      OPEN(In,FILE='provin',STATUS='UNKNOWN')
      OPEN(IOut,FILE='provout',STATUS='UNKNOWN')
      OPEN(IPunch,FILE='provXYZ',STATUS='UNKNOWN')
      Do i=1,100
       SMILES(i)=' '
       TITLE(i:i)=' '
      End Do
      Read(In,*) Title
      Write(IOut,'('' '')')
      Write(IOut,*) Title
      NSMI=1
      FndSp=.False.
   10 Read(In,'(A1)',advance='no',end=20) SMILES(NSMI)
      If(SMILES(NSMI).ne.' ') then
       NSMI=NSMI+1
       FndSP=.true.
       GOTO 10
      ELSEIF(.not.FndSp) then
       GOTO 10
      ENDIF
   20 Write(IOut,*)(SMILES(I),I=1,NSMI)
      Write(IOut,'('' '')')
      Do IAt=1,MxAt
       IAn(IAt)=0
       IRing1(IAt)=0
       IArom(IAt)=0
       IBrnch(IAt)=0
       NPi(IAt)=0
       NBond(IAt)=0
       NH(IAt)=0
       ISter1(IAt)=0
       Isotop(IAt)=0
       IRing(1,IAt)=0
       IRing(2,IAt)=0
       IRing(3,IAt)=0
       If(IAt.le.MxRng) NAtRn1(IAt)=0    
       Do JJ=1,MxBnd
        If(JJ.le.3) ISter2(JJ,IAt)=0
        IBond(JJ,IAt)=0
        BndOrd(JJ,IAt)=0.0D0
       End Do
      End Do
C Analyze SMILES
      CALL SMIF77(IOut,NSMI,SMILES,MxBnd,MxAtRn,MxRng,NAtoms,NRing,
     $  NBrnch,IAn,IArom,NPi,IRing1,IBrnch,NBond,IBond,Charge,BndOrd,
     $  NH,ISter1,ISter2,Isotop,NatRn1,IRing,NatRng,IAtRng,NBTot,
     $  ActivR,IAtomB)
C Set IZ
      CALL IClear(4*NAtoms,IZ)
      CALL SetIZ(IOut,MxBnd,NAtoms,IAn,IBrnch,NH,NBond,IBond,IZ)
      Do IAt=1,NAtoms
       El=IAnEl2(IAn(IAt))
       If(IAt.eq.1) then
        BL(IAt)=0.0D0
        Alpha(IAt)=0.0D0
        Beta(IAt)=0.0D0
        IZ(1,IAt)=0
        IZ(2,IAt)=0
        IZ(3,IAt)=0
        write(IOut,'(A2,I3)') IAnEl2(IAn(IAt)),IAt
       ElseIf(IAt.eq.2) then
        Alpha(IAt)=0.0
        Beta(IAt)=0.0D0
        IZ(2,IAt)=0
        IZ(3,IAt)=0
        write(IOut,'(A2,I3,5X,I3)') IAnEl2(IAn(IAt)),IAt,IZ(1,IAt)
       ElseIf(IAt.eq.3) then
        Beta(IAt)=0.0D0
        IZ(3,IAt)=0
        write(IOut,'(A2,I3,5X,2I3)') IAnEl2(IAn(IAt)),IAt,IZ(1,IAt),
     $    IZ(2,IAt)
       Else
        write(IOut,'(A2,I3,5X,4I3)') El,IAt,(IZ(ii,IAt),ii=1,4) 
        BL(IAt)=1.0d0
        Alpha(IAt)=109.4722D0
        Beta(IAt)=180.0D0
       EndIf
      ENDDO 
C Set Bond Lengths and Valence Angles
      CALL SetBLA(MxBnd,NAtoms,IAn,NPi,IZ,NBond,IBond,BndOrd,BL,Alpha)
C Set Dihedral Angles
      CALL SetDih(IOut,MxBnd,MXatRn,NAtoms,IZ,NBond,IBond,IBrnch,IRing,
     $  NAtRng,IAtRng,Alpha,Beta)
      CALL IClear(NAtoms,LBl)
      CALL IClear(NAtoms,LAlpha)
      CALL IClear(NAtoms,LBeta)
C Temporary
      ToAng=1.0d0
      ToDeg=1.0d0
C Print Z-matrix
      CALL ZMTPRT(IOut,NAtoms,IAN,IZ,LBL,LAlpha,LBeta,BL,ALPHA,BETA,
     $  ToAng,ToDeg)
C transform to Cartesian coordinates
C checking tetrahetral angles (TTest) and 0>teta<180 (TstAng)
C reset fake IZ
      IZ(1,1)=0
      IZ(2,1)=0
      IZ(2,2)=0
      TstAng=.true.
      TTest=.true.
      NZ=NAtoms
      CALL IMove(NAtoms,IAn,IAnZ)
      CALL AClear(3*NAtoms,C)
      CALL AClear(3*NAtoms,CZ)
      CALL AClear(NAtoms,A)
      CALL AClear(NAtoms,B)
      CALL AClear(NAtoms,D)
      CALL Aclear(NAtoms,Alpha1)
      CALL Aclear(NAtoms,Beta1)
C Angles should be in radiants
      pi=4.0d0*ATan(1.0d0)
      ToRad=pi/180.0d0
      Error=.False.
      Do IAt=1,NAtoms
       Alpha(IAt)=Alpha(IAt)*ToRad
       Beta(IAt)=Beta(IAt)*ToRad
      EndDo
      CALL ZtoC(MxAt,NZ,IAnZ,IZ,Bl,Alpha,Beta,TTest,NAtoms,IAn,C,
     $   CZ,A,B,D,Alpha1,Beta1,IOut,Error,TstAng)
      If(Error) Stop
      Write(IOut,'(/,10X,''Cartesian Coords. from SMILES'')')
      Do IAt=1,NAtoms
       Write(IOut,'(I5,2X,A2,3F12.5)') IAt,IAnEl2(IAn(IAt)),
     $  (C(IXYZ,IAt),IXYZ=1,3)
      ENDDO
C Print interatomic distances
      CALL LlinCl(LinScr)
      LinScr(1:22)=' Interatomic Distances (Angstrom)'
      CALL HedPrt(IOut,0,LinScr,Num)
      ScalI=1.0d0
      CALL DisMat(NAtoms,IAN,C,2,5,IOut,Error,0,ScalI)
      write(IOut,'('' '')')
C Write XYZ file
      Write(IPunch,'(I5)') NAtoms
      Write(IPunch,'(A)') Title
      Do IAt=1,NAtoms
      Write(IPunch,'(A2,3F12.5)')IANEl2(IAn(IAt)),(C(IXYZ,IAt),IXYZ=1,3)
      ENDDO
      End
*Deck SMIF77
      SUBROUTINE SMIF77(IOut,NSMI,SMILES,MxBnd,MxAtRn,MxRng,NAtoms,
     $  NRing,NBrnch,IAn,IArom,NPi,IRing1,IBrnch,NBond,IBond,Charge,
     $  BndOrd,NH,ISter1,ISter2,Isotop,NAtRn1,IRing,NAtRng,IAtRng,
     $  NBTot,ActivR,IAtomB)
      Implicit Integer (A-N)
      Dimension IAn(*),NPi(*),IArom(*),IRing1(*),IBrnch(*),NBond(*)
      Dimension NH(*),IBond(MxBnd,*),ISter1(*),ISter2(3,*)
      Dimension Isotop(*),NAtRn1(*),IRing(3,*),NAtRng(*)
      Dimension IAtRng(MxAtRn,*),IAtomB(2,*)
      Character*1 SMILES(NSMI)
      Real*8 Charge(*),BndOrd(MxBnd,*)
      Logical ActivR(MxRng)
C Parse the SMILES and determine Atom properties and bonds
      CALL ParSMI(IOut,NSMI,MxBnd,MxRng,MxAtRn,SMILES,NAtoms,NBrack,
     $  NRing,NBrnch,IAn,IArom,IRing1,IBrnch,NBond,IBond,Charge,BndOrd,
     $  NH,ISter1,NSter2,ISter2,Isotop,NAtRng,IAtRng,ActivR)
C Manage different storage of atoms belonging to rings
      CALL Bridge(Iout,MxAtRn,NAtoms,NRing,IRing,NAtRng,IAtRng)
C Assign Pi bonds and bond orders
      CALL PiBond(IOut,MxBnd,NAtoms,IArom,NPi,NBond,IBond,BndOrd)
C Add Hydrogen Atoms both explicit and implicit
      NHeavy=NAtoms
      CALL AddHyd(IOut,MxBnd,NHeavy,NAtoms,IAn,NPi,NBond,IBond,
     $  BndOrd,NH)
C Print Results
      CALL SmiRes(IOut,MxBnd,MxAtRn,NAtoms,IAn,NPi,IArom,IBrnch,IRing,
     $  NBond,IBond,BndOrd,NH,ISter1,NSter2,ISter2,Isotop)
      RETURN
      END
*Deck ParSMI
      SUBROUTINE ParSMI(IOut,NSMI,MxBnd,MxRng,MxAtRn,SMILES,NAtoms,
     $  NBrack,NRing,NBrnch,IAn,IArom,IRing,IBrnch,NBond,IBond,Charge,
     $  BndOrd,NH,Stereo,NSter2,ISter2,Isotop,NAtRng,IAtRng,ActivR)
      Implicit Integer (A-Z)
      Dimension IAn(*),IArom(*),IRing(3,*),IBrnch(*),NBond(*),NH(*)
      Dimension IBond(MxBnd,*),Stereo(*),ISter2(3,*)
      Dimension Isotop(*)
      Dimension NAtRng(*),IAtRng(MxAtRn,*)
      Real*8 Charge(*),BndOrd(MxBnd,*)
      Logical ActivR(MxRng),Open
      Character*1 SMILES(NSMI)
      Character*2 El
      Character*1 Test,Test1
C Variabili locali
      Logical ActivB(10),ActBrc,UpdRng
      Dimension IFirst(MxRng),ILast(MxRng),IPrev(MxRng),INext(MxRng)
      Dimension ISaved(MxRng),IBefBr(10),IBrFrst(MxRng)
      
      Do IAct=1,10
       ActivR(IAct) = .False.
       ActivB(IAct) = .False.
       IFirst(IAct) = 0
       ILast(IAct) = 0 
       IBefBr(IAct) = 0
      End Do
      NAtoms = 0
      NBrack = 0
      NRing = 0
      NBrnch = 0
      BndTyp = 1.0d0
      IRngMx = 0
      MaxBrn = 0
      NSter2  = 0
      UpdRng = .True.
 
      I = 1
  10  CONTINUE
C       write(IOut,'('' From ParSMI I ='',I3)') I
        IF (I.GT.NSMI) GOTO 100
        Test = SMILES(I)
        
C Gestione dei legami
        IF (Test.eq.'-') then
          BndTyp =1.0D0
          I=I+1
          GOTO 10
        ELSEIF (Test .EQ. '=') THEN
          BndTyp = 2.0d0
          IDBnd  = I
          I=I+1
          GOTO 10
        ELSEIF(Test.eq.'#') then
          BndTyp = 3.0d0
          I = I + 1
          GOTO 10
        ENDIF

C Gestione degli stereocentri
        IF (Test .EQ. '@') THEN
          I = I + 1
          Test1 = SMILES(I)
          IF (Test1 .EQ. '@') THEN
C Stereoxchimica Inversa (@@)
            Stereo(NAtoms) = -1  
            I = I + 1
          ELSE
C Stereochimica Assoluta (@)
            Stereo(NAtoms) = 1 
          ENDIF
          GOTO 10
        ENDIF

C Gestione Isomeria E/Z
        If(Test.eq.'/'.or.Test.eq.'\') THEN
         CALL SetEZ(Test,IAt,NSter2,ISter2)
         I = I + 1
         GOTO 10
        ENDIF

C Gestione della carica ( + o - )
        IF (Test .EQ. '+' .OR. Test .EQ. '-') THEN
         Sign = 1
         IF (Test .EQ. '-') Sign = -1
         I = I + 1
         Test1 = SMILES(I)
         IF (Test1 .GE. '0' .AND. Test1 .LE. '9') THEN
C Correggi Assegnazione Carica
          Charge(Natoms) = Sign * (ICHAR(Test1) - ICHAR('0'))
          I = I + 1
         ELSE
          Charge(NAtoms) = Sign
         ENDIF
        ENDIF

C Gestione delle parentesi quadre [ ]
        IF (Test.EQ.'[') THEN
         If(NBrck.gt.MxBrc) MxBrC=NBrck
         CALL GetBrc(IOut,NAtoms,JAt,NBrck,'open',ActBrc,IBfBrc)
          I = I + 1
C Controlla se c'è un numero di massa (1-2 cifre)
          Test1=SMILES(I)
          MssTmp = 0
          IF (Test1 .GE. '0' .AND. Test1 .LE. '9') THEN
           MssTmp = ICHAR(Test1) - ICHAR('0')
           I = I + 1
           Test1 = SMILES(I)
           IF (Test1 .GE. '0' .AND. Test1 .LE. '9') THEN
            MssTmp = MssTmp * 10 + (ICHAR(Test1) - ICHAR('0'))
            Isotop(NAtoms+1)=MssTmp
            I = I + 1
           ENDIF
C Aggiorna il prossimo carattere
           Test = SMILES(I)  
          ENDIF
          GOTO 10
        Else If(Test.EQ.']') then
         Call GetBrc(IOut,NAtoms,JAt,NBrck,'clos',ActBrc,IBfBrc)
C Gestione della sequenza '(['
         If(SMILES(I-1).eq.')') then
          write(IOut,'('' The sequence )[ is not allowed'')')
          Stop
         EndIf
         I = I + 1
         GOTO 10
        EndIf

C Gestione dei rami ()
        IF (Test .EQ. '(') THEN
          If(NBrnch.gt.MaxBrn) MaxBrn=NBrnch
          Call GetBrn(IOut,NAtoms,JAt,NBrnch,'open',ActivB,IBefBr)
C Gestione della sequenza ']('
          If(SMILES(I-1).eq.']') JAt=IBfBrc+1 
          I = I + 1
          GOTO 10
        ELSE IF (Test .EQ. ')') THEN
          Call GetBrn(IOut,NAtoms,JAt,NBrnch,'clos',ActivB,IBefBr)
          I=I+1
          GOTO 10
        ENDIF

C Gestione degli atomi (lettere maiuscole)
        IF (Test .GE. 'A' .AND. Test .LE. 'Z') THEN
         El(1:1) = Test
         El(2:2) = ' '
         N2Ring=0
         J = I + 1
         IF (J .LE. L) THEN
          Test1 = SMILES(J)
C Gestione degli anelli
          IF (Test1 .GE. '0' .AND. Test1 .LE. '9') THEN
           N0Ring = NRing
           NRing  = ICHAR(Test1) - ICHAR('0')
           If(NRing.gt.IRngMx) IRngMx=NRing
C          Write(IOut,'('' entering GetRng with NRing ='',I2)')NRing 
           CALL GetRng(IOut,MxAtRn,IAt,N0Ring,NRing,NBrnch,ISaved,
     $       ActivR,IFirst,ILast,IPrev,INext,IBrFrst,NAtRng,IAtRng,
     $       Open)
           UpdRng=.false.
           N1Ring=NRing
           N2Ring=0
           If(Open) then
            N2Ring=N0ring
            Open=.false.
           EndIf
           I=I+1
C caso di apertura o chiusura di 2 anelli sullo stesso atomo
           J=I+1
           IF(J.LE.L) THEN
            Test1 = SMILES(J)
            IF (Test1 .GE. '0' .AND. Test1 .LE. '9') THEN
             N0Ring = N1Ring
             NRing = ICHAR(Test1) - ICHAR('0')
             If(NRing.gt.IRngMx) IRngMx=NRing
C            Write(IOut,'('' entering GetRng with NRing ='',I2)')NRing
             CALL GetRng(IOut,MxAtRn,IAt,N0Ring,NRing,NBrnch,ISaved,
     $         ActivR,IFirst,ILast,IPrev,INext,IBrFrst,NAtRng,IAtRng,
     $         Open)
             UpdRng=.false.
             N2Ring=NRing
             N3Ring=0
             If(Open) then
C             N3Ring=N0Ring
              Open=.false.
             EndIf 
             I=I+1
            EndIf
           EndIf
          ElseIF (Test1 .GE. 'a' .AND. Test1 .LE. 'z') THEN
           El(2:2) = Test1
C Avoid problems with normal atoms followed by aromatic atoms
           If(El.eq.'Cc'.or.El.eq.'Cn'.or.El.eq.'Cp') El(2:2)=' '
           If(El.eq.'Nc'.or.El.eq.'Nn'.or.El.eq.'No') El(2:2)=' '
           If(El.eq.'Np'.or.El.eq.'Ns'.or.El.eq.'Pc') El(2:2)=' '
           If(El.eq.'Oc'.or.El.eq.'On'.or.El.eq.'Oo') El(2:2)=' '
           If(El.eq.'Op') El(2:2)=' '
           If(El.eq.'Pn'.or.El.eq.'Po'.or.El.eq.'Pp') El(2:2)=' '
           If(El.eq.'Ps'.or.El.eq.'So'.or.El.eq.'Sp') El(2:2)=' '
           If(El.eq.'Ss') El(2:2)=' '
           If(El(2:2).ne.' ') I = I + 1
          ElseIf(Test1.gt.'1'.and.Test1.le.'9'.and.Test.eq.'H') THEN
           NH(IAt)=NH(IAt) + IChar(Test1) - IChar('0')
           I=I+2
           GOTO 10
          ElseIf(Test.eq.'H') then 
           NH(IAt)=NH(IAt)+1
           I=I+1
           GOTO 10
          EndIf
         EndIf
         IAt=NAtoms+1

         If(UpdRng) then
C         Write(IOut,'('' Entering SetRAt'')')
          CALL SetRAt(MxAtRn,NRing,IAt,NatRng,IAtRng)
         Else
          UpdRng=.true.
         EndIf
         IArom(IAt) = 0
         IRing(1,IAt) = N1Ring
         If(N2Ring.ne.0) then
          IRing(2,IAt)=N2Ring
C         If(N3Ring.ne.0.and.N3Ring.ne.N2Ring) then
C          IRing(3,IAt)=N3Ring
C         Else
           IRing(3,IAt) = 0
C         EndIf
         EndIf
         IBrnch(IAt)= NBrnch
         IAn(IAt)=El2IAn(El) 
         IF(IAt.ne.1.and.IAn(IAt).gt.1) then
          If(JAt.eq.0) JAt=IAt-1
          NBond(IAt) = NBond(IAt) + 1
          NBond(JAt) = NBond(JAt) + 1
          IBond(NBond(IAt), IAt) = JAt
          IBond(NBond(JAt), JAt) = IAt
          BndOrd(NBond(IAt), IAt) = BndTyp
          BndOrd(NBond(JAt), JAt) = BndTyp
         End If 
         JAt=0
         BndTyp=1.0d0
         NAtoms=NAtoms+1
         I = I + 1
         GOTO 10
        END IF

C Gestione degli atomi aromatici (lettere minuscole o maiuscole per C e N)
        IF (Test .GE. 'a' .AND. Test .LE. 'z') THEN
         IF(Test .ne. 'c' .and. Test .ne. 'n' .and. Test .ne. 'o' .and. 
     $    Test .ne. 'p'. and. Test .ne. 's') THEN
          I=I+1
          GOTO 10
         ELSE
          El(1:1) = CHAR(ICHAR(Test) - 32)  ! Converti in maiuscolo (es. 'c' → 'C')
          El(2:2) = ' '
          N2Ring=0
C Gestione degli anelli
          J = I + 1
          IF (J .LE. L) THEN
           Test1 = SMILES(J)
           IF (Test1 .GE. '0' .AND. Test1 .LE. '9') THEN
            N0Ring=NRing
            NRing = ICHAR(Test1) - ICHAR('0')
            If(NRing.gt.IRngMx) IRngMx=NRing
C           Write(IOut,'('' entering GetRng with NRing ='',I2)')NRing
            CALL GetRng(IOut,MxAtRn,IAt,N0Ring,NRing,NBrnch,ISaved,
     $        ActivR,IFirst,ILast,IPrev,INext,IBrFrst,NAtRng,IAtRng,
     $        OPen)
            UpdRng=.False.
            N1Ring=NRing
            N2Ring=0
            If(Open) then
             N2Ring=N0Ring 
             Open=.false.
            EndIf 
            I=I+1
           EndIf
          EndIf
C caso di apertura o chiusura di 2 anelli sullo stesso atomo
          J=I+1
          IF(J.LE.L) THEN
           Test1 = SMILES(J)
           IF (Test1 .GE. '0' .AND. Test1 .LE. '9') THEN
            N0Ring = N1Ring
            NRing = ICHAR(Test1) - ICHAR('0')
            If(NRing.gt.IRngMx) IRngMx=NRing
C           Write(IOut,'('' entering GetRng with NRing ='',I2)')NRing
            CALL GetRng(IOut,MxAtRn,IAt,N0Ring,NRing,NBrnch,ISaved,
     $        ActivR,IFirst,ILast,IPrev,INext,IBrFrst,NAtRng,IAtRng,
     $        Open)
            UpdRng=.false.
            N2Ring=NRing
            N3Ring=0
            If(OPen) then
C            N3Ring=N0Ring
             Open=.False.
            EndIf
            I=I+1
           EndIf
          EndIf
          IAt=NAtoms+1
          If(UpdRng) then
           CALL SetRAt(MxAtRn,NRing,IAt,NatRng,IAtRng)
          Else
           UpdRng=.true.
          EndIf
          IArom(IAt) = 1
          IRing(1,IAt) = N1Ring
          If(N2Ring.ne.0) then
           IRing(2,IAt) = N2Ring
C          If(N3Ring.ne.0.and.N3Ring.ne.N2Ring) then
C           IRing(3,IAt) = N3Ring
C          Else
            IRing(3,IAt) = 0
C          EndIf  
          EndIf
          IBrnch(IAt)= NBrnch
          IAn(IAt)=El2IAn(El) 
          If(IAt.gt.1) then
           If(JAt.eq.0) JAt=IAt-1
           If(IArom(JAt).eq.1) BndTyp=1.5D0
           NBond(IAt) = NBond(IAt) + 1
           NBond(JAt) = NBond(JAt) + 1
           IBond(NBond(IAt),IAt) = JAt
           IBond(NBond(JAt),JAt) = IAt
           BndOrd(NBond(IAt),IAt) = BndTyp
           BndOrd(NBond(JAt),JAt) = BndTyp
          EndIf 
          JAt=0
          BndTyp=1.0d0
          NAtoms=NAtoms+1

          I=I+1
          GOTO 10
         ENDIF
        ENDIF  

 100  CONTINUE
      NRing=IRngMx  
      If(NRing.eq.0) Return 
C Creazione legami fra 1 e ultimo atomo anelli
      write(IOut,'(I2,'' Rings'')') NRing
      Do Irr=1,NRing
       IAt=IFirst(IRr)
       JAt=ILast(IRr)
       IRing(1,IAt)=IRr
       IRing(1,JAt)=IRr 
      Write(IOut,'('' Ring'',I2,'' is closed by atoms:'',2I3)') 
     $   IRr,IAt,Jat
       If(IAt.gt.0.and.JAt.gt.0) then
        NBond(IAt) = NBond(IAt) +1
        NBond(JAt) = NBond(JAt) + 1
        IBond(NBond(IAt),IAt) = JAt
        IBond(NBond(JAt),JAt) = IAt
        BndTyp=1.0d0
        If(IArom(IAt).eq.1.and.IArom(JAt).eq.1) BndTyp=1.5D0
        BndOrd(NBond(IAt),IAt) = BndTyp
        BndOrd(NBond(JAt),JAt) = BndTyp
       Else
        Write(IOut,'(A,I3)') ' Error: Ring Closure Failed for ring',IRr
       ENDIF
      ENDDO
C Elimina Atomi Spuri dagli anelli
      Do IAt=1,NAtoms
       If(NBond(IAt).eq.1) IRing(1,IAt)=0
       write(IOut,'('' Atom'',I3,'' Ring'',3I3)') IAt,(IRing(I,IAt),
     $   I=1,3)
      EndDo
      Call IClear(NRing,NAtRng)
      Do IAt=1,NAtoms
       Do II=1,3
        If(IRing(II,IAt).ne.0) then
         IRI=IRing(II,IAt)
         NAtRng(IRI)=NatRng(IRI)+1
         NRI=NAtRng(IRI)
         IAtRNg(NRI,IRI)=IAt
        EndIf
       EndDo
      EndDo
      Do IR=1,NRing
       Write(IOut,'('' Ring'',I2,'' has the following'',I3,'' Atoms:''
     $   )',advance='no') IR,NAtRng(IR)
       Write(IOut,'(10I3)') (IAtRng(IAt,IR),IAt=1,NAtRng(IR))
      EndDo
      RETURN
      END
*Deck Bridge
      Subroutine Bridge(Iout,MxAtRn,NAtoms,NRing,IRing,NAtRng,IAtRng)
      Implicit None
C Builds array IRING(3,NAtoms) from arrays NAtRng(NRing) and IAtRng(NN,NAtoms)
      Integer IOut,MxAtRn,NAtoms,NRing
      Integer IRing(3,*),NAtRng(*),IAtRng(MxAtRn,*)
C Local
      Integer IR,II,IAt
C No Rings
      If(NRing.eq.0) Return
      Do IR=1,NRing
C      Write(IOut,'('' Ring'',I2,'' has the following'',I3,'' Atoms:''
C    $  )',advance='no') IR,NatRng(IR) 
C      Write(IOut,'(10I3)') (IAtRng(IAt,IR),IAt=1,NAtRng(IR)) 
       Do II=1,NAtRng(IR)
        IAt=IAtRng(II,IR)
        If(IRing(1,IAt).eq.0) then
         IRing(1,IAt)=IR
        ElseIF(IRing(2,IAt).eq.0) then
         IRing(2,IAt)=IR
        ElseIf(IRing(3,IAt).eq.0) then
         IRing(3,IAt)=IR
        Else
         write(IOut,'('' For Atom'',I3,'' IRing3 ='',I3,'' IR ='',I3)') 
     $     IAt,IRing(3,IAt),IR
C        Stop
        EndIf
       EndDo
      EndDo
      Do IAt=1,NAtoms
       If(IRing(1,IAt).eq.0) then
        Cycle
       ElseIf(IRing(2,IAt).eq.0) then
        Write(IOut,'('' Atom'',I3,'' Belongs to Ring :'',I3)') IAt,
     $    IRing(1,IAt)
       ElseIf(IRing(3,IAt).eq.0) then
        Write(IOut,'('' Atom'',I3,'' Belongs to Rings:'',2I3)') IAt,
     $    IRing(1,IAt),IRing(2,IAt)      
       Else 
        Write(IOut,'('' Atom'',I3,'' Belongs to Rings:'',3I3)') IAt,
     $    IRing(1,IAt),IRing(2,IAt),IRing(3,IAt)
       EndIf
      EndDo 
      Return
      End
*Deck SetEZ
      Subroutine SetEZ(EZ,IAt,NSter2,ISter2)
      Implicit None
      Integer IEZ,IAt,NSter2,ISter2(3,*)
      Character*1 EZ
      If(EZ.eq.'/') then
       IEZ=-1
      ElseIF(EZ.eq.'\') then
       IEZ=1
      EndIf
      If(ISter2(1,NSter2).eq.0) then
       NSter2=NSter2+1
       ISter2(1,NSter2)=IAt
       ISter2(3,NSter2)=IEZ
      ElseIf(ISter2(1,NSter2).ne.0) then
       ISter2(2,NSter2)=IAt+1
       If(ISter2(3,NSter2).eq.IEZ) then
        ISter2(3,NSter2)=-1
       Else
        ISter2(3,NSter2)=1
       EndIf 
      EndIf
      Return
      End
*Deck GetBrn
      Subroutine GetBrn(IOut,IAt,JAt,IBrnch,Action,ActivB,IBefBr)
      Implicit None
      Integer IOut,I,IAt,JAt,IBrnch,IBefBr(10)
      Character Action*4
      Logical ActivB(10)
      If(IBrnch.lt.0.or.IBrnch.gt.9) then
       write(IOut,'(I3,'' Branches, but only 1 to 9 are allowed'')')
     $   IBrnch
       Stop
      ElseIf(IAt.le.0) then
       write(IOut,'(I3,'' Atoms, but only 1 to 99 are allowed'')') IAt
      EndIf
      do i=1,10
       If(ActivB(i)) IBrnch=I
      End Do
      If(action.eq.'stay') Return
      If(action.eq.'open') then
        IBrnch=IBrnch+1
        ActivB(IBrnch)=.true.
        IBefBr(IBrnch)=IAt 
      ElseIf(action.eq.'clos') then
       ActivB(IBrnch)=.false.
       JAt=IBefBr(IBrnch)
       IBrnch=IBrnch-1
      End If
      Return
      End
*Deck GetRng
      Subroutine GetRng(IOut,MxAtRn,IAt,IRing0,IRing,IBrnch,ISaved,
     $  ActivR,IFirst,ILast,IPrev,INext,IBrFrst,NAtRng,IAtRng,Open)
      Implicit None
      Integer IOut,MxAtRn,IAt,IRing0,IRing,IS,IBrnch
      Integer IFirst(*),ILast(*),NAtRng(*),IAtRng(MxAtRn,*)
      Integer IPrev(*),INext(*),ISaved(*),IBrFrst(*) 
      Logical ActivR(*),Open
      If(IRing.le.0.or.IRing.gt.9) then
       write(IOut,'(I3,'' Rings, but only 1 to 9 are allowed'')') IRing 
       Stop
      EndIf
      If(ActivR(IRing)) then
C Close ring
       Open=.false.
       ActivR(IRing)=.False.
       ILast(IRing)=IAt+1
       IF(IBrnch.eq.IBrFrst(IRing)) then
        INext(IRing)=ISaved(IRing)
       Else
C Check
        INext(IRing)=IRing+1
        ISaved(IRing)=IRing+1
       EndIf 
       NAtRng(IRing)=NatRng(IRing)+1
       IAtRng(NAtRng(IRing),IRing)=IAt+1
       If(ISaved(IRing).ne.0) then
        IS=ISaved(IRing)
        NAtRng(IS)=NAtRng(IS)+1
        IAtRng(NAtRng(IS),IS)=IAt+1
       EndIf
       Write(IOut,'('' Close Ring'',I3,'' First Atom'',I3,
     $   '' Last Atom'',I3,'' Son Ring'',I2,'' New Active Ring'',I2,
     $   '' Branch'',2I2)') IRing,Ifirst(IRing),ILast(IRing),
     $   IPrev(IRing),ISaved(IRing),IBrnch,IBrFrst(IRing)
       IRing=ISaved(IRing)
      Else
C Open Ring
       Open=.true.
       ActivR(IRing)=.True.
       IFirst(IRing)=IAt+1
       IBrFrst(IRing)=IBrnch 
       If(IRing.eq.1) then
        IFirst(IRing)=IAt
        IPrev(IRing)=0
        ISaved(IRing)=0
       Else
        IPrev(IRing)=IRing0
        ISaved(IRing)=Iring0
       EndIf
       NAtRng(IRing)=1
       IAtRng(NatRng(IRing),IRing)=IAt
       Write(IOut,'('' Open Ring '',I3,'' First Atom'',I3,
     $   '' Father Ring'',I3,'' ISaved='',I3,'' Branch'',I3)') IRing,
     $   Ifirst(IRing),IPrev(IRing),ISaved(IRing),IBrFrst(Iring)
      EndIf
      Return
      End
*Deck SetRAt
      Subroutine SetRAt(MxAtRn,IRing,IAt,NatRng,IAtRng)
      Integer NAtRng(*),IAtRng(MxAtRn,*)
      NAtRng(IRing)=NAtRng(IRing)+1
      IAtRng(NAtRng(IRing),IRing)=IAt 
      Return
      End
*Deck GetBrc
      Subroutine GetBrc(IOut,IAt,JAt,NBrck,Action,ActBrc,IBfBrc)
      Implicit None
      Integer IOut,I,IAt,JAt,NBrck,IBfBrc
      Character Action*4
      Logical ActBrc
      If(IAt.le.0) then
       write(IOut,'(I3,'' Negative Atom Number '')') IAt
      EndIf
      If(action.eq.'open') then
        ActBrc=.true.
        IBfBrc=IAt 
        JAt=IBfBrc
      ElseIf(action.eq.'clos') then
        ActBrc=.false.
        JAt=IBfBrc
      End If
      Return
      End
*Deck GetStr
      Subroutine GetStr(IOut,MxBnd,IAt,JAt,NBond,ISter2,Symbol)
      Implicit None
      Integer IOut, MxBnd, IAt, JAt
      Integer ISter2(MxBnd,*), NBond(*)
      Character Symbol*1
      If (IAt .LE. 0 .OR. JAt .LE. 0) Return
      If (Symbol .EQ. '//') Then
       ISter2(NBond(JAt),JAt) = -1
       ISter2(NBond(IAt),IAt) = -1
      ElseIf (Symbol .EQ. '\') Then
       ISter2(NBond(JAt),JAt) = 1
       ISter2(NBond(IAt),IAt) = 1
      EndIf
      Return
      End
*Deck PiBond
      Subroutine PiBond(IOut,MxBnd,NAtoms,IArom,NPi,NBond,IBond,BndOrd)     
      Integer IOut,MxBnd,NAtoms,IAt,JAt,NBI
      Integer II
      Integer IArom(*),NPi(*),NBond(*),IBond(MxBnd,*)
      Real*8 BndOrd(MxBnd,*)
      Do IAt=1,NAtoms
       NPI(IAt)=0
       NBI=NBond(IAt)
       Do II=1,NBI
        JAt=IBond(II,IAt)
        If(IArom(IAt).eq.1.and.IArom(JAt).eq.1) then
         NPI(IAt)=1
         NPI(JAt)=1
         BndOrd(II,IAt)=1.5D0
        ElseIf(BndOrd(II,IAt).eq.2.0d0) then
         NPI(IAt)=NPI(IAt)+1
        ElseIf(BndOrd(II,IAt).eq.3.0d0) then
         NPI(IAt)=3
        EndIf
       End Do
      End Do 
      Return
      End
*Deck SmiRes
      Subroutine SmiRes(IOut,MxBnd,MxAtRn,NAtoms,IAn,NPi,IArom,IBrnch,
     $  IRing,NBond,IBond,BndOrd,NH,ISter1,NSter2,ISter2,Isotop)
      Implicit Real*8 (A-H,O-Z)
      Dimension NPi(*),IArom(*),NH(*),NBond(*),IBond(MxBnd,*)
      Dimension IBrnch(*),IRing(3,*),ISter1(*),ISter2(3,*)
      Integer IAn(*),Isotop(*)
      Dimension BndOrd(MxBnd,*)
      Character*1 Star,Chir,OpenP,CloseP
      Character*1 CT(MxBnd) 
      Character*2 El,IAnEl2
C Stampa intestazione
      OpenP='('
      CloseP=')' 
      write(IOut,'(/,I5,'' Atoms with * labeling Aromatic Atoms'')') 
     $  NAtoms
      write(IOut,'(/,'' Atom Branch Ring  Pi El. nH Stereo '',
     $  6X,''Bonded Atoms'','' and Bond Orders'')')
C Ciclo per ogni atomo
      Do IAt = 1, NAtoms
       El=IAnEl2(IAn(IAt))
       If(IArom(IAt).ne.0) then
        star = '*'
       Else
        star = ' '
       EndIf
       If(ISter1(IAt).eq.1) then
        Chir='R'
       ElseIf(ISter1(IAt).eq.-1) then
        Chir='S'
       Else
        Chir=' '
       EndIf
       NBI = NBond(IAt)
       If(NBI.le.0.or.NBI.gt.MxBnd) then
        write(IOut,'('' Wrong Number of Bonds('',I2,'') for atom'',
     $    I3)')NBI,IAt
        Stop
C Scrittura dei dati in base al numero di legami (NBI)
       Else
        Do ii=1,NBI
         If(ISter2(ii,IAt).eq.0) then
          CT(ii)=' '
         ElseIf(ISter2(ii,IAt).eq.-1) then
          CT(ii)='E'
         ElseIf(ISter2(ii,IAt).eq.1) then
          CT(ii)='Z'
         End If 
        End Do 
        write(IOut,'(1X,A2,I2,2I5,1X,A1,3X,I2,3X,I2,4X,A1,4X,
     $   4(I3,A1,F8.4,A1,A1))') El,IAt,IBrnch(IAt),
     $   IRing(1,IAt),star,NPi(IAt),NH(IAt),Chir,
     $   (IBond(J,IAt),OpenP,BndOrd(J,IAt),CloseP,CT(J),J=1,NBI) 
       End If 
      End Do
C Fine della stampa
      write(IOut,'('' '')')
      Do ii=1,NAtoms
       If(ISotop(ii).eq.0) CYCLE
       write(IOut,'('' Isotope'',I3,'' for Atom'',I3)')Isotop(ii),ii
      End Do
      write(IOut,'('' '')')
      If(NSter2.gt.0) then
       write(IOut,'(I3,'' E/Z Isomers'')') NSter2
       do IEZ=1,NSter2
        If(ISter2(3,IEZ).eq.-1) then
         write(IOut,'(''Atoms'',I3,'' and'',I3,'' are in Z arrangement''
     $     )') ISter2(1,IEZ),ISter2(2,IEZ)
        ElseIf(ISter2(3,IEZ).eq.1) then
         write(IOut,'(''Atoms'',I3,'' and'',I3,'' are in E arrangement''
     $     )') ISter2(1,IEZ),ISter2(2,IEZ)
        EndIf
       End Do
       write(IOut,'('' '')')
      EndIf
      Return
      End
*Deck AddHyd
      Subroutine AddHyd(IOut,MxBnd,NHeavy,NAtoms,IAn,NPi,NBond,
     $  IBond,BndOrd,NH)
      Implicit Integer (A-Z) 
      Dimension NBond(*),IAn(*),NPi(*),NH(*),IBond(MxBnd,*)
      Real*8 BndOrd(MxBnd,*)
      IHAt = NHeavy
C Add Explicit Hydrogens
      DO IAt=1,NHeavy
       If(NH(IAt).gt.0) then
        write(IOut,'('' Atom'',I3,'' has'',I2,'' explicit hydrogens'')')
     $   IAt,NH(IAt)
        NHI=NH(IAt)
        Do IH=1,NHI
         IHAt=IHAt+1
         IAn(IHAt)=1
         NBond(IAt)=NBond(IAt)+1
         NBI=NBond(IAt)
         IBond(NBI,IAt)=IHAt
         BndOrd(NBI,IAt)=1.0d0
         NBond(IHAt)=1
         IBond(1,IHAt)=IAt
         BndOrd(1,IHAt)=1.0D0
        End Do
       End If
      End Do 
C Add Implicit Hydrogens
      Do 10 IAt=1,NHeavy
       If(IAn(IAt).eq.1) goto 10
       MissH=NFree(IOut,IAt,IAn(IAt),NBond(IAt),NPi(IAt))  
       If(NBond(IAt)+MissH.gt.MxBnd) then
        write(IOut,'('' Too many Bonds for atom'',I3)') IAt
        stop
       End If
       NH(IAt)=NH(IAt)+MissH
       IF(MissH.gt.0) then 
        Do 20 JAt = 1, MissH
         IHAt = IHAt + 1
         IAn(IHAt)=1
         NH(IHAt)=0
         NBond(IAt)=NBond(IAt)+1
         NBond(IHAt)=1
         IBond(NBond(IAt),IAt)=IHAt 
         BndOrd(NBond(IAt),IAt)=1.0d0
         IBond(1,IHAt)=IAt
         BndOrd(1,IHAt)=1.0d0
  20    Continue
       EndIf 
  10  Continue
      NAtoms=IHAt
      Return
      End
*Deck NFree
      Integer Function NFree(IOut,IAt,IZAt,NSigma,NPi)
      Implicit Integer (A-N)
      Dimension IAv(18)
      Data IAv/1,0,1,2,3,4,3,2,1,0,1,2,3,4,3,2,1,0/
      If(IZAt.gt.18) then
       write(Iout,'('' Atom'',I3,'' has a not allowed Atomic Number:'',
     $   I2)')IAt,IZAt 
       Stop
      ElseIf(IZAt.gt.0) then
       NAv=IAv(IZAt)
C      write(IOut,'('' Atom'',I3,'' with atomic number'',I2,
C    $   '' has NAv='',I2,'' Nsigma='',I2,'' NPi='',I2)') 
C    $   IAt,IZAt,NAv,NSigma,NPi
      Else
       write(Iout,'('' Atom'',I3,'' has a not allowed Atomic Number:'',
     $   I2)')IAt,IZAt 
       Stop
      EndIf
      NFree=NAv-NSigma-NPi
      If(NFree.lt.0) then
       write(IOut,'('' Atom'',I3,'' with Atomic Number'',I2,'' has'',I2,
     $ '' sigma and'',I2,'' pi bonds'')') IAt,IZAt,NSigma,NPi  
       Stop
      EndIf
      Return
      End
*Deck IAnEl2
      Character*2 Function IAnEl2(IAn)
      Implicit Integer(A-Z)
C
C     Returns the atomic symbol for a given atomic number.
C
C Input:
C     IAn    : Atomic Number
C Output:
C     IAnEl2 : Atomic Symbol (A2)
      Logical NoCase
C     Local
      Integer i
      Character AllSmb*222
C
        AllSmb = ' ' // 'H HeLiBeB C N O F NeNaMgAlSiP S ClArK CaScTiV '
     $  //'CrMnFeCoNiCuZnGaGeAsSeBrKrRbSrY ZrNbMoTcRuRhPdAgCdInSnSbTeI '
     $  //'XeCsBaLaCePrNdPmSmEuGdTbDyHoErTmYbLuHfTaW ReOsIrPtAuHgTlPbBi'
     $  //'PoAtRnFrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLrRfDbSgBhHsMtDs'
      Ini=2*IAn
      IEnd=Ini+1
      If(Ini.le.0) then
       IAnEl2='  '
      Else
       IAnEl2=AllSmb(Ini:IEnd)
      EndIf
      Return
      End
*Deck El2IAN
      Integer Function El2IAN(AtSymb)
      Implicit Integer(A-Z)
C
C     Returns the atomic number for a given atomic symbol.
C
C Input:
C     AtSymb : Atomic symbol (character)
C     
C     Input
      Character AtSymb*2
C     Local
      Integer i
      Character AllSmb*222, ASymb*2
C     
        AllSmb = ' ' // 'H HeLiBeB C N O F NeNaMgAlSiP S ClArK CaScTiV '
     $  //'CrMnFeCoNiCuZnGaGeAsSeBrKrRbSrY ZrNbMoTcRuRhPdAgCdInSnSbTeI '
     $  //'XeCsBaLaCePrNdPmSmEuGdTbDyHoErTmYbLuHfTaW ReOsIrPtAuHgTlPbBi'
     $  //'PoAtRnFrRaAcThPaU NpPuAmCmBkCfEsFmMdNoLrRfDbSgBhHsMtDs'
        ASymb = AtSymb

      i = Index(AllSmb,ASymb)
      If(i.eq.0) then
        El2IAN = -1
      else
        El2IAN = i/2
      endIf
      Return
      End
*Deck SetIZ
      Subroutine SetIZ(IOut,MxBnd,NAtoms,IAn,IBrnch,NH,NBond,
     $  IBond,IZ)
      Implicit None
      Integer IOut,MxBnd,NAtoms,IAt
      Integer IAn(*),IBrnch(*),NH(*),NBond(*),IBond(MxBnd,*),IZ(4,*)
C Local
      Integer IAnI,NBI,II,IBrI,JAt,NBJ,JJ,IBrJ,KK,NBK,MxBAt,MxSBAt
      Integer MxDBAt,KAt,LAt,LAtOK,ICase
      Integer IUsed(1000)
      Call IClear(NAtoms,IUsed) 
C First Atom
      IAt=1
      IZ(1,IAt)=2
      IZ(2,IAt)=0
      IZ(3,IAt)=0
      IZ(4,IAt)=0
      If(NAtoms.eq.1) Return
C Second Atom
      IAt=2
      IZ(1,IAt)=1
      IZ(2,IAt)=0
      IZ(3,IAt)=0
      IZ(4,IAt)=0
      IUsed(IZ(1,IAt))=1
      If(NAtoms.eq.2) Return
C Third Atom
      IAt=3
      IZ(3,IAt)=0
      IZ(4,IAt)=0
      NBI=NBond(IAt)
      Do II=1,NBI
       If(IBond(II,IAt).eq.1) then
        IZ(1,IAt)=1
        IZ(2,IAt)=2
        IZ(2,1)=0
        IZ(2,2)=3
        GoTo 10
       EndIf
      EndDo
      IZ(1,IAt)=2
      IZ(2,IAt)=1
      IZ(2,1)=3
      IZ(2,2)=0       
   10 Continue
      IUsed(IZ(1,IAt))=1
      If(NAtoms.eq.3) Return
C Fourth Atom
      IAt=4
      NBI=NBond(IAt)
      MxBAt=0
      Do II=1,NBI
       JAt=IBond(II,IAt)
       If(JAt.gt.MxBAt.and.JAt.lt.IAt) MxBAt=JAt 
      EndDo 
      If(MxBAt.eq.0) then
       write(IOut,'('' Wrong Connectivity for Atom'',I3)') IAt
       Stop
      EndIf
      JAt=MxBAt      
      NBJ=NBond(JAt)
      IZ(1,IAt)=JAt
      MxBAt=0
      Do JJ=1,NBJ 
       KAt=IBond(JJ,JAt)
       If(KAt.gt.MxBAt.and.KAt.lt.IAt) MxBAt=KAt
      EndDo 
      If(MxBAt.eq.0) then
       write(IOut,'('' Wrong Connectivity for Atom'',I3)') IAt
       Stop
      EndIf
      KAt=MxBAt 
      IZ(2,IAt)=KAt 
      NBK=NBond(KAt)
      MxBAt=0
      Do kk=1,NBK
       LAt=IBond(KK,KAt)
       If(KAt.gt.MxBAt.and.LAt.ne.JAt.and.LAt.lt.IAt) MxBAt=LAt
      EndDo
      If(MxBAt.ne.0) then
       LAt=MxBAt
       IZ(3,IAt)=LAt 
       IZ(4,IAt)=0
      Else
       NBJ=NBond(JAt)
       LAtOK=0
       Do JJ=1,NBJ 
        LAt=IBond(JJ,JAt)
        IF(LAt.lt.IAt.and.LAt.ne.KAt) LAtOK=LAt
       EndDo 
       If(LAtOK.eq.0) then
        Write(IOut,'('' Bad Bond Pattern for Atom'',I3)') IAt
        Stop
       EndIf 
       IZ(3,IAt)=LAtOK
       IZ(4,IAt)=1
      EndIf
      IUsed(IZ(1,IAt))=1
      If(NAtoms.eq.4) Return
C Other Atoms
      Do 30 IAt=5,NAtoms
       IAnI=IAn(IAt)
       NBI=NBond(IAt)
       IBrI=IBrnch(IAt) 
       MxSBAt=0   
       MxDBat=0
       Do II=1,NBI
        JAt=IBond(II,IAt)
        IBrJ=IBrnch(JAt)
        If(JAt.gt.MxSBAt.and.JAt.lt.IAt.and.IBrI.eq.IBrJ) then
         MxSBAt=JAt
        ElseIf(JAt.gt.MxDBAt.and.JAt.lt.IAt.and.IBrI.ne.IBrJ) then
         MxDBAt=JAt
        EndIf
       EndDo
       MxBAt=Max0(MxSBAt,MxDBAt) 
       If(MxBAt.eq.0) then
        write(IOut,'('' MxBAt=0 for Atom'',I3,'' Program Stops'')')IAt
        STop
       EndIf
C Set IRef=1 for specific cases
C NBond(MxBAt)=2 
C ICase=1
       If(MxBAt.ne.0.and.NBond(MxBAt).eq.2) then
        JAt=MxSBAt
        If(MxSBAt.eq.0) JAt=MxDBAt
        IZ(1,IAt)=JAt
        IZ(2,IAt)=IZ(1,JAt)
        IZ(3,IAt)=IZ(2,JAt)
        IZ(4,IAt)=IUsed(JAt)
        If(IUsed(IZ(1,IAt)).eq.2) IZ(4,IAt)=-1
        ICase=1
        GOTO 25
C IAt,JAt in the Same Branch and IAt not Terminal Atom
C ICase=2
       ElseIf(MxSBAt.gt.0.and.NBond(IAt).ne.1) then
        JAt=MxSBAt
        IZ(1,IAt)=JAt
        IZ(2,IAt)=IZ(1,JAt) 
        IZ(3,IAt)=IZ(2,JAt)
        IZ(4,IAt)=IUsed(JAt)
        If(IUsed(JAt).eq.2) IZ(4,IAt)=-1
        ICase=2
        GOTO 25 
C IAt,JAt in different Branches or IAt terminal Atom
       ElseIf(MxDBAt.gt.0.or.NBond(IAt).eq.1) then
        MxBAt=MxDBAt
        If(MxSBAt.ne.0) MxBAt=MxSBAt
        JAt=MxBAt
        IZ(1,IAt)=JAt 
        KAt=IZ(1,JAt)
        IZ(2,IAt)=KAt
        MxBAt=0
        Do JJ=1,NBJ
         LAt=IBond(JJ,JAt)
      If(IZ(1,LAt).ne.0.and.LAt.lt.IAt.and.LAt.ne.KAt.and.LAt.gt.MxBAt) 
     $    MxBAt=LAt
        EndDo
C Change of Branch or terminal atom and first atom bonded to JAt
        If(MxBAt.eq.0) then
         If(IZ(2,JAt).ne.0) then
          IZ(3,IAt)=IZ(2,JAt) 
          IZ(4,IAt)=IUsed(JAt)
         If(NBond(IAt).eq.1.and.NBond(JAt).gt.2.and.IUsed(JAt).eq.0)then
          IZ(4,IAt)=1
          IUsed(JAt)=1
         EndIf
          If(IUsed(JAt).eq.2) IZ(4,IAt)=-1 
          ICASE=3
          GOTO 25
         Else
          Write(IOut,'('' Wrong IZ for atom'',I3)') IAt
          Stop
         EndIf
C Change of Branch and second atom bonded to JAt
C ICase=4
        ElseIf(MxBAt.gt.0) then
         LAt=MxBAt
         IZ(3,IAt)=LAt 
         IZ(4,IAt)=IUsed(JAt) 
         If(IUsed(JAt).eq.2) IZ(4,IAt)=-1
         ICase=4
         GOTO 25 
        EndIf
       EndIf 
       write(IOut,'('' No Condition fulfilled by Atom'',I3)') IAt
       Stop
   25  IUsed(IZ(1,IAt))=IUsed(IZ(1,IAt))+1
       write(IOut,'('' For Atom'',I3,'' ICase='',I2,''(IZ(4)='',
     $   I2)') IAt,ICase,IZ(4,IAt)
   30 Continue
      Return
      End
*Deck SetBLA
      Subroutine SetBLA(MxBnd,NAtoms,IAn,NPi,IZ,NBond,IBond,BndOrd,BL,
     $  Alpha)
      Implicit Real*8 (A-H,O-Z)
C Set standard bond lengths and valence angles
      Dimension IAn(*),NPi(*),IZ(4,*),NBond(*),IBond(MxBnd,*)
      Dimension BndOrd(MxBnd,*),BL(*),Alpha(*)
      DO IAt = 1, NAtoms
       IAnI=IAn(IAt)
       If(IAt.ge.2) then
        JAt  = IZ(1,IAt)
        do J = 1, NBond(IAt)
         if(IBond(J,IAt).eq.JAt) JRef=J
        enddo
        TestI=BndOrd(JRef,IAt)
        IAnJ=IAn(JAt)
        Val0=RCovCT(IAnI,IAnJ)
        DVal = 0.3D0*log(TestI)
        BL(IAt) = Val0 - DVal
       End If
       if(IAt.ge.3) then
        KAt=IZ(2,IAt)
        IAnK=IAn(KAt)
        do K = 1, NBond(JAt)
         if(IBond(K,JAt).eq.KAt) KRef=K
        enddo
        TestJ=BndOrd(KRef,JAt)
        If(TestI*TestJ.gt.2.8D0) then
         Alpha(IAt)=180.0D0
        ElseIf(TestI*TestJ.gt.1.8D0) then
         Alpha(IAt)=120.0D0
        Else
         Alpha(IAt)=109.4722D0
        EndIf
       EndIf
       If(IAn(IAt).eq.1) then
        NPJ=NPi(JAt)
        If(NPJ.eq.0) then
         Alpha(IAt)=109.4722D0
        ElseIf(NPJ.eq.1) then
         Alpha(IAt)=120.0D0
        ElseIf(NPJ.eq.2) then
         Alpha(IAt)=180.0D0
        EndIf
       EndIf
      End Do
      RETURN
      End
*Deck SetDih
      Subroutine SetDih(IOut,MxBnd,MXatRn,NAtoms,IZ,NBond,IBond,IBrnch,
     $   IRing,NAtRng,IAtRng,Alpha,Beta)
      Implicit None
      Logical Fnd(NAtoms),OKBeta(NAtoms)
      Integer IOut,MxBnd,MxAtRn,NAtoms,IAt,JAt,KAt,LAt,IR,JR,KR,LR,IRI
      Integer IRJ,IRK,IRL,IDiff,ITest,IZI1,IZI2,IZI3,IZI4,IZJ1,IZJ2,IZJ3
      Integer IZJ4,IR1,IR2,IR3,IR4
      Integer IZ(4,*),NBond(*),IBond(MxBnd,*),IBrnch(*)
      Integer IRing(3,*),NAtRng(*),IAtRng(MxAtRn,*)
      Real*8 Alpha(*),Beta(*),SNAt
      Do 10 IAt=1,NAtoms
       JAt=IZ(1,IAt)
       KAt=IZ(2,IAt)
       LAt=IZ(3,IAt)
       IR1=IRing(1,IAt)
       IR2=IRing(1,JAt)
       IR3=IRing(1,KAt)
       IR4=IRing(1,LAt)
       Fnd(IAt)=.False.
       OKBeta(IAt)=.False.
       If(Abs(IZ(4,IAt)).ge.1) then
        If(NBond(JAt).eq.4) then
         Beta(IAt)=109.47122063D0
        ElseIf(NBond(JAt).eq.3) then
         Beta(IAt)=120.0D0
        EndIf
       EndIf
C Manage Cycles
       If(IR1.eq.0.or.IR2.eq.0.or.IR3.eq.0) goto 10
C Valence angle
       SNat=Float(NAtRng(IR1)+NAtRng(IR2)+NAtRng(IR3))/3.0D0
       Alpha(IAt)=180.0*(SNat-2.0d0)/SNat
       If(IR4.eq.0) goto 10
       OKBeta(IAt)=.True.
       IDiff=4
       Do IR=1,3
        IRI=IRing(IR,IAt)
        If(IRI.eq.0) Cycle
        Do JR=1,3
         IRJ=IRing(JR,JAt)
         If(IRJ.eq.0) Cycle
         Do KR=1,3
          IRK=IRing(KR,KAt)
          If(IRK.eq.0) Cycle
          Do LR=1,3
           IRL=IRing(LR,LAt)
           If(IRL.eq.0) Cycle
           ITest=Abs(IRI-IRJ)+ABs(IRJ-IRK)+ABs(IRK-IRL) 
           If(ITest.lt.IDiff) then
            IDiff=ITest         
            IR1=IRI
            IR2=IRJ
            IR3=IRK
            IR4=IRL
           EndIf
          EndDo
         EndDo
        EndDo
       EndDo 
C Dihedral Angle
       IF(IDiff.eq.1.or.IDiff.eq.3) then
        If(IZ(4,IAt).eq.0) Beta(IAt)=180.0D0
       ElseIf(Idiff.eq.0.or.IDiff.eq.2) then
        If(IZ(4,IAt).eq.0) Beta(IAt)=0.0D0
       EndIf 
       write(IOut,'('' Atom'',I3,'' Ring:'',4I2,'' IRTot:'',I2)')
     $  IAt,IR1,IR2,IR3,IR4,IDiff
   10 Continue
      write(IOut,'('' '')')
      IF(NAtoms.lt.5) Return
      Do 20 IAt=5,NAtoms
       If(OKBeta(IAt)) goto 20
       IZI1=IZ(1,IAt)
       IZI2=IZ(2,IAt)
       IZI3=IZ(3,IAt)
       IZI4=IZ(4,IAt)
       Do 30 JAt=4,IAt-1
        IZJ1=IZ(1,JAt)
        IZJ2=IZ(2,JAt)
        IZJ3=IZ(3,JAt)
        IZJ4=IZ(4,JAt)
        If(IZI1.eq.IZJ1.and.IZI2.eq.IZJ2.and.IZI3.eq.IZJ3) then
         If(NBond(IZI1).eq.4) then
          If(Abs(IZI4).eq.1) then
           Beta(IAt)=Alpha(IAt)
           IZ(4,IAt)=-1
          ElseIf(Fnd(JAt)) then
           Beta(IAt)=Beta(JAt)-120.0D0
          Else
           Beta(IAt)=Beta(JAt)+120.0D0
          EndIf
         ElseIf(NBond(IZI1).eq.3) then
          If(Abs(IZI4).eq.1) then
           Beta(IAt)=Alpha(IAt)
           IZ(4,IAt)=-1
          Else
           Beta(IAt)=Beta(JAt)+180.0d0
          EndIf
         EndIf
         Fnd(JAt)=.true.
        EndIf
   30  Continue
       If(IZI1.lt.3.and.IRing(1,IZI2).ne.0.and.IRing(1,IZI3).ne.0) then
        If(IRing(1,IZI2).ne.IRing(1,IZI1)) then
         If(IZ(4,IAt).eq.0) Beta(IAt)=Beta(IAt)+180.0D0
        End If
       EndIf 
       If(Beta(IAt).ge.360.0d0) Beta(IAt)=Beta(IAt)-360.0D0
   20 Continue
      Return
      End
*Deck Aclear
      Subroutine AClear(N,A)
      Implicit Real*8 (A-H,O-Z)
C
C     Clear N elements of A
C
      Dimension A(*)
      Zero =0.0D0
      Do 10 I=1,N
       A(I)=Zero
   10 continue
      Return
      End
*Deck IClear
      Subroutine IClear(N,IA)
      Integer N,IA,I
C
C     Clear N elements of IA
C     
      Dimension IA(*)
      Do 10 I=1,N
       IA(I)=0
  10  Continue
      Return
      End
*Deck RCovCT 
      Function RCovCT(IA,JA)
      Implicit Real*8(A-H,O-Z)
C     
C     Return the bond distance in Angstroms between elements IA and JA
C     using radii from M. Mantina, R. Valero, C. J. Cramer, and
C     D. G. Truhlar, listed in the CRC Handbook.
C         
      Parameter (MaxAn=118)
      Dimension Rad(0:MaxAn)
      Save Rad
      Data Rad/0.0d0,0.32D0,0.37D0,1.30D0,0.99D0,0.84D0,0.75D0,0.71D0,
     $  0.64D0,0.60D0,0.62D0,1.60D0,1.40D0,1.24D0,1.14D0,1.09D0,1.04D0,
     $  1.00D0,1.01D0,2.00D0,1.74D0,1.59D0,1.48D0,1.44D0,1.30D0,1.29D0,
     $  1.24D0,1.18D0,1.17D0,1.22D0,1.20D0,1.23D0,1.20D0,1.20D0,1.18D0,
     $  1.17D0,1.16D0,2.15D0,1.90D0,1.76D0,1.64D0,1.56D0,1.46D0,1.38D0,
     $  1.36D0,1.34D0,1.30D0,1.36D0,1.40D0,1.42D0,1.40D0,1.40D0,1.37D0,
     $  1.36D0,1.36D0,2.38D0,2.06D0,1.94D0,1.84D0,1.90D0,1.88D0,1.86D0,
     $  1.85D0,1.83D0,1.82D0,1.81D0,1.80D0,1.79D0,1.77D0,1.77D0,1.78D0,
     $  1.74D0,1.64D0,1.58D0,1.50D0,1.41D0,1.36D0,1.32D0,1.30D0,1.30D0,
     $  1.32D0,1.44D0,1.45D0,1.50D0,1.42D0,1.48D0,1.46D0,2.42D0,2.11D0,
     $  2.01D0,1.90D0,1.84D0,1.83D0,1.80D0,1.80D0,1.73D0,1.68D0,1.68D0,
     $  1.68D0,1.65D0,1.67D0,1.73D0,1.76D0,1.61D0,1.57D0,1.49D0,1.43D0,
     $  1.41D0,1.34D0,1.29D0,1.28D0,1.21D0,1.22D0,1.36D0,1.43D0,1.62D0,
     $  1.75D0,1.65D0,1.57D0/
C
      IA1 = Min(Max(IA,0),MaxAn)
      JA1 = Min(Max(JA,0),MaxAn)
      RCovCT = Rad(IA1) + Rad(JA1)
      Return
      End
*Deck ZMTPRT
      Subroutine ZMTPRT(IOut,NZ,IANZ,IZ,LBL,LAlpha,LBeta,BL,ALPHA,BETA,
     $  ToAng,ToDeg)
      Implicit Real*8(A-H,O-Z)
C
C     Z-MATRIX PRINTING ROUTINE.
C     CONVERTS FROM INTERNAL (BOHR/RADIAN) UNITS TO EXTERNAL
C     (ANGSTROM/DEGREE) UNITS LOCALLY FOR PRINTING USING ToAng and ToDeg.
C     Modified to avoid FillEl (VB 2025) 
      Logical Coord
      Character*2 IAnEl2,IEl(NZ)
      Parameter (MaxEl=204)
      Dimension IANZ(*), LBL(*), LAlpha(*), LBeta(*)
      Dimension BL(*), ALPHA(*), BETA(*), IZ(4,*)  
      Save Zero, One, F45
      Data Zero/0.0d0/, One/1.0d0/, F45/45.0D0/
 1000 Format(1X,99('-'))
 1010 Format(1X,27X,'Z-MATRIX (ANGSTROMS AND DEGREES)')
 1020 Format('   CD    Cent   Atom    N1       Length/X',
     $  '        N2       Alpha/Y        N3        Beta/Z          J')
 2110 Format(1X,I6,1X,I6,2X,A2)
 2120 Format(1X,I6,9X,A2)
 2210 Format(1X,I6,1X,I6,2X,A2,2X,I6,F11.6,'(',I6,')')
 2220 Format(1X,I6,      9X,A2,2X,I6,F11.6,'(',I6,')')
 2310 Format(1X,I6,1X,I6,2X,A2,2X,I6,F11.6,'(',I6,') ',
     $          I6,1X,F8.3,'(',I6,')')
 2320 Format(1X,I6,      9X,A2,2X,I6,F11.6,'(',I6,') ',
     $          I6,1X,F8.3,'(',I6,')')
 2410 Format(1X,I6,1X,I6,2X,A2,2X,I6,F11.6,'(',I6,') ',
     $          I6,1X,F8.3,'(',I6,') ',I6,1X,F8.3,'(',I6,') ',
     $          I6)
 2420 Format(1X,I6,      9X,A2,2X,I6,F11.6,'(',I6,') ',
     $          I6,1X,F8.3,'(',I6,') ',I6,1X,F8.3,'(',I6,') ',
     $          I6)
 2430 Format(1X,I6,1X,I6,2X,A2,2X,I6,F11.6,7X,F11.6,7X,F11.6)
 2440 Format(1X,I6,1X,6X,2X,A2,2X,I6,F11.6,7X,F11.6,7X,F11.6)
C
C     PRINT THE HEADING.
C
C     TODEG = F45 / ATan(ONE)
      WRITE(IOUT,1000)
      WRITE(IOUT,1010)
      WRITE(IOUT,1020)
      WRITE(IOUT,1000)
C
C     First card.
C
      If(NZ.lt.1) goto 900
      ICard = 1
      Idx   = IAnZ(1)
      IEl(Idx)=IAnEl2(Idx)
      Coord = IZ(1,1).lt.0.or.Bl(1).ne.Zero.or.Alpha(1).ne.Zero.or.
     $  Beta(1).ne.Zero
      PBl = Bl(1) * ToAng 
      PA = Alpha(1) * ToAng
      PB = Beta(1) * ToAng
      ICent = 0
      If(IAnZ(1).ge.0) ICent = 1
      If(IAnZ(1).ge.0.and..not.Coord)
     $  Write(IOut,2110) ICard,ICent,IEl(Idx)
      If(IAnZ(1).lt.0.and..not.Coord) Write(IOut,2120) ICard, IEl(Idx)
      If(IAnZ(1).ge.0.and.Coord)
     $  Write(IOut,2430) ICard, ICent, IEl(Idx), IZ(1,1), PBl, PA, PB
      If(IAnZ(1).lt.0.and.Coord)
     $  Write(IOut,2440) ICard, IEl(Idx), IZ(1,1), PBl, PA, PB
      If(NZ.eq.1) goto 900
C
C     Second card.
C
      NP1   = LBL(2)
      ICard = 2
      Idx   = IAnZ(2)
      IEl(Idx)=IAnEl2(Idx)
      PBl   = BL(2) * ToAng 
      PA = Alpha(2) * ToAng
      PB = Beta(2) * ToAng
      If(IAnZ(2).ge.0) ICent = ICent + 1
      If(IAnZ(2).ge.0.and.IZ(1,2).gt.0)
     $  Write(IOut,2210) ICard,ICent,IEl(Idx),IZ(1,2),PBl,NP1
      If(IAnZ(2).lt.0.and.IZ(1,2).gt.0)
     $  Write(IOut,2220) ICard,IEl(Idx),IZ(1,2),PBl,NP1
      If(IAnZ(2).ge.0.and.IZ(1,2).le.0)
     $  Write(IOut,2430) ICard, ICent, IEl(Idx), IZ(1,2), PBl, PA, PB
      If(IAnZ(2).lt.0.and.IZ(1,2).le.0)
     $  Write(IOut,2440) ICard, IEl(Idx), IZ(1,2), PBl, PA, PB
      If(NZ.eq.2) goto 900
C
C     Third card.
C
      NP1   = LBL(3)
      NP2   = LAlpha(3)
      ICard = 3
      Idx   = IAnZ(3)
      IEl(Idx)=IAnEl2(Idx)
      PBl   = Bl(3) * ToAng
      PA    = Alpha(3) * ToDeg
      If(IZ(1,3).le.0) PA = Alpha(3) * ToAng 
      PB = Beta(3) * Conver
      If(IAnZ(3).ge.0) ICent = ICent + 1
      If(IAnZ(3).ge.0.and.IZ(1,3).gt.0) Write(IOut,2310) ICard,ICent,
     $  IEl(Idx),IZ(1,3),PBl,NP1,IZ(2,3),PA,NP2
      If(IAnZ(3).lt.0.and.IZ(3,1).gt.0) Write(IOut,2320) ICard,
     $  IEl(Idx),IZ(1,3),PBl,NP1,IZ(2,3),PA,NP2
      If(IAnZ(3).ge.0.and.IZ(1,3).le.0)
     $  Write(IOut,2430) ICard, ICent, IEl(Idx), IZ(1,3), PBl, PA, PB
      If(IAnZ(3).lt.0.and.IZ(1,3).le.0)
     $  Write(IOut,2440) ICard, IEl(Idx), IZ(1,3), PBl, PA, PB
      If(NZ.eq.3) goto 900
C
C     Cards 4 through NZ.
C
      Do 500 ICard = 4, NZ
        NP1   = LBL(ICard)
        NP2   = LAlpha(ICard)  
        NP3   = LBeta(ICard)
        Idx   = IAnZ(ICard)
        IEl(Idx)=IAnEl2(Idx)
        PBl   = Bl(ICard) * ToAng
        PA    = Alpha(ICard) * ToDeg
        PB    = Beta(ICard)  * ToDeg
        If(IZ(1,ICard).le.0) PA = Alpha(ICard) * ToAng 
        If(IZ(1,ICard).le.0) PB = Beta(ICard) * ToAng
        If(IAnZ(ICard).ge.0) ICent = ICent + 1
        If(IAnZ(ICard).ge.0.and.IZ(1,ICard).gt.0) Write(IOut,2410)
     $    ICard,ICent,IEl(Idx),IZ(1,ICard),PBl,NP1,IZ(2,ICard),PA,
     $    NP2,IZ(3,ICard),PB,NP3,IZ(4,ICard)
        If(IAnZ(ICard).lt.0.and.IZ(1,ICard).gt.0) Write(IOut,2420)
     $    ICard,IEl(Idx),IZ(1,ICard),PBl,NP1,IZ(2,ICard),PA,
     $    NP2,IZ(3,ICard),PB,NP3,IZ(4,ICard)
        If(IAnZ(ICard).ge.0.and.IZ(1,ICard).le.0)
     $    Write(IOut,2430) ICard, ICent, IEl(Idx), IZ(1,ICard),
     $    PBl, PA, PB
        If(IAnZ(ICard).lt.0.and.IZ(1,ICard).le.0)
     $    Write(IOut,2440) ICard, IEl(Idx), IZ(1,ICard), PBl, PA, PB
  500   Continue
  900 Write(IOut,1000)
      Return
      End
*Deck DisMat
      Subroutine DisMat(NAtoms,IAn,C,IFlag,NColX,IOut,Error,ICrowd,
     $  Conver)
      Implicit Real*8(A-H,O-Z)
C
C     Routine to print a lower triangular matrix of intermolecular
C     distances, given the coordinates.
C
C     NAtoms ... Number of atoms.
C     IAn    ... Integer vector containing the atomic numbers for
C                each of the centers.
C     C      ... Floating point vector containing the coordinates
C                of the NAtoms centers.
C     IFlag  ... Conversion flag.  This variable controls scaling
C                of the distances before printing.
C                1 ... Scale output distances by one.
C                2 ... Scale output distances by Conver.
C                3 ... Scale output distances by one/Conver
C     NColX  ... Number of columns per printed block.
C     IOut   ... Fortran output unit for printing.  If zero, prinitng is
C                suppressed.  If negative, only error messages are printed.
C     Error  ... Logical variable set to true if any pair of atoms is
C                less than 0.5 angstroms (TooClo).
C     ICrowd   ... 0,1 Flag an error on zero distances only
C                  2   Flag an error on distances < 0.5 A
C                  3   Do not flag errors.
C
C     Conver ... A unit conversion factor.
C
      Logical Error, ErrI
      Parameter (NCMax=11, MaxSav=20)
      Real*8 MDCutO
      Character*2 IanEl2
      Dimension C(3,*),S(NCMax),IAn(NAtoms),ISav(MaxSav),JSav(MaxSav)
      Save TooClo, Thresh, Zero, One
      Data TooClo/0.5d0/,Thresh/1.0d-6/,Zero/0.0d0/,One/1.0d0/
 2001 Format(10X,11(5X,I6))
 2003 Format(I6,2X,A2,11F11.6)
 2004 Format(' Small interatomic distances encountered: ',10I6)
 2005 Format(' Small interatomic distances encountered:',
     $  (/,1X,2I5,1PD9.2))
C
      Error = .False.
      If(ICrowd.le.1) then
        Cut = MDCutO(0)
      else if(ICrowd.eq.2) then
        Cut = TooClo
      else
        If(IOut.le.0) Return
        Cut = -Float(1)
        endIf
      Scale = One
      If(IFlag.eq.2) Scale = Conver
      If(IFlag.eq.3) Scale = One / Conver
      If(IOut.le.0) then
        Cut = Cut/Scale
        Do 30 I = 2, NAtoms
          Error = .False.
          If(IAn(I).gt.0) then
            Do 10 J = 1, (I-1)
              R = Sqrt( (C(1,I)-C(1,J))**2 +
     $                  (C(2,I)-C(2,J))**2 +
     $                  (C(3,I)-C(3,J))**2 )
              Error = Error.or.(R.lt.Cut.and.IAn(J).gt.0)
   10         Continue
              If(Error) then
              Do 20 J = 1, (I-1)
                R = Sqrt( (C(1,I)-C(1,J))**2 +
     $                    (C(2,I)-C(2,J))**2 +
     $                    (C(3,I)-C(3,J))**2 )
                If(R.lt.Cut) then
                  If(IOut.ne.0) Write(-IOut,2005) I, J, R*Scale
                  STOP 
                  endIf
   20           Continue
              endIf
            endIf
   30     Continue
      else
        NCol = Min(NColX,NCMax)
        IStart = 1
        Kount = 0
        Do 60 IStart = 1, NAtoms, NCol
          M = 0
          IEnd = Min(IStart+NCol-1,NAtoms)
          Write(IOut,2001) (IR,IR=IStart,IEnd)
          Do 50 I = IStart, NAtoms
            M = M + 1
            IRange = Min(M,NCol)
            Do 40 IR = 1, IRange
              J = IStart + IR - 1
              Temp = (C(1,I)-C(1,J))**2
     $          + (C(2,I)-C(2,J))**2 +(C(3,I)-C(3,J))**2
              S(IR) = Zero
              If(Temp.gt.Thresh) S(IR) = Scale*Sqrt(Temp)
              If(I.ne.J) then
                ErrI = S(IR).lt.Cut.and.IAn(I).gt.0.and.IAn(J).gt.0
                Error = Error.or.ErrI
                If(ErrI.and.Kount.lt.MaxSav) then
                  Kount = Kount + 1
                  ISav(Kount) = I
                  JSav(Kount) = J
                  endIf
                endIf
   40         Continue
            Write(IOut,2003) I,IAnEl2(IAn(I)),(S(IR),IR=1,IRange)
   50       Continue
   60     Continue
        If(Error) Write(IOut,2004) (ISav(I),JSav(I),I=1,Kount)
        endIf
      Return
      End
*Deck HedPrt
      Subroutine HedPrt(IOut,Width,String,Num)
      Implicit Integer(A-Z)
C
C     Print a header for a matrix.  If Width is >0, the
C     header is centered in a line of length Width.  If
C     Width=-1 or -2, then real/imaginary is appended to
C     the label.
C
      Character String*(*), Blank*80, LabRI*20
 1000 Format(A,A,A,I10)
 1010 Format(A,A,A)
 1020 Format(A,A,A,I10,A,':')
 1030 Format(A,A,A,A,':')
C
      Blank = ' '
      LStr = Max(LinEnd(String),1)
      If(Width.le.0) then
        LBl = 1
      else if(Num.gt.0) then
        LBl = (Width-LStr-10)/2 + 1
      else
        LBl = (Width-LStr)/2 + 1
        endIf
      If(Width.eq.-1) then
        LabRI = '(real)'
      else if(Width.eq.-2) then
        LabRI = '(imag)'
      else
        LabRI = ' '
        endIf
      LLabRI = LinEnd(LabRI)
      LBl = Max(LBl,1)
      If(Width.ge.0.and.Num.gt.0) then
        Write(IOut,1000) Blank(1:LBl),String(1:LStr),Blank(1:LBl),Num
      else if(Width.ge.0) then
        Write(IOut,1010) Blank(1:LBl), String(1:LStr), Blank(1:LBl)
      else if(Num.gt.0) then
        Write(IOut,1020) Blank(1:LBl),String(1:LStr),Blank(1:LBl),Num,
     $    LabRI(1:LLabRI)
      else
        Write(IOut,1030) Blank(1:LBl),String(1:LStr),Blank(1:LBl),
     $    LabRI(1:LLabRI)
        endIf
      Return
      End
*Deck LTOutC
      Subroutine LTOutC(IOut,IOpt,NRI,IRI,NSpin,N,A,Key)
      Implicit Real*8(A-H,O-Z)
C
C     Working precision routine to print out the lower triangular part
C     of a matrix:
C
C     IOpt .. Output format:
C             1 ... D13.6
C             2 ... F13.4
C     NRI ... Number of components (1 for real, 2 for complex).
C     IRI ... Which component to print here (1=real, 2=imaginary).
C             Negative for anti-Hermetian matrices.
C     NSpin . -1 for just large component, 1 for just one spin component,
C             -2 for l and s components, 2 for a and b components.
C     N   ... Dimension of matrix.
C     A   ... Array to be printed.
C     Key ... <0: Suppress elements with absolute values less than
C                 10**(-6+Key)
C              1: Print complete matrix.
C
      Parameter (NumCol=5)
      Character*1 LabX(NumCol)
      Dimension A(NRI,*), S(NumCol), IRX(NumCol)
      Save Zero
      Data Zero/0.0d0/
 1010 Format(4X,10(7X,I6,A1))
 1020 Format(I7,A,D13.6,9D14.6)
 1030 Format(I7,A,F13.4,9F14.4)
C
      If(IOpt.lt.1.or.IOpt.gt.2) then
       write(IOut,'(/,''Illegal IOpt in LTOutC.'')')
       Stop
      EndIf
       Call PrtThr(Key,Thresh)
      NDo = IAbs(NSpin)*N
      IRIA = IAbs(IRI)
      Do 100 IStart = 1, NDo, NumCol
       IEnd = Min(IStart+NumCol-1,NDo)
       NCol = IEnd - IStart + 1
       Call LabSpn(NSpin,IStart-1,NCol,IRX,LabX)
       Write(IOut,1010) (IRX(I),LabX(I),I=1,NCol)
       Do 200 IRow = IStart, NDo
        ILim = Min(IRow-IStart+1,NumCol)
        L = (IRow*(IRow-1))/2 + IStart
        Do 10 I = 1, ILim
         S(I) = A(IRIA,L+I-1)
         If(Key.ne.1.and.Abs(S(I)).lt.Thresh) S(I) = Zero
   10   Continue
        If((IRI.eq.-1.or.IRI.eq.2)) then
         ILim1 = Min(ILim,IRow-IStart)
         Do 20 I = 1, ILim1
          S(I) = -S(I)
   20    Continue
        endIf
        Call LabSpn(NSpin,IRow-1,1,IRX,LabX)
        If(IOpt.eq.1) then
         Write(IOut,1020) IRX(1), LabX(1), (S(I),I=1,ILim)
        else
         Write(IOut,1030) IRX(1), LabX(1), (S(I),I=1,ILim)
        endIf
  200  Continue
  100 Continue
      Return
      End
*Deck LTOut
      Subroutine LTOut(IOut,N,A,Key)
      Implicit Real*8(A-H,O-Z)
C
C     Working precision routine to print out the lower triangular part
C     of a symmetric matrix stored in compressed lower triangular form.
C
C        N         Dimension of matrix.
C        A         Array to be printed.
C        KEY    N<=0 ... Suppress elements with absolute values less
C                        than 10**(-6+N)
C                  1 ... Print complete matrix.
C
      Dimension A(*)
C
      Call LTOutC(IOut,1,1,1,1,N,A,Key)
      Return
      End
*Deck LTOutS
      Subroutine LTOutS(IOut,String,Num,N,A,Key)
      Implicit Real*8(A-H,O-Z)
C
C     Print a heading, followed by a lower triangular matrix.
C     If Num is positive, it is also printed.
C
      Character*(*) String
      Dimension A(*)
C
      Call HedPrt(IOut,0,String,Num)
      Call LTOutC(IOut,1,1,1,1,N,A,Key)
      Return
      End
*Deck LabSpn
      Subroutine LabSpn(NSpin,IOff,N,IRX,LabX)
      Implicit Integer(A-Z)
C
C     Generate numbers and labels for a set of functions for printing.
C
      Dimension IRX(*)
      Character*1 LabX(*), LabLS(0:1), LabAB(0:1)
      Save LabLS, LabAB
      Data LabLS/'l','s'/, LabAB/'a','b'/
C
      If(NSpin.lt.0) then
        Do 10 I = 1, N
          IP = IOff + I
          IB = (IP-1) / N
          IRX(I) = IP - N*IB
         LabX(I) = LabLS(IB)
   10   Continue
      else if(NSpin.eq.1) then
        Do 20 I = 1, N
          IRX(I) = IOff + I
          LabX(I) = ' '
   20   Continue
      else
        Do 30 I = 1, N
          IP = IOff + I
          IB = 1 - Mod(IP,2)
          IRX(I) = (IP+1) / 2
          LabX(I) = LabAB(IB)
   30    Continue
        endIf
      Return
      End
*Deck PrtThr
      Subroutine PrtThr(Key,Thresh)
      Implicit Real*8(A-H,O-Z)
C
C     Decode the Key argument to print routines and return the threshold.
C
      Save Zero, Ten
      Data Zero/0.0d0/, Ten/10.0d0/
C
      If(Key.gt.0) then
        Thresh = Zero
      else
        Thresh = Ten**(Key-6)
        endIf
      Return
      End
*Deck AMove
      Subroutine AMove(N,A,B)
      Implicit Real*8(A-H,O-Z)
C
C     Move N words from A to B.  Do not use this routine for
C     overlapping arrays!
C
      Common/IO/in,iout,ipunch
      Dimension A(*), B(*)
C
        If(.not.(N.le.0.or.Loc(B(N)).lt.Loc(A(1)).or.
     $  Loc(B(1)).gt.Loc(A(N)).or.
     $  (Loc(B(1)).eq.Loc(A(1)).and.Loc(B(N)).eq.Loc(A(N))))) then
         write(IOut,'('' Overlap in AMove'')')
         Stop
        endif
      Do 10 I = 1, N
       B(I) = A(I)
   10 continue
      Return
      End
*Deck IMove
      Subroutine IMove(N,A,B)
      Implicit Integer(A-Z)
C
C     Copy N elements from A to B.  Do not use this routine for
C     overlapping arrays!
C
      Dimension A(*), B(*)
       Do 10 I = 1, N
        B(I) = A(I)
   10  Continue  
      Return
      End
*Deck MDCutO
      Function MDCutO(I)
      Implicit Real*8(A-H,O-Z)
      Real*8 MDCutO
C     
C     Return a machine-dependent cutoff.
C     I = 0 ... A value which should be considered small compared to
C               1.0, taking into account round-off error.
C     I = 1 ... A value suitable for a general test for zero values.
C     I = 2 ... A value suitable for thresholding in rotation from
C               local atomic axes.
C     
      If(I.eq.0) then
       MDCutO = 1.d-12
      else if(I.eq.1) then
       MDCutO = 1.d-30
      else if(I.eq.2) then 
       MDCutO = 1.d-6
      endIf
      Return
      End
*Deck LinEnd
      integer function linend(cline)
      implicit real*8(a-h,o-z)
c
c     function which returns the length of a character string,
c     excluding trailing blanks, tabs, nulls, and carriage returns.
c      
      character*(*) cline, blank, null*1, tab*1, cr*1
      parameter (blank=' ')
c     
      null = char(0)
      tab = char(9)
      cr = char(13)
      Do 5 i = len(cline), 1, -1
        if(cline(i:i).ne.blank.and.cline(i:i).ne.cr.and.
     $    cline(i:i).ne.null.and.
     $    cline(i:i).ne.tab) goto 10
    5     continue
      linend = 0
      return
   10 linend = i
      return
      end
*Deck LlinCl
      Subroutine LLinCl(cline)
      implicit real*8(a-h,o-z)
c
c     subroutine for clearing a string
c     
      character*(*) cline, blank
      parameter (blank=' ')
c     
      do 10 i = 1, len(cline), 1
       cline(i:i)=blank
   10  continue
      return
      end
*Deck ZToC
      Subroutine ZToC(MaxNZ,NZ,IAnZ,IZ,Bl,Alph,Bet,TTest,NAtoms,IAn,C,
     $  CZ,A,B,D,Alpha,Beta,IOut,Error,TstAng)
      Implicit Real*8(A-H,O-Z)
C1ZToC
C
C SYNOPSIS:
C   ZToC(MaxNZ,NZ,IAnZ,IZ,Bl,Alph,Bet,TTest,NAtoms,IAn,C,CZ,A,B,D,
C   Alpha,Beta,IOut,Error,TstAng)
C
C DESCRIPTION:
C   This subroutine computes the cartesian coordinates, given the
C   Z-matrix. This routine returns coordinates both with and without
C   the dummy atoms.	
C
C VARIABLES:
C   MaxNZ  ... Maximum number of lines in Z-matrix.
C   NZ     ... Number of lines in the Z-matrix.
C   IAnZ   ... The atomic numbers of the Z-matrix centers.
C   IZ     ... The integer components of the Z-matrix.
C   Bl     ... The bond-lengths from the Z-matrix.
C   Alph   ... The bond-angles from the Z-matrix.
C   Bet    ... The dihedral angles from the Z-matrix.
C   TTest  ... Logical flag to enable testing for tetrahedral angles.
C	       This feature is useful in obtaining exact tetrahedral
C	       angles.	If any are found and this flag is set, the exact
C	       value is used.
C   NAtoms ... Number of atoms (dummies removed), returned.
C   IAn    ... Atomic numbers (dummies removed), returned.
C   C      ... Cartesian Coordinates (dummies removed), returned.
C   CZ     ... Cartesian Coordinates with dummies, returned.
C   A      ... Scratch vector of length NZ.
C   B      ... Scratch vector of length NZ.
C   D      ... Scratch vector of length NZ.
C   Alpha  ... Scratch vector of length NZ.
C   Beta   ... Scratch vector of length NZ.
C   IOut   ... Output unit for error print.  If IOut=0, no print is done.
C   Error  ... A logical variable set to true if ZToC is unable to complete
C              its task.
C   TstAng ... Whether to require angles to be .ge.0.and.lt.180.
C
C2HISTORY
C   Updated error messages - J. L. Sonnenberg 7/2009
C?
      Logical TTest,Error,TstAng,TetRnd,NonLin
      Dimension IAnZ(*),IZ(4,*),Bl(*),Alpha(*),Beta(*),IAN(*),
     $  C(3,*),CZ(3,*),A(*),B(*),D(*),Alph(*),Bet(*),U1(3),U2(3),U3(3),
     $  U4(3),VJ(3),VP(3),V3(3)
      Save Zero,One,Two,CutOff,Four,RTo0
      Data Zero,One,Two,CutOff/0.D0,1.D0,2.D0,1.D-08/
      Data Four,RTo0/4.D0,1.D-12/
 1000 Format(1X,I7,' Z-matrix lines is greater than the maximum of',I7,
     $  ' in subroutine ZToC.')
 1010 Format(' Error on Z-matrix line number',I7,': ',
     $  ' invalid Beta angle type',I5,'.')
 1020 Format(' Error on Z-matrix line number',I7,': ',
     $  ' reference made to an undefined center.')
 1030 Format(' Error on Z-matrix line number',I7,': ',
     $  ' multiple references to a center on the same line.')
 1041 Format(' Error on Z-matrix line number',I7,': ',
     $      ' Atom',I7,' too close to atom',I7)
 1042 Format(' Error on Z-matrix line number',I7,': ',
     $  ' Normalization problem. Are the reference atoms co-linear?')
 1043 Format(' Error on Z-matrix line number',I7,': ',
     $       ' Atoms',3I7,' are colinear')
 1044 Format(' Error on Z-matrix line number',I7,': ',
     $       ' No solution for these values of Alpha and Beta')
 1050 Format(' Error on Z-matrix line number',I7,': ',
     $  ' angle Alpha is outside the valid range of 0 to 180.')
 1060 Format(' Bond length on Z-matrix line number',I7,
     $       ' is not positive.')
 1070 Format(' Error on Z-matrix line number',I7,': ',
     $  ' angle Beta is outside the valid range of 0 to 180.')
 1080 Format(1X,I5,' tetrahedral angles replaced.')
C
C     Check for potential overflow.
C
      Error = NZ.gt.MaxNZ
      If(Error) then
        If(IOut.ne.0) Write(IOut,1000) NZ, MaxNZ
        Return
        endIf
C
C     Check for nonsense in the connectivity.
C
      Do 13 I = 2, NZ
        If(IZ(1,I).le.0.and.IZ(2,I).eq.0.and.IZ(3,I).eq.0.and.
     $     IZ(4,I).eq.0) goto 13
        If(IZ(1,I).ge.I.or.IZ(2,I).ge.I.or.IZ(3,I).ge.I.or.
     $     IZ(1,I).le.0) then
          Error = .True.
          If(IOut.ne.0) Write(IOut,1020) I
          endIf
        If(I.eq.3) then
          If(IZ(2,I).le.0) then
            Error = .True.
            If(IOut.ne.0) Write(IOut,1020) I
            endIf
          If(IZ(1,I).eq.IZ(2,I)) then
            Error = .True.
            If(IOut.ne.0) Write(IOut,1030) I
            endIf
        else if(I.ge.4) then
          If(IAbs(IZ(4,I)).gt.2) then
            Error = .true.
            If(IOut.ne.0) Write(IOut,1010) I, IZ(4,I)
            endIf
          If(IZ(3,I).le.0) then
            Error = .True.
            If(IOut.ne.0) Write(IOut,1020) I
            endIf
          If((IZ(1,I).eq.IZ(2,I)).or.(IZ(1,I).eq.IZ(3,I)).or.
     $       (IZ(2,I).eq.IZ(3,I))) then
            Error = .True.
            If(IOut.ne.0) Write(IOut,1030) I
            endIf
          endIf
   13   Continue
      If(Error) Return
      Pi = Four*ATan(One)
      Call AClear(3*NZ,CZ)
C
C     Move angles to local arrays and optionally test for tetrahedral
C     angles.  If TstAng, test Alpha for out of range 0 to 180 degrees.
C     Test for negative bond lengths.
C
      NumTet = 0
      Do 20 I = 1, NZ
        Alpha(I) = Alph(I)
        Beta(I)  = Bet(I)
        If(IZ(1,I).gt.0) then
          If(I.gt.1.and.Bl(I).le.Zero) then
            Error = .True.
            If(IOut.ne.0) Write(IOut,1060) I
            endIf
          If(I.gt.2.and.TstAng.and.
     $      (Alpha(I).le.Zero.or.Alpha(I).ge.Pi)) then
            Error = .True.
            If(IOut.ne.0) Write(IOut,1050) I
            endIf
          If(TTest) then
            If(TetRnd(Alpha(I))) then
              Alph(I) = Alpha(I)
              NumTet = NumTet + 1
              endIf
            If(TetRnd(Beta(I))) then
              Bet(I)  = Beta(I)
              NumTet = NumTet + 1
              endIf
            endIf
          If(I.gt.3.and.ABs(IZ(4,I)).eq.1.and.
     $      (Beta(I).le.Zero.or.Beta(I).ge.Pi)) then
            Error = .True.
            If(IOut.ne.0) Write(IOut,1070) I
            endIf
          If(I.gt.3.and.Abs(IZ(4,I)).eq.2.and.
     $      (Beta(I).lt.Zero.or.Beta(I).ge.Pi)) then
            Error = .True.
            If(IOut.ne.0) Write(IOut,1070) I
            endIf
          endIf
   20   Continue
      If(NumTet.ne.0.and.IOut.ne.0) Write(IOut,1080) NumTet
      If(Error) Return 
C
C     Atom 1.
C
      If(IZ(1,1).le.0) then
        CZ(1,1) = Bl(1)
        CZ(2,1) = Alpha(1)
        CZ(3,1) = Beta(1)
        endIf
C
C     Atom 2.
C
      If(NZ.ge.2) then
        If(IZ(1,2).le.0) then
          CZ(1,2) = Bl(2)
          CZ(2,2) = Alpha(2)
          CZ(3,2) = Beta(2)
        else
          CZ(1,2) = CZ(1,1)
          CZ(2,2) = CZ(2,1)
          CZ(3,2) = CZ(3,1) + Bl(2)
          endIf
        Call ASUnit(3,CZ(1,2),CZ(1,1),U2)
        endIf
C
C     Atom 3  Put in the plane defined by atoms 1 and 2 and
C     the X-axis.  If this is impossible, use the Z-axis.
C     Be careful of placing atoms until at least one is off-axis.
C
      Do 80 I = 3, NZ
        If(IZ(1,I).le.0) then
          CZ(1,I) = Bl(I)
          CZ(2,I) = Alpha(I)
          CZ(3,I) = Beta(I)
          Call ASUnit(3,CZ(1,I),CZ(1,1),U1)
          NonLin = (One-Abs(SProd(3,U1,U2))).gt.CutOff
        else
          If(IZ(1,I).eq.1) then
            Call ASUnit(3,CZ(1,IZ(1,I)),CZ(1,2),U1)
          else
            Call ASUnit(3,CZ(1,IZ(1,I)),CZ(1,1),U1)
            endIf
          Call Place(CutOff,U1,CZ(1,IZ(1,I)),Bl(I),Alpha(I),CZ(1,I),
     $      NonLin)
          endIf
        If(NonLin) goto 90
   80   Continue
      I = NZ
C
C     General Nth atom for non-linear molecule.
C
   90 K = I + 1
      Do 250 J = K, NZ
        DCAJ = Cos(Alpha(J))
        DSAJ = Sin(Alpha(J))
        DCBJ = Cos(Beta(J))
        DSBJ = Sin(Beta(J))
C
C       Cartesian coordinates
        If(IZ(1,J).le.0) then
          CZ(1,J) = Bl(J)
          CZ(2,J) = Alpha(J)
          CZ(3,J) = Beta(J)
C
C       Dihedral angle
        else if(IZ(4,J).eq.0) then
          Call Vec(CutOff,Error,U1,CZ,IZ(2,J),IZ(3,J))
          If(Error) then
            If(IOut.ne.0) Write(IOut,1041) J,IZ(2,J),IZ(3,J)
            Return
            endIf
          Call Vec(CutOff,Error,U2,CZ,IZ(1,J),IZ(2,J))
          If(Error) then
            If(IOut.ne.0) Write(IOut,1041) J,IZ(1,J),IZ(2,J)
            Return
            endIf
          Call VProd(VP,U1,U2)
          Arg = One - SProd(3,U1,U2)**2
          Error = Arg.lt.Zero
          If(Error) then
            If(IOut.ne.0) Write(IOut,1042) J
            Return
            endIf
          R = Sqrt(Arg)
          Error = R.lt.CutOff
          If(Error) then
            If(IOut.ne.0) Write(IOut,1043) J,(IZ(I,J),I=1,3)
            Return
            endIf
          Call AScale(3,(One/R),VP,U3)
          Call VProd(U4,U3,U2)
          Do 130 I = 1, 3
            VJ(I) = (-U2(I)*DCAJ+U4(I)*DSAJ*DCBJ+U3(I)*DSAJ*DSBJ)*Bl(J)
            CZ(I,J) = VJ(I) + CZ(I,IZ(1,J))
  130     Continue
C
C       Second bond angle specified.
        else if(IAbs(IZ(4,J)).eq.1) then
          Call Vec(CutOff,Error,U1,CZ,IZ(1,J),IZ(3,J))
          If(Error) then
            If(IOut.ne.0) Write(IOut,1041) J,IZ(1,J),IZ(3,J)
            Return
            endIf
          Call Vec(CutOff,Error,U2,CZ,IZ(2,J),IZ(1,J))
          If(Error) then
            If(IOut.ne.0) Write(IOut,1041) J,IZ(2,J),IZ(1,J)
            Return
            endIf
          Zeta = -SProd(3,U1,U2)
          Denom = One - Zeta** 2
          Error = Denom.lt.CutOff
          If(Error) then
            If(IOut.ne.0) Write(IOut,1043) J,(IZ(I,J),I=1,3)
            Return
            endIf
          A(J) = (-DCBJ+Zeta*DCAJ)/Denom
          B(J) = (DCAJ-Zeta*DCBJ)/Denom
          R = Zero
          Gamma = Pi / Two
          If(Abs(Zeta).ge.CutOff) then
            If(Zeta.lt.Zero) R = Pi
            Error = Denom.lt.Zero
            If(Error) then
              If(IOut.ne.0) Write(IOut,1043) J,(IZ(I,J),I=1,3)
              Return
              endIf
            Gamma = ATan(Sqrt(Denom)/Zeta) + R
            endIf
          D(J) = Zero
          If(Abs(Gamma+Alpha(J)+Beta(J)-Two*PI).ge.CutOff) then
            Arg = (One + A(J)*DCBJ - B(J)*DCAJ)  /  Denom
            Error = Arg.lt.Zero
            If(Error) then
              If(IOut.ne.0) Write(IOut,1044) J
              Return
              endIf
            D(J) = Float(IZ(4,J)) * Sqrt(Arg)
            endIf
          Call VProd(V3,U1,U2)
          Do 200 I = 1, 3
            U3(I) = A(J)*U1(I)+B(J)*U2(I)+D(J)*V3(I)
            VJ(I) = Bl(J)*U3(I)
            CZ(I,J) = VJ(I) + CZ(I,IZ(1,J))
  200     Continue
C
C       Angle to plane specified.
        else if(IAbs(IZ(4,J)).eq.2) then
          Call Vec(CutOff,Error,U1,CZ,IZ(1,J),IZ(3,J))
          If(Error) then
            If(IOut.ne.0) Write(IOut,1041) J,IZ(1,J),IZ(3,J)
            Return
            endIf
          Call Vec(CutOff,Error,U2,CZ,IZ(2,J),IZ(1,J))
          If(Error) then
            If(IOut.ne.0) Write(IOut,1041) J,IZ(2,J),IZ(1,J)
            Return
            endIf
          Zeta = -SProd(3,U1,U2)
          Call VProd(V3,U1,U2)
          V3MAG = Sqrt(V3(1)*V3(1)+V3(2)*V3(2)+V3(3)*V3(3))
          Denom = One - Zeta**2
          Error = Abs(Denom).lt.CutOff
          If(Error) then
            If(IOut.ne.0) Write(IOut,1042) J
            Return
            endIf
          A(J) = V3MAG*DCBJ / Denom
          Arg = (One-DCAJ*DCAJ-A(J)*DCBJ*V3MAG) / Denom
          Error = Arg.lt.Zero
          If(Error) then
            If(IOut.ne.0) Write(IOut,1044) J
            Return
            endIf
          B(J) = Sqrt(Arg)
          If(IZ(4,J).ne.2) B(J) = -B(J)
          D(J) = B(J)*Zeta+DCAJ
          Do 240 I = 1, 3
            U3(I) = B(J)*U1(I)+D(J)*U2(I)+A(J)*V3(I)
            VJ(I) = Bl(J)*U3(I)
            CZ(I,J) = VJ(I) + CZ(I,IZ(1,J))
  240     Continue 
          endIf
  250   Continue
C
C     Eliminate dummy atoms, which are characterized by negative atomic
C     numbers.  Ghost atoms (zero atomic number) are retained.  Tidy up
C     the coordinates by zeroing elements less than RTo0.
C
      NAtoms = 0
      Do 280 I = 1, NZ
        If(IAnZ(I).ne.-1) then
          NAtoms = NAtoms + 1
          IAn(NAtoms) = IAnZ(I)
          C(1,NAtoms) = CZ(1,I)
          C(2,NAtoms) = CZ(2,I)
          C(3,NAtoms) = CZ(3,I)
          endIf
  280   Continue
      I = IPopVc(3*NAtoms,C,RTo0,C)
      Return
      End
*Deck SProd
      Function SProd(N,A,B)
      Implicit Real*8(A-H,O-Z)
C     
C     This function returns the scalar product of vectors A and B.
C     Same as SProd but a separate GPU kernel.
C       
      Parameter (Zero=0.0D0)
      Dimension A(*), B(*)
C       
      Prod = Zero
      Do 10 I = 1, N
       Prod = Prod + A(I) * B(I)
   10 Continue
      SProd = Prod
      Return
      End
*Deck VEC
      SUBROUTINE VEC(SMALL,OHOH,U,C,J,K)
      Implicit Real*8(A-H,O-Z)
      LOGICAL OHOH
      DIMENSION C(1),R(3),U(3)
      Save ZERO 
      DATA ZERO/0.0D0/
C  
      R2=ZERO
      JTEMP=(J-1)*3
      KTEMP=(K-1)*3
      DO 10 I=1,3
       R(I)=C(I+JTEMP)-C(I+KTEMP)
       R2=R2+R(I)*R(I)
   10 Continue 
      R2=Sqrt(R2) 
      OHOH = R2 .LT. SMALL
      IF (OHOH) RETURN 
      DO 20 I=1,3
       U(I)=R(I)/R2
   20 Continue
      RETURN
      END
*Deck VProd
      Subroutine VProd(VP,X,Y)
      Implicit Real*8(A-H,O-Z)
C
C     VP = X Cross Y
C     
      Dimension VP(3),X(3),Y(3)
C  
      VP(1)=X(2)*Y(3)-X(3)*Y(2)
      VP(2)=X(3)*Y(1)-X(1)*Y(3)
      VP(3)=X(1)*Y(2)-X(2)*Y(1)
      Return
      End
*Deck AScale
      Subroutine AScale(N,S,A,B)
      Implicit Real*8(A-H,O-Z)
C
C     Vector B = S * A
C     
      Dimension A(*),B(*)
      Do 10 I = 1, N
       B(I) = S*A(I)
   10 continue
      Return
      End
*Deck ASub
      Subroutine ASub(N,A,B,C)
      Implicit Real*8(A-H,O-Z)
C
C     Routine to perform the vector operation
C     C = A - B
C     where C, A, and B are vectors of length N.
C      
      Dimension C(*), A(*), B(*)
        Do 10 I = 1, N
         C(I) = A(I) - B(I)
   10   continue
      Return
      End
*Deck ASUNIT
      SUBROUTINE ASUNIT(N,A,B,C)
      Implicit Real*8(A-H,O-Z)
C
C     THIS SUBROUTINE LOADS C WITH A UNIT VECTOR IN THE A-B DIRECTION.
C     
      DIMENSION A(1), B(1), C(1)
      CALL ASUB(N,A,B,C)
      CALL AUNIT(N,C,C)
      RETURN
      END
*Deck AUnit
      Subroutine AUnit(N,A,B)
      Implicit Real*8(A-H,O-Z)
C
C     Load B with a unit vector in the direction of A.
C     A and B can be the same.
C      
      Dimension A(*), B(*)
      Save Zero, One
      Data Zero/0.0d0/, One/1.0d0/
C
      R = Sqrt(SProd(N,A,A))
      If(R.eq.Zero) then
        Call AClear(N,B)
      else
        Call AScale(N,One/R,A,B)
        endIf
      Return
      End
*Deck TetRnd
      Logical Function TetRnd(Angle)
      Implicit Real*8(A-H,O-Z)
C
C     Test Angle for being near to tetrahedral and round it to the
C     exact value.  The value of this function is .true. if rounding
C     was performed. 
C  
      Save One, Three, F45, Tol
      Data One/1.0d0/, Three/3.0d0/, F45/45.0d0/, Tol/1.d-3/
C
      RefAng = ACos(-One/Three)
      Tol1 = Tol * ATan(One) / F45
      If(Abs(Angle-RefAng).lt.Tol1) then
        TetRnd = .True.
        Angle = RefAng
      else if(Abs(Angle+RefAng).lt.Tol1) then
        TetRnd = .True.
        Angle = -RefAng
      else
        TetRnd = .False.
        endIf
      Return
      End
*Deck IPopVC
      Function IPopVc(N,A,Thresh,B)
      Implicit Real*8(A-H,O-Z)
C
C     Copy vector A to B, set all elements less than Thresh to 0
C     and return the number of zeroes.
C      
      Dimension A(1), B(1)
      Save Zero
      Data Zero/0.0d0/
C
      ICount = 0 
      Do 10 I = 1, N
        B(I) = A(I)
        If(Abs(B(I)).le.Thresh) then
          ICount = ICount + 1
          B(I) = Zero
          endIf
   10   Continue
      IPopVC = ICount
      Return
      End
*Deck Place
      Subroutine Place(CutOff,ZH,CRef,Bl,Alpha,C,NonLin)
      Implicit Real*8(A-H,O-Z)
C
C     Place an atom when all previous atoms are collinear.
C     ZH is a unit vector along the collinear axis, CRef
C     is the atom referenced to the present atom, Bl and
C     Alpha are the bond lengths, and C is the output
C     vector.  The atom is placed in the plane defined
C     by ZH and the X-axis if possible; otherwise the Z-axis
C     is used.
C
      Logical NonLin
      Dimension ZH(3), CRef(3), C(3), CLoc(3), XH(3)
      Save Zero, One
      Data Zero/0.0d0/, One/1.0d0/
C
      CLoc(1) = Bl*Sin(Alpha)
      CLoc(2) = Zero
      CLoc(3) = -Bl*Cos(Alpha)
      NonLin = Abs(CLoc(1)).gt.CutOff
      XH(1) = One - ZH(1)*ZH(1)
      XH(2) = -ZH(1)*ZH(2)
      XH(3) = -ZH(1)*ZH(3)
      R = SProd(3,XH,XH)
      If(R.le.CutOff) then
        XH(1) = -ZH(3)*ZH(1)
        XH(2) = -ZH(3)*ZH(2)
        XH(3) = One - ZH(3)*ZH(3)
        R = SProd(3,XH,XH)
        endIf
      Do 10 I = 1, 3
        C(I) = CLoc(1)*XH(I) + CLoc(3)*ZH(I) + CRef(I)
   10 Continue
      Return
      End
