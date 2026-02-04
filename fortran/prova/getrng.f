      Implicit None
      Integer MxAt,MxRng,MxRnAt,MxBnd,MxBrn
      Parameter(MxAt=100,MxRng=10,MxRnAt=4,MxBnd=4,MxBrn=10)
      Integer IAt,JAt,NRing,NBrnch,NSter2,NTest,NrAct,I,J,IR,IEZ
      Integer In,IOut,IPunch,NSMI,NAtoms,IB,IElIAn,I1,IL,ISign,MssTmp,NBrCk,
      Integer II,NMol,IMol,NHeavy,IndAr,IXYZ,Num,NZ
      Integer NNAtR(MxAt),IRing(MxRnAt,MxAt),IBefBr(MxBrn)
      Integer IRAct(MxRng),IFirst(MxRng),ILast(MxRng)
      Integer IAn(MxAt),IArom(MxAt),IBrnch(MxAt),NH(MxAt),NPI(MxAt)
      Integer NBond(MxAt),IBond(MxBnd,MxAt)
      Integer Isotop(MxAt),ISter1(MxAt),ISter2(3,MxAt)
      Integer NAtRng(MxRng),IAtRng(MxRnAt,MxRng)
      Real*8 BndTyp
      Real*8 Charge(MxAt),BndOrd(MxBnd,MxAt)
      Logical Open,FndSp,FndAt,ActBrc,First,ActivB(MxBrn)
      Logical Error,TstAng,Ttest 
      Character*100 TITLE
      Character*1 Test,Test1,SMILES(100)
      Character*2 El,IAnEl2,El2
      Character*80 LinScr
      Integer IZ(4,MxAt)
      Integer LBl(MxAt),LAlpha(MxAt),LBeta(MxAt)
      Integer IAnZ(MxAt),IAtomB(2,MxAt)
      Real*8 BL(MxAt),Alpha(MxAt),Beta(MxAt)
      Real*8 Pi,ToAng,ToDeg,ToRad,ScalI
      Real*8 C(3,MxAt),CZ(3,MxAt)
      Real*8 A(MxAt),B(MxAt),D(MxAt),Alpha1(MxAt),Beta1(MxAt)
C Define I/O files
      In=5
      Iout=6
      IPunch=7
      OPEN(In,FILE='ringin',STATUS='UNKNOWN')
      OPEN(IOut,FILE='ringout',STATUS='UNKNOWN')
      OPEN(IPunch,FILE='ringXYZ',STATUS='UNKNOWN')
      Read(IN,*) NMol
      Do IMol=1,NMol
C Interpret SMILES
       CALL SMIF77(In,IOut,MxAt,MxBnd,MxRng,MxBrn,MxRnAt,NAtoms,NRing,
     $   NBrnch,IAn,Isotop,IArom,NPI,NH,IBrnch,IRing,ISter1,ISter2,
     $   NBond,IBond,NAtRng,IAtRng,BndOrd,Charge)  
C Set IZ
       CALL IClear(4*NAtoms,IZ)
       write(IOut,'(''   Atom  Branch  NH  NBond          IBond'')')
       Do IAt=1,NAtoms
        Write(IOut,'(8I6)') IAn(IAt),IBrnch(IAt),NH(IAt),NBond(IAt),
     $   (IBond(ii,IAt),ii=1,NBond(IAt))
       EndDo   
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
     $     IZ(2,IAt)
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
       CALL SetDih(IOut,MxBnd,MXRnAt,NAtoms,IZ,NBond,IBond,IBrnch,IRing,
     $  NAtRng,IAtRng,Alpha,Beta)
       CALL IClear(NAtoms,LBl)
       CALL IClear(NAtoms,LAlpha)
       CALL IClear(NAtoms,LBeta)
C Temporary
       ToAng=1.0d0
       ToDeg=1.0d0
C Print Z-matrix
       CALL ZMTPRT(IOut,NAtoms,IAN,IZ,LBL,LAlpha,LBeta,BL,ALPHA,BETA,
     $   ToAng,ToDeg)
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
C Next Molecule
      EndDo
      End
*Deck SMIF77
      Subroutine SMIF77(In,IOut,MxAt,MxBnd,MxRng,MxBrn,MxRnAt,NAtoms,
     $  NRing,NBrnch,IAn,Isotop,IArom,NPI,NH,IBrnch,IRing,ISter1,ISter2,
     $  NBond,IBond,NAtRng,IAtRng,BndOrd,Charge)
      Implicit None
      Integer MxAt,MxRng,MxRnAt,MxBnd,MxBrn
C     Parameter(MxAt=100,MxRng=10,MxRnAt=4,MxBnd=4,MxBrn=10)
      Integer IAt,JAt,NRing,NBrnch,NSter2,NTest,NrAct,I,J,IR,IEZ
      Integer In,IOut,NSMI,NAtoms,IB,IElIAn,I1,IL,ISign,MssTmp,NBrCk,II
      Integer NMol,IMol,NHeavy,IndAr
      Integer NNAtR(MxAt),IRing(MxRnAt,MxAt),IBefBr(MxBrn)
      Integer IRAct(MxRng),IFirst(MxRng),ILast(MxRng)
      Integer IAn(MxAt),IArom(MxAt),IBrnch(MxAt),NH(MxAt),NPI(MxAt)
      Integer NBond(MxAt),IBond(MxBnd,MxAt)
      Integer Isotop(MxAt),ISter1(MxAt),ISter2(3,MxAt)
      Integer NAtRng(MxRng),IAtRng(MxRnAt,MxRng)
      Real*8 BndTyp
      Real*8 Charge(MxAt),BndOrd(MxBnd,MxAt)
      Logical Open,FndSp,FndAt,ActBrc,First,ActivB(MxBrn)
      Character*100 TITLE
      Character*1 Test,Test1,SMILES(100)
      Character*2 El,IAnEl2,El2
      Do i=1,100
       SMILES(i)=' '
       TITLE(i:i)=' '
      End Do
      Read(In,*) Title
      Write(IOut,'('' '')')
      Write(IOut,'(A80)') Title(1:80)
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
      Do I=1,MxRng 
       IFirst(I)=0
       ILAst(I)=0
       IRAct(I)=0
       NNAtr(I)=0
       NAtRng(I)=0
       IBefBr(I)=0
       Do J=1,MxRnAt
        IRing(J,I)=0 
        IAtRng(J,I)=0
       EndDo
      EndDo
      Do I=1,MxBrn
       ActivB(I)=.false.
      EndDo
      Do IAt=1,MxAt
       NBond(IAt)=0
       ISotop(IAt)=0
       IBrnch(IAt)=0
       NH(IAt)=0
       ISter1(IAt)=0
       IArom(IAt)=0
       NPi(IAt)=0
       Do IB=1,MxBnd
        IBond(IB,IAt)=0
       EndDo
       Do IB=1,3
        ISter2(IB,IAt)=0
       EndDo
      EndDo
      NRAct=0
      NRing=0
      NBrCk=0
      NBrnch=0
      NSter2=0
      IndAr=0
      IAt=0
      JAt=0
      I=1
      First=.true.
      ActBrc=.False.
 100  CONTINUE
      IF (I.GT.NSMI) GOTO 200
      FndAt=.false.
      Test = SMILES(I)
C Handle Bonds
        IF (Test.eq.'-') then
          BndTyp =1.0D0
          I=I+1
          GOTO 100
        ELSEIF (Test .EQ. '=') THEN
          BndTyp = 2.0d0
          I=I+1
          GOTO 100
        ELSEIF(Test.eq.'#') then
          BndTyp = 3.0d0
          I = I + 1
          GOTO 100
        ENDIF
C Handle Square Brackets [ ]
        IF (Test.eq.'[') Then
         ActBrc=.True.
         I=I+1
         GOTO 100
        ElseIf(Test.eq.']') Then
         ActBrc=.False.
         I=I+1
         GOTO 100
        EndIf 
C Handle Info within '['
        If(ActBrc) then
C Handle Stereochemistry
         IF (Test.EQ.'@') THEN
          I = I + 1
          Test1 = SMILES(I)
          IF (Test1.EQ.'@') THEN
C S Stereoisomer (@@)
           ISter1(IAt) = -1
           I = I + 1
          ELSE
C R Stereoisomer (@)
           ISter1(IAt) = 1
          ENDIF
          GOTO 100
         ENDIF
C Handle E/Z Isomers
         If(Test.eq.'/'.or.Test.eq.'\') THEN
          CALL SetEZ(Test,IAt,NSter2,ISter2)
          I = I + 1
          GOTO 100
         ENDIF
C Handle Charges ( + o - )
         IF (Test .EQ. '+' .OR. Test .EQ. '-') THEN
          ISign = 1
          IF (Test .EQ. '-') ISign = -1
          I = I + 1
          Test1 = SMILES(I)
          IF (Test1 .GE. '0' .AND. Test1 .LE. '9') THEN
C Account for Sign
           Charge(Natoms) = ISign * (ICHAR(Test1) - ICHAR('0'))
           I = I + 1
          ELSE
           Charge(NAtoms) = ISign
          ENDIF
          I = I + 1
          GOTO 100
         ENDIF
C Check Masses (1-2 Digits)
         MssTmp=0
         If(Test.ge.'0'.and.test.le.'9')then
          MssTmp = ICHAR(Test) - ICHAR('0')
          J=I+1
          Test1 = SMILES(J)
          IF (Test1 .GE. '0' .AND. Test1 .LE. '9') THEN
           MssTmp = MssTmp * 10 + (ICHAR(Test1) - ICHAR('0'))
          ENDIF 
          ISotop(IAt+1)=MssTmp
          I=J+1
          GOTO 100 
         ENDIF
C End of Info within '[]'
        EndIf
C Handle Branches '()'
        IF (Test .EQ. '(') THEN
         Call GetBrn(IOut,MxBrn,IAt,JAt,NBrnch,'open',ActivB,IBefBr)
C Handle the sequence ']('
C        If(SMILES(I-1).eq.']') JAt=IBefBr(NBrnch)
         I = I + 1
         GOTO 100
        ELSE IF (Test .EQ. ')') THEN
         Call GetBrn(IOut,MxBrn,IAt,JAt,NBrnch,'clos',ActivB,IBefBr)
         I=I+1
         GOTO 100
        ENDIF
C Handle Atoms
      If (Test.eq.'H') then
       J=I+1
       Test1=SMILES(J)
       If(Test1.gt.'1'.and.Test1.le.'9'.and.Test.eq.'H') THEN
        NH(IAt)=NH(IAt) + IChar(Test1) - IChar('0')
        I=I+2
        GOTO 100
       Else
        NH(IAt)=NH(IAt)+1
        I=I+1
        GOTO 100
       EndIf
      ElSEIF (Test .GE. 'A' .AND. Test .LE. 'Z') THEN
       FndAt=.true.
       El(1:1)=Test
       El(2:2)=' '
       IAt=IAt+1
       NTest=0
       J=I+1
       Test1 = SMILES(J)
       IF (Test1 .GE. 'a' .AND. Test .LE. 'z') THEN
        El(2:2)=test1
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
        J=I+1
        Test1=SMILES(J)
       EndIf
C Aromatic Atoms (lowercase letters)
      ElseIf(Test .eq. 'c' .or. Test .eq. 'n' .or. Test .eq. 'o' .or.
     $   Test .eq. 'p'. or. Test .eq. 's') THEN
C Convert to capital letters (e.g. 'c' --> 'C')
       FndAt=.true.
       IndAr=1
       IAt=IAt+1
       NTest=0
       El(1:1) = CHAR(ICHAR(Test) - 32) 
       El(2:2) = ' '
       J=I+1
       Test1=SMILES(J)
      EndIf
      If(.not.FndAt) then
       I=I+1
       goto 100
      EndIf
      If(Test1. GT. '0' .AND. TEST1.LE. '9') THEN
       I=I+1
       NTest  = ICHAR(Test1) - ICHAR('0')
      EndIf
  300 Write(IOut,'('' Entering GetRng with IAtom='',I3,'' NTest='',I2,
     $  '' Nring='',I2,'' NRAct='',I2,'' El ='',A2)') IAt,NTest,NRing,
     $  NRAct,El 
      Call GetRng(IOut,MxRNAt,IAt,NTest,NRing,NRAct,NNAtR,IRing,IRAct,
     $ IFirst,ILast)
      If(First) then
       IArom(IAt)=IndAr
       IBrnch(IAt)= NBrnch
       IAn(IAt)=IElIAn(El)
       IF(IAt.ne.1.and.IAn(IAt).gt.1) then
        If(JAt.eq.0) JAt=IAt-1
        NBond(IAt) = NBond(IAt) + 1
        NBond(JAt) = NBond(JAt) + 1
        IBond(NBond(IAt), IAt) = JAt
        IBond(NBond(JAt), JAt) = IAt
        BndOrd(NBond(IAt), IAt) = BndTyp
        BndOrd(NBond(JAt), JAt) = BndTyp
       EndIf 
      EndIf
C Manage the opening and closing of several rings on the same atom
      J=I+1
      Test1=SMILES(j)
      If(Test1. GT. '0' .AND. TEST1.LE. '9') THEN
       First=.false.
       I=I+1
       NTest  = ICHAR(Test1) - ICHAR('0')       
       goto 300
      Else 
       JAt=0
       BndTyp=1.0d0
       IndAr=0
       First=.true.
       I = I + 1
       GOTO 100
      EndIf
  200 CONTINUE
      NAtoms=IAt
      Write(IOut,'(I3,'' Atoms and'',I3,'' Rings Found'')') NAtoms,
     $  NRing
C Add the bond between first and last atom of rings
      If(NRing.gt.0) then
       Do IR=1,Nring
        I1=IFirst(IR)
        IL=ILast(IR)
        NBond(I1)=NBond(I1)+1
        NBond(IL)=NBond(IL)+1
        IBond(NBond(I1),I1)=IL
        IBond(NBond(IL),IL)=I1
        BndOrd(NBond(I1),I1)=1.0D0
        BndOrd(NBond(IL),IL)=1.0D0
       EndDo
      EndIf
C Update NatRng and IAtRng and print results
C Check that monovalent atoms are not in rings
      Do IAt=1,NAtoms
       El2=IAnEl2(IAn(IAt))
       If(NBond(IAt).eq.1.and.NRing.gt.0) then
        Do ii=1,NNAtR(IAt) 
         IRing(ii,IAt)=0
        EndDo 
        NNAtr(IAt)=0
       EndIf
       If(NNAtR(IAt).eq.0) cycle
       Do IR=1,NNAtr(IAt)
        NAtRng(IRing(IR,IAt))=NAtRng(IRing(IR,IAt))+1
        IAtRng(NAtRng(IRing(IR,IAt)),IRing(IR,IAt))=IAt
       EndDo
      EndDo
      IF(NRing.gt.0) then
       Do IR=1,Nring
        write(IOut,'(''  Ring'',I2,'' with'',I2,'' Atoms:'',10I3)')
     $    IR,NAtRng(IR),(IAtRng(I,IR),I=1,NatRng(IR))
       EndDo
      EndIf
C Assign Pi bonds and bond orders
      CALL PiBond(IOut,MxBnd,NAtoms,IArom,NPi,NBond,IBond,BndOrd)
C Add Hydrogen Atoms both explicit and implicit
      NHeavy=NAtoms
      CALL AddHyd(IOut,MxBnd,NHeavy,NAtoms,IAn,NPi,NBond,IBond,
     $  BndOrd,NH)
C Print Results
      CALL SmiRes(IOut,MxBnd,MxRnAt,NAtoms,IAn,NPi,IArom,IBrnch,IRing,
     $  NBond,IBond,BndOrd,NH,ISter1,NSter2,ISter2,Isotop)
      Return
      End
*Deck GetRng
      Subroutine GetRng(IOut,MxRnAt,IAt,NTest,Nring,NRAct,NNAtR,IRing,
     $  IRAct,IFirst,ILAst) 
      Implicit None
C IAt           = Atomo che precede il numero del ciclo
C MxRnAt        = Massimo numero di cicli cui un atomo può partecipare 
C NRing         = Numero Totale di cicli 
C NTest         = Numero d'ordine del ciclo che cambia (0 per non cambiare)
C NNAtR(IAt)    = Numero di cicli cui IAt già partecipa
C IRing(NNAtr(IAt),IAt) = Numero del ciclo NNAtR(IAt) dell'atomo IAt
C NRAct         = Numero di Cicli Attivi
C IRAct(NRAct)  = Numero d'ordine dei cicli attivi      
C IFirst(NRing) = Primo Atomo di un ciclo
C ILAst(NRing)  = Ultimo Atomo di un ciclo
      Integer IOut,MxRnAt,IAt,NTest,NRing,NRAct,IPos,IR,II,JJ,IEl
      Integer IFirst(*),ILast(*),IRAct(*),IRing(MxRnAt,*),NNAtR(*)
      Logical Open,OK
      If(NTest.lt.0.or.NTest.gt.9) then
       write(IOut,'(I3,'' Rings, but only 1 to 9 are allowed'')') NTest 
       Stop
      ElseIf(NTest.eq.0) then
C Stay on the Same Ring
       IF(NRAct.gt.0) then
        NNATR(IAt)=NNAtR(IAt)+1
        IRing(NNAtR(IAt),IAt)=IRAct(NRAct)
       EndIf
       Return
      EndIf
C Check ring
      IPos=0 
      Do IR=1,NRAct
       If(IRAct(IR).eq.NTest) IPOs=IR
      EndDo
      IF(IPos.gt.0) then
C Close an open ring
       NNAtR(IAt)=NNAtR(IAt)+1
       IRing(NNAtR(IAt),IAt)=NTest 
       ILAst(NTest)=IAt
       IRAct(IPos)=0
       If(IPos.gt.1.and.IRAct(IPos-1).ne.0) then
        NNAtR(IAt)=NNAtR(IAt)+1
        IRing(NNAtr(IAt),IAt)=IRAct(IPos-1)
       EndIf
       IF(IPos.ne.NRAct) then
        NNAtR(IAt)=NNAtR(IAt)+1
        IRing(NNAtr(IAt),IAt)=IRAct(NRAct)
       EndIf
       Do IR=IPos+1,NRAct
        IF(IR.gt.1) IRact(IR-1)=IRAct(IR)
       EndDo
C Check that the preceding atom belongs to the ring 
       OK=.False.
       Do II=1,NNAtR(IAt-1)
        If(IRing(ii,IAt-1).eq.NTest) OK=.true.
       EndDo
       If(.not.OK) then
        NNAtR(IAt-1)=NNAtR(IAt-1)+1
        IRing(NNAtR(IAt-1),IAt-1) = NTest
       EndIf
       NRAct=NRAct-1
C Avoid repeating the same ring several times
       If(NNAtR(IAt).ge.2) then
        IEl=0
        do ii=NNAtR(IAt),2,-1
         do jj=1,ii-1
          IF(IRing(ii,IAt).eq.IRing(jj,IAt)) then
           IRing(ii,IAt)=0
           IEl=IEl+1
          EndIf
         EndDo
        EndDo
         NNAtR(IAt)=NNAtR(IAt)-IEl 
       EndIf
      Else
C Open a new ring
       NRing=NRing+1
       IF(NRAct.gt.0) then
        NNAtR(IAt)=NNAtR(IAt)+1
        IRing(NNAtR(IAt),IAt)=IRAct(NRAct)
       EndIf
       NRAct=NRAct+1
       IRAct(NRAct)=NTest 
       IFirst(NTest)=IAt
       NNAtR(IAt)=NNAtR(IAt)+1
       IRing(NNAtR(IAt),IAt)=NTest
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
*Deck IElIAN
      Integer Function IElIAN(AtSymb)
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
       IElIAN = -1
      else
       IElIAN = i/2
      endIf
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
      Subroutine GetBrn(IOut,MxBrn,IAt,JAt,IBrnch,Action,ActivB,IBefBr)
      Implicit None
      Integer IOut,MxBrn,I,IAt,JAt,IBrnch,IBefBr(*)
      Character Action*4
      Logical ActivB(*)
      If(IBrnch.lt.0.or.IBrnch.gt.9) then
       write(IOut,'(I3,'' Branches, but only 1 to 9 are allowed'')')
     $   IBrnch
       Stop
      ElseIf(IAt.lt.0.or.IAt.gt.99) then
       write(IOut,'(I3,'' Atoms, but only 1 to 99 are allowed'')') IAt
      EndIf
      do i=1,MxBrn
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
      Subroutine SmiRes(IOut,MxBnd,MxRnAt,NAtoms,IAn,NPi,IArom,IBrnch,
     $  IRing,NBond,IBond,BndOrd,NH,ISter1,NSter2,ISter2,Isotop)
      Implicit Real*8 (A-H,O-Z)
      Dimension NPi(*),IArom(*),NH(*),NBond(*),IBond(MxBnd,*)
      Dimension IBrnch(*),IRing(MxRnAt,*),ISter1(*),ISter2(3,*)
      Integer IAn(*),Isotop(*)
      Dimension BndOrd(MxBnd,*)
      Character*1 Star,Chir,OpenP,CloseP
      Character*1 CT(MxBnd) 
      Character*2 El,IAnEl2
      OpenP='('
      CloseP=')' 
      write(IOut,'(/,'' Atom Branch Ring  Pi El. nH Stereo '',
     $  6X,''Bonded Atoms'','' and Bond Orders'')')
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
       Else
        write(IOut,'(1X,A2,I2,2I5,1X,A1,3X,I2,3X,I2,4X,A1,4X,
     $   4(I3,A1,F8.4,A1))') El,IAt,IBrnch(IAt),
     $   IRing(1,IAt),star,NPi(IAt),NH(IAt),Chir,
     $   (IBond(J,IAt),OpenP,BndOrd(J,IAt),CloseP,J=1,NBI) 
       End If 
      End Do
C Description of Isotopes and Isomers 
      write(IOut,'('' * labels Aromatic Atoms'')')
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
C     DO IAt=1,NHeavy
C      If(NH(IAt).gt.0) then
C       write(IOut,'('' Atom'',I3,'' has'',I2,'' explicit hydrogens'')')
C    $   IAt,NH(IAt)
C        NHI=NH(IAt)
C       Do IH=1,NHI
C        IHAt=IHAt+1
C        IAn(IHAt)=1
C        NBond(IAt)=NBond(IAt)+1
C        NBI=NBond(IAt)
C        IBond(NBI,IAt)=IHAt
C        BndOrd(NBI,IAt)=1.0d0
C        NBond(IHAt)=1
C        IBond(1,IHAt)=IAt
C        BndOrd(1,IHAt)=1.0D0
C       End Do
C      End If
C     End Do 
C Add Hydrogen Atoms
      Do 10 IAt=1,NHeavy
       If(IAn(IAt).eq.1) goto 10
       MissH=NFree(IOut,IAt,IAn(IAt),NBond(IAt),NPi(IAt))  
       If(MissH.lt.NH(IAt)) MissH=NH(IAt)
       If(NBond(IAt)+MissH.gt.MxBnd) then
        write(IOut,'('' Too many Bonds for atom'',I3)') IAt
        stop
       EndIf 
       write(IOut,'('' Atom'',I3,'' has'',I2,'' explicit and'',I2,
     $    '' implicit hydrogens'')')IAt,NH(IAt),MissH-NH(IAt)
       NH(IAt)=MissH
       IF(NH(IAt).gt.0) then 
        Do 20 JAt = 1, NH(IAt) 
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
C      Stop
      EndIf
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
        If(IUsed(IZ(1,IAt)).eq.3) IZ(4,IAt)=0
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
        If(IUsed(JAt).eq.3) IZ(4,IAt)=0
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
          If(IUsed(JAt).eq.3) IZ(4,IAt)=0 
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
         If(IUsed(JAt).eq.3) IZ(4,Iat)=0
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
