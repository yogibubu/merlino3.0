*Deck DiNa25
C
C New version started 2th November 2024
C This deck contains the driver and the routines
C Aclear,Amove,IClear,LinUpC,PhyFil
C and the function CnvFct
C
      Implicit Real*8 (A-H,O-Z)
C MaxEl  = Number of different elements
C MxScr  = Dimension of Scratch Array
C MxAt   = Maximum Number of Atoms
C MaxNZ  = Maximum Number of ZMatrix Elements
C MxBnd  = Maximum coordination 
C MxBox  = Maximum number of boxes for new connectivity
C MxKwd  = Maximum Number of KeyWords
C LenPhy = Number of Physical Constants
C MxGNIC = Maximum Number of GNICS
C MxTrm  = Maximum Number of Terms in GNICs
C MxAtP  = Maximum Number of Atoms in a Primitive
C MxCyc  = Maximum Number of Cycles
C MxAtCy = Maximum Number of Atoms in a Cycle
C MxPot  = Max Terms in Potential Fitting
      Parameter(MaxEl=200)
      Parameter(MxScr=100000,MxAt=1000,MxBnd=10,MxKwd=50,LenPhy=30)
      Parameter(MxGNIC=1000,MxTrm=15,MxAtP=4,MxCyc=10,MxAtCy=10)
      Parameter(MxFrg=100,MxAtFr=30,MxAtB=100,MaxNZ=1000,MxBox=1000)
      Parameter(MxPot=20)  
      Character*80 FilNam,InFil,OutFil,GauKwd,Title,React,LinScr
      Character*20 StrInp
      Character*100 SMILES
CENZO
      Character*16 Group
      Common/IO/In,IOut,IPunch
      Common/PhyCon/PhyCon(30)
      Common/bic/N2Cyc,N3Cyc,IAt2C(3,20),Iat3C(4,20)
      Common/bic1/NBrL,NBrA,NBrD,IBrL(4,20),IBrA(5,20),IBrD(6,20)
      Integer St2Int,El2IAn,Istart(100)
      Logical Kwd(MxKwd),Loose
      Logical Linear,ImpDih,DoVolt,DoGDV,PrtPic,SyGNIC,DoBPCS,Clean
      Logical DoEck,Do1Dih,DoNorm,DVibRot,InvDst,LConn,DoSySt,Aver
      Logical DoBMat,DoMW,Inv1,DoScan,DoRig,RIgB,RigA,RigL,RigD,RigO
      Logical DoneC,DoB1,DoVCI,DoDVR,DoFit,RdData,DoGNIC,DoZMAt,RdIsot
      Logical WrZMat,RdSMI,DoG16,DoColl,DoGor,DoSpin
      Logical RdB0,RdB0Er,RdVib,RdDvEr,RdEle,RdDeEr
      Logical RdXYZ,RdFChk,TstAng,TTest,Error,DoVMSR,PrtVal
      Dimension IEl(0:MaxEl)
      Dimension IScr(MxScr) 
      Dimension Scr(MxScr)
C Dimensions for Atom Properties (including Coordinates)
      Dimension IAn(MxAt),IFrag(MxAt)
      Dimension ISot(MxAt),MultN(MxAt)
      Dimension IAnZ(MaxNZ),IZ(4,MaxNZ),MapZAt(MaxNZ)
      Dimension LBl(MaxNZ),LAlpha(MaxNZ),LBeta(MaxNZ)
      Dimension BL(MaxNZ),Alpha(MaxNZ),Beta(MaxNZ)
      Dimension QMom(MxAt),GFac(MxAt)
      Dimension PMom(3),RotGHz(3),RTemp(3),DBVib(3),AtMass(MxAt)
      Dimension XYZCM(3),RotMat(3),Rotcm1(3),PMomB(3)
      Dimension C(3,MxAt),CZ(3,MaxNZ),TMom(6)
C Scratch for ZTOC
      Dimension A(MaxNZ),B(MaxNZ),D(MaxNZ),Alpha1(MaxNZ),Beta1(MaxNZ)
C Dimensions for vibrations
      Dimension Frq(3*MxAt)
C Dimension for Fragments
      Dimension IFrsAt(MxFrg),ILstAt(MxFrg),NAtFr(MxFrg)
      Dimension IAtFr(MxAtFr,MxFrg),LConn(MxFrg)
C Dimensions for Connectivity
      Dimension NBond(MxAt),IBond(MxBnd,MxAt),NH(MxAt),IArom(MxAt)
      Dimension BndOrd(MxBnd,MxAt),InCyc(MxAt)
C Dimensions for Stretchings
      Dimension IAtomB(MxAtP,MxTrm,MxGNIC),NTermB(MxGNIC),ITVB(MxGNIC)
      Dimension NEqAt(2,MxGNIC),IFixB(MxGNIC)
      Dimension CoefB(MxTrm,MxGNIC),ValTB(MxGNIC)
C Dimensions for Bendings
      Dimension IAtomA(MxAtP,MxTrm,MxGNIC),NTermA(MxGNIC),ITVA(MxGNIC)
      Dimension IFixA(MxGNIC)
      Dimension CoefA(MxTrm,MxGNIC),ValTA(MxGNIC)
C Dimensions for Linear Bendings
      Dimension IAtomL(MxAtP,MxTrm,MxGNIC),NTermL(MxGNIC),ITVLA(MxGNIC)
      Dimension IFixL(MxGnic)
      Dimension CoefL(MxTrm,MxGNIC),ValTL(MxGNIC)
C Dimensions for Torsions
      Dimension IAtomD(MxAtP,MxTrm,MxGNIC),NTermD(MxGNIC),ITVD(MxGNIC)
      Dimension IPerD(MxGNIC),IFixD(MxGNIC)
      Dimension CoefD(MxTrm,MxGNIC),ValTD(MxGNIC)
C Dimensions for Out-of-Plane Bendings
      Dimension IAtomO(MxAtP,MxTrm,MxGNIC),NTermO(MxGNIC),ITVO(MxGNIC)
      Dimension IFixO(MxGNIC)
      Dimension CoefO(MxTrm,MxGNIC),ValTO(MxGNIC)
C Dimensions for Syntons
      Dimension EAN(MxAt),EANZ(MaxNZ)
C Dimensions for Cycles
      Dimension NatC(MxCyc),ICAt(MxAtCy,MxCyc)
      Dimension IAtCyc(MxAt),IBr(2,MxAt)
C Dimensions for primitives
      Dimension IAtmBR(MxAtP,MxTrm,MxGNIC),IPrimB(MxTrm,MxGNIC)
      Dimension IAtmAR(MxAtP,MxTrm,MxGNIC),IPrimA(MxTrm,MxGNIC)
      Dimension IAtmLR(MxAtP,MxTrm,MxGNIC),IPrimL(MxTrm,MxGNIC)
      Dimension IAtmDR(MxAtP,MxTrm,MxGNIC),IPrimD(MxTrm,MxGNIC)
      Dimension IAtmOR(MxAtP,MxTrm,MxGNIC),IPrimO(MxTrm,MxGNIC)
C Dimension for B and DB matrices
      Dimension BMat(3*MxAtB,MxGNIC)
C Dimensions for BDPCS3
      Dimension R0IJ(MxGNIC),RBIJ(MxGNIC)
C ACN and ACN1 are coordination numbers (from D3 and D4 of Grimme)
C SPIJ is the sum of bond orders
      Dimension ACN(MxAt),ACN1(MxAt),SPIJ(MxAt),Teta0(4)
C Dimensions for Potential Fitting
      Dimension IPot(MxPot),Cof(MxPot)
C Dimension for Gaussian Keywords
      Dimension IDeriv(4)
C Long Formats
 1000 Format(' Atomic Number =',I3,'; Isotope =',I3,
     $  '; Nuc.Spin Multip =',I5,/,' Mass =',F12.5,
     $  '; Quadrup.Mom. =',F10.5,'; Nuc.Magn.Mom. =',F10.5)
C Input/Output Files
      In=1
      IOut=2
      IPunch=7
C Input and Output Files are wayin and wayout
      filnam='topo24'
      ITpFil=0
      OPEN(In,FILE='provin',STATUS='OLD')
      OPEN(IOut,FILE='provout',STATUS='UNKNOWN')
      Rewind(IOut)
C which version of physical constants
      IOpt=2010
C Set physical constants and treshold for linear angles in radiants 
C (171.9 degrees)
      Call PhyFil(IOpt,LenPhy,PhyCon)
      pi=4.0d0*ATan(1.0d0)
      ToDeg=1.80d+2/pi
      TreshL=3.0d0 
C Set Atomic Symbols
      call FillEl(0,MaxEl,IEl)
C Read keywords
      Call LlinCl(LinScr)
      Read(In,'(A80)') LinScr
      Write(IOut,'(A80)') LinScr
      Call FndKwd(IOut,IPrint,ModPCS,IDeriv,LinScr,Kwd,MxKwd)
C     RdSMI=Kwd(1)
      DoG16=Kwd(2)
      DoGDV=Kwd(3)
      Aver=Kwd(5)
      DoEck=Kwd(6)
      ImpDih=Kwd(7)
      DoGNIC=Kwd(8)
      InvDst=Kwd(9)
      DoSySt=Kwd(10)
      Do1Dih=Kwd(11)
      DoNorm=Kwd(12)
      DoBPCS=Kwd(13)
      SyGNIC=Kwd(14)
      DoVolt=Kwd(18)
      DoBMat=Kwd(19)
      DoScan=Kwd(20)
      DoRig=Kwd(21)
      DoZMat=Kwd(22)
      WrZMat=Kwd(23)
      RdIsot=Kwd(24)
      DoFit=Kwd(25)
      DoVCI=Kwd(26)
      DoDVR=Kwd(27)
      DoColl=Kwd(28)
      DoGor=Kwd(29)
      DoSpin=Kwd(30)
      Clean=Kwd(31)
      DVIBRot=Kwd(32)
      Loose=Kwd(33)
      RdXYZ=Kwd(34)
      RdFChk=Kwd(35)
      RdSMI=Kwd(36)
      DoVMSR=Kwd(37)
      If(RdSMI) Loose=.True.
      If(DoG16) ImpDih=.true.
      RigB=.False.
      RigA=.False.
      RigL=.False.
      RigD=.False.
      RigO=.False.
      If(DoBPCS) RigB=.True.
      If(DoScan.and.DoRig) then
       RigB=.True.
       RigA=.True.
       RigL=.True.
       RigO=.True.
      EndIf
C Read Title preceeded and followed by a blank card (Gaussian style)
      Read(In,'(A80)')    LinScr
      Read(In,'(A80)')    Title
      Write(IOut,'(/,A80)') Title
      Read(In,'(A80)')    LinScr
C Fitting and Vibrational programs 
      If(DoFit.or.DoVCI.or.DoDVR) then
       call NewVib(In,Iout,IPrint,DoFit,DoVCI,DoDVR,MxScr,Scr,IScr)
       goto 999
      EndIf
C Read charge and multiplicity
      call LlinCl(StrInp)
      Read(In,'(A20)') Strinp
      call SubStr(StrInp,2,IStart,NValue)
      If(NValue.ne.2) then
       write(IOut,'(A20)') Strinp
       write(IOut,'(''Wrong Charge or Multiplicity'')')
       Stop
      EndIf
      ICharg=St2Int(StrInp(IStart(1):IStart(2)-1),-999)
      If(ICharg.eq.-999) then
       write(IOut,'(''Wrong Charge'')')
       Stop
      EndIf
      Multip=St2Int(StrInp(IStart(2):len(StrInp)),-999)
      If(Multip.eq.-999) then
       write(IOut,'(''Wrong Multiplicity'')')
       Stop
      EndIf
C Read geometry and set connectivity
      If(RdSMI) then
       Read(In,*) SMILES
       write(IOut,*) SMILES
C Interpret SMILES
       call SMIF77(IOut,IPrint,MxAt,MxBnd,NAtoms,IAn,IArom,
     $  InCyc,NBond,IBond,IZ,BndOrd,BL,Alpha,Beta,SMILES)
C Print Atom Informations 
C      call SmiBnd(IOut,MxBnd,NAtoms,IAn,IArom,NBond,IBond,
C    $   BndOrd)
       Call IClear(NAtoms,LBl)
       Call IClear(NAtoms,LAlpha)
       Call IClear(NAtoms,LBeta)
C Temporary
       ToAng=1.0d0
       ToDeg=1.0d0
C Print Z-matrix      
       Call ZPrint(IOut,NAtoms,IAN,IZ,LBL,LAlpha,LBeta,BL,ALPHA,BETA,
     $  ToAng,ToDeg)
C transform to Cartesian coordinates
C checking tetrahetral angles (TTest) and 0>teta<180 (TstAng)
       TstAng=.true.
       TTest=.true.
       NZ=NAtoms
       call AMove(NAtoms,IAn,IAnZ)
       call AClear(3*NAtoms,C)
       call AClear(3*NAtoms,CZ)
       call AClear(NAtoms,A)
       call AClear(NAtoms,B)
       call AClear(NAtoms,D)
       call Aclear(NAtoms,Alpha1)
       call Aclear(NAtoms,Beta1)
C Angles should be in radiants
       ToRad=pi/180.0d0
       Error=.False.
       Do IAt=1,NAtoms
        Alpha(IAt)=Alpha(IAt)*ToRad
        Beta(IAt)=Beta(IAt)*ToRad
       EndDo
       call ZtoC(MaxNZ,NZ,IAnZ,IZ,Bl,Alpha,Beta,TTest,NAtoms,IAn,C,
     $   CZ,A,B,D,Alpha1,Beta1,IOut,Error,TstAng)
       If(Error) Stop 
       Write(IOut,'(/,10X,''Cartesian Coords. from SMILES'')')
       Do IAt=1,NAtoms
        Write(IOut,'(I5,2X,A2,3F12.5)') IAt,IEl(IAn(IAt)),
     $  (C(IXYZ,IAt),IXYZ=1,3)
       End Do
C Print interatomic distances
       Call LlinCl(LinScr) 
       LinScr(1:22)=' Interatomic Distances (Angstrom)'
       Call HedPrt(IOut,0,LinScr,Num)
       ScalI=1.0d0
       Call DisMat(NAtoms,IAN,C,2,5,IOut,Error,0,ScalI)
       write(IOut,'('' '')')
       Stop
      Else 
       call Coord(In,IOut,IPunch,IPrint,MxAt,MaxNZ,MxBnd,MxBox,PhyCon,
     $  KWd,Multip,NAtoms,NFrag,NBond,NH,NZ,IAn,Isot,IFrag,IBond,IANZ,
     $  IZ,MapZAt,Linear,C,CZ,EAN,EANZ,AtMass,TotWt,PMom,RotGHz,RTemp,
     $  MultN,QMom,GFac,NSpec,IScr,Scr,Group)
       NTRot=6*NFrag
       If(Linear) NTRot=5*NFrag 
      EndIf 
C Print Coordinates and Rotational Constants
      if(DoEck) then
       Write(IOut,'(/,10X,''Cartesian Coords.in Eckart Orientation'')')
      else
       Write(IOut,'(/,10X,''Cartesian Coords. in Input Orientation'')')
      endif
      If(NFrag.eq.1) then
       Write(IOut,'(4X,''Atom'',4X,''Pauling EAN  Topolog.EAN'',6X,
     $   ''X'',10X,''Y'',12X,''Z'')')
      Else
       Write(IOut,'(4X,''Atom'',4X,''Pauling EAN  Topolog.EAN'',6X,
     $  ''X'',10X,''Y'',12X,''Z'',5X,''Fragment'')')
      EndIf 
      Do 5 IAt=1,NAtoms
       If(NFrag.eq.1) then
        Write(IOut,'(I5,2X,A2,2(4X,F8.5),3F12.5)') IAt,IEl(IAn(IAt)),
     $   EAn(IAt),EANZ(IAt),(C(IXYZ,IAt),IXYZ=1,3)
       Else
        Write(IOut,'(I5,2X,A2,2(4X,F8.5),3F12.5,2X,I5)') IAt,
     $   IEl(IAn(IAt)),EAn(IAt),EANZ(IAt),(C(IXYZ,IAt),IXYZ=1,3),
     $   IFrag(IAt)
       EndIf
   5  Continue
      write(IOut,'(/,I5,'' Atoms with Total Mass ='',F10.3)') NAtoms,
     $  TotWt
      write(IOut,'('' Equilibrium Rotational Constants (MHz)'',
     $  3(1X,F12.5))') (1.0d+3*RotGHz(i),i=1,3)
      if(DVibRot) write(IOut,'('' Ground State Rotational Constants'',
     $  '' (MHz)  '',3F12.5)') (1.0d+3*RotGHz(i)+DBVib(i),i=1,3)
      write(IOut,'('' '')')
C Print interatomic distances
      Call LlinCl(LinScr)
      LinScr(1:22)=' Interatomic Distances (Angstrom)'
      Call HedPrt(IOut,0,LinScr,Num)
      ScalI=1.0d0
      Call DisMat(NAtoms,IAN,C,2,5,IOut,Error,0,ScalI)
      write(IOut,'('' '')')
C Make partition functions
      call Therm1(In,IOut,IPrint,PhyCon,Multip,NAtoms,Linear,
     $ TotWt,ZPE,RTemp,DBVib,Frq)         
C Make primitive bond lengths and valence angles including linear ones
      NLenR=0
      NAngR=0
      NLAngR=0
      NTotR=0
      call MkBAL(IOut,IPrint,MxBnd,MxAtP,MxTrm,NAtoms,NLenR,NAngR,
     $  NLAngR,IAn,NBond,IBond,NTermB,NTermA,NTermL,IAtmBR,IAtmAR,
     $  IAtmLR,CoefB,CoefA,CoefL,C,TreshL)
      NTotR=NLenR+NAngR+NLAngR
      if(NLenR.ne.0.and.DoBPCS) then
      write(IOut,'(/,'' Bonded Atoms  rDSD Distance BDPCS3 distance'')')
       do 10 il=1,NLenR
        IAt=IAtmBr(1,1,Il)
        JAt=IAtmBr(2,1,Il)
        R0IJ(IL)=Distan(C,IAt,JAt,0)   
        RBIJ(IL)=RBPCS(IAt,JAt,IAn,R0IJ(IL))
        write(IOut,'(2I5,7X,F8.5,7X,F8.5)') IAt,JAt,R0IJ(IL),RBIJ(IL) 
   10  continue
       write(IOut,'('' '')')
      endif 
C Make primitive dihedrals
      NDihR=0
      IPrDih=0
      call MkPrmD(IOut,IPrDih,MxBnd,MxTrm,MxTrm,MxAtP,MxAtP,MxAtCy,
     $  MxCyc,NAtoms,NBond,NLenR,NDihR,NTotR,NCyc,IBond,NTermD,IAtmBR,
     $  IAtmDR,NAtC,ICAt,IAtCyc,CoefD,C,TreshL)
      NTotR=NTotR+NDihR
C Make primitive out-of-plane bendings (use improper dihedrals for G16)
      NOuplR=0
      call MkPrmO(IOut,IPrOup,MxBnd,MxTrm,MxAtP,NAtoms,NBond,
     $  NOuplR,IBond,NTermO,IAtmOR,CoefO,C,ImpDih)
      NTotR=NTotR+NOuPlR
      if(NTotR.eq.0) NTotR=NLenR+NAngR+NLAngR+NDihR+NOuPlR
C Print summary of internal primitives
      If(.not.DoGNIC) then
       Write(IOut,'(/,I5,'' Atoms and'',I5,'' Internal Coordinates'')')
     $  NAtoms,3*NAtoms-NTRot
       Write(IOut,
     $ '(/,14X, '' Stretch.  Bend.  L. Bend. Tors.  Out-Pl. Total'')')
        Write(IOut,'('' Redundant  '',6I8)') NLenR,NAngR,NLAngR,NDihR,
     $  NOuPlR,NTotR
      EndIf
C Compute B matrix (not for G16)
      If(DoBPCS) then
       If(DoG16) then
        write(IOut,'('' BDPCS3 not available for G16'')')
        Stop
       EndIf
C      Call DrvBG(IOut,IPrint,Linear,DoBPCS,DoneC,MxAtP,MxTrm,NAtoms,
C    $   NLenR,NAngR,NLAngR,NOuplR,NDihR,NTotR,IAn,IAtmBR,IAtmAR,IAtmLR,
C    $   IAtmDR,IAtmOR,NTermB,NTermA,NTermL,NTermD,NTermO,C,Atmass,R0IJ,
C    $   RBIJ,CoefB,CoefA,CoefL,CoefD,CoefO,Scr,Iscr) 
      Call PCSGeo(IOut,IPrint,Linear,DoneC,MxAtP,MxTrm,NAtoms,NLenR,
     $  NAngR,NLAngR,NOuplR,NDihR,NTotR,IAn,IAtmBR,IAtmAR,IAtmLR,
     $  IAtmDR,IAtmOR,NTermB,NTermA,NTermL,NTermD,NTermO,C,Atmass,R0IJ,
     $  RBIJ,CoefB,CoefA,CoefL,CoefD,CoefO,Scr,Iscr)
       If(.not.DoneC) then
        write(IOut,'('' Failed Conversion to Cartesian Coordinates'')')
        Stop
       EndIf
       Write(IOut,'(/,'' BDPCS3 Cartesian Coordinates'')')
       Write(IOut,'(3X,''Atom'',2X,''At.Numb. Eff.At.Numb.'',5X,''X'',
     $    10X,''Y'',12X,''Z'')')
       Write(IOut,'(I5,I7,7X,F8.5,3F12.5)')(IAt,IAn(IAt),EAn(IAt),
     $   (C(IXYZ,IAt),IXYZ=1,3),IAt=1,NAtoms)
      Call MofI(IOut,Iprint,.false.,1,NAtoms,C,AtMass,XYZCM,TMom,PMom,
     $  RotMat)
C Compute rotational parameters
       FactA = CnvFct('FactA')
       Fac1  = PhyCon(1)**2
       AvPMom=1.0D0
       AvPMoB=1.0D0
       do 85 I=1,3
        If(PMom(i).le.0.0d0) goto 85
        AvPMom=AvPMom*PMom(i)
        PMomB(i)=PMom(i)/Fac1
        AvPMoB=AvPMoB*PMomB(i)
        Rotcm1(i)=FactA/PMom(i)
   85  Continue
       AvPMom=AvPMom**(1.0D0/3.0D0)
       AvPMoB=AvPMoB**(1.0d0/3.0d0)
       Call RotCon(IOut,IPrint,Linear,NAtoms,PhyCon,PMomB,RotGHz,RTemp)
       Write(IOut,'(/,'' Equilibrium Rotational Constants (MHz)'',
     $   3(1X,F12.5))') (1.0d+3*RotGHz(i),i=1,3)
       if (DVibRot) then
        write (IOut,'('' Ground State Rotational Constants (MHz)'',
     $    3(1X,F12.5))') (1.0d+3*RotGHz(ii)+DBVib(ii),ii=1,3)
       EndIf
       Write(IOut,'('' '')')
C Print interatomic distances
       Call LlinCl(LinScr)
       LinScr(1:22)=' Interatomic Distances (Angstrom)'
       Call HedPrt(IOut,0,LinScr,Num)
       ScalI=1.0d0
       Call DisMat(NAtoms,IAN,C,2,5,IOut,Error,0,ScalI)
       write(IOut,'('' '')')
      EndIf
      If(.not.DoGNIC) then
       goto 999
      Else
       write(IOut,'(/,'' Definition of Generalized Natural Internal'',
     $   '' Coordinates'')')
      EndIf
C Make cycles
      NExpCy=NLenR-NAtoms+NFrag
      IPrCyc=0
      call MkCyc(IOut,IPrCyc,MxBnd,MxAtP,MxTrm,
     $  MxAtCy,MxCyc,NBond,IBond,NLenR,NAngR,NDihR,NExpCy,
     $  IAtmBR,IAtmAR,IAtmDR,NCyc,NatC,ICAt,IAtCyc,IBr,IAn,EAN)
      If(NCyc.gt.0.or.NExpCy.gt.0) write(IOut,'(I3,
     $  '' Rings Found Over'',I3,'' Expected'',/)') NCyc,NExpCy
C Make bond GNICs: retain individual bonds except possibly for terminal 
C Atoms (if DoSysy.eq.true.) and cycles
      NLen=0
      IPrBnd=0
      Call MkGNCB(IOut,IPrBnd,DoSySt,InvDst,MxBnd,MxTrm,MxAtP,NAtoms,
     $  IAn,NBond,NLen,IBond,NTermB,IAtomB,ITVB,IAtCyc,CoefB,C)
C     If(NCyc.gt.0.and.DoSySt) then
C      IPrtCB=0
C Cycles NYI for stretchings
C      Call CySalc(IOut,IprtCB,MxAtP,MxAtCy,MxTrm,NAtoms,NLen,
C    $   NLenR,NatC,ICAt,IAtCyc,IAtmBr,IAtomB,NTermB,ITVB,CoefB,EAn,C)
C     EndIf
C Make angle GNICs
C      NAng=0
      IPrAng=0
      call MkGNCA(IOut,IPrAng,MxBnd,MxGNIC,MxTrm,MxAtP,NAtoms,
     $  NCyc,NBond,NAng,IBond,NTermA,IAtomA,IAn,IAtCyc,ITVA,
     $  CoefA,C,EAN,TreshL)
C Make ring coordinates for valence angles
      If(NCyc.gt.0) then
       IPrtCA=0
       NAng00=NAng
       do 20 icyc=1,NCyc
        NAng0=NAng
        call CyGNA(IOut,IPrtCA,MxAtCy,MxAtP,MxTrm,NBond,NAtC,NAng,ICyc,
     $    ICAt,NTermA,IAtomA,ITVA,CoefA)
        If(NAtC(ICyc).eq.3) then
         write(IOut,'(''   3-Membered Ring ('',3I4,'' ) '',39X,I2)') 
     $     (ICAt(ii,ICyc),ii=1,NAtC(ICyc)),NAng-NAng0 
        ElseIf(NAtC(ICyc).eq.4) then
         write(IOut,'(''   4-Membered Ring ('',4I4,'' ) '',35X,I2)') 
     $     (ICAt(ii,ICyc),ii=1,NAtC(ICyc)),NAng-NAng0
        ElseIf(NAtC(ICyc).eq.5) then
         write(IOut,'(''   5-Membered Ring ('',5I4,'' ) '',31X,I2)') 
     $     (ICAt(ii,ICyc),ii=1,NAtC(ICyc)),NAng-NAng0
        ElseIf(NAtC(ICyc).eq.6) then
         write(IOut,'(''   6-Membered Ring ('',6I4,'' ) '',27X,I2)') 
     $     (ICAt(ii,ICyc),ii=1,NAtC(ICyc)),NAng-NAng0
        ElseIf(NAtC(ICyc).eq.7) then
         write(IOut,'(''   7-Membered Ring ('',7I4,'' ) '',23X,I2)') 
     $     (ICAt(ii,ICyc),ii=1,NAtC(ICyc)),NAng-NAng0
        ElseIf(NAtC(ICyc).eq.8) then
         write(IOut,'(''   8-Membered Ring ('',8I4,'' ) '',19X,I2)')
     $     (ICAt(ii,ICyc),ii=1,NAtC(ICyc)),NAng-NAng0
        EndIf
   20  continue
      EndIf
C Make Linear Angle GNICs: make 2 angles for linear molecules
      NLAng=0
      IPrLAn=0
      call MkGNLA(IOut,IPrLAn,MxBnd,MxGNIC,MxTrm,MxAtP,NAtoms,
     $  NBond,NLAng,Linear,IBond,NTermL,IAtomL,IAn,CoefL,C,TreshL)
C Make dihedral GNICs
      NDih=0
      call MkGNCD(IOut,IPrint,MxBnd,MxTrm,MxTrm,MxAtP,MxAtP,MxAtCy,
     $  Do1Dih,NAtoms,IAn,NBond,NLenR,NDih,NTot,NCyc,IBond,NTermD,
     $  IAtmBR,IAtomD,IBr,NAtC,ICAt,IAtCyc,ITVD,IPerD,NEqAt,CoefD,C,
     $  EAN,TreshL,DoNorm)
C Make ring coordinates for dihedra angles
      If(NCyc.gt.0) then
       NDihCh=NDih
       IPrtCD=0
       do 30 ICyc=1,NCyc
        call CyGND(IOut,IPrtCD,MxAtCy,MxAtP,MxTrm,NAtC,NDih,ICyc,ICAt,
     $    NTermD,IAtomD,ITVD,CoefD)
   30  continue
       write(IOut,'(/,I3,'' Ring coordinates for Dihedral Angles'')')
     $    NDih-NDihCh
      endif
C Make Out-of-Plane GNICs
      NOupl=0
      call MkGNCO(IOut,IPrint,.True.,MxBnd,MxGNIC,MxTrm,MxAtP,NAtoms,
     $  NCyc,NBond,NOuPl,IBond,NTermO,IAtomO,IAn,IAtCyc,CoefO,C,ImpDih)
C
      NTot=NLen+NAng+NLAng+NDih+NOUPl
CENZO Print Information on Torsions
      Indd=0
      do IPuf=1,NLenR
       IScr(indd+1)=IAtmBr(1,1,IPuf)
       Iscr(Indd+2)=IAtmBr(2,1,IPuf)
       Indd=Indd+2
      enddo
      Call DrvTrs(IOut,IPrint,MxBnd,NAtoms,NLenR,IAn,IAtCyc,NBond,IBond,
     $IScr,EAn,C)
CENZO 
C Print results
      Write(IOut,'(/,I5,'' Atoms and'',I5,'' Internal Coordinates'')')
     $  NAtoms,3*NAtoms-NTRot 
      Write(IOut,
     $ '(/,14X, '' Stretch.  Bend.  L. Bend. Tors.  Out-Pl. Total'')')
      Write(IOut,'('' Redundant  '',6I8)') NLenR,NAngR,NLAngR,NDihR,
     $  NOuPlR,NTotR
      Write(IOut,'('' Non Redund.'',6I8,/)') NLen,NAng,NLAng,NDih,
     $  NOuPl,NTot
      NRed=NTot-3*NAtoms+NTRot
      If(NRed.eq.0) then
       write(IOut,'('' All the Redundancies have been Eliminated'',/)')
      else
       write(IOut,'('' WARNING:'',I3,'' Redundancies Still Present'',
     $   /)') NRed 
      endif        
C Possibly prepare the input for VOLT
      NTT=0
C First part of Volt input
      IVlt=0
      if(DoVolt) then
       if(DoG16) then
        write(IOut,'('' No Volt Input for G16'')')
        Stop
       endif
       OPEN(IPunch,FILE='voltin',STATUS='UNKNOWN')
       Rewind(IPunch)
       Write(IPunch,'(/,A80,/)') Title
       do 35 IAt=1,NAtoms
        Write(IPunch,'(I3,3F12.5)') IAn(IAt),(C(ii,IAt),ii=1,3)
   35  Continue
       Write(IPunch,'('' '')')
       NPV=0
       IVlt=IPunch
       if(NLenR.ne.0)call PrtPrm(IVlt,MxAtP,MxTrm,NLenR,1,IAtmBR,NPV)
       if(NAngR.ne.0)call PrtPrm(IVlt,MxAtP,MxTrm,NAngR,2,IAtmAR,NPV)
       if(NLAngR.ne.0)call PrtPrm(IVlt,MxAtP,MxTrm,NLAngR,3,IAtmLR,NPV)
       if(NDihR.ne.0)call PrtPrm(IVlt,MxAtP,MxTrm,NDihR,4,IAtmDR,NPV)
       if(NOuplR.ne.0)call PrtPrm(IVlt,MxAtP,MxTrm,NOuplR,5,IAtmOR,NPV)
       write(IPunch,'('' '')')
      endif
C Bond Lengths
      IType=1
      Ini=1
      NVar=NLen
      IniP=1
      NVarP=NLenR
      If(DoBPCS) write(IOut,'('' Bond Lengths corrected by BDPCS3'',/)')
      IFill=0
      If(RigB) IFill=1
      Do 40 Ir=1,NLen
       IFixB(Ir)=IFill        
   40 Continue
      call OrdRed(IOut,IVlt,IPrint,MxAtP,MxTrm,DoBPCS,IType,InvDst,NVar,
     $  Ini,IniP,NTermB,IAtomB,IPrimB,ITVB,IFixB,IAn,CoefB,ValTB,C,
     $  ImpDih,Clean)
C Valence Angles
      IType=2
      Ini=NLen+1
      NVar=NAng 
      IniP=NLenR+1
      NVarP=NAngR
      IFill=0
      If(RigA) IFill=1
      Do 50 Ir=1,NAng
       IFixA(Ir)=IFill
   50 Continue
      call OrdRed(IOut,IVlt,IPrint,MxAtP,MxTrm,DoBPCS,IType,.False.,
     $ NVar,Ini,IniP,NTermA,IAtomA,IPrimA,ITVA,IFixA,IAn,CoefA,ValTA,C,
     $ ImpDih,Clean)
C Linear Valence Angles
      IType=3
      Ini=NLen+NAng+1
      NVar=NLAng
      IniP=NLenR+NAngR+1
      NVarP=NLAngR
      IFill=0
      If(RigL) IFill=1
      Do 60 Ir=1,NLAng
       IFixL(Ir)=IFill
   60 Continue
      call OrdRed(IOut,IVlt,IPrint,MxAtP,MxTrm,DoBPCS,IType,.False.,
     $ NVar,Ini,IniP,NTermL,IAtomL,IPrimL,ITVLA,IFixL,IAn,CoefL,ValTL,C,
     $ ImpDih,Clean)
C Dihedrals
      Itype=4
      Ini=NLen+NAng+NLAng+1 
      NVar=NDih
      IniP=NLenR+NAngR+NLAngR+1
      NVarP=NDihR
      IFill=0
      If(RigD) IFill=1
      Do 70 Ir=1,NDih
       IFixD(Ir)=IFill
   70 Continue
      call OrdRed(IOut,IVlt,IPrint,MxAtP,MxTrm,DoBPCS,Itype,.False.,
     $ NVar,Ini,IniP,NTermD,IAtomD,IPrimD,ITVD,IFixD,IAn,CoefD,ValTD,C,
     $ ImpDih,Clean)
C Out of Plane
      Itype=5
      Ini=NLen+NAng+NLAng+NDih+1
      NVar=NOuPl
      IniP=NLenR+NAngR+NLAngR+NDihR+1
      NVarP=NOuPlR
      IFill=0
      If(RigO) IFill=1
      Do 80 Ir=1,NOupl
       IFixO(Ir)=IFill
   80 Continue
      call OrdRed(IOut,IVlt,IPrint,MxAtP,MxTrm,DoBPCS,Itype,.False.,
     $ NVar,Ini,IniP,NTermO,IAtomO,IPrimO,ITVO,IFixO,IAn,CoefO,ValTO,C,
     $ ImpDih,Clean)
      if(DoVolt) close(ipunch)
C Gaussian Input
      if(DoG16.or.DoGDV) then
       OPEN(IPunch,FILE='gauin',STATUS='UNKNOWN')
       Rewind(IPunch) 
       call SetGKw(IPunch,DoGNIC,SyGNIC,Loose,DoScan,ModPCS,IDeriv)
       Write(IPunch,'(/,A80,/)') Title
       Write(IPunch,'(2I3)') ICharg,Multip
       do 90 IAt=1,NAtoms
        If(Loose) then
         Write(IPunch,'(I3,3F9.2)') IAn(IAt),(C(ii,IAt),ii=1,3)
        Else
         Write(IPunch,'(I3,3F12.6)') IAn(IAt),(C(ii,IAt),ii=1,3)
        EndIf
   90  continue
       Write(IPunch,'('' '')')
C Bond Lengths
CENZO test
       PrtVal=.False.
CENZO
       Call PrtBnd(IPunch,MxAtP,MxTrm,InvDst,NLen,NTT,NTermB,IAtomB,
     $   ITVB,IFixB,IAn,CoefB,ValtB,C,PrtVal)
C Valence Angles
       Call PrtAng(IPunch,MxAtP,MxTrm,NAng,NTT,NTermA,IAtomA,ITVA,
     $  IFixA,CoefA,ValtA,C,PrtVal)
C Linear Angles
       Call PrtLAn(IPunch,MxAtP,MxTrm,NLAng,NTT,Linear,NTermL,IAtomL,
     $   ITVLA,IFixL,CoefL,ValtL,C,PrtVal)
C Dihedral Angles
       Call PrtDih(IPunch,MxAtP,MxTrm,NDih,NTT,DoScan,NTermD,IAtomD,
     $   ITVD,IPerD,IFixD,CoefD,ValtD,C,Clean,PrtVal)
C Out of Plane Angles(U) or Impr. Dihedrals (Only option in g16)      
       Call PrtOut(IPunch,MxAtP,MxTrm,NOuPl,NTT,NTermO,IAtomO,ITVO,
     $   IFixO,CoefO,ValtO,C,ImpDih,PrtVal)
C Fragments
C     If(NFrag.gt.1) call OrdFrg(IPunch,IPrint,MxAtFr,NAtoms,NFrag,
C    $  IFrag,NAtFr,IFrsAt,IlstAt,IAtFr,LConn)
       If(NFrag.gt.1.and.DoGDV) write(IPunch,'('' Rotor(1)'')')
C rDSD
       If(ModPCS.eq.3) write(IPunch,'(/,''@3F12red.gbs'',/)')
       close(ipunch)
      EndIf
  999 Continue
C Possibly Compute Substitution Structure
      If(DoVMSR) then
C Read NSpec, Errors,DeltaVib,DeltaEl
C The Default is to read only B0 values
      call LlinCl(LinScr)
      read(In,'(A80)') LinScr
      call ParseFlags(LinScr, NSPEC, RdB0Er, RdVib, RdEle)
      RdB0=.true.
      RdDvEr=.false.
      RdDeEr=.false.      
C AtMass
       IAtm=1
C B0
       IB0=IAtm+3*NSpec
C ErrB0
       IErrB0=IB0+3*NSpec
C DBVib
       IDBVib=IErrB0+3*NSpec
C ErrVi
       IErrVi=IDBVib+3*NSpec
C DBEle 
       IDBEle=IErrVi+3*NSpec
C ErrEl
       IErrEl=IDBEle+3*NSpec
C Scratch
       IEnd=IErrEl+3*NSpec
C Compute Substitution Structure
       call VMSRot(In,IOut,IPrint,NAtoms,NSpec,RdB0,RdB0Er,RdVib,
     $  RdDvEr,RdEle,RdDeEr,IAn,AtMass,C,Scr(IAtm),Scr(IB0),Scr(IErrB0),
     $  Scr(IDBVib),Scr(IErrVi),Scr(IDBEle),Scr(IErrEl),Scr(IEnd),IEnd)
      EndIf 
C Possibly Prepare input for MSR Program
      If(RdIsot) then
       OPEN(IPunch,FILE='msrin',STATUS='UNKNOWN')
  788  Read(IPunch,*,end=789) StrInp
       goto 788
  789  backspace(IPunch)     
C Atmass
       IAtm=1
C B0
       IB0=IAtm+3*NSpec 
C IDBVib
       IDBVib=IB0+3*NSpec
C IDBEle
       IDBEle=IDBVib+3*NSpec
C IEnd
       IEnd=IDBEle+3*NSpec
C
       call SemiEx(In,IPunch,IPrint,NAtoms,NSpec,RdB0,RdVib,RdEle,IAn,
     $  Scr(IAtM),Scr(IB0),Scr(IDBVib),Scr(IDBEle))
C      close(IPunch)
      EndIf
C Compute Rates
      Call Rates(IOut,IPrint,DoColl,DoGor,DoSpin)
C
      write(*,'('' Normal Termination of DiNa25'')')
      write(*,'('' DiNa25   output on file xxx.out'')')
      If(RdSMI)  write(*,'('' AVOGADRO input  on file xxx.xyz'')')
      If(DoG16)  write(*,'('' G16-C01  input  on file xxx.gjf'')')
      If(DoGDV)  write(*,'('' GDV-J28  input  on file xxx.gjf'')')
      If(DoVolt) write(*,'('' VOLT     input  on file xxx.vlt'')')
      If(RdIsot) write(*,'('' MSRS     input  on file xxx.msr'')')
      End
C===============================================================
C  DRIVER per il calcolo di energie torsionali
C===============================================================
*Deck DrvTrs
      Subroutine DrvTrs(IOut,IPrint,MxBnd,NAtoms,NLen,IAn,IRing,NBond,
     $  IBond,IAtom,EAn,C)
      INTEGER IOut, IPrint, MxBnd, NAtoms, NLen
      INTEGER IAn(NAtoms), IRing(NAtoms), IAtom(2,NAtoms)
      INTEGER NBond(NAtoms), IBond(MxBnd, NAtoms)
      DOUBLE PRECISION EAn(NAtoms),C(3, NAtoms)
C     Write(IOut,'(/,I5,'' Atoms'')') NAtoms
C     do 10 IAt=1,NAtoms
C      write(IOut,'('' Atom'',I5,'' Ian ='',I3,'' EAn ='',F8.4,
C    $   '' forms'',I2,'' bonds with atoms'',4I3)') IAt,IAn(IAt),
C    $   EAn(IAt),NBond(IAt),(IBond(IB,IAt),IB=1,NBond(IAt))
C 10  Continue 
C     Write(IOut,'('' IRing'')')
C     Write(IOut,'(20I3)') (IRing(I),I=1,NAtoms)
C     Write(IOut,'('' Atoms Involved in'',I5,'' Bonds'')') NLen
C     Write(IOut,'(2I3)') (Iatom(1,ILen),IAtom(2,ILen),ILen=1,NLen)
      CALL SCAN_ALL_BONDS_RIGID(IOut, IPrint,
     &     MxBnd, NAtoms, IAn, EAn, IRing, NLen, IAtom,
     &     NBond, IBond, C)

      Return
      END
*Deck SetIO
      Subroutine SetIO(InFil,OutFil,Ext)
      Character*(*) InFil, OutFil
      Character*3 Ext
      IStart= LineSt(InFil)
      IEnd=IStart+LinEnd(InFil)
      OutFil=InFil(IStart:IEnd)
      OutFil(IEnd:IEnd)='.'
      OutFil(IEnd+1:IEnd+3)=Ext
      Return 
      End 
*Deck FndKwd
      Subroutine FndKwd(IOut,IPrint,ModPCS,IDeriv,CLine,Kwd,MxKwd)
      Implicit Real*8 (A-H,O-Z)
      Character*80  CLine
      Character*8   Test
      Logical Kwd(MxKwd),Error,FndPCS,FndDer
      Dimension Ini(50)
      Dimension IDeriv(4)
      call LClear(MxKwd,Kwd)
      FndPCS=.false.
      FndDer=.false.
      ModPCS=-1
      If(CLine(1:1).eq.'#') then
       CLine(1:1)=' '
      Else
       write(IOut,'('' No Input for DiNa25'')')
       Stop
      EndIf
      If(Cline(2:3).eq.'P ') then
       CLine(2:2)=' '
       IPrint=1
      EndIf
      call SubStr(CLine,MxKwd-1,Ini,NumKwd)
      Ini(NumKwd+1)=len(CLine)
      Error=.False.
      do 10 i=1,NumKwd
       I1=Ini(i)
       If(I.eq.NumKwd) then
        I2=I1+7 
       Else
        I2=Min0(I1+7,Ini(i+1)-1)
       EndIf
C the default is to normalize GNICs
       Kwd(12)=.true.
       Call LinUpC(CLine(I1:I2),Test)
       If(Test(1:5).eq.'PRINT') then
        IPrint=1
       ElseIf(Test(1:6).eq.'SMILES') then
        Kwd(1)=.True.
       ElseIf(Test(1:3).eq.'G16') then
        Kwd(2)=.True.
        Kwd(7)=.True.
       ElseIF(Test(1:3).eq.'GDV') then  
        Kwd(3)=.True.
       ElseIf(Test(1:5).eq.'CUBIC') then
        Kwd(4)=.True. 
       ElseIf(Test(1:7).eq.'NATURAL') then
        Kwd(5)=.True.
       ElseIf(Test(1:6).eq.'ECKART') then
        Kwd(6)=.True.
       ElseIf(Test.eq.'IMPDIH') then
        Kwd(7)=.True.
       ElseIf(Test(1:4).eq.'GNIC') then
        Kwd(8)=.True.
       ElseIf(Test(1:7).eq.'INVDIST') then
        Kwd(9)=.True.
       ElseIf(Test(1:7).eq.'SYMMSTR') then
        Kwd(10)=.True.
       ElseIf(Test(1:6).eq.'ONEDIH') then
        Kwd(11)=.True.
       ElseIf(Test(1:6).eq.'NONORM') then
        Kwd(12)=.False.
       ElseIf(Test(1:6).eq.'BDPCS3') then
        Kwd(13)=.True.
       ElseIf(Test(1:7).eq.'SYMMALL') then
        Kwd(14)=.True.
       ElseIf(Test(1:6).eq.'FINDFR') then
        Kwd(15)=.True.
       ElseIf(Test(1:6).eq.'JOINFR') then
        Kwd(16)=.True.
       ElseIf(Test(1:5).eq.'HBOND') then
        Kwd(17)=.True.
       ElseIf(Test.eq.'VOLT') then
        Kwd(18)=.True.
       ElseIf(Test(1:4).eq.'BMAT') then
        Kwd(19)=.True.
       ElseIf(Test(1:4).eq.'SCAN') then
        Kwd(20)=.True.
       ElseIf(Test(1:5).eq.'RIGID') then
        Kwd(21)=.True.
       ElseIf(Test(1:4).eq.'BLDZ') then
        Kwd(22)=.True. 
       ElseIf(Test(1:4).eq.'WRTZ') then
        Kwd(23)=.True.
       ElseIf(Test(1:6).eq.'SEMIEX') then
        Kwd(24)=.True.
       ElseIf(Test(1:6).eq.'FITPOT') then
        Kwd(25)=.True.
       ElseIf(Test(1:4).eq.'VCI1') then
        Kwd(26)=.True.
       ElseIf(Test(1:4).eq.'DVR1') then
        Kwd(27)=.True.
       ElseIf(Test(1:5).eq.'CLEAN') then 
        Kwd(31)=.true.  
       ElseIf(Test(1:5).eq.'DVIBROT') then
        Kwd(32)=.true.
       ElseIf(Test(1:6).eq.'ISOTOP') then
        Kwd(37)=.true.
       ElseIf(Test(1:5).eq.'LOOSE') then
        Kwd(33)=.true.
       ElseIf(Test(1:5).eq.'RDXYZ') then
        Kwd(34)=.true.
       ElseIf(Test(1:6).eq.'RDFCHK') then
        Kwd(35)=.true.
       ElseIf(Test(1:5).eq.'RDSMI') then
        Kwd(36)=.true.
       ElseIf(Test(1:3).eq.'OPT') then
        IDeriv(1)=1
       ElseIf(Test(1:4).eq.'FREQ') then
        IDeriv(2)=1
       ElseIf(Test(1:4).eq.'CUBIC') then
        IDeriv(3)=1
       ElseIf(Test(1:6).eq.'ANHARM') then
        IDeriv(4)=1
       ElseIf(.not.FndPCS) then
        call SetPCS(Test,ModPCS,FndPCS)
        If(ModPCS.lt.0) Error=.true.
       EndIf
C       If(.not.FndDer) Call SetDer(Test,IDeriv,FndDer)
C       ITDer=IDeriv(1)+IDeriv(2)+IDeriv(3)+IDeriv(4)
C       If(ModPCS.eq.-1.and.ITDer.eq.0) Error=.true.
C       EndIf
   10 continue
C Avoid contemporary SemiEx and substitution Structures
      If(Kwd(24)) Kwd(37)=.False.
      If(Error) then
       write(IOut,'(/,'' Spurious Keywords in the Input '')')
C      write(IOut,'(A80)') CLine
       write(IOut,'('' The Following Keywords are Allowed'',/)')
       do 20 i=1,MxKwd
        Kwd(i)=.true.
   20  continue
      EndIf
      If(Kwd(1))  write(IOut,'('' SMILES    : SMILES by RdKIT'')') 
      If(Kwd(36)) write(IOut,'('' RDSMI     : SMILES by GICForge'')')       
      If(Kwd(34)) write(IOut,'('' RDXYZ     : Coords.from XYZ File'')')
      If(Kwd(35)) write(IOut,'('' RDFCHK    : Coords.from FCHK File'')')
      If(Kwd(2))  write(IOut,'('' G16       : Make G16 Input'')')
      If(Kwd(3))  write(IOut,'('' GDV       : Make GDV Input'')')
      If(Kwd(4))  write(IOut,'('' CUBIC     : Freq=Cubic(GDV) and'',
     $ '' Vibr=Gauconv MSR'')')
      If(Kwd(5))  write(IOut,'('' NATURAL   : Natural Isotopic'',
     $ '' Abundance '')')
      If(Kwd(6))  write(IOut,'('' ECKART    : Enforce Eckart'',
     $  '' Orientation'')')
      If(Kwd(7))  write(IOut,'('' IMPDIH    : Use Improper Dihedrals'', 
     $   '' in place of Out-of-plane Bends'')')
      If(Kwd(8))  write(IOut,'('' GNIC      : Make Generalized'',
     $  ''  Natural Internal Coords.'')')
      If(Kwd(9))  write(IOut,'('' INVDIST   : Inv.Dist. for Stretching''
     $  )')
      If(Kwd(10)) write(IOut,'('' SYMMSTR   : Symmetrize Stretching'')')
      If(Kwd(11)) write(IOut,'('' ONEDIH    : 1 Dihedral per Bond'')')
      If(.not.Kwd(12)) write(IOut,'('' NONORM    : Not Normalize'',
     $  '' Dihedral GNICS'')')
      If(Kwd(31)) write(IOut,'('' CLEAN     : Clean GNIC values'')')
      If(Kwd(33)) write(IOut,'('' LOOSE     : Loose Symmetry for'',
     $  '' Gaussian'')')
      If(Kwd(13)) write(IOut,'('' BDPCS3    : Make BDPCS3 Bond'',
     $  '' Lengths'')')
      If(Kwd(14)) write(IOut,'('' SYMMALL   : Symmetrize GNICs'')')
      If(Kwd(15)) write(IOut,'('' FINDFR    : Find Fragments'')')
      If(Kwd(16)) write(IOut,'('' JOINFR    : Join Fragments'')')
      If(Kwd(17)) write(IOut,'('' HBOND     : Make H-Bonds'')')
      If(Kwd(18)) write(IOut,'('' VOLT      : Make VOLT Input'')')
      If(Kwd(19)) write(IOut,'('' BMAT      : Build B Matrix'')')
      If(Kwd(20)) write(IOut,'('' SCAN      : Scan for Soft DOF'')')
      If(Kwd(21)) write(IOut,'('' RIGID     : Freeze Hard Modes'',
     $  '' in Scan'')')
      If(Kwd(22)) write(IOut,'('' BLDZ      : Build ZMAT from'',
     $  '' Input XYZ'')') 
      If(Kwd(23)) write(IOut,'('' WRTZ      : Write ZMAT for MSR'')')
      If(Kwd(24)) write(IOut,'('' SEMIEX    : Read Isotopes and '',
     $ ''make MSR input (requires ZMAt)'')')
      If(Kwd(32)) write(IOut,'('' DVIBROT   : Read Vibr.Corr. to'', 
     $  '' Rot.Const.'')')
      If(Kwd(37)) write(IOut,'('' ISOTOP    : Read Isotopes and '',
     $  ''Compute Substitution Structure'')')
      If(Kwd(25)) write(IOut,'('' FITPOT    : 1D-Polynomial Fit'')')
      If(Kwd(26)) write(IOut,'('' VCI1      : 1D-vibrational CI'')')
      If(Kwd(27)) write(IOut,'('' DVR1      : 1D-vibrational DVR'')')
      If(Kwd(28)) write(IOut,'('' COLL      : Collision Rate'')')
      If(Kwd(29)) write(IOut,'('' GORIN     : Gorin Barrierless'')')
      If(Kwd(30)) write(IOut,'('' SPINFOR   : Spin Forbidden'')')
      If(ModPCS.eq.0.or.error)write(IOut,'('' PCS0      : PCS0 Model'',
     $  '' (UFF)'')')
      If(ModPCS.eq.1.or.error)write(IOut,'('' PCS1      : PCS1 Model'',
     $  '' (HF3C)'')')
      If(ModPCS.eq.2.or.error)write(IOut,'('' HPCS2     : HPCS2 Model'',
     $  '' (B3LYP-D4/6-31G*)'')')
      If(ModPCS.eq.3.or.error)write(IOut,'('' DPCS3     : DPCS3 Model'',
     $  '' (rev-dsd-PBEP86D4/3F12-)'')')
      If(IDeriv(1).eq.1.or.error)write(IOut,'('' OPT       : Geometry'',
     $  ''  Optimization'')') 
      If(IDeriv(2).eq.1.or.error)write(IOut,'('' FREQ      : Second'',
     $  '' Energy Derivatives'')')
      If(IDeriv(3).eq.1.or.error)write(IOut,'('' CUBIC     : Third'',
     $  '' Energy Derivatives'')')
      If(IDeriv(4).eq.1.or.error)write(IOut,'('' ANHARM    : Fourth'',
     $  '' Energy Derivatives (up to ijkk'')')       
      If(Error) STOP
      Return
      End
*Deck NewVib
      Subroutine NewVib(In,IOut,IPrint,DoFit,DoVCI,DoDVR,MxScr,Scr,
     $  IScr) 
      Implicit Real*8 (A-H,O-Z)  
      Logical DoFit,DoVCI,DoDVR,RdData
      Dimension Scr(MxScr),IScr(MxScr)
      Dimension IPot(20),Cof(20),Delt(20),F(20),G(20),Red(20)
      Common /PhyCon/ PhyCon(30)
C
      Read(In,*) IMode,NPts,NPrm,NSig,NBas,NEig
C Possibly fit a potential
C Scr(IX) = X
      IX = 1
C Scr(IY) = Y
      IY = IX + NPts
C SCr(IS) = Sigma
      IS = IY + NPts
C Scr(IEnd) = Scratch
      IAvail = IS + NPts
C
      If(NPts.eq.0) then
       If(DoFit) then
        write(IOut,'('' NewVib: No Fit with 0 Data Points'')')
        Stop
       EndIf
C Read F,G parameters for VCI
       Read(In,*) (F(I),I=1,NPrm)
       Read(In,*) (G(I),I=1,NPrm)
       Write(IOut,'(/,I3,'' Potential and Kinetic Parameters'')')NPrm
C      Write(IOut,'(6D12.5)') (F(i),I=1,NPrm)
C      Write(IOut,'(6D12.5)') (G(I),I=1,NPrm)
      Else
       If(DoVCI) then
        write(IOut,'('' NewVib: iVCI1 with given points NYI'')')
        Stop
       EndIf 
C Read IPot,Cof parameters for Polynomial fitting
       Read(In,*) (IPot(i),I=1,NPrm)
       Read(In,*) (Cof(i),I=1,NPrm)
C      Write(IOut,'(/,'' Polynomial Fitting with '',I3,'' Terms'')')NPrm
C      Write(IOut,'(6I3)') (IPot(i),I=1,NPrm)
C      Write(IOut,'(6D12.5)') (Cof(I),I=1,NPrm)
      EndIf
      If(NPts.gt.0) then
C Read X,Y,Sigma
       do 10 I=1,NPts
        II=I-1
        If(NSig.gt.0) then
         Read(In,*) Scr(IX+II),Scr(IY+II),Scr(IS+II) 
        Else
         Read(In,*) Scr(IX+II),Scr(IY+II) 
         Scr(IS+II)=1.0d0
        EndIf
   10  continue
      EndIf 
      IPrint=1
      If(DoFit.and.NPts.gt.0) then
       Call NewFit(In,IOut,IPrint,NPts,Scr(IX),Scr(IY),Scr(IS),NPrm,
     $    IPot,Cof,Delt,MxScr,Scr(IAvail),IScr)
       Write(IOut,'(/'' Fitted Polynomial'')')
       Do 100 I=1,NPrm 
        If(IPot(I).eq.2) Freq=SQrt(CnvFct('Fac2AU')*Cof(I)/2.0d0)
        If(I.lt.10) then
         Write(IOut,'(D12.5,'' X**'',I1)') Cof(I),IPot(I)
        Else
         Write(IOut,'(D12.5,'' X**'',I2)') Cof(I),IPot(I)
        EndIf
  100  Continue
C
C Test for the transformation from a.u. to reduced
C
       Call Aclear(NPrm,Red)
       FG=CnvFct('FactG')
       F3=CnvFCt('Fac3AU')
       F4=CnvFct('Fac4AU')
       F4N=F3/(Sqrt(FG)*PhyCon(1))
       Write(IOut,'(/,''Fact4'',2D12.5,/)') F4,F4N
       Do 110 I=1,NPrm
        If(IPot(I).eq.2) then
         Red(2)=Freq
         write(IOut,'('' Half Harmonic Frequency'',F12.5)') Red(2)
        ElseIf(IPot(I).eq.3) then
         SNum3=CnvFct('Fac3AU')
         Den3=6.0d0*Sqrt(Freq**3)
         Red(3)=Cof(I)*SNum3/Den3
         write(IOut,'('' Cubic Anharmonicity'',F12.5)') Red(3)
        ElseIf(IPot(I).eq.4) then
         SNum4=CnvFct('Fac4AU')
         Den4=24.0d0*Sqrt(Freq**4)
         Red(4)=Cof(I)*Snum4/Den4
         write (IOut,'('' Quartic Anharmonicity'',F12.5)') Red(4)
        ElseIf(IPot(I).eq.5) then
         SNum4=CnvFct('Fac4AU')
         SNum5=SNum4/(Sqrt(FG)*PhyCon(1)**2)
         Den4=24.0d0*Sqrt(Freq**4)
         Den5=5.0D0*Den4*SQrt(Freq)
         Red(5)=Cof(I)*Snum5/Den5
         write (IOut,'('' Quintic Anharmonicity'',F12.5)') Red(5)
        ElseIf(IPot(I).eq.6) then
         SNum4=CnvFct('Fac4AU')
         SNum6=SNum4/(FG*PhyCon(1)**4)
         Den4=24.0d0*Sqrt(Freq**4)
         Den6=5.0D0*6.0D0*Den4*Freq
         Red(6)=Cof(I)*Snum6/Den6
         write (IOut,'('' Sextic Anharmonicity'',F12.5)') Red(6)
        EndIf
  110  Continue   
C      Call SetPts(IOut,IPrint,NPts,NBas,NPrm,NEig,Scr(IX),Scr(IY),
C    $   IPot,Cof,Scr(IAVail),Scr(IAVail+NBas))
C      Call AMove(NBas,Scr(IAVail),SCr(IX))
C      Call AMove(NBas,Scr(IAvail+NBas),Scr(IY))
C      NPts=NBas   
      EndIf 
C Variational VCI1
      If(DoVCI) then
C for the moment does not use the polynomial fitting
       Alpha=1.0D0
       IAvail = 1
       Call VCI1D(In,IOut,IPrint,NBas,NEig,Alpha,F,G,MxScr,
     $  Scr(IAvail),IScr)
       Return 
      EndIf
C DVR1
      If(DoDVR) then
       IAvail=IS
       call DVR1D(In,IOut,IPrint,IMode,NPts,NBas,NEig,SCr(IX),
     $   SCr(IY),MxScr,Scr(IAvail),IScr)
      EndIf
      Return
      End
*Deck NewFit
      Subroutine NewFit(In,IOut,IPrint,NPts,X,Y,Sig,NPrm,IPot,Cof0,
     $   Delt0,MxScr,V,IV)
      Implicit Real*8 (A-H,O-Z)
C ---------------------------------------------------------------------
C Input
C NPts = Number of Points
C NPrm = number of parameters to be fitted
C NSig = if 0 no errors
C X(NPts) = abscissas
C Y(NPts) = ordinates
C Sig(NPts) = errors
C IPot(NPrm) = powers of the polynomial terms
C Cof0(NPrm) = values of the fixed parameters (0.0 for free)
C
C Output
C Cof(NPrmf)  = Free Coefficients of the LSQ Polynomial
C Delt(NPrmf) = Estimated Errors of the Free Coefficients
C Cof0(NPrm)  = Fixed + Free Coefficients
C ---------------------------------------------------------------------
      Dimension IPot(*),Cof0(*),Delt0(*),V(*),IV(*)
      Dimension X(*),Y(*),Sig(*)
      Dimension Cof(20),Delt(20),IPotF(20),YFit(500)
      Write(IOut,'(/,'' *** Least Squares Polynomial Fitting ***'')')
      write(IOut,'(I3,'' Points fitted with'',I2,'' Parameters'')')
     $  NPts,NPrm
C
C Partitioning of V
C
      IAF=1
C
      IB=IAf+NPts
C
      IA=IB+NPts
C
      IAtA=IA+NPts*NPrm
C
      IAtAm1=IAtA+NPrm*NPrm
C
      ID=IAtAm1+NPrm*NPrm
C
      IBet=ID+NPrm*NPrm
C
      IPrint=0
C Subtract fixed parameters
      Do 10 I=1,NPts
       YFit(I)=Y(I)-YPoly(X(I),NPrm,IPot,Cof0)
   10 Continue
      NPrmF=0
      Do 20 I=1,NPrm
       If(Cof0(I).ne.0.0d0) go to 20
       NPrmF=NPrmF+1
       IPotF(NPrmF)=IPot(I)
   20 Continue
      Call LnLstq(IOut,IPrint,NPts,NPrmF,ChiSq,X,YFit,Sig,IpotF,V(IAf),
     $  V(IB),V(IA),V(IAtA),V(IAtAm1),V(ID),V(IBet),Cof,Delt,MxScr,IV)
      Write(IOut,'(/,'' Chi Square'',D12.5)') ChiSq
      Write(IOut,'(/,'' Polynomial Coefficients'')')
      IFree=0
      Do 50 I=1,NPrm
       If(Cof0(I).eq.0.0d0) then
        IFree=IFree+1
        Cof0(I)=Cof(IFree)
        Delt0(I)=Delt(IFree)
        Write(IOut,'(D12.5,'' X**'',I1,''('',D10.5,'')'')') Cof0(I),
     $   IPot(I),Delt0(I)
       Else
        Write(IOut,'(D12.5,'' X**'',I1,'' Fixed'')') Cof0(I),IPot(I)
       EndIf
   50 Continue 
      End   
*Deck VCI1D
      Subroutine VCI1D(In,IOut,IPrint,NBas,NEig,Alpha,F,G,MxScr,V,IScr)
      Implicit Real*8 (A-H,O-Z)
C
C Driver of VCI program
C
      Dimension F(*),G(*),V(*),IScr(*)
      If(Alpha.lt.thresh) Alpha=1.0D0
      Write(IOut,'(/,'' *** VCI one-dimensional program ***'')')
      Write(IOUt,'(I5,'' Basis Functions for'',I3,'' Eigenstates'')')
     $   NBas,NEig
      Write(IOut,'(6X,''Alpha ='',F10.2)') Alpha
      write(IOut,'(/,'' F Coefficients (cm-1):'',6F10.2)')
     $  (F(I),I=1,6)
      write(IOut,'(  '' G Coefficients (cm-1):'',6F10.2)')
     $  (G(I),I=1,6)
      Call Var1D(IOut,IPrint,NBas,NEig,Alpha,F,G,V,IScr)
      Return
      End
*Deck DVR1D
      Subroutine DVR1D(In,IOut,IPrint,IMode,NPts,NBas,NEig,XX,YY,
     $  MxScr,V,IScr)
      Implicit Real*8 (A-H,O-Z)
C
C Driver of DVR program
C
      Data EMin/1.0D5/,EMax/-1.0d5/
      Dimension XX(*),YY(*),V(*)
      Dimension IScr(*)
      Write(IOut,'(/,'' *** DVR one-dimensional program ***'')')
      Write(IOut,'('' IMode ='',I2,'' with'',I5,'' Points,'',I5,
     $  '' Basis Functions and'',I5,'' States'')') IMode,NPTS,NPts,NEig
      Write(IOut,'('' Coordinate'',4X,'' Pot. Energy'')')
      EMin=1.0d10
      EMax=-1.0D10
      Do 10 IPts=1,NPts
       If(YY(IPts).lt.EMin) then
        EMin=YY(IPts)
       ElseIf(YY(IPts).gt.EMax) then
        EMax=YY(IPts)
       EndIf  
       Write(IOut,'(D12.5,3X,D12.5)') XX(IPts),YY(IPts)
   10 Continue
      XMin=XX(1)
      XMax=XX(NPts)
      DX=(XMax-XMin)/DFloat(NPts+1)
      EMin=0.0D0
      EMax=EMax-EMin
      Call QVar1D(In,IOut,IPrint,IMode,NPts,NBas,NEig,XMin,XMax,DX,
     $   EMin,EMax,XX,YY,MxScr,V,IScr)    
      Return
      End
*Deck Therm1
      Subroutine Therm1(In,IOut,IPrint,PhyCon,Multip,NAtoms,Linear,
     $ TotWt,ZPE,RTemp,DBVib,Frq)
      Implicit Real*8 (A-H,O-Z)
      Logical Linear 
      Save TRef, PRef
      Data Zero/0.0d0/,One/1.0d0/,Two/2.0d0/,Three/3.0d0/
      Data TRef/2.9815d2/, PRef/1.0d0/
      Dimension PhyCon(*),RTemp(3),DBVib(3),Frq(*)
      P=PRef
      T=TRef
      pi=4.0d0*ATan(1.0d0)
      ToRad=Pi/1.80D+2
      NAt3 = 3*NAtoms
      NVib = NAt3-6
      If(Linear) NVib=NVib+1
C compute electronic partition function
      Boltz = PhyCon(10)
      Avog  = PhyCon(5)
      Gas   = Avog * Boltz
      Degen = Float(Multip)
      QElec = Degen
      SElec = Log(Degen) * Gas
      EElec = Zero
      CElec = Zero
C compute translational partition function
      Call ThrTra(PhyCon,TotWt,P,T,QTran,STran,ETran,CTran)
C compute rotational partition function
      if(NAtoms.eq.1) then
       QRot = one
       ERot = zero
       CRot = zero
       SRot = zero
      Else
       If(Linear) then
        QRot = T/RTemp(1)
        CRot = Gas
        ERot = Gas*T
        SRot = Gas*(Log(QRot)+one)
       Else
        QRot = T/RTemp(1)+T/RTemp(2)+T/RTemp(3)
        CRot = three*gas/two
        ERot = CRot*T
        SRot = GAS*(Log(QROT)+three/two)
       EndIf
      EndIf
      QVib=zero
      QZVib=zero 
      QTot=QTran*QRot*QVib*QElec
      QZTot=QTran*QRot*QZVib*QElec
      write(IOut,'(/,5X,'' Partition Functions at '',F7.2,'' K'')') T
      Write(IOut,'(6X,32(''-''),2X)')
      write(IOut,'(5X,'' Translation         '',D12.5)') QTran
      write(IOut,'(5X,'' Rotation            '',D12.5)') QRot
      write(IOut,'(5X,'' Vibration (Bottom)  '',D12.5)') QVib
      write(IOut,'(5X,'' Vibration (from ZPE)'',D12.5)') QZVib
      write(IOut,'(5X,'' Electronic          '',D12.5)') QElec
      Write(IOut,'(6X,32(''-''),2X)')
      write(IOut,'(5X,'' Total (from Bottom) '',D12.5)') QTot
      write(IOut,'(5X,'' Total (from ZPE)    '',D12.5)') QZTot
      return
      end
*Deck VMSRot
      Subroutine VMSRot(In,IOut,IPrint,NAtoms,NSpec,RdB0,RdB0Er,RdVib,
     $  RdDvEr,RdEle,RdDeEr,IAn,AtMass0,C,AtMass,B0,ErrB0,DBVib,ErrBv,
     $  DBEle,ErrBel,Scr,IEnd)
      Implicit Real*8 (A-H,O-Z)
      Character StrInp*80,Test*5,Plane*2
      Logical RdB0,RdVib,RdEle,RdB0Er,RdDvEr,RdDeEr
      Dimension IAn(*),AtMass0(*),AtMass(*),C(3,*)
      Dimension B0(3,*),DBVib(3,*),DBEle(3,*)
      Dimension ErrB0(3,*),ErrBv(3,*),ErrBel(3,*)
      Dimension Scr(*)
C Local
      Dimension ISot(1000),IStart(30)
      Dimension WT(3)
C for Kraitchman
      Real*8 CISub(3),ErrSub(3)
      Real*8 CKRat(3,100),EKRat(3,100)
      Integer IKRAt(100)
      Logical IsPlan
C Clean initial values
      NKR=0
      Call IClear(3*NSpec,IKRAt)
      Call AClear(3*NSpec,B0)
      Call AClear(3*NSpec,ErrB0) 
      Call AClear(3*NSpec,DBVib)
      Call AClear(3*NSpec,ErrBv) 
      Call AClear(3*NSpec,DBEle)
      Call AClear(3*NSpec,ErrBEl) 
      Call AClear(3*NSpec,CKRAt)
      Call AClear(3*NSpec,EKRat)
C
C Check Planarity
C
      Tol=1.0D-6
      Call FndPl(NAtoms,C,Tol,Plane,IsPlan)
      If(IsPlan) then
       write(IOut,'(/,'' 2D Structure in Plane '',A2,'' with'',I3,
     &   '' Atoms'',/)') Plane,NAtoms
      Else
       write(IOut,'(/,'' 3D Structure with'',I3,'' Atoms'',/)') NAtoms
      EndIf
      Do IAt=1,NAtoms
       write(IOut,'(I3,3F12.5)') IAn(IAt),(C(I,IAt),I=1,3) 
      EndDo
C
      write(IOut,'(/,'' The Following Data are Given for'',I3,
     &   '' Species'')') NSpec
      If(RdVib.and.RdEle) then
       If(RdB0Er) then
         write(IOut,'( '' B0 (with errors), Vibrational and'',
     &     '' Electronic Corrections'')')
       Else
         write(IOut,'( '' B0, Vibrational and'',
     &     '' Electronic Corrections'')')
       EndIf
      ElseIf(RdVib) then
       If(RdB0Er) then
         write(IOut,'( '' B0 (with errors) and'',
     &     '' Vibrational Corrections'')')
       Else
         write(IOut,'( '' B0 and Vibrational Corrections'')')
       EndIf
      Else
       If(RdB0Er) then
         write(IOut,'( '' B0 (with errors)'')')
       Else
         write(IOut,'( '' B0 only'')')
       EndIf
      EndIf
C
      Do 10 ISpec=1,NSpec
       read(In,'(A80)') StrInp
       If(Ispec.eq.1.and.IPrint.gt.0) then
        write(IOut,'('' '')')
       ElseIf(IPrint.gt.0) then
        write(IOut,'('' For Isotopologue'',I3)') ISpec
       EndIf
       call IClear(NAtoms,ISot)
       call AClear(NAtoms,AtMass)
       Call IClear(30,IStart)
       call LlinCl(StrInp)
       Read(In,'(A80)') StrInp
       call SubStr(StrInp,20,IStart,NValue)
       If(StrInp(IStart(1):IStart(1)).eq.'0') goto 25 
       do 20 jj=1,NValue,2        
        IStr1=IStart(jj)
        IEnd1=IStart(jj+1)-1
        IStr2=IEnd1+1
        If(jj.eq.NValue-1) then
         IEnd2=IStr2+3 
        Else        
         IEnd2=IStart(jj+2)-1
        EndIf 
        call ChrNum(StrInp,IStr1,IEnd1,Numb1)
        call ChrNum(StrInp,IStr2,IEnd2,Numb2)
        ISot(Numb1)=Numb2
   20  continue 
   25  continue     
       TotWt=0.0
       Do 30 IAt=1,NAtoms
        RMass = 0.0d0
        ISpin = -1
        QMom1 = 0.0d0
        GFac1 = 0.0d0
        IAnI  = IAn(IAt)
        MNI   = ISot(IAt)
        Call FilMag(MNI,IAnI,JUse,MNO,RMass,ISpin,QMom1,GFac1)
        AtMass(IAt) = RMass
        If(IPrint.gt.0) Write(IOut,'(F12.5)') AtMass(IAt)
        TotWt=TotWt+RMass
   30  Continue
       If(ISpec.eq.1) then
C Salva le masse del parent
         Do 31 IAt=1,NAtoms
           AtMass0(IAt) = AtMass(IAt)
   31    Continue
       EndIf
       If(RdB0) read(In,*)   (B0(ii,ISpec),ii=1,3)
       If(RdB0Er) read(in,*) (ErrB0(ii,ISpec),ii=1,3)
       If(RdVib) read(In,*)  (DBVib(ii,ISpec),ii=1,3)
       If(RdDvEr) read(In,*) (ErrBv(ii,ISpec),ii=1,3)
       If(RdEle) read(In,*)  (DBEle(ii,ISpec),ii=1,3)
       If(RdDeEr) read(in,*) (ErrBel(ii,ISpec),ii=1,3)
C
C Substitution Structure
C For Planar molecules the default correction is IPlCor=2, i.e. the
C only axis j is corrected for planarity defect (see Gordy)
       IPlCor=1
       IPrin1=1
C Kraitchman per specie monosostituita
C NValue=2 vuol dire che nella riga di sostituzione c'era "iat newmass"
C (quindi un solo sito è cambiato)
C
       If (NValue.eq.2 .and. ISpec.gt.1) Then
        Call KSubSt(IOut,IPrin1,NAtoms,ISpec,AtMass0,AtMass,B0,DBVib,
     &       DBEle,ErrB0,ErrBv,ErrBel,Plane,IPlCor,ISub,CISub,ErrSub)
        NKR = NKR + 1
        IKRAt(NKR) = ISUB
        Do INew=1,3
         CKRat(INew,NKR)= CISub(INew)
         EKRAt(INew,NKR)= ErrSub(INew)
        EndDo
       EndIf
   10 Continue
C B0
      write(IOut,'(/,'' B0'')')
      do 40 ISpec=1,NSpec
       write(IOut,'(3F12.5)') (B0(ii,ISpec),ii=1,3)
   40 Continue
C B0 Errors
      if(RdB0Er) then 
       write(IOut,'(/,'' B0 errors'')')
       do 45 ISpec=1,NSpec
        write(IOut,'(3F12.5)') (ErrB0(ii,ISpec),ii=1,3)
   45  Continue
      EndIf
C DBVib
      if(RdVib) then
       write(IOut,'(/,'' DBvib'')')
       do 50 ISpec=1,NSpec
        write(IOut,'(3F12.5)') (DBVib(ii,ISpec),ii=1,3)
  50   Continue
      EndIf
C DBEle
      if(RdEle) then
       write(IOut,'(/,'' DBEle'')')
       do 60 ISpec=1,NSpec
        write(IOut,'(3F12.5)') (DBEle(ii,ISpec),ii=1,3)
  60   Continue
      EndIf
C Substitution Coordinates
      Write(IOut,'(/,'' Substitution Coordinates'')')
      Write(IOut,'('' Atom'',15X,''A'',25X,''B'',25X,''C'')')
      Do 80 I=1,NKR
       Write(IOut,'(I3,3(4X,F10.5,''('',F10.5,'')'')  )') IKRAt(I),
     &   (CKRAT(J,I),EKRAT(J,I),J=1,3)
   80 Continue
      Return
      End
*Deck UpCase
      SUBROUTINE UpCase(Str)
C===================================================================
C  Converte tutti i caratteri alfabetici di Str in maiuscolo.
C  Funziona in Fortran 77 puro usando ICHAR/CHAR.
C  La lunghezza della stringa è determinata con LEN(Str).
C===================================================================
        CHARACTER*(*) Str
        INTEGER i, ia, L
        L = LEN(Str)
        DO 10 i = 1, L
           ia = ICHAR(Str(i:i))
           IF (ia .GE. ICHAR('a') .AND. ia .LE. ICHAR('z')) THEN
C            Converte da minuscolo a maiuscolo sottraendo 32
              Str(i:i) = CHAR(ia - 32)
           ENDIF
   10   CONTINUE
        RETURN
      END
*Deck ParseFlags
      SUBROUTINE ParseFlags(Line, NSPEC, RdB0Er, RdVib, RdEle)
C===================================================================
C  Interpreta una riga contenente:
C    - un numero intero (NSPEC)
C    - zero o più delle parole chiave: 'errors', 'deltavib', 'deltael'
C  La ricerca delle parole chiave è insensibile al maiuscolo/minuscolo.
C===================================================================
        IMPLICIT REAL*8 (A-H,O-Z)
        CHARACTER*(*) Line
        INTEGER NSPEC
        LOGICAL RdB0Er, RdVib, RdEle
        CHARACTER*256 Tmp
C  Inizializza i valori di ritorno
        NSPEC  = 0
        RdB0Er = .FALSE.
        RdVib  = .FALSE.
        RdEle  = .FALSE.
C  Copia la stringa in una variabile temporanea (dimensione adeguata)
        Tmp = Line
C  Converte Tmp in maiuscolo
        CALL UpCase(Tmp)
C  Legge il numero (se presente) dal file interno. Se fallisce, NSPEC resta 0.
        READ(Tmp, *, ERR=100) NSPEC
  100   CONTINUE
C  Controlla le parole chiave in maiuscolo
        IF (INDEX(Tmp, 'ERRORS')   .GT. 0) RdB0Er = .TRUE.
        IF (INDEX(Tmp, 'DELTAVIB') .GT. 0) RdVib  = .TRUE.
        IF (INDEX(Tmp, 'DELTAEL')  .GT. 0) RdEle  = .TRUE.
        RETURN
      END
*Deck FndPl
      SUBROUTINE FndPl(NAtoms,C,tol,Plane,IsPlanar)
C===================================================================
C  Determina automaticamente il piano in cui giace una molecola
C  orientata lungo gli assi principali e centrata nel baricentro.
C  Se tutte le coordinate lungo un asse sono ~0 → molecola planare.
C
C  Input:
C    NAtoms  : numero di atomi
C    C(3,*)  : coordinate [Å] nel sistema degli assi principali
C    tol     : tolleranza per considerare una coordinata come zero
C  Output:
C    Plane   : 'ab', 'ac', 'bc' se planare; stringa vuota altrimenti
C    IsPlanar: .TRUE. se la molecola è planare, .FALSE. altrimenti
C===================================================================
        IMPLICIT REAL*8 (A-H,O-Z)
        INTEGER NAtoms
        REAL*8 C(3,*)
        REAL*8 tol
        CHARACTER*(*) Plane
        LOGICAL IsPlanar
        LOGICAL allZero
C  Inizializza i valori di ritorno
        Plane    = '  '
        IsPlanar = .FALSE.
C  Controlla asse X → piano bc
        allZero = .TRUE.
        DO I = 1, NAtoms
           IF (ABS(C(1,I)) .GT. tol) allZero = .FALSE.
        ENDDO
        IF (allZero) THEN
           Plane    = 'bc'
           IsPlanar = .TRUE.
           RETURN
        ENDIF
C  Controlla asse Y → piano ac
        allZero = .TRUE.
        DO I = 1, NAtoms
           IF (ABS(C(2,I)) .GT. tol) allZero = .FALSE.
        ENDDO
        IF (allZero) THEN
           Plane    = 'ac'
           IsPlanar = .TRUE.
           RETURN
        ENDIF
C  Controlla asse Z → piano ab
        allZero = .TRUE.
        DO I = 1, NAtoms
           IF (ABS(C(3,I)) .GT. tol) allZero = .FALSE.
        ENDDO
        IF (allZero) THEN
           Plane    = 'ab'
           IsPlanar = .TRUE.
           RETURN
        ENDIF
C  Se arriva qui, non è planare
        Plane    = '  '
        IsPlanar = .FALSE.
        RETURN
      END
*Deck KSubSt
      SUBROUTINE KSubSt(IOut,IPrint,NAtoms,ISpec,
     &   AtMass0,AtMass,B0,DBVib,DBEle,
     &   dB0,dDBVib,dDBEle,Plane,Method,
     &   ISub,CISub,ErrSub)
C===================================================================
C Calcola le coordinate di Kraitchman per un isotopologo planare
C includendo:
C  - correzione del difetto inerziale per il piano specificato,
C  - propagazione degli errori delle costanti rotazionali,
C  - aggiunta dell'errore empirico secondo la regola di Costain.
C
C Input:
C   AtMass0, AtMass : masse (amu) della specie di riferimento e sostituita
C   B0, DBVib, DBEle: matrici (3,Nspec) delle costanti rotazionali (MHz)
C   dB0, dDBVib, dDBEle: errori (MHz) associati a B0, DBVib, DBEle
C   Plane : stringa 'ab','ac' o 'bc' che indica il piano della molecola
C   Method: type of planar correction (0=none, 1=only j axis, 2=symmetric)
C Output:
C   CISub(1:3) : |a|,|b|,|c| (Å)
C   ErrSub(1:3): errori totali (Å) sulle coordinate
C===================================================================
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*(*) Plane
      INTEGER IOut,NAtoms,ISpec
      REAL*8 AtMass0(*),AtMass(*)
      REAL*8 B0(3,*),DBVib(3,*),DBEle(3,*)
      REAL*8 dB0(3,*),dDBVib(3,*),dDBEle(3,*)
      REAL*8 CISub(*),ErrSub(*)
C-------------------------------------------------------------------
      REAL*8 I0(3),II(3),DI(3),dI0(3),dII(3),dDI(3)
      REAL*8 DIcorr(3),dDIcorr(3)
      REAL*8 A2,B2,C2,dA2,dB2,dC2
      REAL*8 A,B,C,dA,dB,dC
      REAL*8 M0,M1,DeltaM
      REAL*8 B0E,BIE,dB0E,dBIE
      REAL*8 H,CONV,PI,FACT
      REAL*8 KCONST,delta_def,dDelta_def
      INTEGER I,IA
C-------------------------------------------------------------------
C  Physical Constants (B in MHz → I in amu·Å²)
      H     = 6.62607015D-34
      CONV  = 1.66053906660D-47
      PI    = 4.0D0*DATAN(1.0D0)
      FACT  = H/(8.0D0*PI*PI*CONV*1.0D6)
C  Costant for Costain rule (Å)
      KCONST = 0.0015D0
C-------------------------------------------------------------------
C  1. Find Substituted Atom
      ISub = 0
      NSub = 0
      DO IA=1,NAtoms
         IF (DABS( AtMass(IA) - AtMass0(IA) ).GT.1.0D-8) THEN
            ISub = IA
            NSub = NSub +1 
         ENDIF
      ENDDO
      If(NSub.eq.0.or.NSub.gt.1) then
         WRITE(IOut,'(/,I3,A)') NSub,
     &     ' Isotopic Substitutions: Rs not computed'
         RETURN
      ENDIF
C-------------------------------------------------------------------
C  2. Local Mass Difference (amu)
      M0 = AtMass0(ISub)
      M1 = AtMass(ISub)
      DeltaM = M1 - M0
      IF (DABS(DeltaM).LT.1.0D-12) THEN
         WRITE(IOut,'(/,A,I3,A)') 
     &     ' For Substituted Atom',ISub,' DeltaM ~ 0: Rs not computed'
        RETURN
      ENDIF
C-------------------------------------------------------------------
C  3. Inertia Moments and their Errors
      DO I=1,3
         B0E  = B0(I,1)     + DBVib(I,1)     + DBEle(I,1)
         BIE  = B0(I,ISpec) + DBVib(I,ISpec) + DBEle(I,ISpec)
         dB0E = dB0(I,1)    + dDBVib(I,1)    + dDBEle(I,1)
         dBIE = dB0(I,ISpec)+ dDBVib(I,ISpec)+ dDBEle(I,ISpec)
         IF (B0E.LE.0.0D0) B0E = 1.0D-12
         IF (BIE.LE.0.0D0) BIE = 1.0D-12
         I0(I)   = FACT / B0E
         II(I)   = FACT / BIE
C        Error Propagation: δI = I * δB / B
         dI0(I)  = I0(I) * (dB0E / B0E)
         dII(I)  = II(I) * (dBIE / BIE)
      ENDDO
C  ΔI differences and their errors
      DO I=1,3
         DI(I)  = II(I) - I0(I)
         dDI(I) = DSQRT( dII(I)**2 + dI0(I)**2 )
      ENDDO
C-------------------------------------------------------------------
C  4. Correzione del difetto inerziale e propagazione degli errori
      Call FixDDI(DI, dDI, Plane, Method, DIcorr, dDIcorr)
      IF (Plane .EQ. 'ab' .OR. Plane .EQ. 'AB') THEN
         delta_def = DIcorr(1) + DIcorr(2) - DIcorr(3)
      ELSEIF (Plane .EQ. 'ac' .OR. Plane .EQ. 'AC') THEN
         delta_def = DIcorr(1) + DIcorr(3) - DIcorr(2)
      ELSEIF (Plane .EQ. 'bc' .OR. Plane .EQ. 'BC') THEN
         delta_def = DIcorr(2) + DIcorr(3) - DIcorr(1)
      ENDIF
C-------------------------------------------------------------------
C  5. Coordinate (|a|^2, |b|^2, |c|^2) a seconda del piano
      IF (Plane .EQ. 'ab' .OR. Plane .EQ. 'AB') THEN
         A2  = DIcorr(2) / ABS(DeltaM)   ! |a|^2 = ΔI_b_corr / |Δm|
         B2  = DIcorr(1) / ABS(DeltaM)   ! |b|^2 = ΔI_a_corr / |Δm|
         C2  = 0.0D0
         dA2 = dDIcorr(2) / ABS(DeltaM)
         dB2 = dDIcorr(1) / ABS(DeltaM)
         dC2 = 0.0D0
      ELSEIF (Plane .EQ. 'ac' .OR. Plane .EQ. 'AC') THEN
         A2  = DIcorr(3) / ABS(DeltaM)   ! |a|^2 = ΔI_c_corr / |Δm|
         C2  = DIcorr(1) / ABS(DeltaM)   ! |c|^2 = ΔI_a_corr / |Δm|
         B2  = 0.0D0
         dA2 = dDIcorr(3) / ABS(DeltaM)
         dC2 = dDIcorr(1) / ABS(DeltaM)
         dB2 = 0.0D0
      ELSEIF (Plane .EQ. 'bc' .OR. Plane .EQ. 'BC') THEN
         B2  = DIcorr(3) / ABS(DeltaM)   ! |b|^2 = ΔI_c_corr / |Δm|
         C2  = DIcorr(2) / ABS(DeltaM)   ! |c|^2 = ΔI_b_corr / |Δm|
         A2  = 0.0D0
         dB2 = dDIcorr(3) / ABS(DeltaM)
         dC2 = dDIcorr(2) / ABS(DeltaM)
         dA2 = 0.0D0
      ENDIF
      If(IPrint.gt.0) write(IOut,'('' Plane = '',A2)') plane
C-------------------------------------------------------------------
C  6. Converte in coordinate lineari e propaga l'errore
C     (se A2<0 la coordinata è immaginaria; usiamo il modulo)
      IF (A2.GE.0.D0) THEN
         A  = DSQRT(A2)
         dA = 0.5D0 * dA2 / MAX(A,1.0D-12)
      ELSE
         A  = DSQRT(-A2)
         dA = 0.5D0 * dA2 / MAX(A,1.0D-12)
      ENDIF
      IF (B2.GE.0.D0) THEN
         B  = DSQRT(B2)
         dB = 0.5D0 * dB2 / MAX(B,1.0D-12)
      ELSE
         B  = DSQRT(-B2)
         dB = 0.5D0 * dB2 / MAX(B,1.0D-12)
      ENDIF
      IF (C2.GE.0.D0) THEN
         C  = DSQRT(C2)
         dC = 0.5D0 * dC2 / MAX(C,1.0D-12)
      ELSE
         C  = DSQRT(-C2)
         dC = 0.5D0 * dC2 / MAX(C,1.0D-12)
      ENDIF
C-------------------------------------------------------------------
C  7. Aggiungi in quadratura l'errore di Costain: σ_tot = sqrt(δ^2 + (K/|z|)^2)
      CISub(1) = A
      CISub(2) = B
      CISub(3) = C
      IF (A .NE. 0.D0) THEN
         ErrSub(1) = DSQRT( dA**2 + (KCONST / ABS(A))**2 )
      ELSE
         ErrSub(1) = 0.D0
      ENDIF
      IF (B .NE. 0.D0) THEN
         ErrSub(2) = DSQRT( dB**2 + (KCONST / ABS(B))**2 )
      ELSE
         ErrSub(2) = 0.D0
      ENDIF
      IF (C .NE. 0.D0) THEN
         ErrSub(3) = DSQRT( dC**2 + (KCONST / ABS(C))**2 )
      ELSE
         ErrSub(3) = 0.D0
      ENDIF
C-------------------------------------------------------------------
C  8. Output facoltativo
      IF (IPrint.GT.0) THEN
         WRITE(IOut,'(/,A)') 
     &    ' --- Kraitchman planare con difetto inerziale ed errori ---'
         WRITE(IOut,'(A,I4)') ' Atomo sostituito =', ISub
         WRITE(IOut,'(A,3F14.6)')
     &    ' |a|, |b|, |c|  (Å)  =', CISub(1),CISub(2),CISub(3)
         WRITE(IOut,'(A,3F14.6)')
     &    ' σ(|a|), σ(|b|), σ(|c|) (Å) =', ErrSub(1),ErrSub(2),ErrSub(3)
         WRITE(IOut,'(A,3E16.8)')
     &    ' ΔI_a,ΔI_b,ΔI_c (amu Å^2) =', DI(1),DI(2),DI(3)
         WRITE(IOut,'(A,E16.8)')
     &    ' Difetto inerziale corretto (amu Å^2) =', delta_def
         WRITE(IOut,'(A)')
     &    ' ----------------------------------------------------------'
      ENDIF

      RETURN
      END
*Deck FixDDI
      SUBROUTINE FixDDI(DI, dDI, Plane, Method, DIcorr, dDIcorr)
C===================================================================
C  Corregge le differenze di momento d'inerzia ΔI e propaga gli errori.
C
C  Input:
C    DI(3)   : ΔI_a, ΔI_b, ΔI_c (amu·Å²)
C    dDI(3)  : errori su ΔI_a, ΔI_b, ΔI_c
C    Plane   : 'ab','ac','bc' come spiegato in FixPlanarDI
C    Method  : 0 = nessuna correzione
C               1 = correzione PROSPE (ricalcola solo l'asse j)
C               2 = correzione simmetrica (ripartisce il difetto)
C  Output:
C    DIcorr(3) : ΔI corretti
C    dDIcorr(3): errori corretti
C===================================================================
        REAL*8 DI(3), dDI(3), DIcorr(3), dDIcorr(3)
        CHARACTER*(*) Plane
        INTEGER Method
        INTEGER i, j, k
        REAL*8 delta, ddelta

C  Copia i valori originali
        DO 5 i = 1,3
           DIcorr(i)  = DI(i)
           dDIcorr(i) = dDI(i)
    5   CONTINUE

C  Determina gli indici i,j,k in base al piano
        IF (Plane .EQ. 'ab' .OR. Plane .EQ. 'AB') THEN
           i = 1
           j = 2
           k = 3
        ELSEIF (Plane .EQ. 'ac' .OR. Plane .EQ. 'AC') THEN
           i = 1
           j = 3
           k = 2
        ELSEIF (Plane .EQ. 'bc' .OR. Plane .EQ. 'BC') THEN
           i = 2
           j = 3
           k = 1
        ELSE
C         Piano non valido ⇒ nessuna correzione
           RETURN
        ENDIF

        IF (Method .EQ. 1) THEN
C         Correzione PROSPE: ΔI_j_corr = ΔI_k - ΔI_i
           DIcorr(j)  = DI(k) - DI(i)
           dDIcorr(j) = SQRT( dDI(k)**2 + dDI(i)**2 )
C         Gli altri rimangono invariati
        ELSEIF (Method .EQ. 2) THEN
C         Correzione simmetrica: ripartisce il difetto δ = ΔI_i + ΔI_j - ΔI_k
           delta  = DI(i) + DI(j) - DI(k)
           ddelta = SQRT( dDI(i)**2 + dDI(j)**2 + dDI(k)**2 )
           DIcorr(i)  = DI(i) - 0.5D0*delta
           DIcorr(j)  = DI(j) - 0.5D0*delta
           dDIcorr(i) = SQRT( dDI(i)**2 + (0.5D0*ddelta)**2 )
           dDIcorr(j) = SQRT( dDI(j)**2 + (0.5D0*ddelta)**2 )
C         ΔI_k e dΔI_k restano invariati
        ELSE
C         Method = 0: nessuna correzione
        ENDIF

        RETURN
      END
*Deck SemiEx
      Subroutine SemiEx(In,IOut,IPrint,NAtoms,NSpec,RdExp,RdVib,RdEle,
     $  IAn,AtMass,B0,DBVib,DBEle)
      Implicit Real*8 (A-H,O-Z)
      Character StrInp*80,Test*5
      Logical RdExp,RdVib,RdEle
      Dimension IAn(*),AtMass(*),B0(3,*),DBVib(3,*),DBEle(3,*)
C Local
      Dimension ISot(1000),IStart(30)
      Dimension WT(3)
C
      RdExp=.true.
      RdVib=.true.
      RdEle=.true. 
C
      Call AClear(3*NSpec,B0)
      Call AClear(3*NSpec,DBVib)
      Call AClear(3*NSpec,DBEle)
C
      Do 10 ISpec=1,NSpec
       read(In,'(A80)') StrInp
       If(Ispec.eq.1) then
        write(IOut,'('' '')')
       Else
        write(IOut,'('' For Isotopologue'',I3)') ISpec
       EndIf
       call IClear(NAtoms,ISot)
       call AClear(NAtoms,AtMass)
       Call IClear(30,IStart)
       call LlinCl(StrInp)
       Read(In,'(A80)') StrInp
       call SubStr(StrInp,20,IStart,NValue)
       If(StrInp(IStart(1):IStart(1)).eq.'0') goto 25 
       do 20 jj=1,NValue,2        
        IStr1=IStart(jj)
        IEnd1=IStart(jj+1)-1
        IStr2=IEnd1+1
        If(jj.eq.NValue-1) then
         IEnd2=IStr2+3 
        Else        
         IEnd2=IStart(jj+2)-1
        EndIf 
        call ChrNum(StrInp,IStr1,IEnd1,Numb1)
        call ChrNum(StrInp,IStr2,IEnd2,Numb2)
        ISot(Numb1)=Numb2
   20  continue 
   25  continue     
       Do 30 IAt=1,NAtoms
        RMass = 0.0d0
        ISpin = -1
        QMom1 = 0.0d0
        GFac1 = 0.0d0
        IAnI  = IAn(IAt)
        MNI   = ISot(IAt)
        Call FilMag(MNI,IAnI,JUse,MNO,RMass,ISpin,QMom1,GFac1)
        AtMass(IAt) = RMass
        Write(IOut,'(F12.5)') AtMass(IAt)
   30  Continue
       If(RdExp) read(In,*) (B0(ii,ISpec),ii=1,3)
       If(RdVib) read(In,*) (DBVib(ii,ISpec),ii=1,3)
       If(RdEle) read(In,*) (DBEle(ii,ISpec),ii=1,3)
   10 Continue
C Bexp
      write(IOut,'(/,'' Bexp'')')
      do 40 ISpec=1,NSpec
       write(IOut,'(3F12.5)') (B0(ii,ISpec),ii=1,3)
   40 Continue
C DBVib
      write(IOut,'(/,'' DBvib'')')
      do 50 ISpec=1,NSpec
       write(IOut,'(3F12.5)') (DBVib(ii,ISpec),ii=1,3)
   50 Continue
C DBEle
      write(IOut,'(/,'' DBEle'')')
      do 60 ISpec=1,NSpec
      write(IOut,'(3F12.5)') (DBEle(ii,ISpec),ii=1,3)
   60 Continue
      write(IOut,'(/,'' Weights'')')
      do 70 ISpec=1,NSpec
       Wt(1)=1.0D0
       Wt(2)=1.0D0
       Wt(3)=1.0D0
       write(IOut,'(3F12.5)') (WT(ii),ii=1,3)
   70 Continue
      write(IOut,'('' '')')
      Return
      End
*Deck ChrNum 
      Subroutine ChrNum(String,IStart,IEnd,Number)
      Implicit Real*8 (A-H,O-Z)
C
C     FORM A NUMBER FROM THE Character DIGITS In String.
C
      Character String*80
      Character N*10
      Logical Found
      Save N
      DATA N/'0123456789'/
C
      NTot=IEnd-IStart
      Number=0
      NPow=-1
      do 10 ii=IEnd,IStart,-1 
       if(String(ii:ii).eq.' ') goto 10
       found=.false.
       IBas=0
       do 20 j=1,10
        if(found) goto 20
        if(String(ii:ii).eq.N(j:j)) then
         IBas=j-1
         Npow=NPow+1
         found=.true.
        endif
   20  continue
       Number=Number+IBas*10**NPow  
   10 continue
      return
      end
*Deck SetPCS
      Subroutine SetPCS(String,ModPCS,FndPCS)
      Character String*8
      Logical FndPCS
      ModPCS=-1
      FndPCS=.False.
      If(String(1:4).eq.'PCS0') then
       ModPCS=0
      ElseIf(String(1:4).eq.'PCS1') then
       ModPCS=1
      ElseIf(String(1:5).eq.'HPCS2') then
       ModPCS=2
      ElseIf(String(1:5).eq.'DPCS3') then
       ModPCS=3
      EndIf
      If(ModPCS.ge.0) FndPCS=.true. 
      Return
      End          
*Deck SetDer
      Subroutine SetDer(String,IDeriv,FndDer)
      Character String*8
      Logical FndDer
      Integer IDeriv(4) 
      Call IClear(4,IDeriv)
      ITot=0
      FndDer=.false.
      If(String(1:3).eq.'OPT') then
       IDeriv(1)=1
       ITot=ITot+1
      EndIf
      If(String(1:4).eq.'FREQ') then
       IDeriv(2)=1
       ITot=ITot+1
      ElseIf(String(1:5).eq.'CUBIC') then
       IDERIV(3)=1
       ITot=ITot+1
      ElseIf(String(1:6).eq.'ANHARM') then
       IDeriv(4)=1
       ITot=ITot+1
      EndIf
C     If(ITot.gt.0) FndDer=.true.
      Return
      End
*Deck SetGKw
      Subroutine SetGKw(IOut,DoGNIC,SyGNIC,Loose,DoScan,ModPCS,IDeriv)
      Logical DoGNIC,SyGNIC,Loose,DoScan
      Integer IOut,ModPCS,IDeriv(4) 
      Write(IOut,'(''%Nprocshared=8'')')
      Write(IOut,'(''%Mem=32GB'')')
      Write(IOut,'(''%chk=gicforge.chk'')')
      If(DoGNIC) then
       If(SyGNIC) then
        If(Ideriv(1).eq.1.and.IDeriv(2).eq.0) then
         Write(IOut,'(''#P geom=(readallgic,gicsymm) '')',
     $    advance='no')
        Else
         Write(IOut,'(''#P geom=(readallgic,gicallsymm) '')',
     $    advance='no')
        EndIf
       Else
        Write(IOut,'(''#P geom=readallgic '')',advance='no')
       EndIf
      Else
       Write(IOut,'(''#P'')',advance='no')
      EndIf
      If(Loose) Write(IOut,'(''Symm=loose '')',advance='no')
      If(ModPCS.eq.0) then
       Write(IOut,'(''UFF '')',advance='no')
      ElseIf(ModPCS.eq.1) then
       Write(IOut,'(''HF3C '')',advance='no')
      ElseIf(ModPCS.eq.2) then
       Write(IOut,'(''B3LYP EMPIRICALDISPERSION=GD4 6-31G* '')',
     $   advance='no') 
      ElseIf(ModPCS.eq.3) then
       Write(IOut,'(''revDSDPBEP86D4 gen '')',advance='no')
      EndIf 
      If(IDeriv(1).eq.1.or.DoScan) then 
       If(ModPCS.eq.0) then
        Write(IOut,'(''OPT=nomicro '')',advance='no')
       ElseIf(ModPCS.gt.1.and..not.DoScan) then
        Write(IOut,'(''OPT=calcHFFC '')',advance='no')
       Else
        Write(IOut,'(''OPT '')',advance='no')
       EndIf
      EndIf
      If(IDeriv(2).eq.1) Write(IOut,'(''Freq=IntModes '')',advance='no')
      If(IDeriv(3).eq.1) Write(IOut,'(''Freq=Cubic '')',advance='no')
      If(IDeriv(4).eq.1) Write(IOut,'(''Freq=Anharm '')',advance='no')
      write(IOut,'(''Output=Pickett'')')
      Return
      End
*Deck PMOMG
      Subroutine PMOMG(NAtoms,AtMass,C,PMom1)
      Implicit Real*8 (A-H,O-Z)
C Computes the first derivatives of the diagonal elements of
C the moment of inertia tensor w.r.t. Cartesian coordinates
C I/O
C NAtoms            = Number of Atoms
C C(3,NAtoms)       = Cartesian Coordinates
C AtMass(NAtoms)    = Atomic Masses
C PMom1(3,3*NAtoms) = Cartesian first derivatives of the diagonal
      Dimension AtMass(*),C(3,*),PMom1(3,*)
      Zero=0.0d0
      Two=2.0d0
      Do 10 IAt=1,NAtoms
       ix = 3*(i-1) + 1
       iy = 3*(i-1) + 2
       iz = 3*(i-1) + 3
C XX Component
       PMom1(1,ix) = zero
       PMom1(1,iy) = two*atmass(i)*c(2,i) 
       PMom1(1,iz) = two*atmass(i)*c(3,i)
C YY Component
       PMom1(2,ix) = two*atmass(i)*c(1,i) 
       PMom1(2,iy) = zero
       PMom1(2,iz) = two*atmass(i)*c(3,i)
C ZZ Component
       PMom1(3,ix) = two*atmass(i)*c(1,i)
       PMom1(3,iy) = two*atmass(i)*c(2,i)
       PMom1(3,iz) = zero
   10 Continue
      Return
      End   
*Deck PMOMH
      Subroutine PMOMH(NAtoms,AtMass,C,PMom2)
      Implicit Real*8 (A-H,O-Z)
C Computes the second derivatives of the diagonal elements
C of the moment of inertia tensor w.r.t. Cartesian coordinates 
C I/O
C NAtoms             = Number of Atoms (Nat3 = 3*NAtoms)
C C(3,NAtoms)        = Cartesian Coordinates
C AtMass(NAtoms)     = Atomic Masses
C PMom2(3,NAt3,NAt3) = Second derivatives of the diagonal elements of
C                      the Inertia Moment Tensor w.r.t. Cartesian Coordinates
      Dimension AtMass(*),C(3,*),PMom2(3,3*Natoms,*)
      Zero=0.0d0
      Two=2.0d0
      Do 10 IAt=1,NAtoms
       ix = 3*(i-1) + 1
       iy = 3*(i-1) + 2
       iz = 3*(i-1) + 3
C XX Component
       PMom2(1,ix,ix) = zero
       PMom2(1,iy,iy) = two*atmass(i) 
       PMom2(1,iz,iz) = two*atmass(i)
C YY Component 
       PMom2(2,ix,ix) = two*atmass(i)  
       PMom2(2,iy,iy) = zero
       PMom2(2,iz,iz) = two*atmass(i)  
C ZZ Component
       PMom2(3,ix,ix) = two*atmass(i)
       PMom2(3,iy,iy) = two*atmass(i)
       PMom2(3,iz,iz) = zero
   10 Continue
      Return
      End
* Deck SMIF77
      Subroutine SMIF77(IOut, IPrint, MxAt, MxBnd, NAtoms, IAn, IArom,
     $  InCyc, NBond, IBond, IZ, BndOrd, BL, Alpha, Beta, SMILES)
      IMPLICIT Integer (A-Z)
C     Determine Molecular Topology from a SMILES
C I/O
      Dimension IAn(*), IArom(*), InCyc(*), NBond(*), IBond(MxBnd,*)
      Dimension IZ(4,*),IQAt(MxAt)
      Dimension IBranch(MxAt),IRing(MxAt),IRStrt(10),NAtCyc(10)
      Real*8 BL(*), Alpha(*), Beta(*), BndOrd(MxBnd,*)
      Character*100 SMILES
C Local
      Dimension NPi(MxAt), IStack(MxAt)
      Character*1 ch, nextch
      Character*2 AtSymb, AtomTy(MxAt)
C Initialize Arrays
      LenSMI = LEN_TRIM(SMILES)
      Call IClear(MxAt, NBond)
      Call IClear(MxAt, NPi)
      Call IClear(MxAt, IArom)
      Call IClear(4*MxAt, IZ)
      Call IClear(10, IRStrt)
      Call IClear(10, NatCyc)
      Call IClear(MxAt, IStack)
      Call IClear(MxAt, IBranch)
      Call IClear(MxAt, IRing)
      Call IClear(MxAt, IQAt)
      Call AClear(MxAt, BL)
      Call AClear(MxAt, Alpha)
      Call Aclear(MxAt, Beta)
      IRId = 0
      do IAt = 2, MxAt
         BL(IAt) = 1.0d0
         If (IAt.gt.2) Alpha(IAt) = 109.47122063D0
         If (IAt.gt.3) Beta(IAt) = 180.0D0
      End Do
      NAtoms = 0
      LAstAt = 0
      IBType = 1
      ITop = 0
      IPos = 1
   10 If (IPos.gt.LenSMI) goto 100
      ch = SMILES(IPos:IPos)
C Begin new code
      if (ch.eq.'[') then
         InBracket = 1
         HCount = 0
         LCharg = 0
         AtSymb = '  '
      elseif (ch.eq.']') then
         InBracket = 0
         NAtoms = NAtoms + 1
         AtomTy(NAtoms) = AtSymb
         If (HCount.gt.0) NPi(NAtoms) = HCount
         If (LCharg.ne.0) Call SetCharge(NAtoms,LCharg,IQAt)
         If (LastAt.ne.0) then
            call SetBnd(MxBnd,NAtoms,LastAt,IBType,NBond, IBond,BndOrd)
            call SetIZ(MxBnd, NAtoms, LastAt, NBond, IBond, IZ)
         End If
         LastAt = NAtoms
         IBType = 1
      elseif (InBracket.eq.1) then
         if (ch.eq.'H') then
            nextch = SMILES(IPos+1:IPos+1)
            if (nextch.ge.'0'.and.nextch.le.'9') then
               HCount = IChar(nextch) - ichar('0')
               IPos = IPos + 1
            else
               HCount = 1
            endif
         elseif (ch.eq.'+' .or. ch.eq.'-') then
            LCharg = IChar(ch)
            nextch = SMILES(IPos+1:IPos+1)
            if (nextch.ge.'0'.and.nextch.le.'9') then
               LCharg = LCharg * (IChar(nextch) - ichar('0'))
               IPos = IPos + 1
            endif
         else
            AtSymb(1:1) = ch
            nextch = SMILES(IPos+1:IPos+1)
            if (nextch.ge.'a'.and.nextch.le.'z') then
               AtSymb(2:2) = nextch
               IPos = IPos + 1
            endif
         endif
C End New Code
      elseif (ch.eq.'(') then
         ITop = ITop + 1
         IStack(ITop) = LastAt
         IBranch(NAtoms) = ITop
      elseif (ch.eq.')') then
         LastAt = IStack(ITop)
         ITop = ITop - 1
      elseif (ch.eq.'-') then
         IBType = 1
      elseif (ch.eq.'=') then
         IBType = 2
      elseif (ch.eq.'#') then
         IBType = 3
      elseif (ch.ge.'a'.and.ch.le.'z') then
         AtSymb(1:1) = ch
         AtSymb(2:2) = ' '
         NAtoms = NAtoms + 1
         AtomTy(NAtoms) = AtSymb
         IArom(NAtoms) = 1
         If (LastAt.ne.0) then
            call SetBnd(MxBnd,NAtoms,LastAt,IBType,NBond, IBond,BndOrd)
            call SetIZ(MxBnd, NAtoms, LastAt, NBond, IBond, IZ)
         End If
         LastAt = NAtoms
         IBType = 1
      elseif (ch.ge.'0'.and.ch.le.'9') then
         IDig = IChar(ch) - ichar('0')
         if (IRStrt(IDig).eq.0) then
            IRStrt(IDig) = LastAt
            IRId = IRId + 1
            IRing(LastAt) = IRId
         else
            IRidPrev = IRStrt(IDig)
            call SetBnd(MxBnd,IRidPrev,LastAt,IBType,NBond,IBond,BndOrd)
            NPi(LastAt) = NPi(LastAt) + IBType - 1
            NPi(IRidPrev) = NPi(IRidPrev) + IBType - 1
            IRing(LastAt) = IRing(IRidPrev)
            IRStrt(IDig) = 0
            ! Propagate Ring ID to all connected atoms
            call PropagateRingID(MxAt,MxBnd,NBond,IBond,IRing,
     $       IRing(LastAt),IRidPrev,LastAt)
         endif
         IBType = 1
      elseif (ch.ge.'A'.and.ch.le.'Z') then
         AtSymb(1:1) = ch
         If (IPos.lt.LenSMI) then
            nextch = SMILES(IPos+1:IPos+1)
            if (nextch.ge.'a'.and.nextch.le.'z') then
               AtSymb(2:2) = nextch
               IPos = IPos + 1
            else
               AtSymb(2:2) = ' '
            endif
         else
            AtSymb(2:2) = ' '
         endif
         NAtoms = NAtoms + 1
         AtomTy(NAtoms) = AtSymb
         If (LastAt.ne.0) then
            call SetBnd(MxBnd,NAtoms,LastAt,IBType,NBond,IBond,BndOrd)
            call SetIZ(MxBnd, NAtoms, LastAt, NBond, IBond, IZ)
            NPi(NAtoms) = NPi(NAtoms) + IBType - 1
            NPi(LastAt) = NPi(LastAt) + IBType - 1
         End If
         LastAt = NAtoms
         IBType = 1
      End If
      IPos = IPos + 1
      goto 10
  100 Continue
C Compute number of atoms in rings and Atomic Numbers. Correct Bond Orders
      NRing=0
      Do IAt = 1, NAtoms
         IRAt=IRing(IAt)
         If(IRAt.gt.NRing) NRing=IRAt
         If(IRAt.ne.0) NAtCyc(IRAt)=NatCyc(IRAt)+1
         Call FilIAn(AtomTy(IAt), IAn(IAt))
         If (IArom(IAt).ne.0) then
            NBI = NBond(IAt)
            Do J = 1, NBI
               JAt = IBond(J, IAt)
               If (IArom(JAt).ne.0) BndOrd(J, IAt) = 1.5D0
            end do
         endif
      end do
      If(NRing.gt.0) then
       write(IOut,'('' The Molecule contains'',I2,'' Rings'')') NRing
       do ir=1,NRing
        write(IOut,'('' Ring'',I2,'' has'',I2,'' Atoms'')')IR,NAtCyc(IR)
       end do
      endif
      NHeavy = NAtoms
      ! Add Hydrogens
      Call AddHyd(MxAt,MxBnd,NHeavy,NAtoms,IAn,NPi,NBond,IBond,IZ,
     $  BndOrd,Alpha,AtomTy)
      ! Set Bond Lengths and Valence Angles
      call SetBLA(MxBnd,NAtoms,IAn,IZ,NBond,IBond,BndOrd,BL,Alpha)
      ! Set Dihedral Angles
      call SetDih(NAtoms,IZ,NBond,NAtCyc,IBranch,IRing,Alpha,Beta)
      ! Print Results
      call SMIBnd(IOut,MxBnd,NAtoms,IAn,IArom,IBranch,IRing,
     $  NBond, IBond, BndOrd)
      Return
      End
*Deck SetBnd
      Subroutine SetBnd(MxBnd,IAt,JAt,IBType,NBond,IBond,BndOrd)
      Implicit Integer (A-Z)
C Set Parameters for IAt and JAt (IZ only for IAt)
      Dimension NBond(*),IBond(MxBnd,*)
      Real*8 BndOrd(MxBnd,*)
C Set Bond patterns
      NBond(IAt)=NBond(IAt)+1
      IBond(NBond(IAt),IAt)=JAt
      BndOrd(NBond(IAt),IAt)=Float(IBType)
      NBond(JAt)=NBond(JAt)+1
      IBond(NBond(JAt),JAt)=IAt
      BndOrd(NBond(JAt),JAt)=Float(IBType)  
      Return
      End
*Deck SetIZ
      Subroutine SetIZ(MxBnd,IAt,JAt,NBond,IBond,IZ)
      Dimension NBond(*),IBond(MxBnd,*),IZ(4,*)
      IZ(1,IAt)=JAt
      If(IAt.eq.2) return
      KAt=IZ(1,JAt)
      If(KAt.eq.0) then
       KAt=IBond(2,JAt)
       If(IBond(2,JAt).eq.IAt) KAt=IBond(1,JAt)
      EndIf
      IZ(2,IAt)=KAt
      If(IAt.eq.3) return
C fourth and other atoms
      If(KAt.eq.1) then
       If(NBond(KAt).gt.1) then
        LAt=IBond(1,KAt)
        If(LAt.eq.JAt) LAt=IBond(2,KAt)
        IZ(4,IAt)=0
       Else
        LAt=IBond(2,JAt)
        If(LAt.eq.KAt) LAt=IBond(1,JAt)
        IZ(4,IAt)=1
       EndIf
      Else
       If(NBond(KAt).gt.1) then
        LAt=IZ(1,KAt)
        IZ(4,IAt)=0
       Else
        LAt=IBond(2,JAt)
        If(LAt.eq.KAt) LAt=IBond(1,JAt)
        IZ(4,IAt)=1
       EndIf
      EndIf
      IZ(3,IAt)=LAt
      Return
      End
*Deck SmiBnd
      Subroutine SmiBnd(IOut,MxBnd,NAtoms,IAn,IArom,IBranch,IRing,
     $  NBond,IBond,BndOrd)
      Implicit Real*8 (A-H,O-Z)
      PARAMETER(MaxEl=200)
      Dimension IAn(*),IArom(*),NBond(*),IBond(MxBnd,*)
      Dimension IBranch(*),IRing(*)
      Dimension BndOrd(MxBnd,*)
      Character Star*1
      Dimension IEl(0:MaxEl)
      Call FillEl(0,MaxEl,IEl)
      write(IOut,'(/,'' Aromatic Atoms are marked by * '')')
      write(IOut,'('' Atom Branch Ring'',10X,''Bonded Atoms'',10X,
     $  ''Bond Orders'')')
      DO IAt = 1, NAtoms
       NBI = NBond(IAt)
       IAnI=IAn(IAt)
       If(IArom(IAt).ne.0) then
        star='*'
       Else
        star=' '
       EndIf
       If(NBI.eq.1) then
         write(IOut,'(1X,A2,I2,2I5,1X,A1,7X,I3,12X,F8.4)')
     $     IEL(IAnI),IAt,IBranch(IAt),IRing(IAt),Star,
     $     (IBond(J,IAt),J=1,NBI),(BndOrd(J,IAt),J=1,NBI)
       ElseIf(NBI.eq.2) then
         write(IOut,'(1X,A2,I2,2I5,1X,A1,7X,2I3,9X,2F8.4)')
     $     IEL(IAnI),IAt,IBranch(IAt),IRing(IAt),Star,
     $     (IBond(J,IAt),J=1,NBI),(BndOrd(J,IAt),J=1,NBI)
       ElseIf(NBI.eq.3) then
         write(IOut,'(1X,A2,I2,2I5,1X,A1,7X,3I3,6X,3F8.4)')
     $     IEL(IAnI),IAt,IBranch(IAt),IRing(IAt),Star,
     $     (IBond(J,IAt),J=1,NBI),(BndOrd(J,IAt),J=1,NBI)
       ElseIf(NBI.eq.4) then
         write(IOut,'(1X,A2,I2,2I5,1X,A1,7X,4I3,3X,4F8.4)')
     $     IEl(IAnI),IAt,IBranch(IAt),IRing(IAt),Star,
     $     (IBond(J,IAt),J=1,NBI),(BndOrd(J,IAt),J=1,NBI)
       Else
         write(IOut,'('' Wrong Number of Bonds ('',I1,'' ) for Atom'',
     $     I3)') NBI,IAt
          stop
         EndIf
      END DO
      write(IOut,'('' '')')
      Return
      End 
*Deck Valence
      FUNCTION Valence(AtomSy)
        CHARACTER*2 AtomSy
        INTEGER Valence,IArom
        IF (AtomSy .EQ. 'H ') THEN
           Valence = 1    
        ELSE IF (AtomSy(1:1) .EQ. 'c') THEN
           Valence = 3
        ELSE IF (AtomSy .EQ. 'C ') THEN  
           Valence = 4
        ELSE IF (AtomSy(1:1) .EQ. 'n') THEN
           Valence = 2
        ELSE IF (AtomSy .EQ. 'N ') THEN  
           Valence = 3
        ELSE IF (AtomSy .EQ. 'O ') THEN
           Valence = 2
        ELSE IF (AtomSy .EQ. 'F ') THEN
           Valence = 1
        ELSE IF (AtomSy .EQ. 'Cl') THEN
           Valence = 1
        ELSE IF (AtomSy .EQ. 'Br') THEN
           Valence = 1
        ELSE IF (AtomSy .EQ. 'S ') THEN
           Valence = 2
        ELSE
           Valence = 0
        END IF
      END FUNCTION
*Deck AddHyd
      Subroutine AddHyd(MxAt,MxBnd,NHeavy,NAtoms,IAn,NPi,NBond,IBond,
     $  IZ,BndOrd,Alpha,AtomTy)
      Implicit Real*8 (A-H,O-Z)
C Add Hydrogens
      Dimension BndOrd(MxBnd,*),Alpha(*)
      Dimension NBond(*),IAn(*),NPi(*),IBond(MxBnd,*),IZ(4,*)
      CHARACTER*2 AtomTy(MxAt)
      Integer Valence
      IHAt = NHeavy
      DO 10 IAt=1,NHeavy
       NSig = Valence(AtomTy(IAt))-NPi(IAt)
       NBNH=NBond(IAt)
       MissH=NSig-NBNH
       IF(MissH.gt.0) then 
        Do 20 JAt = 1, MissH
         IHAt = IHAt + 1
         IAn(IHAt)=1 
         NBond(IAt)=NBond(IAt)+1
         NBond(IHAt)=1
         IBond(NBond(IAt),IAt)=IHAt 
         BndOrd(NBond(IAt),IAt)=1.0d0
         IBond(1,IHAt)=IAt
         BndOrd(1,IHAt)=1.0d0
         call SetIZ(MxBnd,IHAt,IAt,NBond,IBond,IZ)
         If(NSig.eq.3) Alpha(IHAt)=120.0d0
  20    Continue
       EndIf 
  10  Continue
      NAtoms=IHAt
      Return
      End 
*Deck SetBLA
      Subroutine SetBLA(MxBnd,NAtoms,IAn,IZ,NBond,IBond,BndOrd,BL,Alpha)
      Implicit Real*8 (A-H,O-Z)
C Set standard bond lengths and valence angles
      Dimension IAn(*),IZ(4,*),NBond(*),IBond(MxBnd,*)
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
        EndIf
       EndIf
      End Do
      RETURN
      End
*Deck SetDih
      Subroutine SetDih(NAtoms,IZ,NBond,NAtCyc,IBranch,IRing,Alpha,Beta)
      Implicit Real*8 (A-H,O-Z)
      Logical Fnd(100)
      Dimension IZ(4,*),NBond(*),IBranch(*),IRing(*)
      Dimension Alpha(*),Beta(*),NAtCyc(*)
      Do 10 IAt=1,NAtoms
       Fnd(IAt)=.False.
       If(IZ(4,IAt).eq.1) then
        If(NBond(IZ(1,IAt)).eq.4) then
         Beta(IAt)=109.47122063D0
        ElseIf(NBond(IZ(1,IAt)).eq.3) then
         Beta(IAt)=120.0D0
        EndIf
       EndIf
       If(IRing(IAt).eq.0) goto 10
       IR0=IRing(IAt)
       IR1=IRing(IZ(1,IAt))
       IR2=IRing(IZ(2,IAt))
       IR3=IRing(IZ(3,IAt))
       If(IR0.eq.IR1.and.IR1.eq.IR2) then
C Set valence angles for cycles
        SNat=Float(NAtCyc(IR1)) 
        Alpha(IAt)=180.0*(SNat-2.0d0)/SNat
C Set dihedral angles for cycles
        If(IR2.eq.IR3) Beta(IAt)=0.0D0  
       EndIf 
   10 Continue   
      Do 20 IAt=5,NAtoms
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
           IZ(4,IAt)=-1
          ElseIf(Fnd(JAt)) then
           Beta(IAt)=Beta(JAt)-120.0D0
          Else
           Beta(IAt)=Beta(JAt)+120.0D0
          EndIf
         ElseIf(NBond(IZI1).eq.3) then
          If(Abs(IZI4).eq.1) then
           IZ(4,IAt)=-1
          Else
           Beta(IAt)=Beta(JAt)+180.0d0
          EndIf
         EndIf
         Fnd(JAt)=.true.
        EndIf
   30  Continue
       If(Beta(IAt).ge.360.0d0) Beta(IAt)=Beta(IAt)-360.0D0       
   20 Continue
      Return
      End
*Deck PropagateRing
      Recursive Subroutine PropagateRingID(MxAt,MxBnd,NBond,IBond,IRing, 
     $  CurrentRingID, StartAt, EndAt)
      IMPLICIT Integer (A-Z)
      Integer MxAt, NBond(*), IBond(MxBnd,*),IRing(*),StartAt
      Integer EndAt, I, J
! Propagate the ring ID to all atoms connected to StartAt or EndAt
      Do I = 1, MxAt
       If (IRing(I) .eq. 0) then
! Check if atom is connected to StartAt or EndAt
        Do J = 1, NBond(I)
         If (IBond(J, I) .eq. StartAt .or. IBond(J, I) .eq. EndAt) then
          IRing(I) = CurrentRingID
! Recursively propagate to connected atoms
          Call PropagateRingID(MxAt,MxBnd,NBond,IBond,IRing, 
     $     CurrentRingID,StartAt,I)
         End If
        End Do
       End If
      End Do
      Return
      End
*Deck SetCharge
      Subroutine SetCharge(IAt,LCharg,IQAt)
      IMPLICIT Integer (A-Z)
      Dimension IQAt(*)
      IQAt(IAt) = LCharg
      Return
      End
