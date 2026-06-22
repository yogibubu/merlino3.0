*Deck GICPrune
C=======================================================================
C  GICPRUNE
C
C  Final type-local redundancy pruning for GICForge.
C
C  The reduction is intentionally block-local: bends, linear bends,
C  torsions, and out-of-plane coordinates are pruned independently.
C  Stretchings are kept as primitive bond coordinates and are not
C  pruned here.
C
C  Method: modified Gram-Schmidt over the B rows of one block.  A coordinate is
C  retained if its B row increases the numerical rank of its own type block.
C=======================================================================

      Subroutine PruneGICBlocks(IOut,IPrint,MxAtP,MxTrm,NAtoms,NLen,
     $ NAng,NLAng,NOupl,NDih,NTermB,NTermA,NTermL,NTermD,NTermO,
     $ IAtomB,IAtomA,IAtomL,IAtomD,IAtomO,IPrimB,IPrimA,IPrimL,
     $ IPrimD,IPrimO,ITVB,ITVA,ITVLA,ITVD,ITVO,IFixB,IFixA,IFixL,
     $ IFixD,IFixO,CoefB,CoefA,CoefL,CoefD,CoefO,ValTB,ValTA,ValTL,
     $ ValTD,ValTO,C,DoBMat,BMat,Scr)
      Implicit Real*8 (A-H,O-Z)
      Integer MxAtP,MxTrm,NAtoms,NLen,NAng,NLAng,NOupl,NDih
      Dimension NTermB(*),NTermA(*),NTermL(*),NTermD(*),NTermO(*)
      Dimension IAtomB(MxAtP,MxTrm,*),IAtomA(MxAtP,MxTrm,*)
      Dimension IAtomL(MxAtP,MxTrm,*),IAtomD(MxAtP,MxTrm,*)
      Dimension IAtomO(MxAtP,MxTrm,*)
      Dimension IPrimB(MxTrm,*),IPrimA(MxTrm,*),IPrimL(MxTrm,*)
      Dimension IPrimD(MxTrm,*),IPrimO(MxTrm,*)
      Dimension ITVB(*),ITVA(*),ITVLA(*),ITVD(*),ITVO(*)
      Dimension IFixB(*),IFixA(*),IFixL(*),IFixD(*),IFixO(*)
      Dimension CoefB(MxTrm,*),CoefA(MxTrm,*),CoefL(MxTrm,*)
      Dimension CoefD(MxTrm,*),CoefO(MxTrm,*)
      Dimension ValTB(*),ValTA(*),ValTL(*),ValTD(*),ValTO(*)
      Dimension C(3,*),BMat(3*NAtoms,*),Scr(*)
      Logical Keep(1000)
      Logical DoB1,DoBMat
      Character*16 Label

      If(NAtoms.le.0) return
      NTot=NLen+NAng+NLAng+NDih+NOupl
      If(NTot.le.0) return

      DoB1=.False.
      Call MkBNew(IOut,0,DoB1,MxAtP,MxTrm,NAtoms,NLen,NAng,NLAng,
     $ NOupl,NDih,IAtomB,IAtomA,IAtomL,IAtomD,IAtomO,NTermB,NTermA,
     $ NTermL,NTermD,NTermO,CoefB,CoefA,CoefL,CoefD,CoefO,C,BMat)

      Write(IOut,'(/,'' Type-local residual GIC redundancy pruning'')')
      Write(IOut,'(''   Method: modified Gram-Schmidt on B rows;'',
     $ '' coordinate types are not mixed.'')')

      Write(IOut,'(''   Stretch: kept all '',I5,
     $ '' primitive coordinates.'')') NLen

      Label='Bend'
      IOff=NLen
      Call PruneOneBlock(IOut,Label,NAtoms,NAng,IOff,BMat,Keep,Scr,
     $ NKeep)
      If(NKeep.lt.NAng) Call PackGICBlock(MxAtP,MxTrm,NAng,Keep,
     $ NTermA,IAtomA,IPrimA,ITVA,IFixA,CoefA,ValTA)

      Label='Linear bend'
      IOff=NLen+NAng
      Call PruneOneBlock(IOut,Label,NAtoms,NLAng,IOff,BMat,Keep,Scr,
     $ NKeep)
      If(NKeep.lt.NLAng) Call PackGICBlock(MxAtP,MxTrm,NLAng,Keep,
     $ NTermL,IAtomL,IPrimL,ITVLA,IFixL,CoefL,ValTL)

      Label='Torsion'
      IOff=NLen+NAng+NLAng
      Call PruneOneBlock(IOut,Label,NAtoms,NDih,IOff,BMat,Keep,Scr,
     $ NKeep)
      If(NKeep.lt.NDih) Call PackGICBlock(MxAtP,MxTrm,NDih,Keep,
     $ NTermD,IAtomD,IPrimD,ITVD,IFixD,CoefD,ValTD)

      Label='Out-of-plane'
      IOff=NLen+NAng+NLAng+NDih
      Call PruneOneBlock(IOut,Label,NAtoms,NOupl,IOff,BMat,Keep,Scr,
     $ NKeep)
      If(NKeep.lt.NOupl) Call PackGICBlock(MxAtP,MxTrm,NOupl,Keep,
     $ NTermO,IAtomO,IPrimO,ITVO,IFixO,CoefO,ValTO)

      Write(IOut,'(''   Final active GIC counts:'')')
      Write(IOut,'(''     Stretch='',I5,'' Bend='',I5,'' Linear='',I5,
     $ '' Torsion='',I5,'' Out-of-plane='',I5)') NLen,NAng,NLAng,NDih,
     $ NOupl
      If(DoBMat) then
       NTot=NLen+NAng+NLAng+NDih+NOupl
       Call MkBNew(IOut,0,DoB1,MxAtP,MxTrm,NAtoms,NLen,NAng,NLAng,
     $ NOupl,NDih,IAtomB,IAtomA,IAtomL,IAtomD,IAtomO,NTermB,NTermA,
     $ NTermL,NTermD,NTermO,CoefB,CoefA,CoefL,CoefD,CoefO,C,BMat)
       Call WriteGICBMat(NAtoms,NTot,BMat)
       Write(IOut,'(''   Machine-readable final B matrix: bmat.out'')')
      EndIf
      Return
      End

*Deck WriteGICBMat
      Subroutine WriteGICBMat(NAtoms,NInt,BMat)
      Implicit Real*8 (A-H,O-Z)
      Dimension BMat(3*NAtoms,*)
      NCart=3*NAtoms
      Open(77,File='bmat.out',Status='Unknown')
      Rewind(77)
      Write(77,'(A)') '# merlino.gicforge.bmatrix.v1'
      Write(77,'(A)') '# row col value; rows are final GICs'
      Write(77,'(2I8)') NInt,NCart
      Do 20 IInt=1,NInt
       Do 10 ICart=1,NCart
        Write(77,'(2I8,1X,D24.16)') IInt,ICart,BMat(ICart,IInt)
   10  Continue
   20 Continue
      Close(77)
      Return
      End

*Deck PruneOneBlock
      Subroutine PruneOneBlock(IOut,Label,NAtoms,NVar,IOff,BMat,Keep,
     $ Scr,NKeep)
      Implicit Real*8 (A-H,O-Z)
      Integer IOut,NAtoms,NVar,IOff,NKeep
      Character*(*) Label
      Dimension BMat(3*NAtoms,*),Scr(3*NAtoms,*)
      Logical Keep(*)
      Logical AnyRem
      Real*8 Norm,Norm0,Dot,TAbs,TRel
      Data TAbs/1.0D-10/, TRel/1.0D-08/

      NCart=3*NAtoms
      NKeep=0
      If(NVar.le.0) return

      Do 10 I=1,NVar
       Keep(I)=.False.
   10 Continue

      Do 100 I=1,NVar
       Row=IOff+I
       Norm0=0.0D0
       Do 110 K=1,NCart
        Scr(K,NKeep+1)=BMat(K,Row)
        Norm0=Norm0+Scr(K,NKeep+1)*Scr(K,NKeep+1)
  110  Continue
       Norm0=DSqrt(Norm0)
       If(Norm0.le.TAbs) go to 100

       Do 130 J=1,NKeep
        Dot=0.0D0
        Do 120 K=1,NCart
         Dot=Dot+Scr(K,NKeep+1)*Scr(K,J)
  120   Continue
        Do 125 K=1,NCart
         Scr(K,NKeep+1)=Scr(K,NKeep+1)-Dot*Scr(K,J)
  125   Continue
  130  Continue

       Norm=0.0D0
       Do 140 K=1,NCart
        Norm=Norm+Scr(K,NKeep+1)*Scr(K,NKeep+1)
  140  Continue
       Norm=DSqrt(Norm)
       If(Norm.gt.TAbs.and.Norm.gt.TRel*Norm0) then
        NKeep=NKeep+1
        Keep(I)=.True.
        Do 150 K=1,NCart
         Scr(K,NKeep)=Scr(K,NKeep)/Norm
  150   Continue
       EndIf
  100 Continue

      If(NKeep.lt.NVar) then
       Write(IOut,'(''   '',A,'': kept '',I5,'' of '',I5,
     $ '' coordinates; removed '',I5,'' residual redundant.'')')
     $ Label,NKeep,NVar,NVar-NKeep
       Write(IOut,'(''     Removed local indices:'')')
       AnyRem=.False.
       Do 210 I=1,NVar
        If(.not.Keep(I)) then
         Write(IOut,'(I6)',advance='no') I
         AnyRem=.True.
        EndIf
  210  Continue
       If(AnyRem) Write(IOut,'('' '')')
      Else
       Write(IOut,'(''   '',A,'': kept all '',I5,
     $ '' coordinates.'')') Label,NVar
      EndIf
      Return
      End

*Deck PackGICBlock
      Subroutine PackGICBlock(MxAtP,MxTrm,NVar,Keep,NTerm,IAtom,IPrim,
     $ ITPV,IFixG,Coef,ValTot)
      Implicit Real*8 (A-H,O-Z)
      Integer MxAtP,MxTrm,NVar
      Integer I,J,T,A
      Dimension NTerm(*),IAtom(MxAtP,MxTrm,*),IPrim(MxTrm,*),ITPV(*)
      Dimension IFixG(*),Coef(MxTrm,*),ValTot(*)
      Logical Keep(*)

      J=0
      Do 100 I=1,NVar
       If(.not.Keep(I)) go to 100
       J=J+1
       If(J.ne.I) then
        NTerm(J)=NTerm(I)
        ITPV(J)=ITPV(I)
        IFixG(J)=IFixG(I)
        ValTot(J)=ValTot(I)
        Do 30 T=1,MxTrm
         IPrim(T,J)=IPrim(T,I)
         Coef(T,J)=Coef(T,I)
         Do 20 A=1,MxAtP
          IAtom(A,T,J)=IAtom(A,T,I)
   20    Continue
   30   Continue
       EndIf
  100 Continue
      NVar=J
      Return
      End
