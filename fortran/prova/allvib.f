C -----------  Least Squares Routines ----------
*Deck LnLstq
      Subroutine LnLstq(IOut,IPrint,NPts,NPrm,ChiSq,X,Y,Sig,IPot,
     $  AFunc,B,A,AtA,AtAm1,D,Beta,Cof,Delt,MxScr,IScr)
      Implicit Real*8 (A-H,O-Z)
      Logical SVD
C Given X,Y, and Sig values at NPts points build the design matrix for
C LSQ polynomial fitting of degree NPrm-1
C Afunc has dimension NPrm, B has dimension NPts and U(NPts,NPrm)
      Dimension X(*),Y(*),Sig(*),IPot(*),Cof(*),Delt(*)
      Dimension Afunc(*),B(*),A(NPts,NPrm),AtA(NPrm,NPrm)
      Dimension AtAm1(NPrm,NPrm),D(NPrm,NPrm),Beta(NPrm)
      Dimension IScr(*)
      Do 10 I=1,NPts
       XX=X(I)
       call FPoly(XX,AFunc,NPrm,IPot)
       Tmp=1.0D0/Sig(i)
       Do 20 J=1,NPrm
        A(i,j)=AFunc(j)*Tmp
   20  Continue
       B(i)=Y(i)*Tmp
   10 Continue
      Do 30 I=1,NPrm
       Do 40 J=1,NPrm
        AtA(I,J)=0.0d0
        Do 50 K=1,NPts
         AtA(I,J)=AtA(I,J)+A(K,I)*A(K,J) 
   50   Continue
   40  Continue
   30 Continue   
      Do 60 I=1,NPrm
       Beta(i)=0.0d0
       Do 70 J=1,NPts
        Beta(i)=Beta(i)+A(J,I)*B(J)
   70  Continue
   60 Continue
      SVD=.False.
      Call AMove(NPrm*NPrm,AtA,AtAm1)
      Call MkAAm1(IOut,IPrint,SVD,NPrm,AtAm1,D,IScr)
      Do 80 i=1,NPrm
       Cof(i)=0.0d0
       Do 85 j=1,NPrm
        Cof(i)=Cof(i)+AtAm1(i,j)*Beta(j)
   85  Continue
   80 Continue
      ChiSq=0.0D0
      Do 90 IPts=1,NPts
       DY = Y(IPts) - YPoly(X(IPts),NPrm,IPot,Cof)
       Chisq = ChiSq + DY*DY
   90 Continue
      ChiSq=ChiSq/DFloat(NPts-NPrm)
      If (IPrint.gt.0) Write(IOut,'(/,'' Chi Square ='',D12.5)') ChiSq
      Write(IOut,'(/,'' Variance-Covariance Matrix'')')
      Write(IOut,'(6X,''IPot'',2X,6(I5,7X))')(IPot(I),I=1,NPrm)
      Do 100 I=1,NPrm
       Delt(i)=SQrt(ChiSq*AtAM1(i,i))
       write(IOut,'(6X,I2,2X,6D12.4)') IPot(I),(ChiSq*AtAm1(I,J),
     $   J=1,NPrm)
  100 Continue
      If(IPrint.eq.0) Return
      write(IOut,'(/,'' Fitted Parameters '')')
      Do 110 i=1,NPrm
       Write(IOut,'(F12.5,'' X**'',I1,''('',D12.5,'')'')') 
     $   Cof(I),IPot(I),Delt(i)
  110 Continue
      Return
      End
*Deck MkAAm1
      Subroutine MkAAm1(IOut,IPrint,SVD,NPts,ATA,D,IScr)
      Implicit Real*8 (A-H,O-Z)
      Logical Inv1,SVD
      Dimension ATA(*),D(*),IScr(*)
      If(SVD) then
C      Call GenInv(IOut,IPrint,NPts,NPts,ATA,IScr,D)
       write(IOut,'('' No GenInv yet'')')
       stop
      Else
       InIS=1
       InIAD1=INIS+2*NPts
       InIAd2=InIAd1+NPts
       If(.not.Inv1(ATA,NPts,IScr(InIS),IScr(InIAD1),IScr(InIAD2),D,
     $   NPts,Det)) then
        Write(IOut,'('' Inversion of ATA Matrix Failed'')')
        Stop
       EndIf
      EndIf
      If(IPrint.eq.0) Return
      Write(IOut,'(/,'' (A+A)-1 Matrix'')')
      Do 150 IPts=1,NPts
       Write(IOut,'('' Data Points:'',I5)') IPts
       Ini=NPts*(IR-1)+1
       IEnd=Ini+NPts-1
       Write(IOut,'(6F10.5)') (ATA(ii),ii=ini,iend)
  150 Continue
      Return
      End
*Deck FPoly
      Subroutine FPoly(X,P,NPrm,IPot)
      Implicit Real*8 (A-H,O-Z)
      Dimension P(*),IPot(*)
      Do 10 I=1,NPrm
       P(i)=X**(IPot(i))
   10 Continue
      Return
      End
*Deck YPoly
      Real*8 Function YPoly(X,NPrm,IPot,Cof)
      Implicit Real*8 (A-H,O-Z)
      Dimension IPot(*),Cof(*)
      YPoly=0.0d0
      Do 10 I=1,NPrm
       YPoly=YPoly+Cof(i)*X**(IPot(i))
   10 Continue
      Return
      End
*Deck SetPts
      Subroutine SetPts(IOut,IPrint,NPts,NBas,NPrm,NEig,X,Y,IPot,Cof,
     $  XDVR,YDVR)
      Implicit Real*8 (A-H,O-Z)
      Dimension X(*),Y(*),XDVR(*),YDVR(*)
      Dimension IPot(*),Cof(*)
      IEMin=0
      IEMax=0
      EMin=1.0d10
      EMax=-1.0D10
      Write(IOut,'(/,I5,'' Original Points'')') NPts
      Do 10 IPts=1,NPts
       If(Y(IPts).lt.EMin) then
        EMin=Y(IPts)
        IEMin=IPts
       ElseIf(Y(IPts).gt.EMax) then
        EMax=Y(IPts)
        IEMax=IPts
       EndIf
       Write(IOut,'(D12.5,3X,D12.5)') X(IPts),Y(IPts)
   10 Continue
      Write(IOut,'(/,''EMin'',D12.5,''('',I3,'')'')')EMin,IEMin
      Write(IOut,'(''EMax'',D12.5,''('',I3,'')'')')EMax,IEMax
      DX1=Abs(X(IEMin)-X(IEmin-1))
      DX2=ABs(X(IEMin)-X(IEMin+1))
      DX=(DX1+DX2)/DFloat(2)
      FK=(Y(IEMin-1)+Y(IEMin+1)-2.0D0*Y(IEMin))/(DX*DX)
      EMax=DFLoat(NEig/3)*SQrt(FK) 
      write(IOut,'('' FK,Emax'',2D12.5)') FK,EMax
      XMax=SQrt(EMax/FK)/2.0d0-0.5d0
      XMin=-XMax
      DX=(XMax-XMin)/(NBas-1)
      Write(IOut,'('' XMin,XMax'',2D12.5)') XMin,XMax
      Write(IOut,'(I3,'' IPot,Cof'')') NPrm
      Write(IOut,'(6I3)') (IPot(I),I=1,NPrm)
      Write(IOut,'(6D12.5)') (Cof(I),I=1,NPrm)
      Write(IOut,'(/,I5,'' Fitted Points'')') NBas 
      do 20 I=1,NBas
       XDVr(I)=XMin+DFloat(I-1)*DX
       YDVr(I)=YPoly(XDVr(I),NPrm,IPot,Cof)
       Write(IOut,'(2D12.5)') XDVr(I),YDVr(I)
   20 Continue
      Return
      End
C ------------  Start of VCI1 Routines ---------------
*Deck Var1D 
      Subroutine Var1D(IOut,IPrint,NBas,NEig,Alpha,F,G,V,IScr)
      Implicit Real*8 (A-H,O-Z)
C Input
      Dimension F(*),G(*),V(*)
      Dimension IScr(*)
C Split scratch for diagonalization
C V(I1) H  matrix (lower triangular)
      I1=1
C V(I2) Eigenvalues
      I2=I1+NBas*(NBas+1)/2
C V(I3)=EVec
      I3=I2+NBas
C V(I4)=WA or Scratch
      I4=I3+NBas*NBas
C V(I5)=Scr
      I5=I4+6*NBas
      I6=I5+NBas*(NBas+1)/2+1
C
      Thresh=1.0D-05
      NMax=4
      IJ=0
      Call Aclear(NBas*(NBas+1)/2,V(I1))
      Do 10 I=1,NBas
       NI=I-1
       Do 20 J=1,I  
        IJ=IJ+1
        NJ=J-1
        V(IJ)=HIJ(IOut,IPrint,NI,NJ,NMax,Alpha,F,G)
        If(DAbs(V(IJ)).lt.Thresh) V(IJ)=0.0D0
   20  Continue
   10 Continue
      Call HQRII1(IOut,NBas,1,NBas,-1,V(I1),V(I2),NBas,V(I3),.true.,
     $  IErr,IScr,V(I4),V(I5),I6-I5)
      If(IErr.ne.0) then
       write(IOut,'(/,'' After HQRII1: IErr ='',I3)') IErr
       Stop
      EndIf
      Harm0=F(2)/2.0D0
      HarmI=F(2) 
      X0=F(4)/6.4D+01-7.0D0*F(3)*F(3)/(6.4D+01*9.0D0*F(2)) 
      XI=F(4)/1.6D+01-5.0D0*F(3)*F(3)/(4.8D+01*F(2))
      Pert0=Harm0+X0+XI/4.0D0  
      Do 40 IBas=1,NEig
        Ini=(IBas-1)*NBas+I3
       If(IBas.eq.1) then
        IF(DAbs(F(3)).lt.DAbs(F(2)).and. DAbs(F(4)).lt.DAbs(F(2))) then
         write(IOut,'(/,'' Ground State ZPE (cm-1):'',6X,''Harmonic'',
     $    F10.2,5X,''VPT2'',F10.2,5X,''Variational'',F10.2)') 
     $    DABs(Harm0),Pert0,V(I2)
        Else
         write(IOut,'(/,'' Ground State ZPE (cm-1):'',6X,''Harmonic'',
     $    F10.2,5X,''Variational'',F10.2)') DAbs(Harm0),V(I2)
        EndIf
        call Clean(IOut,NBas,MaxTrm,V(Ini))
       Else
        SNI=DFloat(IBas-1)
        PertI=SNI*F(2)+(SNI*SNI+SNI)*XI
        IF(DAbs(F(3)).lt.DAbs(F(2)).and. DAbs(F(4)).lt.DAbs(F(2))) then
         write(IOut,'(/,'' State'',I3,'' Omega(I)-ZPE (cm-1):'',
     $    1X,''Harmonic'',F10.2,5X,''VPT2'',F10.2,5X,''Variational'',
     $    F10.2)') IBas-1,DABs(SNI*HarmI),PertI,(V(I2+IBas-1)-V(I2))
        Else
         write(IOut,'(/,'' State'',I3,'' Omega(I)-ZPE (cm-1):'',
     $    1X,''Harmonic'',F10.2,5X,''Variational'',
     $    F10.2)') IBas-1,DABs(SNI*HarmI),(V(I2+IBas-1)-V(I2))
        EndIf 
        call Clean(IOut,NBas,MaxTrm,V(Ini))
       EndIf
   40 Continue
      Return
      End
*Deck HIJ 
      Real*8 Function HIJ(IOut,IPrint,NI,NJ,NMax,Alpha,f,g)
      Implicit Real*8 (A-H,O-Z) 
C     +-------------------------------------------------------------------+
C        POTENTIAL OPERATORS (Hf2,Hf3,etc. WITH COEFFICIENTS)
C      + KINETIC OPERATORS   (Hg2,Hg3,etc. where the number is the power of P +
C                             the power of Q WITH COEFFICIENTS, i.e. 1/2P**2)
C     +-------------------------------------------------------------------+
C     WARNING: Values OK for a potential in Phi.
C              HOWEVER, if a K potential is available, the factors 2!,
C              3!... (collected in DenF2,DenF3, etc.) must be removed
C              In the same vein the factors DenD2, etc. must be removed
C              alpha allows to use basis functions with an exponent different 
C              from that optimized for the quadratic term
C
      Dimension f(*),g(*)
      Save Zero,Pt5,F1,F2,F3,F4,F5,F6,F8,Ten,F12,F15,F16
      Data Zero/0.0D0/,Pt5/5.0D-01/,F1/1.0D0/,F2/2.0D0/,
     $  F3/3.0D0/,F4/4.0D0/,F5/5.0D0/,F6/6.0D0/,F8/8.0D0/,
     $  Ten/1.0D1/,F12/1.2D1/,F15/1.5D1/,F16/1.6D1/
      If(Alpha.lt.1.0d-01) Alpha=F1
      HIJ=Zero
      IDelt=IAbs(NI-NJ)
      If(IDelt.gt.NMax) Return
      NMin=Min0(NI,NJ)
      Div=DFloat(NMin) 
C -- Diagonal Terms --
      If(IDelt.eq.0) then
       DenD2=F2
       DenF2=F2
       DenD4=F4
       DenF4=F2*F3*F4
       DenD6=F4*F2
       DenF6=F2*F3*F4*F5*F6
       Hf2=f(2)*(F2*div+F1)/(alpha**2*DenD2*DenF2)
       Hg2=g(2)*(F2*div+F1)*alpha**2/(DenD2*DenF2)
       Hf4=f(4)*(div*div+div+pt5)/(F16*alpha**4)
       Hg4=g(4)*(div*div+div+F1+pt5)*alpha**4/F8
       SNumf6=F5*(F3+F8*div+F6*div*div+F4*div*div*div)
       Hf6=f(6)*SNumf6/(DenD6*DenF6*alpha**6)
       Hg6=g(6)*F6*alpha**4/(DenD6*DenF6*alpha**6)
       HIJ=Hf2+Hg2+Hf4+Hg4+Hf6+Hg6 
       Return
C -- i,i+1 Terms  --
      ElseIf(IDelt.eq.1) then
       DenD1=DSqrt(F2)
       DenF1=F1
       Den1=DenD1*DenF1
       SNum1=DSQrt(div+F1)
       Hf1 = f(1)*SNum1/Den1
       DenD3=F2*DSqrt(F2)
       DenF3=F2*F3
       Den3=DenD3*DenF3*alpha**3
       SNumf3=F3*DSQrt((div+F1)**3)
       Hf3 = f(3)*SNumf3/Den3
       Hg3 = g(3)*Snumf3*alpha**4/Den3
       DenD5= F4*DSQrt(F2)
       DenF5 = F2*F3*F4*F5
       Den5=DenD5*DenF5*alpha**5
       SNumf5=F5*(F2*div*div+F4*div+F3)*DSQrt(Div+F1)
       Hf5 = f(5)*SNumf5/Den5
       Hg5 = g(5)*F4*alpha**4*SNumf5/Den5
       Hij = Hf1+Hf3+Hg3+Hf5+Hg5
       Return
C -- i,i+2 TERMS --
      ElseIf(IDelt.eq.2) then
       DenD2=F2
       DenF2=F2
       Den2=DenD2*DenF2*alpha**2
       SNumf2=DSQrt((div+F1)*(div+F2))
       Hf2 = f(2)*SNumf2/Den2
       Hg2 = g(2)*SNumf2*alpha**4/Den2
       DenD4=F4
       DenF4=F2*F3*F4
       den4=DenD4*DenF4*alpha**4
       Snumf4=(F4*div+F6)*DSQrt((div+F1)*(div+F2))
       Hf4 = f(4)*SNumf4/Den4
       Hg4 = g(4)*SNumf4*alpha**4/Den4
       DenD6=F4*F2
       DenF6=F2*F3*F4*F5*F6
       Den6=DenD6*DenF6*alpha**6
       SNumf6=F15*(div*div+F3*div+F3)*DSQrt((div+F1)*(div+F2))
       Hf6 = f(6)*SNumF6/Den6
       Hg6 = g(6)*F2*SNumF6*alpha**4/Den6
       HIJ = Hf2-Hg2+Hf4+Hg4+Hf6
       Return
C --i,i+3 Terms --
      ElseIf(IDelt.eq.3) then
       DenD3=F2*DSQrt(F2)
       DenF3=F2*F3
       Den3=DenD3*DenF3*alpha**3
       SNumf3=DSQrt((div+F1)*(div+F2)*(div+F3))
       Hf3 = f(3)*SNumf3/Den3
       Hg3 = g(3)*F3*SNumf3*alpha**4/Den3
       DenD5=F4*DSQrt(F2)
       DenF5=F2*F3*F4*F5
       Den5=DenD5*DenF5*alpha**5
       SNumf5=F5*(div+F2)*DSQRt((div+F1)*(div+F2)*(div+F3))
       Hf5 = f(5)*SNumf5/Den5
       Hg5 = g(5)*F4*SNumf5*alpha**4/Den5
       HIJ = Hf3-Hg3+Hf5-Hg5 
       Return
C --i,i+4 Terms
      ElseIf(IDelt.eq.4) then
       DenD4=F4
       DenF4=F2*F3*F4
       Den4=DenD4*DenF4*alpha**4
       SNumf4=DSQrt((div+F1)*(div+F2)*(div+F3)*(div+F4))
       Hf4 = f(4)*SNumf4/Den4
       Hg4 = g(4)*F6*SNumf4*alpha**4/Den4
       DenD6=F4*F2
       DenF6=F2*F3*F4*F5*F6
       Den6=DenD6*DenF6*alpha6
       SNumf6=F3*(F2*div+F5)*DSQrt((div+F1)*(div+F2)*(div+F3)*(div+F4))
       Hf6 = f(6)*SNumf6/Den6
       Hg6 = g(6)*F10*alpha**4*SNumf6/Den6
       HIJ = Hf4-Hg4+HF6-Hg6
       Return
C --i,i+5 Terms
      ElseIf(IDelt.eq.5) then
       DenD5=F4*DSQrt(F2)
       DenF5=F2*F3*F4*F5
       Den5=DenD5*DenF5*alpha**5
       SNumf5=DSQrt((div+F1)*(div+F2)*(div+F3)*(div+F4)*(div+F5))
       Hf5 = f(5)*SNumf5/Den5
       Hg5 = g(5)*F4*F15*SNumf5*alpha**5/Den5
       HIJ = Hf5-Hg5
       Return
C --i,i+6 Terms
      ElseIf(IDelt.eq.6) then
       DenD6=F4*F2
       DenF6=F2*F3*F4*F5*F6
       Den6=DenD6*DenF6*alpha**6
       SNumf5=DSQrt((div+F1)*(div+F2)*(div+F3)*(div+F4)*(div+F5))
       SNumf6=SNumf5*DSQrt(div+F6)
       Hf6 = f(6)*SNumf6/Den6
       Hg6 = g(6)*F2*F15*SNumf6*alpha**4/Den6
       HIJ = Hf6-Hg6 
       Return
      Else
       HIJ = 0.0D0
      EndIf
      Return
      End
*Deck VFP
      Real*8 Function VFP(N,DNv)
C DNV is the double precision representation of the vibrational quantum number
      Real*8 DNv
      VFP=1.0D0
      Do 10 I=1,N
       VFP=VFP*(DNV+DFloat(I))
   10 Continue
      VFP=DSqrt(VFP)
      Return
      End
*Deck Clean
      Subroutine Clean(IOut,N,MaxTrm,V)
      Implicit Real*8 (A-H,O-Z)
      Dimension V(*)
      Thresh=1.0D-01
      I1=0
      Do 10 I=1,N
       If(DAbs(V(I)).gt.Thresh) then
        If(I.le.10) then
         Write(IOut,'(2X,F5.2,''|'',I1,''>'')',advance='no') V(I),I-1
        Else
         Write(IOut,'(2X,F5.2,''|'',I2,''>'')',advance='no') V(I),I-1
        EndIf
       EndIf
   10 Continue     
      Write(IOut,'('' '')')
      Return
      End
C ---------------- Start of DVR1 routines ------------
*Deck QVar1D 
      Subroutine QVar1D(In,IOut,IPrint,IMode,NPts,NBas,NEig,XMin,XMax,
     $  DX,EMin,EMax,X,EPot,MxScr,V,IScr)
      Implicit Real*8 (A-H,O-Z)
C If(ITyp.eq.0) -XMin,XMax (NPts points including XMin and XMax)
C If(ITyp.eq.1) -inf,inf X is not needed (uses DX)
C If(ITyp.eq.2) 0,inf X is not needed (uses DX)
C If(ITyp.eq.3) periodic X is not needed (first point = last point?)
      Data JR/0/
      Dimension X(*),EPot(*),V(*)
      Dimension IScr(*)
      au2cm1 = CnvFct('au2cm1')
      XMin=X(1)
      XMax=X(NBas)
      DX=(XMax-XMin)/DFloat(NBas+1)
      EMin=0.0D0
      EMax=EMax-EMin
      Write(IOut,'(/,'' XMin='',D12.5,'' XMax='',D12.5,'' DX='',D12.5,
     $  '' Delta(E)max='',F8.1,'' cm-1'')')XMin,XMax,DX,
     $  (EMax-EMin)*au2cm1
C
C Split scratch for diagonalization
C V(I1) H  matrix (lower triangular)
      I1=1
C V(I2) Eigenvalues
      I2=I1+NBas*(NBas+1)/2
C V(I3)=EVec
      I3=I2+NBas
C V(I4)=WA or Scratch
      I4=I3+NBas*NBas
C V(I5)=Scr
      I5=I4+6*NBas
      I6=I5+NBas*(NBas+1)/2+1
C
      Call AClear(I6,V)
      Thresh=1.0D-05
      NMax=4
      IJ=0
      Call Aclear(NBas*(NBas+1)/2,V(I1))
C Set points
C      Call SetPts(NPts,IMode,JR,DX,XMin,XMax,EMin,Epot)
C Build DVR Hamiltonian
      Call HamDVR(In,IOut,NBas,IMode,JR,DX,XMin,XMax,EMin,EPot,V(I1))
C Diagonalize Hamiltonian
      Call HQRII1(IOut,NBas,1,NBas,-1,V(I1),V(I2),NBas,V(I3),.true.,
     $  IErr,IScr,V(I4),V(I5),I6-I5)
      If(IErr.ne.0) then
       write(IOut,'(/,'' After HQRII1: IErr ='',I3)') IErr
       STOP
      EndIf
C
C Write Eigenvalues and EigenVectors
C
      Do 20 IBas=1,NEig
       Ini=(IBas-1)*NBas+I3
       If(IBas.eq.1) then
        write(IOut,'(/,'' Ground State ZPE (cm-1):  '',F8.1)')
     $   V(I2+IBas-1)*au2cm1
       Else 
        write(IOut,'(/,'' State'',I3,'' Nu(I)-ZPE (cm-1):'',F8.1)')
     $   IBas-1,(V(I2+IBas-1)-V(I2))*au2cm1
       EndIf
       call PrtDVR(IOut,NBas,MaxTrm,V(Ini))
   20 Continue
C Analyze Results
C     Call DVRAnl(IOut,IPrint,MaxPts,NPts,NEner,NEVec,JR,JMax,
C    $  VScal,IVType,XGrid,YGrid,EVal,EVec,E00,E01,We,WeXe,Be,Alpe,De)
      Return
      End
*Deck PrtDVR
      Subroutine PrtDVR(IOut,N,MaxTrm,V)
      Implicit Real*8 (A-H,O-Z)
      Dimension V(*)
      Thresh=1.0D-01
      II=-1
      Do 10 I=1,N
       IM=I-1 
       If(DAbs(V(I)).gt.Thresh) then
        II=II+1
        ITst=II/6
        IF((II-ITst*6).le.0.and.II.ne.0) write(IOut,'('' '')')
        If(I.le.10) then
         Write(IOut,'(2X,F5.2,''|'',I1,''> '')',advance='no') V(I),IM 
        Else
         Write(IOut,'(2X,F5.2,''|'',I2,''>'')',advance='no') V(I),IM
        EndIf
       EndIf
   10 Continue
      Write(IOut,'('' '')')
      Return
      End
*Deck HamDVR 
      Subroutine HamDVR(In,IOut,NPts,Ityp,JR,DX,XMin,XMax,EMin,EPot,
     $  HDVR)
      Implicit Real*8 (A-H,O-Z)
C     Input
C     -----
C     NPts   :: Number of data points
C     ITyp   :: Type of coordinate 
C               0 | equispaced DVR in general interval [XMin,XMax]
C               1 | cartesian coordinate (from -inf to +inf)
C               2 | stretching coordinate (from 0 to +inf)
C               3 | periodic variable
C     JR     :: Rotational level
C     DX     :: Interval between F2 abcissas
C     XMin   :: Minimum value in abcissa
C     XMax   :: Maximum value in abcissa
C     EMin   :: Minimum value of the potential energy
C     EPot   :: (NPts) Potential energy array
C
C     Output
C     ------
C     HDVR   :: (NPts*(NPts+1)/2) Hamiltonian matrix
C
      Dimension EPot(*),HDVR(*)
      Save F1, F2, F3, F4
      Data F1/1.0D0/, F2/2.0D0/, F3/3.0D0/, F4/4.0D0/
      Pi    = F4*ATan(F1)
      Pi2   = Pi**2
      amu2au = F1/CnvFct('au2amu')
      If(ITyp.eq.0) then
C Equispaced DVR in general interval [XMin,XMax]
       NTot = NPts + 1
       XX   = XMax - XMin
       Den  = F4*amu2au*XX**2
       Fac  = Pi/(F2*DFloat(NTot))
       SNum1i = DFloat(2*NTot**2+1)/F3
       Fac2   = (JR+1)*JR
       WW     = F2*amu2au
       IJ=0
       Do 100 i = 1, NPts
        Do 110 j = 1, i
         IJ=IJ+1
         If(i.eq.j) then
          GS2    = Sin(DFloat(2*i)*Fac)
          SNum2i = F1/(GS2**2)
          Tii    = Pi2*(SNum1i-SNum2i)/Den
          Xi     = (XMin + XX*i)/DFloat(NTot)
          Xi2    = Xi**2
          Fac1   = WW*XI2
          HDVR(ij) = Tii + EPot(i) - EMin + (Fac2/fac1)
         Else
          GS1j   = Sin(DFloat(i-j)*Fac)
          GS2j   = Sin(DFloat(i+j)*Fac)
          SNum1j = F1/(GS1j**2)
          SNum2j = F1/(GS2j**2)
          Tij    = (-F1)**(i-j)*Pi2*(SNum1j-SNum2j)/Den
          HDVR(ij) = Tij
         EndIf
  110   Continue
  100  Continue
      ElseIf(ITyp.eq.1) then
C Cartesian coordinate [-infinite,infinite]
       Den = F2*amu2au*DX**2
       Tii = Pi2/(F3*Den)
       IJ=0
       Do 200 i = 1, NPts
        Do 210 j = 1, i-1
         IJ=IJ+1
         If(i.eq.j) then
          HDVR(IJ) = Tii + EPot(i) - EMin
         Else 
          Tij = F2*(-F1)**(i-j)/(Den*DFloat((j-i)*(j-i)))
          HDVR(ij) = Tij
         EndIf 
  210   Continue
  200  Continue
      ElseIf(ITyp.eq.2) then
C Stretching coordinate (from 0 to +infinite)
       Den = F2*amu2au*DX**2
       Tii = Pi2/(F3*Den)
       IJ=0
       Do 300 i = 1, NPts
        Do 310 j = 1, i-1
         IJ=IJ+1
         If(i.eq.j) then
          HDVR(ij) = Tii + EPot(i) - EMin - F1/DFloat(2*i**2)
         Else 
          Tij1 = F2/DFloat((i-j)*(i-j))
          Tij2 = F2/DFloat((i+j)*(i+j))
          Tij  = (((-F1)**(i-j))*(TIJ1-TIJ2))/Den
          HDVR(ij) = Tij
         EndIf
  310   Continue
  300  Continue
      Elseif(ITyp.eq.3) then
C Periodic Variable
       Den = F2*amu2au
       NN  = (NPts-1)/2
       Tii = DFloat(NN*(NN+1))/(F3*Den)
       IJ=0
       Do 400 i = 1, NPts
        Do 410 j = 1, i-1
         IJ=IJ+1
         If(I.eq.J) then
          HDVR(ij) = Tii  + EPot(i) - EMin
         Else
          Angle = Pi*DFloat(i-j)/DFloat(NPts)
          CNum  = Cos(Angle)
          SDen  = Sin(Angle)
          Tij   = ((-F1)**(i-j))*CNum/(F2*Den*SDen**2)
          HDVR(ij) = Tij
         EndIf
  410   Continue
  400  Continue
      Else
       Write(IOut,'(''Unsupported value for IMode in V1DVR'')')
       Stop
      EndIf
      Return
      End
*Deck DVRAnl 
      Subroutine DVRAnl(IOut,IPrint,MaxPts,NPts,NEner,NEVec,JR,JMax,
     $  VScal,IVType,XGrid,YGrid,EVal,EVec,E00,E01,We,WeXe,Be,Alpe,De)
      Implicit Real*8(A-H,O-Z)
C
C     Anharmonic 1D Treatment: Output
C     ===============================
C
C     Description
C     -----------
C     Print detailed information on the results of the 1D treatment.
C
C     Input
C     -----
C     MaxPts :: Maximum number of points considered in treatment
C     NPts   :: Number of data points actually used
C     NEner  :: Number of energy levels treated
C     NEVec  :: Number of computed eigenfunctions
C     JR     :: Current rotational state of interest
C     JMax   :: Highest rotational state to be treated
C     VScal  :: Scaling factor to apply to the eigenvectors
C     IVType :: (NPts) Type of vibrational state
C               -> See [AnD1SrtV] for details
C     XGrid  :: (NPts) Computed/fitted X data points
C     YGrid  :: (NPts) Computed/fitted Y data points
C     EVal   :: (NEner) Vibrational energies
C     EVec   :: (MaxPts,NEVec) Eigenvec. of the vibrot Hamiltonian
C
C     In-Out
C     ------
C O   E00    :: Wavenumbers of the fundamental rotational state (0,0)
C O   E01    :: Wavenumbers of the first rotational state (0,1)
C O   We     :: Vibrational constant (in cm^-1^)
C O   WeXe   :: Anharmonic correction to vibrational constant (in cm^-1^)
C O   Be     :: Rotational constant (in cm^-1^)
C O   Alpe   :: Vibro-rotational alpha matrix (in cm^-1^)
C O   De     :: Anharmonic correction to rotational constant (in cm^-1^)
C
CEND
C
C     Dimensions
      Integer MaxPts, NEner, NEVec, NPts
C     Input
      Integer IVType(*), IOut, IPrint, JR, JMax
      Real*8 EVec(MaxPts,*), EVal(*), XGrid(*), YGrid(*), VScal
C     External
      Real*8 CnvFct 
C     InOut
      Real*8 Alpe, Be, De, E00, E01, We, WeXe
C     Local
      Integer i, i1, i2, j, JCont, k, NCol, Nk
      Real*8 au2cm1, D0100, D0200, D10, D1101, D20, Eh2KCa, Eh2kJM,
     $  F1, F2, F3, F4, F12, F24, F10M6
      Character SymLb(3)*1, FmtO*80
      Save F2, F3, F4, F12, F24, F10M6
      Save SymLb
      Data F1/1.0D0/, F2/2.0D0/, F3/3.0D0/, F4/4.0D0/, F12/12.0D0/,
     $  F24/24.0D0/, F10M6/1.0D-6/
      Data SymLb/'/','+','-'/
C
 7010 Format(/,2X,76('-'),/,2X,'Vib.State(Symm)',3X,'Freq(cm^-1)',10X,
     $ 'E(a.u.)',6X,'E(Kcal/mol)',4X,'E(kJ/mol)',/,3X,'(J=',I4,')',/,2X,
     $  76('-'))
 7015 Format(5X,I3,'(',A1,')',3X,F10.3,'(',F10.3,')',2X,F10.6,
     $  2(5X,F10.3))
 7020 Format(' Vibrational Constants (cm^-1) from 0-2 levs.')
 7021 Format(' Rotational Constants (cm^-1) from J = 0,1')
 7022 Format(' Rotational Constants (cm^-1) from J = 0,1,2')
 7025 Format(' We     =',F11.5,/,' WeXe   =',F11.5)
 7026 Format(' Be     =',F11.5,/,' Alphae =',F11.5,:,' De     =',F11.5)
 7030 Format(' ** Values of Psi')
 7031 Format('(/,7X,''X'',7X,''E(kJ/mol)'',',I1,'(8X,I4))')
c 7035 Format('(2F12.5,4X,',I1,'F12.5)')
 7035 Format('(F12.5,E12.5,4X,',I1,'F12.5)')
 7900 Format(/,' WARNING: Energy levels too close. ',
     $  'Eigenvectors scaled by 1.')
C
      au2cm1 = CnvFct('au2cm1')
      Eh2kJM = CnvFct('au2kJM')
      Eh2kCa = CnvFct('au2kcaM')
      JCont  = 1
  100 If(VScal.lt.F10M6) then
        VScal = Eh2kJM*(EVal(JCont+1)-EVal(JCont))/F2
        JCont = JCont + 1
        If(JCont.lt.NEVec) then
          goto 100
        else
          Write(IOut,7900)
          VScal = F1
          endIf
        endIf
      If(JR.eq.0.or.IPrint.ne.0) then
        Write(IOut,7010) JR
        Write(IOut,7015) 0, SymLb(IVType(1)), EVal(1)*au2cm1,
     $    (EVal(1)-EVal(1))*au2cm1, EVal(1), EVal(1)*Eh2kCa,
     $    EVal(1)*Eh2kJM
        Do 200 i = 2, NEner
          Write(IOut,7015) i-1, SymLb(IVType(i)), EVal(i)*au2cm1,
     $      (EVal(i)-EVal(i-1))*au2cm1, EVal(i), EVal(i)*Eh2kCa,
     $      EVal(i)*Eh2kJM
  200     Continue
        endIf
      If(JR.eq.0) then
        D10  = (EVal(2)-EVal(1))*au2cm1
        D20  = (EVal(3)-EVal(1))*au2cm1
        WeXe = (D20-F2*D10)/F2
        We   = D10 - F2*WeXe
        E00  = EVal(1)
        If(JMax.eq.0) then
          Write(IOut,7020)
          Write(IOut,7025) We, WeXe
          endIf
      else if(JR.eq.1.and.JMax.lt.2) then
        D1101 = (EVal(2) - EVal(1))*au2cm1
        Alpe  = (We-D1101)/2 + WeXe
        D0100 = (EVal(1)-E00)*au2cm1
        Be    = (D0100+Alpe)/2
        Write(IOut,7020)
        Write(IOut,7025) We, WeXe
        Write(IOut,7021)
        Write(IOut,7026) Be, Alpe
        E01 = EVal(1)
      else if(JR.eq.1.and.JMax.ge.2) then
        D1101 = (EVal(2)-EVal(1))*au2cm1
        Alpe  = (We-D1101)/F2 + WeXe
        E01   = EVal(1)
      else if(JR.eq.2) then
        D0100 = (E01-E00)*au2cm1
        D0200 = (EVal(1)-E00)*au2cm1
        Be    = (Alpe/F2) + (F3*D0100/F4) - (D0200/F12)
        De    = (F3*D0100-D0200)/F24
        Write(IOut,7020)
        Write(IOut,7025) We, WeXe
        Write(IOut,7022)
        Write(IOut,7026) Be, Alpe, De
        endIf
      If(IPrint.gt.0) then
        Write(IOut,7030)
        NCOL = 4
        Do 300 i = 1, NEVec, NCol
          i1 = i
          i2 = Min(NEVec,i+NCol-1)
          Nk = i2-i1+1
          Write(FmtO,7031) Nk
          Write(IOut,FmtO) (k, k=i1,i2)
          Write(FmtO,7035) Nk
          Do 310 j = 1, NPts
            Write(IOut,FmtO) XGrid(j), YGrid(j)*Eh2kJM,
     $        (EVec(j,k)*VScal+EVal(k)*Eh2kJM,k=i1,i2)
c            Write(IOut,FmtO) XGrid(j), YGrid(j)*Eh2kJM, (EVec(j,k),
c     $        k=i1,i2)
  310       Continue
  300     Continue
        endIf
      Return
      End


