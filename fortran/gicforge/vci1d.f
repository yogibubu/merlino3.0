*Deck VCI1D 
      Subroutine VCI1D(In,IOut,IPrint,MxScr,V,IScr)
      Implicit Real*8 (A-H,O-Z)
C
C Driver of VCI program
C
      Dimension V(*),IScr(*)
      Dimension F(10),G(10)
      pi=4.0d0*ATan(1.0d0)
      ToDeg=1.80d+2/pi
      Thresh=1.0D-5
C Read Data
      Read(In,*) NBas,NEig
      Read(In,*) Alpha
      Read(In,*) (F(I),I=1,6)
      Read(In,*) (G(I),I=1,6)
C
      If(Alpha.lt.thresh) Alpha=1.0D0
      Write(IOut,'(//,'' *** VCI one-dimensional program ***'')')
      Write(IOUt,'(I5,'' Basis Functions for'',I3,'' Eigenstates'')')
     $   NBas,NEig
      Write(IOut,'(6X,''Alpha ='',F10.2)') Alpha
      write(IOut,'(/,'' F Coefficients (cm-1):'',6F10.2)')
     $  (F(I),I=1,6)
      write(IOut,'(  '' G Coefficients (cm-1):'',6F10.2)')
     $  (G(I),I=1,6)
      IPrint=0
      Call Var1D(IOut,IPrint,NBas,NEig,Alpha,F,G,V,IScr)
      Return 
      End 
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

