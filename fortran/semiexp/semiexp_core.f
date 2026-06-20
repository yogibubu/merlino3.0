      Subroutine M4SEBondB(NAtom,IAt,JAt,XYZ,BRow,Info)
C     Analytic Wilson B row for a bond distance in Angstrom.
      Integer NAtom,IAt,JAt,Info
      Double Precision XYZ(3,NAtom),BRow(3,NAtom)
      Double Precision DX,DY,DZ,R
      Integer I
      Info=0
      Do 10 I=1,NAtom
         BRow(1,I)=0.0D0
         BRow(2,I)=0.0D0
         BRow(3,I)=0.0D0
10    Continue
      DX=XYZ(1,IAt)-XYZ(1,JAt)
      DY=XYZ(2,IAt)-XYZ(2,JAt)
      DZ=XYZ(3,IAt)-XYZ(3,JAt)
      R=DSQRT(DX*DX+DY*DY+DZ*DZ)
      If(R.le.1.0D-12) Then
         Info=1
         Return
      End If
      BRow(1,IAt)= DX/R
      BRow(2,IAt)= DY/R
      BRow(3,IAt)= DZ/R
      BRow(1,JAt)=-DX/R
      BRow(2,JAt)=-DY/R
      BRow(3,JAt)=-DZ/R
      Return
      End

      Subroutine M4SEKraitchman(A0,B0,C0,A1,B1,C1,DMass,PMass,
     $                          Ref,Coord,Info)
C     Single-substitution Kraitchman coordinates in Angstrom.
C     Rotational constants are MHz.  DMass is the isotope mass increment
C     and PMass is the parent total mass, both in amu.  The substitution
C     mass is DMass*PMass/(PMass+DMass), matching the Python SE solver.
      Integer Info,I
      Double Precision A0,B0,C0,A1,B1,C1,DMass,PMass
      Double Precision Ref(3),Coord(3),Conv,Mu
      Double Precision I0(3),I1(3),DI(3),X2(3)
      Info=0
      Conv=505379.006D0
      Do 10 I=1,3
         Coord(I)=0.0D0
10    Continue
      If(A0.le.0.0D0.or.B0.le.0.0D0.or.C0.le.0.0D0) Then
         Info=1
         Return
      End If
      If(A1.le.0.0D0.or.B1.le.0.0D0.or.C1.le.0.0D0) Then
         Info=1
         Return
      End If
      If(DMass.le.0.0D0.or.PMass.le.0.0D0) Then
         Info=2
         Return
      End If
      Mu=DMass*PMass/(PMass+DMass)
      If(Mu.le.0.0D0) Then
         Info=2
         Return
      End If
      I0(1)=Conv/A0
      I0(2)=Conv/B0
      I0(3)=Conv/C0
      I1(1)=Conv/A1
      I1(2)=Conv/B1
      I1(3)=Conv/C1
      Do 20 I=1,3
         DI(I)=I1(I)-I0(I)
20    Continue
      X2(1)=(DI(2)+DI(3)-DI(1))/(2.0D0*Mu)
      X2(2)=(DI(1)+DI(3)-DI(2))/(2.0D0*Mu)
      X2(3)=(DI(1)+DI(2)-DI(3))/(2.0D0*Mu)
      Do 30 I=1,3
         If(X2(I).gt.0.0D0) Coord(I)=DSQRT(X2(I))
         If(Ref(I).lt.0.0D0) Coord(I)=-Coord(I)
30    Continue
      Return
      End

      Subroutine M4SEClassNormalEq(NObs,NPar,NClass,ClassMap,J,Res,W,
     $                             Damp,DQ,Cov,Hess,Info)
C     Weighted least-squares with parameter classes.
C     ClassMap(I)=0 freezes parameter I. Positive values identify shared
C     active classes; all parameters in the same class receive one step.
      Integer NObs,NPar,NClass,ClassMap(NPar),Info
      Double Precision J(NObs,NPar),Res(NObs),W(NObs),Damp
      Double Precision DQ(NPar),Cov(NPar,NPar),Hess(NPar,NPar)
      Double Precision JR(500,100),DQR(100),CovR(100,100),HessR(100,100)
      Integer I,K,L,C1,C2
      Info=0
      If(NObs.gt.500.or.NPar.gt.100.or.NClass.gt.100) Then
         Info=2
         Return
      End If
      Do 20 I=1,NPar
         DQ(I)=0.0D0
         Do 10 K=1,NPar
            Cov(I,K)=0.0D0
            Hess(I,K)=0.0D0
10       Continue
20    Continue
      Do 40 K=1,NObs
         Do 30 C1=1,NClass
            JR(K,C1)=0.0D0
30       Continue
40    Continue
      Do 70 I=1,NPar
         C1=ClassMap(I)
         If(C1.gt.0.and.C1.le.NClass) Then
            Do 50 K=1,NObs
               JR(K,C1)=JR(K,C1)+J(K,I)
50          Continue
         Else If(C1.lt.0.or.C1.gt.NClass) Then
            Info=3
            Return
         End If
70    Continue
      Call M4SENormalEq(NObs,NClass,JR,Res,W,Damp,DQR,CovR,HessR,
     $                  Info)
      If(Info.ne.0) Return
      Do 100 I=1,NPar
         C1=ClassMap(I)
         If(C1.gt.0) DQ(I)=DQR(C1)
         Do 90 L=1,NPar
            C2=ClassMap(L)
            If(C1.gt.0.and.C2.gt.0) Then
               Cov(I,L)=CovR(C1,C2)
               Hess(I,L)=HessR(C1,C2)
            End If
90       Continue
100   Continue
      Return
      End

      Subroutine M4SEAngleB(NAtom,IAt,JAt,KAt,XYZ,BRow,Info)
C     Analytic Wilson B row for angle I-J-K in radians.
      Integer NAtom,IAt,JAt,KAt,Info
      Double Precision XYZ(3,NAtom),BRow(3,NAtom)
      Double Precision A(3),B(3),U(3),V(3),GI(3),GK(3)
      Double Precision RA,RB,C,S
      Integer I,M
      Info=0
      Do 10 I=1,NAtom
         Do 20 M=1,3
            BRow(M,I)=0.0D0
20       Continue
10    Continue
      RA=0.0D0
      RB=0.0D0
      Do 30 M=1,3
         A(M)=XYZ(M,IAt)-XYZ(M,JAt)
         B(M)=XYZ(M,KAt)-XYZ(M,JAt)
         RA=RA+A(M)*A(M)
         RB=RB+B(M)*B(M)
30    Continue
      RA=DSQRT(RA)
      RB=DSQRT(RB)
      If(RA.le.1.0D-12.or.RB.le.1.0D-12) Then
         Info=1
         Return
      End If
      C=0.0D0
      Do 40 M=1,3
         U(M)=A(M)/RA
         V(M)=B(M)/RB
         C=C+U(M)*V(M)
40    Continue
      If(C.gt.1.0D0) C=1.0D0
      If(C.lt.-1.0D0) C=-1.0D0
      S=DSQRT(DMAX1(1.0D0-C*C,1.0D-16))
      Do 50 M=1,3
         GI(M)=(-V(M)+C*U(M))/(S*RA)
         GK(M)=(-U(M)+C*V(M))/(S*RB)
         BRow(M,IAt)=GI(M)
         BRow(M,KAt)=GK(M)
         BRow(M,JAt)=-(GI(M)+GK(M))
50    Continue
      Return
      End

      Subroutine M4SERotConst(NAtom,Mass,XYZ,ABC,PMom,Info)
C     Principal moments (amu A**2) and rotational constants (MHz).
C     The constant matches h/(8*pi*pi*c*I) expressed for I in amu A**2.
      Integer NAtom,Info
      Double Precision Mass(NAtom),XYZ(3,NAtom),ABC(3),PMom(3)
      Double Precision COM(3),Tot,X,Y,Z,Ixx,Iyy,Izz,Ixy,Ixz,Iyz
      Double Precision A(3,3),Work(3,3),Conv
      Integer I,M
      Info=0
      Conv=505379.006D0
      Tot=0.0D0
      Do 10 M=1,3
         COM(M)=0.0D0
10    Continue
      Do 30 I=1,NAtom
         Tot=Tot+Mass(I)
         Do 20 M=1,3
            COM(M)=COM(M)+Mass(I)*XYZ(M,I)
20       Continue
30    Continue
      If(Tot.le.0.0D0) Then
         Info=1
         Return
      End If
      Do 40 M=1,3
         COM(M)=COM(M)/Tot
40    Continue
      Ixx=0.0D0
      Iyy=0.0D0
      Izz=0.0D0
      Ixy=0.0D0
      Ixz=0.0D0
      Iyz=0.0D0
      Do 50 I=1,NAtom
         X=XYZ(1,I)-COM(1)
         Y=XYZ(2,I)-COM(2)
         Z=XYZ(3,I)-COM(3)
         Ixx=Ixx+Mass(I)*(Y*Y+Z*Z)
         Iyy=Iyy+Mass(I)*(X*X+Z*Z)
         Izz=Izz+Mass(I)*(X*X+Y*Y)
         Ixy=Ixy-Mass(I)*X*Y
         Ixz=Ixz-Mass(I)*X*Z
         Iyz=Iyz-Mass(I)*Y*Z
50    Continue
      A(1,1)=Ixx
      A(1,2)=Ixy
      A(1,3)=Ixz
      A(2,1)=Ixy
      A(2,2)=Iyy
      A(2,3)=Iyz
      A(3,1)=Ixz
      A(3,2)=Iyz
      A(3,3)=Izz
      Call M4SEJacobi3(A,PMom,Work,Info)
      If(Info.ne.0) Return
      Call M4SESort3(PMom)
      Do 60 I=1,3
         If(PMom(I).le.1.0D-14) Then
            ABC(4-I)=0.0D0
         Else
            ABC(4-I)=Conv/PMom(I)
         End If
60    Continue
      Return
      End

      Subroutine M4SENormalEq(NObs,NPar,J,Res,W,Damp,DQ,Cov,Hess,
     $                        Info)
C     Weighted least-squares normal equations for semiexperimental fits.
      Integer NObs,NPar,Info
      Double Precision J(NObs,NPar),Res(NObs),W(NObs),Damp
      Double Precision DQ(NPar),Cov(NPar,NPar),Hess(NPar,NPar)
      Double Precision A(100,100),B(100),Scale,Sigma2
      Integer I,JCol,K,L,Dof
      Info=0
      If(NPar.gt.100) Then
         Info=2
         Return
      End If
      Do 20 I=1,NPar
         B(I)=0.0D0
         Do 10 JCol=1,NPar
            A(I,JCol)=0.0D0
            Hess(I,JCol)=0.0D0
10       Continue
20    Continue
      Do 50 K=1,NObs
         Scale=W(K)
         Do 40 I=1,NPar
            B(I)=B(I)+J(K,I)*Res(K)*Scale
            Do 30 JCol=1,NPar
               A(I,JCol)=A(I,JCol)+J(K,I)*J(K,JCol)*Scale
30          Continue
40       Continue
50    Continue
      Do 70 I=1,NPar
         Do 60 JCol=1,NPar
            Hess(I,JCol)=2.0D0*A(I,JCol)
            Cov(I,JCol)=A(I,JCol)
60       Continue
         A(I,I)=A(I,I)+Damp
70    Continue
      Call M4SEGauss(NPar,A,B,DQ,Info)
      If(Info.ne.0) Return
      Call M4SEInvert(NPar,Cov,Info)
      If(Info.ne.0) Return
      Sigma2=0.0D0
      Do 80 K=1,NObs
         Sigma2=Sigma2+W(K)*Res(K)*Res(K)
80    Continue
      Dof=NObs-NPar
      If(Dof.lt.1) Dof=1
      Sigma2=Sigma2/DBLE(Dof)
      Do 100 I=1,NPar
         Do 90 JCol=1,NPar
            Cov(I,JCol)=Sigma2*Cov(I,JCol)
90       Continue
100   Continue
      Return
      End

      Subroutine M4SEGauss(N,A,B,X,Info)
      Integer N,Info
      Double Precision A(100,100),B(100),X(100)
      Double Precision Pivot,Factor,Sum,Temp
      Integer I,J,K,IP
      Info=0
      Do 60 K=1,N
         IP=K
         Pivot=DABS(A(K,K))
         Do 10 I=K+1,N
            If(DABS(A(I,K)).gt.Pivot) Then
               Pivot=DABS(A(I,K))
               IP=I
            End If
10       Continue
         If(Pivot.le.1.0D-20) Then
            Info=1
            Return
         End If
         If(IP.ne.K) Then
            Do 20 J=K,N
               Temp=A(K,J)
               A(K,J)=A(IP,J)
               A(IP,J)=Temp
20          Continue
            Temp=B(K)
            B(K)=B(IP)
            B(IP)=Temp
         End If
         Do 50 I=K+1,N
            Factor=A(I,K)/A(K,K)
            A(I,K)=0.0D0
            Do 30 J=K+1,N
               A(I,J)=A(I,J)-Factor*A(K,J)
30          Continue
            B(I)=B(I)-Factor*B(K)
50       Continue
60    Continue
      Do 90 I=N,1,-1
         Sum=B(I)
         Do 70 J=I+1,N
            Sum=Sum-A(I,J)*X(J)
70       Continue
         X(I)=Sum/A(I,I)
90    Continue
      Return
      End

      Subroutine M4SEInvert(N,A,Info)
      Integer N,Info
      Double Precision A(100,100),Aug(100,200),Pivot,Factor,Temp
      Integer I,J,K,IP
      Info=0
      Do 20 I=1,N
         Do 10 J=1,N
            Aug(I,J)=A(I,J)
            Aug(I,N+J)=0.0D0
10       Continue
         Aug(I,N+I)=1.0D0
20    Continue
      Do 80 K=1,N
         IP=K
         Pivot=DABS(Aug(K,K))
         Do 30 I=K+1,N
            If(DABS(Aug(I,K)).gt.Pivot) Then
               Pivot=DABS(Aug(I,K))
               IP=I
            End If
30       Continue
         If(Pivot.le.1.0D-20) Then
            Info=1
            Return
         End If
         If(IP.ne.K) Then
            Do 40 J=1,2*N
               Temp=Aug(K,J)
               Aug(K,J)=Aug(IP,J)
               Aug(IP,J)=Temp
40          Continue
         End If
         Pivot=Aug(K,K)
         Do 50 J=1,2*N
            Aug(K,J)=Aug(K,J)/Pivot
50       Continue
         Do 70 I=1,N
            If(I.ne.K) Then
               Factor=Aug(I,K)
               Do 60 J=1,2*N
                  Aug(I,J)=Aug(I,J)-Factor*Aug(K,J)
60             Continue
            End If
70       Continue
80    Continue
      Do 100 I=1,N
         Do 90 J=1,N
            A(I,J)=Aug(I,N+J)
90       Continue
100   Continue
      Return
      End

      Subroutine M4SEJacobi3(A,D,V,Info)
      Integer Info
      Double Precision A(3,3),D(3),V(3,3)
      Double Precision App,Aqq,Apq,Phi,C,S,Aip,Aiq,Vip,Viq
      Integer I,J,Iter,P,Q
      Info=0
      Do 20 I=1,3
         Do 10 J=1,3
            V(I,J)=0.0D0
10       Continue
         V(I,I)=1.0D0
20    Continue
      Do 100 Iter=1,60
         P=1
         Q=2
         If(DABS(A(1,3)).gt.DABS(A(P,Q))) Then
            P=1
            Q=3
         End If
         If(DABS(A(2,3)).gt.DABS(A(P,Q))) Then
            P=2
            Q=3
         End If
         If(DABS(A(P,Q)).lt.1.0D-14) Goto 120
         App=A(P,P)
         Aqq=A(Q,Q)
         Apq=A(P,Q)
         Phi=0.5D0*DATAN2(2.0D0*Apq,Aqq-App)
         C=DCOS(Phi)
         S=DSIN(Phi)
         Do 40 I=1,3
            If(I.ne.P.and.I.ne.Q) Then
               Aip=A(I,P)
               Aiq=A(I,Q)
               A(I,P)=C*Aip-S*Aiq
               A(P,I)=A(I,P)
               A(I,Q)=S*Aip+C*Aiq
               A(Q,I)=A(I,Q)
            End If
40       Continue
         A(P,P)=C*C*App-2.0D0*S*C*Apq+S*S*Aqq
         A(Q,Q)=S*S*App+2.0D0*S*C*Apq+C*C*Aqq
         A(P,Q)=0.0D0
         A(Q,P)=0.0D0
         Do 50 I=1,3
            Vip=V(I,P)
            Viq=V(I,Q)
            V(I,P)=C*Vip-S*Viq
            V(I,Q)=S*Vip+C*Viq
50       Continue
100   Continue
      Info=1
      Return
120   Continue
      D(1)=A(1,1)
      D(2)=A(2,2)
      D(3)=A(3,3)
      Return
      End

      Subroutine M4SESort3(X)
      Double Precision X(3),T
      If(X(1).gt.X(2)) Then
         T=X(1)
         X(1)=X(2)
         X(2)=T
      End If
      If(X(2).gt.X(3)) Then
         T=X(2)
         X(2)=X(3)
         X(3)=T
      End If
      If(X(1).gt.X(2)) Then
         T=X(1)
         X(1)=X(2)
         X(2)=T
      End If
      Return
      End
