C     Merlino4 independent dense small-space VCI helpers.
C     Fixed-form Fortran77. Large spaces will use a Davidson driver later.

      Subroutine M4VCIBasis(NMode,MaxQ,MaxSt,NState,Basis,Info)
      Integer NMode,MaxQ,MaxSt,NState,Basis(MaxSt,NMode),Info
      Integer State(64),I,K,Sum,Done
      If(NMode.gt.64) Then
         Info=2
         Return
      End If
      NState=0
      Info=0
      Do 10 I=1,NMode
         State(I)=0
10    Continue
      Done=0
20    Continue
      Sum=0
      Do 30 I=1,NMode
         Sum=Sum+State(I)
30    Continue
      If(Sum.le.MaxQ) Then
         NState=NState+1
         If(NState.gt.MaxSt) Then
            Info=1
            Return
         End If
         Do 40 I=1,NMode
            Basis(NState,I)=State(I)
40    Continue
      End If
      K=NMode
50    Continue
      If(K.le.0) Then
         Done=1
      Else If(State(K).lt.MaxQ) Then
         State(K)=State(K)+1
         Do 60 I=K+1,NMode
            State(I)=0
60       Continue
      Else
         K=K-1
         Go To 50
      End If
      If(Done.eq.0) Go To 20
      Return
      End

      Subroutine M4VCIHarm(NMode,NState,Basis,Freq,H)
      Integer NMode,NState,Basis(NState,NMode)
      Double Precision Freq(NMode),H(NState,NState),E
      Integer I,J,K
      Do 20 I=1,NState
         Do 10 J=1,NState
            H(I,J)=0.0D0
10       Continue
20    Continue
      Do 40 I=1,NState
         E=0.0D0
         Do 30 K=1,NMode
            E=E+Freq(K)*(Dble(Basis(I,K))+0.5D0)
30       Continue
         H(I,I)=E
40    Continue
      Return
      End

      Double Precision Function M4X4Diag(N)
      Integer N
C     <n|x**4|n> for x=(a+a+)/sqrt(2)
      M4X4Diag=0.75D0*(2.0D0*Dble(N)*Dble(N)+
     $         2.0D0*Dble(N)+1.0D0)
      Return
      End
