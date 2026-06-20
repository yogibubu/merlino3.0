*Deck DNICFq
      Subroutine DNICFq(IOut,IPrint,NVar,FMat,GMat,NImag,EVal,EVec,V,
     $  MDV)
      Implicit Real*8(A-H,O-Z)
c
C     DiNaUtility - Internal Coordinates: Harmonic Frequency calculations
C     ===================================================================
C
C     Description
C     -----------
C     Given as input the Wilson G matrix and the Hessian matrix in internal
C     coordinates (F matrix) and compute the normal modes in internal coordinates
C     with the GF analysis by using the Miyazawa method (JCP 29,246(1958))
C
C     Input
C     -----
C     NVar :: integer
C         Number of internal coordinates.
C         NOTE: `NAt3` for complete non-redundant internal coordinates.
C     FMat :: real*8, dimension(NAt3,NAt3)
C         Hessian (F) matrix in internal coordinates.
C     GMat :: real*8, dimension(NAt3TT)
C         Wilson G matrix. Stored as lower triangular.
C
C     Output
C     ------
C     NImag :: integer
C         Number of imaginary frequencies.
C     EVal :: real*8, dimension(NVar)
C         Eigenvalues of the GF matrix (in cm^-1^).
C     EVec :: real*8, dimension(NVar,NVar)
C         Eigenvectors of the GF matrix.
C
CEND
C
C     Dimensions
      Integer MDV, NVar
      Common /PhyCon/ PhyCon(30)
C     Input
      Integer IOut, IPrint
      Real*8 FMat(*), GMat(NVar,*), DinFac, PhyCon
C     Output
      Real*8 EVal(NVar), EVec(NVar,*)
      Integer NImag
C     Local
      Real*8 V(*), amu2au, ToCm1, X, One, Ten, Zero
      Integer i
C     Indexes
      Integer IEnd, INorm, IScr, IScr2, IScr3, IScrM, IW
C     Parameters
      Save Zero, One, Ten
      Data Zero/0.0D0/, One/1.0D0/, Ten/1.0D1/
C
C
C     Initialization
C     ==============
C
C     Parameters
C     ----------
      ToCm1  = DinFac('au2cm1')**2
      amu2au = PhyCon(2)/PhyCon(12)
C
C     Memory Allocation
C     -----------------
      IW    = 1
      IScrM = IW    + NVar**2
      IScr  = IScrM + NVar**2
      IScr2 = IScr  + NVar**2
      IScr3 = IScr2 + NVar**2
      INorm = IScr3 + NVar**2
      IEnd  = INorm + NVar
      Call TstCor(IEnd,MDV,'FreqIn init')
      Call AClear(IEnd-1,V)
C
C     Debugging: Input Data
C     ---------------------
      If(IPrint.gt.1) then
        Call OutMtS(IOut,'G Matrix',0,1,GMat,NVar,NVar,NVar,NVar)
        Call LTOutS(IOut,'F Matrix',0,NVar,FMat,1)
        endIf
C
C
C     GF Analysis
C     ===========
C     - calculate the square root of the Wilson G matrix
C     - diagonalize the (G^1/2*F*G^1/2) matrix (using specific
C       subroutines for the symmetric matrices)
C     - scale back the eigenvectors by the G^-1/2 matrix
      Call AMove(NVar**2,GMat,V(IW))
      Call RootMt(V(IW),V(IScrM),V(IScr),V(IScr2),NVar,NVar,0)
C
      Call Square(FMat,V(IScr2),NVar,NVar,0)
      Call XGEMM(1,'N','N',NVar,NVar,NVar,One,V(IW),NVar,V(IScr2),NVar,
     $  Zero,V(IScr),NVar)
      Call XGEMM(1,'N','N',NVar,NVar,NVar,One,V(IScr),NVar,V(IW),NVar,
     $  Zero,V(IScr2),NVar)
      Call XSyLT1(.false.,NVar,NVar,V(IScr2),V(IScr))
      Call DiagDN(IOut,IPrint,1,V(IScr),V(IScrM),EVal,NVar,V(IScr2),
     $  NVar,V(IScr3),NVar,.false.)
C     Scale back the eigenvectors
      Call AMove(NVar**2,GMat,V(IW))
      Call RootMt(V(IW),V(IScr3),V(IScr),V(IScr2),NVar,NVar,0)
      Call XGEMM(1,'N','N',NVar,NVar,NVar,One,V(IW),NVar,V(IScrM),NVar,
     $  Zero,V(IScr2),NVar)
C     Scale the frequencies to cm^-1 (and divide to take into account
C     the units of the mass). Scale back also the normal modes.
      NImag = 0
      Do 100 i = 1, NVar
        X = Sqrt(Abs(EVal(i))*ToCm1/amu2au)
        If(EVal(i).lt.Zero) X = -X
        If(X.lt.-Ten) NImag = NImag + 1
        EVal(i) = X
  100   Continue
      Call AScale(NVar**2,One/sqrt(amu2au),V(IScr2),V(IScr2))
C
C
C     PED Analysis
C     ============
      If(IPrint.ge.1) then
        Call BldPED(IOut,NVar,EVal,V(IScr2),V(IScrM),V(IEnd),MDV-IEnd+1)
        endIf
C
C
C     Save EVec
C     =========
      Call AMove(NVar**2,V(IScr2),EVec)
      Return
      End
*Deck DNICGF
      Subroutine DNICGF(IOut,IPrint,NAtoms,NVib,AtMass,Crd,FX,FFX,BMat,
     $  DBMat,GMat,FMat,V,MDV)
      Implicit Real*8(A-H,O-Z)
C
C     DiNaUtility - Internal Coordinates: G/F Matrices Builder
C     ========================================================
C
C     Description
C     -----------
C     Given the B matrix, builds the F and G matrix in non-redundant
C       internal coordinates. The dimension of the matrices is NAt3*NAt3
C       since the six translation and rotation vectors are included.
C     Please note that the contributions due to the residual cartesian
C       gradient are eliminated.
C
C     Input
C     -----
C     AtMass :: real*8, dimension(NAtoms)
C         Atomic masses
C     Crd :: real*8, dimensions(3,NAtoms)
C         Atomic Cartesian coordinates
C     FX :: real*8, dimension(NAt3)
C         Cartesian forces
C     FFX :: real*8, dimension(NAt3TT)
C         Cartesian force constants in lower triangular format
C     BMat :: real*8, dimension(NVib,NAt3)
C         Non-redundant Wilson B matrix
C     DBMat :: real*8, dimension(NVib,NAt3,NAt3)
C         Non-redundant Wilson B matrix deriv.
C
C     Output
C     ------
C     GMat :: real*8, dimensions(NAt3,NAt3)
C         Wilson G matrix, including translations/rotations.
C     FMat :: real*8, dimensions(NAt3,NAt3)
C         Hessian (F) matrix in internal coordinates.
C
CEND
C
C     Dimensions
      Integer MDV, NAt3, NAt3Sq, NAtoms, NVib
C     External
      Integer InToWP
C     Input
      Integer IOut, IPrint
      Real*8 DBMat(NVib,3*NAtoms,*), BMat(NVib,*), Crd(3,*), AtMass(*),
     $  FX(*), FFX(*)
C     Output
      Real*8 GMat(3*NAtoms,*), FMat(*)
C     Local
      Integer i, ii, IOff1, IOff2, IOff3, Info, j, jj, k
      Real*8 V(*), One, X, Zero
C     Indexes
      Integer IBInv, IBMat, IGInt, IEnd, IIPiv, IProj, IScr, IScr2,
     $  IScr3, IScrG
      Save One, Zero
      Data One/1.0d0/, Zero/0.0d0/
C
 9000 Format(' ERROR: Null diagonal element from the LU decomposition ',
     $  'of B:',I5)
 9001 Format(' ERROR: Internal error in LU decomposition of B (arg. ',
     $  I2,')')
C
C     ================
C      INITIALIZATION
C     ================
      NAt3   = 3*NAtoms
      NAt3Sq = NAt3**2
C     =========
C      STORAGE
C     =========
      IBMat = 1
      IBInv = IBMat + NAt3Sq
      IScr  = IBInv + NAt3Sq
      IScr2 = IScr  + NAt3Sq
      IScr3 = IScr2 + NAt3Sq
      IScrG = IScr3 + NAt3Sq
      IGInt = IScrG + NAt3Sq
      IProj = IGInt + NAt3
      IIPiv = IProj + NAt3Sq
      IEnd  = IIPiv + InToWp(NAt3)
      Call TstCor(IEnd,MDV,'IntFG - init')
      Call AClear(IEnd,V)
C     =============================================
C      BUILD SQUARE B MATRIX WITH TRANSL. AND ROT.
C     =============================================
      Do 100 i = 1, NVib
        Do 110 j = 1, NAt3
  110     V(IBMat+NAt3*(j-1)+i-1) = BMat(i,j)
  100   Continue
      Call DNTRBM(NAtoms,.False.,AtMass,Crd,V(IScr))
C
      Do 200 i = 1, NAtoms
        ii = 3*(i-1)
        X  = AtMass(i)
        Do 210 j = 1, 6
          IOff1 = NAt3*(j-1) + ii
          IOff2 = NAt3*(NVib+j-1) + ii
          V(IScr2+IOff1)   = V(IScr+IOff1)*X
          V(IScr2+IOff1+1) = V(IScr+IOff1+1)*X
          V(IScr2+IOff1+2) = V(IScr+IOff1+2)*X
          V(IBInv+IOff2)   = V(IScr+IOff1)
          V(IBInv+IOff2+1) = V(IScr+IOff1+1)
          V(IBInv+IOff2+2) = V(IScr+IOff1+2)
  210     Continue
  200   Continue
      Call XGEMM(1,'T','N',6,6,NAt3,One,V(IScr2),NAt3,V(IScr),NAt3,
     $  Zero,V(IScr3),6)
      Call GGetRF(6,6,V(IScr3),6,V(IIPiv),Info)
      Call GGetRI(6,V(IScr3),6,V(IIPiv),V(IScr),NAt3**2,Info)
      Call XGEMM(1,'N','T',6,NAt3,6,One,V(IScr3),6,V(IScr2),NAt3,Zero,
     $  V(IScr),6)
      Do 300 i = 1, 6
        Do 310 j = 1, NAt3
  310     V(IBMat+NAt3*(j-1)+NVib+i-1) = V(IScr+6*(j-1)+i-1)
  300   Continue
C     ================================================
C      BUILDS THE GENERALIZED INVERSE OF B (NVib,NAt3)
C     ================================================
      Do 400 i = 1, NAtoms
        ii = 3*(i-1)
        X  = AtMass(i)
        Do 410 j = 1, NVib
          IOff1 = NAt3*(j-1) + ii
          V(IScr+IOff1)   = BMat(j,ii+1)/X
          V(IScr+IOff1+1) = BMat(j,ii+2)/X
          V(IScr+IOff1+2) = BMat(j,ii+3)/X
  410     Continue
  400   Continue
      Call XGEMM(1,'N','N',NVib,NVib,NAt3,One,BMat,NVib,V(IScr),NAt3,
     $  Zero,V(IScr2),NVib)
      Call GGetRF(NVib,NVib,V(IScr2),NVib,V(IIPiv),Info)
      If(IPrint.gt.5) then
        Call HedPrt(IOut,0,
     $    'Diagonal Elements of the LU decomposition of B',0)
        Write(IOut,'(5D15.6)') (V(ISCr2+NVib*i+i),i=0,NVib-1)
        endIf
      Call GGetRI(NVib,V(IScr2),NVib,V(IIPiv),V(IScr3),NVib**2,Info)
      If(Info.gt.0) then
        Write(IOut,9000) Info
        Call Lnk1E(0)
      else if(Info.lt.0) then
        Write(IOut,9001) -Info
        Call Lnk1E(0)
        endIf
      Call XGEMM(1,'T','N',NAt3,NVib,NVib,One,BMat,NVib,V(IScr2),NVib,
     $  Zero,V(IScr),NAt3)
      Do 500 i = 1, NAtoms
        ii = 3*(i-1)
        X  = AtMass(i)
        Do 510 j = 1, NVib
          IOff1 = NAt3*(j-1) + ii
          V(IScr+IOff1)   = V(IScr+IOff1)/X
          V(IScr+IOff1+1) = V(IScr+IOff1+1)/X
          V(IScr+IOff1+2) = V(IScr+IOff1+2)/X
  510     Continue
  500   Continue
      Call AMove(NAt3*NVib,V(IScr),V(IBInv))
C     ===================================================
C      BUILDS PROJECTOR ON THE PURE VIBRATIONAL SUBSPACE
C     ===================================================
C     The pseudoinverse of B is used to build the projector onto the
C     pure vibrational subspace (ref. JCP, 132, 184101)
      Call XGEMM(1,'N','T',NVib,NVib,NAt3,One,BMat,NVib,BMat,NVib,Zero,
     $  V(IScr),NVib)
      Call GGetRF(NVib,NVib,V(IScr),NVib,V(IIPiv),Info)
      Call GGetRI(NVib,V(IScr),NVib,V(IIPiv),V(IScr2),NVib**2,Info)
      If(Info.gt.0) then
        Write(IOut,9000) Info
        Call Lnk1E(0)
      else if(Info.lt.0) then
        Write(IOut,9001) -Info
        Call Lnk1E(0)
       endIf
      Call XGEMM(1,'T','N',NAt3,NVib,NVib,One,BMat,NVib,V(IScr),NVib,
     $ Zero,V(IScr2),NAt3)
      Call XGEMM(1,'N','N',NAt3,NAt3,NVib,One,V(IScr2),NAt3,BMat,NVib,
     $  Zero,V(IProj),NAt3)
C     ================
C      BUILD F MATRIX
C     ================
C     1) Put F in square form and compute the Hessian in cartesian
C     coordinates (projecting out contributions due to tr. and rot.)
      Call Square(FFX,V(IScr),NAt3,NAt3,0)
      If(IPrint.ge.3)
     $  Call OutMtS(IOut,'Cartesian Force Constants matrix',0,1,V(IScr),
     $    NAt3,NAt3,NAt3,NAt3)
C     Compute the gradient in internal coordinates
      Call GGEMV('T',NAt3,NAt3,One,V(IProj),NAt3,FX,1,Zero,V(IScr2),1)
      Call GGEMV('T',NVib,NAt3,One,V(IBInv),NAt3,V(IScr2),1,Zero,
     $  V(IGInt),1)
      Call AClear(NAt3**2,V(IScrG))
      Do 600 i = 1, NAt3
        IOff1 = NAt3*(i-1)
        Do 610 j = 1, NAt3
          IOff2 = IOff1 + j - 1
          X = Zero
          Do 620 k = 1, NVib
  620       X = X + V(IGInt+k-1)*DBMat(k,i,j)
          V(IScrG+IOff2) = V(IScrG+IOff2) + X
  610     Continue
  600   Continue
      If(IPrint.ge.6)
     $  Call OutMtS(IOut,'Gradient contrib. to F matrix',0,1,V(IScrG),
     $    NAt3,NAt3,NAt3,NAt3)
      Call ASub(NAt3**2,V(IScr),V(IScrG),V(IScr))
      Call XGEMM(1,'T','N',NAt3,NAt3,NAt3,One,V(IProj),NAt3,V(IScr),
     $  NAt3,Zero,V(IScr2),NAt3)
      Call XGEMM(1,'N','N',NAt3,NAt3,NAt3,One,V(IScr2),NAt3,V(IProj),
     $  NAt3,Zero,V(IScr),NAt3)
      Call XGEMM(1,'T','N',NAt3,NAt3,NAt3,One,V(IBInv),NAt3,V(IScr),
     $  NAt3,Zero,V(IScr2),NAt3)
      Call XGEMM(1,'N','N',NAt3,NAt3,NAt3,One,V(IScr2),NAt3,V(IBInv),
     $  NAt3,Zero,V(IScr),NAt3)
      Call XSyLT1(.false.,NAt3,NAt3,V(IScr),FMat)
C     ================
C      BUILD G MATRIX
C     ================
      Do 700 j = 1, NAtoms
        jj = 3*(j-1)
        IOff1 = NAt3*jj
        IOff2 = NAt3*(jj+1)
        IOff3 = NAt3*(jj+2)
        X = AtMass(j)
        Do 710 i = 1, NAt3
          V(IScr+IOff1+i-1) = V(IBMat+IOff1+i-1)/X
          V(IScr+IOff2+i-1) = V(IBMat+IOff2+i-1)/X
          V(IScr+IOff3+i-1) = V(IBMat+IOff3+i-1)/X
  710     Continue
  700   Continue
      Call XGEMM(1,'N','T',NAt3,NAt3,NAt3,One,V(IScr),NAt3,V(IBMat),
     $  NAt3,Zero,GMat,NAt3)
      Return
      End
