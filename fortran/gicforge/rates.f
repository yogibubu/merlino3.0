*Deck Ratesd
      Subroutine Rates(IOut,IPrint,DoColl,Dogor,DoSpin) 
      Implicit Real*8 (A-H,O-Z)
      Common/PhyCon/PhyCon(30)
      Dimension TotWTF(2),Rvdw(2)
      Logical DoColl,DoGor,DoSpin
 2600 Format(/,' Gorin Reaction Rate')
 2610 Format(' T =',F7.2,5X,'KGorin =',D12.5,' cm^3/sec')
 2700 Format(/,' Spin-Forbidden Reaction Rate')
 2710 Format(' En.cross.point(cm-1) =',F10.2,5X,'Spin-Orbit Coupling',
     $ '(cm-1) =',F10.2,/,' Gradient diff.(Ha/Ang)',D10.5,5X,
     $ 'Red.Mass(uma)',F12.5)
 2720 Format(' T(K) =',F10.2,5X,/,' Part.Func.MECP.Compl. =',
     $  D12.5,5X,'Prod.Part.Func.React. =',D12.5)
 2730 Format(' Landau-Zener Prob.    =',D12.5,5X,
     $  'Rate Const.(cm^3/sec) =',D12.5)
      T=2.9815D2
C     Collision rate
      If(DoColl) then
       write(IOut,'(/,'' Collision Reaction Rate'')')
       TotWtF(1)=3.20D1
       TotWtF(2)=5.0D1
       Rvdw(1)=1.8D0
       Rvdw(2)=2.45D0
       ToKg   = PhyCon(2)
       Avog   = PhyCon(5)
       Boltz  = PhyCon(10)
       Call Collis(IOut,IPrint,ToKg,Avog,Boltz,T,TotWtF(1),TotWtF(2),
     $   Rvdw(1),Rvdw(2),SKColT)
       Write(IOut,'('' T ='',F7.2,5X,''Collision ='',D12.5,
     $   '' cm^3/sec'')') T,SKColT
      EndIf
C Gorin Rate
      If(DoGor) then
C FMA and FMB are the masses of the two fragments
C C6 is the r-6 factor
       FMA=1.0d0
       FMB=1.0d0
       C6=1.0d0 
       write(IOut,2600)
       call Gorin (IOut,PhyCon,T,FMA,FMB,C6,GorKT)
       write(IOut,2610) T,GorKT
       EndIf
C Spin forbidden rate
      If(DoSpin) then
       write(IOut,2700)
       UmaRdM=1.331D1
       HSOcm1=4.48D1
       DGHA1=1.318D-1
       EMecm1=1.419D3
       Qmecp=1.19406D15
       QaQb=1.15586D16
       Call NAdTST(IOut,IPrint,T,UmaRdM,HSOcm1,DGHA1,Emecm1,Qmecp,QaQb,
     $   PLZ,sksfcm)
       write(IOut,2710) Emecm1,HSOcm1,DGHA1,UmaRdM
       write(IOut,2720) T,Qmecp,QaQb
       write(IOut,2730) PLZ,sksfcm
      EndIf
      Return
      End
*Deck Collis
      Subroutine Collis(IOut,IPrint,ToKg,Avog,Boltz,TUser,FMa,FMb,Ra,Rb,
     $  SKColT)
      Implicit Real*8 (A-H,O-Z)
C
C     Contrary to standard rA and rB are entered in Angstrom
C
      Save Pt5,One,Four,Eight,TStd
      Data Pt5,One,Four,Eight /5.0D-1,1.0D0,4.0D0,8.0D0/
      Data TStd/2.9815D2/
  100 Format(/,' From Collis: r(Ang):',2D15.5,' mass (uma)',
     $  2D15.5)
  200 Format(' T =',F7.2,5X,'Kcollision =',D12.5,' cm^3/sec')
      if(IPrint.gt.1) write(IOut,100)RA,RB,FMA,FMB
      pi     = Four * ATan(One)
      AngCm  = 1.0D-8
      T = TUser
      IF(T.lt.1.0d-7) T=TStd
      rAcm=rA*AngCM
      FMAKg=FMA*ToKg
      rBcm=rB*AngCM
      FMBKg=FMB*ToKg
      fac1=pi*(rAcm+rBcm)**2
      fac2a=Eight*Boltz/pi
      RedMIn=(One/FMAKg+One/FMBKg)
      fac2=sqrt(fac2a*RedmIn)
      fac3=exp(pt5)
      skcoll=fac1*fac2*fac3
      sKcolT=sKcoll*sqrt(T)
      if(IPrint.gt.1) write(IOut,200)T,skcolT
      Return
      End
*Deck Gorin
      Subroutine Gorin (IOut,PhyCon,TUser,FMA,FMB,C6,GorKT)
      Implicit Real*8 (A-H,O-Z)
      logical error
      dimension PhyCon(*)
  100 Format(/,' Gorin, masses:',2D10.5,5X,'C6 =',D10.5,5X,'T =',F10.3)
  200 Format(' T =',F7.2,5X,'KGorin =',D12.5,' cm^3/sec')
      Save zero,one,two,three,six,eleven,TStd,T0C
      Data zero,one,two,three,six/0.0d0,1.0d0,2.0d0,3.0d0,6.0d0/
      Data Eleven,TStd,T0C/1.1d1,2.9815d2,2.7315d2/ 
      TOANG  = PhyCon(1)
      TOKG   = PhyCon(2)
      BOLTZ  = PhyCon(10)
      PLANCK = PhyCon(4)
      AVOG   = PhyCon(5)
      VolMol = PhyCon(13)
C
C     Compute the GAS constant,pi, conversion from Bohr to cm
C
      Gas   = Avog * Boltz
      Pi    = Four * ATan(One)
      Angcm = 1.0D-8
      Bohrcm= Angcm*ToAng
C
      If(TUser.eq.Zero) T=TStd
      redmas=FMA*FMB/(FMA+FMB)
      write(IOut,100)FMA,FMB,C6,T
      elsxth=eleven/six
      fac1=two**elsxth
      argam=two/three 
      call LogGam(argam,gln,error)
      if(error) then
       write(IOut,'('' From Gorin: error in Gamma'')')
       stop
      endif
      fac2=exp(argam)
      onthrd=one/three
      fac3=C6**onthrd
      fac4=sqrt(pi/redmas)
      onsxth=one/six
      fac5=(Boltz*T)**onsxth
      GorKT=fac1*fac2*fac3*fac4*fac5
      write(IOut,200)T,GorKT
      Return
      End
*Deck NAdTST
      Subroutine NAdTST(IOut,IPrint,TUser,AmuM,HSOcm1,DGHA1,Emecm1,
     $  Qmecp,Qab,PLZ,sksfcm)
      Implicit Real*8 (A-H,O-Z)
      Save One,Two,Three,Four,TStd
      Data One,Two,Three,Four /1.0D0,2.0D0,3.0D0,4.0D0/
      Data TStd /2.9815D2/
  100 Format(/,' From NAdTST')
  200 Format(' En.cross.point(cm-1) =',F10.2,5X,'Spin-Orbit Coupling',
     $ '(cm-1) =',F10.2,/,' Gradient diff.(Ha/Ang)',D10.5,5X,
     $ 'Red.Mass(uma)',F12.5)
  300 Format(' T(K) =',F10.2,5X,' Part.Func.MECP.Compl. =',
     $  D12.5,5X,'Prod.Part.Func.React. =',D12.5)
  400 Format(' Landau-Zener Prob.=',D12.5,5X,'Rate Const.(cm^3/sec)=',
     $  D12.5,5X,'at ',F6.2,' K')
      pi     = Four * ATan(One)
      AngCm  = 1.0D-8
      T = TUser
      IF(T.Eq.Zero) T=TStd
      RKBJK=1.38064852D-23
      cm1ToJ=1.98630D-23
      ToKg= 1.66054D-27
      hJs=6.62607004D-34
      HA1Jm1=4.359425D-8
      fatmPa=1.01325D5
      fm3cm3=1.0D6
      if(IPrint.gt.1) then
       write(IOut,100)
       write(IOut,200) Emecm1,HSOcm1,DGHA1,AmuM 
       write(IOut,300) T,Qmecp,QaQb
      endif
      RedM = AmuM*ToKg 
      vel=sqrt(Three*RKBJK*T/RedM)
      HsoJ=Hsocm1*cm1ToJ
      DGJm1=DgHA1*HA1Jm1
      twopi=two*pi
      sav=twopi*HsoJ/dsqrt(hJs)
      term=sav*sav/(vel*DgJm1)
      fac=exp(-term)
      PLZ = one-fac
      EmeJ=Emecm1*cm1toJ
      fac1=RKBJK*T*Qmecp/(hJs*Qab)
      fac2=exp(-EmeJ/(RKBJK*T))
      sksf=PLZ*fac1*fac2
      sksfcm=sksf*RKBJK*fm3cm3/fatmPa
      If(IPrint.gt.1) write(IOut,400) PLZ,sksfcm,T
      return
      end
*Deck LogGam
      Subroutine LogGam(xx,gln,error)
      Real*8 x,xx,y,tmp,ser,gln,cof(14)
      Integer j
      Logical error
      Save cof
      Data cof /57.1562356658629235,-59.5979603554754912,
     $ 14.1360979747417471,-0.491913816097620199,
     $ .339946499848118887e-4, .465236289270485756e-4,
     $ -.983744753048795646e-4, .158088703224912494e-3,
     $ -.210264441724104883e-3, .217439618115212643e-3,
     $ -.164318106536763890e-3, .844182239838527433e-4,
     $ -.261908384015814087e-4,.368991826595316234e-5/
      error=.true.
      if(xx.lt.0.0d0) return
      x=xx
      y=x
      tmp=x+5.2421875000000
      tmp=(x+0.5)*log(tmp)-tmp
      ser=0.999999999999997092
      gln=0.0D0
      Do j=1,14
      ser=ser+cof(j)/(y+j)
      enddo
      gln=tmp+log(2.5066282746310005*ser/x)
      error=.false.
      return
      end

