subroutine fractions(vectf,vectfrac)
use system
use const
 implicit none
integer i
real*16 vectf(4),vectfrac(6)
real*16 xphiA,xphiEO
real*16 xmphiA,xmphiB,xmphiNa,xmphiCl
real*16 aa, bb, cc,auxA,auxB

xmphiA=vectf(1)
xmphiB=vectf(2)
xmphiNa=vectf(3) 
xmphiCl=vectf(4)

vectfrac=0.

auxA = MB*xmphiB/(Ma*xmphiA)
auxB=(K0A+xmphiNa*vsol)*(K0B+xmphiCl*vsol)
auxB=auxB/(Ma*xmphiA*vab*K0D*K0A*K0B)
!auxB = xmphiNa*vsol*Ma*xmphiA*vab/K0D
!auxB = auxB*K0A/(K0A+xmphiNa*vsol)
!auxB = AuxB*K0EO/(K0EO+xmphiNa*vsol)

aa=auxA
!bb=-1.-auxA-1./auxB
bb=-1.-auxA-auxB
cc=1.

vectfrac(1) = (-bb - SQRT(bb**2 - 4.0*aa*cc))/(2.0*auxA) ! fB_aspol
vectfrac(2) = vectfrac(1)*auxA                        ! fA_aspol

!if (abs(vectfrac(1)).gt.1. .or. abs(vectfrac(2)).gt.1.)then
!   print*,'frac',vectfrac(1),vectfrac(2),xmphiA,xmphiEo,xmphiNa
   !stop
!endif

vectfrac(3) = (1.0 - vectfrac(1))*(xmphiCl*vsol/(K0B+xmphiCl*vsol)) ! fB_asion
vectfrac(4) = (1.0 - vectfrac(2))*(xmphiNA*vsol/(K0A+xmphiNa*vsol)) ! fA_asion

!print*, xmphiCl, xmphiNA
!print*, K0B, K0A
!print*, vectfrac(1) , vectfrac(2) 
!print*, vectfrac(3) , vectfrac(4) 
!stop

!if (abs(vectfrac(3)).gt.1. .or. abs(vectfrac(4)).gt.1.)then
!   print*,'frac',vectfrac(3),vectfrac(4)
   !stop
!endif

vectfrac(5) = (1.0 - vectfrac(1)-vectfrac(3)) ! fEO_unas
vectfrac(6) = (1.0 - vectfrac(2)-vectfrac(4)) ! fA_unas

!if (abs(vectfrac(5)).gt.1. .or. abs(vectfrac(6)).gt.1.)then
  !print*,'frac',vectfrac(5),vectfrac(6)
  !stop
!endif

!print*, K0A, xmphiNa*vsol*vectfrac(6)/vectfrac(4)
!print*, K0EO, xmphiNa*vsol*vectfrac(5)/vectfrac(3)
!print*, K0D, xmphiNa*vsol*vectfrac(6)*mA*xmphiA*vab*vectfrac(5)/vectfrac(1)
!stop

end subroutine

