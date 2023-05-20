subroutine fkfun(x,f,ier)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! User provided routine for kinsol
! x is the input vector
! f is the output vector, kinsol will change x in order to get f = 0
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!   xmNaalphatot = xmaddedNaCl + xmAalpha*MA ! Na+ = Na+_added + PolA * MA, total (including polymer assoc)


use solver
use system
use const
use results
implicit none
 
real*8 xmsalt
real*8 Penality,testKa,testkeo,testkd
real*8 testneuta,testneutb,testpcka,testpckb
integer*4 ier
real*16 vectfalpha(4),vectfbeta(4),vectfrac(6)
real*8 x(7),f(7)
real*8 potA,potB,potNa,free_ener
real*8 muAalpha,muAbeta,muBalpha,muBbeta,fealpha,febeta
real*8 potquimA,elib,potquimB,potquimNa,muNaalpha,muNabeta
real*8 neutralalpha,neutralbeta,elecneualpha,elecneubeta
real*8 penalityA,penalityNa,penalityB
real*8 diffNaalpha
integer i

xmNaalpha=exp(x(7))     ! xmNa en alfa
xmBalpha =xmAalpha/ratioalpha 
xmClalpha=exp(x(2))

xmAbeta=exp(x(3))        !xmA en beta
xmNabeta=exp(x(4))        !xmNa en beta
xmBbeta=exp(x(5))
xmClbeta=exp(x(6))

print*,'xm',xmAalpha,xmnaalpha,xmBalpha,xmClalpha,xmAbeta,xmNabeta,xmBbeta,xmClbeta
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FRACTIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

vectfalpha(1)=xmAalpha
vectfalpha(2)=xmBalpha
vectfalpha(3)=xmNaalpha
vectfalpha(4)=xmClalpha

vectfbeta(1)=xmAbeta
vectfbeta(2)=xmBbeta
vectfbeta(3)=xmNabeta
vectfbeta(4)=xmClbeta
call fractions(vectfalpha,vectfrac)

!vectfractions-> (1) fEO_aspol,(2)fA_aspol,(3) fEO_asion (4)fA_Asion,(5)fEO_unas,f(6)fA_unas
fB_aspol_alpha=vectfrac(1)
fA_aspol_alpha=vectfrac(2)

fB_asion_alpha=vectfrac(3)
fA_asion_alpha=vectfrac(4)

fB_unas_alpha=vectfrac(5)
fA_unas_alpha=vectfrac(6)

print*, 'ALFA', 'fB_aspol:', vectfrac(1), 'fA_aspol:', vectfrac(2)
print*, 'ALFA', 'fB_asion:', vectfrac(3), 'fA_asion:', vectfrac(4)
print*, 'ALFA', 'fB_unas:', vectfrac(5), 'fA_unas:', vectfrac(6)

!stop

print*,'ALFA, test', fB_aspol_alpha*xmBalpha*Mb, fA_aspol_alpha*xmAalpha*Ma

!stop

!testeo fracciones
!testKa=-log10(vsol*(Na/1.0d24)*xmNaalpha*vsol*(1.-fA_asion_alpha-fA_aspol_alpha)/fA_asion_alpha )-pKa
!testkeo=-log10(vsol*(Na/1.0d24)*xmNaalpha*vsol*(1.-fEO_asion_alpha-fEO_aspol_alpha)/fEO_asion_alpha )-pKEO
!testkd=-log10((Na/1.0d24)*xmnaalpha*vsol*fA_unas_alpha*Ma*xmAalpha*vab*(1.-fEo_asion_alpha-fEO_aspol_alpha)/fEO_aspol_alpha) -pKd
!print*,'testalpha',testKa,testkeo,testkd
!stop

call fractions(vectfbeta,vectfrac)

!vectfractions-> (1) fEO_aspol,(2)fA_aspol,(3) fEO_asion (4)fA_Asion,(5)fEO_unas,f(6)fA_unas

fB_aspol_beta=vectfrac(1)
fA_aspol_beta=vectfrac(2)
fB_asion_beta=vectfrac(3)
fA_asion_beta=vectfrac(4)
fB_unas_beta=vectfrac(5)
fA_unas_beta=vectfrac(6)


print*, 'BETA', 'fB_aspol:', vectfrac(1), 'fA_aspol:', vectfrac(2)
print*, 'BETA', 'fB_asion:', vectfrac(3), 'fA_asion:', vectfrac(4)
print*, 'BETA', 'fB_unas:', vectfrac(5), 'fA_unas:', vectfrac(6)

!stop

!testeo fracciones
!testKa=-log10(vsol*(Na/1.0d24)*xmNaalpha*vsol*(1.-fA_asion_alpha-fA_aspol_alpha)/fA_asion_alpha )-pKa
!testkeo=-log10(vsol*(Na/1.0d24)*xmNaalpha*vsol*(1.-feO_asion_alpha-fEo_aspol_alpha)/feO_asion_alpha )-pKEO
!testkd=-log10((Na/1.0d24)*xmNabeta*vsol*fA_unas_beta*Ma*xmAbeta*vab*(1-fEo_asion_beta-fEO_aspol_beta)/fEO_aspol_beta) -pKd
!print*,'testbeta',testKa,testkeo,testkd
!stop
!print*,'beta',vectfrac
!!!!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  AUXILIARY CALC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!xmClalpha=-fA_unas_alpha*Ma*xmAalpha + fB_asion_alpha*Meo*xmEOalpha + xmNaalpha  ! xmCl en alpha
!xmClbeta=-fA_unas_beta*Ma*xmAbeta + fB_asion_beta*MB*xmEObeta + xmNabeta  ! xmCl en alpha


xSolventalpha=1. -Ma*vpol*vsol*xmAalpha -Mb*vpol*vsol*xmBalpha&
-fa_asion_alpha*Ma*xmAalpha*vpos*vsol-fB_asion_alpha*Mb*xmBalpha*vneg*vsol&
-xmNaalpha*vpos*vsol -xmClalpha*vneg*vsol

xSolventbeta=1. - Ma*vpol*vsol*xmAbeta  -Mb*vpol*vsol*xmBbeta&
-fA_asion_beta*Ma*xmAbeta*vpos*vsol-fB_asion_beta*Mb*xmBbeta*vneg*vsol&
-xmNabeta*vpos*vsol -xmClbeta*vneg*vsol

xmSolventalpha=xSolventalpha/vsol
xmSolventbeta =xSolventbeta/vsol

packconst=(1./vsol)*(log(xSolventalpha)-log(xSolventbeta) ) ! betapi 
neutralconst=log(xmClbeta*vsol)+vneg*vsol*packconst-log(xmClalpha*vsol) ! phi

print*,'Solven_alpha, Solven_beta,packconst,neutralconst',xSolventalpha,xSolventbeta,packconst,neutralconst
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111
!  CHEMICAL POTS AND FREE ENERGY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


potquimNa=0.
call muNa(potquimNa)
potNa=potquimNa

potquimA=0.
call muA(potquimA)
potA = potquimA

potquimB=0.
call muB(potquimB)
potB = potquimB

elib=0.
call fe(elib)
free_ener=elib

elecneualpha=0.
call electroneutroalpha(elecneualpha)
neutralalpha=elecneualpha

elecneubeta=0.
call electroneutrobeta(elecneubeta)
neutralbeta=elecneubeta


xmaddedNaCl = xmNaalpha + fA_asion_alpha*MA*xmAalpha - xmAalpha*MA

Penality=abs(xmNabeta-xmNaalpha)/(xmNabeta*0.5+xmNaalpha*0.5)
Penality=Penality+abs(xmAalpha-xmAbeta)/(xmAalpha*0.5+xmAbeta*0.5)
Penality=Penality+abs(xmBalpha-xmBbeta)/(xmBalpha*0.5+xmBbeta*0.5)
Penality=Penality+abs(xmClalpha-xmClbeta)/(xmClalpha*0.5+xmClbeta*0.5)

!Penality=Penality+abs(ratioBAalpha-ratioEOAbeta)/(ratioeoaalpha*0.5+ratioeoabeta*0.5)

f(1)=-free_ener/Penality
f(2)=-potNa/Penality
f(3)=-potA/Penality
f(4)=-potB/Penality
f(5)=-neutralalpha/Penality
f(6)=-neutralbeta/Penality
f(7)=exp(x(1))-xmBalpha

iter = iter + 1
norma = 0.0

do i = 1, 7
  print*,'f',i,f(i)
  norma = norma +(f(i))**2    
enddo

!testneuta=xmNaalpha -xmClalpha-Ma*xmAalpha*fA_unas_alpha+fEO_asion_alpha*Meo*xmEOalpha
!testneutb=xmNabeta -xmClbeta-Ma*xmAbeta*fA_unas_beta+fEO_asion_beta*Meo*xmEObeta

!testpcka=1.-Ma*xmAalpha*vpol*vsol-xmSolventalpha*vsol -Meo*xmEOalpha*vpol*vsol-xmNaalpha*vpos*vsol-xmClalpha*vneg*vsol
!testpcka=testpcka -vneg*vsol*(Ma*xmAalpha*(fA_asion_alpha+fA_aspol_alpha)+Meo*xmEOalpha*fEO_asion_alpha)
 
!testpckb=1.-Ma*xmAbeta*vpol*vsol-xmSolventbeta*vsol -Meo*xmEObeta*vpol*vsol-xmNabeta*vpos*vsol-xmClbeta*vneg*vsol
!testpckb=testpckb -vneg*vsol*(Ma*xmAbeta*(fA_asion_beta+fA_aspol_beta)+Meo*xmEObeta*fEO_asion_beta)
!print*,'testneutra',testneuta,testneutb,testpcka,testpckb

print*,norma,penality !if (Norma.lt.1E-5)then

ier = 0.0
return
end subroutine
