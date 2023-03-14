!  ###############################################################################!     
  !     Pss/PDMA Molecular Theory Program 
  !    
  !###############################################################################
!use pks
use system
use const
use solver
use results
implicit none
integer i, j, iii
real*8 logxmAalpha
real*8 logratioalpha
integer k, kk, kkk
integer kkkk, kkkkk

print*, ' PSS/PDDMA GIT Version: ', _VERSION

call readinput ! read input from file
!xmBAlphainitial=10**(xmBalphainitial)
!xmClalphainitial=
!!xmNabetainitial =
!xmAbetainitial=10**(xmAbetainitial)
!xmBbetainitial=10**(xmBbetainitial)
!xmClbetainitial=

call allocation

vab=1.
vpol=vpolcero/vsol

vneg=4./3.*pi*rsal**3/vsol !volume of anion in units of vsol
vpos=4./3.*pi*rsal**3/vsol !volume of cation in units of vsol
yes=0 ! es para  chequear si encuentra o no xalpha, xbeta


do k=1,pKDp
do kk=1,pKAp
do kkk=1,pKBp
do kkkk = 1,Map 
do kkkkk = 1,Mbp 


pKD = pKDi + (pKDf-pKDi)*float(k-1)/float(pKDp)
pKA = pKAi + (pKAf-pKAi)*float(kk-1)/float(pKAp)
pKB = pKBi + (pKBf-pKBi)*float(kkk-1)/float(pKBp)
Ma = Mai + (Maf-Mai)*float(kkkk-1)/float(Map)
Mb = Mbi + (Mbf-Mbi)*float(kkkkk-1)/float(Mbp)

print*,'pkD ,pKa, pKB, Ma, Mab',pKd ,pKa,pKB, Mb, Mb

KD=10**(-pKD)
KA=10**(-pKA)
KB=10**(-pKB)

K0A = (KA)*(vsol/Na*1.0d24) ! thermodynamic constants
K0B = (KB)*(vsol/Na*1.0d24)
K0D = (KD)*(Na/1.0d24)


do j=1, npasosratioalpha  ! loop over xmtot alpha

  logratioalpha = (logratioalphaf-logratioalphai)*float(j-1)/float(npasosratioalpha) + logratioalphai
  ratioalpha= 10**(logratioalpha) !10**  !Segunda variable que fijamos  xmpoltotalalpha


 do i = 1, npasosxmAalpha ! loop over ratio_alpha

  logxmAalpha = logxmAalphai  + (logxmAalphaf-logxmAalphai) &
  /float(npasosxmAalpha)*float(i-1)  !Na

  xmAalpha = 10**(logxmAalpha)

      iter=0
      call solve

  enddo ! j

enddo ! i

enddo !k
enddo !kk
enddo !kkk
enddo !kkkk
enddo !kkkkk


! SAVE RESULTS TO FILE

open (unit=3,file='csal_poltot_mol_alpha.txt',status='replace')

do iii=1,yes
   write (3,*) arraympoltot(1,iii), arraymcsal(1,iii)
end do

open (unit=4,file='csal_poltot_mol_beta.txt',status='replace')

do iii=1,yes
   write (4,*) arraympoltot(2,iii), arraymcsal(2,iii)
end do

open (unit=40,file='cpoltot_ratioBA_mol_alpha.txt',status='replace')

do iii=1,yes
   write (40,*) arrayratioBA(1,iii), arraympoltot(1,iii)
end do

open (unit=30,file='cpoltot_ratioBA_mol_beta.txt',status='replace')

do iii=1,yes
   write (30,*) arrayratioBA(2,iii), arraympoltot(2,iii)
end do

open (unit=400,file='csal_ratioBA_mol_alpha.txt',status='replace')

do iii=1,yes
   write (400,*) arrayratioBA(1,iii), arraymcsal(1,iii)
end do

open (unit=300,file='csal_ratioBA_mol_beta.txt',status='replace')

do iii=1,yes
   write (300,*) arrayratioBA(2,iii), arraymcsal(2,iii)
end do

open (unit=600,file='polA_polB_alpha.txt',status='replace')

do iii=1,yes
   write (600,*) arraymA(1,iii), arraymB(1,iii)
end do
 

open (unit=500,file='polA_polB_beta.txt',status='replace')

do iii=1,yes
   write (500,*) arraymA(2,iii), arraymB(2,iii)
end do

open (unit=600,file='polA_addedNa_alpha.txt',status='replace')

do iii=1,yes
   write (600,*) arraymA(1,iii), arrayaddedNaCl(iii)
end do



call endall     ! clean up and terminate
end 



subroutine endall
 stop
end subroutine
