subroutine readinput

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!  This routine reads variables from fort.8
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!use pks
use system
use solver
!use kai

implicit none
integer i
real*8 xmNaalphaini
character*30 basura1
character*30 basura2
character*30 basura3
character*30 basura4
character*30 basura5
character*30 basura6
character*30 basura7
character*30 basura8
character*30 basura9
character*30 basura10
character*30 basura11
character*30 basura12

! read starts here, not that read is performed sequentially! 

read(8,*), basura1
read(8,*), Mai, Maf, Map    ! Ma  polA

read(8, *), basura2
read(8, *), Mbi, Mbf, Mbp  !   MEO  polEO

read(8, *), basura3! vp
read(8, *), vpolcero !

read(8, *), basura4! rsal
read(8, *), rsal !

read(8, *), basura5
read(8, *), logxmAalphai, logxmAalphaf, npasosxmAalpha  !  scan total monomer density in molar

read(8, *), basura6
read(8, *), logxmNaalphai, logxmNaalphaf, npasosxmNaalpha   ! scan ratio monomer density
!
read(8, *), basura7
read(8, *), pKDi, pKDf, pKDp     ! polymer-polymer attraction strenght in kBT
!
read(8, *), basura8
read(8, *), pKAi, pKAf, pKAp      ! polymer-polymer attraction strenght in kBT
!
read(8, *), basura9
read(8, *), pkBi, pKBf, pKBp     ! polymer-polymer attraction strenght in kBT

read(8, *), basura10
read(8, *), chi  ! cutoff for porr sv interaction in lattice sites

read(8,*), basura11
read(8,*), justone ! solves only one point

read(8, *), basura12
read(8, *), xmNaalphaini,xmNaalphainitial,xmBalphainitial, xmClalphainitial, &
            xmAbetainitial, xmNabetainitial,xmBbetainitial,xmClbetainitial ! read initial guess in the same order as output, ratio and xmtot in alpha are not solved for, so the result is stores in basura

! SAVE INVERTED
! read starts here, not that read is performed sequentially! 

write(9,*), basura1
write(9,*), Mai, Maf, Map    ! Ma  polA

write(9, *), basura2
write(9, *), Mbi, Mbf, Mbp  !   MEO  polEO

write(9, *), basura3! vp
write(9, *), vpolcero !

write(9, *), basura4! rsal
write(9, *), rsal !

write(9, *), basura5
write(9, *), log10(xmBalphainitial), logxmAalphaf, npasosxmAalpha  !  scan total monomer density in molar

write(9, *), basura6
write(9, *), logxmNaalphai, logxmNaalphaf, npasosxmNaalpha   ! scan ratio monomer density
!
write(9, *), basura7
write(9, *), pKDi, pKDf, pKDp     ! polymer-polymer attraction strenght in kBT
!
write(9, *), basura8
write(9, *), pKAi, pKAf, pKAp      ! polymer-polymer attraction strenght in kBT
!
write(9, *), basura9
write(9, *), pkBi, pKBf, pKBp     ! polymer-polymer attraction strenght in kBT

write(9, *), basura10
write(9, *), chi  ! cutoff for porr sv interaction in lattice sites

write(9,*), basura11
write(9,*), justone ! solves only one point

write(9, *), basura12
write(9, *), xmBalphainitial, xmClalphainitial, xmNaalphaini,xmNaalphainitial &
             ,xmBbetainitial,xmClbetainitial, xmAbetainitial, xmNabetainitial ! read initial guess in the same order as output, ratio and xmtot in alpha are not solved for, so the result is stores in basura
close(9)

end subroutine
