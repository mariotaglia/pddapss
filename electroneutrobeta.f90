
subroutine electroneutrobeta(electrob)
use const
use system
use results
implicit none
real*8 electrob
!!PSSPDDA
electrob= -fA_unas_beta*Ma*xmAbeta +fB_unas_beta*Mb*xmBbeta+xmNabeta-xmClbeta

end  subroutine 
