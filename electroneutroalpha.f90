
subroutine electroneutroalpha(electroa)
use const
use system
use results
implicit none
real*8 electroa
!!PSSPDDA
electroa= -fA_unas_alpha*Ma*xmAalpha +fB_unas_alpha*Mb*xmBalpha+xmNaalpha-xmClalpha

end  subroutine 
