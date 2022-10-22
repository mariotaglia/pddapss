subroutine fe(elib)

use results
use system
use const
implicit none
real*8 Free_energy,elib 
!!checkear xmsolvent
!!PSSPDDA
Free_Energy = 0.0
Free_Energy=Free_Energy-xmAalpha-xmBalpha-xmSolventalpha-xmNaalpha-xmClalpha
Free_Energy=Free_Energy + 0.5*chi*(xmAalpha*Ma+xmBalpha*Mb)**2+xmAalpha*Ma*fA_aspol_alpha
Free_Energy=Free_Energy+xmAbeta+xmBbeta+xmSolventbeta+xmNabeta+xmClbeta
Free_Energy=Free_Energy+packconst-0.5*chi*(Ma*xmAbeta+MB*xmBbeta)**2-xmAbeta*Ma*fA_aspol_beta

elib=Free_Energy

return 
end subroutine

