subroutine padres_conocidos(sol)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Comprueba si cumple las paternidades !
! conocidas con anterioridad           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use variable
implicit none
integer sol,i,j,padre
dimension sol(nind,7)

penal(1)=0
do i=1,nind
  if(sol(i,7).le.nvirt) cycle ! es individuo virtual
  if(sum(fijo_p(sol(i,7)-nvirt,:)).eq.0) cycle ! no se saben padres
  do j=1,2
!print*,j,fijo_p(sol(i,7)-nvirt,j)
    if(fijo_p(sol(i,7)-nvirt,j).eq.0) cycle ! este parental no se sabe
    if(fijo_p(sol(i,7)-nvirt,j).eq.-9) then
	padre=sol(i,j+1)
	if(padre.eq.0) cycle ! puede ser fundador directamente
!print*,i,sol(i,7),fijo_p(sol(i,7)-nvirt,:)
!print*,i,sol(i,7),sol(i,7)-nvirt,j,padre
!print*,sol(padre,2:3),sol(sol(padre,2),7),sol(sol(padre,3),7)
      if(sum(sol(padre,2:3)).eq.0) cycle ! el padre no es nacido en cautividad
!print*,'no son padres fundadores'
	  if((sol(sol(padre,2),7).gt.nvirt).or.(sol(sol(padre,3),7).gt.nvirt)) then
!print*,'es hijo de no fundadores'
	    penal(1)=penal(1)+1! penaliza porque es hijo de uno con padres conocidos
!pause
		cycle
	  endif
	  cycle
	endif
	if(sol(i,j+1).eq.0) then
!print*,'no asigna',i,j,sol(i,7)-nvirt,fijo_p(sol(i,7)-nvirt,j),sol(sol(i,j+1),7)
!pause
	  penal(1)=penal(1)+1 ! penaliza porque no le asigna parental
	  cycle
	endif
	if(fijo_p(sol(i,7)-nvirt,j).ne.sol(sol(i,j+1),7)) penal(1)=penal(1)+1 ! penaliza si no coincide
!print*,'diferente',i,j,sol(i,7)-nvirt,fijo_p(sol(i,7)-nvirt,j),sol(sol(i,j+1),7)
!print*,penal(1)
!pause
  enddo
enddo
!print*,penal(1)
!pause
end subroutine