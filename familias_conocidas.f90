subroutine familias_conocidas(sol)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Comprueba si cumple las familias !
! FS conocidas con anterioridad    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use variable
implicit none
integer sol,i,j,k,l,n
dimension sol(nind,7)

penal(2)=0
do n=1,2
  if(sum(fijo_fam(:,n)).eq.0) cycle! no hay familias de parental n (padre = 1 / madre = 2)
  do k=1,fijo_fam(0,n)
    do i=1,nind
	  if(sol(i,7).le.nvirt) cycle! es individuo virtual
	  if(fijo_fam(sol(i,7)-nvirt,n).ne.k) cycle! no es de la familia k
	  if(sol(i,n+1).eq.0) then
	    penal(2)=penal(2)+1! penaliza porque no le asigna parental
		cycle
	  endif
	  do l=i+1,nind
	    if(sol(l,7).le.nvirt) cycle! es individuo virtual
        if(fijo_fam(sol(l,7)-nvirt,n).ne.k) cycle! no es de la familia k
		if(sol(l,n+1).ne.sol(i,n+1)) penal(2)=penal(2)+1! penaliza porque no coinciden
	  enddo
	enddo
  enddo
enddo

end subroutine