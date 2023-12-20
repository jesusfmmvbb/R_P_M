!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Renumera los individuos para que     !
! vayan de 1 al total correlativamente !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine recodifica(sol)
use variable
implicit none
integer a1,sol,i,j,k
dimension sol(nind,7),a1(nind,3)

a1=sol(:,1:3)

do i=1,nind
  parental:do j=2,3
!print*,i,j
  if(a1(i,j).eq.0) cycle parental ! No tiene padres conocidos que recodificar
    do k=1,nind
!print*,i,j,k
!print*,a1(i,j),a1(k,1)
	  if(a1(i,j).eq.a1(k,1)) then
	    sol(i,j)=k
		cycle parental
	  endif
	enddo
  enddo parental
enddo

do i=1,nind
  sol(i,1)=i
enddo

end subroutine