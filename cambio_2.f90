subroutine cambio_2(sol)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Genera un cambio en la        !
! solucion actual del annealing !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use variable
implicit none
integer sol,pos,tipo,cont,n,cambia,s,buc
real u
dimension sol(nind,7),pos(2),cont(2),cambia(4)

! Decide el tipo de cambio

call mother(u)
n=min(int(u*2)+1,2)


!print*,'cambia=',n

select case (n)
  case (1) ! Cambia el parental de un individuo (sin moverlo)
    do
      call mother(u)
      pos(1)=min(int(u*nind)+1,nind) ! Elije el individuo
      call mother(u)
      s=min(int(u*2)+1,2) ! Elije el sexo del parental
!	  if(count(sol(1:pos(1)-1,4)==s).gt.1) exit ! Acepta si hay algún individuo del sexo requerido por "encima" en la genealogia
exit
	enddo
!print*,'Primero',pos(1),s
!write(50,*) 'Primero',pos(1),s
	do
      call mother(u)
!      pos(2)=min(int(u*(pos(1)-1)),pos(1)-1) ! Elije el nuevo parental (anterior a el)
      pos(2)=min(int(u*(nind+1)),nind) ! Elije el nuevo parental entre todos los posibles
!print*,'Segundo',pos(2)
	  if(pos(2).eq.0) exit ! Parental desconocido
!write(50,*)'Segundo',pos(2),sol(pos(2),4),s
      if(pos(2).eq.pos(1)) cycle ! Rechaza si es el mismo individuo
	  if(sol(pos(2),4).ne.s) cycle ! Rechaza si no es del sexo adecuado
      if(sol(pos(2),1).eq.sol(pos(1),s+1)) cycle ! Rechaza si es el mismo parental


      buc=0
      call ancestro(sol(:,2:3),pos(1),pos(2),sol(pos(1),4),buc) ! Comprueba los bucles
      if(buc.eq.0) exit

	enddo
    if(pos(2).eq.0) then
      sol(pos(1),s+1)=0
!write(50,*)'Se queda con 0'
	else
	  sol(pos(1),s+1)=sol(pos(2),1)
!write(50,*)'Se queda con',pos(2),sol(pos(2),4),s
	endif


  case (2) ! Intercambia la posición de dos individuos en la genealogia
    call mother(u)
    pos(1)=min(int(u*nind)+1,nind) ! Elije el primer individuo
	s=sol(pos(1),4)
!print*,pos(1),s
	do
      call mother(u)
	  pos(2)=min(int(u*nind)+1,nind) ! Elije el segundo individuo
	  if(pos(2).eq.pos(1)) cycle ! Repite porque son el mismo
	  if(sol(pos(2),4).eq.s) exit ! Acepta si son del mismo sexo
	enddo
!print*,pos(1),pos(2),s

! Comprueba si hay padres reciprocos
    cont=0
	if(sol(pos(1),1+s).eq.sol(pos(2),1)) cont(1)=1
	if(sol(pos(2),1+s).eq.sol(pos(1),1)) cont(2)=1
	if(product(cont).ne.0) then
	  print*,'*************************'
	  print*,'Error. Reciprocal parents'
	  print*,'*************************'
	  pause
	  stop
	endif

! Cambia los individuos elegidos
    cambia=sol(pos(1),4:7)
	sol(pos(1),4:7)=sol(pos(2),4:7)
	sol(pos(2),4:7)=cambia
endselect

endsubroutine

recursive subroutine ancestro(ped,act,desc,s,buc)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Comprueba que no se producen  !
! bucles al cambiar el parental !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use variable
implicit none
integer ped,buc,act,desc,i,s,temp
dimension ped(nind,2)

!print*,act,desc,ped(desc,:),buc
!pause

if(buc.eq.1) return
if((ped(desc,1).eq.act).or.(ped(desc,2).eq.act)) then
  buc=1
  return
else
  do i=1,2
    if(ped(desc,i).eq.0) cycle
	temp=ped(desc,i)
	call ancestro(ped,act,temp,s,buc)
  enddo
endif

endsubroutine