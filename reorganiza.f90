!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reorganiza genealogia para que los !
! padres esten antes que los hijos   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine reorganiza(sol)
use variable
implicit none
integer pos,cont,sol,i,j,cambia,reales,rehacer
dimension sol(nind,7),cambia(7)

orden:do
  rehacer=0
  do i=1,nind ! Desde el principio de la genealogia
    pos=i
    reales=count(sol(pos,2:3)/=0) ! Cuenta los padres que se le concocen
    if(reales.eq.0) cycle ! Si no tiene padres pasa al siguiente
    do
      cont=0
! Mira si sus padres están antes que el !
	  do j=2,3
	    cont=cont+min(1,count(sol(1:pos-1,1)==sol(pos,j)))
	  enddo
	  select case (reales-cont) ! Los padres faltantes
	    case (0) ! Como estan pasa al siguiente
!print*,pos,reales-cont,sol(pos,1:3),sol(pos,7)
!write(50,*) pos,reales-cont,sol(pos,1:3),sol(pos,7)
	      exit 
	    case (1,2) !Como no estan "baja" el individuo hasta que esten
!print*,pos,reales-cont,sol(pos,1:3),sol(pos,7)
!write(50,*) pos,reales-cont,sol(pos,1:3),sol(pos,7)
	      cambia=sol(pos,:)
		  sol(pos,:)=sol(pos+1,:)
		  sol(pos+1,:)=cambia
		  pos=pos+1
		  rehacer=1
	    case default ! Por si acaso
	      print*, '****************************************'
	      print*, 'Error in parents known but not detected',reales,cont
	      print*, '****************************************'
		  pause
		  stop
	  endselect
    enddo
	if(rehacer.eq.1) cycle orden ! Si ha habido algún cambio vuelve a comprobar
  enddo
  exit
enddo orden
!close(50)
end subroutine