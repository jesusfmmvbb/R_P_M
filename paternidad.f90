subroutine paternidad(sol)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Comprueba coherencia en paternidades !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use variable

implicit none
integer sol,i,j,k,k1,k2,s,desc,pad,hijo,conpad,cont,loc,o,p,h,exceso
dimension sol(nind,7),desc(2),pad(2,2),hijo(2),conpad(2,2)

fallo(1)=0

do i=1,nind
desc=0
exceso=0
  if(sol(i,6).eq.0) cycle ! Si no tiene genotipo salta al siguiente
  do s=2,3
    if(sol(i,s).eq.0) then ! Cuenta padres desconocidos
	  desc(s-1)=desc(s-1)+1 
	  cycle
	endif
	if(sol(sol(i,s),6).eq.0) desc(s-1)=desc(s-1)+1 ! Parental no genotipado
  enddo

!print*,desc
!print*,i,count(desc/=0)
  select case (count(desc/=0)) ! Cuenta padres no disponibles 
    case (2) ! No se conocen los padres (o no genotipados) y se pasa al siguiente
	  cycle
	case (1) ! Se conoce genotipo de un solo padre
	  do s=2,3
	    if(sol(i,s).eq.0) cycle ! Parental desconocido
		if(sol(sol(i,s),6).eq.0) cycle ! Parental no genotipado
		do loc=1,nloci,salto
		  cont=0
		  do k1=1,2
			if(marcadores(k1,loc,sol(i,7)-nvirt).eq.0) cont=cont+1
		    do k2=1,2
			  if(marcadores(k2,loc,sol(sol(i,s),7)-nvirt).eq.0) cont=cont+1
			  if(marcadores(k1,loc,sol(i,7)-nvirt).eq.marcadores(k2,loc,sol(sol(i,s),7)-nvirt)) cont=cont+1
			enddo
		  enddo
		  if(cont.eq.0) exceso=exceso+1 ! No coincide ninguno de los alelos de padre e hijo
!		  if(cont.eq.0) write(55,*) loc,sol(i,7),sol(sol(i,2),7),sol(sol(i,3),7) ! No coincide ninguno de los alelos de padre e hijo
!print*,i,sol(i,s),sol(sol(i,s),7)
!if(cont.eq.0) pause
		enddo
	  enddo
	case (0) ! Se conocen los genotipos de los dos padres
	  do loc=1,nloci,salto
        do k=1,2
          hijo(k)=marcadores(k,loc,sol(i,7)-nvirt)
	    enddo
	    if(product(hijo).eq.0) cycle
        do s=2,3
          do k=1,2
	        pad(s-1,k)=marcadores(k,loc,sol(sol(i,s),7)-nvirt)
		  enddo
		enddo
    	if(product(pad).eq.0) cycle
		conpad=0
		do o=1,2
		  do p=1,2
		    do h=1,2
			  if(hijo(o).eq.pad(p,h)) conpad(p,o)=conpad(p,o)+1
			enddo
		  enddo
		enddo
		if((sum(conpad(:,1))*sum(conpad(:,2)).eq.0).or.&
		   (sum(conpad(1,:))*sum(conpad(2,:)).eq.0)) exceso=exceso+1
!		if((sum(conpad(:,1))*sum(conpad(:,2)).eq.0).or.&
!		   (sum(conpad(1,:))*sum(conpad(2,:)).eq.0)) write(55,*) loc,sol(i,7),sol(sol(i,2),7),sol(sol(i,3),7)
	  enddo
  end select
!  fallo(1)=fallo(1)+max(exceso-tolerancia(1),0)
  if(exceso.gt.tolerancia(1)) fallo(1)=fallo(1)+1
enddo
end subroutine