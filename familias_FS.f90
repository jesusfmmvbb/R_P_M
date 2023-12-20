subroutine familias_FS(sol)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Comprueba que las familias de FS sean compatibles !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use variable
implicit none
integer i,j,m,n,hom,alelos,genotipos,acum,sumalelos,sol,loc,exceso
allocatable alelos(:),genotipos(:,:)
dimension sol(nind,7)
allocate (alelos(maxval(nalel)),genotipos(maxval(nalel),maxval(nalel)))

fallo(2)=0
do n=1,nind
  exceso=0
  if(sol(n,6).eq.0) cycle ! Individuo sin genotipo
  if(product(sol(n,2:3)).eq.0) cycle ! Al menos un padre desconocido
!print*,'Primero',n
  loci:do loc=1,nloci,salto
    if(nalel(loc).eq.2) cycle loci ! No lo hace para bialelicos
    alelos=0
    genotipos=0
	if(product(marcadores(:,loc,sol(n,7)-nvirt)).eq.0) cycle ! Fallo de genotipado para el locus
    do m=n,nind
      if(product(sol(m,2:3)).eq.0) cycle ! Al menos un padre desconocido
      if(sol(m,6).eq.0) cycle ! Individuo sin genotipo
	  if(product(marcadores(:,loc,sol(m,7)-nvirt)).eq.0) cycle ! Fallo de genotipado para el locus
      if((sol(n,2).eq.sol(m,2)).and.(sol(n,3).eq.sol(m,3))) then ! Si tienen los mismos padres (son FS)
!print*,'Otros',m
        genotipos(marcadores(1,loc,sol(m,7)-nvirt),marcadores(2,loc,sol(m,7)-nvirt))=1
        genotipos(marcadores(2,loc,sol(m,7)-nvirt),marcadores(1,loc,sol(m,7)-nvirt))=1
        do hom=1,2
          alelos(marcadores(hom,loc,sol(m,7)-nvirt))=1
        enddo
      endif
    enddo
! Comprueba si hay mas de cuatro genotipos *
	acum=0
	do i=1,nalel(loc)
	  acum=acum+sum(genotipos(i,i:nalel(loc)))
	enddo
	if(acum.gt.4) then
	  exceso=exceso+1
!write(45,*) n,loc,'>cuatro genotipos'
	  cycle loci
	endif
! Comprueba si hay un alelo en heterocigosis con mas de dos alelos *
    do i=1,nalel(loc)
      if((sum(genotipos(i,1:i-1))+sum(genotipos(i,i+1:))).gt.2) then
        exceso=exceso+1
!write(45,*) n,loc,'tres heterocigotos'
	    cycle loci
	  endif
    enddo
    sumalelos=sum(alelos)
! Clasifica segun número de alelos *
    select case(sumalelos)
      case(5:)
	    exceso=exceso+1
!write(45,*) n,loc,sumalelos
	    cycle loci
      case(4)
        acum=0
        do i=1,nalel(loc)
          acum=acum+genotipos(i,i)
        enddo
        if(acum.ne.0) then        
	      exceso=exceso+1
!write(45,*) n,loc,sumalelos
 	      cycle loci
		endif
      case(3)
        acum=0
        do i=1,nalel(loc)
          acum=acum+genotipos(i,i)
        enddo
        if(acum.gt.1) then
          exceso=exceso+1
!write(45,*) n,loc,sumalelos
	      cycle loci
	    endif
    end select
  enddo loci
!  fallo(2)=fallo(2)+max(exceso-tolerancia(2),0)
  if(exceso.gt.tolerancia(2)) fallo(2)=fallo(2)+1
enddo
end subroutine