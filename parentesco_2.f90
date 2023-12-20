subroutine parentesco_2(sol)
! Calcula el parentesco genealogico entre los genotipados
use variable
implicit none
integer::i,n,j,sol,cont,aux1
real ft
dimension ft(nind,nind),aux1(n_genot,2),sol(nind,7)


! Copia posicion y codigo de los genotipados en vector auxiliar
cont=0
do i=1,nind
!print('(7i4)'),sol(i,:)
  if(sol(i,6).eq.0) cycle
  cont=cont+1
  aux1(cont,1)=sol(i,1)
  aux1(cont,2)=sol(i,7)
enddo


! Ordena por codigo
call ordena(aux1,n_genot,2,2)


! Calcula la matriz de parentesco de todos (genotipados y virtuales)
ft=0
do i=1,nind
  if(product(sol(i,2:3)).eq.0) then
    ft(i,i)=.5
  else
    ft(i,i)=.5*(1+ft(sol(i,2),sol(i,3)))
  endif
  do j=i+1,nind
    do n=1,2
	  if(sol(j,n+1).eq.0) cycle
	  ft(i,j)=ft(i,j)+ft(i,sol(j,n+1))
	enddo
	ft(i,j)=ft(i,j)/2
	ft(j,i)=ft(i,j)
  enddo
enddo

! Copia los parentescos a la matriz genealogica
!cont=nind
do i=1,n_genot
  do j=i,n_genot
	ftf(i,j)=ft(aux1(i,1),aux1(j,1)) ! Parentescos cruzados
	ftf(j,i)=ftf(i,j)
  enddo
enddo

end subroutine

subroutine ordena(ix,nf,nc,nx)
! Ordena una matriz por los valores de una de sus columnas !
implicit none
integer i,j,nf,nc,nx,ix,ic
dimension ix(nf,nc),ic(nc)
do i=1,nf-1
  do j=i+1,nf
    if(ix(j,nx).lt.ix(i,nx)) then
      ic=ix(i,:)
      ix(i,:)=ix(j,:)
      ix(j,:)=ic
     endif
  enddo
enddo
end subroutine
