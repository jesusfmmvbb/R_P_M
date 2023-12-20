subroutine edad(sol)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Comprueba coherencia de edades !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use variable
implicit none
integer::i,n,j,sol,pos,gen,buc
dimension sol(nind,7)

fallo(3)=0
gen=0

!do i=1,nind
!  print('(7i3)'),sol(i,:)
!enddo
!pause

do i=1,nind
!print*,'individuo',i,gen
  if(sol(i,5).lt.0) cycle ! No se sabe edad del individuo
  buc=0
  call error_edad(sol(:,2:5),i,i,buc,gen)
  fallo(3)=fallo(3)+buc
!print*,i,fallo(3)
!print*,'fallos',fallo(3)
enddo
!pause

end subroutine


recursive subroutine error_edad(ped,act,anc,buc,gen)
use variable
implicit none
integer i,j,ped,buc,act,anc,temp,gen,tempgen
dimension ped(nind,4)

!print*,'act',act,'anc',anc,'buc',buc,'gen',gen

tempgen=gen+1
do i=1,2
if(buc.eq.1) return
  temp=ped(anc,i)
!print*,i,temp
  if(temp.eq.0) cycle !No se sabe el padre
!print*,temp,ped(temp,4),ped(act,4),puber,tempgen
  if(ped(temp,4).lt.0) then !El padre no tiene edad
    elseif((ped(temp,4)-ped(act,4)).lt.(puber*tempgen).or.(ped(temp,4)-ped(act,4)).gt.(senes*tempgen)) then
      buc=1
!print*,'buc'
	  return
  endif
  call error_edad(ped,act,temp,buc,tempgen)
enddo

end subroutine