module variable
implicit none
integer marcadores,solopt,nact,nind,nvirt,nm,nh,nmr,clases_virt,id_ban,&
        nhr,ngen,nloci,nalel,par,fallo,falla,tolerancia,peso,aux,aux2,n_banned,&
	datos,fijo_p,fijo_fam,penal,n_genot,cont1,cont2,banned,maxproh,no_repr,salto
real fest,ftf,eneopt,puber,senes,pondera,fact
character*40 input
dimension fallo(3),tolerancia(3),peso(3),datos(5),penal(4),pondera(4),n_banned(2)
allocatable marcadores(:,:,:),fest(:,:),solopt(:,:),ftf(:,:),nalel(:),id_ban(:),&
            par(:,:),aux(:,:),aux2(:,:),fijo_p(:,:),fijo_fam(:,:),banned(:,:,:),no_repr(:),fact(:,:)
end module

program multi_gen_simul
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Busca la genealogía que da la maxima correlacion !
! entre la matriz de parentesco genealogico y la   !
! estimada a partir del parentesco molecular.      !
! Una o varias generaciones comprobando la         !
! compatibilidad de familias de hermanos, la       !
! paternidad y la edad de todos los ancestros.     !
! Permite relajar compatibilidad para tener en     !
! cuenta errores de genotipado o de datado         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

use variable
implicit none
integer nmin,l,j,k1,i,g,ant,n,k,k2,parmol,comp,iseed,m
real time1,time2,u
character*40 coanc
character*1 biallelic
allocatable comp(:,:)

print*,'Data files?'
read(*,*) input

input=adjustl(input)

open(UNIT=55,FILE=trim(input)//'_param')

! Lee semilla para números aleatorios *
read(55,*) iseed
!print*,'semilla',iseed
! Lee número de individuos reales (genotipados o no) y virtuales *
read(55,*) nact,nvirt
nind=nact+nvirt
!print*,'Individuos',nact,nvirt
! N. loci y de alelos por locus *
read(55,*) nloci,salto,biallelic
salto=nloci/salto
allocate(nalel(nloci))
select case (biallelic)
 case ('S')
   nalel=2
 case ('N')
   read(55,*) (nalel(n),n=1,nloci)
 case default
   print*,'******************************'
   print*,'Wrong code for the type of markers',biallelic
   print*,'******************************'
   pause
   stop
endselect

nmin=int(maxval(nalel)/2)+1
!print*,'loci',nloci,salto,'alelos',nalel(1:3)

! Minima y máxima edad reproductiva *
read(55,*) puber,senes
!print*,'edad',puber,senes
! N. de fallos tolerados para paternidad, familias de FS y edades *
read(55,*) tolerancia
!print*,'toleranca',tolerancia
! Importancia relativa de cada tipo de fallos *
read(55,*) peso
!print*,'peso',peso
! Importancia relativa de cada tipo de fallos datos conocidos a priori*
read(55,*) pondera
!print*,'pondera',pondera
if(nvirt.ne.0) then
  ! Lee numero de clases de edad para virtuales *
  read(55,*) clases_virt
  allocate(par(0:clases_virt,3))

  ! Lee numero de machos y hembras en virtuales *
  do i=0,clases_virt
    read(55,*) par(i,:)
  enddo

  if(sum(par(:,1:2)).ne.nvirt) then
    print*,'******************************'
    print*,'Check number of virtual males and females'
    do i=0,clases_virt
      print*,par(i,:)
      print*
      print*,nvirt
    enddo
    print*,'******************************'
    pause
    stop
  endif
endif

! Lee número de generaciones posibles *
read(55,*) ngen
!print*,'generaciones',ngen
! Calcula parentesco molecular o lo importa *
read(55,*) datos(1)

! Se conocen paternidades o no *
read(55,*) datos(2)

! Se conocen familias FS o no *
read(55,*) datos(3)

! Hay paternidades prohibidas o no *
read(55,*) datos(4)

! Hay individuos que no pueden ser parentales *
read(55,*) datos(5)
!print*,'datos',datos

allocate (solopt(nind,7))

! Rellena los datos de los individuos *
solopt=0
if(nvirt.ne.0) then
  ! Primero los de los virtuales *
  l=0
  do i=0,clases_virt
    do j=1,2
      do k=1,par(i,j)
        l=l+1
	    solopt(l,1)=l
	    solopt(l,4)=j
	    solopt(l,5)=par(i,3)
	    solopt(l,6)=0
	    solopt(l,7)=l
	  enddo
    enddo
  enddo
endif
! Despues lee los de los reales *
do i=nvirt+1,nind
  solopt(i,1)=i
  read(55,*) solopt(i,4:6)
!print*,'datos',i,solopt(i,4:6)
  solopt(i,7)=i
enddo
close(55)

open(UNIT=55,FILE=trim(input)//'_genot')

n_genot=sum(solopt(:,6))

! Datos de los marcadores *
allocate(marcadores(2,nloci,nact))
do i=1,nact
  read(55,*) ((marcadores(k,l,i),k=1,2),l=1,nloci)
  select case(solopt(i+nvirt,6))
    case (0)
      if(sum(marcadores(:,:,i)).ne.0) then
	    print*,'******************************'
        print*,'Code 0 for a genotyped individual',i
	    print*,'******************************'
	    pause
	    stop
	  endif
	case(1)
	case default
	  print*,'******************************'
      print*,'Wrong code for the genotype status of individual',i,solopt(i+nvirt,6)
	  print*,'******************************'
	  pause
	  stop
  endselect
!print*,'marcadores',i,marcadores(:,1,i)
enddo
!print*,'marcadores',sum(marcadores)
!pause
close(55)


allocate (fest(n_genot,n_genot))

select case (datos(1))
case(0)! Calcula el parentesco molecular *
print*,'calculating molecular coancestry'
    allocate (comp(n_genot,n_genot))
  comp=nloci
  fest=0
  falla=0
  cont1=0
  do j=1,nact
    if(solopt(j+nvirt,6).eq.0) cycle
	cont1=cont1+1
	cont2=cont1-1
    do i=j,nact
      if(solopt(i+nvirt,6).eq.0) cycle
	  cont2=cont2+1
      loci:do l=1,nloci
        parmol=0
        do k1=1,2
          if(marcadores(k1,l,j).eq.0) then
            comp(cont2,cont1)=comp(cont2,cont1)-1
            comp(cont1,cont2)=comp(cont1,cont2)-1
            parmol=0
	        cycle loci
	      endif
          do k2=1,2
            if(marcadores(k2,l,i).eq.0) then
            comp(cont2,cont1)=comp(cont2,cont1)-1
            comp(cont1,cont2)=comp(cont1,cont2)-1
	          parmol=0
	          cycle loci
            endif
            if(marcadores(k1,l,j).eq.marcadores(k2,l,i)) parmol=parmol+1
          enddo
        enddo
        fest(cont2,cont1)=fest(cont2,cont1)+parmol
      enddo loci
    enddo
  enddo
  do j=1,n_genot
    do i=j,n_genot
	  if(comp(i,j).eq.0) then
	    print*,i,j,' no common loci'
        falla=falla+1
	    cycle
      endif
      fest(i,j)=fest(i,j)/(comp(i,j)*4)
      fest(j,i)=fest(i,j)
	enddo
  enddo
  deallocate (comp)


! Escribe matriz de parentesco molecular *
  open(UNIT=55,FILE=trim(input)//'_IBS_coanc')
  do i=1,n_genot
    write(55,'(10000f8.5)') fest(i,:)
  enddo
  close(55)
print*,'molecular coancestry done'

  case(1)! Lee matriz de parentesco estimado de un archivo *

open(UNIT=55,FILE=trim(input)//'_coanc')
  do i=1,n_genot
    read(55,*) fest(i,:)
  enddo
  close(55)



  case default
  print*,'*******************************************************'
  print*,'Error in the data source (coancestry matrix file)',datos(1)
  print*,'*******************************************************'
  pause
  stop
end select

allocate (ftf(n_genot,n_genot))

penal=0

select case (datos(2))
  case(0)

  case(1)! Existe un archivo con paternidades conocidas
  allocate(fijo_p(nact,2))
  open(UNIT=55,FILE=trim(input)//'_known_parentage')
  do i=1,nact
    read(55,*) n,fijo_p(i,:)
    do j=1,2
      if((fijo_p(i,j).ne.0).and.(fijo_p(i,j).ne.-9)) fijo_p(i,j)=fijo_p(i,j)+nvirt
    enddo
  enddo
  close(55)

  case default
  print*,'**************************************************'
  print*,'Error in the data source (known paternities)',datos(2)
  print*,'**************************************************'
  pause
  stop
end select

select case (datos(3))
  case(0)

  case(1)! Existe un archivo con familias conocidas
  allocate(fijo_fam(0:nact,2))
  fijo_fam=0
  open(UNIT=55,FILE=trim(input)//'_known_families')
  do n=1,2
    read(55,*) k
    fijo_fam(0,n)=k
    do i=1,k
      read(55,*) l
  	  do j=1,l
	    read(55,*) m
	    fijo_fam(m,n)=i
	  enddo
	  read(55,*)
	enddo
  enddo
  close(55)

!do i=0,nact
!  print*,fijo_fam(i,:)
!enddo
!pause

  case default
  print*,'**************************************************'
  print*,'Error in the data source (known families)',datos(3)
  print*,'**************************************************'
  pause
  stop
end select

select case (datos(4))
  case(0)

  case(1)! Existe un archivo con incompatibilidades
  open(UNIT=55,FILE=trim(input)//'_banned_parentage')
  read(55,*) n_banned(2)
  read(55,*) maxproh
  allocate(banned(2,0:maxproh,n_banned(1)),id_ban(n_banned(1)))
  do i=1,n_banned(1)
    read(55,*) id_ban(i),banned(1,0,i),banned(1,1:banned(1,0,i),i),banned(2,0,i),banned(2,1:banned(2,0,i),i)
  enddo
  close(55)


  case default
  print*,'**************************************************'
  print*,'Error in the data source (banned paternities)',datos(4)
  print*,'**************************************************'
  pause
  stop
end select

select case (datos(5))
  case(0)

  case(1)! Existe un archivo con individuos que no se reprodujeron
  open(UNIT=55,FILE=trim(input)//'_no_reproduction')
  read(55,*) n_banned(2)
  allocate(no_repr(n_banned(2)))
  do i=1,n_banned(2)
    read(55,*) no_repr(i)
  enddo
  close(55)


  case default
  print*,'**************************************************'
  print*,'Error in the data source (no reproduced individuals)',datos(5)
  print*,'**************************************************'
  pause
  stop
end select


call cpu_time(time1)



! Inicializa generador de numeros aleatorios *
call setup(iseed)

! Busca la genealogia por annealing *
call annealing()

call cpu_time(time2)

open(UNIT=55,FILE=trim(input)//'_est_geneal.txt')
write(55,*) 'No. real individuals (genotyped or not)',nact
write(55,*)
write(55,*) 'No. genotyped individuals',n_genot
write(55,*)
write(55,*) 'No. virtual sires and dams'
do i=0,clases_virt
  write(55,*) par(i,:)
enddo
write(55,*)
write(55,*) 'No. initial generations',ngen
write(55,*)
write(55,*) 'No. loci',nloci
write(55,*)
write(55,*) 'Max. No. alelles/loci',maxval(nalel)
write(55,*)
write(55,*) 'No. allowed mismatches (parentage, FS(HS) families, age)'
write(55,*)  tolerancia
write(55,*) 'Calculations took ',time2-time1,' seconds'
write(55,*)

write(55,*) 'Correlation between input and solution coancestry matrix'
write(55,*) -eneopt
!do i=1,nind
!  write(55,*) solopt(i,1:3)
!enddo
write(55,*)
aux=0
do i=1,nind
  do j=1,2
    if(solopt(i,j+1).eq.0)cycle
    aux(i,j)=solopt(solopt(i,j+1),7)
  enddo
enddo

allocate(aux2(nind,3))
do i=1,nind
  aux2(i,1)=solopt(i,7)
  aux2(i,2:)=aux(i,:)
enddo

call ordena(aux2,nind,3,1)

do i=1,nind
  write(55,*) aux2(i,:)
enddo

deallocate(aux2)

close(55)

allocate(fact(nact,nact))
call parentesco(solopt)
open(UNIT=55,FILE=trim(input)//'_est_coanc.txt')
write(55,*) 'No. real individuals (genotyped or not)',nact
write(55,*)
write(55,*) 'No. genotyped individuals',n_genot
write(55,*)
write(55,*) 'No. virtual sires and dams'
do i=0,clases_virt
  write(55,*) par(i,:)
enddo
write(55,*)
write(55,*) 'No. initial generations',ngen
write(55,*)
write(55,*) 'No. loci',nloci
write(55,*)
write(55,*) 'Max. No. alelles/loci',maxval(nalel)
write(55,*)
write(55,*) 'No. allowed mismatches (parentage, FS(HS) families, age)'
write(55,*)  tolerancia
write(55,*) 'Calculations took ',time2-time1,' seconds'
write(55,*)

write(55,*) 'Correlation between input and solution coancestry matrix'
write(55,*) -eneopt

do i=1,nact
  write(55,'(10000f8.5)') fact(i,:)
enddo
close(55)

if(n_genot.ne.nact) then
  call parentesco_2(solopt)
  open(UNIT=55,FILE=trim(input)//'_est_coanc_genot.txt')
  do i=1,n_genot
    write(55,'(1000f8.5)') ftf(i,:)
  enddo
  close(55)
endif

print*,'****************************'
print*,'Normal ending of the program'
print*,'****************************'
!pause
!stop


end program

!*************
! SUBRUTINAS *
!*************

! Minimiza la funcion usando simulated annealing *
subroutine annealing()
use variable
implicit none
real t,kt,k,eneac,eneal,u,delta,omega
integer sol,sola,ch,s,rang,l,nmax,g,n,m,i,j,niv,rep,ivig,tot,prim,antes,ind,cambios,&
        pad,mad,ant,copia,acum,camb,loc,alelos,genotipos,hom,sumalelos,irep,iniv,tam
character*2 inicial
allocatable alelos(:),genotipos(:,:)
dimension sol(nind,7),sola(nind,7)
allocate (alelos(maxval(nalel)),genotipos(maxval(nalel),maxval(nalel)))


! Lee parametros del annealing *
open(UNIT=55,FILE='anneal_param.txt')
read(55,*) iniv
read(55,*) irep
read(55,*) cambios
read(55,*) t
read(55,*) k
read(55,*) inicial
close(55)

sol=solopt
print*,'starts simulated annealing'

select case(inicial)
 case('N') ! Solucion inicial al azar
   tam=int(nind/ngen)
   sol(1:tam+mod(nind,ngen),2:3)=0 ! Los 'tam' primeros (más los redondeos) son fundadores
   prim=tam+mod(nind,ngen)+1
   if(count(sol(1:tam,4)==1)*count(sol(1:tam,4)==2).ne.0) then
     do i=1,ngen-1
!print*,'generacion',i
       do j=prim,prim+tam-1
!print*,'individuo',j
         do s=1,2 ! Para padre y madre
           do
	     call mother(u)
	     n=min(int(u*tam)+prim-tam,prim-1)
             if(sol(n,4).eq.s) exit ! es del sexo requerido
	   enddo
	   sol(j,s+1)=n
	 enddo
       enddo
       prim=prim+tam
     enddo
   endif

 case('S') ! Lee la solucion inicial de un archivo
   open(UNIT=55,FILE=trim(input)//'_initial')
   do i=1,nind
     read(55,*) n,sol(i,2:3)
   enddo
   close(55)

 case default
   print*,'******************************'
   print*,'Wrong code for the type of markers',inicial
   print*,'******************************'
   pause
   stop
end select
print*,'initial solution'
!do i=1,nind
!  print('(9i5)'),sol(i,:)!,sol(sol(i,2),5)-sol(i,5),sol(sol(i,3),5)-sol(i,5)
!enddo
!print*

!print*,'marcadores',sum(marcadores)

! Calcula matriz parentesco genealogico solucion inicial *
!print*,'entra en parentesco geneal'
call parentesco_2(sol)

!print*,'parentesco inicial'
!do i=1,nind
!  print('(100f6.3)'),ftf(i,:)
!enddo
!print*

!print*,'sale de parentesco'
!print*,'marcadores',sum(marcadores)

! Comprueba cuantas incompatibilidades padres-hijo existen
!print*,'entra en paternidad'
      call paternidad(sol)

!print*,'sale de paternidad'
!print*,'marcadores',sum(marcadores)


! Comprueba cuantas incompatibilidades de familias_FS existen
!print*,'entra en familias_FS'
      call familias_FS(sol)

!print*,'sale de familias_FS'
!print*,'marcadores',sum(marcadores)


! Comprueba cuantas incompatibilidades de edades existen
!print*,'entra en edad'
      call edad(sol)

!print*,'sale de edad'

! Comprueba cuantas incompatibilidades de padres conocidos existen
	  if(datos(2).eq.1) call padres_conocidos(sol)

! Comprueba cuantas incompatibilidades de familias conocidas existen
	  if(datos(3).eq.1) call familias_conocidas(sol)

! Comprueba cuantas incompatibilidades de padres prohibidos existen
	  if(datos(4).eq.1) call prohibidos(sol)

! Comprueba cuantas incompatibilidades de individuos que no se reprodujeron existen
	  if(datos(5).eq.1) call prohibidos_2(sol)

! Calcula valor de estado inicial *
      call energia(eneac)
print*,fallo,'| ',penal
print*,-eneac
!pause

!print*,'sale de energia',eneac,fallo

      solopt=sol
      eneopt=eneac

allocate(aux(nind,2))

      kt=t/k
	ivig=irep-1

print*,'Optimisation starts'
print*,'        Level   Accep. sol.  Changes  Current sol.   Best solution'
print*,'        =====   ===========  =======  ============   ============='

! Bucle de niveles *
      do niv=1,iniv
open(UNIT=55,FILE=trim(input)//'_est_geneal_partial.txt')
aux=0
do i=1,nind
  do j=1,2
    if(solopt(i,j+1).eq.0)cycle
    aux(i,j)=solopt(solopt(i,j+1),7)
  enddo
enddo
do i=1,nind
  write(55,*) solopt(i,7),aux(i,:)
enddo

close(55)

	  camb=1+(cambios*ivig/irep)
        kt=kt*k
        ivig=0
! Bucle de repeticiones *
        do rep=1,irep
!print*,'replica',rep
            ! Genera solucion alternativa *
! Copia la solucion actual en la alternativa *
          sola=sol
	do l=1,camb
!print*,'entra en cambio'
!write(50,*) rep,l
	  call cambio_2(sola)
!print*,rep,l
	enddo
!print*,'nueva solucion'
!do i=1,nind
!  write(55,'(7i4)') sola(i,:)
!  print('(7i4)'),sola(i,:)
!enddo
!close(55)
!print*

! Reordena para que padres precedan a hijos *
      call reorganiza(sola)

!pause
!print*,'nueva solucion reorganiza'
!do i=1,nind
!  print('(7i4)'),sola(i,:)
!enddo
!print*

! Recodifica desde 1 a nind *
      call recodifica(sola)

!print*,'nueva solucion recodifica'
!do i=1,nind
!  print('(7i4)'),sola(i,:)
!enddo
!pause

close(50)

! Calcula matriz parentesco genealogico solucion alternativa *
      call parentesco_2(sola)

!print*,'sale de parentesco'

! Comprueba cuantas incompatibilidades padres-hijo existen
      if(peso(1).ne.0) call paternidad(sola)

!print*,'sale de paternidad'

! Comprueba cuantas incompatibilidades de familias_FS existen
      if(peso(2).ne.0) call familias_FS(sola)

!print*,'sale de familias_FS'

! Comprueba cuantas incompatibilidades de edades existen
	  if(peso(3).ne.0) call edad(sola)

!print*,'sale de edad'

! Comprueba cuantas incompatibilidades de padres conocidos existen
	  if(datos(2).eq.1) call padres_conocidos(sola)

!print*,'sale de padres conocidos'

! Comprueba cuantas incompatibilidades de familias conocidas existen
	  if(datos(3).eq.1) call familias_conocidas(sola)

!print*,'sale de familias conocidas'

! Comprueba cuantas incompatibilidades de padres prohibidos existen
	  if(datos(4).eq.1) call prohibidos(sol)

! Comprueba cuantas incompatibilidades de individuos que no se reprodujeron existen
	  if(datos(5).eq.1) call prohibidos_2(sol)

! Calcula valor de estado alternativo *
      call energia(eneal)

!print*,'sale de energia',eneal,fallo

! Anneling *
          ch=0
          delta=amax1(eneal-eneac,0.)
          omega=exp(-delta/kt)
          if(omega.ge.1) then
            ch=1
            if(eneal.lt.eneopt) then
              eneopt=eneal
              solopt=sola
!print*,fallo,'|',penal
            endif
          else
            call mother(u)
            if(u.lt.omega) then
              ch=1
            endif
          endif
          if(ch.eq.1) then
!print*,rep,eneal,eneac,omega,u,ch
!pause
            sol=sola
	      eneac=eneal
            ivig=ivig+1
          endif
        enddo
      print*,niv,ivig,1+(cambios*ivig/irep),-eneac,-eneopt
!print*,fallo,'|',penal
!pause
        if(ivig.eq.0) then
          return
        endif
      enddo
      print*,'Limit of steps/temperatures reached'
      end

      subroutine energia(corre)
! Calcula valor de la funcion objetivo *
      use variable
	implicit none
      integer i,usa
      real corre,sum1,sum2,sum1cuad,sum2cuad,sumprod

	usa=-(falla/2)+(n_genot+(n_genot*n_genot))/2

    sum1=0
	sum2=0
	sum1cuad=0
	sum2cuad=0
	sumprod=0

	do i=1,nact
      sum1=sum1+sum(ftf(i,i:))
	  sum2=sum2+sum(fest(i,i:))
	  sum1cuad=sum1cuad+sum(ftf(i,i:)*ftf(i,i:))
	  sum2cuad=sum2cuad+sum(fest(i,i:)*fest(i,i:))
	  sumprod=sumprod+sum(ftf(i,i:)*fest(i,i:))
	enddo

if((usa*sum1cuad-sum1*sum1)*(usa*sum2cuad-sum2*sum2).eq.0) then
corre=0
return

!print*,'suma de cuadrados 1',sum1cuad
!print*,'suma al cuadrado 1',sum1*sum1
!print*,'suma de cuadrados 2',sum2cuad
!print*,'suma al cuadrado 2',sum2*sum2
!print*,'producto cruzado',sumprod
!print*,usa
!do i=1,nact
!  print'(100f7.4)',fest(i,:)
!enddo
!print*
!do i=1,nact
!  print'(100f7.4)',ftf(i,:)
!enddo
!pause
endif

	corre=-(usa*sumprod-sum1*sum2)/&
           sqrt((usa*sum1cuad-sum1*sum1)*(usa*sum2cuad-sum2*sum2))


!    do i=1,3
!	  fallo(i)=max(fallo(i)-tolerancia(i),0)
!	enddo

	corre=corre+dot_product(fallo,peso)+sum(pondera*penal)

!print*,fallo,'|',penal
!print*,corre

      end subroutine
! Rutina para numeros aleatorios uniformes (0 - 1) *
      subroutine setup(iseed)
      Common /seeds/xseed(8),yseed(8),i16,i32
      integer iseed,xseed,yseed,carry
      i16=65535
      i32=2147483647
      n=iseed
      do 10 i=1,8
         k=mod(n,1000)
         carry=mod(n/1000,1000)
         n=672*k+carry
         xseed(i)=mod(n,i16)
         n=672*carry+k
         yseed(i)=mod(n,i16)
  10  continue
      end

      subroutine mother(u)
      Common /seeds/x(8),y(8),i16,i32
      integer x,y,carx,cary
      real u,w,t
      k=mod(x(8),21475)*100000+y(8)+mod(x(8),18112)
       if(k.lt.0)k=-k
      k=mod(k,i32)
      L=mod(y(8),21475)*100000+x(8)+mod(y(8),18112)
       if(L.lt.0)L=-L
      L=mod(L,i32)
      t=i32
      carx=mod(L,i16)
      cary=mod(k,i16)
       n=x(1)*12013 + x(2)*1066 + x(3)*1215 + x(4)*1492&
        + x(5)*1776 + x(6)*1812 + x(7)*1860 + x(8)*1941 + carx
       m=y(1)*9272 + y(2)*7777 + y(3)*6666 + y(4)*5555&
        + y(5)*4444 + y(6)*3333 + y(7)*2222 + y(8)*1111 + cary
      do 10 i=1,7
       x(i)=x(i+1)
  10   y(i)=y(i+1)

      x(8)=mod(n,i16)
      y(8)=mod(m,i16)
      k=mod(n,21475)*100000+mod(m,83648)
       if(k.lt.0)k=-k
      k=mod(k,i32)
      w=k
       u=w/t
      end

subroutine prohibidos(sol)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Comprueba si aparecen paternidades prohibidas !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use variable
implicit none
integer sol,i,j,k,l
dimension sol(nind,7)

!print*,'entra en prohibidos'

penal(3)=0
do j=1,n_banned(1)
  do i=1,nind
    if(id_ban(j).eq.sol(i,7)-nvirt) exit ! es individuo que tiene incompatibilidades
!print*,i,sol(i,7),sol(i,7)-nvirt,banned(:,0,sol(i,7)-nvirt)
    do k=1,2
      if(banned(k,0,j).eq.0) cycle ! no hay prohibidos para este parental
	  if(sol(i,k+1).eq.0) cycle ! aparece como fundador
      do l=1,banned(k,0,j)
!print*,sol(i,7),sol(i,j+1),sol(sol(i,j+1),7),sol(i,7)-nvirt,sol(sol(i,j+1),7)-nvirt,banned(j,k,sol(i,7)-nvirt)
	    if(sol(sol(i,k+1),7)-nvirt.eq.banned(k,l,j)) penal(3)=penal(3)+1 ! penaliza si coincide
!print*,penal(3)
      enddo
	enddo
  enddo
!pause
enddo

!print*,penal(3)
!print*,'sale de prohibidos'

end subroutine

subroutine prohibidos_2(sol)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Comprueba si aparecen como padres !
! individuos prohibidos             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use variable
implicit none
integer sol,i,j
dimension sol(nind,7)

!print*,'entra en prohibidos'

penal(4)=0
do i=1,n_banned(2)
  do j=1,nind
    if(no_repr(i).eq.sol(j,7)-nvirt) exit ! este inidividuo no puede ser parental
  enddo
  penal(4)=penal(4)+count(sol(:,2:3)==sol(j,1)) ! penaliza si aparece como parental

!pause
enddo

!print*,penal(4)
!print*,'sale de prohibidos'

end subroutine

subroutine parentesco(sol)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calcula el parentesco genealogico entre los reales !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use variable
implicit none
integer::i,n,j,sol,cont,aux1
real ft
dimension ft(nind,nind),aux1(nact,2),sol(nind,7)


! Copia posicion y codigo de los reales en vector auxiliar
cont=0
do i=1,nind
!print('(7i4)'),sol(i,:)
  if(sol(i,7).le.nvirt) cycle ! Es un individuo virtual
  cont=cont+1
  aux1(cont,1)=sol(i,1)
  aux1(cont,2)=sol(i,7)
enddo


! Ordena por codigo
call ordena(aux1,nact,2,2)


! Calcula la matriz de parentesco de todos (reales y virtuales)
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
do i=1,nact
  do j=i,nact
	fact(i,j)=ft(aux1(i,1),aux1(j,1)) ! Parentescos cruzados
	fact(j,i)=fact(i,j)
  enddo
enddo

end subroutine
