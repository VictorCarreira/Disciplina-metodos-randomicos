!     PROGRAMA - modelagem magnetica
!     VARIAVEIS UTILIZADAS
!	xp,yp e zp = COORDENADAS DO PONTO DE OBSERVAÇÃO (metros)
!	xq,yq,zq = coordenadas do centro da esfera (metros)
!	mi = inclinação da magnetização (graus, positivo abaixo da horizontal) 
!	md = declinação da magnetização (graus, positiva para Leste, a partir do Norte verdadeiro)
!	mom = momento de dipolo magnetico (A.m²)
!     bx,by,bz,t = elementos do campo B gerados delo dipolo - dado de saida (nT)
!	eixo x, na direção Norte
!	eixu y, na direção Leste
!	eixo z, para baixo 
	

!	DEFINICAO DO TIPO DAS VARIAVEIS
	implicit real*8(a-h,o-z)

        real*8 perfil(10000),anom(10000),anomc(10000)
 	real*8 m,mi,md,mom,inp(3)
	character*30 alfa

      open(1,file='input.txt') 
      open(2,file='anomalia.txt') 
      open(3,file='fit.txt') 

	do i=1,3
	read(1,222) alfa,xnum
	inp(i)=xnum
	end do
	close(1)
222	format(A7,f10.4)

	ig=1
	do while (.true.)
	read(2,*,end=111) a1,a2   
	perfil(ig)=a1
	anom(ig)=a2
	ig=ig+1
	end do
111	continue
	close(2)

	np=ig-1

	zp=-1d0
	yq=0d0

!!  caracteristicas do campo geomagnético local
	mi=-34d0
	md=0d0
	F=23500d0

	xq=inp(1)
	zq=inp(2)
	mom=inp(3)

!    calculo das componentes do campo regional
	
	call dircos(mi,md,amx,amy,amz)

	fx=F*amx
	fy=F*amy
	fz=F*amz


!	 ----cálculo da amomalia de campo total através das subrotinas dipole e dircos

	do i=1,np
	xp=perfil(i)
!	yp=yp
	call dipole(xq,yq,zq,mi,md,mom,xp,yp,zp,bx,by,bz)
	bxx=bx+fx
	byy=by+fy
	bzz=bz+fz
	bt=sqrt(bxx**2+byy**2+bzz**2)

	anomc(i)=bt-F

	write(3,11) perfil(i),anomc(i),anom(i)

	end do

11	format(3(ES14.6E3,2x))


	write(6,*) "numero de pontos=",np

!  cálculo do Rrms

	soma=0d0
	do i=1,np
	soma=soma+(anomc(i)-anom(i))**2
	end do

	write(6,*) "============="

	Rrms=dsqrt(soma/dfloat(np))

	write(6,*) "Rrms=", Rrms , 'nT'
	
	stop
	end
	
!!!!!!! SUBROTINAS  !!!!!!!!!!

	subroutine dipole(xq,yq,zq,mi,md,moment,xp,yp,zp,bx,by,bz)

!     VARIAVEIS UTILIZADAS
!	xp,yp e zp = COORDENADAS DO PONTO DE OBSERVAÇÃO (metros)
!	xq,yq,zq = coordenadas do centro da esfera (metros)
!	mi = inclinação da magnetização (graus, positivo abaixo da horizontal) 
!	md = declinação da magnetização (graus, positiva para Leste, a partir do Norte verdadeiro)
!	moment = momento de dipolo magnetico (A.m²)
!	bx,by,bz,t = elementos do campo B gerados delo dipolo - dado de saida (nT)
!	eixo x, na direção Norte
!	eixu y, na direção Leste
!	eixo z, para baixo 

	implicit real*8 (a-h,o-z)
     	real*8 mi,md,m,mx,my,mz,moment
	data t2nt/1.d9/,cm/1.d-7/
	pi=2.d0*acos(0.d0)
	call dircos(mi,md,mx,my,mz)
	rx=xp-xq
	ry=yp-yq
	rz=zp-zq
	r2=rx**2+ry**2+rz**2
	r=sqrt(r2)
	if(r.eq.0.)pause 'DIPOLE: bad argument detected.'
	r5=r**5
	dot=rx*mx+ry*my+rz*mz
	bx=cm*moment*(3.d0*dot*rx-r2*mx)/r5
	by=cm*moment*(3.d0*dot*ry-r2*my)/r5
	bz=cm*moment*(3.d0*dot*rz-r2*mz)/r5
	bx=bx*t2nt
	by=by*t2nt
	bz=bz*t2nt
	return
	end	
!!!!!!!!!!!!!!!!!!

    	subroutine dircos(incl,decl,a,b,c)
	implicit real*8 (a-h,o-z) 
      real*8 incl
	pi=2.d0*acos(0.d0)
	d2rad=pi/180.d0
	xincl=incl*d2rad
	xdecl=decl*d2rad
	a=cos(xincl)*cos(xdecl)
	b=cos(xincl)*sin(xdecl)
	c=sin(xincl)
	return
	end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 

    
