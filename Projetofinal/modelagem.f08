PROGRAM modmag
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!Programa de modelagem direta magnética                  !
!Criador: Victor Carreira                                !
!Este programa é baseado em Blakely(1995) e Cosme(2018)  !
!atualizado para o FORTRAN 2008.                         !
!Para execução: gfortran modelagem.f08 -o make           !
!Para compilação: ./make                                 !
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!--------------------------------------------------------!
! Para usar compilação com flags utilize:                !
!-> gfortran -fbounds-check -fbacktrace -Wall -Wextra    !
!-pedantic modelagem.f08 -o nomedoexecutável             !
!--------------------------------------------------------!


   !+++++++++++++++++++++++++++++++++++++++++++++++++++!
   !-----------------Lista de Variáveis----------------!
   !+++++++++++++++++++++++++++++++++++++++++++++++++++!
   !xp,yp,zp: coordenadas do pto de obs                !
   !xq,yq,zq: coordenadas do centro da esfera          !
   !mi: inclinação da magnetização (graus,             !
   !positivo abaixo da linha horizontal)               !
   !md: declinação da magnetização (graus,             !
   !positiva para Leste, a partir do Norte verdadeiro) !
   !mom: momento de dipolo magnetico (A.m²)            !
   !bx,by,bz,t: elementos do campo B gerados pelo      ! 
   !dipolo. Saída em (nT).                             !
   !eixo x, na direção Norte                           !
   !eixo y, na direção Leste                           !
   !eixo z, para baixo                                 !
   !+++++++++++++++++++++++++++++++++++++++++++++++++++!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!DEFINIÇÃO DE VARIÁVEIS GLOBAIS!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IMPLICIT NONE
  INTEGER,PARAMETER::SP = SELECTED_INT_KIND(r=8)
  INTEGER,PARAMETER::DP = SELECTED_REAL_KIND(8,10)
  INTEGER(KIND=SP):: i,ig,np
  REAL(KIND=DP),ALLOCATABLE,DIMENSION(:):: perfil,anom,anomc,inp
  REAL(KIND=DP):: mi,md,mom,a1,a2,zp,yq,F,xq,zq,amx,amy,amz,soma,rms!,m
  REAL(KIND=DP)::ti,tf,tt,bt,bx,bxx,by,byy,bz,bzz,fx,fy,fz,xp,yp
  REAL(KIND=DP),PARAMETER::az=0.0!,ra=1.0

  CALL cpu_time(ti)

  ALLOCATE(perfil(1000),anom(1000),anomc(1000),inp(3))

  !Digite os dados de entrada no modelo no terminal:
  WRITE(*,'(a,$)')'xq:';READ*,inp(1)!Posição em x
  WRITE(*,'(a,$)')'zq:'; READ*,inp(2)!Profundidade em z
  WRITE(*,'(a,$)')'Momentum:'; READ*,inp(3)!Momento magnético

  PRINT*,"Vetor de entradas:",inp

  OPEN(2,FILE='anomalia.txt')!Dados reais (D.obs)
  OPEN(3,FILE='fit.txt')!Ajuste

   ig=1
    DO WHILE(.TRUE.)
     READ(2,FMT=*,END=111)a1,a2
      perfil(ig)=a1
      anom(ig)=a2
      ig=ig+1
    END DO
   111 CONTINUE
   CLOSE(2)

   np=ig-1

   zp=-1!Altura do magnetômero
   yq=0!Problema bi-dimensional não há dimensão y no ponto de medida.
   yp=0!Problema bi-dimensional não há dimensão y no corpo anômalo.

!Características do Campo Geomagnético local (Rio de Janeiro)
   mi=-34
   md=0!-22.54 !(IGRF-2014. obtido em http://extranet.on.br/jlkm/magdec/index.html)
   F=23500

!Caracteristicas da Fonte Magnética
   xq=inp(1)
   zq=inp(2)
   mom=inp(3)

!Cálculo das componentes do campo Regional (Planeta Terra)

  CALL dircos(mi,md,az,amx,amy,amz)!no caso 2D considera-se o azimute 0

   fx=F*amx
   fy=F*amy
   fz=F*amz

!Cálculo da anomalia do campo total
 
   DO i=1,np
     xp=perfil(i)
     CALL dipole(xp,yp,zp,mi,md,mom,xp,yp,zp,bx,by,bz)!O problema do NaN está aqui
       bxx=bx+fx
       byy=by+fy
       bzz=bz+fz
       bt=SQRT(bxx**2+byy**2+bzz**2)
        anomc(i)=bt-F
     WRITE(3,FMT=11) perfil(i),anomc(i),anom(i)
   END DO 
     WRITE(6,FMT=*)'======== AVALIAÇÃO ======='
     WRITE(6,FMT=*)'Número de pontos=',np

!Cálculo do Resíduo:

  soma=0
   DO i=1,np
      soma=soma+(anomc(i)-anom(i))**2
      !print*,soma
   END DO 

  rms=DSQRT(soma/DFLOAT(np))

    WRITE(6,FMT=*)'Resíduo (RMS)',rms,"nT" !Escreve na tela o RMS

CALL cpu_time(tf)
    
   tt=tf-ti

   WRITE(6,FMT=*)'Tempo de máquina=',tt,"segundos"


!!!!!!!!!!!!!!!!!!!!!FORMATOS UTILIZADOS!!!!!!!!!!!!!!!!!!!!!!!

11 FORMAT(3(ES14.6E3,2X))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!----------------------------FIM--------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
 CONTAINS

 SUBROUTINE dircos(incl,decl,azim,a,b,c)
!This subroutine computes direction cosines from inclination and
!declination.

!INPUT PARAMETERS:
!incl: inclination in degrees positive below horizontal.
!decl: declination in degrees positive east of true north.
!azim: azimuth of x axis in degrees positive east of north.

!OUTPUT PARAMETERS:
!a,b,c: the three direction cosines.

 IMPLICIT NONE
    REAL(KIND=DP), INTENT(IN)::incl,decl,azim
    REAL(KIND=DP), INTENT(OUT)::a,b,c
    REAL(KIND=DP):: xincl,xdecl,pi,d2rad
    REAL(KIND=DP), PARAMETER:: xazim=0.0 !No caso onde o corpo não está inclinado
    
    !d2rad=0.017453293

    pi=2.0*ACOS(0.0)
    d2rad=pi/180.0
   
     xincl=incl*d2rad
     xdecl=decl*d2rad
     !xazim=azim*d2rad
     a=COS(xincl)*COS(xdecl-xazim)
     b=COS(xincl)*SIN(xdecl-xazim)
     c=SIN(xincl)

END SUBROUTINE dircos

!-----------------------------------------------------------------

SUBROUTINE dipole(xq,yq,zq,mi,md,moment,xp,yp,zp,bx,by,bz)
!This subroutine computes the three components of magnetic induction
!caused a uniformly magnetized sphere. X axis is north, Z axis is
!donw. It is a adaptation of Blakely "dipole" subroutine to input 
!momentum. 

!INPUT PARAMETERS:
! Observation point located at (xp,yp,zp). Shpere centered at (xq,yq,zq).
!Dipole Momentum of sphere defined by "moment" (A.m²), inclination "mi",
!and declination "md". Units of distance irrelevant but must be consistent.
!All angles in degrees. Intensity of magnetization in A/m. Requires
!subroutine DIRCOS.

!OUTPUT PARAMETERS:
! The three components of magnetic induction (bx,by,bz) in units of nT.
IMPLICIT NONE
    REAL(KIND=DP), INTENT(IN)::xq,yq,zq,mi,md,xp,yp,zp,moment!,m,ra
    REAL(KIND=DP), INTENT(OUT)::bx,by,bz
    REAL(KIND=DP):: mx,my,mz,r,rx,ry,rz,r2,r5,dot
    REAL(KIND=DP), PARAMETER::t2nt=1.0E9,cm=1.0E-7,az=0.0!,pi=3.14159265
    

    print*

    CALL dircos(mi,md,az,mx,my,mz)
    !print*,rx,ry,rz
     rx=xp-xq
     ry=yp-yq
     rz=zp-zq
     r2=rx**2+ry**2+rz**2
     r=SQRT(r2)
     !IF(r.eq.0.0)PAUSE 'Dipole: Bad argument detected!'
     r5=r**5
     dot=rx*mx+ry*my+rz*mz
     !moment=4.0*pi*(ra**3)*m/3.0 utilizado quando a entrada é o magnetização(m) e o volume(ra)
     bx=cm*moment*(3.0*dot*rx-r2*mx)/r5
     by=cm*moment*(3.0*dot*ry-r2*my)/r5
     bz=cm*moment*(3.0*dot*rz-r2*mz)/r5
     bx=bx*t2nt
     by=by*t2nt
     bz=bz*t2nt

END SUBROUTINE dipole



END PROGRAM modmag
