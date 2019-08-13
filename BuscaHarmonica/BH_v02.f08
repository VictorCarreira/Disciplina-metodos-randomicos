PROGRAM BuscaHarmonica
! Ajuste de uma reta via busca harmonica
!
!                 f(x)=ax+b            
!
!
!     f(x)|
!         |          /  * 
!         |        */
!         |        /
!         |      */
!         |   *  /
!         |     /
!         |    /*     
!         | * /    
!         |  /   
!         | /  
!       b_|/ * 
!         |____________________________________ 
!                                              x   
!
!*-dados (xd,yd)
!/-ptos da reta: dado calculado (xc,yc)
!a,b-parâmetros do modelo
!O melhor modelo é aquele que minimiza a fução data misfit (phi):
!                  Φ=Σ(yd-yc)²
!                  Φ=Σ(yd-axd-b)²
!Entra (xd[i],yd[i]) e sai (a,b).


    !-----------------Sumário das Variáveis-----------------!
    !HS: Harmony memory                                     !
    !HMS: Harmony memory size (number of measured points)   !
    !HMCR: Harmony consideration rate                       !
    !PAR: Pitch adjusting rate                              !
    !bw : termination criterion                             !
    !MaxItr: Maximal iteraction                             !
    !NVAR: Number of parameters                             !
    !-------------------------------------------------------!

!     O HMRC pode ser definido como a taxa de probabilidade de se selecionar um 
!componente da memória de harmonia atual, Harmony memobry - HM. Daí a diferença
!'1-HMCR' é a probabilidade de se gerar indivíduos aleatoriamente. Por exemplo, 
!HMCR=0,70 indica uma probabilidade de 70% de se selecionar um novo indivíduo
!a partir de HM. Ao passo que PAR, por sua vez, pode ser definido como a probabi
!lidade de um candidato presente em HM ser modificado, através de uma iteração. 
!Por exemplo, PAR=0,1 representa uma probabilidade de 10% de um indivíduo adjacen
!te a HM ser selecionado através do ajuste de passo. Consequentemente, '1-PAR' in
!dica a probabilidade de não selecionar um indivíduo adjacente, distanciado do passo.
!Os parâmetros HMCR e PAR do algoritmo controlam a geração de soluções candidatas
!, assim como a velocidade de convergência.



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!DEFINIÇÃO DE VARIÁVEIS GLOBAIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
 IMPLICIT NONE
   INTEGER,PARAMETER::SP = SELECTED_INT_KIND(r=4)
   INTEGER,PARAMETER::DP = SELECTED_REAL_KIND(8,8)
   INTEGER(KIND=DP)::NVAR,NG,NH,HMS
   INTEGER(KIND=DP)::currentIteration,MaxItr,i,j,index
   REAL(KIND=DP)::PARmin,PARmax,bwmin,bwmax,HMCR,PAR,newFitness
   REAL(KIND=DP)::coef,RES,RAN,pvbRan,ti,tf,tt
   INTEGER(KIND=DP)::pp
   REAL(KIND=DP), ALLOCATABLE, DIMENSION(:)::fit,NCHV,BestGen,BW,gx
   REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:)::HM,PVB,dado

 
OPEN(1,FILE='entrada.txt')
OPEN(2,FILE='saida.txt')



!Dimensiona os parâmetros da Busca:
 NVAR=2 !número de parâmetros da matriz HM
 NG=1
 NH=0
 MaxItr=5000
 HMS=5 ! número de linhas da matriz HM
 HMCR=0.6
 PARmin=0.4
 PARmax=0.85
 bwmin=0.01
 bwmax=4

 CALL cpu_time(ti)

 ALLOCATE(HM(1:HMS,1:NVAR),PVB(1:NVAR,1:2),dado(6,2))
 ALLOCATE(fit(1:HMS),NCHV(1:NVAR),BestGen(1:NVAR),BW(1:NVAR),gx(1:NG))

 !ZERANDO VARIÁVEIS:
 HM=0.0d0
 PVB=0.0d0
 dado=0.0d0
 fit=0.0d0
 NCHV=0.0d0
 BestGen=0.0d0
 BW=0.0d0
 gx=0.0d0


 DO i=1,6
  READ(1,FMT=*)dado(i,1),dado(i,2) 
 END DO  


 CALL INITIALIZE(HMS,NVAR,fit,PVB,HM)
 



print*,'Matriz HM(inicial)=' 
print*,HM(1,1),HM(1,2),phi(dado,HM(1,1),HM(1,2))
print*,HM(2,1),HM(2,2),phi(dado,HM(2,1),HM(2,2))
print*,HM(3,1),HM(3,2),phi(dado,HM(3,1),HM(3,2))
print*,HM(4,1),HM(4,2),phi(dado,HM(4,1),HM(4,2))
print*,HM(5,1),HM(5,2),phi(dado,HM(5,1),HM(5,2))
print*,"========================================"



 currentIteration  = 0

DO WHILE(currentIteration<MaxItr)
   PAR=(PARmax-PARmin)/(MaxItr)*currentIteration+PARmin
   coef=LOG(bwmin/bwmax)/MaxItr
  DO pp =1,NVAR
     BW(pp)=bwmax*EXP(coef*currentIteration)
  END DO
    DO i=1,NVAR
       CALL RANDOM_NUMBER(RAN)
       IF( ran < HMCR ) THEN         
         index = randint(1,HMS)
         NCHV(i) = HM(index,i)
         CALL RANDOM_NUMBER(pvbRan)
       IF( pvbRan < PAR) THEN
         CALL RANDOM_NUMBER(pvbRan)
         RES = NCHV(i)
       IF( pvbRan < 0.5) THEN
         CALL RANDOM_NUMBER(pvbRan)
         RES =RES+  pvbRan * BW(i)
       IF( RES < PVB(i,2)) THEN
         NCHV(i) = RES
       END IF
       ELSE
         CALL RANDOM_NUMBER(pvbRan)
         RES =RES- pvbRan * BW(i)
       IF( RES > PVB(i,1)) THEN
         NCHV(i) = RES
       END IF
       END IF
       END IF    
       ELSE
        NCHV(i) = randval( PVB(i,1), PVB(i,2) )
       END IF       
   newFitness = Phi(dado,HM(i,1),HM(i,2))
   CALL UpdateHM(currentIteration,HMS,BestGen,fit,NCHV, HM,  newFitness )


  


   currentIteration=currentIteration+1
  END DO

END DO




print*,'Matriz HM=' 
print*,HM(1,1),HM(1,2),phi(dado,HM(1,1),HM(1,2))
print*,HM(2,1),HM(2,2),phi(dado,HM(2,1),HM(2,2))
print*,HM(3,1),HM(3,2),phi(dado,HM(3,1),HM(3,2))
print*,HM(4,1),HM(4,2),phi(dado,HM(4,1),HM(4,2))
print*,HM(5,1),HM(5,2),phi(dado,HM(5,1),HM(5,2))

pause



 CALL cpu_time(tf)

 tt=tf-ti

 WRITE(6,FMT=*)'Tempo de máquina',tt,'segundos'


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!--------------------------------FIM------------------------------------------!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



CONTAINS
!---------------------------------------------------------
 FUNCTION Phi(dado,a,b)
 !Função de ajuste de dados. Datamisfit
  INTEGER,PARAMETER::SP = SELECTED_INT_KIND(r=4)
  INTEGER,PARAMETER::DP = SELECTED_REAL_KIND(8,8)
  INTEGER(KIND=SP):: i
  REAL(KIND=DP)::phi
  REAL(KIND=DP),INTENT(IN)::a,b
  REAL(KIND=DP),DIMENSION(:,:),INTENT(IN)::dado
 
  phi=0.0
  DO i=1,6
    phi=phi+(dado(i,2)-a*dado(i,1)-b)**2
  END DO 

 END FUNCTION Phi
!---------------------------------------------------------
 SUBROUTINE INITIALIZE(HMS,NVAR,fit,PVB,HM)
 INTEGER(8),INTENT(IN)::HMS,NVAR
 REAL(8),INTENT(INOUT),DIMENSION(:)::fit
 REAL(8),INTENT(INOUT),DIMENSION(:,:)::PVB,HM
 INTEGER::i,j,z
 REAL(8), DIMENSION(:,:), ALLOCATABLE::SUBS
 REAL(8)::a,b
 

 ALLOCATE(SUBS(1:NVAR,1:HMS))

 ! Procedimento para inicializar HM randomicamente
 PVB(1,1) = 0.5 !amin
 PVB(1,2) = 0.1  !bmin
 PVB(2,1) = 5.0  !amax
 PVB(2,2) = 3.0  !bmax

!print*,'Matriz PVB=' 
!print*,PVB(1,1),PVB(1,2)
!print*,PVB(2,1),PVB(2,2)



DO i=1,HMS !5
  DO j=1,NVAR !2
    HM(i,j)=randval(PVB(1,j),PVB(2,j))
  END DO
END DO


!pause
  !DO z=1,NVAR
  !  SUBS(z,i)=HM(i,z)
  !END DO

 DO i=1,HMS
    a=HM(i,1)
    b=HM(i,2)
    fit(i) = Phi(dado,a,b)
   ! print*,'fit(i)=',fit(i)
 END DO

 END SUBROUTINE INITIALIZE
!----------------------------------------------------------------------
SUBROUTINE UpdateHM(currentIteration,HMS,BestGen,fit,NCHV, HM, NewFit )

  REAL(8),INTENT(INOUT)::NewFit
  REAL(8),INTENT(INOUT),DIMENSION(:):: BestGen,fit,NCHV
  REAL(8),INTENT(INOUT),DIMENSION(:,:)::HM
  INTEGER(8),INTENT(IN):: currentIteration,HMS
  REAL(8)::BestFit,WorstFit
  INTEGER(8)::BestIndex,WorstIndex,i,j

IF(currentIteration==0) THEN
  BestFit=fit(1)
   DO i = 1,HMS
      IF( fit(i) < BestFit ) THEN
        BestFit = fit(i)
        BestIndex =i
      END IF
   END DO
  WorstFit=fit(1)
  WorstIndex =1
   DO i = 1,HMS
     IF( fit(i) > WorstFit ) THEN 
       WorstFit = fit(i)
       WorstIndex =i
     END IF
   END DO  
   ELSE
     IF (NewFit< WorstFit) THEN
      IF( NewFit < BestFit ) THEN
        DO j=1,NVAR
          HM(WorstIndex,j)=NCHV(j)
        END DO
         BestFit=NewFit
         BestGen=NCHV
         fit(WorstIndex)=NewFit
         BestIndex=WorstIndex
        ELSE
         DO j=1,NVAR 
           HM(WorstIndex,j)=NCHV(j)
         END DO
          fit(WorstIndex)=NewFit
          WorstFit=fit(1)
          WorstIndex =1
         DO i = 1,HMS
          IF( fit(i) > WorstFit ) THEN
           WorstFit = fit(i)
           WorstIndex =i
          END IF
         END DO
      END IF
    END IF 
END IF

IF(CURRENTITERATION/1000*1000==CURRENTITERATION) THEN
 PRINT*,BESTFIT,HM(:,1),HM(:,2)
 WRITE(2,FMT=*) BESTFIT,HM(:,1),HM(:,2)
END IF

END SUBROUTINE UpdateHM 

! function
!------------------------------------------
FUNCTION randval(Maxv,Minv)

 REAL(8)::randval,Maxv,Minv
 !REAL(8)::randval,RAND,Maxv,Minv
 !CALL RANDOM_NUMBER(RAND)

 randval=RAND()*(Maxv-Minv)+Minv

 END FUNCTION randval
!------------------------------------------
FUNCTION randint(Maxv,Minv)

 REAL(8)::RAND
 INTEGER(8)::Minv
 INTEGER(4)::randint,Maxv

 CALL RANDOM_NUMBER(RAND)

 randint=INT(RAND*(Maxv-Minv)+Minv + 0.5)
 

 END FUNCTION randint
!------------------------------------------
END PROGRAM BuscaHarmonica


