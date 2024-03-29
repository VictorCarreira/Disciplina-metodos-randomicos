PROGRAM BuscaHarmonica

    !-----------------Sumário das Variáveis-----------------!
    !HS: Harmony memory                                     !
    !HMS: Harmony memory size (number of measured points)   !
    !HMCR: Harmony consideration rate                       !
    !PAR: Pitch adjusting rate                              !
    !bw : termination criterion                             !
    !MaxItr: Maximal iteraction                             !
    !NVAR: Number of parameters                             !
    !-------------------------------------------------------!


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!DEFINIÇÃO DE VARIÁVEIS GLOBAIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
 IMPLICIT NONE
   INTEGER,PARAMETER::SP = SELECTED_INT_KIND(r=12)
   INTEGER,PARAMETER::DP = SELECTED_REAL_KIND(12,14)
   INTEGER::NVAR,NG,NH,HMS
   INTEGER::currentIteration,MaxItr,i,j,index
   REAL(KIND=DP)::PARmin,PARmax,bwmin,bwmax,HMCR,PAR,newFitness,BestFit,WorstFit
   REAL(KIND=DP)::coef,RES,RAN,pvbRan,ti,tf,tt
   INTEGER(KIND=DP)::pp
   REAL(KIND=DP), ALLOCATABLE, DIMENSION(:)::fit,NCHV,BestGen,BW,gx
   REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:)::HM,PVB,dado

 
OPEN(1,FILE='entrada.txt')
OPEN(2,FILE='saida.txt')



 !Dimensionando o problema:
 NVAR=10 !Dimensão horizontal
 NG=1
 NH=0
 MaxItr=10000
 HMS=5
 HMCR=0.6
 PARmin=0.4
 PARmax=0.85
 bwmin=4
 bwmax=0.01

 CALL cpu_time(ti)

 ALLOCATE(HM(1:HMS,1:NVAR),PVB(1:NVAR,1:2),dado(6,2))
 ALLOCATE(fit(1:HMS),NCHV(1:NVAR),BestGen(1:NVAR),BW(1:NVAR),gx(1:NG))


 DO i=1,6
  READ(1,FMT=*)dado(i,1),dado(i,2) 
 END DO  


 CALL initialize()

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
   newFitness = Phi(dado,1.3,0.2)
   CALL UpdateHM( newFitness )
   currentIteration=currentIteration+1
  END DO
END DO


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
  REAL(8),DIMENSION(:,:),INTENT(IN)::dado
  REAL(4)::phi,a,b

  phi=0
  DO i=1,6
    phi=phi+(dado(i,2)-a*dado(i,1)+b)**2
  END DO 

 END FUNCTION Phi
!---------------------------------------------------------
 SUBROUTINE INITIALIZE()

 INTEGER::i,j,z
 REAL(8), DIMENSION(:), ALLOCATABLE::SUBS

 ALLOCATE(SUBS(1:NVAR))

 ! Tamanho das variáveis
 PVB(1,1) = 1.0 
 PVB(1,2) = 4.0
 PVB(2,1) = 5.0
 PVB(2,2) = 10.0

DO i=1,HMS
  DO j=1,NVAR
    HM(i,j)=randval(PVB(j,1),PVB(j,2))
  END DO 

  DO z=1,NVAR
    SUBS(z)=HM(i,z)
  END DO
    fit(i) = Phi(SUBS,1,2)
END DO

 END SUBROUTINE INITIALIZE
!----------------------------------------------------------------------
SUBROUTINE UpdateHM( NewFit )

  REAL(8)::NewFit
  INTEGER::BestIndex,WorstIndex,i,j

IF(currentIteration==0) THEN

BestFit=fit(1)
DO i = 1,HMS

IF( fit(i) < BestFit ) THEN

BestFit = fit(i)

BestIndex =i

end IF

end DO

WorstFit=fit(1)
WorstIndex =1
DO i = 1,HMS

if( fit(i) > WorstFit ) THEN

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

 END IF ! main if

END IF

IF(CURRENTITERATION/1000*1000==CURRENTITERATION) THEN

 PRINT*,BESTFIT,BestGen(1),BestGen(2)
 WRITE(2,FMT=*) BESTFIT,BestGen(1), BestGen(2)

END IF

END SUBROUTINE UpdateHM ! function
!------------------------------------------
FUNCTION randval(Maxv,Minv)

 REAL(8)::randval,RAND,Maxv,Minv

 CALL RANDOM_NUMBER(RAND)

 randval=RAND*(Maxv-Minv)+Minv

 END FUNCTION randval
!------------------------------------------
FUNCTION randint(Maxv,Minv)

 REAL(8)::RAND
 INTEGER::randint,Maxv,Minv

 CALL RANDOM_NUMBER(RAND)

 randint=INT(RAND*(Maxv-Minv)+Minv + 0.5)

 END FUNCTION randint
!------------------------------------------
END PROGRAM BuscaHarmonica


