PROGRAM BuscaHarmonica

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!DEFINIÇÃO DE VARIÁVEIS GLOBAIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
 IMPLICIT NONE
   INTEGER,PARAMETER::SP = SELECTED_INT_KIND(r=12)
   INTEGER,PARAMETER::DP = SELECTED_REAL_KIND(12,14)
   INTEGER::NVAR,NG,NH,HMS
   INTEGER::currentIteration,MaxItr,I,index
   REAL(KIND=DP)::PARmin,PARmax,bwmin,bwmax,HMCR,PAR,newFitness,BestFit,WorstFit
   REAL(KIND=DP)::coef,RES,RAN,pvbRan,ti,tf,tt
   INTEGER(KIND=DP)::pp
   REAL(KIND=DP), ALLOCATABLE, DIMENSION(:)::fit,NCHV,BestGen,BW,gx
   REAL(KIND=DP), ALLOCATABLE, DIMENSION(:,:)::HM,PVB



 !Entradas
 NVAR=10
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

 ALLOCATE(HM(1:HMS,1:NVAR),PVB(1:NVAR,1:2))
 ALLOCATE(fit(1:HMS),NCHV(1:NVAR),BestGen(1:NVAR),BW(1:NVAR),gx(1:NG))

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
   newFitness = Fitness(NCHV)
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
 FUNCTION Fitness(sol)

  REAL(8),DIMENSION(1:),INTENT(IN)::sol
  REAL(8)::Fitness

  Fitness=-(sol(1)+sol(2))
!  Fitness= -(EXP(1/2*(sol(1)+sol(2)-25)**2)+sin(4*sol(1)-3*sol(2))**4+1/2*(2*sol(1)+sol(2)-10)**2)

 END FUNCTION Fitness
!---------------------------------------------------------
 SUBROUTINE INITIALIZE()

 INTEGER::i,j,z
 REAL(8), DIMENSION(:), ALLOCATABLE::SUBS

 ALLOCATE(SUBS(1:NVAR))

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
    fit(i) = Fitness(SUBS)
END DO

 END SUBROUTINE INITIALIZE
!----------------------------------------------------------------------
 SUBROUTINE UpdateHM( NewFit )

  REAL(8)::NewFit
  INTEGER::BestIndex,WorstIndex,I,J

IF(currentIteration==0) THEN

BestFit=fit(1)
DO I = 1,HMS

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

DO J=1,NVAR

HM(WorstIndex,J)=NCHV(J)

END DO

BestFit=NewFit
BestGen=NCHV
fit(WorstIndex)=NewFit
BestIndex=WorstIndex

ELSE

DO J=1,NVAR

  HM(WorstIndex,J)=NCHV(J)

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

END IF

end SUBROUTINE UpdateHM ! function
!------------------------------------------
function randval(Maxv,Minv)

 REAL(8)::randval,RAND,Maxv,Minv

 CALL RANDOM_NUMBER(RAND)

 randval=RAND*(Maxv-Minv)+Minv

end function randval
!------------------------------------------
function randint(Maxv,Minv)

 REAL(8)::RAND
 INTEGER::randint,Maxv,Minv

 CALL RANDOM_NUMBER(RAND)

 randint=INT(RAND*(Maxv-Minv)+Minv + 0.5)

end function randint
!------------------------------------------
END PROGRAM BuscaHarmonica


