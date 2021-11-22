PROGRAM integral
IMPLICIT NONE
REAL*8,EXTERNAL::primitiva,integranda
REAL*8::a,b,esatto,integrale
INTEGER::i
PRINT *, 'DAMMI GLI ESTREMI DI INTEGRAZIONE'
READ(*,*) a,b
esatto=primitiva(b)-primitiva(a)

CALL quadgauss(integranda,a,b,integrale)
!!!!!!!serve un print dei risultati

END PROGRAM integral
!questo metodo è comodo perchè con solo la stima di 4 punti è possibile risolvere integrali molto complessi risparmiando nel numero di cicli necessari
SUBROUTINE QUADGAUSS(f,a,b,risultato)
  IMPLICIT NONE
  REAL*8,EXTERNAL:: f
  REAL*8::a,b,risultato,summa
  REAL*8:: c(0:1),x(0:1) ,xd(0:1)
  INTEGER::i


  c(0)=1.
  c(1)=1.
  x(0)=-1./SQRT(3.)
  x(1)=1./SQRT(3.)
  DO i=0,1
    x(i)=0.5*((b-a)*xd(i))
  END DO
  summa=0.
  DO i=0,1
    summa=summa+c(i)*f(x(i))
  END DO
  risultato=summa*(b-a)*0.5


END SUBROUTINE  QUADGAUSS




REAL*8 FUNCTION integranda(x)
IMPLICIT NONE
REAL*8::x
integranda=1./(SQRT(9-x**2))
END FUNCTION integranda

REAL*8 FUNCTION primitiva(x)
IMPLICIT NONE
REAL*8::x
primitiva=ASIN(x/3)
END FUNCTION primitiva
