PROGRAM deleting
!QUESTO PROGRAMMA FA L'ELIMINAZIONE DI GAUS PER I SISTEMI LINEARI, LI RISOLVE E VERIFICA LA SOLUZIONE
IMPLICIT NONE
INTEGER, PARAMETER::nd=3
REAL*8::a(nd,nd), c(nd), x(nd), a_copy(nd,nd),sum, c_copy(nd)
INTEGER::i,j

OPEN(21,file='matrice_elim.dat')
DO i=1,nd
  READ(21,*) (a(i,j),j=1,nd), c(i) 
END DO
PRINT *,'MATRICE INIZIALE'

DO i=1,nd
  WRITE(*,11) (a(i,j),j=1,nd)
  11 FORMAT(f6.3,x,f6.3,x,f6.3)
END DO

PRINT*, 'TERMINI NOTI'
DO i=1,nd
  WRITE(*,12) c(i)
  12 FORMAT(f6.3,x,f6.3,x,f6.3)
END DO
a_copy=a
c_copy=c

CALL gauss(a,c,nd,x)



print *, 'le soluzioni sono'
DO i=1,nd
  !READ(*,*)
  WRITE(*,77) i, x(i)
  77 format(x,i2,f10.4)
END DO
!a questo punto voglio verificare se il mio programma funziona con a_copy e c_copy
print *, 'verifica '
do i=1,nd
  sum=0
  do j=1,nd
    sum=sum+a_copy(i,j)*x(j)
  end do
  WRITE(*,10) i, sum, c_copy(i)
  10 FORMAT(i2,x,f9.5,x,f9.5)
  !READ(*,*)
end do

END PROGRAM deleting

SUBROUTINE gauss(a,c,nd,x)
  IMPLICIT NONE
  INTEGER::nd
  REAL*8::a(nd,nd), c(nd), x(nd),fakt,sum
  INTEGER::i,j,k


  print *, 'elimino le variabili '
  DO i=1,nd-1 !SCEGLIE LA VARIABILE DA ELIMIUNARE
    DO j=i+1,nd !SCEGLIE LA RIGA DA CUI ELIMINARE
      fakt=a(j,i)/a(i,i)
      DO k=1,nd
        a(j,k)=a(j,k)-a(i,k)*fakt
      END DO
      c(j)=c(j)-c(i)*fakt        
    END DO
  END DO
!  READ(*,*)

  print *, 'ricavo le soluzioni '
  x(nd)=c(nd)/a(nd,nd)
  DO i=nd-1,1,-1
    sum=0
    DO j=i+1,nd
      sum=sum+a(i,j)*x(j)
    END DO
    x(i)=(c(i)-sum)/a(i,i)
  END DO
!  READ(*,*)

end subroutine gauss
