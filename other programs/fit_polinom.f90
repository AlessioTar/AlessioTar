PROGRAM fit_polinomiale
IMPLICIT NONE
REAL*8, ALLOCATABLE:: m(:),t(:),a(:)
INTEGER::i,r, grado, dim_a,j
REAL*8::summa,valore,s_r

r=0
OPEN(21, file='m_t.dat')
DO
  READ(21,*, END=10)
  r=r+1
END DO
10 WRITE(*,*) 'HO LETTO', r, 'RIGHE'
ALLOCATE(m(r),t(r))
REWIND(21)
DO i=1,r
  READ(21,*) m(i),t(i)
  !WRITE(*,*) i, m(i), t(i)
END DO

m=LOG10(m)
t=LOG10(t)

24 WRITE(*,*) 'DI CHE GRADO VUOI IL POLINOMIO?'
WRITE(*,*) 'INSERIRE 0 FA CHIUDERE IL PROGRAMMA'
READ(*,*) grado

IF(grado==0)THEN
  go to 25
END IF


dim_a=grado+1
ALLOCATE(a(dim_a))

CALL fit(m,t,r,a,dim_a)
WRITE(*,*) 'i coefficienti del fit M-T sono', (a(i), i=1,dim_a)
WRITE(*,*) '     '
summa=0.
DO i=1,r
   valore=0.
   DO j=1,dim_a
      valore=valore+a(j)*(m(i)**(j-1))
   END DO
   summa=summa+(t(i)-valore)**2
END DO
s_r=summa
WRITE(*,11)'al grado', grado, 'lo scarto quadratico vale', s_r
11 format(a,i3,x,a,e15.8)

!CALL fit(t,m,r,a,dim_a)
!WRITE(*,*) 'i coefficienti del fit T-M sono', (a(i), i=1,dim_a)
DEALLOCATE(a)
write(*,*) '    '
write(*,*) '    '
write(*,*) '    '
write(*,*) '    '

go to 24
25 continue
END PROGRAM fit_polinomiale


SUBROUTINE fit(x,y,dim,a,dim2)
  IMPLICIT NONE
  INTEGER::dim,i,j,k,dim2
  REAL*8::mat(dim2,dim2),c(dim2),somma
  REAL*8::x(dim),y(dim),a(dim2)
  mat=0
  c=0
  DO i=1,dim2
    DO j=1,dim2
      somma=0.
      DO k=1,dim
        somma=somma+x(k)**(i+j-2)
      END DO
      mat(i,j)=somma
    END DO
  END DO

  DO i=1,dim2
    somma=0.
    DO k=1,dim
      somma=somma+y(k)*x(k)**(i-1)
    END DO
    c(i)=somma
  END DO

!SERVE PER STAMPARE LA MATRICE DELLE SOMMATORIE
!DO i=1,2
! WRITE(*,*) (mat(i,j), j=1,2), c(i)
!END DO

CALL gauss(mat,c,dim2,a)

END SUBROUTINE fit

SUBROUTINE gauss(a,c,nd,x)
  IMPLICIT NONE
  INTEGER::nd
  REAL*8::a(nd,nd), c(nd), x(nd),fakt,sum
  INTEGER::i,j,k


  ! 'elimino le variabili '
  DO i=1,nd-1 !SCEGLIE LA VARIABILE DA ELIMIUNARE
    DO j=i+1,nd !SCEGLIE LA RIGA DA CUI ELIMINARE
      fakt=a(j,i)/a(i,i)
      DO k=1,nd
        a(j,k)=a(j,k)-a(i,k)*fakt
      END DO
      c(j)=c(j)-c(i)*fakt        !QUESTO è FUORI DA CICLO DI PRIMA PERCHè C() è COMPOSTO SOLO DA UNA COLONNA QUINDI NON SERVE CHE SIA DENTRO IL DO POICHè FAREBBE TRE VOLTER LA STESSA COSA
    END DO
  END DO
!  READ(*,*)

  !ricavo le soluzioni
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
