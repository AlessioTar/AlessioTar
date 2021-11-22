PROGRAM sorting
  IMPLICIT NONE
  INTEGER::conta, nd, i, i_col, i_opt

  REAL*8, ALLOCATABLE:: col1(:), col2(:)

  OPEN(21,file='spline.dat')

  conta=0
  DO
     READ(21,*, END=66)
     conta=conta+1
  END DO
66 PRINT *, "ho letto ", conta, " righe"
  nd=conta
  REWIND(21)
  ALLOCATE(col1(nd), col2(nd))
  DO i=1,nd
     READ(21,*) col1(i), col2(i)
     WRITE(*,*) col1(i), col2(i)
  END DO

  55 PRINT *, "scegli il vettore su cui fare l'ordinamento"
  PRINT *, " scrivi 1 se colonna 1"
  PRINT *, " scrivi 2 se colonna 2"
  !READ(*,*) i_col

  

  555 PRINT *, "scegli il tipo di ordinamento"
  PRINT *, " scrivi 1 se crescente"
  PRINT *, " scrivi 2 se decrescente"
  !READ(*,*) i_opt


  CALL sort_colonne(col1,col2,nd,1,1)
  DO i=1,nd
     WRITE(*,77) i, col1(i), col2(i)
     77 FORMAT(x,i4,2(x,f10.4))
  END DO
read(*,*)
END PROGRAM sorting

SUBROUTINE sort_colonne(vect1, vect2, nd, i_col, i_opt)
  IMPLICIT NONE
  INTEGER:: nd,i_opt, i_col, i, posizione,j
  REAL*8:: valore_rif, provv
  REAL*8:: vect1(nd), vect2(nd)
  REAL*8:: aux1(nd), aux2(nd)

  IF(i_col==1) THEN
     aux1=vect1
     aux2=vect2
  ELSE
     aux1=vect2
     aux2=vect1
  END IF

  DO i=1,nd-1
     valore_rif=aux1(i)
     posizione=i
     DO j=i+1,nd
        SELECT CASE(i_opt)
        CASE(1)
           IF(aux1(j)<valore_rif) THEN
              valore_rif=aux1(j)
              posizione=j
           END IF
        CASE(2)
            IF(aux1(j)>valore_rif) THEN
              valore_rif=aux1(j)
              posizione=j
           END IF
        END SELECT
     END DO
     aux1(posizione)=aux1(i)
     aux1(i)=valore_rif
     provv=aux2(posizione)
     aux2(posizione)=aux2(i)
     aux2(i)=provv
  END DO

  SELECT CASE(i_col)
  CASE(1)
     vect1=aux1
     vect2=aux2
  CASE(2)
     vect1=aux2
     vect2=aux1
  END SELECT

END SUBROUTINE  sort_colonne
