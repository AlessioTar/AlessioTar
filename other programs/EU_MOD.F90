MODULE RHS
  IMPLICIT NONE
  CONTAINS
  SUBROUTINE dxdy(neq,y,yp)          !definisce la mia equazione diff
    REAL*8::y(neq),yp(neq)
    INTEGER::neq
    yp(1)=y(2)
    yp(2)=-y(1)
  END SUBROUTINE dxdy
END MODULE RHS

MODULE ODE_SOLVER
  USE RHS
  IMPLICIT NONE
  CONTAINS
  SUBROUTINE save_result(fname,npoint,neq,x,y)
    IMPLICIT NONE
    CHARACTER(len=*)::fname
    REAL*8::x(neq,npoint),y(neq,npoint)
    INTEGER::npoint,i,neq
    OPEN(40,file=fname)

    DO i=0, npoint
      WRITE(40,'(3(1pe20.12))') x(i), y(i)
    END DO

    CLOSE(40)   !salva un file .txt con i risultati
  END SUBROUTINE save_result

  SUBROUTINE fwd_euler(neq,h,yold,ynew)            !applica eulero in avanti
    IMPLICIT NONE
    REAL*8::h,yold(neq),ynew(neq),yp(neq)
    INTEGER::neq
    CALL dxdy(neq,yold,yp)
    ynew(i)=yold(i)+h*yp(i)
  END SUBROUTINE


END MODULE ODE_SOLVER





PROGRAM eu
  USE ODE_SOLVER
  USE RHS
  IMPLICIT NONE
  INTEGER::i,neq 
  INTEGER, PARAMETER::npunti=50
  REAL*8,DIMENSION(0:npunti):: x,y

  REAL*8::xmin,xmax,h
  !inizializzo la griglia
  xmin=0.d0
  xmax=4.d0*atan(1.0d0)
  h=(xmax-xmin)/npunti
  !inizializzo x
  DO i=0,npunti
    x(i)=xmin+i*h
  END DO
  !inizializzo y
  y(0)=1.d0
  DO i=1,npunti
    CALL fwd_euler(neq,h,y(:,i-1),y(:,i))
    PRINT*, x(i),y(i)

  END DO

  CALL save_result('eulero_mod.txt',npunti,neq,x,y)




END PROGRAM eu
