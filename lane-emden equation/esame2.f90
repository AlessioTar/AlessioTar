MODULE sistema
  IMPLICIT NONE
  CONTAINS

  SUBROUTINE dxdy(x,y1,y2,neq,n,yp)

    REAL*8,INTENT(IN)::x,y1,y2
    INTEGER,INTENT(IN)::neq,n

    REAL*8,INTENT(OUT)::yp(neq)
    !inizializzo il sistema da risolvere
    !x rappresenta la variabile di integrazione xi
    !y1 è theta, y2 è eta
    !yp(1) corrisponde al rhs della prima riga del sistema e yp(2) alla seconda

    yp(1)=y2               !dTheta/dXi=Eta
    yp(2)=-y1**n-(2/x)*y2  !dEta/dXi=-Theta**n-2*eta/Xi

  END SUBROUTINE dxdy
END MODULE sistema

SUBROUTINE RK4(neq,h,xold,y1old,y2old,y1_ris,y2_ris,n)
  USE sistema
  IMPLICIT NONE

    REAL*8, INTENT(IN) :: y1old,y2old,h,xold
    REAL*8, INTENT(OUT) :: y1_ris,y2_ris
    REAL*8 :: k1(2), k2(2), k3(2), k4(2),y1new,y2new,xnew

    INTEGER :: i,neq,n

    !primo stadio
    CALL dxdy(xold, y1old,y2old,neq,n,k1)

    !secondo stadio
    y1new=y1old+0.5*h*k1(1)
    y2new=y2old+0.5*h*k1(2)

    xnew=xold+0.5d0*h
    CALL dxdy(xnew, y1new,y2new,neq,n,k2)

    !terzo stadio
    y1new=y1old+0.5*h*k2(1)
    y2new=y2old+0.5*h*k2(2)

    xnew=xold+0.5d0*h
    CALL dxdy(xnew, y1new,y2new,neq,n,k3)

    !quarto stadio
    y1new=y1old+h*k3(1)
    y2new=y2old+h*k3(2)
    xnew=xold+h
    CALL dxdy(xnew, y1new,y2new, neq,n,k4)

    !formula finale che fornisce il risultato atteso
    y1_ris=y1old+h*(k1(1)+2.0d0*(k2(1)+k3(1))+k4(1))/6.0d0
    y2_ris=y2old+h*(k1(2)+2.0d0*(k2(2)+k3(2))+k4(2))/6.0d0
END SUBROUTINE RK4

SUBROUTINE zeri(left,right,theta,eta,n_val,xi_z)
  IMPLICIT NONE
  REAL*8::left,right,theta,eta,h,med,test,etaris,int
  REAL*8,INTENT(OUT)::xi_z
  INTEGER::n_val,conta=0
  int=ABS(right-left)
    DO WHILE(int>0.000000000000001d0)   !finchè la mia coordinata verticale non è vicina a 0 per 7 cifre significative continuo
    med=0.5*(right+left)
    conta=conta+1
    CALL rk4(2,int,left,theta,eta,test,etaris,n_val)
    IF(test==0.d0) THEN
      exit
    ELSEIF(test<0.d0) THEN
       right=med
       !theta=test
     ELSE
       left=med
       theta=test
!CAMBIARE ANCHE LA ETA!
     END IF
     !READ(*,*)
     int=ABS(right-left)
     PRINT*, int
   END DO
   print*, left,right
   xi_z=(left+right)/2.
END SUBROUTINE zeri


REAL *8 FUNCTION theta_n0(xi)
	IMPLICIT NONE
	REAL*8::xi
	theta_n0 = 1.d0 - (xi**2)/6.d0
END FUNCTION theta_n0

REAL *8 FUNCTION eta_n0(xi)
	IMPLICIT NONE
	REAL*8::xi
	eta_n0=-xi/3.d0
END FUNCTION eta_n0

REAL *8 FUNCTION theta_n1(xi)
	IMPLICIT NONE
	REAL*8::xi
	theta_n1=SIN(xi)/xi
END FUNCTION theta_n1

REAL *8 FUNCTION eta_n1(xi)
	IMPLICIT NONE
	REAL*8::xi
	eta_n1=(xi*COS(xi)-SIN(xi))/(xi**2)
END FUNCTION eta_n1

REAL *8 FUNCTION mass_profile(xi,et)
	IMPLICIT NONE
	REAL*8::xi,et
  mass_profile=-(xi**2)*et
END FUNCTION mass_profile



PROGRAM laneemden
  use sistema
  IMPLICIT NONE
  !inizializzo le condizioni su cui integrare
  REAL*8::xi_min,xi_max,h,h1,intervallo,xi_zero(3,2),xi_aux,theta_aux,eta_aux,an(3)
  REAL*8,ALLOCATABLE::theta(:),xi(:),eta(:),er_theta0(:),er_theta1(:),er_eta0(:),er_eta1(:),er_mas0(:),er_mas1(:)
  INTEGER::np,n(3),i,j,k,n_val,conta,l
  REAL*8,EXTERNAL::theta_n0,theta_n1,eta_n0,eta_n1,mass_profile

  an(1)=(sqrt(6.d0))
  an(2)=4*ATAN(1.d0)
  an(3)=0.d0      !soluzione analitica ignota, non è possibile stimare un errore
  n(1)=0
  n(2)=1
  n(3)=3
  xi_zero(:,:)=0.d0

  xi_min=0.d0
  xi_max=10.d0
  open(40,file='ris_n0.txt')
  open(41,file='ris_n1.txt')
  open(43,file='ris_n3.txt')
  open(80,file='mass_profile_n0.txt')
  open(81,file='mass_profile_n1.txt')
  open(83,file='mass_profile_n3.txt')


  DO j=1,2
    IF(j==1) THEN
      np=200
    ELSE
      np=2000
    END IF
    ALLOCATE(theta(np),eta(np),xi(np),er_theta0(np),er_theta1(np),er_eta0(np),er_eta1(np),er_mas0(np),er_mas1(np))
    h=(xi_max-xi_min)/np
    !inizializzo xi
    xi(1)=xi_min
    DO i=2,np
      xi(i)=xi(i-1)+h
    END DO
    !inizializzo eta e theta con le condizioni iniziali note e quelle fornite dallo sviluppo di taylor
    DO k=1,3
      n_val=n(k)
      PRINT *, 'CALCOLATO CON N UGUALE A',n_val, np,'punti'
      eta(1)=0.d0
      theta(1)=1.d0
      eta(2)=-(xi(2))/3.d0+(n_val*(xi(1))**3)/30.d0
      theta(2)=1.d0-((xi(2))**2)/6.d0+(n_val*(xi(2))**4)/120.d0
if(np==200)then
if(n_val==0)then
      WRITE(40,*)'theta         eta          phi'
      WRITE(40,*)xi(1),theta(1),eta(1),mass_profile(xi(1),eta(1))
      WRITE(40,*)xi(2),theta(2),eta(2),mass_profile(xi(2),eta(2))

elseif(n_val==1)then
      WRITE(40,*)'theta         eta          phi'
      WRITE(41,*)xi(1),theta(1),eta(1),mass_profile(xi(1),eta(1))
      WRITE(41,*)xi(2),theta(2),eta(2),mass_profile(xi(2),eta_n1(xi(2)))

elseif(n_val==3)then
      WRITE(43,*)xi(1),theta(1),eta(1),mass_profile(xi(1),eta(1))
      WRITE(43,*)xi(2),theta(2),eta(2),mass_profile(xi(2),eta(2))


end if
end if
      conta=0
      DO i=3,np
        CALL rk4(2,h,xi(i-1),theta(i-1),eta(i-1),theta(i),eta(i),n_val)
        if(np==200)then
        if(n_val==0)then
          WRITE(40,*)xi(i),theta(i),eta(i),mass_profile(xi(i),eta(i))
        elseif(n_val==1)then
          WRITE(41,*)xi(i),theta(i),eta(i),mass_profile(xi(i),eta(i))
        elseif(n_val==3)then
          WRITE(43,*)xi(i),theta(i),eta(i),mass_profile(xi(i),eta(i))
        end if
        end if
        !PRINT *,i, 'xi',xi(i),theta(i),'theta',eta(i),'eta'
        xi_aux=xi(i-1)
        theta_aux=theta(i-1)
        eta_aux=eta(i-1)

        IF(conta==0)THEN
          IF(theta(i)*theta(i-1)<0.) THEN
            CALL zeri(xi(i-1),xi(i),theta_aux,eta_aux,n_val,xi_zero(k,j))
            conta=1
          END IF
        END IF

      END DO
      open(50,file='errori.txt')
      open(51,file='errori 10.txt')
      DO l=1,np
        IF(n_val==0)THEN
          er_theta0(l)=ABS(theta(l)-theta_n0(xi(l)))
          er_eta0(l)=ABS(eta(l)-eta_n0(xi(l)))
          er_mas0(l)=ABS(mass_profile(xi(l),eta(l))-mass_profile(xi(l),eta_n0(xi(l))))
        ELSEIF(n_val==1)THEN
          er_theta1(l)=ABS(theta(l)-theta_n1(xi(l)))
          er_eta1(l)=ABS(eta(l)-eta_n1(xi(l)))
          er_mas1(l)=ABS(mass_profile(xi(l),eta(l))-mass_profile(xi(l),eta_n1(xi(l))))


          if(np==2000)then
            WRITE(51,*) er_theta0(l),er_eta0(l),er_theta1(l),er_eta1(l),er_mas0(l),er_mas1(l)
          elseif(np==200)then
            WRITE(50,*) er_theta0(l),er_eta0(l),er_theta1(l),er_eta1(l),er_mas0(l),er_mas1(l)
          end if
        ELSEIF(n_val==3)then
          exit
        END IF
      END DO
      an(3)=xi_zero(k,j)     !soluzione analitica ignota, non è possibile stimare un errore
      PRINT*,'xi zero', xi_zero(k,j),'errore', ABS(xi_zero(k,j)-an(k))     !nel terzo caso apparirà zero perchè non è noto il risultato analitico
      print *, '    '

    END DO
    DEALLOCATE(theta,eta,xi,er_theta0,er_theta1,er_eta0,er_eta1,er_mas0,er_mas1)

  END DO
  close(40)
  close(41)
  close(43)
  close(50)
  close(51)

  WRITE (*,*) "PREMERE INVIO PER CONCLUDERE L'ESECUZIONE"
  READ(*,*)
END PROGRAM laneemden
