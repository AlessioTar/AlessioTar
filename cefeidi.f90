MODULE dati
  IMPLICIT NONE
  REAL*8:: om,oa
  SAVE
END MODULE dati

PROGRAM cefeidi
  use dati
  IMPLICIT NONE
  !variabili riguardanti i 19 cataloghi di galassie
  INTEGER::ngal,nval,i,j,bin,k
  REAL*8::val,ris,incr=0.01,h0,errh0,num,den
  REAL*8::pmed,magmed=0.d0,mv,distanza,errpe,errmapp,errmv,errdist
  REAL*8,ALLOCATABLE::gg(:),mag(:),matrice(:,:),der(:),c(:),curva(:),vel(:),errvel(:),h(:),errh(:)
  REAL*8,ALLOCATABLE::tot_p(:),tot_erp(:),tot_mapp(:),tot_ermap(:),tot_mv(:),tot_ermv(:),tot_dist(:),tot_erdist(:)

  CHARACTER(LEN=16)::filename
  CHARACTER(LEN=7),ALLOCATABLE::nome(:)

  !variabili riguardanti la stima delle costanti c1 e c2
  REAL*8,ALLOCATABLE::pe(:),magass(:)
  REAL*8::cof(2)
  INTEGER::ncat

  !variabili riguardanti la stima di OmegaM e OmegaA
  REAL*8, EXTERNAL:: integranda,integranda2
  REAL*8:: x1,x2, inte1, inte2, risultato

!tramite un catalogo noto stimo le costanti c1 e c2 che mi serviranno
!in seguito per calcolare la magnitudine assolata delle cefeidi in esame
  OPEN(13,file='ceph_catalog.txt',status='old')
  ncat=0
  DO
     READ(13,*, END=88)
     ncat=ncat+1
  END DO
  88 CONTINUE
  REWIND(13)
  ALLOCATE(pe(ncat),magass(ncat))
  DO i=1,ncat
    READ(13,*) pe(i),magass(i)
  END DO
  CLOSE(13)
  CALL fit(LOG10(pe),magass,ncat,cof)


!leggo i nomi delle galassie dal file gal_vel
  OPEN(11,file='gal_vel.txt',status='old')
  ngal=0
  DO
     READ(11,*, END=66)
     ngal=ngal+1
  END DO
  66 continue
  REWIND(11)
  ALLOCATE(nome(ngal),vel(ngal),errvel(ngal))
  ALLOCATE(tot_p(ngal),tot_erp(ngal),tot_mapp(ngal),tot_ermap(ngal),&
          &tot_mv(ngal),tot_ermv(ngal),tot_dist(ngal),tot_erdist(ngal))
  DO i=1,ngal
    READ(11,*) nome(i),vel(i),errvel(i)
  END DO
  CLOSE(11)

  ALLOCATE(h(ngal),errh(ngal))
  DO i=1,ngal
    filename(1:5)='ceph_'
    filename(6:12)=nome(i)
    filename(13:16)='.txt'
    OPEN(12,file=filename,status='old')
    nval=0
    DO
       READ(12,*, END=77)
       nval=nval+1
    END DO
    77 CONTINUE
    REWIND(12)
    ALLOCATE(gg(nval),mag(nval))
    DO j=1,nval
      READ(12,*) gg(j),mag(j)
    END DO
    CLOSE(12)
    CALL sort_colonne(gg,mag,nval)                         !odrino i dati in base al giorno giuliano
    DO j=1,nval
    END DO
    ALLOCATE(matrice(nval,nval),c(nval))
    CALL mat_gen(gg,mag,nval,matrice,c)                    !genero la matrice tridiagonale
    ALLOCATE(der(nval))
    CALL gauss(matrice,nval,c,der)                         !col metodo di gauss risolvo il sistema associato alla matrice
    val=MINVAL(gg)
    bin=INT((MAXVAL(gg)-MINVAL(gg))/incr)

    ALLOCATE(curva(bin))
    magmed=0.
    DO k=1,bin
      CALL interp(val,gg,mag,der,nval,ris)                 !trovo il valore della magnitidine corrispondente a un certo giorno giuliano
      curva(k)=ris

      val=val+incr                                         !incremento di 0.05 alla il giorno giuliano così da ricostruire la curva di luce
    END DO

    CALL periodi_magnitudini(curva,bin,incr,pmed,errpe,magmed,errmapp)


    !a questo punto sono in grado di calcolare la magnitudine assoluta, quindi di conseguenza le distanze delle rispettive galassie
    mv=cof(2)*LOG10(pmed)+cof(1)
    errmv=(ABS(cof(2))/(LOG(10.)*pmed))*(errpe/pmed)

    distanza=((magmed-mv)*0.2)+1.d0
    distanza=(10**distanza)/10**6                          !così la distanza è espressa in Mpc
    errdist=ABS(errmapp+errmv)*2d0*LOG(10.d0)*10.d0**((magmed-mv)*0.2d0)
    errdist=errdist/10.d0**6                               !normalizzo l'errore in Mpc
    !calcolo la costante di hubble per ogni galassia e in seguito facendone la media pesata sugli errori ricavo una stima per il valore reale
    h(i)=vel(i)/distanza
    errh(i)=SQRT((errvel(i)**2/distanza**2)+((vel(i)*errdist)**2/distanza**4))

    tot_p(i)=pmed
    tot_erp(i)=errpe
    tot_mapp(i)=magmed
    tot_ermap(i)=errmapp
    tot_mv(i)=mv
    tot_ermv(i)=errmv
    tot_dist(i)=distanza
    tot_erdist(i)=errdist

!SCOMMENTARE QUESTA REGIONE DI CODICE PER LEGGEE I DATI DELLE SINGOLE STELLA SULLO SCHERMO, ALTRIMENTI APRIRE IL FILE output_data.txt
!__________________________________________________________________________________________________
    WRITE(*,*) i
    WRITE(*,*) 'dati relativi alla galassia ',filename
    WRITE (*,80)  'periodo medio               ', pmed,'+/-',ABS(errpe),' giorni giuliani'
    WRITE (*,80)  'magnitudine apparente media ',magmed,'+/-',ABS(errmapp),' mag'
    80 format(x,a,f8.4,a,f8.5,a)
    WRITE (*,80)  'magnitudine assoluta        ',mv,'+/-',ABS(errmv),' mag'
    WRITE (*,90) filename(1:12),' dista ',distanza,'+/-',ABS(errdist), ' Mpc'
    90 format (x,a,a,9(x),f8.4,a,f8.4,a)
    WRITE (*, 100) 'stima di H0                ',h(i),'+/-',ABS(errh(i)),'km/s/Mpc'
    100 format(x,a,x,f8.4,a,f8.4,x,a)
    WRITE(*,*) '   '
!__________________________________________________________________________________________________
    DEALLOCATE(gg,mag,matrice,c,der,curva)
    !READ(*,*)                                 scommentare per vedere i dati di ogni galassia comparire uno alla volta
  END DO

  num=0.
  den=0.
  DO i=1,ngal
   num=num+(h(i)/(errh(i)**2))
   den=den+(1./errh(i)**2)
 END DO
 h0=num/den

 errh0=SQRT(1./den)

!__________________________________________________________________________________________________
 WRITE(*,*) '   '
 WRITE(*,*) 'STIMA GENERALE DELLA COSTANTE DI HUBBLE'
 WRITE(*,*) h0,'+/-',errh0
!__________________________________________________________________________________________________
 CALL save_txt_data(ngal,nome,tot_p,tot_erp,tot_mapp,tot_ermap,tot_mv,tot_ermv,tot_dist,tot_erdist,h,errh,h0,errh0)

 !ora con la stima di H0, avendo nota l'età dell'universo,
 !ricavo dei possibili valori per la funzione E(z)
 x1=0.d0                                                   !estremi di integrazione x1 e x2
 x2=1.d0
 oa=0.d0
 om=1.d0
OPEN(15,file='omega_val.txt')
 WRITE(15,*)'  '
 WRITE(15,*) "l'eta' nota dell'universo e' 13.82 Gyr"
 WRITE(15,*) "valori di OmegaA ed OmegaM compatibili con l'eta' dell'universo ottenuti con il vincolo della somma pari a 1"
 WRITE(15,*) "stima eta'  OmegaA  OmegaM  OmegaTOT"

 do i=1,100
   CALL quadgauss(integranda,x1,x2,inte1)
   CALL quadgauss(integranda2,x1,x2,inte2)
   risultato=978.d0*(inte1+inte2)/h0
   if(risultato>=13.69.and.risultato<=13.95) then       !i valori accettati di discostano dal valore reale di +/- 1%
     WRITE(15,120) risultato,oa,om,oa+om
   end if
   oa=oa+0.01
   om=om-0.01
 end do
 WRITE(15,*) '  '
 WRITE(15,*) "valori di OmegaA ed OmegaM compatibili con l'eta' dell'universo ottenuti senza vincolo di somma pari a 1"
 WRITE(15,*) "stima eta'  OmegaA  OmegaM  OmegaTOT"
 oa=0.d0
 do i=1,101
   om=0.d0
   do j=1,101
     CALL quadgauss(integranda,x1,x2,inte1)
     CALL quadgauss(integranda2,x1,x2,inte2)
     risultato=978.d0*(inte1+inte2)/h0                     !la moltiplicazione per 978 è dovuta alla conversione delle unità di misura

     if(risultato>=13.69d0.and.risultato<=13.95d0) then    !i valori accettati di discostano dal valore reale di +/- 1%
      WRITE(15,120) risultato,oa,om,oa+om
     end if
     om=om+0.01

   end do
   oa=oa+0.01
 end do
 120 format(3x,f5.2,5x,f5.2,3x,f5.2,3x,f5.2)

CLOSE(15)
WRITE(*,*) '  '
WRITE(*,*) '  '
WRITE(*,*) 'TUTTI I DATI SONO STATI SALVATI NEI FILE output_data.txt omega_val.txt   '

OPEN(103,file="grafdata.txt")
!WRITE(103,223)ngal
!WRITE(103,223)nome
WRITE(103,223)tot_p
WRITE(103,223)tot_erp
WRITE(103,223)tot_mapp
WRITE(103,223)tot_ermap
WRITE(103,223)tot_mv
WRITE(103,223)tot_ermv
WRITE(103,223)tot_dist
WRITE(103,223)tot_erdist
WRITE(103,223)h
WRITE(103,223)errh
223 format(19(','f9.5))
CLOSE(103)
write(*,*) "PREMERE INVIO PER CONCLUDERE L'ESECUZIONE"
READ(*,*)
END PROGRAM cefeidi
!======================================================
!======================================================
SUBROUTINE sort_colonne(vect1, vect2, nd)
  IMPLICIT NONE
  INTEGER:: nd, i, posizione,j
  REAL*8:: valore_rif, provv
  REAL*8:: vect1(nd), vect2(nd)
  REAL*8:: aux1(nd), aux2(nd)
  aux1=vect1
  aux2=vect2
  DO i=1,nd-1
     valore_rif=aux1(i)
     posizione=i
     DO j=i+1,nd
       IF(aux1(j)<valore_rif) THEN
         valore_rif=aux1(j)
         posizione=j
       END IF
     END DO
     aux1(posizione)=aux1(i)
     aux1(i)=valore_rif
     provv=aux2(posizione)
     aux2(posizione)=aux2(i)
     aux2(i)=provv
  END DO
  vect1=aux1
  vect2=aux2
  DO i=1,nd-1
    IF(vect1(i)==vect1(i+1)) THEN
       vect1(i)=vect1(i)+0.0000001  !evita il caso ci siano due valori esattamente uguali nella colonna secondo cui ordino
    END IF
  END DO
END SUBROUTINE  sort_colonne
!------------------------------------------------------
SUBROUTINE mat_gen(x,fdx,n,mat,c)
  IMPLICIT NONE
  INTEGER::i,n
  REAL*8::x(n),fdx(n),mat(n,n),c(n),e(n),g(n),r(n)
  mat=0.d0
  DO i=1,n
    IF(i==1 .or. i==n) THEN
    mat(i,i)=1.d0
    c(i)=0.d0
  ELSE
    e(i)=(x(i)-x(i-1))
    mat(i,i-1)=e(i)
    g(i)=(x(i+1)-x(i))
    mat(i,i+1)=g(i)
    r(i)=(x(i+1)-x(i-1))
    mat(i,i)=2*r(i)
    c(i)=(6.d0/(x(i+1)-x(i)))*(fdx(i+1)-fdx(i))&
       &+(6.d0/(x(i)-x(i-1)))*(fdx(i-1)-fdx(i))
  END IF
END Do
END SUBROUTINE mat_gen
!-----------------------------------------------------
SUBROUTINE gauss(mat,n,c,der)   !der è il vettore che contiene i risultati
  REAL*8::fact,mat(n,n),c(n),der(n),sum
  INTEGER::n,i,j,k
  DO i=1,n-1
    DO j=i+1,n
      fact=mat(j,i)/mat(i,i)
      DO k=1,n
        mat(j,k)=mat(j,k)-mat(j-1,k)*fact
      END DO
      c(j)=c(j)-c(j-1)*fact
    END DO
  END DO
  !trovo i valori delle derivate seconde nei punti noti
  der(n)=c(n)/mat(n,n)
  DO i=n-1,1,-1
    sum=0.d0
    DO j=i+1,n
      sum=sum+mat(i,j)*der(j)
    END DO
    der(i)=(c(i)-sum)/mat(i,i)

  END DO
END SUBROUTINE gauss
!------------------------------------------------------
SUBROUTINE interp(val,x,fdx,der,n,fun)
  IMPLICIT NONE
  REAL*8::val,x(n),fun,f1,f2,f3,f4,f0,der(n),fdx(n)
  INTEGER::i,n
  DO i=2,n
    IF(val>=x(i-1).and.val<=x(i)) THEN
      f0=(x(i)-x(i-1))
      f1=(der(i-1)*(x(i)-val)**3)/(6.d0*f0)
      f2=(der(i)*(val-x(i-1))**3)/(6.d0*f0)
      f3=((fdx(i-1)/f0)-(der(i-1)*f0)/6.d0)*(x(i)-val)
      f4=((fdx(i)/f0)-(der(i)*f0)/6.d0)*(val-x(i-1))
      fun=f1+f2+f3+f4
    END IF
  END DO
END SUBROUTINE interp
!------------------------------------------------------
SUBROUTINE fit(x,y,nd,c)
  IMPLICIT NONE
INTEGER::i,nd
REAL*8:: x(nd),y(nd),c(2),a,b,delta


delta=nd*sum(x**2)-(sum(x))**2
a=(sum(x**2)*sum(y)-sum(x)*sum(x*y))/delta
b=(nd*sum(x*y)-sum(x)*sum(y))/delta
c(1)=a
c(2)=b

END SUBROUTINE fit
!------------------------------------------------------
SUBROUTINE periodi_magnitudini(lc,bin,incr,pmed,errpe,magmed,errmapp)
  IMPLICIT NONE
  !variabili d'uscita
  REAL*8::pmed,errpe,magmed,errmapp,incr
  !variabili interne
  INTEGER::count,test,j,k,l,bin
  REAL*8::massimi(5),minimi(5),erm(6),p(6),m(6),aux1,aux,magrif,lc(bin)
  !inizializzazione
  massimi(:)=0.d0
  minimi(:)=0.d0
  erm(:)=0.d0
  p(:)=0.d0
  m(:)=0.d0
  count=0
  test=0
  !nel ciclo sottostante conto le volte che la curva di luce supera il riferimento
  !durante il fronte d'onda discendente, salvando il corrispondente giorno giuliano
  !del primo e dell'ultimo attraversamento sono in grado di calcolare il periodo medio di ciascuna stella
  magrif=sum(lc)/size(lc)
  DO k=1,bin-1
    IF(lc(k)>magrif.and.lc(k+1)<magrif) THEN
      count=count+1
      massimi(count)=k
    ELSEIF(lc(k)<magrif.and.lc(k+1)>magrif) THEN
      test=test+1
      minimi(test)=k
    END IF
  END DO

  do l=1,count-1                                           !riempio i vettori p ed m con i dati calcolati tramite gli attraversamenti dall'alto
    aux=0.d0
    do k=int(massimi(l)),int(massimi(l+1))
      aux=aux+lc(k)
    end do
    p(l)=(massimi(l+1)-massimi(l))*incr
    m(l)=aux/(massimi(l+1)-massimi(l))
    aux=0.d0
    do k=int(massimi(l)),int(massimi(l+1))
      aux1=aux1+(lc(k)-m(l))**2
    end do
    erm(l)=SQRT(aux1/(massimi(l+1)-massimi(l)))
  end do

  do l=1,test-1                                            !riempio la parte finale di p ed m con i dati calcolati tramite gli attraversamenti dal basso
    aux=0.d0
    do k=int(minimi(l)),int(minimi(l+1))
      aux=aux+lc(k)
    end do
    p(l+count-1)=(minimi(l+1)-minimi(l))*incr
    m(l+count-1)=aux/(minimi(l+1)-minimi(l))
    aux1=0.d0
    do k=int(minimi(l)),int(minimi(l+1))
      aux1=aux1+(lc(k)-m(l))**2
    end do
    erm(l+count-1)=SQRT(aux1/(minimi(l+1)-minimi(l)))
  end do

  pmed=sum(p)/real(count+test-2)
  magmed=sum(m)/real(count+test-2)
  aux=0.d0
  aux1=0.d0
  do j=1,count+test-2
    aux=aux+(p(j)-pmed)**2
    aux1=aux1+(m(j)-magmed)**2
  end do
  errpe=SQRT(aux/real(count+test-3))                       !divido per il numero di valori -1 per avitare di sottostimare l'errore
  errmapp=SQRT(aux1/real(count+test-3))


END SUBROUTINE periodi_magnitudini
!------------------------------------------------------
SUBROUTINE quadgauss(f,a,b,res) !sempre usata con npunti=4
  use dati
  IMPLICIT NONE
  REAL*8, EXTERNAL:: f
  REAL*8:: a,b,res,summa
  REAL*8, ALLOCATABLE:: c(:), x(:), xd(:)
  INTEGER::i

  ALLOCATE(x(0:3),c(0:3), xd(0:3))
  xd(0)=SQRT(3./7.-2./7.*SQRT(6./5.))
  c(0)=(18.+SQRT(30.))/36.
  c(1)=c(0)
  xd(1)=-xd(0)
  c(2)=(18.-SQRT(30.))/36.
  c(3)=c(2)
  xd(2)=SQRT(3./7.+2./7.*SQRT(6./5.))
  xd(3)=-xd(2)
  DO i=0,3
     x(i)=0.5*((b+a)+(b-a)*xd(i))
  END DO
  summa=0.
  DO i=0,3
     summa=summa+c(i)*f(x(i))
  END DO
  res=summa*(b-a)*0.5
END SUBROUTINE quadgauss
!------------------------------------------------------
REAL*8 FUNCTION integranda(z)
  use dati
  IMPLICIT NONE
  REAL*8::z
  integranda=1./((1+z)*SQRT(om*((1+z)**3)+((1-om-oa)*((1+z)**2))+oa))
END FUNCTION
!------------------------------------------------------
REAL*8 FUNCTION integranda2(t)
    use dati
  IMPLICIT NONE
  REAL*8::t
  REAL*8, EXTERNAL:: integranda
  integranda2=integranda(1./t)/(t**2)
END FUNCTION integranda2
!------------------------------------------------------
SUBROUTINE save_txt_data(ngal,nome,vet1,vet2,vet3,vet4,vet5,vet6,vet7,vet8,vet9,vet10,val,erval)
  IMPLICIT NONE
  INTEGER::i,ngal
  REAL*8::vet1(ngal),vet2(ngal),vet3(ngal),vet4(ngal),vet5(ngal),vet6(ngal)
  REAL*8::vet7(ngal),vet8(ngal),vet9(ngal),vet10(ngal),val,erval
  CHARACTER(len=*)::nome(ngal)

  OPEN(14,file='output_data.txt')
  WRITE (14,*)'Nome      periodo    er periodo    mag app      er mag     mag ass     er mag       dist       er dist&
              &       h0         er h0'
  WRITE(14,*)'  '
  DO i=1,ngal
    WRITE(14,130) nome(i),vet1(i),vet2(i),vet3(i),vet4(i),vet5(i),vet6(i),vet7(i),vet8(i),vet9(i),vet10(i)
  end do
WRITE(14,*)'  '
  WRITE(14,*)'  stima finale della costante di Hubble ---> ',val,' +/- ',erval
  130 format(a,10(x,x,x,f9.5))

CLOSE(14)
END SUBROUTINE save_txt_data
