PROGRAM istogramma
  INTEGER::i,r,n_bin,bin,opt
  INTEGER, ALLOCATABLE:: contabin(:)
  REAL*8, ALLOCATABLE::colonna1(:),colonna2(:),col1_copy(:),col2_copy(:),bin_lim(:)
  REAL*8::massimo,minimo,bin_dim

r=0
OPEN(21, file='m_t.dat')
DO
  READ(21,*, END=10)
  r=r+1
END DO
10 WRITE(*,*) 'HO LETTO', r, 'RIGHE'
ALLOCATE(colonna1(r),colonna2(r))
REWIND(21)
DO i=1,r
  READ(21,*) colonna1(i),colonna2(i)
END DO
WRITE(*,*) "seleziona colonna 1 o colonna 2 du cui fare l'istogramma"
READ(*,*) opt
WRITE(*,*) "seleziona seleziona numero di bin dell'istogramma"
READ(*,*) n_bin
ALLOCATE(contabin(nbin),bin_lim(nbin+1))
contabin=0
bin_lim=0

SELECT CASE(opt)
CASE (1)
  massimo=MAXVAL(colonna1)*1.0000001
  WRITE(*,*)'valore massimo', massimo
  minimo=MINVAL(colonna1)*0.9999999
  WRITE(*,*) 'valore minimo', minimo
  bin_dim=(massimo-minimo)/REAL(n_bin)
  WRITE(*,*)'dimensione bin', bin_dim
  WRITE(*,*)


  DO i=1,nbin+1
    IF(i==1)THEN
      bin_lim(i)=REAL(minimo)
      cycle
    END IF
    bin_lim(i)=bin_lim(i-1)+bin_dim
  END DO
WRITE(*,*) 'LIMITI DEI BIN', (bin_lim(i),i=1,n_bin+1)
WRITE(*,*)'  DATO         CONTEGGIO DEI RISPETTIVI DATI APPARTENENTI AI BIN '
  DO i=1,r
    bin=((colonna1(i)-minimo)/bin_dim)+1
    contabin(bin)=contabin(bin)+1
    WRITE(*,81)colonna1(i)
    WRITE(*,*)'      ', (contabin(j), j=1,n_bin)
  END DO


CASE (2)
  massimo=MAXVAL(colonna2)*1.0000001
  WRITE(*,*)'valore massimo', massimo
  minimo=MINVAL(colonna2)*0.9999999
  WRITE(*,*) 'valore minimo', minimo
  bin_dim=(massimo-minimo)/REAL(n_bin)
  WRITE(*,*)'dimensione bin', bin_dim
  WRITE(*,*)
  DO i=1,nbin+1
    IF(i==1)THEN
      bin_lim(i)=REAL(minimo)
      cycle
    END IF
    bin_lim(i)=bin_lim(i-1)+bin_dim
  END DO

  WRITE(*,*) 'LIMITI DEI BIN', (bin_lim(i),i=1,n_bin+1)
  WRITE(*,*)'  DATO         CONTEGGIO DEI RISPETTIVI DATI APPARTENENTI AI BIN '

  DO i=1,r
    bin=((colonna2(i)-minimo)/bin_dim)+1
    contabin(bin)=contabin(bin)+1
    WRITE(*,81)colonna2(i)
    WRITE(*,*) '      ', (contabin(j), j=1,n_bin)
  END DO


END SELECT

81 format(f9.4,x)
END PROGRAM istogramma
