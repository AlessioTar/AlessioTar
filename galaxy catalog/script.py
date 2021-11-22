
import numpy as np
import os
from astropy.io import fits
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.stats as spt


'''definisco tutte le directory che userò nel programma partendo da quella home'''
dir_home=os.getcwd()
if not os.path.exists (dir_home):
    os.mkdir(dir_home)
dir_out=dir_home + '\\output\\'
if not os.path.exists (dir_out):
    os.mkdir(dir_out)
dir_out_4=dir_out + '\\step_4\\'
if not os.path.exists (dir_out_4):
    os.mkdir(dir_out_4)
    dir_data=dir_home+'\\data\\'
if not os.path.exists (dir_data):
    os.mkdir(dir_data)
dir_plots=dir_home + '\\plots\\'
if not os.path.exists (dir_plots):
    os.mkdir(dir_plots)
dir_trend=dir_plots + '\\trend\\'
if not os.path.exists (dir_trend):
    os.mkdir(dir_trend)
dir_hist=dir_plots + '\\histograms\\'
if not os.path.exists (dir_hist):
    os.mkdir(dir_hist)
    
'''zona in cui definisco tutte le funzioni che andrò a utilizzare'''

def ManualMean(array):
    somma=0
    for i in array:
        somma+=i
    mean=somma/len(array)
    return mean

def ManualStd(x):
    somma=0
    for i in x:
        somma+=(i-ManualMean(x))**2
    denominatore=len(x)-1
    std=np.sqrt(somma/denominatore)
    return std
        
    
def gauss(bins,mean,sigma):   
    x=np.zeros(len(bins)-1)
    for i in range(len(x)):
        x[i]=(bins[i]+bins [i+1])/2
    return [x,1/(sigma*np.sqrt(2*np.pi))
            *np.exp(-(x-mean)**2/(2*sigma**2))]

def MAD(val):
    mad=1.4826*np.median(np.absolute(val-(np.median(val))))
    return mad

def check_if_gaussian (x):
    sigma=x.std()
    media=x.mean()
    y=x[np.logical_and(media-sigma<x,x<media+sigma)]
    y2=x[np.logical_and(media-sigma<x,x<media+sigma)]
    if len(y)>0.682*len(x) and len(y2)>0.954*len(x):
        return True
    else:
        return False

def save_fig(nomefig):
    plt.savefig(nomefig,dpi=360)

def values_redshift(grandezza):  #questa funzione plotta delle grandezze in funzione del redshift
    plt.figure(redshift_funct_label[j])
    plt.scatter(redshift,grandezza,c='r',s=3)
    coeff=np.polyfit(redshift,grandezza,1)
    y_fit=np.polyval(coeff,redshift)
    plt.plot(redshift,y_fit, label='best fit grado 1', color='k')
    plt.xlabel('redshift')
    plt.ylabel(redshift_funct_label[j])
    plt.legend(loc='upper right',fontsize='small')

def check_trend(grandezza):
    pearson_coefficient,p_val=spt.pearsonr(redshift,grandezza)
    if 0.7<pearson_coefficient or pearson_coefficient<-0.7:
        return True
    else:
        return False
    



def save_data_step2(array,nome):
    data=[nome,check_if_gaussian(array),'%3.6f'%array.mean(),'%3.6f'%array.std(),'%3.6f'%np.median(array),'%3.6f'%MAD(array)]
    
    return data

def residui_step4(x,y,teoric):
    ax2=plt.subplot(gs[0,0])
    x_inf=x[y<teoric]
    x_sup=x[y>teoric]
    plt.hist(x_inf,bins=30,density=True,alpha=0.4,label='valori x<della rel teorica')
    plt.hist(x_sup,bins=30,density=True,alpha=0.4,color='r',label='valori x>della rel teorica')
    plt.legend(loc='best',fontsize='small')

    ax3=plt.subplot(gs[1,1])
    y_inf=y[y<teoric]
    y_sup=y[y>teoric]
    plt.hist(y_inf,bins=30,density=True,alpha=0.4,orientation='horizontal',label='valori y<della rel teorica')
    plt.hist(y_sup,bins=30,density=True,alpha=0.4,color='r',orientation='horizontal',label='valori y>della rel teorica')
    plt.legend(loc='best',fontsize='small')

    return x_inf,x_sup,y_inf,y_sup


def percentuale(x, tot):
    return (float(len(x))/float(tot))*100




    
###STEP 1
os.chdir(dir_data)
filename='data_SDSS_Info.fit'
hdul=fits.open(filename)
dati=hdul[1]
dati.header
columns=dati.columns

righe=hdul[1].data
righe=righe[righe['ID']==1]             #
righe=righe[righe['lgm_tot_p50']>0]     #faccio lo slicing dei dati outliners
righe=righe[righe['sfr_tot_p50']>-8000] #
righe=righe[righe['oiii_5007_flux_err']>-1000] #


galassia=righe['specobjid']
redshift=righe['z']
magU=righe['petroMag_u']
errmagU=righe['petroMagErr_u']
magG=righe['petroMag_g']
errmagG=righe['petroMagErr_g']
magR=righe['petroMag_r']
errmagR=righe['petroMagErr_r']
magI=righe['petroMag_i']
errmagI=righe['petroMagErr_i']
magZ=righe['petroMag_z']
errmagZ=righe['petroMagErr_z']
aflux=righe['h_alpha_flux']
erraflux=righe['h_alpha_flux_err']
bflux=righe['h_beta_flux']
errbflux=righe['h_beta_flux_err']
o3=righe['oiii_5007_flux']
erro3=righe['oiii_5007_flux_err']
n2=righe['nii_6584_flux']
errn2=righe['nii_6584_flux_err']
mass_50=righe['lgm_tot_p50']
mass_16=righe['lgm_tot_p16']
mass_84=righe['lgm_tot_p84']
sfr_50=righe['sfr_tot_p50']
sfr_16=righe['sfr_tot_p16']
sfr_84=righe['sfr_tot_p84']
absMagU=righe['absMagU']
absMagG=righe['absMagG']
absMagR=righe['absMagR']
absMagI=righe['absMagI']
absMagZ=righe['absMagZ']

print 'ho letto e rinominato le colonne con successo'

###STEP 2
label=[['redshift'],['magU    ','magR    ','magZ    ','magI    ','magG    '],['aflux   ','bflux   ','o3      ','n2      '],['mass_50 ','mass_16 ','mass_84 '],['sfr_50  ','sfr_16  ','sfr_84  '],['absMagU ','absMagR ','absMagZ ','absMagI ','absMagG ']]
nomi=[[redshift],[magU,magR,magZ,magI,magG],[aflux,bflux, o3, n2],[mass_50,mass_16,mass_84],[sfr_50,sfr_16,sfr_84],[absMagU,absMagR,absMagZ,absMagI,absMagG]]
#nomi=[redshift,magU,errmagU,magR,errmagR,magZ,errmagZ,magI,errmagI,magG,errmagG,aflux,erraflux,bflux,errbflux,o3,erro3,n2,errn2,mass_50,mass_16,mass_84,sfr_50,sfr_16,sfr_84,absMagU,absMagR,absMagZ,absMagI,absMagG]
figure=['redshift','magnitude','flux','mass','star_formation_ratio','absolute_magnitude']
errori=[[errmagU,errmagR,errmagZ,errmagI,errmagG],[erraflux,errbflux],erro3,errn2]



#plotto e salvo le mie figure che mostrano la gaussiana sovrapposta all'istogramma nella prima riga della figura e i residui nella seconda


for j in range(len(nomi)):
    f=plt.figure(figure[j],figsize=(18,10))
    for i in range(len(nomi[j])):
        plt.figure()
        count,bins,ign=plt.hist(nomi[j][i],130,density=True)
        plt.close()
        
        x_gaus,y_gaus=gauss(bins,nomi[j][i].mean(),nomi[j][i].std())
        
        gs=gridspec.GridSpec(2,len(nomi[j]))
        
        ax1=plt.subplot(gs[0,i])
        
        plt.hist(nomi[j][i],130,histtype='bar',density=True, label=label[j][i])
        plt.plot(x_gaus,y_gaus,ls='--',c='k',label='gaussian_fit')
        
#se la distrubuzione dei dati è gaussiana sul grafico apparirà una retta verticale denominata mean, se non è gaussiana apparirà median
        
        if check_if_gaussian(nomi[j][i])==True: 
            plt.axvline(nomi[j][i].mean(), label='mean')
            plt.axvspan(nomi[j][i].mean()-nomi[j][i].std(),nomi[j][i].mean()+nomi[j][i].std(),alpha=0.3,color='k')
            
        else:
            plt.axvline(np.median(nomi[j][i]), label='median')
            plt.axvspan(np.median(nomi[j][i])-MAD(nomi[j][i]),np.median(nomi[j][i])+MAD(nomi[j][i]),alpha=0.3,color='k')

        plt.legend(loc='best',fontsize='small')
        
        ax2=plt.subplot(gs[1,i])
        
        residui=count-y_gaus
        plt.scatter(x_gaus,residui,s=5,label='residui')
        plt.axhline(0,c='k')
        plt.axhspan(-np.std(residui),np.std(residui),alpha=0.3)
        plt.legend(loc='best',fontsize='small')

    os.chdir(dir_hist)
    save_fig(figure[j])
    plt.close(f)

print 'ho salvato le immagini .png di tutti gli istogrammi nella cartella "histograms" con successo'

#salvo in un file .txt i dati statistici dei miei dati

os.chdir(dir_out)
t=0
sav_dati=np.zeros(22,dtype=object)
for j in range(len(nomi)):        
    for i in range(len(nomi[j])):
        sav_dati[t]=save_data_step2(nomi[j][i],label[j][i])
        t=t+1

sav_dati[21]=['absMagG ',check_if_gaussian(absMagG),ManualMean(absMagG),ManualStd(absMagG), 'dati calcolati con le funzioni definite manualmente']

'''dato che per un caso devo stimare la mean e la std con delle funzioni scritte da me aggiungo un altra riga al mio
   array sav_dati in cui sono contenuti questi due valori. per far si che siano più facilmente confrontabili ho scelto
   di fare questa operazione su absMagG
'''

fileout_data='valori_step_2.txt'
fil=open(fileout_data, 'w')
np.savetxt(fil,sav_dati,header='nome   is gaussian    mean        std          median       mad    ',fmt='%s', delimiter='  ,      .  ',newline='\n',comments='#')
fil.close()

print 'ho salvato nella cartella "output" i dati statistici dello step 2'
###STEP 3
redshift_funct_label=['magG','aflux','mass_50','sfr_50','absMagG']
redshift_funct=[[magG],[aflux],[mass_50],[sfr_50],[absMagG]]

for j in range(len(redshift_funct)):        
    for i in range(len(redshift_funct[j])):
        values_redshift(redshift_funct[j][i])     #plotta tutti i grafici scatter con anche i fit di grado 1
        #check_trend(redshift_funct[j][i])         #verifica se ci sono dei trend
        filename=redshift_funct_label[j]+'__trend'
        os.chdir(dir_trend)
        save_fig(filename)

                
        plt.close()
os.chdir(dir_data)
print 'ho plottato le grandezze in funzione della redshift e verificato eventuali trend'
'''siccome l'unico trend che trovo è nella massa studio le proprietà solo al 50esimo percentile'''

for j in range(len(redshift_funct)):
    if check_trend(redshift_funct[j][0])==True:
        print redshift_funct_label[j]+' presenta un trend con il redshift'
        

        valore_1=redshift_funct[j][0][np.logical_and(0<redshift,redshift<=0.025)]
        valore_2=redshift_funct[j][0][np.logical_and(0.025<redshift,redshift<=0.05)]
        valore_3=redshift_funct[j][0][np.logical_and(0.05<redshift,redshift<0.075)]
        valore_4=redshift_funct[j][0][np.logical_and(0.075<redshift,redshift<0.1)]
        lista_valori=[valore_1,valore_2,valore_3,valore_4]
        print 'dividendolo noto che: '
        print 'il primo blocco contiene il '+str(  len(valore_1)*100.0/len(redshift_funct[j][0])) +'% dei valori '+'la cui media, std e mediana valgono rispettivamente: '+str(np.mean(valore_1))+'    '+str(np.std(valore_1))+'    '+str(np.median(valore_1))
        print 'il secondo blocco contiene il '+str(len(valore_2)*100.0/len(redshift_funct[j][0])) +'% dei valori '+'la cui media, std e mediana valgono rispettivamente: '+str(np.mean(valore_2))+'    '+str(np.std(valore_2))+'    '+str(np.median(valore_2))
        print 'il terzo blocco contiene il '+ str( len(valore_3)*100.0/len(redshift_funct[j][0])) +'% dei valori '+'la cui media, std e mediana valgono rispettivamente: '+str(np.mean(valore_3))+'    '+str(np.std(valore_3))+'    '+str(np.median(valore_3))
        print 'il quarto blocco contiene il '+str( len(valore_4)*100.0/len(redshift_funct[j][0])) +'% dei valori '+'la cui media, std e mediana valgono rispettivamente: '+str(np.mean(valore_4))+'    '+str(np.std(valore_4))+'    '+str(np.median(valore_4))
        for i in range(4):

            gs1=gridspec.GridSpec(2,4)
            plt.figure()
            count,bins,ign=plt.hist(lista_valori[i],bins=60,density=True)
            plt.close()

            plt.figure(200,figsize=(18,10))
            ax1=plt.subplot(gs[0,i])
            plt.hist(lista_valori[i],bins=40,density=True, label=redshift_funct_label[j])
            x_gauss,y_gauss=gauss(bins,lista_valori[i].mean(),lista_valori[i].std())
            plt.plot(x_gauss,y_gauss, color='k', label='gaussian_fit', ls='--')
            if check_if_gaussian(lista_valori[i])==True:
                plt.axvline(lista_valori[i].mean, label='mean', color='k', zorder=1)
            else:
                plt.axvline(np.median(lista_valori[i]), label='median', color='k', zorder=1)

            plt.legend(loc='best',fontsize='small')

            ax2=plt.subplot(gs[1,i])
            residui=count-y_gauss
            plt.scatter(x_gauss,residui)
            plt.axhline(0, color='k')
            plt.axhspan(-np.std(residui),np.std(residui), alpha=0.3, color='r')
            
        os.chdir(dir_plots)
        save_fig('divided trend redshift')
        plt.close(200)

###STEP 4

#BPT
redshift_bpt=redshift[np.logical_and(np.logical_and(bflux>0,o3>0),np.logical_and(aflux>0,n2>0))]

o3_bpt=o3[np.logical_and(np.logical_and(bflux>0,o3>0),np.logical_and(aflux>0,n2>0))]
bflux_bpt=bflux[np.logical_and(np.logical_and(bflux>0,o3>0),np.logical_and(aflux>0,n2>0))]
n2_bpt=n2[np.logical_and(np.logical_and(bflux>0,o3>0),np.logical_and(aflux>0,n2>0))]
aflux_bpt=aflux[np.logical_and(np.logical_and(bflux>0,o3>0),np.logical_and(aflux>0,n2>0))]
lins=np.linspace(-4,-0.05,835)

kauf=(0.61/(lins-0.05))+1.3

log_n2=np.log(n2_bpt/aflux_bpt)
log_o3=np.log(o3_bpt/bflux_bpt)
              
plt.figure('bpt',figsize=(18,10))
gs=gridspec.GridSpec(2,2)
ax1=plt.subplot(gs[1,0])
plt.scatter(log_n2,log_o3, c=redshift_bpt,cmap='viridis', s=3,label='bpt')
plt.plot(lins, kauf, ls='--', color='k',label='relazione teorica')
plt.xlim(-4,2)
plt.ylim(-4,4)
plt.xlabel('log_n2')
plt.ylabel('log_o3')
plt.legend(loc='best',fontsize='small')
p,o,i,u=residui_step4(log_n2,log_o3,kauf)
os.chdir(dir_plots)
save_fig('bpt_diagram')
dati_bpt=[['bpt_inf',check_if_gaussian(p), percentuale(p,len(o3_bpt)),np.mean(p),np.median(p),np.std(p),MAD(p)],['bpt_sup',check_if_gaussian(o), percentuale(o,len(o3_bpt)),np.mean(o),np.median(o),np.std(o),MAD(o)]]


print 'ho plottato e salvato il grafico BPT'




#color-mass
gs=gridspec.GridSpec(2,2)
plt.figure('color-mass')
teor_funct=-0.495+(0.25*mass_50)
u_r=magU-magR
ax1=plt.subplot(gs[1,0])
plt.xlabel('mass_50')
plt.ylabel('U-R')
plt.scatter(mass_50,u_r,c=mass_50,cmap='YlOrRd',s=5)
plt.plot(mass_50,teor_funct, ls='--',color='k')
a,s,d,f=residui_step4(mass_50,u_r,teor_funct)
save_fig('color-mass')
dati_color=[['color_inf',check_if_gaussian(a), percentuale(a,len(mass_50)),np.mean(a),np.median(a),np.std(a),MAD(a)],['color_sup',check_if_gaussian(s), percentuale(s,len(mass_50)),np.mean(s),np.median(s),np.std(s),MAD(s)]]
#plt.show()
print 'ho plottato e salvato il grafico color-mass'

#sfr-mass
sfr_teor=(0.76*mass_50)-8.64

gs=gridspec.GridSpec(2,2)
plt.figure('sfr-mass')
ax1=plt.subplot(gs[1,0])
plt.xlabel('mass_50')
plt.ylabel('sfr_50')
plt.plot(mass_50,sfr_teor,label='relazione teorica',c='k')
plt.scatter(mass_50,sfr_50,c=mass_50,cmap='RdYlGn_r',s=5)

q,w,e,r=residui_step4(mass_50,sfr_50,sfr_teor)
save_fig('sfr-mass')
print 'ho plottato e salvato il grafico sfr-mass'
plt.close()
dati_sfr=[['sfr_inf',check_if_gaussian(q), percentuale(q,len(sfr_50)),np.mean(q),np.median(q),np.std(q),MAD(q)],['sfr_sup',check_if_gaussian(w), percentuale(w,len(sfr_50)),np.mean(w),np.median(w),np.std(w),MAD(w)]]


os.chdir(dir_out_4)
fil_bpt=open('dati_bpt.txt', 'w')
np.savetxt(fil_bpt,dati_bpt,    header='nome   gaussiano   percentuale     media        mediana      std           mad',fmt='%s', delimiter=' , ',newline='\n',comments='#')
fil_bpt.close()

fil_color=open('dati_color.txt', 'w')
np.savetxt(fil_color,dati_color,header='nome   gaussiano   percentuale     media        mediana      std           mad' ,fmt='%s', delimiter=' , ',newline='\n',comments='#')
fil_color.close()

fil_sfr=open('dati_sfr.txt', 'w')
np.savetxt(fil_sfr,dati_sfr,    header='nome   gaussiano   percentuale     media        mediana      std           mad' ,fmt='%s', delimiter=' , ' ,newline='\n',comments='#')
fil_sfr.close()
















