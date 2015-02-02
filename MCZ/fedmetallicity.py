##############################################################################
##Calculate metalicity, code originally in IDL written by Lisa
##new calculation based on the most recent version of the .pro file.
##
##I ignored the section she had of not calculating the lines that didn't
##meet some mass condition.
##inputs:
## measured - flux data, must be the format returned by readfile()
## num - also returned by readfile()
## outfilename - the name of the file the results will be appended to
## red_corr - True by default
## disp - if True prints the results, default False
## saveres - if True appends the results onto outfilename, True by default
##############################################################################
import numpy as np
from pylab import hist,show

##all the lines that go like
##if S_mass[i] > low_lim and S_mass[i] < 14.0 and all_lines[i] != 0.0:
##are ignored, replaced by (if not G) where G is set to False

G=True#False

##list of metallicity methods, in order calculated
Zs=["KD02comb_updated", #always
    "KD02_N2O2", #Halpha, Hbeta,  [OII]3727, [NII]6584
    "KD03_N2Ha",#Halpha, Hbeta,  [OII]3727, [NII]6584
    "KD03new_R23", #Hbeta,  [OII]3727, [OIII]5007, [OIII]4959 
    "M91", #Hbeta,  [OII]3727,  [OIII]5007, [OIII]4959 
    "Z94", #Hbeta,  [OII]3727, [OIII]5007, [OIII]4959 
    "PP04_N2",   #Halpha, [NII]6584
    "PP04_O3N2", #Halpha, Hbeta,  [OIII]5007, [NII]6584
    "Pi01",  #Hbeta,  [OII]3727,  [OIII]5007
    "D02", #Halpha, [NII]
    "E(B-V)", #Halpha, Hbeta]
    "C01_R23","C01"]

def get_keys():
    return Zs

##############################################################################
##fz_roots function as used in the IDL code  FED:reference the code here!
##############################################################################

def calculation(diags,measured,num,(bsmeas,bserr),outfilename='blah.txt',red_corr=True,disp=False,saveres=False): 

    raw_lines={}

    for k in measured.iterkeys():
        #kills all non-finite terms 
        measured[k][np.where(np.isfinite(measured[k][:])==False)]=0.0 
        raw_lines[k]=measured[k]
    
    
    raw_lines['[OIII]4959']=raw_lines['[OIII]5007']/3.

    diags.setHab(raw_lines['Ha'],raw_lines['Hb'])
    

    #if Ha or Hb is zero, cannot do red correction
    if red_corr and diags.hasHa and diags.hasHb:
        with np.errstate(invalid='ignore'):
            #print 'extinction correction ',i
            diags.calcEB_V()
    elif red_corr:
        response=raw_input("WARNING: reddening correction cannot be done without both H_alpha and H_beta measurement!! Continuing without reddening correction? [Y/n]\n").lower()
        assert(not (response.startswith('n'))),"please fix the input file to include Halpha and Hbeta measurements"
        red_corr=False
        diags.EB_V=np.ones(len(raw_lines['Ha']))*1e-5
    else:
        diags.EB_V=np.ones(len(raw_lines['Ha']))*1e-5
        

    #####setting up lines first
    diags.setOlines(raw_lines['[OII]3727'], raw_lines['[OIII]5007'], raw_lines['[OI]6300'], raw_lines['[OIII]4959'])
    diags.setNII(raw_lines['[NII]6584'])
    diags.setS(raw_lines['[SII]6717'],raw_lines['[SII]6731'],raw_lines['[SIII]9069'],raw_lines['[SIII]9532']),

    diags.calcNIIOII()
    diags.calcNIISII()
    diags.calcSIIHa()

    diags.calcR23()
    #diags.calcS23()

    diags.calcZ94_Z()
    diags.calcD02_Z()
    diags.calcPP04_N2_Z()

    diags.Pmethod()

    diags.calcPP04_O3N2_Z()
    diags.calcKD02_N2O2_Z()
    diags.calcKD03_NHa()
    #diags.calcC01_ZR23()
    diags.calcM91()

'''

            ''
            else:
                if (raw_lines['[NII]6584'][i] != 0.0) and (raw_lines['[OII]3727'][i] != 0.0) :
                    logNIIOII[i]=np.log10(raw_lines['[NII]6584'][i]/raw_lines['[OII]3727'][i])

                if (raw_lines['[OIII]49595007'][i] != 0.0) and (raw_lines['[OII]3727'][i] != 0.0) :    
                    logR23[i]=np.log10((raw_lines['[OII]3727'][i] + raw_lines['[OIII]49595007'][i])/raw_lines['Hb'][i])
                    R23[i]=((raw_lines['[OII]3727'][i] + raw_lines['[OIII]49595007'][i])/raw_lines['Hb'][i])
                    logOIIIOII[i]=np.log10(raw_lines['[OIII]49595007'][i]/raw_lines['[OII]3727'][i])     

                if (SII67176731_raw[i] != 0.0) and (SIII90699532_raw[i] != 0.0) :    
                    logS23[i]=np.log10((SII67176731_raw[i] + SIII90699532_raw[i])/raw_lines['Hb'][i])
                    logSIIISII[i]=np.log10(SIII90699532_raw[i]/SII67176731_raw[i])
                  
                if (raw_lines['[OIII]5007'][i] != 0.0) and (raw_lines['[OII]3727'][i] != 0.0) :
                    logOIII5007OII[i]=np.log10(raw_lines['[OIII]5007'][i]/raw_lines['[OII]3727'][i])
                    logOIIOIII5007[i]=np.log10(raw_lines['[OII]3727'][i]/raw_lines['[OIII]5007'][i])
                    logOIIOIII5007Hb[i]=np.log10((raw_lines['[OII]3727'][i] + raw_lines['[OIII]5007'][i])/raw_lines['Hb'][i])
                    OIII5007OII[i]=10**(logOIII5007OII[i])
                    OIIOIII5007[i]=10**(logOIIOIII5007[i])
                    # ratios for other diagnostics - slightly different ratios needed
                    #         OIIIOII[i]=1.347*OIII5007OII[i]
                    #         logOIIIOII[i]=np.log10(OIIIOII[i])
                    # R23 but without [OIII]4959
                    R23_5007[i]=(1./OIII5007OII[i] + 1.)/(1./OIII5007OII[i] + 1.347)*R23[i] 
                    OIIOIII5007Hb[i]=R23_5007[i]
                    logOIIOIII5007Hb[i]=np.log10(OIIOIII5007Hb[i])

                if (raw_lines['[NII]6584'][i] != 0.0) and (SII67176731_raw[i] != 0.0) :  
                    logNIISII[i]=np.log10(raw_lines['[NII]6584'][i]/SII67176731_raw[i])

                if (raw_lines['[NII]6584'][i] != 0.0) and (raw_lines['Ha'][i] != 0.0) :  
                    logNIIHa[i]=np.log10(raw_lines['[NII]6584'][i]/raw_lines['Ha'][i])

                if (SII67176731_raw[i] != 0.0) and (raw_lines['Ha'][i] != 0.0) :  
                    logSIIHa[i]=np.log10(SII67176731_raw[i]/raw_lines['Ha'][i])

                if (raw_lines['[OI]6300'][i] != 0.0) and (raw_lines['Ha'][i] != 0.0) :  
                    logOIHa[i]=np.log10(raw_lines['[OI]6300'][i]/raw_lines['Ha'][i])

                if (raw_lines['[OIII]5007'][i] != 0.0) and (raw_lines['Hb'][i] != 0.0) :
                    logOIIIHb[i]=np.log10(raw_lines['[OIII]5007'][i]/raw_lines['Hb'][i])

                if (raw_lines['[OIII]5007'][i] != 0.0) and (raw_lines['[OIII]4959'][i] != 0.0) and (raw_lines['Hb'][i] != 0.0) :
                    logOIII49595007Hb[i]=np.log10((raw_lines['[OIII]5007'][i]+raw_lines['[OIII]4959'][i])/raw_lines['Hb'][i])

                if (raw_lines['[OII]3727'][i] != 0.0) and (raw_lines['Hb'][i] != 0.0) :
                    logOIIHb[i]=np.log10(raw_lines['[OII]3727'][i]/raw_lines['Hb'][i])
         
            ''

    N2=logNIIHa
    D02_Z=9.12+0.73*N2

    #### Pettini & Pagel diagnostics - Pettini & Pagel 2004, MNRAS, 348, L59

    #print 'PP04'
    PP04_N2_Z=9.37 + 2.03* N2 + 1.26* N2**2 + 0.32* N2**3

    O3N2=logOIIIHb-logNIIHa

    PP04_O3N2_Z=8.73 - 0.32*O3N2

    #print 'setting initial guess'

    # NEW Initial Guess - Nov 1 2006
    Z_init_guess=np.zeros(num)+8.7 # upper branch if nothing else
    Z_branch=np.zeros(num)
    # [NII]/Ha if no [NII]/[OII]

    Z_init_guess[(logNIIHa < -1.3)&(raw_lines['[NII]6584'] != 0.0)&(raw_lines['Ha'] != 0.0)]=8.2
    Z_init_guess[(logNIIHa < -1.1)&(logNIIHa >= -1.3)&(raw_lines['[NII]6584'] != 0.0)&(raw_lines['Ha'] != 0.0)]=8.4
    Z_init_guess[(logNIIHa >= -1.1) & (raw_lines['[NII]6584'] != 0.0) & (raw_lines['Ha'] != 0.0)]=8.7


    NIIOII_lines=raw_lines['[NII]6584']*raw_lines['Ha']*raw_lines['Hb']*raw_lines['[OII]3727']

    # [NII]/[OII] if at all possible

    Z_init_guess[(logNIIOII < -1.2)&(NIIOII_lines != 0.0)]=8.2  #  1.2 using low-Z gals,
    Z_init_guess[(logNIIOII >= -1.2)&(NIIOII_lines != 0.0)]=8.7  # 1.1 using HIi regions


    # #### P-method #####
    #

    #print 'P01'

    P_abund_up=np.zeros(num)
    P_abund_low=np.zeros(num)
    #P = 10**logOIIIOII/(1+10**logOIIIOII)
    R3=10**logOIII49595007Hb
    R2=10**logOIIHb
    P = R3/(R2+R3)
    P_R23=R2+R3

    for j in range(num):
        if G :
            #print 'P01 ',j
            P_abund_up[j]=(P_R23[j]+726.1+842.2*P[j]+337.5*P[j]**2)/(85.96+82.76*P[j]+43.98*P[j]**2+1.793*P_R23[j])
            P_abund_low[j]=(P_R23[j]+106.4+106.8*P[j]-3.40*P[j]**2)/(17.72+6.60*P[j]+6.95*P[j]**2-0.302*P_R23[j])


    Pi01_Z=P_abund_up.copy()
    Pi01_Z[(Z_init_guess < 8.4)]=P_abund_low[(Z_init_guess < 8.4)].copy()
    Pi01_Z[(Z_init_guess >= 8.4)]=P_abund_up[(Z_init_guess >= 8.4)].copy()

    # P-method 2001 upper branch

    #print 'P01 old'
    P_abund_old=np.zeros(num)
    Pi01_Z_old=np.zeros(num)
    P = 10**logOIIIOII/(1+10**logOIIIOII)#
    P_R23=10**logR23
    for j in range(num) :
        #print 'P01 old ',j
        if  G :
            P_abund_old[j]=(P_R23[j]+54.2+59.45*P[j]+7.31*P[j]**2)/(6.07+6.71*P[j]+0.371*P[j]**2+0.243*P_R23[j])

    Pi01_Z_old[(Z_init_guess >= 8.4)]=P_abund_old[(Z_init_guess >= 8.4)].copy()

    # P-method lower branch




    # ## calculating Kewley & Dopita (2002) (KD02) estimates of abundance ##

    # KD02 [NII]/[OII] estimate (can be used for whole log(O/H)+12 range, but rms 
    # scatter increases to 0.11 rms for log(O/H)+12 < 8.6#  rms = 0.04 for
    # log(O/H)+12 > 8.6
    # uses equation (4) from paper:

    NIIOII_roots=np.zeros((4,num),dtype=complex)
    KD02_NIIOII_Z=np.zeros(num)
    NIIOII_coef=np.array([1106.8660,-532.15451,96.373260,-7.8106123,0.23928247])# q=2e7 line (approx average)

    for i in range(num):
        if G :
            #print '[NII]/[OII]  ',i
     
            if (raw_lines['[NII]6584'][i] != 0.0) and (raw_lines['[OII]3727'][i] != 0.0) :
                if (raw_lines['Ha'][i] != 0.0) and (raw_lines['Hb'][i] != 0.0) :
                    
                    NIIOII_coef[0]=1106.8860-logNIIOII[i]
                    NIIOII_roots[:,i]=fz_roots(NIIOII_coef)          # finding roots for == (4)
                    for k in range(4) :                                         # root must be real and
                        if (NIIOII_roots[k,i].imag) == 0.0 :        # between 7.5 and 9.4
                            if (abs(NIIOII_roots[k,i]) >= 7.5) and (abs(NIIOII_roots[k,i]) <= 9.4) :
                                KD02_NIIOII_Z[i]=abs(NIIOII_roots[k,i])

         

    B_NIIOII=logNIIOII.copy()


    # calculating [NII]/Ha abundance estimates using [OIII]/[OII] also


    nite=5  # number of iteations+1
    Z_new_NIIHa=np.zeros((nite,num))
    Z_new_NIIHa[0,:]=KD02_NIIOII_Z.copy()  # was 8.6
    logq=np.zeros((nite,num))
    Zmax=np.zeros((nite,num))
    noOIIIOII=0.
    logq_NIIOII=np.zeros(num)

    for j in range(num) :
        if G :
            #print '[NII]/Ha',j
            for i in range(nite-1):
       # ionization parameter
                if (raw_lines['[NII]6584'][j] != 0.0) :
                    if (raw_lines['Ha'][j] > 0.0) and (raw_lines['Hb'][j] > 0.0) :
                        if (raw_lines['[OIII]5007'][j] != 0.0) and (raw_lines['[OII]3727'][j] != 0.0) :
                            logq[i,j]=(32.81 + 0.0*logOIIIOII[j]-1.153*logOIIIOII[j]**2 + Z_new_NIIHa[i,j]*(-3.396 -0.025*logOIIIOII[j] + 0.1444*logOIIIOII[j]**2))/(4.603-0.3119*logOIIIOII[j] -0.163*logOIIIOII[j]**2+Z_new_NIIHa[i,j]*(-0.48 + 0.0271*logOIIIOII[j]+ 0.02037*logOIIIOII[j]**2))

    # calculating logq using the [NII]/[OII] metallicities for comparison

                            logq_NIIOII[j]=(32.81 + 0.0*logOIIIOII[j]-1.153*logOIIIOII[j]**2 + KD02_NIIOII_Z[j]*(-3.396 -0.025*logOIIIOII[j] + 0.1444*logOIIIOII[j]**2))/(4.603-0.3119*logOIIIOII[j] -0.163*logOIIIOII[j]**2+KD02_NIIOII_Z[j]*(-0.48 + 0.0271*logOIIIOII[j]+ 0.02037*logOIIIOII[j]**2))

                        else :
                            logq[i,j]=7.37177
     #       no_OIIIOII_yes_HaHb_flag(noOIIIOII)=j
                            noOIIIOII+=1
          
                    else:
                        logq[i,j]=7.37177
                        logNIIHa[j]=np.log10(raw_lines['[NII]6584'][j]/raw_lines['Ha'][j])

                    Z_new_NIIHa[i+1,j]=(7.04 + 5.28*logNIIHa[j]+6.28*logNIIHa[j]**2+2.37*logNIIHa[j]**3)-logq[i,j]*(-2.44-2.01*logNIIHa[j]-0.325*logNIIHa[j]**2+0.128*logNIIHa[j]**3)+10**(logNIIHa[j]-0.2)*logq[i,j]*(-3.16+4.65*logNIIHa[j])
       


    NIIHalogq=logq.copy()
    KD03_NIIHa_abund=Z_new_NIIHa[nite-1,:].copy()


    # calculating upper and lower metallicities for objects without
    # Hb  and for objects without [OIII] and/or [OII]

    Z_new_NIIHa_up=np.zeros(num)
    logq_up=np.zeros(num)
    Zmax_up=np.zeros(num)
    Z_new_NIIHa_low=np.zeros(num)
    logq_low=np.zeros(num)
    Zmax_low=np.zeros(num)
    KD03_NIIHa_abund_low=np.zeros(num)
    KD03_NIIHa_abund_up=np.zeros(num)
    Hb_up_flag=np.zeros(num)
    flagnoHb=0
    Hb_up_ID=np.zeros(100)

    for i in range(num) :
        if G :
            #print '[NII]/Ha, q',i
            if (raw_lines['[NII]6584'][i] != 0.0) :
                if (raw_lines['Ha'][i] > 0.0) and (raw_lines['Hb'][i] == 0.0) :
                    logNIIHa[i]=np.log10(raw_lines['[NII]6584'][i]/raw_lines['Ha'][i])
    # ionization parameter
                    logq_low[i]=6.9
                    logq_up[i]=8.38
           
                    Z_new_NIIHa_low[i]=(7.04 + 5.28*logNIIHa[i]+6.28*logNIIHa[i]**2+2.37*logNIIHa[i]**3)-logq_low[i]*(-2.44-2.01*logNIIHa[i]-0.325*logNIIHa[i]**2+0.128*logNIIHa[i]**3)+10**(logNIIHa[i]-0.2)*logq_low[i]*(-3.16+4.65*logNIIHa[i])

                    Z_new_NIIHa_up[i]=(7.04 + 5.28*logNIIHa[i]+6.28*logNIIHa[i]**2+2.37*logNIIHa[i]**3)-logq_up[i]*(-2.44-2.01*logNIIHa[i]-0.325*logNIIHa[i]**2+0.128*logNIIHa[i]**3)+10**(logNIIHa[i]-0.2)*logq_up[i]*(-3.16+4.65*logNIIHa[i])
                    KD03_NIIHa_abund_low[i]= Z_new_NIIHa_low[i].copy()
                    KD03_NIIHa_abund_up[i]= Z_new_NIIHa_up[i].copy()

                    Hb_up_flag[i]=1
                    Hb_up_ID[flagnoHb]=i
                    flagnoHb+=1
        
        
            if (raw_lines['[NII]6584'][i] != 0.0) :
                if (raw_lines['Ha'][i] > 0.0) and (raw_lines['Hb'][i] > 0.0) :
                    if (raw_lines['[OIII]5007'][i] == 0.0) or (raw_lines['[OII]3727'][i] == 0.0) :
                        logq_low[i]=6.9
                        logq_up[i]=8.38

                        Z_new_NIIHa_low[i]=(7.04 + 5.28*logNIIHa[i]+6.28*logNIIHa[i]**2+2.37*logNIIHa[i]**3)-logq_low[i]*(-2.44-2.01*logNIIHa[i]-0.325*logNIIHa[i]**2+0.128*logNIIHa[i]**3)+10**(logNIIHa[i]-0.2)*logq_low[i]*(-3.16+4.65*logNIIHa[i])

                        Z_new_NIIHa_up[i]=(7.04 + 5.28*logNIIHa[i]+6.28*logNIIHa[i]**2+2.37*logNIIHa[i]**3)-logq_up[i]*(-2.44-2.01*logNIIHa[i]-0.325*logNIIHa[i]**2+0.128*logNIIHa[i]**3)+10**(logNIIHa[i]-0.2)*logq_up[i]*(-3.16+4.65*logNIIHa[i])
                        KD03_NIIHa_abund_low[i]= Z_new_NIIHa_low[i].copy()
                        KD03_NIIHa_abund_up[i]= Z_new_NIIHa_up[i].copy()

                      




    B_R23=logR23.copy()

    # ## calculating z from Kobulnicky parameterization of Zaritzky et al. (1994) ##




    # ## NEW (May 5 2004)  Initial guess using (when available) 
    # 1. [NII]/[OII] 
    # 2. [NII]/Ha 
    # 3. assume upper branch

    # old initial guess (pre-Nov 1 2006)

    #Z_init_guess=np.zeros(num)
    #Z_init_guess=KD03_NIIHa_abund 
    #for j in range(num) :
    #    if S_mass[j] < low_lim then goto, endloop5
    #     if S_mass[j] > 14.0 then goto, endloop5
    #     if all_lines[j] == 0.0 then goto, endloop5
    #     print 'R23 initial guess',j
    #      if (raw_lines['[NII]6584'][j] != 0.0) and (raw_lines['[OII]3727'][j] != 0.0) :
    #      if (raw_lines['Ha'][j] != 0.0) and (raw_lines['Hb'][j] != 0.0) :
    #      if (KD02_NIIOII_Z[j] != 0.0) :
    #           Z_init_guess[j]=KD02_NIIOII_Z[j]
    #      
    #      
    #           
    #   if (raw_lines['[NII]6584'][j] == 0.0) or (raw_lines['Ha'][j] == 0.0) :
    #          Z_init_guess[j]=8.7
    #      
    #endloop5:
    #end


    # ## calculating Charlot 01 R23 calibration: (case F) ##

    x2=OIIOIII5007/1.5
    x3=(10**logOIIIHb)/2.
    C01_ZR23=np.zeros(num)

    for i in range(num) : 
        if OIIOIII5007[i] < 0.8 :
            C01_ZR23[i]=np.log10(3.78e-4 * x2[i]**0.17 * x3[i]**(-0.44))+12.0    
         
        if OIIOIII5007[i] >= 0.8 :
            C01_ZR23[i]=np.log10(3.96e-4 * x3[i]**(-0.46))+12.0   

        
    # ## calculating Charlot 01 calibration: (case A) based on [NII]/[SII]##

    NIISII=10**logNIISII
    x2=OIIOIII5007/1.5
    x6=NIISII/0.85

    C01_Z=np.log10(5.09e-4*(x2**0.17)*(x6**1.17))+12




    # #### New ionization parameter and metallicity diagnostics #######
    # NEW R23 diagnostics from Kobulnicky & Kewley 
    # See NEW_fitfin.coefs

    nite=5  # number of iteations+1
    Z_new=np.zeros((nite,num))
    Z_new_low=np.zeros((nite,num))
    Z_new_up=np.zeros((nite,num))

    # trying out just using a single [NII]/Ha value

    #####KK04:

    #Z_init_KK04=fltarr(Snum)
    #Z_init_KK04(where(logNIIHa > -1.0))=8.9
    #Z_init_KK04(where(logNIIHa <= -1.0))=8.3
    #Z_new(0,*)=Z_init_KK04

    Z_new[0,:]=Z_init.copy()
    logq=np.zeros((nite,num))
    Zmax=np.zeros((nite,num))

    for i in range(nite-1) :

    # ionization parameter
        logq[i,:]=(32.81 -1.153*logOIIIOII**2 + Z_new[i,:]*(-3.396 -0.025*logOIIIOII + 0.1444*logOIIIOII**2))/(4.603  -0.3119*logOIIIOII -0.163*logOIIIOII**2+Z_new[i,:]*(-0.48 +0.0271*logOIIIOII+ 0.02037*logOIIIOII**2))

    # maximum of R23 curve:
        for j in range(num) :
            if G :
                #print 'new R23',j

                if (logq[i,j] >= 8.0) and (logq[i,j] < 8.3) :
                    Zmax[i,j]=8.4

                if (logq[i,j] >= 7.6) and (logq[i,j] < 8.0) :
                    Zmax[i,j]=8.4  # was 8.5

                if (logq[i,j] >= 7.15) and (logq[i,j] < 7.6) :
                    Zmax[i,j]=8.4  # was 8.6

                if (logq[i,j] >= 6.7) and (logq[i,j] < 7.15) :
                    Zmax[i,j]=8.4  # was 8.7

                
                Z_new_low[i+1,j]=(9.40 + 4.65*logR23[j]-3.17*logR23[j]**2)-logq[i,j]*(0.272+0.547*logR23[j]-0.513*logR23[j]**2)

                Z_new_up[i+1,j]=(9.72  -0.777*logR23[j]-0.951*logR23[j]**2 -0.072*logR23[j]**3-0.811*logR23[j]**4)-logq[i,j]*(0.0737  -0.0713*logR23[j]  -0.141 *logR23[j]**2+0.0373*logR23[j]**3 -0.058*logR23[j]**4)

                if Z_new[i,j] <= Zmax[i,j] :
                    Z_new[i+1,j]=(9.40 + 4.65*logR23[j]-3.17*logR23[j]**2)-logq[i,j]*(0.272+0.547*logR23[j]-0.513*logR23[j]**2)
                else:

                    Z_new[i+1,j]=(9.72  -0.777*logR23[j]-0.951*logR23[j]**2 -0.072*logR23[j]**3-0.811*logR23[j]**4)-logq[i,j]*(0.0737-0.0713*logR23[j]  -0.141 *logR23[j]**2+0.0373*logR23[j]**3 -0.058*logR23[j]**4)

    OIIIOII_R23_logq=logq.copy()

    KD03new_abund=np.zeros(num)
    KD03new_abund_low=np.zeros(num)
    KD03new_abund_up=np.zeros(num)
    logq_OIIIOII_final=np.zeros(num)
    KD03new_abund_R23=np.zeros(num)
    KD03new_abund_R23_low=np.zeros(num)
    KD03new_abund_R23_up=np.zeros(num)

    for j in range(num) :
        if G :
            #print 'R23 ',j
            if (raw_lines['[OIII]5007'][j] != 0.0) and (raw_lines['[OII]3727'][j] != 0.0) :
                if (raw_lines['Hb'][j] != 0.0) :
                    logq_OIIIOII_final[j]=logq[nite-1,j]
                    KD03new_abund[j]=Z_new[nite-1,j]
                    KD03new_abund_low[j]=Z_new_low[nite-1,j]
                    KD03new_abund_up[j]=Z_new_up[nite-1,j]
                    KD03new_abund_R23_low[j]=Z_new_low[nite-1,j]
                    KD03new_abund_R23_up[j]=Z_new_up[nite-1,j]
                    KD03new_abund_R23[j]=KD03new_abund[j]


    #2014 FED: changed moc value to None, not 0!
    KD03new_abund[(KD03new_abund_low > KD03new_abund_up)]=None
    KD03new_abund_R23[(KD03new_abund_low > KD03new_abund_up)]=None


    # calculating ionization parameter for those objects with [NII]/[OII],Ha & Hb
    # and [NII]/[OII] metallicity is > 8.4
    logq_final=np.zeros(num)

    for i in range(num) :
        if G :
            #print 'ionization parameter',i
            if (raw_lines['[NII]6584'][i] != 0.0) and (raw_lines['[OII]3727'][i] != 0.0) :
                if (raw_lines['Hb'][i] != 0.0) and (raw_lines['Ha'][i] != 0.0) :
                    if KD02_NIIOII_Z[i] > 8.4 :
                        logq_final[i]=(32.81 + 0.0*logOIIIOII[i]-1.153*logOIIIOII[i]**2 +KD02_NIIOII_Z[i]*(-3.396 -0.025*logOIIIOII[i] + 0.1444*logOIIIOII[i]**2))/(4.603  -0.3119*logOIIIOII[i] -0.163*logOIIIOII[i]**2+KD02_NIIOII_Z[i]*(-0.48 +0.0271*logOIIIOII[i]+ 0.02037*logOIIIOII[i]**2))
                    else:
                        logq_final[i]=logq[nite-1,i]
        
       
       
            if (raw_lines['[NII]6584'][i] == 0.0) and (raw_lines['[OII]3727'][i] != 0.0) :
                if (raw_lines['[OIII]5007'][i] != 0.0) :
                    if (raw_lines['Hb'][i] != 0.0) and (raw_lines['Ha'][i] != 0.0) :
                        logq_final[i]=logq[nite-1,i]
                        

    # ### KD02 [NII]/[OII] estimate ###
    # (can be used for for log(O/H)+12 > 8.6 only)
    # uses equation (5) from paper, this should be identical 
    # to the estimate above for abundances log(O/H)+12 > 8.6

    #KD02_NIIOII2_Z=np.zeros(num)

    KD02_NIIOII2_Z=np.log10(8.511e-4*(1.54020+1.26602*logNIIOII+0.167977*logNIIOII**2))+12.


    # if [NII]/[OII] after extinction correction is less than -1.5, then check the data.
    # if it is only slightly less than 1.5, then this can be a result of either noisy
    # data, inaccurate fluxes or extinction correction, or a higher ionization parameter
    # than modelled.  For these cases, the average of the M91,Z94 and C01 should be used.

    # KD02 R23 estimate (not reliable for  8.4 < log(O/H)+12 < 8.8)
    # uses [NII]/[OII] estimate as initial guess - this can be changed below

    KD02_R23_Z=np.zeros(num)

    # initializing:


    Zi=np.array([0.05,0.1,0.2,0.5,1.0,1.5,2.0,3.0])  # model grid abundances in solar units
    ZiOH=np.log10(Zi*8.511e-4)+12            # log(O/H)+12 units
    Zstep=np.array([0.025,0.075,0.15,0.35,0.75,1.25,1.75,2.5,3.5]) # middle of model grid abundances
    ZstepOH=np.log10(Zstep*8.511e-4)+12
    qstep=np.log10([3.5e6,7.5e6,1.5e7,3e7,6e7,1.16e8,2.25e8,3.25e8]) # model grid ionization parameters
    n_ite=3                      # number of iteations for abundance determination
    tol=1.0e-2                    # tolerance for convergance 
    R23_roots=np.zeros((4,num),dtype=complex)   # all possible roots of R23 diagnostic
    q_roots=np.zeros((3,num),dtype=complex)     # possible roots of q diagnostic
    q_R23=np.zeros((num,n_ite+1)) 
    q=np.zeros((num,n_ite+1))        # actual q value
    OIIIOII_coef=np.zeros((4,8))      # coefficients from model grid fits
    R23_coef=np.zeros((5,7))          # coefficients from model grid fits
    R23_Z=np.zeros((num,n_ite+1))    # Z value for each iteation

    R23_Z[:,0]=KD02_NIIOII_Z.copy()          # use NIIOII abundance as initial estimate
                                      # may be changed


    # occasionally, for noisy data or badly fluxed [OII].[OIII] or Hb, 
    # or for high ionization parameter galaxies, R23 is slightly higher
    # than the curves in our models - this will result in all complex roots of
    # the R23 curve unless a slightly lower R23 is used.  These should
    # be checked individually to make sure that it is just noise etc in the
    # data that is causing the problem, rather than wrong fluxes input.
    # the R23 ratio should be close to 0.95, not much more than 1.0 for the
    # data to be physical.

    #openw,1,'../../Dats/R23_check.txt'
    #for i=0,num-1 :
    #   if logR23[i] > 0.95 :
    #      logR23[i]=0.95
    #      printf,1,i,name[i],logR23[i]
    #   
    #end
    #close,1

    for i in range(num) :
        if G :
            #print 'old R23',i
            if (logOIII5007OII[i] != 0.0) and (logR23[i] != 0.0) :
      # coefficients from KD02 paper:

                R23_coef[:,0]=[-3267.93,1611.04,-298.187,24.5508,-0.758310] # q=5e6
                R23_coef[:,1]=[-3727.42,1827.45,-336.340,27.5367,-0.845876] # q=1e7
                R23_coef[:,2]=[-4282.30,2090.55,-383.039,31.2159,-0.954473] # q=2e7
                R23_coef[:,3]=[-4745.18,2309.42,-421.778,34.2598,-1.04411]  # q=4e7
                R23_coef[:,4]=[-4516.46,2199.09,-401.868,32.6686,-0.996645] # q=8e7
                R23_coef[:,5]=[-3509.63,1718.64,-316.057,25.8717,-0.795242] # q=1.5e8
                R23_coef[:,6]=[-1550.53,784.262,-149.245,12.6618,-0.403774] # q=3e8

                OIIIOII_coef[:,0]=[-36.9772,10.2838 ,-0.957421,0.0328614] #z=0.05 
                OIIIOII_coef[:,1]=[-74.2814,24.6206,-2.79194,0.110773]    # z=0.1
                OIIIOII_coef[:,2]=[-36.7948,10.0581,-0.914212,0.0300472]  # z=0.2
                OIIIOII_coef[:,3]=[-81.1880,27.5082,-3.19126,0.128252]    # z=0.5
                OIIIOII_coef[:,4]=[-52.6367,16.0880,-1.67443,0.0608004]   # z=1.0
                OIIIOII_coef[:,5]=[-86.8674,28.0455,-3.01747,0.108311]    # z=1.5
                OIIIOII_coef[:,6]=[-24.4044,2.51913,0.452486,-0.0491711]  # z=2.0
                OIIIOII_coef[:,7]=[49.4728,-27.4711,4.50304,-0.232228]    # z=3.0

                OIIIOII_coef[0,:]=OIIIOII_coef[0,:]-logOIII5007OII[i]
                
                R23_coef[0,:]=R23_coef[0,:]-logR23[i]

                for ite in range(1,n_ite+1) :                 # iteate if tolerance level not met
                    if abs(R23_Z[i,ite]-R23_Z[i,ite-1]) >= tol :

      #   calculate ionization parameter using [OIII]/[OII] with
      #   [NII]/[OII] abundance for the first iteation, and the R23
      #   abundance in consecutive iteations

                        for j in range(8):   
                            if R23_Z[i,ite-1] > ZstepOH[j] :
                                if R23_Z[i,ite-1] <= ZstepOH[j+1] :
                                    q_roots[:,i]=fz_roots(OIIIOII_coef[:,j])      


      #   q must be between 3.5e6 and 3.5e8 cm/s because that is the range it
      #   is defined over by the model grids, and it must be real.

                        for k in range(3) :
                            if (q_roots[k,i].imag) == 0.0 :
                                if (q_roots[k,i].real) >= 6.54407 :   #log units (q >= 1e6 cm/s) 
                                    if (q_roots[k,i].real) <= 8.30103 :   #log units (q <= 2e8 cm/s)
                                        q[i,ite]=(q_roots[k,i].real)
                                        q_R23[i,ite]=(q_roots[k,i].real)
            
       
      #   calculate abundance using ionization parameter:
                        R23_qstepno=0
                        
                        for j in range(7) :   
                            if q[i,ite] > qstep[j] :
                                if q[i,ite] <= qstep[j+1] :
                                    R23_roots[:,i]=fz_roots(R23_coef[:,j])
                                    R23_qstepno=j
          
     
      #   There will be four roots, two complex ones, and two real ones.
      #   use previous R23 value (or [NII]/[OII] if first iteation) 
      #   and q to find which real root to use (ie. which side of R23 curve 
      #   to use).  

                        nn=0
                        #   Rmax=[1.04967,1.06497,1.06684,1.06329,1.03844,0.991261,0.91655]
                        Smax=np.array([8.69020,8.65819,8.61317,8.58916,8.49012,8.44109,8.35907])

                        for k in range(4) :
                            if (R23_roots[k,i].imag) == 0.0 :
                                if (R23_Z[i,ite-1] >= Smax[R23_qstepno]) and ((R23_roots[k,i].real) >= Smax[R23_qstepno]) :
                                    R23_Z[i,ite]=R23_roots[k,i].real

                                if (R23_Z[i,ite-1] <= Smax[R23_qstepno]) and ((R23_roots[k,i].real) <= Smax[R23_qstepno]) :
                                    R23_Z[i,ite]=(R23_roots[k,i].real)



    # around maximum of R23 sometimes the R23 ratio will be slightly higher than
    # that available for the closest q value.  This will depend on noise addded to
    # data.  If this happens, step up in ionization parameter to find abundance
    # using that one instead. Around local maximum, the actual ionization parameter
    # used is not significant compared to the errors associated with the lack of
    # abundance sensitivity of the R23 ratio in this region.

                        while (R23_Z[i,ite] == 0.0) and (R23_qstepno <= 5) :
                            R23_roots[:,i]=fz_roots(R23_coef[:,R23_qstepno+1])
                            for k in range(4):
                                if (R23_roots[k,i].imag) == 0.0 :
                                    if (R23_Z[i,ite-1] >= Smax[R23_qstepno]) and ((R23_roots[k,i].real) >= Smax[R23_qstepno]) :
                                        R23_Z[i,ite]=R23_roots[k,i].real

                                    if (R23_Z[i,ite-1] <= Smax[R23_qstepno]) and ((R23_roots[k,i].real) <= Smax[R23_qstepno]) :
                                        R23_Z[i,ite]=(R23_roots[k,i].real)


      
                            R23_qstepno+=1
        

                    else:
                        R23_Z[i,ite]=R23_Z[i,ite-1]  
                        q[i,ite]=q[i,ite-1]
     
    KD02_R23_Z=R23_Z[:,n_ite].copy()
    q_OIIIOII=q[:,n_ite].copy()

    #  ### Combined \R23\ method outlined in KD02 paper Section 7. ###
    #  ie for objects with only [OII], [OIII], Hb available

    KD02_R23comb_Z=np.zeros(num)

    for i in range(num) :
        if G :
            #print 'Combined R23',i
            if Z94_Z[i] >= 9.0 : 
                KD02_R23comb_Z[i]=(KD02_R23_Z[i]+M91_Z[i]+Z94_Z[i])/3.  # my technique averaged with M91 and Z94

            if (KD02_R23comb_Z[i] <= 9.0) and (KD02_R23comb_Z[i] >= 8.5) :
                KD02_R23comb_Z[i]=(M91_Z[i]+Z94_Z[i])/2.                   # average of M91 and Z94

            if (Z94_Z[i] <= 9.0) and (Z94_Z[i] >= 8.5) :
                KD02_R23comb_Z[i]=(M91_Z[i]+Z94_Z[i])/2.                    # average of M91 and Z94

            if KD02_R23comb_Z[i] <= 8.5 :
                KD02_R23comb_Z[i]=KD02_R23_Z[i]                        # my technique

            if Z94_Z[i] <= 8.5 :
                KD02_R23comb_Z[i]=KD02_R23_Z[i]                        # my technique


    # ### [NII]/[SII] method outlined in KD02 paper ###
    # this method produces a systematic shift of 0.2 dex in log(O/H)+12
    # compared with the average of M91, Z94, and C01.  We believe this
    # is a result of inaccurate abundances or depletion factors, which are known 
    # problems in sulfur modelling.  Initial guess of [NII]/[OII] used
    # can be changed.  Initial guess is not critical except for high
    # ionization parameters.  ionization parameter diagnostic is [OIII]/[OII]

    KD02_NIISII_Z=np.zeros(num)


    n_ite=3                      # number of iteations for abundance determination
    tol=1.0e-2                    # tolerance for convergance 
    NIISII_roots=np.zeros((4,num),dtype=complex)   # all possible roots of NIISII diagnostic
    q_roots=np.zeros((3,num),dtype=complex)     # possible roots of q diagnostic
    q=np.zeros((num,n_ite+1))        # actual q value
    q_NIISII=np.zeros((num,n_ite+1))        # actual q value
    OIIIOII_coef=np.zeros((4,8))      # coefficients from model grid fits
    NIISII_coef=np.zeros((5,7))          # coefficients from model grid fits
    NIISII_Z=np.zeros((num,n_ite+1))    # Z value for each iteation

    # initializing:

    NIISII_Z[:,0]=KD02_NIIOII_Z.copy()  # use [NII]/[OII] abundance as initial estimate
                                      # may be changed

    for i in range(num) :
        if G :
            #print '[NII]/[SII]',i
            if (logOIII5007OII[i] != 0.0) and (logNIISII[i] != 0.0) :
      # coefficients from KD02 paper:

                NIISII_coef[:,0]=[-1042.47,521.076,-97.1578,8.00058,-0.245356]
                NIISII_coef[:,1]=[-1879.46,918.362,-167.764,13.5700,-0.409872]
                NIISII_coef[:,2]=[-2027.82,988.218,-180.097,14.5377,-0.438345]
                NIISII_coef[:,3]=[-2080.31,1012.26,-184.215,14.8502,-0.447182]
                NIISII_coef[:,4]=[-2162.93,1048.97,-190.260,15.2859,-0.458717]
                NIISII_coef[:,5]=[-2368.56,1141.97,-205.908,16.4451,-0.490553]
                NIISII_coef[:,6]=[-2910.63,1392.18,-249.012,19.7280,-0.583763]

                OIIIOII_coef[:,0]=[-36.9772,10.2838 ,-0.957421,0.0328614] #z=0.05 
                OIIIOII_coef[:,1]=[-74.2814,24.6206,-2.79194,0.110773]    # z=0.1
                OIIIOII_coef[:,2]=[-36.7948,10.0581,-0.914212,0.0300472]  # z=0.2
                OIIIOII_coef[:,3]=[-81.1880,27.5082,-3.19126,0.128252]    # z=0.5
                OIIIOII_coef[:,4]=[-52.6367,16.0880,-1.67443,0.0608004]   # z=1.0
                OIIIOII_coef[:,5]=[-86.8674,28.0455,-3.01747,0.108311]    # z=1.5
                OIIIOII_coef[:,6]=[-24.4044,2.51913,0.452486,-0.0491711]  # z=2.0
                OIIIOII_coef[:,7]=[49.4728,-27.4711,4.50304,-0.232228]    # z=3.0

                OIIIOII_coef[0,:]=OIIIOII_coef[0,:]-logOIII5007OII[i]
                NIISII_coef[0,:]=NIISII_coef[0,:]-logNIISII[i]

                for ite in range(1, n_ite+1):                 # iteate if tolerance level not met
                    if abs(NIISII_Z[i,ite]-NIISII_Z[i,ite-1]) >= tol :

      #   calculate ionization parameter using [OIII]/[OII] with
      #   [NII]/[OII] abundance for the first iteation, and the [NII]/[SII]
      #   abundance in consecutive iteations

                        for j in range(8) :   
                            if NIISII_Z[i,ite-1] > ZstepOH[j] :
                                if NIISII_Z[i,ite-1] <= ZstepOH[j+1] :
                                    q_roots[:,i]=fz_roots(OIIIOII_coef[:,j])

          

      #   q must be between 3.5e6 and 3.5e8 cm/s because that is the range it
      #   is defined over by the model grids, and it must be real.

                        for k in range(3) :
                            if (q_roots[k,i].imag) == 0.0 :
                                if (q_roots[k,i].real) >= 6.54407 :   #log units (q >= 1e6 cm/s) 
                                    if (q_roots[k,i].real) <= 8.30103 :   #log units (q <= 2e8 cm/s)
                                        q[i,ite]=(q_roots[k,i].real)
                                        q_NIISII[i,ite]=(q_roots[k,i].real)


      #   calculate abundance using ionization parameter:

                        for j in range(7) :   
                            if q[i,ite] > qstep[j] :
                                if q[i,ite] <= qstep[j+1] :
                                    NIISII_roots[:,i]=fz_roots(NIISII_coef[:,j])
                                    q_NIISII_used=j
              
     
      #   There will be four roots, two complex ones, and two real ones.
      #   use previous NIISII value (or [NII]/[OII] if first iteation) 
      #   and q to find which real root to use (ie. which side of R23 curve 
      #   to use).  

                        nn=0

                        for k in range(4) :
                            if (NIISII_roots[k,i].imag) == 0.0 :
                                if ((NIISII_roots[k,i].real) >= 8.0) and ((NIISII_roots[k,i].real) <= 9.35) :
                                    NIISII_Z[i,ite]=(NIISII_roots[k,i].real)


                        if NIISII_Z[i,ite] == 0.0 :
                            NIISII_roots[:,i]=fz_roots(NIISII_coef[:,q_NIISII_used+1])
                            for k in range(4) :
                                if (NIISII_roots[k,i].imag) == 0.0 :
                                    if ((NIISII_roots[k,i].real) >= 8.0) and ((NIISII_roots[k,i].real) <= 9.35) :
                                        NIISII_Z[i,ite]=(NIISII_roots[k,i].real)





                    else:
                        NIISII_Z[i,ite]=NIISII_Z[i,ite-1]  
                        q[i,ite]=q[i,ite-1]

    KD02_NIISII_Z=NIISII_Z[:,n_ite].copy()




    # KD01 combined method (uses [NII], [OII], [OIII], [SII]):
    # uses KD02 [NII]/[OII] method if [NII]/[OII] gives 8.6 < log(O/H)+12
    # uses average of M91 and Z94 if 8.5 < log(O/H)+12 < 8.6
    # uses average of C01 and KD02 if  log(O/H)+12 < 8.5
    # Also calculates comparison average

    KD02_comb_Z=np.zeros(num)
    KD02C01_ave=np.zeros(num)
    M91Z94C01_ave=np.zeros(num)
    M91Z94_ave=np.zeros(num)
    name=np.zeros(num)

    for i in range(num) :
        name[i]=i
        if G :
            #print 'KD01 combined',i
            if (M91_Z[i] != 0.0) and (Z94_Z[i] != 0.0) :
                M91Z94_ave[i]=(M91_Z[i]+Z94_Z[i])/2.
                if (C01_Z[i] != 0.0) :
                    M91Z94C01_ave[i]=(M91_Z[i]+Z94_Z[i]+C01_Z[i])/3.

      
            if (KD02_R23_Z[i] != 0.0) and (C01_Z[i] != 0.0) :
                KD02C01_ave[i]=(KD02_R23_Z[i]+C01_Z[i])/2
      
            if KD02_NIIOII_Z[i] <= 8.6 :
                if M91Z94_ave[i] >= 8.5 :
                    KD02_comb_Z[i]=M91Z94_ave[i]   # average of M91 and Z94
                else:
                    KD02_comb_Z[i]=KD02C01_ave[i]
         
            else:
                KD02_comb_Z[i]=KD02_NIIOII_Z[i]
     


    #-----------------------------------------
    # ### NEW May 25 - NEW combined method ###
    #-----------------------------------------

    # if [NII]/[OII] abundance available and [NII]/Ha abundance < 8.4, then 
    # use R23. 

    KD_comb_NEW=np.zeros(num)

    for j in range(num) :
        if G :
            #print 'NEW combined',j
            if (Z_init_guess[j] > 8.4) :
                KD_comb_NEW[j]=KD02_NIIOII_Z[j]
            else:
                if (KD03new_abund_R23[j] > 0.0) and (M91_Z[j] > 0.0 ):
                    KD_comb_NEW[j]=(KD03new_abund_R23[j]+M91_Z[j])/2.
                else:
                    KD_comb_NEW[j]=KD03_NIIHa_abund[j] 



    ###############################################
    ################ PRINTING OUT #################
    ###############################################
    if saveres:
        wf=open(outfilename,'a')
    if disp:
        print 'ID\tKD CombNEW\tKD02[NII]/[OII]\tR23 NEW Comb\tM91\tZ94'
        print '\tKD02[NII]/Ha\tD02\tPP04_N2\tPP04_O3N2'
    for i in range(num):
        line1=str(i+1)+'\t'+str(KD_comb_NEW[i])+'\t'+str(KD02_NIIOII_Z[i])+'\t'+str(KD03new_abund_R23[i])+'\t'+str(M91_Z[i])+'\t'+str(Z94_Z[i])
        line2='\t'+str(KD03_NIIHa_abund[i])+'\t'+str(D02_Z[i])+'\t'+str(PP04_N2_Z[i])+'\t'+str(PP04_O3N2_Z[i])+'\t'+str(Pi01_Z[i])
        if disp:
            print line1 
            print line2
        if saveres:
            wf.write(line1+line2)
            wf.write('\n')
    if saveres:
        wf.close()

    res={Zs[0]:KD_comb_NEW, Zs[1]:KD02_NIIOII_Z, Zs[2]:KD03new_abund_R23,
         Zs[3]:M91_Z, Zs[4]:Z94_Z, Zs[5]:KD03_NIIHa_abund, Zs[6]:D02_Z,
         Zs[7]:PP04_N2_Z, Zs[8]:PP04_O3N2_Z, Zs[9]:Pi01_Z,Zs[10]:EB_V}
    return res,red_corr
    #return [KD_comb_NEW,KD02_NIIOII_Z,KD03new_abund_R23,M91_Z,Z94_Z,KD03_NIIHa_abund,D02_Z,PP04_N2_Z,PP04_O3N2_Z,Pi01_Z]

'''

