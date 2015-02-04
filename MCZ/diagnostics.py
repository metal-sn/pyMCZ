import numpy as np
import pylab as pl
niter=5  # number of iteations+1 for KD02 methods

k_Ha=2.535  # CCM Rv=3.1
k_Hb=3.609  # CCM Rv=3.1

k_O1=2.66146  # CCM Rv=3.1
k_O2=4.771 # CCM Rv=3.1
k_O35007=3.341 # CCM Rv=3.1
k_O34959=3.384 # CCM Rv=3.1
k_O3=(k_O35007+k_O34959)/2.

k_N2=2.44336 # CCM Rv=3.1
k_S2=2.38089 # CCM Rv=3.1




class diagnostics:
    def __init__(self,num):
        self.nm=num

        self.Ha,self.Hb=None,None
        self.hasHa,self.hasHb=False,False

        #metallicity diagnostics to be returned
        self.mds={
            'E(B-V)':None,
            #
            ###Kobulnicky parameterization of Zaritzky et al. (1994) 
            ###Z94 = Zaritsky, D., Kennicutt, R. C., & Huchra, J. P. 1994, 
            ###ApJ, 420, 87
            'Z94':None,
            #
            ###McGaugh, S.S., 1991, ApJ, 380, 140'
            'M91':None,         
            #
            ### D02_Z Denicolo [NII]/Ha diagnostic Denicolo, Terlevich & Terlevich (2002),  MNRAS, 330, 69
            'D02':None,
            #
            ### PP04_N2_Z, PP04_O3N2_Z Pettini & Pagel diagnostics - Pettini & Pagel (2004), MNRAS, 348, L59
            'PP04_N2':None,
            'PP04_O3N2':None,
            ##Pilyugin method.  Based on [OIII],[OII], Hbeta 
            ##calibrated from Te method'
            'Pi01':None,
            #
            ### Kewley & Dopita (2002) (KD02) estimates of abundance ##
            'KD02_N2O2':None, 
            #
            #KD02 = Kewley, L. J., & Dopita, M. A., 2002, ApJ, submitted '
            'KD02comb_updated':None, 
            #
            #Kewley, L. J., & Dopita, M. A., 2003 
            'KD03new_R23':None,
            'KD03_N2Ha':None,
            #
            #C01 = Charlot, S., & Longhetti, M., 2001, MNRAS, 323, 887'
            'C01_R23':None,
            'C01':None
        }


        #other lines calculated and repeatedly used
        self.R23=None
        self.logR23=None
        self.logN2O2=None
        self.logO3O2=None
        self.N2=None
        self.O23727=None
        self.N26584=None
        self.O35007=None
        self.N2O2_coef0=1106.8660

    def printme(self):
        try: print "\nHa", np.mean(self.Ha),self.Ha
        except:pass
        try:print "\nHb", np.mean(self.Hb),self.Hb
        except:pass
        try:print  "\nE(B-V)",np.mean(self.mds['E(B-V)']),self.mds['E(B-V)']
        except:pass 
        try:       print  "\nZ94_Z",np.mean(self.mds['Z94']),self.mds['Z94']
        except:pass 
        try:       print  "\nR23",np.mean(self.R23),self.R23
        except:pass 
        try:       print "\nlog(R23)", np.mean(self.logR23),self.logR23
        except:pass 
        try:        print  "\nlog([NII][OII])",np.mean(self.logN2O2),self.logN2O2
        except:pass
        try:        print  "\nlog([OIII][OII])",np.mean(self.logO3O2),self.logO3O2
        except:pass

    

    def fz_roots(self):        
        rts=np.zeros((self.nm,4),dtype=complex)
        self.N2O2_coef[np.where(~np.isfinite(self.N2O2_coef))]=0.0
        for i in range(self.nm):

            rts[i]= np.roots(self.N2O2_coef[i][::-1])#::-1][0])
        if rts.size==0:
            print 'fz_roots failed'
        return rts

    def dustcorrect(self,l1,l2,flux=False):
        if not flux:
            return 0.4*self.mds['E(B-V)']*(l1-l2)
        return 10**(0.4*self.mds['E(B-V)']*(l1-l2))        

    def setHab(self,Ha,Hb):
        self.Ha=Ha
        self.Hb=Hb
        if (self.Ha > 0).any():
            self.hasHa=True
        if (self.Hb > 0).any():
            self.hasHb=True

    def setOlines(self, O23727, O35007, O16300, O34959):
        self.O23727 =O23727
        self.O35007=O35007
        if (self.O35007>0).any() : self.hasO35007=True
        if (self.O23727>0).any() : self.hasO23727=True
        if self.hasO23727 and self.hasO35007:
            O35007toO23727=self.O35007/self.O23727
            O23727toO35007=self.O23727/self.O35007
            self.logO35007O2=np.log10( O35007toO23727 )+self.dustcorrect(k_O3,k_O2)
            self.logO2O35007=np.log10( O23727toO35007 )+self.dustcorrect(k_O2,k_O3)
            #self.logO2O35007Hb=np.log10((self.O23727+self.O35007)/self.Hb)
            # ratios for other diagnostics - slightly different ratios needed
            self.O35007O2=10**(self.logO35007O2)
            self.O2O35007=10**(self.logO2O35007)
            self.hasO35007O2 =True
            if self.hasHb:
                self.logO2O35007Hb=np.log10((self.O23727/self.Hb)* self.dustcorrect(k_O2,k_Hb,flux=True))+ (self.O35007/self.Hb)*self.dustcorrect(k_O35007,k_Hb,flux=True)

        if self.hasHb and self.hasO23727:
            self.logO2Hb=np.log10(self.O23727/self.Hb)+self.dustcorrect(k_O2,k_Hb)#0.4*self.mds['E(B-V)']*(k_O2-k_Hb) 

        # never used
        #if self.hasHa:
            #logO1Ha=np.log10(O16300/self.Ha)+self.dustcorrect(k_O1,k_Ha)


        if self.hasO35007:
            if self.hasHb :
                self.logO3Hb=np.log10(self.O35007/self.Hb)+self.dustcorrect(k_O35007,k_Hb)#0.4*self.mds['E(B-V)']*(k_O35007-k_Hb)
                self.hasO3Hb=True
                #if self.hasO23727:
                #self.logO2O35007Hb=np.log10((self.O23727+self.O35007)/self.Hb)
            
            if (O34959>0).any():
                self.logO349595007Hb=np.log10(10**(np.log10(self.O35007/self.Hb)+self.dustcorrect(k_O35007,k_Hb))+10**(np.log10(O34959/self.Hb)+self.dustcorrect(k_O34959,k_Hb)))
                self.O34959p5007=O34959 + self.O35007
                if self.hasO23727:
                    self.logO3O2=np.log10((self.O34959p5007)/self.O23727)+self.dustcorrect(k_O3,k_O2)        
                    self.hasO3O2=True

        #O3O2[i]=1.347*O35007O2[i]
        #logO3O2[i]=np.log10(O3O2[i])

    def setNII(self,N26584):
        if (N26584>0).any():
            self.N26584=N26584
            self.hasN2=True
            
    def setS(self,S26717,S26731=[0],S39069=[0],S39532=[0]):
        if (S26717>0).any():
            self.S26717=S26717
            self.hasS26717=True
        if (S26731>0).any():
            self.S26731=S26731
            self.hasS26731=True
        if (S39069>0).any():
            self.S39069=S39069
            self.hasS39069=True
        if (S39532>0).any():
            self.S39532=S39532
            self.hasS39532=True

    def calcEB_V(self):
        #logHaHb=np.log10(Ha/Hb)
        self.mds['E(B-V)']=np.log10(2.86/(self.Ha/self.Hb))/(0.4*(k_Ha-k_Hb)) # E(B-V)
        #print self.mds['E(B-V)']
        self.mds['E(B-V)'][self.mds['E(B-V)']<=0]=0.00001

    def calcNIIOII(self):
        if self.hasN2 and self.hasO23727 :
            self.hasN2O2=True
            self.logN2O2=np.log10(self.N26584/self.O23727)+self.dustcorrect(k_N2,k_O2) 
            self.N2O2_coef=np.array([[0,-532.15451,96.373260,-7.8106123,0.23928247]]*self.nm).T# q=2e7 line (approx average)
            self.N2O2_coef[0]=self.N2O2_coef[0]+self.N2O2_coef0-self.logN2O2
            self.N2O2_coef=self.N2O2_coef.T
            # finding roots for == (4)
            self.N2O2_roots=np.array([self.fz_roots()])[0]          

    def calcR23(self):
        #R23 NEW Comb, [NII]/Ha: KK04 = Kobulnicky & Kewley, 2004, submitted'
        if  self.hasO3O2 and self.hasO23727 and self.hasHb:
            self.R23=((self.O23727/self.Hb)*self.dustcorrect(k_O2,k_Hb, flux=True) + (self.O34959p5007/self.Hb)*self.dustcorrect(k_O3,k_Hb, flux=True) )
            self.logR23=np.log10(self.R23)

            if self.hasO35007O2:
                # R23 but without [O3]4959
                self.R23_5007=(1./self.O35007O2 + 1.)/(1./self.O35007O2 + 1.347)*self.R23  
        
    def calcS23(self):
        #the original code here uses S267176731, which is however set to 6717 as default
        if  self.hasS26717 :
            if self.hasS39069 and self.hasHb:
                self.logS23=np.log10((self.S26717/self.Hb)*self.dustcorrect(k_S2,k_Hb,flux=True) + (S39069/self.Hb)*self.dustcorrect(k_S3,k_Hb,flux=True))                                 
            self.logS3S2=np.log10(S39069/self.S26717)+self.dustcorrect(k_S3,k_S2)


    def calcNIISII(self):
        if self.hasS26717 and self.hasN2:
            self.N2S2=self.N26584/self.S26717+self.dustcorrect(k_N2,k_S2,flux=True)#0.4*self.mds['E(B-V)']*(k_N2-k_S2) 
            self.logN2S2=np.log10(self.N26584/self.S26717)+self.dustcorrect(k_N2,k_S2)#0.4*self.mds['E(B-V)']*(k_N2-k_S2) 
            self.hasN2S2=True
        else: 
            print "WARNING: needs SII6717 and NII6584 to calculate NIISII: did you run setN2() and setS?"

        if self.hasHa and self.hasN2:
            self.logN2Ha = np.log10(self.N26584/self.Ha)+self.dustcorrect(k_N2,k_Ha)#0.4*self.mds['E(B-V)']*(k_N2-k_Ha) 
            self.hasN2Ha=True
        else: 
            print "WARNING: needs NII6584 and Ha to calculate NIIHa: did you run setHab()?"

    def calcSIIHa(self):
        if self.hasS26717 and self.hasHa:
            self.logS2Ha=np.log10(self.S26717/self.Ha)+self.dustcorrect(k_S2,k_Ha)#0.4*self.mds['E(B-V)']*(k_S2-k_Ha) 
        else: 
            print "WARNING: needs SII6717 and Ha to calculate SIIHa: did you run setHab() and setS()?"



#######################these are the metallicity diagnostics##################

    def calcD02_Z(self):
        if self.hasN2Ha:
            self.mds['D02'] = 9.12+0.73*self.logN2Ha
        else:
            print "WARNING: need N2Ha to do this. did you run calclogNIISI"

    def calcPP04_N2_Z(self):
        if self.hasN2Ha:
            self.mds['PP04_N2']= 9.37 + 2.03*self.logN2Ha + 1.26*self.logN2Ha**2 + 0.32*self.logN2Ha**3
        else:
            print "WARNING: need N2Ha to do this. did you run calclogNIISI"

    def calcZ94_Z(self):
        ### calculating z from Kobulnicky parameterization of Zaritzky et al. (1994) 
        if self.logR23==None:
            print "Must first calculate R23"
        else:
            self.mds['Z94']=9.265-0.33*self.logR23-0.202*self.logR23**2-0.207*self.logR23**3-0.333*self.logR23**4         
            self.mds['Z94'][(self.logR23 > 0.9)]=None
    
    def Pmethod(self):
        # #### P-method #####
        #make sure you run calclogOOs and calclogO2Hb first
        R3=10**self.logO349595007Hb
        R2=10**self.logO2Hb
        P = R3/(R2+R3)
        P_R23=R2+R3
        
        # NEW Initial Guess - Nov 1 2006
        self.Z_init_guess=np.zeros(self.nm)+8.7 # upper branch if nothing else
        # [N2]/Ha if no [N2]/[O2]
        if self.hasHa:
            self.Z_init_guess[(self.logN2Ha < -1.3)&(self.N26584 != 0.0)&(self.hasHa != 0)]=8.2
            self.Z_init_guess[(self.logN2Ha < -1.1)&(self.N26584 != 0.0)&(self.hasHa != 0.0)&(self.logN2Ha >= -1.3)]=8.4
            self.Z_init_guess[(self.logN2Ha >=-1.1)&(self.N26584 != 0.0)&(self.hasHa != 0.0)]=8.7
            N2O2_lines=self.N26584*self.Ha*self.Hb*self.O23727

            if self.hasN2O2:
                # [N2]/[O2] if at all possible
                self.Z_init_guess[(self.logN2O2 < -1.2)&(N2O2_lines != 0.0)]=8.2  #  1.2 using low-Z gals,
                self.Z_init_guess[(self.logN2O2 >= -1.2)&(N2O2_lines != 0.0)]=8.7  # 1.1 using HIi regions
        
        P_abund_up=(P_R23+726.1+842.2*P+337.5*P**2)/(85.96+82.76*P+43.98*P**2+1.793*P_R23)
        P_abund_low=(P_R23+106.4+106.8*P-3.40*P**2)/(17.72+6.60*P+6.95*P**2-0.302*P_R23)

        self.mds['Pi01']=P_abund_up.copy()
        self.mds['Pi01'][(self.Z_init_guess < 8.4)]=P_abund_low[(self.Z_init_guess < 8.4)].copy()
        self.mds['Pi01'][(self.Z_init_guess >= 8.4)]=P_abund_up[(self.Z_init_guess >= 8.4)].copy()

    def calcPi01_Z_old(self):
        # P-method 2001 upper branch (derecated and commented out)
        # available but deprecated
        if self.hasO2O2:
            P = 10**self.logO3O2/(1+10**self.logO3O2)
            P_R23=10**self.logR23
            P_abund_old=(P_R23+54.2+59.45*P+7.31*P**2)/(6.07+6.71*P+0.371*P**2+0.243*P_R23)
            self.mds['Pi01_old']=np.zeros(self.nm)
            self.mds['Pi01_old'][Z_init_guess >= 8.4]=P_abund_old[Z_init_guess >= 8.4]
        else:
            print "WARNING: need OIIIOII to calculate Pi01_Z_old, did you set them up with  setOlines()?"
        
    def calcC01_ZR23(self):
        # Charlot 01 R23 calibration: (case F) ##        
        # available but deprecated
        if self.hasO35007O2 and self.hasO3Hb:
            x2=self.O2O35007/1.5
            x3=(10**self.logO3Hb)/2.
            self.mds['C01_R23']=np.zeros(self.nm)        
            self.mds['C01_R23'][self.O2O35007<0.8]=np.log10(3.78e-4 * (self.O2O35007[self.O2O35007<0.8]/1.5)**0.17 * x3[self.O2O35007<0.8]**(-0.44))+12.0    
         
            self.mds['C01_R23'][ self.O2O35007 >= 0.8]=np.log10(3.96e-4 * x3[self.O2O35007 >= 0.8]**(-0.46))+12.0   
        else:
            print "WARNING: need [OIII]5700, [OII]3727, and Ha to calculate calcC01_ZR23, did you set them up with  setOlines()?"        

        # Charlot 01 calibration: (case A) based on [N2]/[SII]##
        # available but deprecated
        if self.hasN2S2:
            self.mds['C01']=np.log10(5.09e-4*((self.O2O35007/1.5)**0.17)*(((10**self.logN2S2)/0.85)**1.17))+12
        else:
            print "WARNING: need [OIII]5700, [OII]3727, and Ha to calculate calcC01_ZR23, did you set them up with  setOlines() and ?"        



    def calcPP04_O3N2_Z(self):
        if self.hasO3Hb and self.hasN2Ha:
            self.mds['PP04_O3N2']=8.73 - 0.32*(self.logO3Hb-self.logN2Ha)
        else:
            print "WARNING: need O3Hb and N2Ha to calculate PP04_O3N2, did you set them up?"

    def calcM91(self):
        # ## calculating M91 calibration using [N2O2] as 
        # initial estimate of abundance ##:
        # this initial estimate can be changed by replacing 
        # OH_init by another guess, eg C01_Z 
        # NOTE: occasionally the M91 'upper branch' will give a metallicity
        # that is lower than the 'lower branch'.  Happens for very high R23
        # values.  Probably because the R23 value is higher than the R23 point
        # at which the upper & lower branches intersect.  What to do?  if R23
        # is higher than the intersection (calculate the intersection), then
        # the metallicity is likely to be around the R23 maximum = 8.4
        
        # Z_init=Z94_Z        
        # Z_init=KD02N2O2_Z

        self.mds['M91']=np.zeros(self.nm)
        M91_Z_low=12.0-4.944+0.767*self.logR23+0.602*self.logR23**2-self.logO3O2*(0.29+0.332*self.logR23-0.331*self.logR23**2)
        M91_Z_up=12.0-2.939-0.2*self.logR23+-0.237*self.logR23**2-0.305*self.logR23**3-0.0283*self.logR23**4-self.logO3O2*(0.0047-0.0221*self.logR23-0.102*self.logR23**2-0.0817*self.logR23**3-0.00717*self.logR23**4)
        indx,=np.where((np.abs(self.logO3O2)>0) & (np.abs(self.logR23)>0) & (self.Z_init_guess < 8.4))

        self.mds['M91'][indx]=12.0-4.944+0.767*self.logR23[indx]+0.602*self.logR23[indx]**2-self.logO3O2[indx]*(0.29+0.332*self.logR23[indx]-0.331*self.logR23[indx]**2)
        indx,=np.where((np.abs(self.logO3O2)>0) & (np.abs(self.logR23)>0) & (self.Z_init_guess >= 8.4))
        self.mds['M91'][indx]=12.0-2.939-0.2*self.logR23[indx]-0.237*self.logR23[indx]**2-0.305*self.logR23[indx]**3-0.0283*self.logR23[indx]**4-self.logO3O2[indx]*(0.0047-0.0221*self.logR23[indx]-0.102*self.logR23[indx]**2-0.0817*self.logR23[indx]**3-0.00717*self.logR23[indx]**4)

        #2014 FED: changed wrong values to None
        self.mds['M91'][(M91_Z_up < M91_Z_low)]=None

    def calcKD02_N2O2_Z(self):
        ### Kewley & Dopita (2002) (KD02) estimates of abundance ##
        ##  KD02
        # KD02 [N2]/[O2] estimate (can be used for whole log(O/H)+12 range, 
        # but rms scatter increases to 0.11 rms for log(O/H)+12 < 8.6
        # rms = 0.04 for
        # log(O/H)+12 > 8.6
        # uses equation (4) from paper:
        self.mds['KD02_N2O2']=np.zeros(self.nm)
        if self.hasN2 and self.hasO23727 and self.hasHa and self.hasHb:         
            for i in range(self.nm):
                k=3
                while k>=0:
                    if not (self.N2O2_roots[i][k].imag) == 0.0 :        
                        k=k-1
                    else:
                        # between 7.5 and 9.4
                        if (abs(self.N2O2_roots[i][k]) >= 7.5) and (abs(self.N2O2_roots[i][k]) <= 9.4) :
                            self.mds['KD02_N2O2'][i]=abs(self.N2O2_roots[i][k]) 
                            k=0
                        k=k-1
        else:
            print "WARNING: need NII6584 and OII3727 and Ha and Hb to calculate this. did you run setO() setHab() and setNII()?"

    def calcKD03_NHa(self):
        # calculating [N2]/Ha abundance estimates using [O3]/[O2] also
        if self.mds['KD02_N2O2'] ==None:
            self.calcKD02_N2O2_Z()
            if self.mds['KD02_N2O2'] ==None:
                print "WARNING: without KD02_N2O2_Z cannot calculate KD03_NHa"
                return -1
        ##why iterating???Z_new_N2Ha=np.zeros((niter,self.nm))
        Z_new_N2Ha=self.mds['KD02_N2O2'].copy()  # was 8.6

        '''
        ##why iterating???
        for i in range(niter-1):
            # ionization parameter
            if (self.N26584 != 0.0).any() :
                if self.hasHa and self.hasHb:
                    if self.hasO3O2 :
                        # calculating logq using the [N2]/[O2] 
                        #metallicities for comparison
                        #used to include 0.0*logO3O2(j)
                        logq[i]=(32.81 -1.153*self.logO3O2**2 + 
                                 KD02_N2O2_Z*(-3.396 -0.025*self.logO3O2 + 0.1444*self.logO3O2**2))/(4.603-0.3119*self.logO3O2 -0.163*self.logO3O2**2+KD02_N2O2_Z*(-0.48 + 0.0271*self.logO3O2+ 0.02037*self.logO3O2**2))
                    else :
                        logq[i]=([7.37177]*self.nm)
                        #       no_O3O2_yes_HaHb_flag(noO3O2)=j
                else:
                    logq[i]=([7.37177]*self.nm)
#                    logN2Ha=np.log10(self.N26584/self.Ha)
                Z_new_N2Ha[i+1]=(7.04 + 5.28*self.logN2Ha+6.28*self.logN2Ha**2+2.37*self.logN2Ha**3)-logq[i]*(-2.44-2.01*self.logN2Ha-0.325*self.logN2Ha**2+0.128*self.logN2Ha**3)+10**(self.logN2Ha-0.2)*logq[i]*(-3.16+4.65*self.logN2Ha)
       
                
                
        #N2Halogq=logq
        self.mds['KD03_N2Ha']=Z_new_N2Ha[niter-1,:]
        '''

        if self.hasN2 and self.hasHa:
            if self.hasO3O2 :
                # calculating logq using the [N2]/[O2] 
                #metallicities for comparison
                #used to include 0.0*logO3O2(j)
                self.logq=(32.81 -1.153*self.logO3O2**2 + 
                       self.mds['KD02_N2O2']*(-3.396 -0.025*self.logO3O2 + 0.1444*self.logO3O2**2))/(4.603-0.3119*self.logO3O2 -0.163*self.logO3O2**2+ self.mds['KD02_N2O2']*(-0.48 + 0.0271*self.logO3O2+ 0.02037*self.logO3O2**2))
            else :
                self.logq=7.37177
                #       no_O3O2_yes_HaHb_flag(noO3O2)=j
            Z_new_N2Ha=(7.04 + 5.28*self.logN2Ha+6.28*self.logN2Ha**2+2.37*self.logN2Ha**3)-self.logq*(-2.44-2.01*self.logN2Ha-0.325*self.logN2Ha**2+0.128*self.logN2Ha**3)+10**(self.logN2Ha-0.2)*self.logq*(-3.16+4.65*self.logN2Ha)
            
                
            self.mds['KD03_N2Ha']=Z_new_N2Ha

        else:
            print "WARNING: need NII6584  and Ha to calculate this. did you run  setHab() and setNII()?"

    def calclims():
        # calculating upper and lower metallicities for objects without
        # Hb  and for objects without [O3] and/or [O2]
        Hb_up_ID=np.zeros(100)
        logq_low=6.9
        logq_up=8.38
        if self.hasN2 and self.hasHa:
            logN2Ha=np.log10(self.N26584/self.Ha)
            # ionization parameter
            Z_new_N2Ha_low=(7.04 + 5.28*logN2Ha+6.28*logN2Ha**2+2.37*logN2Ha**3)-logq_low*(-2.44-2.01*logN2Ha-0.325*logN2Ha**2+0.128*logN2Ha**3)+10**(logN2Ha-0.2)*logq_low*(-3.16+4.65*logN2Ha)
            Z_new_N2Ha_up= (7.04 + 5.28*logN2Ha+6.28*logN2Ha**2+2.37*logN2Ha**3)-logq_up *(-2.44-2.01*logN2Ha-0.325*logN2Ha**2+0.128*logN2Ha**3)+10**(logN2Ha-0.2)*logq_up *(-3.16+4.65*logN2Ha)
            '''
            this is identical to  the paragraph above!!
            if self.hasO35007 and self.hasO23727 :
                Z_new_N2Ha_low=(7.04 + 5.28*logN2Ha+6.28*logN2Ha**2+2.37*logN2Ha**3)-logq_low*(-2.44-2.01*logN2Ha-0.325*logN2Ha**2+0.128*logN2Ha**3)+10**(logN2Ha-0.2)*logq_low*(-3.16+4.65*logN2Ha)
                
                Z_new_N2Ha_up=(7.04 + 5.28*logN2Ha+6.28*logN2Ha**2+2.37*logN2Ha**3)-logq_up*(-2.44-2.01*logN2Ha-0.325*logN2Ha**2+0.128*logN2Ha**3)+10**(logN2Ha-0.2)*logq_up*(-3.16+4.65*logN2Ha)
            '''



    def calcKK04(self):
        if self.mds['KD02_N2O2']=None:
            # ### KD02 [NII]/[OII] estimate ###
            # (can be used for for log(O/H)+12 > 8.6 only)
            # uses equation (5) from paper, this should be identical 
            # to the estimate above for abundances log(O/H)+12 > 8.6
            
            self.mds['KD02_N2O2']=np.log10(8.511e-4*(1.54020+1.26602*self.logN2O2+0.167977*self.logN2O2**2))+12.


        # ionization parameter
        
        logq_final=np.zeros(self.nm)
        if self.hasN2 and self.hasO23727 and self.hasHb and self.hasHa:
            logq_final=(32.81 + 0.0*self.logO3O2-1.153*self.logO3O2**2 +self.mds['KD02_N2O2']*(-3.396 -0.025*self.logO3O2 + 0.1444*self.logO3O2**2))/(4.603  -0.3119*self.logO3O2 -0.163*self.logO3O2**2+self.mds['KD02_N2O2']*(-0.48 +0.0271*self.logO3O2+ 0.02037*self.logO3O2**2))
            logq_final[self.mds['KD02_N2O2']<=8.4]=self.logq[self.mds['KD02_N2O2']<=8.4]

        if  not self.hasN2 and self.hasO23727 and self.hasO35007_raw and self.hasHb and self.hasHa:
            logq_final=logq

