import numpy as np
import pylab as pl
import sys
import scipy.stats as stats
import numpy.polynomial.polynomial as nppoly
from fedmetallicity import get_keys
niter=5  # number of iteations+1 for KD02 methods

k_Ha=2.535  # CCM Rv=3.1
k_Hb=3.609  # CCM Rv=3.1

#k_O1=2.66146   # CCM Rv=3.1
k_O2=4.771     # CCM Rv=3.1
k_O35007=3.341 # CCM Rv=3.1
k_O34959=3.384 # CCM Rv=3.1
k_O3=(k_O35007+k_O34959)/2.

k_N2=2.44336 # CCM Rv=3.1
k_S2=2.38089 # CCM Rv=3.1

MINMASS= 6.0
MAXMASS=14.0
DUSTCORRECT=True

class diagnostics:
    def __init__(self,num):
        self.nm=num
        self.Ha=None
        self.Hb=None
        self.hasHa,self.hasHb=False,False
        self.hasO2,self.hasO3=False,False
        self.hasS2,self.hasN2=False,False

        self.hasO3Hb=False
        self.hasO3O2=False

        self.hasN2O2=False
        self.hasN2S2=False

        self.hasS26731=False
        self.hasS39532=False
        self.hasS39069=False

        self.N2O2_roots=None
        #other lines calculated and repeatedly used
        self.R23=None
        #self.R23_5007=None
        self.N2=None
        self.N2S2=None
        self.O23727=None
        self.O35007=None
        self.N26584=None
        self.S26717=None
        self.S26731=None
        self.S39069=None
        self.S39532=None

        self.O34959p5007=None
        self.O35007O2=None
        self.O2O35007=None

        self.logR23=None
        self.logN2O2=None
        self.logO3O2=None
        self.logS23=None
        self.logS3S2=None

        self.logO2Hb=None
        self.logO3Hb=None
        self.logN2Ha=None
        self.logS2Ha=None
        self.logO3O2=None
        self.logO35007O2=None
        self.logO2O35007=None
        self.logO349595007Hb=None

        self.logq=None
        self.Z_init_guess=None
        self.N2O2_coef0=1106.8660

        self.OIII_OII=None
        self.OIII_Hb=None
        self.OIII_SII=None

        self.NII_OII=None
        self.NII_SII=None



        #metallicity diagnostics to be returned
        
        self.mds={}
        for Z in  get_keys():
            self.mds[Z]=None
        '''    'E(B-V)':None,
            'logR23':None,
            'Z94':None,
            'M91':None,         
            'D02':None,
            'PP04_N2':None,
            'PP04_O3N2':None,
            'Pi01':None,
            'KD02_N2O2':None, 
            'KD02comb_updated':None, 
            'KD02new_R23':None,
            'KD02_N2Ha':None,
            'C01_R23':None,
            'C01':None
        }
        '''



    def printme(self, verbose=False):
        try: 
            print "\nHa", np.mean(self.Ha)
            if verbose: print self.Ha
        except:pass
        try:
            print "\nHb", np.mean(self.Hb)
            if verbose: print self.Hb
        except:pass
        try:
            print  "\nO2",np.mean(self.O23727)
            if verbose: print self.O23727
        except:pass 
        try:
            print  "\nO3",np.mean(self.O35007)
            if verbose: print self.O35007
        except:pass 
        try:
            print  "\nO34959",np.mean(self.O34959)
            if verbose: print self.O34959
        except:pass 
        try:       
            print  "\nZ94",np.mean(self.mds['Z94'])
            if verbose: print self.mds['Z94']
        except:pass 
        try:       
            print  "\nR23",np.mean(self.R23)
            if verbose: print self.R23
        except:pass 
        try:       
            print "\nlog(R23)", np.mean(self.logR23)
            if verbose: print self.logR23
        except:pass 
        try:        
            print  "\nlog([NII][OII])",stat.nanmean(self.logN2O2)
            if verbose: print self.logN2O2
        except:pass
        try:        
            print  "\nlog([OIII][OII])",stat.nanmean(self.logO3O2)
            if verbose: 
                print self.logO3O2
        except:pass
        for k in self.mds.iterkeys():
            print "\n",k,
            try: print stats.nanmean(self.mds[k]), np.stdev(self.mds[k])
            except: 
                if verbose: print self.mds[k]
    

    def checkminimumreq(self,red_corr,ignoredust, Smass):
        if Smass<MINMASS or Smass > MAXMASS :
            print "will not calculate for this mass: %f MSun"%Smass
            return -1
        if red_corr and not ignoredust:
            if not self.hasHa :
                return -1
            if not self.hasHb :
                return -1
            #if not self.hasO2 :
            #    return -1
            #if not self.hasO3 :
            #    return -1
            #if not self.hasS2 :
            #    return -1
        if not self.hasN2 and not (self.hasO2 and self.hasO3):
            return -1



    def fz_roots(self,coef): 
        if len(coef.shape)==1:
            coef[~(np.isfinite(coef))]=0.0
            rts= np.roots(coef[::-1])
            if rts.size==0:
                print 'fz_roots failed'
                rts=np.zeros(coef.size-1)
            return rts

        else:
            rts=np.zeros((coef.shape[0],coef.shape[1]-1),dtype=complex)
            coef[~(np.isfinite(coef))]=0.0

            for i in range(coef.shape[0]):
                rts[i]= np.roots(coef[i][::-1])#::-1][0])
            return rts

    def setdustcorrect(self):
        global DUSTCORRECT
        DUSTCORRECT=True
    def unsetdustcorrect(self):
        global DUSTCORRECT
        DUSTCORRECT=False

    def dustcorrect(self,l1,l2,flux=False):
        global DUSTCORRECT
        if DUSTCORRECT:
            if not flux:
                return 0.4*self.mds['E(B-V)']*(l1-l2)
            return 10**(0.4*self.mds['E(B-V)']*(l1-l2))        
        else:
            return 1.0

    def setHab(self,Ha,Hb):
        self.Ha=Ha
        self.Hb=Hb
        if sum(self.Ha > 0):
            self.hasHa=True
        if sum(self.Hb > 0):
            self.hasHb=True

    def setOlines(self, O23727, O35007, O16300, O34959):
        self.O23727 = O23727
        self.O35007 = O35007

        if sum(self.O35007>0) : self.hasO3=True
        if sum(self.O23727>0) : self.hasO2=True

        if self.hasO2 and self.hasO3:
            self.logO35007O2=np.log10( self.O35007/self.O23727 )+self.dustcorrect(k_O3,k_O2)
            self.logO2O35007=np.log10( self.O23727/self.O35007 )+self.dustcorrect(k_O2,k_O3)

            #self.logO2O35007Hb=np.log10((self.O23727+self.O35007)/self.Hb)
            # ratios for other diagnostics - slightly different ratios needed
            self.O35007O2=10**(self.logO35007O2)
            self.O2O35007=10**(self.logO2O35007)
            if self.hasHb:
                self.logO2O35007Hb=np.log10((self.O23727/self.Hb)* self.dustcorrect(k_O2,k_Hb,flux=True))+ (self.O35007/self.Hb)*self.dustcorrect(k_O35007,k_Hb,flux=True)

        else: 
            print "WARNING: needs O lines and  and Ha/b: did you run setHab()?"
        if self.hasHb :
            if self.hasO2:
                self.logO2Hb=np.log10(self.O23727/self.Hb)+self.dustcorrect(k_O2,k_Hb)#0.4*self.mds['E(B-V)']*(k_O2-k_Hb) 
            if self.hasO3:
                self.logO3Hb=np.log10(self.O35007/self.Hb)+self.dustcorrect(k_O35007,k_Hb)#0.4*self.mds['E(B-V)']*(k_O2-k_Hb) 
                self.hasO3Hb=True
                if not O34959 == None and sum(O34959>0)>0:
                    self.logO349595007Hb=np.log10(10**(np.log10(self.O35007/self.Hb)+self.dustcorrect(k_O35007,k_Hb))+10**(np.log10(O34959/self.Hb)+self.dustcorrect(k_O34959,k_Hb)))
                    self.O34959p5007=O34959 + self.O35007
                    self.logO3O2=np.log10((self.O34959p5007)/self.O23727)+self.dustcorrect(k_O3,k_O2)
                    self.hasO3O2=True
        # never used
        #if self.hasHa:
            #logO1Ha=np.log10(O16300/self.Ha)+self.dustcorrect(k_O1,k_Ha)
        if self.hasO2 and self.hasO3:
            self.OIII_OII=np.log10(self.O35007/self.O23727+self.dustcorrect(k_O35007,k_O2,flux=True))
        if self.hasHb:
            self.OIII_Hb=np.log10(self.O35007/self.Hb+self.dustcorrect(k_O35007,k_Hb, flux=True))

    def setNII(self,N26584):
        if not N26584==None and sum(N26584>0):
            self.N26584=N26584
            self.hasN2=True
            if self.hasHa :
                self.logN2Ha=np.log10(self.N26584/self.Ha)#+self.dustcorrect(k_N2,k_Ha,flux=True) 
                #lines are very close: no dust correction
                #Note: no dust correction cause the lies are really close!
            else: 
                print "WARNING: needs NII6584 and Ha to calculate NIIHa: did you run setHab()?"
            if self.hasS2 and self.hasS26731 and self.hasN2:
                self.NII_SII=np.log10(self.N26584/(self.S26717+self.S26731))#+self.dustcorrect(k_N2,k_S2,flux=True) 
                #lines are very close: no dust correction
            if self.hasO2 and self.hasN2:
                self.NII_OII=np.log10(self.N26584/self.O23727+self.dustcorrect(k_N2,k_O2,flux=True) )

    def setSII(self,S26717,S26731,S39069,S39532):
        if not S26717==None and sum(S26717>0)>0:
            self.S26717=S26717
            self.hasS2=True

            if self.hasHa:
                self.logS2Ha=np.log10(self.S26717/self.Ha)+self.dustcorrect(k_S2,k_Ha)               
            else: 
                print "WARNING: needs SII6717 and Ha to calculate SIIHa: did you run setHab() and setS()?"
        if not S26731==None and sum(S26731>1e-9)>0:
            self.S26731=S26731
            self.hasS26731=True
        if not S39069==None and sum(S39069>1e-9)>0:
            self.S39069=S39069
            self.hasS39069=True
        if not S39532==None and sum(S39532>1e-9)>0:
            self.S39532=S39532
            self.hasS39532=True
        if self.hasS2 :
            
            if self.hasN2 and self.NII_SII==None and self.hasS26731:
                self.NII_SII=np.log10(self.N26584/(self.S26717+self.S26731))#+self.dustcorrect(k_N2,k_O2,flux=True) 
                    #lines are very close: no dust correction            
            if self.hasO3  and self.OIII_SII==None and self.hasS26731:
                self.OIII_SII=np.log10(self.O35007/(self.S26717+self.S26731)+self.dustcorrect(k_O3,k_S2,flux=True) )


    #@profile
    def calcEB_V(self):
        #logHaHb=np.log10(Ha/Hb)
        self.mds['E(B-V)']=np.log10(2.86/(self.Ha/self.Hb))/(0.4*(k_Ha-k_Hb)) # E(B-V)
        #print self.mds['E(B-V)']
        self.mds['E(B-V)'][self.mds['E(B-V)']<=0]=0.00001
        

    #@profile
    def calcNIISII(self):
        if self.hasS2 and self.hasN2:
            self.N2S2=self.N26584/self.S26717+self.dustcorrect(k_N2,k_S2,flux=True)#0.4*self.mds['E(B-V)']*(k_N2-k_S2) 
            self.logN2S2=np.log10(self.N26584/self.S26717)+self.dustcorrect(k_N2,k_S2)#0.4*self.mds['E(B-V)']*(k_N2-k_S2) 
            self.hasN2S2=True
        else: 
            print "WARNING: needs SII6717 and NII6584 to calculate NIISII: did you run setN2() and setS?"

    #@profile
    def calcNIIOII(self):
        if self.hasN2 and self.hasO2:
            self.logN2O2=np.log10(self.N26584/self.O23727)+self.dustcorrect(k_N2,k_O2) 
            self.hasN2O2=True
            N2O2_coef=np.array([[self.N2O2_coef0,-532.15451,96.373260,-7.8106123,0.23928247]]*self.nm).T# q=2e7 line (approx average)
            N2O2_coef[0]-=self.logN2O2
            N2O2_coef=N2O2_coef.T
            # finding roots for == (4)
            self.N2O2_roots=np.array([self.fz_roots(N2O2_coef)])[0]          

    #@profile
    def calcR23(self):
        #R23 NEW Comb, [NII]/Ha: KK04 = Kobulnicky & Kewley, 2004, submitted'
        if  self.hasO3   and self.hasO2 and self.hasHb:
            self.R23=((self.O23727/self.Hb)*self.dustcorrect(k_O2,k_Hb, flux=True) + (self.O34959p5007/self.Hb)*self.dustcorrect(k_O3,k_Hb, flux=True) )
            self.logR23=np.log10(self.R23)
            #self.R23_5007=(1./self.O35007O2 + 1.)/(1./self.O35007O2 + 1.347)*self.R23  
            self.mds['logR23']=self.logR23
        else:
            print "WARNING: need O3, O2, Hb"
            
    #@profile
    def calcS23(self):
        #the original code here uses S267176731, 
        #which is however set to 6717 as default
        #Vilchez & Esteban (1996)
        if  self.hasS2 :
            if self.hasS39069 and self.hasHb:
                self.logS23=np.log10((self.S26717/self.Hb)*self.dustcorrect(k_S2,k_Hb,flux=True) + (S39069/self.Hb)*self.dustcorrect(k_S3,k_Hb,flux=True))                                 
            self.logS3S2=np.log10(S39069/self.S26717)+self.dustcorrect(k_S3,k_S2)


    def initialguess(self):
        # Initial Guess - appearing in LK code as of Nov 2006
        # upper branch: if no lines are available, metallicity is set to 8.7        
        self.Z_init_guess=np.zeros(self.nm)+8.7 
        # use [N2]/Ha 
        if self.hasHa and self.hasN2:
            self.Z_init_guess[(self.logN2Ha < -1.3)&(self.N26584 != 0.0)]=8.2
            self.Z_init_guess[(self.logN2Ha < -1.1)&(self.logN2Ha >= -1.3)&(self.N26584 != 0.0)]=8.4
            #A1 KE08
            self.Z_init_guess[(self.logN2Ha >=-1.1)&(self.N26584 != 0.0)]=8.7
        #use [N2]/[O2]
        if self.hasN2 and self.hasO2:            
            if self.hasHb:
                N2O2=self.N26584*self.Ha*self.Hb*self.O23727
            if not self.hasN2O2:
                print "WARNING: must calculate logN2O2 first"
                calcNIIOII()
            self.Z_init_guess[(self.logN2O2 < -1.2)&(N2O2 != 0.0)]=8.2  
            # at logN2O2<-1.2 using low-Z gals, A1 KE08
            self.Z_init_guess[(self.logN2O2 >=-1.2)&(N2O2 != 0.0)]=8.7  
            # at logN2O2>-1.2  using HII regions

#######################these are the metallicity diagnostics##################
    #@profile
    def calcpyqz(self, plot=False):
        import pyqz

        if not self.NII_SII==None :
            if not self.OIII_SII ==None:
                self.mds['pyqzN2S2_O3S2']=pyqz.get_qz(20,'z',np.atleast_1d([self.NII_SII]),np.atleast_1d([self.OIII_SII]),'NII/SII','OIII/SII', method='default', plot=plot, n_plot = False, savefig=False )[0].T
            if  not self.OIII_Hb ==None:
                self.mds['pyqzN2S2_O3Hb']=pyqz.get_qz(20,'z',np.atleast_1d([self.NII_SII]),np.atleast_1d([self.OIII_Hb]),'NII/SII','OIII/Hb', method='default', plot=plot, n_plot = False, savefig=False )[0].T
            if  not self.OIII_OII ==None:
                self.mds['pyqzN2S2_O3O2']=pyqz.get_qz(20,'z',np.atleast_1d([self.NII_SII]),np.atleast_1d([self.OIII_OII]),'NII/SII','OIII/OII', method='default', plot=plot, n_plot = False, savefig=False )[0].T

        if not self.NII_OII==None :
            if not self.OIII_SII ==None:
                self.mds['pyqzN2O2_O3S2']=pyqz.get_qz(20,'z',np.atleast_1d([self.NII_OII]),np.atleast_1d([self.OIII_SII]),'NII/OII','OIII/SII', method='default', plot=plot, n_plot = False, savefig=False )[0].T
            if  not self.OIII_Hb ==None:
                self.mds['pyqzN2O2_O3Hb']=pyqz.get_qz(20,'z',np.atleast_1d([self.NII_OII]),np.atleast_1d([self.OIII_Hb]),'NII/OII','OIII/Hb', method='default', plot=plot, n_plot = False, savefig=False )[0].T
            if  not self.OIII_OII ==None:
                self.mds['pyqzN2O2_O3O2']=pyqz.get_qz(20,'z',np.atleast_1d([self.NII_OII]),np.atleast_1d([self.OIII_OII]),'NII/OII','OIII/OII', method='default', plot=plot, n_plot = False, savefig=False )[0].T

        if not self.logN2Ha==None :
            if  not self.OIII_Hb ==None:
                self.mds['pyqzN2Ha_O3Hb']=pyqz.get_qz(20,'z',np.atleast_1d([self.logN2Ha]),np.atleast_1d([self.OIII_Hb]),'NII/Ha','OIII/Hb', method='default', plot=plot, n_plot = False, savefig=False )[0].T
            if  not self.OIII_OII ==None:
                self.mds['pyqzN2Ha_O3O2']=pyqz.get_qz(20,'z',np.atleast_1d([self.logN2Ha]),np.atleast_1d([self.OIII_OII]),'NII/Ha','OIII/OII', method='default', plot=plot, n_plot = False, savefig=False )[0].T


        
    #@profile
    def calcD02(self):
        # [NII]/Ha Denicolo, Terlevich & Terlevich (2002), MNRAS, 330, 69
        #FED:added uncertainties
        e1=np.random.normal(0,0.05,self.nm)
        e2=np.random.normal(0,0.1,self.nm)
        if self.hasN2 and self.hasHa:
            self.mds['D02'] = 9.12+e1+(0.73+e2)*self.logN2Ha
        else:
            print "WARNING: need N2Ha to do this. did you run setHab and setNII"
        
    #@profile
    def calcPP04(self):
        ### PP04_N2_Z, PP04_O3N2_Z Pettini & Pagel diagnostics - Pettini & Pagel (2004), MNRAS, 348, L59
        # [NII]/Ha Pettini & Pagel (2004), MNRAS, 348, L59
        if self.hasN2 and self.hasHa:
            self.mds['PP04_N2']= nppoly.polyval(self.logN2Ha,[9.37, 2.03, 1.26, 0.32])
            #FED: restricting the range as per paper
            index=(self.logN2Ha>-1.5)*(self.logN2Ha<-0.3)
            self.mds['PP04_N2'][~index]=float('NaN')
            #            print self.logN2Ha
            #            print "here", self.mds['PP04_N2']
            if self.hasO3Hb :
                self.mds['PP04_O3N2']=8.73 - 0.32*(self.logO3Hb-self.logN2Ha)
            else:
                print "WARNING: need O3Hb for PP04_O3N2"
        else:
            print "WARNING: need N2Ha to do this. did you run setHab and setNII"


    #@profile
    def calcZ94(self):
        ### calculating z from Kobulnicky,Kennicutt,Pizagno (1998)
        ### parameterization of Zaritzky et al. (1994) 
        ###Z94 = Zaritsky, D., Kennicutt, R. C., & Huchra, J. P. 1994, 
        ###ApJ, 420, 87
        ### only valid on the upper branch of R23 (KE08 A2.4)

        if self.logR23==None:
            print "WARNING: Must first calculate R23"
            self.calcR23()
            if self.logR23==None:
                print "WARNING: Cannot compute this without R23"
        if not self.logR23==None:
            self.mds['Z94']=nppoly.polyval(self.logR23, [9.265,-0.33,-0.202,-0.207,-0.333])
            self.mds['Z94'][(self.logR23 > 0.9)]=None
            ## 0.9 is a conservative constraint to make sure that we are 
            ## only using the upper branch (i.e. 12+log(O/H)>8.4

    def Pi05(self):
        # #### P-method #####
        ##Pilyugin 2001 method.  Based on [OIII],[OII], Hbeta 
        ##calibrated from Te method
        # make sure you run setOlines() first
        if self.logR23==None:
            print "WARNING: Must first calculate R23"
            self.calcR23()
        if self.logR23==None:
            print "WARNING: Cannot compute this without R23"
        else:
            #R3=10**self.logO349595007Hb
            #R2=10**self.logO2Hb
            #P = R3/(R2+R3)
            P=10**self.logO349595007Hb/self.R23
            Psq=P*P
            #P_R23=R2+R3
            #P_R23=self.R23
        
            P_abund_up =(self.R23+726.1+842.2*P+337.5*Psq)/(85.96+82.76*P+43.98*Psq +1.793*self.R23)
            P_abund_low=(self.R23+106.4+106.8*P- 3.40*Psq)/(17.72+ 6.60*P+ 6.95*Psq -0.302*self.R23)

            if self.Z_init_guess==None:
                self.initialguess()
            self.mds['Pi05']=P_abund_up
            self.mds['Pi05'][self.Z_init_guess <  8.4]=P_abund_low[self.Z_init_guess <  8.4]



    #@profile
    def calcPi01_old(self):
        # P-method 2001 upper branch (derecated and commented out)
        # Pilyugin 2001
        # available but deprecated
        if self.Z_init_guess==None:
            self.initialguess()
        if self.hasO3O2 and self.hasO3  and self.hasO2:
            P = 10**self.logO3O2/(1+10**self.logO3O2)
            if self.logR23==None:
                print "WARNING: Must first calculate R23"
                self.calcR23()
                if self.logR23==None:
                    print "WARNING: Cannot compute this without R23"
            if  self.hasHb:
                R3=10**self.logO349595007Hb
                R2=10**self.logO2Hb
                P = R3/(R2+R3)
                P_R23=R2+R3
                P_R23=10**self.logR23
                P_abund_old=(P_R23+54.2+59.45*P+7.31*P**2)/(6.07+6.71*P+0.371*P**2+0.243*P_R23)
            self.mds['Pi01_old']=np.zeros(self.nm)+float('NaN')
            self.mds['Pi01_old'][self.Z_init_guess >= 8.4]=P_abund_old[self.Z_init_guess >= 8.4]
        else:
            print "WARNING: need OIIIOII to calculate Pi01_Z_old, did you set them up with  setOlines()?"
        
    #@profile
    def calcC01_ZR23(self):
        # C01 = Charlot, S., & Longhetti, M., 2001, MNRAS, 323, 887
        # Charlot 01 R23 calibration: (case F) ##        
        # available but deprecated
        if self.hasO3 and self.hasO2 and self.hasO3Hb:
            x2=self.O2O35007/1.5
            x3=(10**self.logO3Hb)*0.5
            self.mds['C01_R23']=np.zeros(self.nm)+float('NaN')        
            self.mds['C01_R23'][self.O2O35007<0.8]=np.log10(3.78e-4 * (x2[self.O2O35007<0.8])**0.17 * x3[self.O2O35007<0.8]**(-0.44))+12.0    
         
            self.mds['C01_R23'][ self.O2O35007 >= 0.8]=np.log10(3.96e-4 * x3[self.O2O35007 >= 0.8]**(-0.46))+12.0   
        else:
            print "WARNING: need [OIII]5700, [OII]3727, and Ha to calculate calcC01_ZR23, did you set them up with  setOlines()?"        

        # Charlot 01 calibration: (case A) based on [N2]/[SII]##
        # available but deprecated
        if not self.hasN2S2:
            print "WARNING: trying to calculate logNIISII"
            self.calcNIISII()
        if self.hasN2S2 and self.hasO3 and self.hasO2 and self.hasO3Hb:
            self.mds['C01']=np.log10(5.09e-4*(x2**0.17)*((self.N2S2/0.85)**1.17))+12
        else:
            print "WARNING: needs [NII]6584, [SII]6717, [OIII]5700, [OII]3727, and Ha to calculate calcC01_ZR23, did you set them up with  setOlines() and ?"        


    #@profile
    def calcM91(self):
        # ## calculating McGaugh (1991)
        # McGaugh, S.S., 1991, ApJ, 380, 140'
        # M91 calibration using [N2O2] as 
        # initial estimate of abundance:
        # this initial estimate can be 
        # changed by replacing 
        # OH_init by another guess, eg C01_Z 
        # NOTE: occasionally the M91 
        # 'upper branch' will give a metallicity
        # that is lower than the 'lower branch'.  
        # Happens for very high R23 values.  
        # If R23 is higher than the intersection 
        # (calculate the intersection), then
        # the metallicity is likely to be around 
        # the R23 maximum = 8.4
        
        if self.logR23==None:
            print "WARNING: Must first calculate R23"
            self.calcR23()
        if self.logR23==None:
            print "WARNING: Cannot compute this without R23"
        else:
            if self.Z_init_guess==None:
                self.initialguess()
                
            self.mds['M91']=np.zeros(self.nm)+float('NaN')
            M91_Z_low=nppoly.polyval(self.logR23,[12.0-4.944,0.767,0.602])-self.logO3O2*nppoly.polyval(self.logR23,[0.29,0.332,-0.331])
            M91_Z_up=nppoly.polyval(self.logR23,[12.0-2.939,-0.2,-0.237,-0.305,-0.0283])-self.logO3O2*nppoly.polyval(self.logR23,[0.0047,-0.0221,-0.102,-0.0817,-0.00717])

            indx=(np.abs(self.logO3O2)>0) * (np.abs(self.logR23)>0) * (self.Z_init_guess < 8.4)
            self.mds['M91'][indx]==M91_Z_low[indx]

            indx=(np.abs(self.logO3O2)>0) * (np.abs(self.logR23)>0) * (self.Z_init_guess >= 8.4)
            self.mds['M91'][indx]=M91_Z_up[indx]

            #2014 FED: changed wrong values to None
            self.mds['M91'][(M91_Z_up < M91_Z_low)]=float('NaN')

    #@profile
    def calcKD02_N2O2(self):
        ##  Kewley & Dopita (2002) estimates of abundance 
        ##  KD02
        # KD02 [N2]/[O2] estimate (can be used for whole log(O/H)+12 range, 
        # but rms scatter increases to 0.11 rms for log(O/H)+12 < 8.6
        # rms = 0.04 for
        # log(O/H)+12 > 8.6
        # uses equation (4) from KD02 paper
        # FED: i vectorized the hell out of this function!!! 
        # from a 7 dimensional if/for loop to 1 if and 1 for :D
        if self.hasN2 and self.hasO2 and self.hasHa and self.hasHb:         
            self.mds['KD02_N2O2']=np.zeros(self.nm)+float('NaN')
            if not self.hasN2O2:
                print "WARNING: must calculate logN2O2 first"
                self.calcNIIOII()

                if  self.N2O2_root == None:
                    print "WARNING:  cannot calculate N2O2"
                    return -1

            roots=self.N2O2_roots.T
            for k in range(4):
                indx=(abs(roots[k]) >= 7.5) * (abs(roots[k]) <= 9.4) * (roots[k][:].imag ==  0.0 )
                self.mds['KD02_N2O2'][indx]=abs(roots[k][indx]) 
        else:
            print "WARNING: need NII6584 and OII3727 and Ha and Hb to calculate this. did you run setO() setHab() and setNII()?"
        return 1



    #@profile
    def calcKD02_N2Ha(self):
        # calculating [N2]/Ha abundance estimates using [O3]/[O2] also
        if self.mds['KD02_N2O2'] == None:
            self.calcKD02_N2O2()
            if self.mds['KD02_N2O2'] == None:
                print "WARNING: without KD02_N2O2 cannot calculate KD02_N2Ha"
                return -1

        Z_new_N2Ha=self.mds['KD02_N2O2'].copy()  # was 8.6

        if self.hasN2 and self.hasHa:
            self.logq_save=np.zeros(self.nm)
            convergence,ii=100,0
            while convergence>1e-3 and ii<10:
                ii+=1
                if self.hasO3O2 :        
                    # calculating logq using the [N2]/[O2] 
                    #metallicities for comparison
                    logO3O2sq=self.logO3O2**2 
                    self.logq=(32.81 -1.153*logO3O2sq + Z_new_N2Ha*(-3.396 -0.025*self.logO3O2 + 0.1444*logO3O2sq))/(4.603-0.3119*self.logO3O2 -0.163*logO3O2sq+ Z_new_N2Ha*(-0.48 + 0.0271*self.logO3O2+ 0.02037*logO3O2sq)) 
                else:        
                    self.logq=7.37177
                Z_new_N2Ha=nppoly.polyval(self.logN2Ha,[7.04, 5.28,6.28,2.37])-self.logq*nppoly.polyval(self.logN2Ha,[-2.44,-2.01,-0.325,+0.128])+10**(self.logN2Ha-0.2)*self.logq*(-3.16+4.65*self.logN2Ha)
                convergence=np.abs(self.logq-self.logq_save).mean()
                self.logq_save=self.logq.copy()
            self.mds['KD02_N2Ha']=Z_new_N2Ha
        else:
            print "WARNING: need NII6584  and Ha to calculate this. did you run  setHab() and setNII()?"




    #@profile
    def calcKD02R23(self):
        #Kewley, L. J., & Dopita, M. A., 2003 
        # calculating upper and lower metallicities for objects without
        # Hb  and for objects without [O3] and/or [O2]
        Hb_up_ID=np.zeros(100)
        if self.hasN2 and self.hasHa:
            logq_lims=[6.9,8.38]
            #logN2Ha=np.log10(self.N26584/self.Ha) CHECK!! why remove dust correction??

            Z_new_N2Ha_lims= np.atleast_2d([1.0,1.0]).T*nppoly.polyval(self.logN2Ha,[7.04, 5.28,6.28,2.37])-np.atleast_2d( logq_lims).T*nppoly.polyval(self.logN2Ha,[-2.44,-2.01,-0.325,0.128])+np.atleast_2d(logq_lims).T*(10**(self.logN2Ha-0.2)*(-3.16+4.65*self.logN2Ha))

            # #### New ionization parameter and metallicity diagnostics #######
            # NEW R23 diagnostics from Kobulnicky & Kewley 
            # See NEW_fitfin.coefs

            # trying out just using a single [NII]/Ha value
        
        Zmax=np.zeros(self.nm)

        # ionization parameter
        if not self.hasO3O2:
            logq=np.zeros(self.nm)
        else:
            if self.Z_init_guess==None:
                self.initialguess()
            Z_new=self.Z_init_guess.copy()
            if self.logR23==None:
                print "WARNING: Must first calculate R23"
                self.calcR23()
            if self.logR23==None:
                print "WARNING: Cannot compute this without R23" 
            else:
                for ii in range(4):
                    #3 iterations are typically enought to achieve convergence KE08 A2.3
                    logq=(32.81 -1.153*self.logO3O2**2 + Z_new*(-3.396 - 0.025*self.logO3O2 + 0.1444*self.logO3O2**2))/(4.603 - 0.3119*self.logO3O2 - 0.163*self.logO3O2**2+Z_new*(-0.48 + 0.0271*self.logO3O2 + 0.02037*self.logO3O2**2))
                    Zmax[(logq >= 6.7) * (logq < 8.3)]=8.4
                    # maximum of R23 curve:               
                    Z_new=nppoly.polyval(self.logR23,[9.72, -0.777,-0.951,-0.072,-0.811])-logq*nppoly.polyval(self.logR23,[0.0737,  -0.0713, -0.141, 0.0373, -0.058])
                    Z_new_lims=[nppoly.polyval(self.logR23,[9.40, 4.65,-3.17])-logq*nppoly.polyval(self.logR23,[0.272,0.547,-0.513]),
                                nppoly.polyval(self.logR23,[9.72, -0.777,-0.951,-0.072,-0.811])-logq*nppoly.polyval(self.logR23,[0.0737,  -0.0713, -0.141, 0.0373, -0.058])]

                    indx=self.Z_init_guess<=Zmax
                    Z_new[indx]=nppoly.polyval(self.logR23[indx], [9.40 ,4.65,-3.17])-logq[indx]*nppoly.polyval(self.logR23[indx],[0.272,+0.547,-0.513])

                #if (self.hasHb):
                    #2014 FED: changed moc value to None, not 0! and removed the if Hb
                Z_new[(Z_new_lims[0]>Z_new_lims[1])]=None
                self.mds['KD02_R23']=Z_new



    #@profile
    def calcKDcombined(self):
        #KD02comb_updated  Kewley, L. J., & Dopita, M. A., 2002, ApJ, submitted '

        # ### KD02 [NII]/[OII] estimate ###
        # (can be used for for log(O/H)+12 > 8.6 only)
        # uses equation (5) from paper, this should be identical 
        # to the estimate above for abundances log(O/H)+12 > 8.6
        if not self.hasN2O2:
            self.calcNIIOII()

        #alternative way to calculate KD02_N2O2, but we forego it for now
        #if not self.logN2O2==None:
        #    self.mds['KD02_N2O2']=np.log10(8.511e-4*(1.54020+1.26602*self.logN2O2+0.167977*self.logN2O2**2))+12.
        #else: self.mds['KD02_N2O2']=np.zeros(self.nm)+float('NaN')
        if self.mds['KD02_N2Ha']==None:
            self.calcKD02_N2Ha()
        # ionization parameter        
        logq_final=np.zeros(self.nm)
        if self.hasN2 and self.hasO2 and self.hasHb and self.hasHa and self.hasO3O2:
            logq_final=(32.81 + 0.0*self.logO3O2-1.153*self.logO3O2**2 +self.mds['KD02_N2O2']*(-3.396 -0.025*self.logO3O2 + 0.1444*self.logO3O2**2))/(4.603  -0.3119*self.logO3O2 -0.163*self.logO3O2**2+self.mds['KD02_N2O2']*(-0.48 +0.0271*self.logO3O2+ 0.02037*self.logO3O2**2))

            logq_final[self.mds['KD02_N2O2']<=8.4]=self.logq[self.mds['KD02_N2O2']<=8.4]



        if  not self.hasN2 and self.hasO2 and self.hasO3 and self.hasHb and self.hasHa:
            logq_final=self.logq

        # if [NII]/[OII] after extinction correction is less than -1.5, then check the data.
        # if it is only slightly less than 1.5, then this can be a result of either noisy
        # data, inaccurate fluxes or extinction correction, or a higher ionization parameter
        # than modelled.  For these cases, the average of the M91,Z94 and C01 should be used.
        
        # KD02 R23 estimate (not reliable for  8.4 < log(O/H)+12 < 8.8)
        # uses [NII]/[OII] estimate as initial guess - this can be changed below
        
        # initializing:
        Zi=np.array([0.05,0.1,0.2,0.5,1.0,1.5,2.0,3.0])      # model grid abundances in solar units
        ZiOH=np.log10(Zi*8.511e-4)+12    # log(O/H)+12 units
        Zstep=np.array([0.025,0.075,0.15,0.35,0.75,1.25,1.75,2.5,3.5]) #middle of model grid abundances
        ZstepOH=np.log10(Zstep*8.511e-4)+12
        qstep=np.log10([3.5e6,7.5e6,1.5e7,3e7,6e7,1.16e8,2.25e8,3.25e8]) # model grid ionization parameters
        n_ite=3                          # number of iteations for abundance determination
        tol=1.0e-2                       # tolerance for convergance 
        R23_roots=np.zeros((4,self.nm),dtype=complex)   # all possible roots of R23 diagnostic
        q_roots=np.zeros((3,self.nm),dtype=complex)     # possible roots of q diagnostic
        q=np.zeros((self.nm,n_ite+1))        # actual q value
        O3O2_coef=np.zeros((4,8))        # coefficients from model grid fits
        R23_coef=np.zeros((5,7))         # coefficients from model grid fits
        R23_Z=np.zeros((self.nm,n_ite+1))    # Z value for each iteation
        if not self.mds['KD02_N2O2'] ==None:  R23_Z[:,0]=self.mds['KD02_N2O2'].copy()  # use NIIOII abundance as initial estimate

        # occasionally, for noisy data or badly fluxed [OII].[OIII] or Hb, 
        # or for high ionization parameter galaxies, R23 is slightly higher
        # than the curves in our models - this will result in all complex roots of
        # the R23 curve unless a slightly lower R23 is used.  These should
        # be checked individually to make sure that it is just noise etc in the
        # data that is causing the problem, rather than wrong fluxes input.
        # the R23 ratio should be close to 0.95, not much more than 1.0 for the
        # data to be physical.


        if self.logR23 == None:
            self.calcR23()
        if not self.hasO3 or not self.hasO2 or self.logR23==None:            
            print "WARNING:  Must first calculate R23 and O350072" 
            
        else:
            R23c0=[-3267,-3727.42,-4282.30,-4745.18,-4516.46,-3509.63,-1550.53]
            R23_coef[:,0]=[-3267.93,1611.04,-298.187,24.5508,-0.758310] # q=5e6
            R23_coef[:,1]=[-3727.42,1827.45,-336.340,27.5367,-0.845876] # q=1e7
            R23_coef[:,2]=[-4282.30,2090.55,-383.039,31.2159,-0.954473] # q=2e7
            R23_coef[:,3]=[-4745.18,2309.42,-421.778,34.2598,-1.04411]  # q=4e7
            R23_coef[:,4]=[-4516.46,2199.09,-401.868,32.6686,-0.996645] # q=8e7
            R23_coef[:,5]=[-3509.63,1718.64,-316.057,25.8717,-0.795242] # q=1.5e8
            R23_coef[:,6]=[-1550.53,784.262,-149.245,12.6618,-0.403774] # q=3e8

            O3O2c0=[-36.9772,-74.2814,-36.7948,-81.1880,-52.6367,-86.8674,-24.4044,49.4728]
            O3O2_coef[:,0]=[-36.9772,10.2838,-0.957421,0.0328614] #z=0.05 
            O3O2_coef[:,1]=[-74.2814,24.6206,-2.79194,0.110773]    # z=0.1
            O3O2_coef[:,2]=[-36.7948,10.0581,-0.914212,0.0300472]  # z=0.2
            O3O2_coef[:,3]=[-81.1880,27.5082,-3.19126,0.128252]    # z=0.5
            O3O2_coef[:,4]=[-52.6367,16.0880,-1.67443,0.0608004]   # z=1.0
            O3O2_coef[:,5]=[-86.8674,28.0455,-3.01747,0.108311]    # z=1.5
            O3O2_coef[:,6]=[-24.4044,2.51913,0.452486,-0.0491711]  # z=2.0
            O3O2_coef[:,7]=[49.4728,-27.4711,4.50304,-0.232228]    # z=3.0
            
            R23_coefi =np.zeros((self.nm,5,7))+R23_coef
            O3O2_coefi=np.zeros((self.nm,4,8))+O3O2_coef

            # coefficients from KD02 paper:
            R23_coefi[:,0,:]= (np.ones((self.nm,1))*R23c0)- (np.ones((7,1))*self.logR23).T
            O3O2_coefi[:,0,:]=(np.ones((self.nm,1))*O3O2c0)-(np.ones((8,1))*self.logO35007O2).T
            # coefficients from KD02 paper:            
            for ite in range(1,n_ite+1) : 
                # iteate if tolerance level not met
                indx=( abs(R23_Z[:,ite]-R23_Z[:,ite-1]) > tol)
                        #   calculate ionization parameter using [OIII]/[OII] with
                        #   [NII]/[OII] abundance for the first iteation, and the R23
                        #   abundance in consecutive iteations
                for j in range(8):  
                    indx1=( R23_Z[:,ite-1] > ZstepOH[j] )*(R23_Z[:,ite-1] <= ZstepOH[j+1] )*indx
                    #                                indx2= R23_Z[i,ite-1] <= ZstepOH[j+1] :
                    if not sum(indx1): continue
                    q_roots[:,indx1]=self.fz_roots(O3O2_coefi[indx1,:,j]).T  

                    for i,ii in enumerate(indx1) :
                        if not ii:
                            continue
                        #q must be between 3.5e6 and 3.5e8 cm/s 
                        #because that is the range it
                        #is defined over by the model grids, and it must be real.

                        for k in range(3) :
                            if (q_roots[k,i].imag) == 0.0 :
                                if (q_roots[k,i].real) >= 6.54407 :   #log units (q >= 1e6 cm/s) 
                                    if (q_roots[k,i].real) <= 8.30103 :   #log units (q <= 2e8 cm/s)
                                        q[i,ite]=(q_roots[k,i].real)
                                        
                        #   calculate abundance using ionization parameter:
                        R23_qstepno=0
                        
                        for j in range(7) :   
                            if q[i,ite] > qstep[j] :
                                if q[i,ite] <= qstep[j+1] :
                                    R23_roots[:,i]=self.fz_roots(R23_coefi[i,:,j])
                                    R23_qstepno=j
     
                                    #   There will be four roots, two complex ones, 
                                    #   and two real ones.
                                    #   use previous R23 value (or [NII]/[OII] if first iteation) 
                                    #   and q to find which real root to use 
                                    #   (ie. which side of R23 curve to use).  
                                    
                        #    Rmax=[1.04967,1.06497,1.06684,1.06329,1.03844,0.991261,0.91655]
                        Smax=np.array([8.69020,8.65819,8.61317,8.58916,8.49012,8.44109,8.35907])
                            
                        for k in range(4) :
                            if (R23_roots[k,i].imag) == 0.0 :
                                if (R23_Z[i,ite-1] >= Smax[R23_qstepno] and R23_roots[k,i].real >= Smax[R23_qstepno]) or (R23_Z[i,ite-1] <= Smax[R23_qstepno] and R23_roots[k,i].real <= Smax[R23_qstepno]):
                                    R23_Z[i,ite]=R23_roots[k,i].real
                                   
                                    # around maximum of R23 sometimes the R23 ratio 
                                    # will be slightly higher than
                                    # that available for the closest q value.  
                                    # This will depend on noise addded to data.  
                                    # If this happens, step up in ionization parameter 
                                    # to find abundance using that one instead. 
                                    # Around local maximum, the actual ionization parameter
                                    # used is not significant compared to the errors 
                                    # associated with the lack of
                                    # abundance sensitivity of the R23 ratio in this region.

                        while R23_Z[i,ite] == 0.0 and R23_qstepno <= 5 :
                            R23_roots[:,i]=self.fz_roots(R23_coef[:,R23_qstepno+1])
                            for k in range(4):
                                if (R23_roots[k,i].imag) == 0.0 :
                                    if (R23_Z[i,ite-1] >= Smax[R23_qstepno] and R23_roots[k,i].real >= Smax[R23_qstepno]) or (R23_Z[i,ite-1] <= Smax[R23_qstepno] and R23_roots[k,i].real <= Smax[R23_qstepno]):
                                        R23_Z[i,ite]=R23_roots[k,i].real
                            R23_qstepno+=1
                            
                            
                R23_Z[indx*(-1),ite]=R23_Z[indx*(-1),ite-1]  
                q[indx*(-1),ite]=q[indx*(-1),ite-1]


        KD02_R23_Z=R23_Z[:,n_ite]

        #  ### Combined \R23\ method outlined in KD02 paper Section 7. ###
        #  ie for objects with only [OII], [OIII], Hb available
        if not self.hasHa and not self.hasHb:
            print "WARNING: need Halpha and Hbeta for this. did you run setHab()?"
        # KD01 combined method (uses [NII], [OII], [OIII], [SII]):
        # uses KD02 [NII]/[OII] method if [NII]/[OII] gives 8.6 < log(O/H)+12
        # uses average of M91 and Z94 if 8.5 < log(O/H)+12 < 8.6
        # uses average of C01 and KD02 if  log(O/H)+12 < 8.5
        # Also calculates comparison average

        KD02C01_ave=np.zeros(self.nm)
        M91Z94C01_ave=np.zeros(self.nm)
        M91Z94_ave=np.zeros(self.nm)
        
        if self.mds['M91'] == None:
            print "WARNING:  Must first calculate M91"
            self.calcM91()
        if self.mds['Z94'] == None:
            print "WARNING:  Must first calculate Z94"
            self.calcZ94()

        if KD02_R23_Z == None or self.mds['Z94']==None or self.mds['M91']==None:
            print "WARNING:  cannot calculate KD02comb_R23 because  KD02_R23, M91, or Z94 failed"
        else:
            self.mds['KD02comb_R23']=np.zeros(self.nm)+float('NaN')

            # LK02 averaged with M91 and Z94            
            indx=self.mds['Z94']>=9.0
            self.mds['KD02comb_R23'][indx]=(KD02_R23_Z[indx]+self.mds['M91'][indx]+self.mds['Z94'][indx])/3.  
            
            # average of M91 and Z94
            indx= (self.mds['KD02comb_R23'] <= 9.0) * (self.mds['KD02comb_R23'] >= 8.5)
            self.mds['KD02comb_R23'][indx]=0.5*(self.mds['M91'][indx]+self.mds['Z94'][indx])                  
        
            # average of M91 and Z94
            indx=(self.mds['Z94'] <= 9.0) * (self.mds['Z94'] >= 8.5)
            self.mds['KD02comb_R23'][indx]=0.5*(self.mds['M91'][indx]+self.mds['Z94'][indx])                 
            
            indx= self.mds['KD02comb_R23'] <= 8.5 
            self.mds['KD02comb_R23'][indx]=KD02_R23_Z[indx]                        
            
            indx= self.mds['Z94'] <= 8.5 
            self.mds['KD02comb_R23'][indx]=KD02_R23_Z[indx]                        

            #KD01 combined
            indx=(np.abs(self.mds['M91'])>0) * (np.abs(self.mds['Z94'])>0)
            M91Z94_ave[indx]=0.5*(self.mds['M91'][indx]+self.mds['Z94'][indx])
            
            #indx =(np.abs(self.mds['C01'])>0) *( np.abs(self.mds['M91'])>0) * (np.abs(self.mds['Z94'])>0)
            #M91Z94C01_ave[indx]=(self.mds['M91'][indx]+self.mds['Z94'][indx]+self.mds['C01'][indx])/3.
            if not  self.mds['C01']==None:            
                indx=(np.abs(KD02_R23_Z)> 0.0) * (np.abs(self.mds['C01'])>0)
                KD02C01_ave[indx]=0.5*(KD02_R23_Z[indx]+self.mds['C01'][indx])

                
        

        # ### [NII]/[SII] method outlined in KD02 paper ###
        # this method produces a systematic shift of 0.2 dex in log(O/H)+12
        # compared with the average of M91, Z94, and C01.  We believe this
        # is a result of inaccurate abundances or depletion factors, which are known 
        # problems in sulfur modelling.  Initial guess of [NII]/[OII] used
        # can be changed.  Initial guess is not critical except for high
        # ionization parameters.  ionization parameter diagnostic is [OIII]/[OII]
        
        KD02_N2S2_Z=np.zeros(self.nm)
        
        
        N2S2_roots=np.zeros((4,self.nm),dtype=complex)   # all possible roots of NIISII diagnostic
        q_roots=np.zeros((3,self.nm),dtype=complex)     # possible roots of q diagnostic
        q=np.zeros((self.nm,n_ite+1))        # actual q value
        N2S2_coef=np.zeros((5,7))          # coefficients from model grid fits
        N2S2_Z=np.zeros((self.nm,n_ite+1))    # Z value for each iteation
        
        # initializing:
        
        N2S2_Z[:,0]=self.mds['KD02_N2O2'].copy()  # use [NII]/[OII] abundance as initial estimate

        if self.hasO3 and self.hasO2 and self.hasN2S2:
            # coefficients from KD02 paper:
            N2S2c0=[-1042.47,-1879.46,-2027.82,-2080.31,-2162.93,-2368.56,-2910.63]
            N2S2_coef[:,0]=[-1042.47,521.076,-97.1578,8.00058,-0.245356]
            N2S2_coef[:,1]=[-1879.46,918.362,-167.764,13.5700,-0.409872]
            N2S2_coef[:,2]=[-2027.82,988.218,-180.097,14.5377,-0.438345]
            N2S2_coef[:,3]=[-2080.31,1012.26,-184.215,14.8502,-0.447182]
            N2S2_coef[:,4]=[-2162.93,1048.97,-190.260,15.2859,-0.458717]
            N2S2_coef[:,5]=[-2368.56,1141.97,-205.908,16.4451,-0.490553]
            N2S2_coef[:,6]=[-2910.63,1392.18,-249.012,19.7280,-0.583763]


            N2S2_coefi =np.zeros((self.nm,5,7))+R23_coef
            N2S2_coefi[:,0,:]=(np.ones((self.nm,1))*N2S2c0)-(np.ones((7,1))*self.logN2S2).T

            for ite in range(1, n_ite+1):                 # iteate if tolerance level not met
                indx = (abs(N2S2_Z[:,ite]-N2S2_Z[:,ite-1]) >= tol )
                        #   calculate ionization parameter using [OIII]/[OII] with
                        #   [NII]/[OII] abundance for the first iteation, and the [NII]/[SII]
                        #   abundance in consecutive iteations
                for j in range(8) :   
                    indx1=(N2S2_Z[:,ite-1]>ZstepOH[j] )*(N2S2_Z[:,ite-1]<=ZstepOH[j+1] )*indx
                        #   q must be between 3.5e6 and 3.5e8 cm/s because that is 
                        #   the range it
                        #   is defined over by the model grids, and it must be real.
                    if not sum(indx1): continue
                    q_roots[:,indx1]=self.fz_roots(O3O2_coefi[indx1,:,j]).T  
#                    print j,q_roots
                    for i,ii in enumerate(indx1) :
                        if not ii:
                            continue
 
                        for k in range(3) :
                            if (q_roots[k,i].imag) == 0.0 :
                                if (q_roots[k,i].real) >= 6.54407 :   #log units (q >= 1e6 cm/s) 
                                    if (q_roots[k,i].real) <= 8.30103 :   #log units (q <= 2e8 cm/s)
                                        q[i,ite]=(q_roots[k,i].real)

                        #   calculate abundance using ionization parameter:
                        N2S2_qstepno=0
                        for j in range(7) :   
                            if q[i,ite] > qstep[j] and q[i,ite] <= qstep[j+1] :
                                    N2S2_roots[:,i]=self.fz_roots(N2S2_coefi[i,:,j])
                                    N2S2_qstepno=j
                   
                                    #   There will be four roots, two complex ones, 
                                    #   and two real ones.
                                    #   use previous NIISII value 
                                    #   (or [NII]/[OII] if first iteation) 
                                    #   and q to find which real root to use 
                                    #   (ie. which side of R23 curve 
                                    #   to use).  

                        for k in range(4) :
                            if (N2S2_roots[k,i].imag) == 0.0 and (N2S2_roots[k,i].real) >= 8.0 and (N2S2_roots[k,i].real) <= 9.35 :
                                    N2S2_Z[i,ite]=(N2S2_roots[k,i].real)


                        if N2S2_Z[i,ite] == 0.0 :
                            N2S2_roots[:,i]=self.fz_roots(N2S2_coef[:, N2S2_qstepno+1])
                            for k in range(4) :
                                if (N2S2_roots[k,i].imag) == 0.0 and ((N2S2_roots[k,i].real) >= 8.0) and ((N2S2_roots[k,i].real) <= 9.35) :
                                        N2S2_Z[i,ite]=(N2S2_roots[k,i].real)

                    
                N2S2_Z[indx*(-1),ite]=N2S2_Z[indx*(-1),ite-1]  
                q[indx*(-1),ite]=q[indx*(-1),ite-1]
        
        KD02_N2S2_Z=N2S2_Z[:,n_ite]
        if not self.mds['KD02_N2O2']==None:
            self.mds['KD02comb']=self.mds['KD02_N2O2'].copy()
        else:
            self.mds['KD02comb']=np.zeros(self.nm)+float('NaN')
        if KD02_R23_Z == None or self.mds['Z94']==None or self.mds['M91']==None:
            print "WARNING:  cannot calculate KD02comb_R23 because  KD02_R23, M91, or Z94 failed"
        else:
            indx= (self.mds['KD02_N2O2'] <= 8.6 ) * (M91Z94_ave >= 8.5 )
            self.mds['KD02comb'][indx]=M91Z94_ave[indx]   # average of M91 and Z94
            indx= (self.mds['KD02_N2O2'] <= 8.6 ) * (M91Z94_ave < 8.5 )
            self.mds['KD02comb'][indx]=KD02C01_ave[indx]
        
            #-----------------------------------------
            # ### combined method ###
            #-----------------------------------------
            
            # if [NII]/[OII] abundance available and [NII]/Ha abundance < 8.4, then 
            # use R23. 
            
            self.mds['KD02comb_updated']=np.zeros(self.nm)+float('NaN')
            indx=self.Z_init_guess > 8.4
            #self.mds['KD02_R23']=KD02_R23_Z#np.zeros(self.nm)+float('NaN')
            self.mds['KD02comb_updated'][indx]=self.mds['KD02_N2O2'][indx].copy()

            indx=(self.mds['KD02_R23'] > 0.0) * (self.mds['M91'] > 0.0 ) * (self.Z_init_guess <= 8.4)
            self.mds['KD02comb_updated'][indx]=0.5*(self.mds['KD02_R23'][indx].copy()+self.mds['M91'][indx].copy())
            indx=(self.mds['KD02_R23'] <= 0.0) * (self.mds['M91'] <= 0.0 ) * (self.Z_init_guess <= 8.4)
            if not self.mds['KD02_N2Ha']==None:
                self.mds['KD02comb_updated'][indx]=self.mds['KD02_N2Ha'][indx].copy()
            
