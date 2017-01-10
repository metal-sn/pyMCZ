from __future__ import print_function
import numpy as np
#import sys
import scipy.stats as  stats
import numpy.polynomial.polynomial as nppoly
from metallicity import get_keys, printsafemulti

niter = 5  # number of iteations+1 for KD02 methods

k_Ha = 2.535  # CCM Rv=3.1
k_Hb = 3.609  # CCM Rv=3.1

#k_O1=2.661   # CCM Rv=3.1
k_O2 = 4.771     # CCM Rv=3.1
k_O35007 = 3.341  # CCM Rv=3.1
k_O34959 = 3.384  # CCM Rv=3.1
k_O3 = (k_O35007 + k_O34959) / 2.

k_N2 = 2.443  # CCM Rv=3.1
k_S2 = 2.381  # CCM Rv=3.1

k_S3 = 1  # guess for CCM Rv=3.1

global DUSTCORRECT
DUSTCORRECT = True

'''
R23_coef=np.zeros((5,7))         # coefficients from model grid fits
R23c0=[-3267,-3727.42,-4282.30,-4745.18,-4516.46,-3509.63,-1550.53]
R23_coef[:,0]=[-3267.93,1611.04,-298.187,24.5508,-0.758310] # q=5e6
R23_coef[:,1]=[-3727.42,1827.45,-336.340,27.5367,-0.845876] # q=1e7
R23_coef[:,2]=[-4282.30,2090.55,-383.039,31.2159,-0.954473] # q=2e7
R23_coef[:,3]=[-4745.18,2309.42,-421.778,34.2598,-1.04411]  # q=4e7
R23_coef[:,4]=[-4516.46,2199.09,-401.868,32.6686,-0.996645] # q=8e7
R23_coef[:,5]=[-3509.63,1718.64,-316.057,25.8717,-0.795242] # q=1.5e8
R23_coef[:,6]=[-1550.53,784.262,-149.245,12.6618,-0.403774] # q=3e8


N2S2_coef=np.zeros((5,7))          # coefficients from model grid fits
N2S2c0=[-1042.47,-1879.46,-2027.82,-2080.31,-2162.93,-2368.56,-2910.63]
N2S2_coef[:,0]=[-1042.47,521.076,-97.1578,8.00058,-0.245356]
N2S2_coef[:,1]=[-1879.46,918.362,-167.764,13.5700,-0.409872]
N2S2_coef[:,2]=[-2027.82,988.218,-180.097,14.5377,-0.438345]
N2S2_coef[:,3]=[-2080.31,1012.26,-184.215,14.8502,-0.447182]
N2S2_coef[:,4]=[-2162.93,1048.97,-190.260,15.2859,-0.458717]
N2S2_coef[:,5]=[-2368.56,1141.97,-205.908,16.4451,-0.490553]
N2S2_coef[:,6]=[-2910.63,1392.18,-249.012,19.7280,-0.583763]

O3O2_coef=np.zeros((4,8))        # coefficients from model grid fits
O3O2c0=[-36.9772,-74.2814,-36.7948,-81.1880,-52.6367,-86.8674,-24.4044,49.4728]
O3O2_coef[:,0]=[-36.9772,10.2838,-0.957421,0.0328614] #z=0.05 
O3O2_coef[:,1]=[-74.2814,24.6206,-2.79194,0.110773]    # z=0.1
O3O2_coef[:,2]=[-36.7948,10.0581,-0.914212,0.0300472]  # z=0.2
O3O2_coef[:,3]=[-81.1880,27.5082,-3.19126,0.128252]    # z=0.5
O3O2_coef[:,4]=[-52.6367,16.0880,-1.67443,0.0608004]   # z=1.0
O3O2_coef[:,5]=[-86.8674,28.0455,-3.01747,0.108311]    # z=1.5
O3O2_coef[:,6]=[-24.4044,2.51913,0.452486,-0.0491711]  # z=2.0
O3O2_coef[:,7]=[49.4728,-27.4711,4.50304,-0.232228]    # z=3.0
'''

M08_coefs = {'R23': [0.7462, -0.7149, -0.9401, -0.6154, -0.2524],
           'N2Ha': [-0.7732, 1.2357, -0.2811, -0.7201, -0.3330],
           'O3Hb': [0.1549, -1.5031, -0.9790, -0.0297],
           'O3O2': [-0.2839, -1.3881, -0.3172],
           'O2Hb': [0.5603, 0.0450, -1.8017, -1.8434, -0.6549],
           'O3N2': [0.4520, -2.6096, -0.7170, 0.1347]}

#this is to check the Maiolino coefficients and find the split maximum of the cirves with degeneracy
'''
import pylab as pl

x=np.arange(7.0,9.5,0.1)

for k in M08_coefs.iterkeys():
    print k,max(nppoly.polyval(x-8.69,M08_coefs[k])),x[nppoly.polyval(x-8.69,M08_coefs[k])==max(nppoly.polyval(x-8.69,M08_coefs[k]))]
    pl.plot(x,nppoly.polyval(x-8.69,M08_coefs[k]), label=k+' max:%.1f'%x[nppoly.polyval(x-8.69,M08_coefs[k])==max(nppoly.polyval(x-8.69,M08_coefs[k]))])

print nppoly.polyval(8.4-8.69,M08_coefs['N2Ha'])
pl.ylabel("log R")
pl.xlabel("12+log(O/H)")
pl.legend(loc=3)
pl.show()
'''


class diagnostics:
    def __init__(self, num, logf, nps):
        self.nm = num
        self.Ha = None
        self.Hb = None

        self.hasHa, self.hasHb = False, False
        self.hasO2, self.hasO3 = False, False
        self.hasS2, self.hasN2 = False, False

        self.hasO3Hb = False
        self.hasO3O2 = False

        self.hasN2O2 = False
        self.hasN2S2 = False

        self.hasS26731 = False
        self.hasS39532 = False
        self.hasS39069 = False
        self.hasS2Hb = False

        self.N2O2_roots = None
        #other lines calculated and repeatedly used
        self.P = None
        self.R2 = None
        self.R3 = None
        self.R23 = None
        self.S2Hb = None
        self.N2 = None
        self.N2S2 = None
        self.O23727 = None
        self.O35007 = None
        self.N26584 = None
        self.S26717 = None
        self.S26731 = None
        self.S39069 = None
        self.S39532 = None

        self.O34959p5007 = None
        self.O35007O2 = None
        self.O2O35007 = None

        self.logR23 = None
        self.logN2O2 = None
        self.logN2S2 = None
        self.logO3O2 = None
        self.logS23 = None
        #self.logS3S2=None

        self.logO2Hb = None
        self.logO3Hb = None
        self.logN2Ha = None
        self.logS2Ha = None

        self.logO35007O2 = None
        self.logO2O35007 = None

        self.logO3O2sq = None
        self.logq = None
        self.Z_init_guess = None
        self.N2O2_coef0 = 1106.8660

        self.OIII_OII = None
        self.OIII_Hb = None
        self.OIII_SII = None

        self.NII_OII = None
        self.NII_SII = None
        #metallicity diagnostics to be returned

        self.mds = {}
        for Z in  get_keys():
            self.mds[Z] = None

        #setting output file
        self.logf = logf
        self.nps = nps

    def printme(self, verbose=False):
        try:
            print ("\nHa", np.mean(self.Ha))
            if verbose:
                print (self.Ha)
        except (IndexError, TypeError):
            pass
        try:
            print ("\nHb", np.mean(self.Hb))
            if verbose:
                print (self.Hb)
        except (IndexError, TypeError):
            pass
        try:
            print ("\nO2", np.mean(self.O23727))
            if verbose:
                print (self.O23727)
        except (IndexError, TypeError):
            pass
        try:
            print ("\nO3", np.mean(self.O35007))
            if verbose:
                print (self.O35007)
        except (IndexError, TypeError):
            pass
        try:
            print ("\nO34959", np.mean(self.O34959))
            if verbose:
                print (self.O34959)
        except (IndexError, TypeError):
            pass
        try:
            print ("\nZ94", np.mean(self.mds['Z94']))
            if verbose:
                print (self.mds['Z94'])
        except (IndexError, TypeError):
            pass
        try:
            print ("\nR23", np.mean(self.R23))
            if verbose:
                print (self.R23)
        except (IndexError, TypeError):
            pass
        try:
            print ("\nlog(R23)", np.mean(self.logR23))
            if verbose:
                print (self.logR23)
        except (TypeError, IndexError):
            pass
        try:
            print ("\nlog([NII][OII])", stats.nanmean(self.logN2O2))
            if verbose:
                print (self.logN2O2)
        except (TypeError, IndexError):
            pass
        try:
            print ("\nlog([OIII][OII])", stats.nanmean(self.logO3O2))
            if verbose:
                print (self.logO3O2)
        except (TypeError, IndexError):
            pass
        for k in self.mds.iterkeys():
            print ("\n", k)
            try:
                print (stats.nanmean(self.mds[k]), np.stdev(self.mds[k]))
            except (IndexError, TypeError):
                if verbose:
                    print (self.mds[k])
    
    def checkminimumreq(self, red_corr, ignoredust):
        if red_corr and not ignoredust:
            if not self.hasHa:
                return -1
            if not self.hasHb:
                return -1
            #if not self.hasO2 :
            #    return -1
            #if not self.hasO3 :
            #    return -1
            #if not self.hasS2 :
            #    return -1
        if not self.hasN2 and not (self.hasO2 and self.hasO3):
            return -1

    def fz_roots(self, coef):
        if len(coef.shape) == 1:
            coef[~(np.isfinite(coef))] = 0.0
            rts = np.roots(coef[::-1])
            if rts.size == 0:
                printsafemulti('WARNING: fz_roots failed', self.logf, self.nps)
                rts = np.zeros(coef.size - 1)
            return rts

        else:
            rts = np.zeros((coef.shape[0], coef.shape[1] - 1), dtype=complex)
            coef[~(np.isfinite(coef))] = 0.0

            for i in range(coef.shape[0]):
                rts[i] = np.roots(coef[i][::-1])  # ::-1][0])
            return rts

    def setdustcorrect(self):
        global DUSTCORRECT
        DUSTCORRECT = True

    def unsetdustcorrect(self):
        global DUSTCORRECT
        DUSTCORRECT = False

    def dustcorrect(self, l1, l2, flux=False):
        #global DUSTCORRECT
        if DUSTCORRECT:
            if not flux:
                return 0.4 * self.mds['E(B-V)'] * (l1 - l2)
            return 10 ** (0.4 * self.mds['E(B-V)'] * (l1 - l2))
        else:
            if not flux:
                return 0
            return 1.0

    def setHab(self, Ha, Hb):
        self.Ha = Ha
        self.Hb = Hb
        if sum(self.Ha > 0):
            self.hasHa = True
        if sum(self.Hb > 0):
            self.hasHb = True

    def setOlines(self, O23727, O35007, O16300, O34959):
        self.O23727 = O23727
        self.O35007 = O35007
        self.O16300 = O16300

        if sum(self.O35007 > 0):
            self.hasO3 = True
        if sum(self.O23727 > 0):
            self.hasO2 = True

        if self.hasO2 and self.hasO3:
            self.O35007O2 = (self.O35007 / self.O23727) * self.dustcorrect(k_O3, k_O2, flux=True)
            self.O2O35007 = (self.O23727 / self.O35007) * self.dustcorrect(k_O2, k_O3, flux=True)

            self.logO35007O2 = np.log10(self.O35007O2)
            self.logO2O35007 = np.log10(self.O2O35007)

            #self.logO2O35007Hb=np.log10((self.O23727+self.O35007)/self.Hb)
            # ratios for other diagnostics - slightly different ratios needed
            if self.hasHb:
                self.logO2O35007Hb = np.log10((self.O23727 / self.Hb) * self.dustcorrect(k_O2, k_Hb, flux=True)) + \
                (self.O35007 / self.Hb) * self.dustcorrect(k_O35007, k_Hb, flux=True)

        else:
            printsafemulti("WARNING: needs O lines and  and Ha/b: did you run setHab()?", self.logf, self.nps)
        if self.hasHb:
            if self.hasO2:
                self.logO2Hb = np.log10(self.O23727 / self.Hb) + self.dustcorrect(k_O2, k_Hb)  # 0.4*self.mds['E(B-V)']*(k_O2-k_Hb) 
            if self.hasO3:
                self.O3Hb = (self.O35007 / self.Hb) + self.dustcorrect(k_O35007, k_Hb, flux=True)  # 0.4*self.mds['E(B-V)']*(k_O2-k_Hb)
                self.logO3Hb = np.log10(self.O3Hb)
                self.hasO3Hb = True


        if self.hasO2 and self.hasO3:
            self.OIII_OII = np.log10(self.O35007 / self.O23727 + self.dustcorrect(k_O35007, k_O2, flux=True))
            if O34959  is not None and sum(O34959 > 0) > 0:
                self.O34959p5007 = (O34959 + self.O35007)
                self.logO3O2 = np.log10((self.O34959p5007) / self.O23727) + self.dustcorrect(k_O3, k_O2)
                #this is useful when we get logq
                self.hasO3O2 = True
        if self.hasHb:
            self.OIII_Hb = np.log10(self.O35007 / self.Hb + self.dustcorrect(k_O35007, k_Hb, flux=True))

    def setNII(self, N26584):
        if N26584 is not None and sum(N26584 > 0):
            self.N26584 = N26584
            self.hasN2 = True
            if self.hasHa:
                self.logN2Ha = np.log10(self.N26584 / self.Ha)  # +self.dustcorrect(k_N2,k_Ha,flux=True) 
                #lines are very close: no dust correction
                #Note: no dust correction cause the lies are really close!
            else:
                printsafemulti("WARNING: needs NII6584 and Ha to calculate NIIHa: did you run setHab()?", self.logf, self.nps)
            if self.hasS2 and self.hasS26731 and self.hasN2:
                self.NII_SII = np.log10(self.N26584 / (self.S26717 + self.S26731))  # +self.dustcorrect(k_N2,k_S2,flux=True) 
                #lines are very close: no dust correction
            if self.hasO2 and self.hasN2:
                self.NII_OII = np.log10(self.N26584 / self.O23727 + self.dustcorrect(k_N2, k_O2, flux=True))

    def setSII(self, S26717, S26731, S39069, S39532):
        if S26717 is not None and sum(S26717 > 0) > 0:
            self.S26717 = S26717
            self.hasS2 = True

            if self.hasHa:
                self.logS2Ha = np.log10(self.S26717 / self.Ha) + self.dustcorrect(k_S2, k_Ha)
            else:
                printsafemulti("WARNING: needs SII6717 and Ha to calculate SIIHa: did you run setHab() and setS()?", self.logf, self.nps)
        if S26731 is not None and sum(S26731 > 1e-9) > 0:
            self.S26731 = S26731
            self.hasS26731 = True
        if S39069 is not None and sum(S39069 > 1e-9) > 0:
            self.S39069 = S39069
            self.hasS39069 = True
        if S39532 is not None and sum(S39532 > 1e-9) > 0:
            self.S39532 = S39532
            self.hasS39532 = True
        if self.hasS2:
            if self.hasN2 and self.NII_SII is None and self.hasS26731:
                self.NII_SII = np.log10(self.N26584 / (self.S26717 + self.S26731))  # +self.dustcorrect(k_N2,k_O2,flux=True) 
                #lines are very close: no dust correction
            if self.hasO3  and self.OIII_SII is None and self.hasS26731:
                self.OIII_SII = np.log10(self.O35007 / (self.S26717 + self.S26731) + self.dustcorrect(k_O3, k_S2, flux=True))

    #@profile
    def calcEB_V(self):
        printsafemulti("calculating E(B-V)", self.logf, self.nps)
        self.mds['E(B-V)'] = np.log10(2.86 * self.Hb / self.Ha) / (0.4 * (k_Ha - k_Hb))  # E(B-V)
        self.mds['E(B-V)'][self.mds['E(B-V)'] <= 0] = 1e-5

    #@profile
    def calcNIISII(self):
        if self.hasS2 and self.hasN2:
            self.N2S2 = self.N26584 / self.S26717 + self.dustcorrect(k_N2, k_S2, flux=True)  # 0.4*self.mds['E(B-V)']*(k_N2-k_S2)
            self.logN2S2 = np.log10(self.N26584 / self.S26717) + self.dustcorrect(k_N2, k_S2)  # 0.4*self.mds['E(B-V)']*(k_N2-k_S2)
            self.hasN2S2 = True
        else:
            printsafemulti("WARNING: needs SII6717 and NII6584 to calculate NIISII: did you run setN2() and setS?", self.logf, self.nps)

    #@profile
    def calcNIIOII(self):
        if self.hasN2 and self.hasO2:
            self.logN2O2 = np.log10(self.N26584 / self.O23727) + self.dustcorrect(k_N2, k_O2)
            self.hasN2O2 = True
        if not self.hasN2O2 or np.mean(self.logN2O2) < 1.2:

            try:
                printsafemulti('''WARNING: the KD02 and KK04 (+M08) methods should only be used for  log([NII]6564/[OII]3727) >1.2, 
                the mean log([NII]6564/[OII]3727)= %f''' % np.mean(self.logN2O2), self.logf, self.nps)
            except TypeError:
                printsafemulti('''WARNING: the KD02 and KK04 (+M08) methods 
                should only be used for  log([NII]6564/[OII]3727) >1.2, 
                the mean log([NII]6564/[OII]3727)= %s''' % self.logN2O2, self.logf, self.nps)

        if not self.hasN2O2:
            self.N2O2_roots = np.zeros(self.nm) + float('NaN')
        else:
            N2O2_coef = np.array([[self.N2O2_coef0, -532.15451, 96.373260, -7.8106123, 0.23928247]] * self.nm).T  # q=2e7 line (approx average)
            N2O2_coef[0] -= self.logN2O2
            N2O2_coef = N2O2_coef.T
            # finding roots for == (4)
            self.N2O2_roots = np.array([self.fz_roots(N2O2_coef)])[0]

    #@profile
    def calcR23(self):
        printsafemulti("calculating R23", self.logf, self.nps)

        #R23 NEW Comb, [NII]/Ha: KK04 = Kobulnicky & Kewley, 2004, submitted'
        if  self.hasO3   and self.hasO2 and self.hasHb:
            self.R2 = (self.O23727 / self.Hb) * self.dustcorrect(k_O2, k_Hb, flux=True)
            self.R3 = (self.O34959p5007 / self.Hb) * self.dustcorrect(k_O3, k_Hb, flux=True)
            self.R23 = self.R2 + self.R3
            self.logR23 = np.log10(self.R23)
            self.mds['logR23'] = self.logR23
            #note that values of logR23 > 0.95 are unphysical.
            #you may choose to uncomment the line below
            #self.logR23[self.logR23>0.95]=0.95
        else:
            printsafemulti("WARNING: need O3, O2, Hb", self.logf, self.nps)

    #@profile
    def calcS23(self):
        printsafemulti("calculating S23", self.logf, self.nps)
        #the original code here uses S267176731,
        #which is however set to 6717 as default
        #Vilchez & Esteban (1996)
        if  self.hasS2:
            if self.hasS39069 and self.hasHb:
                self.logS23 = np.log10((self.S26717 / self.Hb) *
                                       self.dustcorrect(k_S2, k_Hb, flux=True) +
                                       (self.S39069 / self.Hb) *
                                       self.dustcorrect(k_S3, k_Hb, flux=True))

            #self.logS3S2=np.log10(S39069/self.S26717)+self.dustcorrect(k_S3,k_S2)

            ##@profile
    def calclogq(self, Z):
        if not self.hasO3O2:
            printsafemulti("WARNING: needs O3,O2,Hb to calculate logq properly.", self.logf, self.nps)
            return -1
        if self.logO3O2sq is None:
            self.logO3O2sq = self.logO3O2 ** 2
        return (32.81 - 1.153 * self.logO3O2sq + Z * (-3.396 - 0.025 * self.logO3O2 + 0.1444 * self.logO3O2sq)) / (4.603 - 0.3119 * self.logO3O2 - \
                0.163 * self.logO3O2sq + Z * (-0.48 + 0.0271 * self.logO3O2 + 0.02037 * self.logO3O2sq))

    ##@profile
    def initialguess(self):
        # Initial Guess - appearing in LK code as of Nov 2006
        # upper branch: if no lines are available, metallicity is set to 8.7
        self.Z_init_guess = np.zeros(self.nm) + 8.7
        # use [N2]/Ha
        if self.hasHa and self.hasN2:
            self.Z_init_guess[(self.logN2Ha < -1.3) & (self.N26584 != 0.0)] = 8.2
            self.Z_init_guess[(self.logN2Ha < -1.1) & (self.logN2Ha >= -1.3) & (self.N26584 != 0.0)] = 8.4
            #A1 KE08
            self.Z_init_guess[(self.logN2Ha >= -1.1) & (self.N26584 != 0.0)] = 8.7
        #use [N2]/[O2]
        if self.hasN2 and self.hasO2:
            N2O2 = np.zeros(self.nm) + float('nan')
            if self.hasHb:
                ###FED CHECK THIS!
                N2O2 = self.N26584 * self.Ha * self.Hb * self.O23727
            if not self.hasN2O2:
                printsafemulti("WARNING: must calculate logN2O2 first", self.logf, self.nps)
                self.calcNIIOII()
            self.Z_init_guess[(self.logN2O2 < -1.2) & (N2O2 != 0.0)] = 8.2
            # at logN2O2<-1.2 using low-Z gals, A1 KE08
            self.Z_init_guess[(self.logN2O2 >= -1.2) & (N2O2 != 0.0)] = 8.7

            # at logN2O2>-1.2  using HII regions



            #######################these are the metallicity diagnostics##################
            #@profile
    def calcpyqz(self, plot=False, allD13=False):
        printsafemulti("calculating D13", self.logf, self.nps)

        # initializing variable pyqz to avoid style issues
        # (pyqz not defined is reported as error by Landscape.io w/ import in func
        pyqz = None
        try:
            import pyqz
        except ImportError:
            return -1

        #check version of pyqz
        from distutils.version import StrictVersion
        oldpyqz = False
        if StrictVersion(pyqz.__version__) <= StrictVersion('0.5.0'):
            oldpyqz = True

        if self.NII_SII is not None and allD13:
            if self.OIII_SII  is not None:

                if oldpyqz:
                    self.mds['D13_N2S2_O3S2'] = pyqz.get_qz(20, 'z', np.atleast_1d([self.NII_SII]), \
                                        np.atleast_1d([self.OIII_SII]), 'NII/SII', 'OIII/SII', \
                                        method='default', plot=plot, n_plot=False, savefig=False)[0].T
                else:
                    #pyqz.get_grid_fn(Pk=5.0,calibs='GCZO', kappa =20, struct='pp')
                    self.mds['D13_N2S2_O3S2'] = pyqz.interp_qz('Tot[O]+12', [np.atleast_1d([self.NII_SII]), \
                                                                          np.atleast_1d([self.OIII_SII])], \
                                                             '[NII]/[SII]+;[OIII]/[SII]+', \
                                                             show_plot=plot, n_plot=False, \
                                                             save_plot=False, verbose=False)[0].T

            if  self.OIII_Hb  is not None:
                if oldpyqz:
                    self.mds['D13_N2S2_O3Hb'] = pyqz.get_qz(20, 'z', np.atleast_1d([self.NII_SII]), \
                                        np.atleast_1d([self.OIII_Hb]), 'NII/SII', 'OIII/Hb', \
                                        method='default', plot=plot, n_plot=False, savefig=False)[0].T
                else:
                    self.mds['D13_N2S2_O3SHb'] = pyqz.interp_qz('Tot[O]+12', [np.atleast_1d([self.NII_SII]), \
                                                                          np.atleast_1d([self.OIII_Hb])], \
                                                             '[NII]/[SII]+;[OIII]/Hb', \
                                                             show_plot=plot, n_plot=False, \
                                                             save_plot=False, verbose=False)[0].T



            if  self.OIII_OII  is not None:
                if oldpyqz:
                    self.mds['D13_N2S2_O3O2'] = pyqz.get_qz(20, 'z', np.atleast_1d([self.NII_SII]), \
                                        np.atleast_1d([self.OIII_OII]), 'NII/SII', 'OIII/OII', \
                                        method='default', plot=plot, n_plot=False, savefig=False)[0].T
                else:
                    self.mds['D13_N2S2_O3O2'] = pyqz.interp_qz('Tot[O]+12', [np.atleast_1d([self.NII_SII]), \
                                                                          np.atleast_1d([self.OIII_OII])], \
                                                             '[NII]/[SII]+;[OIII]/[OII]+', \
                                                             show_plot=plot, n_plot=False, \
                                                             save_plot=False, verbose=False)[0].T

                    
        if self.NII_OII is not None and allD13:
            if self.OIII_SII  is not None:
                if oldpyqz:
                    self.mds['D13_N2O2_O3S2'] = pyqz.get_qz(20, 'z', np.atleast_1d([self.NII_OII]), \
                                        np.atleast_1d([self.OIII_SII]), 'NII/OII', 'OIII/SII', \
                                        method='default', plot=plot, n_plot=False, savefig=False)[0].T
                else:
                    self.mds['D13_N2O2_O3S2'] = pyqz.interp_qz('Tot[O]+12', [np.atleast_1d([self.NII_OII]), \
                                                                          np.atleast_1d([self.OIII_SII])], \
                                                             '[NII]/[OII]+;[OIII]/[SII]+', \
                                                             show_plot=plot, n_plot=False, \
                                                             save_plot=False, verbose=False)[0].T


            if  self.OIII_Hb  is not None:
                if oldpyqz:
                    self.mds['D13_N2O2_O3Hb'] = pyqz.get_qz(20, 'z', np.atleast_1d([self.NII_OII]), \
                                        np.atleast_1d([self.OIII_Hb]), 'NII/OII', 'OIII/Hb', \
                                        method='default', plot=plot, n_plot=False, savefig=False)[0].T
                else:
                    self.mds['D13_N2O2_O3Hb'] = pyqz.interp_qz('Tot[O]+12', [np.atleast_1d([self.NII_OII]), \
                                                                          np.atleast_1d([self.OIII_Hb])], \
                                                             '[NII]/[OII]+;[OIII]/Hb', \
                                                             show_plot=plot, n_plot=False, \
                                                             save_plot=False, verbose=False)[0].T

            if  self.OIII_OII  is not None:
                if oldpyqz:
                    self.mds['D13_N2O2_O3O2'] = pyqz.get_qz(20, 'z', np.atleast_1d([self.NII_OII]), \
                                        np.atleast_1d([self.OIII_OII]), 'NII/OII', 'OIII/OII', \
                                        method='default', plot=plot, n_plot=False, savefig=False)[0].T
                else:
                    self.mds['D13_N2O2_O3O2'] = pyqz.interp_qz('Tot[O]+12', [np.atleast_1d([self.NII_OII]), \
                                                                          np.atleast_1d([self.OIII_OII])], \
                                                             '[NII]/[OII]+;[OIII]/[OII]+', \
                                                             show_plot=plot, n_plot=False, \
                                                             save_plot=False, verbose=False)[0].T

        if self.logN2Ha is not None:
            if  self.OIII_Hb  is not None:
                if oldpyqz:
                    self.mds['D13_N2Ha_O3Hb'] = pyqz.get_qz(20, 'z', np.atleast_1d([self.logN2Ha]), \
                                        np.atleast_1d([self.OIII_Hb]), 'NII/Ha', 'OIII/Hb', \
                                        method='default', plot=plot, n_plot=False, savefig=False)[0].T
                else:
                    self.mds['D13_N2Ha_O3Hb'] = pyqz.interp_qz('Tot[O]+12', [np.atleast_1d([self.logN2Ha]), \
                                                                          np.atleast_1d([self.OIII_Hb])], \
                                                             '[NII]/Ha;[OIII]/Hb', \
                                                             show_plot=plot, n_plot=False, \
                                                             save_plot=False, verbose=False)[0].T

            if  self.OIII_OII  is not None:
                if oldpyqz:
                    self.mds['D13_N2Ha_O3O2'] = pyqz.get_qz(20, 'z', np.atleast_1d([self.logN2Ha]), \
                                        np.atleast_1d([self.OIII_OII]), 'NII/Ha', 'OIII/OII', \
                                        method='default', plot=plot, n_plot=False, savefig=False)[0].T

                else:
                    self.mds['D13_N2Ha_O3O2'] = pyqz.interp_qz('Tot[O]+12', [np.atleast_1d([self.logN2Ha]), \
                                                                          np.atleast_1d([self.OIII_Hb])], \
                                                             '[NII]/Ha;[OIII]/[OII]+', \
                                                             show_plot=plot, n_plot=False, \
                                                             save_plot=False, verbose=False)[0].T

    #@profile
    def calcDP00(self):
        # Diaz, A. I., & Perez-Montero, E. 2000, MNRAS, 312, 130
        # As per KD02: DP00 diagnostic systematically underestimates the
        # abundance relative to the comparison abundance.
        # A term is added to improve the fit according to KD02 Eq. 6
        # AVAILABLE BUT DEPRECATED
        printsafemulti("calculating DP00", self.logf, self.nps)

        if self.logS23  is None:
            self.calcS23()
            if self.logS23 is None:
                printsafemulti("WARNING: Cannot compute this without S23", self.logf, self.nps)
                return -1
        self.mds['DP00'] = 1.53 * self.logS23 + 8.27 + 1.0 / (2.0 - 9.0 * self.logS23 ** 3)

    #@profile
    def calcD02(self):
        # [NII]/Ha Denicolo, Terlevich & Terlevich (2002), MNRAS, 330, 69
        #FED:added uncertainties
        printsafemulti("calculating D02", self.logf, self.nps)

        e1 = np.random.normal(0, 0.05, self.nm)
        e2 = np.random.normal(0, 0.1, self.nm)
        if self.hasN2 and self.hasHa:
            self.mds['D02'] = 9.12 + e1 + (0.73 + e2) * self.logN2Ha
        else:
            printsafemulti("WARNING: need N2Ha to do this. did you run setHab and setNII", self.logf, self.nps)

    #@profile
    def calcPP04(self):
        ### PP04_N2_Z, PP04_O3N2_Z Pettini & Pagel diagnostics -
        ### Pettini & Pagel (2004), MNRAS, 348, L59
        # [NII]/Ha Pettini & Pagel (2004), MNRAS, 348, L59
        #discriminating lower and upper branch using  [NII]/[OII] or  [NII]/Ha
        printsafemulti("calculating PP04", self.logf, self.nps)
        if self.hasN2 and self.hasHa:
            self.mds['PP04_N2Ha'] = nppoly.polyval(self.logN2Ha, [9.37, 2.03, 1.26, 0.32])

            #FED: restricting the range as per paper
            index = (self.logN2Ha > -2.5) * (self.logN2Ha < -0.3)
            self.mds['PP04_N2Ha'][~index] = float('NaN')
            if self.hasO3Hb:
                self.mds['PP04_O3N2'] = 8.73 - 0.32 * (self.logO3Hb - self.logN2Ha)
                index = (self.logO3Hb > 2)
                self.mds['PP04_O3N2'][index] = float('NaN')
            else:
                printsafemulti("WARNING: need O3Hb for PP04_O3N2", self.logf, self.nps)
        else:
            printsafemulti("WARNING: need N2Ha to do this. did you run setHab and setNII", self.logf, self.nps)

    #@profile
    def calcZ94(self):
        ### calculating z from Kobulnicky,Kennicutt,Pizagno (1998)
        ### parameterization of Zaritzky et al. (1994)
        ###Z94 = Zaritsky, D., Kennicutt, R. C., & Huchra, J. P. 1994,
        ###ApJ, 420, 87
        ### only valid on the upper branch of R23 (KE08 A2.4)

        printsafemulti("calculating Z94", self.logf, self.nps)
        if self.logR23 is None:
            printsafemulti("WARNING: Must first calculate R23", self.logf, self.nps)
            self.calcR23()
            if self.logR23 is None:
                printsafemulti("WARNING: Cannot compute this without R23", self.logf, self.nps)
                return -1
        self.mds['Z94'] = nppoly.polyval(self.logR23, [9.265, -0.33, -0.202, -0.207, -0.333])
        self.mds['Z94'][(self.logR23 > 0.9)] = None
        ## 0.9 is a conservative constraint to make sure that we are
        ## only using the upper branch (i.e. 12+log(O/H)>8.4)

    def calcP(self):
        if self.P is None:
            if self.logR23 is None:
                printsafemulti("WARNING: Must first calculate R23", self.logf, self.nps)
                self.calcR23()
                if self.logR23 is None:
                    printsafemulti("WARNING: Cannot compute this without R23", self.logf, self.nps)
                    return -1
            #R3=10**self.logO349595007Hb
            #R2=10**self.logO2Hb
            #P = R3/(R2+R3)
            self.P = self.R3 / self.R23

    #@profile
    def calcP05(self):
        # #### P-method #####
        ##Pilyugin+ 2005 method.  Based on [OIII],[OII], Hbeta
        ##calibrated from Te method
        # make sure you run setOlines() first
        printsafemulti("calculating P05", self.logf, self.nps)

        if self.calcP() == -1:
            return -1
        if self.Z_init_guess is None:
            self.initialguess()

        Psq = self.P * self.P

        P_abund_up = (self.R23 + 726.1 + 842.2 * self.P + 337.5 * Psq) / (85.96 + 82.76 * self.P + 43.98 * Psq + 1.793 * self.R23)
        P_abund_low = (self.R23 + 106.4 + 106.8 * self.P - 3.40 * Psq) / (17.72 + 6.60 * self.P + 6.95 * Psq - 0.302 * self.R23)

        self.mds['P05'] = P_abund_up
        self.mds['P05'][self.Z_init_guess < 8.4] = P_abund_low[self.Z_init_guess < 8.4]

    #@profile
    def calcP10(self):
        # #### P-method #####
        ##Pilyugin+ 2010 method.
        ##calibrated from Te method
        # need Hb
        #The Astrophysical Journal, Volume 720, Issue 2, pp. 1738-1751 (2010).
        #Published in Sep 2010

        printsafemulti("calculating P10", self.logf, self.nps)

        if not self.hasHb:
            printsafemulti("this method needs Hb", self.logf, self.nps)
            return -1
        self.mds['P10_ONS'] = np.zeros(self.nm) + float('NaN')
        self.mds['P10_ON'] = np.zeros(self.nm) + float('NaN')
        #P10N2=np.zeros(self.nm)+float('NaN')
        #P10S2=np.zeros(self.nm)+float('NaN')
        P10logR3 = np.zeros(self.nm) + float('NaN')
        P10logR2 = np.zeros(self.nm) + float('NaN')
        P10logN2 = np.zeros(self.nm) + float('NaN')
        P10logS2 = np.zeros(self.nm) + float('NaN')

        self.calcP()
        if self.R2 is not None:
            P10logR2 = np.log(self.R2)

        if self.R3 is not None:
            P10logR3 = np.log(self.R3)

        if self.hasN2:
            #the ratio of N26548 and N26548 is N26584/N26548 = 3
            #independent on physical conditions
            #The Physics and Dynamics of Planetary Nebulae
            # By Grigor A. Gurzadyan
            P10logN2 = np.log((self.N26584 * 1.33) / self.Hb) + self.dustcorrect(k_N2, k_Hb)

        if self.hasS2 and self.hasS26731:
            self.S2Hb = ((self.S26717 + self.S26731) / self.Hb) + self.dustcorrect(k_S2, k_Hb, flux=True)
            self.hasS2Hb = True
            P10logS2 = np.log10(self.S2Hb)

        P10logN2S2 = P10logN2 - P10logS2
        P10logN2R2 = P10logN2 - P10logR2
        P10logS2R2 = P10logS2 - P10logR2

        coefsONS0 = np.array([8.277, 0.657, -0.399, -0.061, 0.005])
        coefsONS1 = np.array([8.816, -0.733, 0.454, 0.710, -0.337])
        coefsONS2 = np.array([8.774, -1.855, 1.517, 0.304, 0.328])

        vsONS = np.array([np.ones(self.nm), self.P, P10logR3, P10logN2R2, P10logS2R2]).T

        coefsON0 = np.array([8.606, -0.105, -0.410, -0.150])
        coefsON1 = np.array([8.642, 0.077, 0.411, 0.601])
        coefsON2 = np.array([8.013, 0.905, 0.602, 0.751])

        vsON = np.array([np.ones(self.nm), P10logR3, P10logR2, P10logN2R2]).T

        indx = P10logN2 > -0.1
        if self.P is not None:
            self.mds['P10_ONS'][indx] = np.dot(vsONS[indx], coefsONS0)
        self.mds['P10_ON'][indx] = np.dot(vsON[indx], coefsON0)

        indx = (P10logN2 < -0.1) * (P10logN2S2 > -0.25)
        if self.P is not None:
            self.mds['P10_ONS'][indx] = np.dot(vsONS[indx], coefsONS1)
        self.mds['P10_ON'][indx] = np.dot(vsON[indx], coefsON1)

        indx = (P10logN2 < -0.1) * (P10logN2S2 < -0.25)
        if self.P is not None:
            self.mds['P10_ONS'][indx] = np.dot(vsONS[indx], coefsONS2)
        self.mds['P10_ON'][indx] = np.dot(vsON[indx], coefsON2)

        indx = ~((self.mds['P10_ONS'] > 7.1) * (self.mds['P10_ON'] > 7.1) * (self.mds['P10_ONS'] < 9.4) * (self.mds['P10_ON'] < 9.4))
        if self.P is not None:
            self.mds['P10_ONS'][indx] = float('NaN')
        self.mds['P10_ON'][indx] = float('NaN')

    #@profile
    def calcP01(self):
        # P-method 2001 upper branch (derecated and commented out)
        # Pilyugin 2001
        # available but deprecated
        printsafemulti("calculating old P05", self.logf, self.nps)

        if self.Z_init_guess is None:
            self.initialguess()
        if self.hasO3O2 and self.hasO3  and self.hasO2:
            P = 10 ** self.logO3O2 / (1 + 10 ** self.logO3O2)
            if self.logR23 is None:
                printsafemulti("WARNING: Must first calculate R23", self.logf, self.nps)
                self.calcR23()
                if self.logR23 is None:
                    printsafemulti("WARNING: Cannot compute this without R23", self.logf, self.nps)
                    return -1
            Psq = P ** 2
            P_abund_old = (self.R23 + 54.2 + 59.45 * P + 7.31 * Psq) / (6.07 + 6.71 * P + 0.371 * Psq + 0.243 * self.R23)
            self.mds['P01'] = np.zeros(self.nm) + float('NaN')
            self.mds['P01'][self.Z_init_guess >= 8.4] = P_abund_old[self.Z_init_guess >= 8.4]
        else:
            printsafemulti("WARNING: need OIIIOII to calculate P01, did you set them up with  setOlines()?", self.logf, self.nps)

    #@profile
    def calcC01_ZR23(self):
        # C01 = Charlot, S., & Longhetti, M., 2001, MNRAS, 323, 887
        # Charlot 01 R23 calibration: (case F) ##
        # available but deprecated
        printsafemulti("calculating C01", self.logf, self.nps)

        if self.hasO3 and self.hasO2 and self.hasO3Hb:
            x2 = self.O2O35007 / 1.5
            x3 = (10 ** self.logO3Hb) * 0.5
            self.mds['C01_R23'] = np.zeros(self.nm) + float('NaN')
            self.mds['C01_R23'][self.O2O35007 < 0.8] = np.log10(3.78e-4 * (x2[self.O2O35007 < 0.8]) ** 0.17 * x3[self.O2O35007 < 0.8] ** (-0.44)) + 12.0

            self.mds['C01_R23'][self.O2O35007 >= 0.8] = np.log10(3.96e-4 * x3[self.O2O35007 >= 0.8] ** (-0.46)) + 12.0
        else:
            printsafemulti('''WARNING: need [OIII]5700, [OII]3727, and Ha to calculate calcC01_ZR23, 
did you set them up with  setOlines()?''', self.logf, self.nps)

        # Charlot 01 calibration: (case A) based on [N2]/[SII]##
        # available but deprecated
        if not self.hasN2S2:
            printsafemulti("WARNING: trying to calculate logNIISII", self.logf, self.nps)
            self.calcNIISII()
        if self.hasN2S2 and self.hasO3 and self.hasO2 and self.hasO3Hb:
            self.mds['C01_N2S2'] = np.log10(5.09e-4 * (x2 ** 0.17) * ((self.N2S2 / 0.85) ** 1.17)) + 12
        else:
            printsafemulti('''WARNING: needs [NII]6584, [SII]6717, [OIII]5700, [OII]3727, and Ha to calculate calcC01_ZR23, 
did you set them up with  setOlines() and ?''', self.logf, self.nps)

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

        printsafemulti("calculating M91", self.logf, self.nps)
        self.mds['M91'] = np.zeros(self.nm) + float('NaN')

        if self.logR23 is None:
            printsafemulti("WARNING: Must first calculate R23", self.logf, self.nps)
            self.calcR23()
            if self.logR23 is None:
                printsafemulti("WARNING: Cannot compute this without R23", self.logf, self.nps)
                return -1

        if self.Z_init_guess is None:
            self.initialguess()

        M91_Z_low = nppoly.polyval(self.logR23, [12.0 - 4.944, 0.767, 0.602]) - \
                   self.logO3O2 * nppoly.polyval(self.logR23, [0.29, 0.332, -0.331])
        M91_Z_up = nppoly.polyval(self.logR23, [12.0 - 2.939, -0.2, -0.237, -0.305, -0.0283]) - \
                   self.logO3O2 * nppoly.polyval(self.logR23, [0.0047, -0.0221, -0.102, -0.0817, -0.00717])

        indx = (np.abs(self.logO3O2) > 0) * (np.abs(self.logR23) > 0) * (self.Z_init_guess < 8.4)
        self.mds['M91'][indx] = M91_Z_low[indx]
        indx = (np.abs(self.logO3O2) > 0) * (np.abs(self.logR23) > 0) * (self.Z_init_guess >= 8.4)
        self.mds['M91'][indx] = M91_Z_up[indx]
        self.mds['M91'][(M91_Z_up < M91_Z_low)] = float('NaN')

    #@profile
    def calcM13(self):
        #Marino+ 2013
        printsafemulti("calculating M13", self.logf, self.nps)

        if not self.hasHa  or not self.hasN2:
            printsafemulti("WARNING: need O3, N2, Ha and Hb, ",
                           "or at least N2 and Ha", self.logf, self.nps)
            return -1
        else:
            e1 = np.random.normal(0, 0.027, self.nm)
            e2 = np.random.normal(0, 0.024, self.nm)
            self.mds["M13_N2"] = 8.743 + e1 + (0.462 + e2) * self.logN2Ha
            if   self.hasHb and self.hasO3:
                e1 = np.random.normal(0, 0.012, self.nm)
                e2 = np.random.normal(0, 0.012, self.nm)
                O3N2 = self.logO3Hb - self.logN2Ha
                self.mds["M13_O3N2"] = 8.533 + e1 - (0.214 + e1) * O3N2
                index = (O3N2 > 1.7) 
                self.mds["M13_O3N2"][index] = float('NaN')
                index = (O3N2 < -1.1)
                self.mds["M13_O3N2"][index] = float('NaN')

    #@profile
    def calcM08(self, allM08=False):
        #Maiolino+ 2008
        #Astronomy and Astrophysics, Volume 488, Issue 2, 2008, pp.463-479
        #Published in Sep 2008
        printsafemulti("calculating M08", self.logf, self.nps)
        highZ = None
        if self.logO35007O2 is not None:
            self.mds['M08_O3O2'] = np.zeros(self.nm) + float('NaN')
            coefs = np.array([M08_coefs['O3O2']] * self.nm).T
            coefs[0] = coefs[0] - self.logO35007O2
            sols = np.array([self.fz_roots(coefs.T)])[0] + 8.69
            indx = ((sols.real >= 7.1) * (sols.real <= 9.4) * (sols.imag == 0)).cumsum(1).cumsum(1) == 1
            #the two cumsum assure that if the condition for the ith element
            #of indx is [False, False] then after the first cumsum(1) is [0,0]
            #[False, True] is [0,1]
            #[True, True] is [1,2]
            #but (here is the kicker) [True, False] is [1,1].
            #Because i want only one solution
            #(i'll settle for the first one occurring) [1,1] is ambiguous.
            #The second cumsum(1) makes
            #[0,0]->[0,0], [0,1]->[0,1], [1,2]->[1,3] and finally [1,1]->[1,2]

            self.mds['M08_O3O2'][(indx.sum(1)) > 0] = sols[indx].real
            highZ = np.median(self.logO35007O2) < 0
        if self.logN2Ha is not None:
            self.mds['M08_N2Ha'] = np.zeros(self.nm) + float('NaN')
            coefs = np.array([M08_coefs['N2Ha']] * self.nm).T
            coefs[0] = coefs[0] - self.logN2Ha
            sols = np.array([self.fz_roots(coefs.T)])[0] + 8.69
            indx = ((sols.real >= 7.1) * (sols.real <= 9.4) * (sols.imag == 0)).cumsum(1).cumsum(1) == 1
            self.mds['M08_N2Ha'][(indx.sum(1)) > 0] = sols[indx].real
            if highZ is None:
                highZ = np.median(self.logN2Ha) > -1.3

        if self.logR23 is None:
            printsafemulti("WARNING: Must first calculate R23", self.logf, self.nps)
            self.calcR23()
        if self.logR23 is None:
            printsafemulti("WARNING: Cannot compute M08_R23 without R23", self.logf, self.nps)
        else:
            self.mds['M08_R23'] = np.zeros(self.nm) + float('NaN')
            coefs = np.array([M08_coefs['R23']] * self.nm).T
            coefs[0] = coefs[0] - self.logR23
            sols = np.array([self.fz_roots(coefs.T)])[0] + 8.69
            if highZ is True:
                indx = ((sols.real >= 7.1) * (sols.real <= 9.4) * (sols.imag == 0) * (sols.real >= 8.0)).cumsum(1).cumsum(1) == 1
                self.mds['M08_R23'][(indx.sum(1)) > 0] = sols[indx].real
            elif highZ is False:
                indx = ((sols.real >= 7.1) * (sols.real <= 9.4) * (sols.imag == 0) * (sols.real <= 8.0)).cumsum(1).cumsum(1) == 1
                self.mds['M08_R23'][(indx.sum(1)) > 0] = sols[indx].real
        if not allM08:
            return
        else:
            printsafemulti("calculating other M08s", self.logf, self.nps)

        if self.logO3Hb is not None:

            self.mds['M08_O3Hb'] = np.zeros(self.nm) + float('NaN')
            coefs = np.array([M08_coefs['O3Hb']] * self.nm).T
            coefs[0] = coefs[0] - self.logO3Hb
            sols = np.array([self.fz_roots(coefs.T)])[0] + 8.69
            if highZ is True:
                indx = ((sols.real >= 7.1) * (sols.real <= 9.4) * (sols.imag == 0) * (sols.real >= 7.9)).cumsum(1).cumsum(1) == 1
                self.mds['M08_O3Hb'][(indx.sum(1)) > 0] = sols[indx].real
            elif highZ is False:
                indx = ((sols.real >= 7.1) * (sols.real <= 9.4) * (sols.imag == 0) * (sols.real <= 7.9)).cumsum(1).cumsum(1) == 1
                self.mds['M08_O3Hb'][(indx.sum(1)) > 0] = sols[indx].real

        if self.logO2Hb is not None:
            self.mds['M08_O2Hb'] = np.zeros(self.nm) + float('NaN')
            coefs = np.array([M08_coefs['O2Hb']] * self.nm).T
            coefs[0] = coefs[0] - self.logO2Hb
            sols = np.array([self.fz_roots(coefs.T)])[0] + 8.69
            if highZ is True:
                indx = ((sols.real >= 7.1) * (sols.real <= 9.4) * (sols.imag == 0) * (sols.real >= 8.7)).cumsum(1).cumsum(1) == 1
                self.mds['M08_O2Hb'][(indx.sum(1)) > 0] = sols[indx].real
            elif highZ is False:
                indx = ((sols.real >= 7.1) * (sols.real <= 9.4) * (sols.imag == 0) * (sols.real <= 8.7)).cumsum(1).cumsum(1) == 1
                self.mds['M08_O2Hb'][(indx.sum(1)) > 0] = sols[indx].real



        if self.hasO3  and self.hasN2:
            self.mds['M08_O3N2'] = np.zeros(self.nm) + float('NaN')
            coefs = np.array([M08_coefs['O3N2']] * self.nm).T
            coefs[0] = coefs[0] - np.log(self.O35007 / self.N26584) * self.dustcorrect(k_O35007, k_N2)
            sols = np.array([self.fz_roots(coefs.T)])[0] + 8.69
            indx = ((sols.real >= 7.1) * (sols.real <= 9.4) * (sols.imag == 0)).cumsum(1).cumsum(1) == 1
            self.mds['M08_O3N2'][(indx.sum(1)) > 0] = sols[indx].real

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
        # from a 7 dimensional if/for loop to 1 if and 1 for :-D
        #vectorizing makes fed happy ...

        printsafemulti("calculating KD02_N2O2", self.logf, self.nps)

        if self.hasN2 and self.hasO2 and self.hasHa and self.hasHb:
            self.mds['KD02_N2O2'] = np.zeros(self.nm) + float('NaN')
            if not self.hasN2O2:
                printsafemulti("WARNING: must calculate logN2O2 first", self.logf, self.nps)
                self.calcNIIOII()
            if  not self.hasN2O2 or self.N2O2_roots  is None or sum(np.isnan(self.N2O2_roots.flatten())) == len(self.N2O2_roots.flatten()):
                printsafemulti("WARNING:  cannot calculate N2O2", self.logf, self.nps)
                return -1
            roots = self.N2O2_roots.T
            for k in range(4):
                indx = (abs(roots[k]) >= 7.5) * (abs(roots[k]) <= 9.4) * (roots[k][:].imag == 0.0)
                self.mds['KD02_N2O2'][indx] = abs(roots[k][indx])
        else:
            printsafemulti("WARNING: need NII6584 and OII3727 and Ha and Hb to calculate this. did you run setO() setHab() and setNII()?", self.logf, self.nps)
        return 1

    #@profile
    def calcKK04_N2Ha(self):
        # calculating [N2]/Ha abundance estimates using [O3]/[O2] also
        printsafemulti("calculating KK04_N2Ha", self.logf, self.nps)

        if self.mds['KD02_N2O2']  is None:
            self.calcKD02_N2O2()
        if self.mds['KD02_N2O2']  is None or sum(np.isnan(self.mds['KD02_N2O2'])) == self.nm:
            printsafemulti("WARNING: without KD02_N2O2 cannot calculate KK04_N2Ha properly, but we will do our best...", self.logf, self.nps)
            Z_new_N2Ha = np.zeros(self.nm) + 8.6
        else:
            Z_new_N2Ha = self.mds['KD02_N2O2'].copy()  # was 8.6

        if self.hasN2 and self.hasHa:
            logq_save = np.zeros(self.nm)
            convergence, tol, ii = 100, 1.0e-3, 0
            if self.hasO3O2:
                # calculating logq using the [N2]/[O2]
                # metallicities for comparison
                while convergence > tol and ii < 100:
                    ii += 1
                    self.logq = self.calclogq(Z_new_N2Ha)
                    Z_new_N2Ha = nppoly.polyval(self.logN2Ha, [7.04, 5.28, 6.28, 2.37]) - \
                                self.logq * nppoly.polyval(self.logN2Ha, [-2.44, -2.01, -0.325, +0.128]) + \
                                10 ** (self.logN2Ha - 0.2) * self.logq * (-3.16 + 4.65 * self.logN2Ha)
                    convergence = np.abs(self.logq - logq_save).mean()
                    logq_save = self.logq.copy()
                if ii >= 100:
                    printsafemulti("WARNING: loop did not converge", self.logf, self.nps)
                    Z_new_N2Ha = np.zeros(self.nm) + float('NaN')
            else:
                self.logq = 7.37177 * np.ones(self.nm)
                Z_new_N2Ha = nppoly.polyval(self.logN2Ha, [7.04, 5.28, 6.28, 2.37]) - \
                            self.logq * nppoly.polyval(self.logN2Ha, [-2.44, -2.01, -0.325, +0.128]) + \
                            10 ** (self.logN2Ha - 0.2) * self.logq * (-3.16 + 4.65 * self.logN2Ha)
            self.mds['KK04_N2Ha'] = Z_new_N2Ha
            indx = self.logN2Ha > 0.8
            self.mds['KK04_N2Ha'][indx] = float('NaN')
        else:
            printsafemulti("WARNING: need NII6584  and Ha to calculate this. did you run  setHab() and setNII()?", self.logf, self.nps)

    #@profile
    def calcKK04_R23(self):
        # Kobulnicky & Kewley 2004
        # calculating upper and lower metallicities for objects without
        # Hb  and for objects without O3 and/or O2

        printsafemulti("calculating KK04_R23", self.logf, self.nps)
        #this is in the original code but not used :(
        #if self.hasN2 and self.hasHa:
        #logq_lims=[6.9,8.38]
        #logN2Ha=np.log10(self.N26584/self.Ha) CHECK!! why remove dust correction??
        #Z_new_N2Ha_lims= np.atleast_2d([1.0,1.0]).T*nppoly.polyval(self.logN2Ha,[7.04, 5.28,6.28,2.37])-
        #np.atleast_2d( logq_lims).T*nppoly.polyval(self.logN2Ha,[-2.44,-2.01,-0.325,0.128])+
        #np.atleast_2d(logq_lims).T*(10**(self.logN2Ha-0.2)*(-3.16+4.65*self.logN2Ha))
        # R23 diagnostics from Kobulnicky & Kewley 2004

        Zmax = np.zeros(self.nm)
        # ionization parameter form logR23
        if not self.hasO3O2:
            logq = np.zeros(self.nm)
        else:
            if self.Z_init_guess is None:
                self.initialguess()
            Z_new = self.Z_init_guess.copy()
            if self.logR23 is None:
                printsafemulti("WARNING: Must first calculate R23", self.logf, self.nps)
                self.calcR23()
            if self.logR23 is None:
                printsafemulti("WARNING: Cannot compute this without R23", self.logf, self.nps)
            else:
                logqold, convergence, ii = np.zeros(self.nm) + 100, 100, 0
                tol = 1e-4
                #3 iterations are typically enought to achieve convergence KE08 A2.3
                while convergence > tol and ii < 100:
                    Zmax = Zmax * 0.0
                    ii += 1
                    logq = self.calclogq(Z_new)
                    Zmax[(logq >= 6.7) * (logq < 8.3)] = 8.4
                    # maximum of R23 curve:
                    Z_new = nppoly.polyval(self.logR23, [9.72, -0.777, -0.951, -0.072, -0.811]) - \
                           logq * nppoly.polyval(self.logR23, [0.0737, -0.0713, -0.141, 0.0373, -0.058])
                    indx = self.Z_init_guess <= Zmax
                    Z_new[indx] = nppoly.polyval(self.logR23[indx], [9.40, 4.65, -3.17]) - \
                                 logq[indx] * nppoly.polyval(self.logR23[indx], [0.272, 0.547, -0.513])
                    convergence = np.abs((logqold - logq).mean())
                    logqold = logq.copy()
                if ii >= 100:
                    printsafemulti("WARNING: loop did not converge", self.logf, self.nps)
                    Z_new = np.zeros(self.nm) + float('NaN')
                Z_new_lims = [nppoly.polyval(self.logR23, [9.40, 4.65, -3.17]) - \
                            logq * nppoly.polyval(self.logR23, [0.272, 0.547, -0.513]),
                            nppoly.polyval(self.logR23, [9.72, -0.777, -0.951, -0.072, -0.811]) - \
                            logq * nppoly.polyval(self.logR23, [0.0737, -0.0713, -0.141, 0.0373, -0.058])]
                Z_new[(Z_new_lims[0] > Z_new_lims[1])] = None
                self.mds['KK04_R23'] = Z_new

    #@profile
    def calcKDcombined(self):
        # KD02comb  Kewley, L. J., & Dopita, M. A., 2002, ApJ
        # updated in KE08
        # ### KD02 [NII]/[OII] estimate ###
        # (can be used for log(O/H)+12 > 8.6 only)

        printsafemulti("calculating KD_combined", self.logf, self.nps)

        #We first use the
        #[N ii]/[O ii] ratio to determine whether it lies on the upper
        #or lower R23 branch

        if self.mds['KD02_N2O2']  is None:
            self.calcKD02_N2O2()
        if self.mds['KK04_N2Ha'] is None:
            self.calcKK04_N2Ha()
        if self.logR23  is None:
            self.calcR23()
        if self.mds['M91']  is None:
            printsafemulti("WARNING:  Must first calculate M91", self.logf, self.nps)
            self.calcM91()
#        if self.mds['Z94']  is None:
#            printsafemulti(  "WARNING:  Must first calculate Z94",self.logf,self.nps)
#            self.calcZ94()
        if self.mds['KK04_R23']  is None:
            printsafemulti("WARNING:  Must first calculate KK04_R23", self.logf, self.nps)
            self.calcKK04_R23()
        if not self.hasHa and not self.hasHb:
            printsafemulti("WARNING: need Ha and Hb for this. did you run setHab()?", self.logf, self.nps)

        #alternative way to calculate KD02_N2O2, stated in the paper KD02,
        #valid in high Z regimes (Z>8.4)
        #but we forego it
        #if not self.logN2O2 is None:
        #    self.mds['KD02_N2O2']=np.log10(8.511e-4*(1.54020+1.26602*self.logN2O2+0.167977*self.logN2O2**2))+12.
        #else: self.mds['KD02_N2O2']=np.zeros(self.nm)+float('NaN')

        # ionization parameter
        # calculate an initial ionization parameter by assuming
        # a nominal lower branch [12 + log (O/H ) = 8.2]
        # or upper branch [12 + log (O/H ) = 8.7] metallicity using
        # equation (13) from KK04
        logq = np.zeros(self.nm)
        if self.hasN2 and self.hasO2 and self.hasHb and self.hasHa and self.hasO3O2:
            logq = self.calclogq(self.mds['KD02_N2O2'])
            logq[self.mds['KD02_N2O2'] >= 8.4] = self.logq[self.mds['KD02_N2O2'] >= 8.4]
        else:
            if self.Z_init_guess is None:
                self.initialguess()
            logq = self.calclogq(self.Z_init_guess)
        #FED: CHECK: the paragraph below makes sense in words but i dont see whereit ie enforced.
        # if log([NII]/[OII]) after extinction correction is <-1.5, then check the data.
        # if it is only slightly less than 1.5, then this can be a result of either noisy
        # data, inaccurate fluxes or extinction correction, or a higher ionization parameter
        # than modelled.
        # For these cases, the average of the M91,Z94 and C01 should be used.

        # KD02 R23 estimate (not reliable for  8.4 < log(O/H)+12 < 8.8)
        # uses [NII]/[OII] estimate as initial guess - this can be changed below

        self.mds['KD02comb'] = np.zeros(self.nm) + float('NaN')

        indx_ig = self.Z_init_guess > 8.4
        if self.mds['KD02_N2O2'] is not None:
            self.mds['KD02comb'][indx_ig] = self.mds['KD02_N2O2'][indx_ig].copy()
        if self.mds['KK04_N2Ha'] is not None:
            self.mds['KD02comb'][~indx_ig] = self.mds['KK04_N2Ha'][~indx_ig].copy()
        if self.mds['KK04_R23'] is not None and self.mds['M91'] is not None:
            # if [NII]/[OII] abundance available
            # and [NII]/Ha abundance < 8.4, then use R23.
            indx = (~np.isnan(self.mds['KK04_R23'])) * (~np.isnan(self.mds['M91'])) * (~indx_ig)
            self.mds['KD02comb'][indx] = 0.5 * (self.mds['KK04_R23'][indx].copy() + self.mds['M91'][indx].copy())

        else:
            printsafemulti("WARNING:  cannot calculate KK04comb because  KK04_R23 or M91, failed", self.logf, self.nps)

#######################these are the metallicity diagnostics##################
#@profile
    def calcPM14(self):
        # Perez-Montero 2014
        # (can be used for for log(O/H)+12 > 8.6 only)
        import os
        from subprocess import Popen, PIPE, STDOUT
        from StringIO import StringIO

        printsafemulti("calculating HIICHI", self.logf, self.nps)
        fin_hii_chi = open(os.getenv('HIICHI_DIR') + '/in.tmp', 'w')

        if not self.hasHb:
            printsafemulti("cannot calculate HIICHI without Hbeta", self.logf, self.nps)
            return -1

        ratios = np.zeros((5, self.nm))

        if self.R2 is not None:
            ratios[0] = self.R2
        elif self.hasO2:
            ratios[0] = ((self.O23727 / self.Hb) * self.dustcorrect(k_O2, k_Hb, flux=True))
        else:
            ratios[0] = np.array(['0 '] * self.nm)

        #we will never have 4363...
        ratios[1] = np.zeros(self.nm)

        if self.hasO3Hb:
            ratios[2] = self.O3Hb
        elif self.hasO3:
            ratios[2] = ((self.O35007 / self.Hb) + self.dustcorrect(k_O35007, k_Hb, flux=True))  # 0.4*self.mds['E(B-V)']*(k_O2-k_Hb) 
        else:
            ratios[2] = np.zeros(self.nm)

        if self.hasN2:
            ratios[3] = self.N26584 / self.Hb
        else:
            ratios[3] = np.zeros(self.nm)

        if self.hasS2Hb:
            ratios[4] = self.S2Hb
        elif self.hasS2 and self.hasS26731:
            ratios[4] = (((self.S26717 + self.S26731) / self.Hb) + self.dustcorrect(k_S2, k_Hb, flux=True))
        else:
            ratios[4] = np.zeros(self.nm)

        
        for ni in range(self.nm):
            fin_hii_chi.write('%f %f %f %f %f\n' % (ratios[0][ni], ratios[1][ni], ratios[2][ni], ratios[3][ni], ratios[4][ni]))
        fin_hii_chi.close()
        os.system("ln -s %s/C13*dat . " % os.getenv('HIICHI_DIR'))
        print ("\n\n\n\n\n")
        #os.system("python %s/HII-CHI-mistry_v01.2.py in.tmp"%os.getenv('HIICHI_DIR'))
        #os.system("python %s/HII-CHI-mistry_v01.2.py %s/in.tmp"%(os.getenv('HIICHI_DIR'),os.getenv('HIICHI_DIR')))
        p = Popen(['python', '%s/HII-CHI-mistry_v01.2.py' % os.getenv('HIICHI_DIR'), '%s/in.tmp' % os.getenv('HIICHI_DIR')], stdout=PIPE, stdin=PIPE, stderr=STDOUT)
        out, err = p.communicate(input='%s/in.tmp' % os.getenv('HIICHI_DIR'))
        print ("\n\n\n\n\n")
        out = StringIO(out)
        #       for l in enumerate(out):
        #           if l[0].isdigit():
        #               break
        out = out.readlines()[12:]
        self.mds['PM14'] = np.zeros((self.nm))
        self.mds['PM14err'] = np.zeros((self.nm))
        for i, l in enumerate(out):
            self.mds['PM14'][i], self.mds['PM14err'][i] = map(float, l.split()[3:5])
        #data =  np.loadtxt(out,  skiprows=12, usecols=(3,4))#, dtype=[('lOH','f'),('elOH','f')], delimiter=",", unpack = True)
        print (self.mds)
        os.system('rm -r C13*dat')

    #@profile
    def calcD16(self):
        #Dopita+ 2016
        printsafemulti("calculating D16", self.logf, self.nps)
        
        if not self.hasHa  or not self.hasN2 or not self.hasS2:
            printsafemulti("WARNING: need N2, Ha and SII, ",
                           self.logf, self.nps)
            return -1
        y = self.logN2S2 + 0.264 * self.logN2Ha
        self.mds["D16"] = 8.77 + y - 0.45 * pow(y + 0.3, 5)
        index = (y < -1.)
        self.mds["D16"][index] = float('NaN')                
        index = (y > 0.5)
        self.mds["D16"][index] = float('NaN')                

        
