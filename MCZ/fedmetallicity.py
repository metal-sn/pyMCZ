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


    diags.calcKK04()
