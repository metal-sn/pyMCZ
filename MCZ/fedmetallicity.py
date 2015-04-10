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
import sys,os
import numpy as np
from pylab import hist,show


IGNOREDUST=False
MP=True

##list of metallicity methods, in order calculated
Zs=["E(B-V)", #Halpha, Hbeta
    "logR23", #Hbeta,  [OII]3727, [OIII]5007, [OIII]4959

    "D02",    #Halpha, [NII]6584
    "Z94",    #Hbeta,  [OII]3727, [OIII]5007, ([OIII]4959 )
    "M91",    #Hbeta,  [OII]3727, [OIII]5007, ([OIII]4959 )
    "C01_N2S2",    #[OII]3727, [OIII]5007, [NII]6584, [SII]6717
    "C01_R23",#Hbeta,  [OII]3727, [OIII]5007, ([OIII]4959 )

    "P05",   #Hbeta,  [OII]3727, [OIII]5007, ([OIII]4959 )
    "P01",   #available but deprecated
    "PP04_N2Ha",#Halpha, [NII]6584
    "PP04_O3N2",   #Halpha, Hbeta,[OIII]5007, [NII]6584
    "DP00",   #S23 available but deprecated
    "P10_ONS","P10_ON","M13_O3N2","M13_N2",
    "M08_R23","M08_N2Ha","M08_O3Hb","M08_O2Hb","M08_O3O2","M08_O3N2",
    "D13_N2S2_O3S2",    "D13_N2S2_O3Hb",    "D13_N2S2_O3O2",    "D13_N2O2_O3S2",    "D13_N2O2_O3Hb",    "D13_N2O2_O3O2",    "D13_N2Ha_O3Hb",    "D13_N2Ha_O3O2",

    "KD02_N2O2",   #Halpha, Hbeta,  [OII]3727, [NII]6584
    "KD02_N2S2",
    "KK04_N2Ha",   #Halpha, Hbeta,  [OII]3727, [NII]6584
    "KK04_R23", #Hbeta,  [OII]3727, [OIII]5007, ([OIII]4959 )
    "KD02comb_update","KK04comb"]
#'KD02_N2O2', 'KD03new_R23', 'M91', 'KD03_N2Ha'


def get_keys():
    return Zs
 

##############################################################################
##fz_roots function as used in the IDL code  FED:reference the code here!
##############################################################################

#@profile
def calculation(diags,measured,num,(bsmeas,bserr),mds,nps,outfilename='blah.txt',dust_corr=True,disp=False,saveres=False, verbose=False, scales=None): 

    global IGNOREDUST
    diags.setdustcorrect()
    raw_lines={}
    raw_lines['[OIII]5007']=np.array([float('NaN')])
    raw_lines['Hb']=np.array([float('NaN')])
    for k in measured.iterkeys():
        #kills all non-finite terms 
        measured[k][~(np.isfinite(measured[k][:]))]=0.0 
        raw_lines[k]=measured[k]



    ######FED why this????????????
    raw_lines['[OIII]4959']=raw_lines['[OIII]5007']/3.
    raw_lines['[OIII]49595007']=raw_lines['[OIII]4959']+raw_lines['[OIII]5007']
    diags.setHab(raw_lines['Ha'],raw_lines['Hb'])
    
    print nps
    #if Ha or Hb is zero, cannot do red correction
    if dust_corr and diags.hasHa and diags.hasHb:
        with np.errstate(invalid='ignore'):
            diags.calcEB_V()
    elif dust_corr and not IGNOREDUST:
        if nps>1:
            response=raw_input("WARNING: reddening correction cannot be done without both H_alpha and H_beta measurement!! Continuing without reddening correction? [Y/n]\n").lower()
            assert(not (response.startswith('n'))),"please fix the input file to include Halpha and Hbeta measurements"
        else: print "WARNING: reddening correction cannot be done without both H_alpha and H_beta measurement!!"

        IGNOREDUST=True
        dust_corr=False
        diags.mds['E(B-V)']=np.ones(len(raw_lines['Ha']))*1e-5
    else:
        diags.unsetdustcorrect()
        diags.mds['E(B-V)']=np.ones(len(raw_lines['Ha']))*1e-5

    for k in ['[OII]3727','[OIII]5007','[OI]6300','[OIII]4959',
              '[SII]6717','[SII]6731','[SIII]9069','[SIII]9532'
              '[OII]3727','[OIII]5007','[OI]6300','[OIII]4959',
              '[NII]6584','[SIII]9532']:
        if k not in raw_lines or len(raw_lines[k])==1:
            raw_lines[k]=np.array([0.]*num)

    diags.setOlines(raw_lines['[OII]3727'], raw_lines['[OIII]5007'], raw_lines['[OI]6300'], raw_lines['[OIII]4959'])

    diags.setNII(raw_lines['[NII]6584'])
    diags.setSII(raw_lines['[SII]6717'],raw_lines['[SII]6731'],raw_lines['[SIII]9069'],raw_lines['[SIII]9532'])

    if diags.checkminimumreq(dust_corr,IGNOREDUST) == -1:
        return -1
        
    diags.calcNIIOII()
    diags.calcNIISII()


    diags.calcR23()
    #diags.calcS23()

    diags.initialguess()
    mds=mds.split(',')
    #needs N2 and Ha
    #diags.printme()
    if verbose: print "calculating diagnostics: ",mds
    if 'all' in mds:
         diags.calcD02()
         if   os.getenv("PYQZ_DIR"):
              cmd_folder = os.getenv("PYQZ_DIR")+'/'
              if cmd_folder not in sys.path:
                   sys.path.insert(0, cmd_folder)
              import pyqz
              diags.calcpyqz()
         else:
              print '''set path to pyqz as environmental variable :
export PYQZ_DIR="your/path/where/pyqz/resides/ in bash, for example, if you want this diagnostic. '''


         diags.calcZ94()
         diags.calcM91()

         diags.calcPP04()

         diags.calcP05()
         diags.calcP10()

         diags.calcM08()
         diags.calcM13()

         diags.calcKD02_N2O2()
         diags.calcKK04_N2Ha()
         
         diags.calcKK04R23()
         diags.calcKDcombined()

    if 'DP00' in mds:
        diags.calcDP00()
    if 'P01' in mds:
        diags.calcP01()
         
    if 'D02' in mds:
         diags.calcD02()
    if 'D13' in mds:
         if   os.getenv("PYQZ_DIR"):
              cmd_folder = os.getenv("PYQZ_DIR")
              if cmd_folder not in sys.path:
                   sys.path.insert(0, cmd_folder)
              import pyqz
              #diags.calcpyqz()
              #in order to see the original pyqz plots
              #call pyqz with option plot=True by
              #using the commented line below instead
              diags.calcpyqz(plot=True)
         else:
              print '''set path to pyqz as environmental variable 
PYQZ_DIR if you want this diagnostic. '''

    if 'PP04' in mds:
        diags.calcPP04()
    if 'Z94' in mds:
        diags.calcZ94()
    if 'M91' in mds:
        diags.calcM91()
    if 'P10' in mds:
        diags.calcP10()
    if 'M13' in mds:
        diags.calcM13()
    if 'M08' in mds:
        diags.calcM08()
    if 'P05' in mds:
        diags.calcP05()
    if 'C01' in mds:
        diags.calcC01_ZR23()
    if 'KD02' in mds :
        diags.calcKD02_N2O2()
        diags.calcKK04_N2Ha()       
        diags.calcKK04R23()
        diags.calcKDcombined()
