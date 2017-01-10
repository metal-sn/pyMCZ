##############################################################################
## Calculates oxygen abundance (here called metalicity) based on strong emission lines, 
## based on code originally written in IDL by Lisa Kewley (Kewley & Ellison 2008). Outputs 
## oxygen abundance in many different diagnostics (see Bianco et al. 2016).
##
##new calculation based on the most recent version of the .pro file.
##
##inputs:
## measured - flux data, must be the format returned by readfile()
## num - number of spectra for which to calculate metallicity, also returned by readfile()
## outfilename - the name of the file the results will be appended to
## red_corr - reddening correction flag - True by default
## disp - if True prints the results, default False
##############################################################################

from __future__ import print_function
import sys
import os
import numpy as np

IGNOREDUST = False
MP = True
FIXNEGATIVES = True  # set to true if no negative flux measurements should be allowed. all negative flux measurements are set to 0

##list of metallicity methods, in order calculated
Zs = ["E(B-V)",  # based on Halpha, Hbeta
    "logR23",  # Hbeta,  [OII]3727, [OIII]5007, [OIII]4959

    "D02",  # Halpha, [NII]6584
    "Z94",  # Hbeta,  [OII]3727, [OIII]5007, ([OIII]4959 )
    "M91",  # Hbeta,  [OII]3727, [OIII]5007, ([OIII]4959 )
    "C01_N2S2",  # [OII]3727, [OIII]5007, [NII]6584, [SII]6717
    "C01_R23",  # Hbeta,  [OII]3727, [OIII]5007, ([OIII]4959 )

    "P05",  # Hbeta,  [OII]3727, [OIII]5007, ([OIII]4959 )
    "P01",  # available but deprecated
    "PP04_N2Ha",  # Halpha, [NII]6584
    "PP04_O3N2",  # Halpha, Hbeta,[OIII]5007, [NII]6584
    "DP00",  # S23 available but deprecated
    "P10_ONS", "P10_ON",
    "M08_R23", "M08_N2Ha", "M08_O3Hb", "M08_O2Hb", "M08_O3O2", "M08_O3N2",
    "M13_O3N2", "M13_N2",
    "D13_N2S2_O3S2", "D13_N2S2_O3Hb",
    "D13_N2S2_O3O2", "D13_N2O2_O3S2",
    "D13_N2O2_O3Hb", "D13_N2O2_O3O2",
    "D13_N2Ha_O3Hb", "D13_N2Ha_O3O2",
    "KD02_N2O2",  # Halpha, Hbeta,  [OII]3727, [NII]6584
    "KD02_N2S2",
    "KK04_N2Ha",  # Halpha, Hbeta,  [OII]3727, [NII]6584
    "KK04_R23",  # Hbeta,  [OII]3727, [OIII]5007, ([OIII]4959 )
    "KD02comb",
    "PM14",
    "D16"]  # ,"KK04comb"]
#'KD02_N2O2', 'KD03new_R23', 'M91', 'KD03_N2Ha'

Zserr = ['PM14err']  # ,"KK04comb"]
#'KD02_N2O2', 'KD03new_R23', 'M91', 'KD03_N2Ha'


def get_keys():
    return Zs


def get_errkeys():
    return Zserr


def printsafemulti(string, logf, nps):
    #this is needed because dealing with a log output with multiprocessing
    #is painful. but it introduces a bunch of if checks.
    if nps == 1:
        logf.write(string + "\n")
        #print >> logf, string
    else:
        print (string)


##############################################################################
##fz_roots function as used in the IDL code 
##############################################################################

#@profile
def calculation(mscales, measured, num, mds, nps, logf, dust_corr=True, disp=False, verbose=False):

    global IGNOREDUST
    mscales.setdustcorrect()
    raw_lines = {}
    raw_lines['[OIII]5007'] = np.array([float('NaN')])
    raw_lines['Hb'] = np.array([float('NaN')])
    raw_lines['Hz'] = np.array([float('NaN')])
    for k in measured.iterkeys():
        #kills all non-finite terms
        measured[k][~(np.isfinite(measured[k][:]))] = 0.0
        if FIXNEGATIVES:
            measured[k][measured[k] < 0] = 0.0
        raw_lines[k] = measured[k]

    ######we trust QM better than we trust the measurement of the [OIII]4959
    ######which is typical low S/N so we set it to [OIII]5007/3.
    ######change this only if youre spectra are very high SNR
    raw_lines['[OIII]4959'] = raw_lines['[OIII]5007'] / 3.
    raw_lines['[OIII]49595007'] = raw_lines['[OIII]4959'] + raw_lines['[OIII]5007']
    mscales.setHab(raw_lines['Ha'], raw_lines['Hb'])

    #if Ha or Hb is zero, cannot do red correction
    if dust_corr and mscales.hasHa and mscales.hasHb:
        with np.errstate(invalid='ignore'):
            mscales.calcEB_V()
    elif dust_corr and not IGNOREDUST:

        if nps > 1:
            print ('''WARNING: reddening correction cannot be done 
            without both H_alpha and H_beta measurement!!''')

        else:
            response = raw_input('''WARNING: reddening correction cannot be done without both H_alpha and H_beta measurement!! 
            Continuing without reddening correction? [Y/n]\n''').lower()
            assert(not (response.startswith('n'))), "please fix the input file to include Ha and Hb measurements"

        IGNOREDUST = True
        dust_corr = False
        mscales.mds['E(B-V)'] = np.ones(len(raw_lines['Ha'])) * 1e-5
    else:
        mscales.unsetdustcorrect()
        mscales.mds['E(B-V)'] = np.ones(len(raw_lines['Ha'])) * 1e-5

    for k in ['[OII]3727', '[OIII]5007', '[OI]6300', '[OIII]4959',
              '[SII]6717', '[SII]6731', '[SIII]9069', '[SIII]9532'
              '[OII]3727', '[OIII]5007', '[OI]6300', '[OIII]4959',
              '[NII]6584', '[SIII]9532']:
        if k not in raw_lines or (len(raw_lines[k]) == 1 and np.isnan(raw_lines[k][0])):
            raw_lines[k] = np.array([0.] * num)

    mscales.setOlines(raw_lines['[OII]3727'], raw_lines['[OIII]5007'], raw_lines['[OI]6300'], raw_lines['[OIII]4959'])
    mscales.setSII(raw_lines['[SII]6717'], raw_lines['[SII]6731'], raw_lines['[SIII]9069'], raw_lines['[SIII]9532'])
    mscales.setNII(raw_lines['[NII]6584'])

    if mscales.checkminimumreq(dust_corr, IGNOREDUST) == -1:
        return -1

    mscales.calcNIIOII()
    mscales.calcNIISII()

    mscales.calcR23()
    #mscales.calcS23()

    mscales.initialguess()
    mds = mds.split(',')


    #mscales.printme()
    if verbose:
        print ("calculating metallicity diagnostic scales: ", mds)
    if 'all' in mds:
        mscales.calcD02()
        if   os.getenv("PYQZ_DIR"):
            cmd_folder = os.getenv("PYQZ_DIR") + '/'
            if cmd_folder not in sys.path:
                sys.path.insert(0, cmd_folder)
            mscales.calcpyqz()
        else:
            printsafemulti('''WARNING: CANNOT CALCULATE pyqz: 
            set path to pyqz as environmental variable :
            export PYQZ_DIR="your/path/where/pyqz/resides/ in bash, for example, if you want this scale. \n''', logf, nps)

        mscales.calcZ94()
        mscales.calcM91()

        mscales.calcPP04()

        #mscales.calcP05()
        mscales.calcP10()

        mscales.calcM08()
        mscales.calcM13()

        mscales.calcKD02_N2O2()
        mscales.calcKK04_N2Ha()

        mscales.calcKK04_R23()
        mscales.calcKDcombined()

    if 'DP00' in mds:
        mscales.calcDP00()
    if 'P01' in mds:
        mscales.calcP01()

    if 'D02' in mds:
        mscales.calcD02()
    if 'D13' in mds:
        if   os.getenv("PYQZ_DIR"):
            cmd_folder = os.getenv("PYQZ_DIR")
            if cmd_folder not in sys.path:
                sys.path.insert(0, cmd_folder)
            #mscales.calcpyqz()
            #in order to see the original pyqz plots
            #call pyqz with option plot=True by
            #using the commented line below instead
            mscales.calcpyqz(plot=disp)
        else:
            printsafemulti('''WARNING: CANNOT CALCULATE pyqz: 
            set path to pyqz as environmental variable 
            PYQZ_DIR if you want this scale. ''', logf, nps)

    if 'D13all' in mds:
        if   os.getenv("PYQZ_DIR"):
            cmd_folder = os.getenv("PYQZ_DIR")
            if cmd_folder not in sys.path:
                sys.path.insert(0, cmd_folder)
            #mscales.calcpyqz()
            #in order to see the original pyqz plots
            #call pyqz with option plot=True by
            #using the commented line below instead
            mscales.calcpyqz(plot=disp, allD13=True)
        else:
            printsafemulti('''set path to pyqz as environmental variable 
PYQZ_DIR if you want this scale. ''', logf, nps)

    if 'PM14' in mds:
        if   os.getenv("HIICHI_DIR"):
            cmd_folder = os.getenv("HIICHI_DIR") + '/'
            if cmd_folder not in sys.path:
                sys.path.insert(0, cmd_folder)
            mscales.calcPM14()
    if 'PP04' in mds:
        mscales.calcPP04()
    if 'Z94' in mds:
        mscales.calcZ94()
    if 'M91' in mds:
        mscales.calcM91()
    if 'P10' in mds:
        mscales.calcP10()
    if 'M13' in mds:
        mscales.calcM13()
    if 'M08all' in mds:
        mscales.calcM08(allM08=True)
    elif 'M08' in mds:
        mscales.calcM08()
    if 'P05' in mds:
        mscales.calcP05()
    if 'C01' in mds:
        mscales.calcC01_ZR23()
    if 'KD02' in mds:
        mscales.calcKD02_N2O2()
        mscales.calcKK04_N2Ha()
        mscales.calcKK04_R23()
        mscales.calcKDcombined()
    if 'D16' in mds:
        mscales.calcD16()
