import os,sys,argparse,warnings

import numpy as np
import scipy.stats as stats
from scipy.special import gammaln
from scipy import optimize
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter


#modules of this package
import pylabsetup

#import metallicity_save2 as metallicity
import fedmetallicity as metallicity

PROFILING = True
PROFILING = False


alllines=['[OII]3727','Hb','[OIII]4959','[OIII]5007','[OI]6300','Ha','[NII]6584','[SII]6717','[SII]6731','[SIII]9069','[SIII]9532']
morelines=['E(B-V)','dE(B-V)','scale_blue','d scale_blue']



OLD=False

#pickle may not be installed
NOPICKLE=False
try:
    import pprint, pickle
except:
    NOPICKLE=True


CLOBBER=False
VERBOSE=False
UNPICKLE=False
ASCIIOUTPUT=False
RUNSIM=True
BINMODE='k'
binning={'bb':'Bayesian blocks','k':"Knuth's rule",'d':"Doane's formula",'s':r'$\sqrt{N}$','t':r'$2 N^{1/3}$', 'kd':'Kernel Density'}

def is_number(s):
    if not type(s) is np.string_:
        try:
            tmp=float(s)
            return True
        except :
            return False
    return False
        

def getknuth(m,data,N):
    m=int(m)
    bins=np.linspace(min(data),max(data), int(m) + 1)
    nk,bins=np.histogram(data,bins)
    return -(N*np.log(m) + gammaln(0.5*m) - m*gammaln(0.5) -
             gammaln(N+0.5*m)+np.sum(gammaln(nk+0.5)))

def knuthn(data, maxM=None):
    assert data.ndim==1, "data must be 1D array to calculate Knuth's number of bins"
    N=data.size
    if not maxM:
        maxM=5*np.sqrt(N)
    m0=2.0*N**(1./3.)
    mkall= optimize.fmin(getknuth,m0, args=(data,N), disp=VERBOSE, maxiter=30)#[0]
    mk=mkall[0]
    if mk>maxM or mk<0.3*np.sqrt(N):
        mk=m0
        return mk, 't'
    return mk, 0
    '''

    m0=2.0*(N**(1.0/3.0))
    mkall= optimize.fmin(getknuth,m0, args=(data,N), disp=VERBOSE, maxiter=30)
    mk=mkall[0]
    if mk>maxM:
        mk=m0
    return mk 
    '''

##############################################################################
##Reads the flux file and returns it as an array.
##Ignores non-numeric lines
##Returns  (flux array,num)
##############################################################################
def readfile(filename):
    noheader=1
    findex=-1
    f=open(filename,'r')
    l0=f.readline()
    l1=f.readline().split()
    if l0.startswith('#') or l0.startswith(';'):
        header=l0.strip().replace(";",'').replace("#",'').split(',');
        header[0]=header[0].replace(' ','')
        header=header[:len(l1)]
    else:
        noheader=0
        header=['galnum']+alllines+['flag']+morelines
        header=header[:len(l1)]


    #print header
    formats=['i']+['f']*(len(header)-1)
    if 'flag' in header:
        findex=header.index('flag')
        formats[findex]='S10'
    #cols=tuple([i for i in range(len(header)) if not i==findex])
    #print formats
    
    bstruct={}
    for i,k in enumerate(header):
        bstruct[k]=[i,0]
    print "file header",header
    b = np.loadtxt(filename,skiprows=noheader, dtype={'names':header,'formats':formats})
    #usecols=cols, unpack=True)
    
    for i,k in enumerate(header):
        if not k=='flag' and is_number(b[k][0]):
            bstruct[k][1]=np.count_nonzero(b[k])+sum(np.isnan(b[k]))
    j=len(b['galnum'])
    return b,j,bstruct

##############################################################################
##The input format generator
##############################################################################
def input_format(filename,path):
    p = os.path.join(path,"input") 
    assert os.path.isdir(p), "bad data directory %s"%p
    if os.path.isfile(os.path.join(p,filename+'_err.txt')):
        if os.path.isfile(os.path.join(p,filename+'_meas.txt')):
            return ingest_data(filename,path=p)            
    print "Unable to find _meas and _err files ",filename+'_meas.txt',filename+'_err.txt',"in directory ",p
    return -1

def ingest_data(filename,path):
    ###Initialize###
    measfile=os.path.join(path,filename+"_meas.txt")
    errfile =os.path.join(path,filename+"_err.txt")
    
    ###read the max, meas, min flux files###    
    meas,nm, bsmeas=readfile(measfile)
    err, nn, bserr =readfile(errfile)
    return (filename, meas, err, nm, path, (bsmeas,bserr))

##############################################################################
##sets which metallicity scales can be calculated based on the available lines
#############################################################################

def setscales(bss):
    Zs= metallicity.get_keys()
    scales={}
    #set all to true for now
    for s in Zs:
        scales[s]=True
    return scales

##############################################################################
##returns appropriate bin size for the number of data
##mode 'k' calculates this based on Knuth's rule
##mode 'd' calculates this based on Doane's formula 
##mode 's' calculates this based on sqrt of number of data
##mode 't' calculates this based on 2*n**1/3 (default)
##############################################################################
def getbinsize(n,data,):
    if BINMODE=='d':
        g1=stats.mstats.moment(data,moment=3)
        s1=np.sqrt(float(n)/6.0)
        #s1=1.0/np.sqrt(6.*(n-2.)/((n+1.)*(n+3.)))
        k=1+np.log(n)+np.log(1+(g1*s1)),0
    elif BINMODE=='s':
        k=np.sqrt(n),0
    elif BINMODE=='t':
        k=2.*n**(1./3.),0
    else:
        #from astroML.plotting import hist as amlhist
        #distrib=amlhist(data, bins='knuth', normed=True)
        k= knuthn(data)
        #distrib=amlhist(data, bins='knuth', normed=True)
    return k

##############################################################################        
##estimating error starting at the peak.
##(almost) symmetric - make it go same l and r
##not well behaved in case of multiple peaks
##DEPRECATED
##############################################################################        
def err_est1(count,prob=0.68):
    peak=np.argmax(count)
    total=np.sum(count)
    temp=0
    l=0
    r=0
    while temp<total*prob:
        if count[peak-l]>=count[peak+r]:
            temp+=count[peak-l]
            l+=1
        else:
            temp+=count[peak+r]
            r+=1

    return peak-l,peak-r,(total-temp)/total,0,0

##############################################################################
##Check if hist files need to be replaced
##############################################################################
def checkhist(snname,Zs,nsample,i,path):
    global CLOBBER
    
    name='%s_n%d_%s_%d'%((snname,nsample,Zs,i+1))
    outdir=os.path.join(path,'hist')
    outfile=os.path.join(outdir,name+".pdf")
    if os.path.isfile(outfile) and not CLOBBER:
        replace=raw_input("replacing existing image files, starting with: %s ? [Y/n]\n"%outfile).lower()
        assert(not (replace.startswith('n'))),"save your existing output directory under another name first"
        CLOBBER =True

##############################################################################
##Save the result as histogram as name
## delog - if true de-logs the data. False by default
##############################################################################
def savehist(data,snname,Zs,nsample,i,path,nmeas,delog=False, verbose=False, fs=18):
    global BINMODE
    name='%s_n%d_%s_%d'%((snname,nsample,Zs,i+1))
    outdir=os.path.join(path,'hist')
    outfile=os.path.join(outdir,name+".pdf")
    plt.clf()
    ###de-log###
    if delog:
        with np.errstate(invalid='ignore'):
            data=np.power(10,np.absolute(data-12))
        
    ####kill outliers###
    data=data[np.isfinite(data)]
    #if max(data)-min(data)>0.0001:
    #    data,ignore,ignore=stats.sigmaclip(data,high=5.0,low=5.0)
    n=data.shape[0]
    if not n>0:
        if verbose:print "data must be an actual distribution (n>0 elements!, %s)"%Zs
        return "-1,-1"
    #if not max(data)-min(data)>0.1:
    #    if verbose:print "the data must be in a distribution, not all the same!"
    #    return "-1,-1"
    if data.shape[0]<=0 or np.sum(data)<=0:
        print '{0:15} {1:20} {2:>13d}   {3:>7d}   {4:>7d} '.format(name.split('_')[0],Zs,-1,-1,-1)
        return "-1, -1, -1"    
    if 1:
#    try:
        ###find C.I.###
        median,pc16,pc84=np.percentile(data,[50,16,84])
        std=np.std(data)
        left=pc16
        right=pc84
        maxleft=median-std*5
        maxright=median+std*5
        if "%2f"%maxright=="%2f"%maxleft:
            maxleft=median-1
            maxright=median+1
        if round(right,5)==round(left,5) and round(left,5)==round(median,5):
            print '{0:15} {1:20} {2:>13.3f} - {3:>7.3f} + {4:>7.3f} (no distribution)'.format(name.split('_')[0],Zs,median,0,0 )

            return "-2,-2"
        ######histogram######
        ##if sklearn is available use it to get Kernel Density
        if BINMODE=='kd':
            try:
                from sklearn.neighbors import KernelDensity
            except:
                print '''sklearn is not available, 
                thus we cannot compute kernel density. 
                switching to bayesoan blocks'''
                BINMODE='bb'
        if BINMODE=='kd':
            bw=(data.max()-data.min())/4.
            if bw >0:
                kde = KernelDensity(kernel='tophat', bandwidth=bw).fit(data[:, np.newaxis])
                bins=np.linspace(maxleft,maxright,1000)[:, np.newaxis]
                log_dens = kde.score_samples(bins)
                dens=np.exp(log_dens)
                #print dens
                plt.fill(bins[:,0], dens/dens.max(), fc='#AAAAFF')
            numbin,bm=getbinsize(data.shape[0],data)        
            distrib=np.histogram(data, bins=numbin, density=True)            
            ###make hist###
            counts, bins=distrib[0],distrib[1]
            widths=np.diff(bins)
            countsnorm=counts/np.max(counts)
#            plt.bar(bins[:-1],countsnorm,widths,color=['gray'], alpha=0.3)
        ###find appropriate bin size###
        ##if astroML is available use it to get Bayesian blocks
        else:
            if BINMODE=='bb' :
                try:
                    from astroML.plotting import hist as amlhist
                    if BINMODE=='bb':
                        distrib=amlhist(data, bins='blocks', normed=True)
                    plt.clf()
                except:
                    print "bayesian blocks for histogram requires astroML to be installed"
                    print "defaulting to Knuth's rule "
                    ##otherwise 
                    numbin,bm=getbinsize(data.shape[0],data)        
                    distrib=np.histogram(data, numbin, density=True)
            else:
                numbin,bm=getbinsize(data.shape[0],data)        
                distrib=np.histogram(data, numbin, density=True)            
            ###make hist###
            counts, bins=distrib[0],distrib[1]
            widths=np.diff(bins)
            countsnorm=counts/np.max(counts)

        ###plot hist###
        plt.bar(bins[:-1],countsnorm,widths,color=['gray'])
        plt.minorticks_on()
        plt.gca().xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
        plt.xlim(maxleft,maxright)
        plt.ylim(0,1.15)
        plt.yticks(np.arange(0.2,1.1,0.2 ), [ "%.1f"%x for x in np.arange(0.2,1.1,0.2)])  
        plt.axvspan(left,right,color='DarkOrange',alpha=0.4)
        plt.axvline(x=median,linewidth=2,color='white',ls='--')
        st='%s '%(snname)
        plt.annotate(st, xy=(0.13, 0.6), xycoords='axes fraction',size=fs,fontweight='bold')
        st='%s '%(Zs.replace('_',' '))
        plt.annotate(st, xy=(0.62, 0.93), xycoords='axes fraction',fontsize=fs,fontweight='bold')
        st='measurement %d of %d\n\nmedian: %.3f\n16th Percentile: %.3f\n84th Percentile: %.3f'%(i+1,nmeas,round(median,3),round(left,3),round(right,3))
        plt.annotate(st, xy=(0.62, 0.65), xycoords='axes fraction',fontsize=fs)
        st='MC sample size %d\nhistogram rule: %s'%(nsample,binning[BINMODE])   
        if bm:
            st='MC sample size %d\nhistogram rule: %s'%(nsample,binning[bm])
        plt.annotate(st, xy=(0.62, 0.55), xycoords='axes fraction',fontsize=fs-5)
        if delog:
            plt.xlabel('O/H')
        elif "E(B-V)" in Zs:
            plt.xlabel('E(B-V) [mag]')
        elif "logR23" in Zs:
            plt.xlabel('logR23')
        else:
            plt.xlabel('12+log(O/H)')
        plt.ylabel('relative counts')
        plt.savefig(outfile,format='pdf')
        ###print out the confidence interval###
        print '{0:15} {1:20} {2:>13.3f} - {3:>7.3f} + {4:>7.3f}'.format(snname, Zs, round(median,3), round(median-left,3), round(right-median,3))
        return "%f\t %f\t %f"%(round(median,3), round(median-left,3), round(right-median,3))

#    except (OverflowError,AttributeError,ValueError):
#        if VERBOSE: print data
#        print name, 'had infinities'
#        return "-2, -2"


##############################################################################
## The main function. takes the flux and its error as input. 
##  filename - a string 'filename' common to the three flux files
##  flux - np array of the fluxes
##  err - the flux errors, must be the same dimension as flux
##  nsample - the number of samples the code will generate. Default is 100
##  errmode - determines which method to choose the bin size.
##      mode 'k' calculates this based on Knuth's rule
##      mode 'd' calculates this based on Doane's formula
##      mode 's' calculates this based on sqrt of number of data
##      mode 't' calculates this based on 2*n**1/3 (default)
##############################################################################
def run((name, flux, err, nm, path, bss), nsample,smass,mds,delog=False, unpickle=False, dust_corr=True, verbose=False, fs=18):
    global RUNSIM,BINMODE
    assert(len(flux[0])== len(err[0])), "flux and err must be same dimensions" 
    assert(len(flux['galnum'])== nm), "flux and err must be of declaired size" 
    assert(len(err['galnum'])== nm), "flux and err must be same dimensions" 
    
    newnsample=int(nsample+0.1*nsample)
    p=os.path.join(path,'..')

    ###retrieve the metallicity keys
    Zs= metallicity.get_keys()

    ###make necessary paths
    if not os.path.exists(os.path.join(p,'output','%s'%name)):
        os.makedirs(os.path.join(p,'output','%s'%name))
    if not os.path.exists(os.path.join(p,'output','%s'%name,'hist')):
        os.makedirs(os.path.join(p,'output','%s'%name,'hist'))
    binp=os.path.join(p,'output','%s'%name)
    picklefile=os.path.join(binp,'%s_n%d.pkl'%(name,nsample))

    if VERBOSE: print "output files will be stored in ",binp

    if not CLOBBER:
        for key in Zs:
            for i in range(nm):
                checkhist(name,key,nsample,i,binp)

    if unpickle:
        RUNSIM=False
        if not os.path.isfile(picklefile):
            print "missing pickled file for this simulation: name, nsample.\nrun the MonteCarlo? Ctr-C to exit, Return to continue?\n"
            raw_input()
            RUNSIM=True
        else:
            pklfile = open(picklefile, 'rb')
            res=pickle.load(pklfile)

    if RUNSIM:
        ###Sample 'nsample' points from a gaussian centered on 0 with std 1
        mu=0
        sigma=1
        sample=np.random.normal(mu,sigma,newnsample)
        
        ###Start calculation###
        ## the flux to be feed to the calculation will be
        ## flux + error*i
        ## where i is the sampled gaussian    
        if VERBOSE: print "Starting iteration"
        
        #initialize the dictionary
        res={}
        for key in Zs:
            res[key]=[[] for i in range(nm)]

        #do the iterations
        temp={}
        delkeys=[]
        for k in bss[0].iterkeys():
            if k=='flag' or k=='galnum' or bss[0][k][1]==0 :#or bss[1][k][1]==bss[0][k][1]:
                delkeys.append(k)
        for k in delkeys:
                del bss[0][k]
                del bss[1][k]

        scales=setscales(bss[0])
        import diagnostics as dd
        for i in range(nm):
            diags=dd.diagnostics(newnsample)
            print "\n\nreading in measurements ",i+1
            #for i in range(newnsample):
            fluxi={}#np.zeros((len(bss[0]),nm),float)
            for j,k in enumerate(bss[0].iterkeys()):
                print '{0:15} '.format(k),
                #print bss[0][k][1], bss[1][k][1]
                #print_options.set_float_precision(2)
                print '{0:0.2} +/- {1:0.2}'.format(flux[k][i],err[k][i])
                fluxi[k]=flux[k][i]*np.ones(len(sample))+err[k][i]*sample
                warnings.filterwarnings("ignore")
            success=metallicity.calculation(diags,fluxi,nm,bss,smass,mds,disp=VERBOSE, dust_corr=dust_corr,verbose=VERBOSE)
            if success==-1:
                print "MINIMUM REQUIRED LINES: '[OII]3727','[OIII]5007','[NII]6584','[SII]6717','Ha','Hb' and 6.0<Smass<14 MSun"

#            diags.printme()
#            s=key+"\t "+savehist(t,'test','EB_V',100,i,binp,nm,delog=delog)+'\n'
#            plt.hist(t)
#            plt.show()
#            for k in diags.mds.iterkeys():
#                if k in ['KD_comb_NEW']: print '*',
#                print '{0:20}'.format(k),
#                if not diags.mds[k]==None:
#                    print ' {0:4} {1:4}'.format( stats.nanmean(diags.mds[k]),stats.nanstd(diags.mds[k])),
#                print ""
            for key in diags.mds.iterkeys():
                res[key][i]=diags.mds[key]
                if res[key][i]==None:
                    res[key][i]=[float('NaN')]*len(sample)
        for key in diags.mds.iterkeys():
            res[key]=np.array(res[key]).T
        #recast the result into np.array
        ##        for key in Zs:
        ##           res[key]=np.array(res[key])
        
        if VERBOSE: print "Iteration Complete"
    
        #"I CAN PICKLE THIS!"
        #pickle this realization
        if not NOPICKLE:
            pickle.dump(res,open(picklefile,'wb'))
            

    from matplotlib.font_manager import findfont, FontProperties
    
    print findfont(FontProperties())
    if 'Time' not in  findfont(FontProperties()):
        fs=15
    ###Bin the results and save###
    print '{0:15} {1:20} {2:>13} - {3:>7} + {4:>7} {5:11} {6:>7}'.format("SN","diagnostic", "metallicity","34%", "34%", "(sample size:",'%d)'%nsample)
    #return -1
    for i in range(nm):
        if ASCIIOUTPUT:
            fi=open(os.path.join(binp,'%s_n%d_%d.txt'%(name,nsample,i+1)),'w')
            fi.write("%s\t Median Oxygen abundance (12+log(O/H))\t 16th percentile\t 84th percentile\n"%name)
        
        print "\n\nmeasurement %d-------------------------------------------------------------"%(i+1)
        for key in Zs:
            if len(res[key].shape)>1 and sum(sum(~np.isnan(res[key])))>0:
                s=key+"\t "+savehist(res[key][:,i],name,key,nsample,i,binp,nm,delog=delog, verbose=verbose, fs=fs)+'\n'
                if ASCIIOUTPUT:
                    fi.write(s)

        if ASCIIOUTPUT:
            fi.close()
        
    
        if VERBOSE: print "uncertainty calculation complete"

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('name', metavar='<name>', type=str, help="the SN file name (root of the _min,_max file names")
    parser.add_argument('nsample', metavar='N', type=int, help="number of iterations, minimum 100")
    parser.add_argument('--clobber',default=False, action='store_true', help="replace existing output")
    parser.add_argument('--delog',default=False, action='store_true', help="result in natural, not log space. default is log space")
    parser.add_argument('--binmode', default='k', type=str, choices=['d','s','k','t','bb','kd'], help="method to determine bin size {d: Duanes formula, s: n^1/2, t: 2*n**1/3(default), k: Knuth's rule, bb: Bayesian blocks, kd: Kernel Density}")
    parser.add_argument('--path',   default=None, type=str, help="input/output path (must contain the input _max.txt and _min.txt files in a subdirectory sn_data)")
    parser.add_argument('--unpickle',   default=False, action='store_true', help="read the pickled realization instead of making a new one")

    parser.add_argument('--verbose',default=False, action='store_true', help="verbose mode")
    parser.add_argument('--mass',default=10, type=float,help="stellar mass, which can be validated")
    parser.add_argument('--nodust',default=False, action='store_true', help=" dont do dust corrections (default is to do it)")
    parser.add_argument('--asciiout',default=False, action='store_true', help=" write distribution in an ascii output (default is not to)")
    parser.add_argument('--md',default='all', type =str, help=" metallivity diagnostic to calculate. default is 'all', options are: D02, Z94, M91,C01, Pi01, PP04, pyqz, KD02, KD02comb")
    args=parser.parse_args()

    global CLOBBER
    global VERBOSE
    global BINMODE
    global ASCIIOUTPUT
    CLOBBER=args.clobber
    VERBOSE=args.verbose
    BINMODE=args.binmode

    ASCIIOUTPUT=args.asciiout
    if args.unpickle and NOPICKLE:
        args.unpickle = False
        print "cannot use pickle on this machine, wont save and won't read saved realizations. Ctr-C to exit, Return to continue?\n"
        raw_input();

    if args.path:
        path=args.path
    else:
        assert (os.getenv("MCMetdata"))," the _max, _min (and _med) data must live in a folder named sn_data. pass a path to the sn_data folder, or set up the environmental variable MCMetdata pointing to the path where sn_data lives "
        path=os.getenv("MCMetdata")
    assert(os.path.isdir(path)),"pass a path or set up the environmental variable MCMetdata pointing to the path where the _min _max _med files live"

    if args.nsample>=100:
        fi=input_format(args.name, path=path)
        if fi!=-1:
            run(fi,args.nsample,args.mass,args.md, delog=args.delog, unpickle=args.unpickle, dust_corr=(not args.nodust), verbose=VERBOSE)
    else:
        print "nsample must be at least 100"
    

if __name__ == "__main__":
    if PROFILING:
        import cProfile
        cProfile.run("main()")
    else:
        main()
    #files=['sn2006ss','ptf10eqi-z']
    #filename=files[1]
    #nsample=10000
    

