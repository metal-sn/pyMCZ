import os,sys,argparse,warnings

import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter

#modules of this package
import pylabsetup
import metallicity

#pickle may not be installed
NOPICKLE=False
try:
    import pprint, pickle
except:
    NOPICKLE=True


CLOBBER=False
VERBOSE=False
UNPICKLE=False
RUNSIM=True
BINMODE='t'
binning={'bb':'Bayesian blocks','d':"Doane's formula",'s':r'$\sqrt{N}$','t':r'$2 N^{1/3}$'}

##############################################################################
##Reads the flux file and returns it as an array.
##Ignores non-numeric lines
##Returns  (flux array,num)
##############################################################################
def readfile(filename):
    a=open(filename,'r')
    b=np.array([])
    i=0
    j=0
    
    for line in a:
        err_count=0
        strs=np.array(line.rstrip('\n').split())
        i=strs.size
        for word in range(i):
            try:
                strs[word]=strs[word].astype(np.float)
            except ValueError:
                err_count+=1
                strs[word]=-1
        strs=strs.astype(np.float)
        if err_count<i/2-1:    
            b=np.append(b,strs,axis=0)
            j+=1
    b.resize(j,i)
    b=np.transpose(b)
    a.close()

    return b,j

##############################################################################
##returns appropriate bin size for the number of data
##mode 'd' calculates this based on Doane's formula 
##mode 's' calculates this based on sqrt of number of data
##mode 't' calculates this based on 2*n**1/3 (default)
##############################################################################
def getbinsize(n,data,):
    if BINMODE=='d':
        g1=stats.mstats.moment(data,moment=3)
        s1=np.sqrt(float(n)/6.0)
        #s1=1.0/np.sqrt(6.*(n-2.)/((n+1.)*(n+3.)))
        k=1+np.log(n)+np.log(1+(g1*s1))
    elif BINMODE=='s':
        k=np.sqrt(n)
    elif BINMODE=='t':
        k=2.*n**(1./3.)
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
##Save the result as histogram as name
## delog - if true de-logs the data. False by default
##############################################################################
def savehist(data,snname,Zs,nsample,i,path,nmeas,delog=False):
    
    global CLOBBER
    
    name='%s_n%d_%s_%d'%((snname,nsample,Zs,i+1))
    outdir=os.path.join(path,'hist')
    outfile=os.path.join(outdir,name+".pdf")
    if os.path.isfile(outfile) and not CLOBBER:
        replace=raw_input("replacing existing image files, starting with: %s ? [Y/n]\n"%outfile).lower()
        assert(not (replace.startswith('n'))),"save your existing output directory under another name first"
        CLOBBER =True

    plt.clf()
    ###de-log###
    if delog:
        with np.errstate(invalid='ignore'):
            data=np.power(10,np.absolute(data-12))
        
    ####kill outliers###
    data=data[np.isfinite(data)]
    data,ignore,ignore=stats.sigmaclip(data,high=5.0,low=5.0)
    n=data.shape[0]
    
    if data.shape[0]<=0 or np.sum(data)<=0:
        print name,'is blank' ##if no data
        return "-1, -1, -1"    
    try:
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

        ######histogram######
        ###find appropriate bin size###
        ##if astroML is available use it to get Bayesian blocks
        if BINMODE=='bb':
            try:
                from astroML.plotting import hist as amlhist
                distrib=amlhist(data, bins='blocks', normed=True)
                plt.clf()
            except:
                print "bayesian blocks for histogram requires astroML to be installed"
                print "defaulting to 2*n**1/3 "
                ##otherwise 
                numbin=getbinsize(data.shape[0],data)        
                distrib=np.histogram(data, numbin, density=True)
        else:
            numbin=getbinsize(data.shape[0],data)        
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
        plt.annotate(st, xy=(0.13, 0.6), xycoords='axes fraction',fontsize=18,fontweight='bold')
        st='%s '%(Zs)
        plt.annotate(st, xy=(0.62, 0.93), xycoords='axes fraction',fontsize=18,fontweight='bold')
        st='measurement %d of %d\n\nmedian: %.3f\n16th Percentile: %.3f\n84th Percentile: %.3f'%(i+1,nmeas,round(median,3),round(left,3),round(right,3))
        plt.annotate(st, xy=(0.62, 0.65), xycoords='axes fraction',fontsize=18)
        st='MC sample size %d\nhistogram rule: %s'%(nsample,binning[BINMODE])
        plt.annotate(st, xy=(0.62, 0.55), xycoords='axes fraction',fontsize=13)
        if delog:
            plt.xlabel('O/H')
        elif Zs == "E(B-V)":
            plt.xlabel('E(B-V) [mag]')
        else:
            plt.xlabel('12+log(O/H)')
        plt.ylabel('relative counts')
        plt.savefig(outfile,format='pdf')
        
        ###print out the confidence interval###
        print '{0:15} {1:20} {2:>13.3f} - {3:>7.3f} + {4:>7.3f}'.format(snname, Zs, round(median,3), round(median-left,3), round(right-median,3))

        return "%f\t %f\t %f"%(round(median,3), round(median-left,3), round(right-median,3))

    except (OverflowError,AttributeError,ValueError):
        if VERBOSE: print data
        print name, 'had infinities'
        return "-2, -2"

##############################################################################
## The main function. takes the flux and its error as input. 
##  filename - a string 'filename' common to the three flux files
##  flux - np array of the fluxes
##  err - the flux errors, must be the same dimension as flux
##  nsample - the number of samples the code will generate. Default is 100
##  errmode - determines which method to choose the bin size.
##      mode 'd' calculates this based on Doane's formula
##      mode 's' calculates this based on sqrt of number of data
##      mode 't' calculates this based on 2*n**1/3 (default)
##############################################################################
def run((name, flux, err, path), nsample,delog=False, unpickle=False):
    global RUNSIM
    assert(flux.shape == err.shape), "flux and err must be same dimensions" 

    newnsample=int(nsample+0.1*nsample)
    p=os.path.join(path,'..')
    nm = flux.shape[1]

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
            res[key]=[]
            
        #do the iterations
        for i in range(newnsample):
            temp=flux+err*sample[i]
            warnings.filterwarnings("ignore")
            t=metallicity.calculation(temp,nm,disp=VERBOSE)
            for key in Zs:
                res[key].append(t[key])
            
        #recast the result into np.array
        for key in Zs:
            res[key]=np.array(res[key])
        
        if VERBOSE: print "Iteration Complete"
    
        #"I CAN PICKLE THIS!"
        #pickle this realization
        if not NOPICKLE:
            pickle.dump(res,open(picklefile,'wb'))

    ###Bin the results and save###
    print '{0:15} {1:20} {2:>13} - {3:>7} + {4:>7} {5:11} {6:>7}'.format("SN","diagnostic", "metallicity","34%", "34%", "(sample size:",'%d)'%nsample)
    for i in range(nm):
        fi=open(os.path.join(binp,'%s_n%d_%d.txt'%(name,nsample,i+1)),'w')
        fi.write("%s\t Median Oxygen abundance (12+log(O/H))\t 16th percentile\t 84th percentile\n"%name)
        
        print "\n\nmeasurement %d-------------------------------------------------------------"%(i+1)
        
        for key in Zs:
            s=key+"\t "+savehist(res[key][:,i],name,key,nsample,i,binp,nm,delog=delog)+'\n'
            fi.write(s)
        fi.close()
    
        if VERBOSE: print "uncertainty calculation complete"

##############################################################################
##The input format generator
##############################################################################
def input_format(filename,path):
    p = os.path.join(path,"sn_data") 
    assert os.path.isdir(p), "bad data directory %s"%p
    if os.path.isfile(os.path.join(p,filename+'_max.txt')) and os.path.isfile(os.path.join(p,filename+'_min.txt')):
        if os.path.isfile(os.path.join(p,filename+'_med.txt')):
            return in_mmm(filename,path=p)
        return in_mm(filename,path=p)
    
    print "Unable to find _min and _max files ",filename+'_max.txt',filename+'_min.txt',"in directory ",p
    return -1

def in_mmm(filename,path):
#    p=os.path.abspath('..')
#    p+='\\sn_data\\'
    
    ###Initialize###
    maxfile=os.path.join(path,filename+"_max.txt")
    medfile=os.path.join(path,filename+"_med.txt")
    minfile=os.path.join(path,filename+"_min.txt")
    
    ###read the max, med, min flux files###    
    maxf,nm=readfile(maxfile)
    medf,num=readfile(medfile)
    minf,nn=readfile(minfile)

    ###calculate the flux error as 1/2 [(max-med)+(med-min)]###
    err=0.5*((maxf-medf)+(medf-minf))
    return (filename, medf, err, path)

def in_mm(filename, path):

    ###Initialize###
    maxfile=os.path.join(path,filename+"_max.txt")
    minfile=os.path.join(path,filename+"_min.txt")
    
    ###read the max, med, min flux files###    
    maxf,nm=readfile(maxfile)
    minf,nn=readfile(minfile)
    ###calculate the flux error as 1/2 [(max-min)]###
    err=0.5*(maxf - minf)
    medf= minf+err
    return (filename, medf, err, path)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('name', metavar='<name>', type=str, help="the SN file name (root of the _min,_max file names")
    parser.add_argument('nsample', metavar='N', type=int, help="number of iterations, minimum 100")
    parser.add_argument('--clobber',default=False, action='store_true', help="replace existing output")
    parser.add_argument('--delog',default=False, action='store_true', help="result in natural, not log space. default is log space")
    parser.add_argument('--binmode', default='t', type=str, choices=['d','s','t','bb'], help='method to determine bin size {d: Duanes formula, s: n^1/2, t: 2*n**1/3(default), Bayesian blocks}')
    parser.add_argument('--path',   default=None, type=str, help="input/output path (must contain the input _max.txt and _min.txt files in a subdirectory sn_data)")
    parser.add_argument('--unpickle',   default=False, action='store_true', help="read the pickled realization instead of making a new one")

    parser.add_argument('--verbose',default=False, action='store_true', help="verbose mode")

    args=parser.parse_args()

    global CLOBBER
    global VERBOSE
    global BINMODE
    CLOBBER=args.clobber
    VERBOSE=args.verbose
    BINMODE=args.binmode

    if args.unpickle and NOPICKLE:
        args.unpickle = False
        print "cannot use pickle on this machine, wont save and wont read saved realizations. Ctr-C to exit, Return to continue?\n"
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
            run(fi,args.nsample,delog=args.delog, unpickle=args.unpickle)
    else:
        print "nsample must be at least 100"
    
if __name__ == "__main__":
    main()
    #files=['sn2006ss','ptf10eqi-z']
    #filename=files[1]
    #nsample=10000
    

