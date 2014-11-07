import numpy as np
import matplotlib.pyplot as plt
import scipy.stats.mstats as ssm
import metallicity
import os
from matplotlib.ticker import FuncFormatter

import sys, argparse

CLOBBER=False
VERBOSE=False
BINMODE='t'

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
##a simple gaussian
##############################################################################
def gaussian(x,m,s):
    return 1./(s*np.sqrt(2*np.pi)) * np.exp(-(x-m)**2 /(2 *s**2))

##############################################################################
##returns appropriate bin size for the number of data
##mode 'd' calculates this based on Doane's formula 
##mode 's' calculates this based on sqrt of number of data
##mode 't' calculates this based on 2*n**1/3 (default)
##############################################################################
def getbinsize(n,data,):
    if BINMODE=='d':
        g1=ssm.moment(data,moment=3)
        s1=np.sqrt(6.*(n-2.)/((n+1.)*(n+3.)))
        k=1+np.log2(n)+np.log2(1+np.abs(g1/s1))
    elif BINMODE=='s':
        k=np.sqrt(n)
    elif BINMODE=='t':
        k=2.*n**(1./3.)
    return k
##############################################################################        
##almost (symmetric) method
##estimating error starting at the peak.
##almost (symmetric) - make it go same l and r
##not well behaved in case of multiple peaks
##code not up to date(oct 21 2013)
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
##(central) method
#estmating error starting from the ends.
##input: 'count' returned from plt.hist, 'prob' default 68%
##output: the index of left boundary, right boundary, confidence,
#fraction to be shifted to left, fraction to be shifted from the right
##############################################################################
def err_est(count,prob=0.68):
    total=np.sum(count)
    l=0
    fl=0.
    r=len(count)-1
    fr=0.
    ltemp=count[l]
    rtemp=count[r]
    
    thresh=total*(1.-prob)/2.
    while ltemp<=thresh:
        l+=1
        ltemp+=count[l]
        if ltemp>thresh:
            ltemp-=count[l]
            fl=(thresh-ltemp)/count[l]   
            l-=1
            break
    while rtemp<=thresh:
        r-=1
        rtemp+=count[r]
        if rtemp>thresh:
            rtemp-=count[r]
            fr=(thresh-rtemp)/count[r]
            r+=1
            break
    return l,r,(total-2*thresh)/total,fl,1-fr    

##############################################################################
##Save the result as histogram as name
## delog - if true de-logs the data. False by default
##############################################################################
def savehist(data,filename,Zs,nsample,i,path,delog=False):
    
    global CLOBBER
    name='%s_n%d_%s_i%d'%((filename,nsample,Zs,i))

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
        
    data=data[np.isfinite(data)]
    #    data=np.sort(data)
    
    ####kill outliers###

    mean=np.mean(data)
    std=np.std(data)
    data=data[(data>mean-5*std)*(data>mean-5*std)][:nsample]

    n=data.shape[0]
    
    if data.shape[0]<=0:
        print name,'is blank' ##if no data
        return "-1, -1, -1"    
    try:
        median=np.median(data)
        ######find fit and save hist######
        ###find appropriate bin size###
        numbin=getbinsize(data.shape[0],data)
        
        ###make hist###
        count, bins, ignored = plt.hist(data, numbin, normed=1.0,color=['steelblue'])

        #####FED
        to_unity = lambda y, pos:  "%.2f"%(y / float(max(count)))
        plt.gca().yaxis.set_major_formatter(FuncFormatter(to_unity))

        ###find error###
        l,r,t,fl,fr=err_est(count)
        y=np.zeros(len(bins))
        left=bins[l]+(bins[l+1]-bins[l])*fl
        right=bins[r]-(bins[r]-bins[r-1])*fr
        if right<median or left>median:
            median = (left+right)/2
        ###plot hist###
        plt.plot(bins,y)
        plt.xlim(8.6,9.3)
        plt.axvspan(left,right,color='red',alpha=0.4)
        st='%s i=%d\nn=%d\nconfidence: %.2f\nmedian: %.3f\n16th Percentile: %.3f\n84th Percentile: %.3f'%(Zs,i,n,t,round(median,3),round(left,3),round(right,3))
        plt.annotate(st, xy=(0.60, 0.70), xycoords='axes fraction',fontsize=15)
        if delog:
            plt.xlabel('O/H',fontsize=18)
        else:
            plt.xlabel('12+log(O/H)',fontsize=18)
        plt.ylabel('counts',fontsize=18)
        plt.axvline(x=median,linewidth=2,color='black',ls='--')
        plt.savefig(outfile,clobber=False,format='pdf')
        
        ###print out the confidence interval###
#        print name, ':\t%f +- %f'%((left+right)/2.,(right-left)/2.)

        print '{0:40} {1:>13.3f} - {2:>7.3f} + {3:>7.3f}'.format(name, round(median,3), round(median-left,3), round(right-median,3))

        return "%f, %f, %f"%(round(median,3), round(median-left,3), round(right-median,3))

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
def run((filename, flux, err, path), nsample,delog=False):
    ###flux and err must be same dimensions
    newnsample=int(nsample+0.1*nsample)
    p=os.path.join(path,'..')
    if flux.shape != err.shape:
        print "flux and err must be of same dimensions"
        return
    nm = flux.shape[1]

    ###retrieve the metallicity keys
    Zs= metallicity.get_keys()

    ###make necessary paths
    if not os.path.exists(os.path.join(p,'bins','%s'%filename)):
        os.makedirs(os.path.join(p,'bins','%s'%filename))
    if not os.path.exists(os.path.join(p,'bins','%s'%filename,'hist')):
        os.makedirs(os.path.join(p,'bins','%s'%filename,'hist'))
    binp=os.path.join(p,'bins','%s'%filename)

    if VERBOSE: print "output files will be stored in ",binp

    ###Sample 'nsample' points from a gaussian###
    ##a gaussian centered on 0 with std 1
    mu=0
    sigma=1
    sample=np.random.normal(mu,sigma,newnsample)
    
    ###save this sampled gaussians into a png file###
    count, bins, ignored = plt.hist(sample, 40,normed=1)
    plt.plot(bins,gaussian(bins,mu,sigma))
    plt.title("Sampled")
    st="n=%d"%newnsample
    plt.annotate(st, xy=(0.70, 0.85), xycoords='axes fraction')

    plt.savefig(os.path.join(binp,'%s_n%d_sample.png'%(filename,newnsample)))
    plt.clf()

    
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
        t=metallicity.calculation(temp,nm,disp=VERBOSE)
        for key in Zs:
            res[key].append(t[key])
            
    #make the result a np.array
    for key in Zs:
        res[key]=np.array(res[key])
        
    if VERBOSE: print "Iteration Complete"
    
    ###Bin the results and save the uncertainty###
    print '{0:40} {1:>13} - {2:>7} + {3:>7}'.format("diagnostic", "metallicity","left", "right")
    for i in range(nm):
        fi=open(os.path.join(binp,'%s_n%d_i%d.csv'%(filename,nsample,i)),'w')
        fi.write("%s, Median (Z), Left, Right\n"%filename)

        print "measurement %d-------------------------------------------------------------"%i

        for key in Zs:
            s=key+", "+savehist(res[key][:,i],filename,key,nsample,i,binp,delog=delog)+'\n'
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
    parser.add_argument('filename', metavar='<filename>', type=str, help="the common filename")
    parser.add_argument('nsample', metavar='N', type=int, help="number of iterations")
    parser.add_argument('--clobber',default=False, action='store_true', help="replace existing output")
    parser.add_argument('--delog',default=False, action='store_true', help="result in natural, not log space. default is log space")
    parser.add_argument('--binmode', default='t', type=str, choices=['d','s','t'], help='method to determine bin size {d: Duanes formula, s: n^1/2, t: 2*n**1/3(default)}')
    parser.add_argument('--path',   default=None, type=str, help="input/output path (must contain the input _max.txt and _min.txt files in a subdirectory sn_data)")
    parser.add_argument('--verbose',default=False, action='store_true', help="verbose mode")

    args=parser.parse_args()

    global CLOBBER
    global VERBOSE
    global BINMODE
    CLOBBER=args.clobber
    VERBOSE=args.verbose
    BINMODE=args.binmode
    if args.path:
        path=args.path
    else:
        assert (os.getenv("MCMetdata"))," pass a path or set up the environmental variable MCMetdata pointing to the path where the _min _max _med files live"
        path=os.getenv("MCMetdata")
    assert(os.path.isdir(path)),"pass a path or set up the environmental variable MCMetdata pointing to the path where the _min _max _med files live"



    if args.nsample>0:
        fi=input_format(args.filename, path=path)
        if fi!=-1:
            run(fi,args.nsample,delog=args.delog)

    else:
        print "nsample must be positive number"
    
if __name__ == "__main__":
    main()
    #files=['sn2006ss','ptf10eqi-z']
    #filename=files[1]
    #nsample=10000
    

