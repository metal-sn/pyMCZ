import numpy as np
import matplotlib.pyplot as plt
import scipy.stats.mstats as ssm
import metallicity
import os

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
##mode 'd' calculates this based on Doane's formula (default)
##mode 's' calculates this based on sqrt of number of data
##mode 't' calculates this based on 2*n**1/3
##############################################################################
def getbinsize(n,data,mode='t'):
    if mode=='d':
        g1=ssm.moment(data,moment=3)
        s1=np.sqrt(6.*(n-2.)/((n+1.)*(n+3.)))
        k=1+np.log2(n)+np.log2(1+np.abs(g1/s1))
    elif mode=='s':
        k=np.sqrt(n)
    elif mode=='t':
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
def savehist(data,filename,Zs,nsample,i,delog=False,mode='t'):
    p=os.path.abspath('..')
    name='%s_n%d_%s_i%d'%((filename,nsample,Zs,i))
    outfile=p+'\\bins\\%s\\hist\\%s.png'%(filename,name)
    plt.clf()

    ###de-log###
    if delog==True:
        data=np.power(10,np.absolute(data-12))
        
    data=np.sort(data)
    data=data[np.where(np.isfinite(data)==True)]
    
    ####kill outliers###
    mean=np.mean(data)
    std=np.std(data)
    data=data[np.where(data>mean-5*std)]
    data=data[np.where(data<mean+5*std)]
    n=data.shape[0]
    if data.shape[0]<=0:
        print name,'is blank' ##if no data
        return "-1, -1"
    
    try:
        ######find fit and save hist######
        ###find appropriate bin size###
        numbin=getbinsize(data.shape[0],data,mode)
        
        ###make hist###
        count, bins, ignored = plt.hist(data, numbin, normed=1.0)

        ###find error###
        l,r,t,fl,fr=err_est(count)
        y=np.zeros(len(bins))
        left=bins[l]+(bins[l+1]-bins[l])*fl
        right=bins[r]-(bins[r]-bins[r-1])*fr
        
        ###plot hist###
        plt.plot(bins,y)
        plt.axvspan(left,right,color='red',alpha=0.4)
        st='n=%d\nconfidence: %f\nleft: %f\nright: %f'%(n,t,left,right)
        plt.annotate(st, xy=(0.70, 0.80), xycoords='axes fraction')
        plt.title(name)
        if delog==True:
            plt.xlabel('O/H')
        else:
            plt.xlabel('12+log(O/H)')
        plt.ylabel('counts')
        plt.savefig(outfile)
        
        ###print out the confidence interval###
        print name, ':\t%f +- %f'%((left+right)/2.,(right-left)/2.)
        return "%f, %f"%((left+right)/2.,(right-left)/2.)

    except (OverflowError,AttributeError,ValueError):
        print data
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
def run((filename, flux, err), nsample,binmode='t'):
    p=os.path.abspath('..')
    ###flux and err must be same dimensions
    if flux.shape != err.shape:
        print "flux and err must be of same dimensions"
        return
    nm = flux.shape[1]

    ###retrieve the metallicity keys
    Zs= metallicity.get_keys()

    ###make necessary paths
    if not os.path.exists(p+'\\bins\\%s'%filename):
        os.makedirs(p+'\\bins\\%s'%filename)
        os.makedirs(p+'\\bins\\%s\\hist'%filename)
    binp=p+'\\bins\\%s\\'%filename
    
    ###Sample 'nsample' points from a gaussian###
    ##a gaussian centered on 0 with std 1
    mu=0
    sigma=1
    sample=np.random.normal(mu,sigma,nsample)
    
    ###save this sampled gaussians into a png file###
    count, bins, ignored = plt.hist(sample, 40,normed=1)
    plt.plot(bins,gaussian(bins,mu,sigma))
    plt.title("Sampled")
    st="n=%d"%nsample
    plt.annotate(st, xy=(0.70, 0.85), xycoords='axes fraction')
    plt.savefig(binp+filename+'_n%d_sample.png'%nsample)
    plt.clf()

    
    ###Start calculation###
    ## the flux to be feed to the calculation will be
    ## flux + error*i
    ## where i is the sampled gaussian    
    print "Starting iteration"

    #initialize the dictionary
    res={}
    for key in Zs:
        res[key]=[]

    #do the iterations
    for i in range(len(sample)):
        temp=flux+err*sample[i]
        t=metallicity.calculation(temp,nm)
        for key in Zs:
            res[key].append(t[key])
            
    #make the result a np.array
    for key in Zs:
        res[key]=np.array(res[key])
        
    print "Iteration Complete"
    
    ###Bin the results and save the uncertainty###
    for i in range(nm):
        fi=open(binp+'%s_n%d_i%d.csv'%(filename,nsample,i),'w')
        fi.write("%s, Metalicity, Uncertainty\n"%filename)

        for key in Zs:
            s=key+", "+savehist(res[key][:,i],filename,key,nsample,i,binmode)+'\n'
            fi.write(s)
        print "---------------------------------------------"
        fi.close()
    
    print "uncertainty calculation complete"

##############################################################################
##The input format generator
##############################################################################
def input_format(filename):
    p=os.path.abspath('..')
    
    if os.path.isfile(p+'\\sn_data\\%s_max.txt'%filename):
        if os.path.isfile(p+'\\sn_data\\%s_min.txt'%filename):
            if os.path.isfile(p+'\\sn_data\\%s_med.txt'%filename):
                return in_mmm(filename)
            else:
                return in_mm(filename)
    print "Unable to find _min and _max files in directory sn_data"
    return -1

def in_mmm(filename):
    p=os.path.abspath('..')
    p+='\\sn_data\\'
    ###Initialize###
    maxfile=p+filename+"_max.txt"
    medfile=p+filename+"_med.txt"
    minfile=p+filename+"_min.txt"
    
    ###read the max, med, min flux files###    
    maxf,nm=readfile(maxfile)
    medf,num=readfile(medfile)
    minf,nn=readfile(minfile)

    ###calculate the flux error as 1/2 [(max-med)+(med-min)]###
    err=0.5*((maxf-medf)+(medf-minf))
    return (filename, medf, err)

def in_mm(filename):
    p=os.path.abspath('..')
    p+='\\sn_data\\'
    ###Initialize###
    maxfile=p+filename+"_max.txt"
    minfile=p+filename+"_min.txt"
    
    ###read the max, med, min flux files###    
    maxf,nm=readfile(maxfile)
    minf,nn=readfile(minfile)
    ###calculate the flux error as 1/2 [(max-min)]###
    err=0.5*(maxf - minf)
    medf= minf+err
    return (filename, medf, err)
##############################################################################
##############################################################################    
#### How to Use this code:
##############################################################################    
##############################################################################    
####example) if you have:
####     sn2006ss_max.txt
####     sn2006ss_med.txt
####     sn2006ss_min.txt
####
####     and you want to sample 5000 points
####     do:
####
#### filename="sn2006ss"
#### nsample=5000
####
##############################################################################
####Just edit this part to use
##############################################################################

files=['sn2006ss','ptf10eqi-z']
filename=files[1]

nsample=10000
fi=input_format(filename)
if fi!=-1:
    run(fi,nsample)
