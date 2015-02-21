import os,sys
import numpy as np
import pylab as pl
import fedMCZ_err
from scipy import stats, interpolate

#pickle may not be installed
NOPICKLE=False
try:
    import pprint, pickle
except:
    NOPICKLE=True

def fitdistrib(picklefile):
    if NOPICKLE:
        print "you must install pickle to read in the distributions and fit them"
        return -1

    assert( os.path.isfile(picklefile)), "missing pickled file %s"%picklefile
    pklfile = open(picklefile, 'rb')
    res=pickle.load(pklfile)

    testingdiags=['E(B-V)','D02','Z94','KD02comb_updated']
    ebvs=res['E(B-V)'].T
    nm=len(ebvs)
    Ndata= len(ebvs[0])
    assert( Ndata>0), "something is wrong with your distribution"


    invalids=[sum(np.isnan(res['D02'].T[i]))+sum(np.isnan(res['KD02comb_updated'].T[i]))+sum(np.isnan(res['Z94'].T[i])) for i in range(nm)]
    myi, = np.where(invalids==min(invalids))
    
    try: myi=myi[1] #in case there are more than one solutions
    except: 
        try:myi=myip[0]
        except:pass
    print "myi",myi
    assert( Ndata-invalids[myi]>0),  "something is wrong with your distribution"

    xold=np.zeros(10)
    fig = pl.figure(figsize=(15, 15))
    for i,d in enumerate(testingdiags):
        ax = fig.add_subplot(2,2,i+1)
        print d

        for f in [0.1,0.25,0.5,0.75,1]:
            n0=int(f*Ndata)
            x= np.random.choice(res[d].T[myi],n0,replace=False)
            print res[d].T[myi]
            x=x[x>0]
            
            x.sort()
            print x
            n=len(x)
            Px_cuml = np.linspace(0, 1, n)
            
            # set up an interpolation of the inverse cumulative distribution
            tck = interpolate.splrep(Px_cuml, x)
            
            # sample evenly along the cumulative distribution, and interpolate
            Px_cuml_sample = np.linspace(0, 1, 10 * n)
            x_sample = interpolate.splev(Px_cuml_sample, tck)
            indices = np.linspace(0, n - 1, 20).astype(int)
            
            ax.plot(x[indices], Px_cuml[indices], 'o', lw=0, label="%d"%n0)
            ax.plot(x, Px_cuml, '-k')
            ax.set_title('Cumulative Distribution for %s'%d)
            ax.set_xlabel('$x$')
            ax.set_ylabel('$p(<x)$')
            D, p = stats.ks_2samp(x, xold)
            print "KS test for %s n = %d D = %.2g; p = %.2g" % (d,n0,D, p)
            xold=x.copy()
        ax.legend(scatterpoints=1,loc=4)
    pl.savefig(picklefile.replace('.pkl','_testcomplete.png'))
    pl.show()
    
