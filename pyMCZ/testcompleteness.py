import os
import sys
import numpy as np
import pylab as pl

from scipy import stats  # , interpolate
import pylabsetup

#pickle may not be installed
NOPICKLE = False
try:
    import pickle
except ImportError:
    NOPICKLE = True


def fitdistrib(picklefile, scales=None):
    from matplotlib.font_manager import findfont, FontProperties
    from matplotlib.ticker import FormatStrFormatter
    majorFormatter = FormatStrFormatter('%.2f')

    fs = 20  # FontProperties.get_size(FontProperties())
    if 'Time' not in  findfont(FontProperties()):
        fs = 20
    print "FONT: %s, %d" % (findfont(FontProperties()), fs)
    fs = fs - 2
    params = {'axes.labelsize': fs,
              'xtick.labelsize': fs - 3,
              'ytick.labelsize': fs - 3,
              'legend.fontsize': fs - 2,
              'legend.handletextpad': 0.2
          }
    pl.rcParams.update(params)


    if NOPICKLE:
        print "you must install pickle to read in the distributions and fit them"
        return -1
    assert(os.path.isfile(picklefile)), "missing pickled file %s" % picklefile
    pklfile = open(picklefile, 'rb')
    res = pickle.load(pklfile)
    print scales
    if scales:
        testingdiags = scales.split(',')
    else:
        testingdiags = ['E(B-V)', 'D02', 'M13_N2', 'KD02comb']

    print "\n\n###testing with scales: ", testingdiags,
    print "###\n\n"
    ebvs = res['E(B-V)'].T
    nm = len(ebvs)
    Ndata = len(ebvs[0])
    assert(Ndata > 0), "something is wrong with your distribution"
    invalids = [sum(np.isnan(res[testingdiags[1]].T[i])) + sum(np.isnan(res[testingdiags[2]].T[i])) + sum(np.isnan(res[testingdiags[3]].T[i])) for i in range(nm)]
    myi, = np.where(invalids == min(invalids))
    try:
        myi = myi[1]  # in case there are more than one solutions
    except IndexError:
        try:
            myi = myi[0]
        except IndexError:
            pass
    xold = np.zeros(10)
    assert(Ndata - invalids[myi] > 0), "something is wrong with your distribution. Perhaps one of the default scale has all invalid values? try select different scales"
    # -- plotting utils
    fwid = 10.  # fig width
    rat = 1.   # fig aspect ratio
    offx = 0.08  # x offset in figure coords
    offy = 0.05  # y offset in figure coords
    psz = 0.4  # size of the plot in figure width coords
    nr = 2    # number of row plots
    nc = 2    # number of col plots

    fig, axx = pl.subplots(nr, nc, figsize=(fwid, fwid / rat))
    print fig

    for i, d in enumerate(testingdiags):
        ii, jj = int(i / 2), int((i + 1) / 2) - int(i / 2)
        ax = axx[ii][jj]
        ax.set_position([jj / float(nc) + offx, 1.0 - (ii + 1) / float(nr) + offy, psz, psz * rat])
        for f in [0.1, 0.25, 0.5, 0.75, 1]:
            n0 = int(f * Ndata)
            x = np.random.choice(res[d].T[myi], n0, replace=False)
            #print res[d].T[myi]
            x = x[x > 0]

            x.sort()
            #print x
            n = len(x)

            Px_cuml = np.linspace(0, 1, n)
            # set up an interpolation of the inverse cumulative distribution
            #tck = interpolate.splrep(Px_cuml, x)

            # sample evenly along the cumulative distribution, and interpolate
            #Px_cuml_sample = np.linspace(0, 1, 10 * n)

            #x_sample = interpolate.splev(Px_cuml_sample, tck)

            indices = np.linspace(0, n - 1, 20).astype(int)
            ax.set_ylim(0, 1.05)
            ax.plot(x[indices], Px_cuml[indices], 'o', lw=0, label="%d" % n0)
            ax.plot(x, Px_cuml, '-k')
            maxleft, maxright = min(x), max(x)
            D, p = stats.ks_2samp(x, xold)
            print "KS test for %s n = %d D = %.2g; p = %.2g" % (d, n0, D, p)
            xold = x.copy()
            lims = ax.set_xlim((maxleft, maxright))
        axratio = (lims[1] - lims[0]) / 1.05
        ax.set_aspect(aspect=axratio)

        ax.set_title('%s Cumulative Distribution' % d.replace('_', ' '), fontsize=fs)
        ax.set_xlabel('$x$')
        ax.set_ylabel('$p(<x)$')
        ax.xaxis.set_major_formatter(majorFormatter)

        ax.legend(scatterpoints=1, loc=4)
        xticks = ax.get_xticks()
        dx = xticks[-1] - xticks[-2]
        xticks = xticks[(xticks < maxright) * (xticks > maxleft)]
        while (maxright - xticks[-1]) < 0.25 * dx:
            xticks = xticks[:-1]
    pl.xticks(xticks, ['%.2f' % s for s in xticks])
    pl.savefig(picklefile.replace('.pkl', '_testcomplete.pdf'))
    pl.show()

if __name__ == "__main__":
    if len(sys.argv) > 2:
        print '''only argument allowed: name of the pickle file containing the distribution
        default is '../output/exampledata/exampledata_n2000.pkl'
       '''
        sys.exit()
    infile = '../output/exampledata/exampledata_n2000.pkl'
    if len(sys.argv) > 1:
        infile = sys.argv[1]
    fitdistrib(infile)
