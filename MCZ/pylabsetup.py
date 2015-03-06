import matplotlib as mpl
import pylab as plt
from pylab import rc


rc('axes', linewidth=1.2)
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['font.size'] = 18.
mpl.rcParams['font.family'] = 'serif'
#mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = [  'Times New Roman', 'Times','Palatino', 'Charter', 'serif']
mpl.rcParams['font.sans-serif'] = ['Helvetica']
mpl.rcParams['axes.labelsize'] = 18
mpl.rcParams['xtick.labelsize'] = 18.
mpl.rcParams['ytick.labelsize'] = 18.
mpl.rcParams['xtick.major.size']= 10.
mpl.rcParams['xtick.minor.size']= 5.
mpl.rcParams['ytick.major.size']= 10.
mpl.rcParams['ytick.minor.size']= 5.

fontsize=20
#mpl.rc('axes',  titlesize=fontsize)
#mpl.rc('axes',  labelsize=fontsize)
#mpl.rc('xtick', labelsize=fontsize)
#mpl.rc('ytick', labelsize=fontsize)
#mpl.rc('font', size=fontsize, family='serif', serif='Utopia',
#              style='normal', variant='normal',
#              stretch='normal', weight='normal')
mpl.rc('font',**{'family':'sans-serif', 'serif':['Bitstream Vera Serif'], 
                 'sans-serif':['Helvetica'], 'size':20, 
                 'weight':'normal'})
mpl.rc('axes',**{'labelweight':'bold', 'linewidth':1})
mpl.rc('ytick',**{'major.pad':8, 'color':'k'})
mpl.rc('xtick',**{'major.pad':8, 'color':'k'})
mpl.rc('mathtext',**{'default':'regular','fontset':'cm', 
                     'bf':'monospace:bold'})
params = {'legend.fontsize': 20,
          'legend.linewidth': 1,
          'legend.numpoints':1,
          'legend.handletextpad':1
      }

plt.rcParams.update(params)   
plt.minorticks_on()



