import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

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
        if err_count<i/2:    
            b=np.append(b,strs,axis=0)
            j+=1
    b.resize(j,i)
    b=np.transpose(b)
    a.close()

    return b,j
def f(x,a,b):
    return a*x+b

RNAME="ptf10eqi"
ZNAME="ptf10eqi-z"

###open everything
p=os.path.abspath('..')

#open radius file
r,n=readfile(p+"\\r_data\\"+RNAME+"_rr.txt")

#open Z % Z_err file
z=ZNAME
Zs=[]
for i in range(n):
    temp=np.genfromtxt(p+"\\bins\\%s\\%s_n10000_i%d.csv"%(z,z,i),delimiter=',',dtype=(str, float),autostrip=True)
    Zs.append(temp)
Zs=np.array(Zs)

###make things nice to plot
#Z: Zs[i][j][1], Z_err: Zs[i][j][2] pos: r[1,i] pos_err: r[2:i]
#j in range(1,Zs.shape[1])

x=r[1,:]
x_e=r[2,:]
for j in range(1,Zs.shape[1]):
    y=Zs[:,j,1].astype(float)
    y_e=Zs[:,j,2].astype(float)

    if np.sum(y)<=0:
        #if the metalicity was not calculated, the default value was set to -1
        print Zs[0,j,0]+"was not calcutated"
    
    else:
        #format the plt
        plt.clf()
        plt.title(Zs[0,j,0])
        plt.xlabel("Position")
        plt.ylabel("Metalicity")

        #plot the data points
        plt.plot(x,y,marker='o',linestyle='none')
        plt.errorbar(x,y,xerr=x_e, yerr=y_e, fmt=None)

        #fit from y err
        si=y_e
        py, pcov= curve_fit(f,x,y,sigma=si)
        resy = y - f(x,py[0],py[1])
        fresy=sum(resy**2)
        plt.plot(x,f(x,py[0],py[1]),color='red',label='y_err')

        #fit from x err
        ##using y_eff = x_error * dy/dx
        ##where the slope previously calculated with y_err was used for dy/dx
        si=x_e * py[0]
        px, pcov= curve_fit(f,x,y,sigma=si)
        resx = y - f(x,px[0],px[1])
        fresx=sum(resx**2)
        plt.plot(x,f(x,px[0],px[1]),color='blue', label='x_err')

        plt.legend()
        plt.savefig("%s_plot.png"%Zs[0,j,0])
        print Zs[0,j,0]+"\ty_fres:%f\tx_fres:%f"%(fresy,fresx)

        
