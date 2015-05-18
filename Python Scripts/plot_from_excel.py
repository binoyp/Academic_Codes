# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 21:19:21 2015

@author: Binoy Pilakkat
"""
from xlwings import Range,Workbook,Sheet
import matplotlib.pyplot as plt
import numpy as np

def multiplot(col,data):
    fig = plt.figure()
#    fig.suptitle('Grid Convergence')
    nc = len(col)
    for i in range(1,nc):
        ax = fig.add_subplot(nc-1,1,i)
#        ax = fig.add_subplot(3,2,i)
        ax.plot(data[:,0],data[:,i],'^:r')
#        ax.set_title(col[0],fontsize=10)
        ax.set_xlabel(col[0],fontsize=8)
        ax.set_ylabel(col[i])
        ax.grid(1)
    plt.tight_layout()
    
def multiplot_custom(col,data):
    fig = plt.figure()
#    fig.suptitle('Grid Convergence')
    nc = len(col)
    print nc
    count =1
    for i in range(1,nc,2):
        ax = fig.add_subplot((nc-1)/2,1,count)
        ax.plot(data[:,0],data[:,i],'^:b',label=col[i])
#        ax.set_title(col[0],fontsize=10)
        ax.set_xlabel(col[0],fontsize=8)
        ax.set_ylabel(col[i])
#        ax.set_ylabel(r'$C_W$')
        ax.legend(loc=5,prop={'size':9})
        ax.grid(1)
        
        ax2 = ax.twinx()
        ax2.plot(data[:,0],data[:,i+1],'D:r',label="Relative Variation")
        ax2.set_ylabel('%')
        count +=1
        for tl in ax2.get_yticklabels():
            tl.set_color('r')
        ax2.legend(prop={'size':9})
    plt.tight_layout()    
    
def singleplot(col,data):
    styles = ['-', '--', '-.', ':']
    markers = list('+ox^psDv')
#    fig = plt.figure()
#    fig.suptitle('Grid Convergence')
    nc = len(col)
    fig,ax = plt.subplots()
    for i in range(1,nc):
        s = styles[i % len(styles)]
        m = markers[i % len(markers)]
#        ax = fig.add_subplot(nc-1,1,i)
        ax.plot(data[:,0],data[:,i],marker=m, linewidth=2, linestyle=s,label= col[i])
#        ax.set_title(col[0],fontsize=10)
        ax.set_xlabel(col[0],fontsize=8)
#        ax.set_ylabel(col[i])
        ax.grid(1)
    ax.legend(loc=1)
    plt.tight_layout()   
    
    
def singleplot_custom(col,data,xc,yc):
    styles = ['-', '--', '-.', ':']
    markers = list('+ox^psDv')

    fig,ax = plt.subplots()

    s = styles[0]
    m = 's'
    mcount = 0
    print col

#        ax = fig.add_subplot(nc-1,1,i)
    l1 =ax.plot(data[:,xc],data[:,yc],marker='o', linewidth=2,label= col[yc])
#        ax.set_title(col[0],fontsize=10)
    ax2 = ax.twinx()
    l2= ax2.plot(data[:,xc],data[:,yc+1],'D:r',label=col[yc+1])
    ax2.set_ylabel('time(s)')
    ax.set_xlabel('Grid Size',fontsize=12)
    mcount += 1
#        ax.set_ylabel(col[i])
    ax.grid(1)
    ax.legend(loc=2)
    ax2.legend(loc=1)
    plt.tight_layout() 
     
#wb = Workbook(r'E:\Academic\ECN\CFD\Resistance Report\q2.xlsx')
wb = Workbook(r'E:\Academic\ECN\Python Scripts\book.xlsx')
cols = Range(1,"A1:D1").value

data = Range(1,'A2',asarray=1).table.value

data =data.astype('float')
#multiplot(cols,data)
singleplot_custom(cols,data,0,2)