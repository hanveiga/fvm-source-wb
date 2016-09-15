import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import math
import sys
import os

shape = [int(sys.argv[2]),int(sys.argv[2])]
fig=plt.figure()
ax = fig.add_subplot(121, projection='3d')
ax2 = fig.add_subplot(122)
folder_path = sys.argv[1]

def animate(filename):
        #m=np.loadtxt(sys.argv[1]+str(n).zfill(5)+".dat")
        m=np.loadtxt(folder_path+'/'+filename)
        X = np.array(m[:,0]).reshape(shape,order='C')
        Y = np.array(m[:,1]).reshape(shape,order='C')
        F = np.array(m[:,2]).reshape(shape,order='C')
        ax.cla()
	ax.plot_surface(X,Y,F,alpha=0.7,color="cyan",rstride=1, cstride=1,linewidth=0)
	ax2.cla()
        mid_pt = int(int(sys.argv[2])/2.)
        rho = np.array(m[:,2]).reshape(shape,order='C')
        velx = np.array(m[:,3]).reshape(shape,order='C')
        vely = np.array(m[:,4]).reshape(shape,order='C')
        E =  np.array(m[:,5]).reshape(shape,order='C') 
        ax2.plot(X[:,mid_pt],rho[mid_pt,:],"b")
        ax2.plot(X[:,mid_pt],velx[mid_pt,:],"r")
	ax2.plot(X[:,mid_pt],vely[mid_pt,:],"g")
	ax2.plot(X[:,mid_pt],E[mid_pt,:],"y")
	ax.set_zlim3d(-1,4)
	ax2.set_ylim(-2,2)
	#plt.title("n = "+str(n))
	#fig.canvas.draw_idle()  

list_of_names = [f for f in os.listdir(folder_path) if f != '.DS_Store']
list_of_names.sort()
list_of_names.insert(0, list_of_names.pop(-1))
#print list_of_names
anim = animation.FuncAnimation(fig, animate, list_of_names, blit=True)
anim.save('sw.mp4', writer = 'ffmpeg',fps=5)
