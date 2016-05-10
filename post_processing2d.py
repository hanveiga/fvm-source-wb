import numpy as np
from matplotlib import pyplot as plt
import sys
from mpl_toolkits.mplot3d import Axes3D

def plot_curves(curve1, curve2, size_x, picture_name):
    """ function eats .dat files """
    size_y=size_x
    shape = [int(size_x),int(size_y)]
    data_1 = np.loadtxt(curve1)
    x_1 = data_1[:,0].reshape(shape,order='C')
    y_1 = data_1[:,1].reshape(shape,order='C')
    f_1 = data_1[:,2].reshape(shape,order='C')
    # plt.plot(x_end,rho_end, label='curve 1')
    a = plt.contour(x_1, y_1, f_1,30, label='curve 1')
    print np.max(f_1)
    data_2 = np.loadtxt(curve2)
    x_2 = data_2[:,0].reshape(shape,order='C')
    y_2 = data_2[:,1].reshape(shape,order='C')
    rho_2 = data_2[:,2].reshape(shape,order='C')
    b = plt.contour(x_2,y_2,rho_2,30, label='curve 2')
    # b = plt.contour(xcoords, ycoords, step)
    print np.max(rho_2)
    plt.clabel(a, inline=1, fontsize=10)
    plt.clabel(b, inline=1, fontsize=10)
    plt.legend()
    
    cbar = plt.colorbar(b)
    # ymax = max(np.max(rho_end1),np.max(rho_end))
    # ymin = min(np.min(rho_end),np.min(rho_end1))

    # xmax = max(np.max(x_end1),np.max(x_end))
    # xmin = min(np.min(x_end),np.min(x_end1))
    
    #print xmax, xmin, ymax, ymin
    #plt.axis([0,1/6.,0,1])

    plt.axes().set_aspect('equal')
    plt.savefig(str(picture_name)+'.png', dpi=300) 
  
   

if __name__=='__main__':
    plot_curves(sys.argv[1],sys.argv[2],sys.argv[3], sys.argv[4])
