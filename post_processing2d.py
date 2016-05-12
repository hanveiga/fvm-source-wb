import numpy as np
from matplotlib import pyplot as plt
import sys
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as manimation

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
    b = plt.contourf(x_2,y_2,rho_2,30, label='curve 2')
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

def plot_image(curve1,curve2,size_x,picture_name):
    size_y=size_x
    shape = [int(size_x),int(size_y)]
    data_1 = np.loadtxt(curve1)
    x_1 = data_1[:,0].reshape(shape,order='C')
    y_1 = data_1[:,1].reshape(shape,order='C')
    f_1 = data_1[:,2].reshape(shape,order='C')
    # plt.plot(x_end,rho_end, label='curve 1')

    plt.imshow(f_1,interpolation='none')

    plt.axes().set_aspect('equal')
    plt.savefig(str(picture_name)+'.png') 

def get_image(snapfile,size_x):
    size_y=size_x
    shape = [int(size_x),int(size_y)]
    data_1 = np.loadtxt(snapfile)
    x_1 = data_1[:,0].reshape(shape,order='C')
    y_1 = data_1[:,1].reshape(shape,order='C')
    f_1 = data_1[:,2].reshape(shape,order='C')

    #a = plt.imshow(f_1,interpolation='sinc')

    return f_1
    
def make_movie(folder_path, size_x):
    # list files in folder alphabetically
    list_of_snapshots = [f for f in os.listdir(folder_path) if isfolder(f)]
    # pass file to plot_image
    movie = []
    for snap in list_of_snapshots:
        movie.append(get_image(snap,size_x))
    # output film

    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='Movie Test', artist='Matplotlib',
            comment='Movie support!')
    writer = FFMpegWriter(fps=15, metadata=metadata)

    with writer.saving(fig, name+".mp4", 100):
        for step in movie:
            plt.imshow(step)
            writer.grab_frame()
            plt.clf()


if __name__=='__main__':
    #plot_curves(sys.argv[1],sys.argv[2],sys.argv[3], sys.argv[4])
    plot_image(sys.argv[1],sys.argv[2],sys.argv[3], sys.argv[5])
