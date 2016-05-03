import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
class Field():

    def __init__(self,t,nsteps=8,dt=np.pi*1e-5,nx=384,ny=128+6):

        self.ymin = np.loadtxt('domain_y.dat')
        self.xmin = np.loadtxt('domain_x.dat')
        self.xmed = .5*(self.xmin[1:]+self.xmin[:-1])
        self.ymed = .5*(self.ymin[1:]+self.ymin[:-1])
        self.nx = nx
        self.ny = ny
        self.dt = dt
        self.t = t
        self.dens = [None]*nsteps
        self.vx = [None]*nsteps
        self.vy = [None]*nsteps

        s = nx*ny
        for i in range(nsteps):
            dat = np.fromfile('substep_%d_%d.dat'%(i,t))
            self.dens[i] = dat[:s].reshape(ny,nx)
            self.vy[i] = dat[s:2*s].reshape(ny,nx)
            self.vx[i] = dat[2*s:3*s].reshape(ny,nx)
        return

class Snap():

    def __init__(self,directory='temp_files/',nx=384,ny=128+6):

        dat = np.loadtxt('param_file.txt')
        self.alpha = dat[3]
        self.mp = dat[4]
        self.q = dat[5]
        self.omf = dat[6]
        self.flaringindex = dat[7]
        self.soft = dat[9]
        self.nx = nx
        self.ny = ny
        s = nx*ny
        dat = np.fromfile(directory+'output.dat')
        self.ymin = dat[:(ny+1)]
        dat=dat[(ny+1):]
        self.xmin = dat[:(nx+1)]
        dat=dat[(nx+1):]
        self.dens = dat[:s].reshape(ny,nx)
        self.vy = dat[s:2*s].reshape(ny,nx)
        self.vx= dat[2*s:3*s].reshape(ny,nx)
        self.pot= dat[3*s:4*s].reshape(ny,nx)

        self.ymed = .5*(self.ymin[1:]+self.ymin[:-1])
        self.xmed = .5*(self.xmin[1:] + self.xmin[:-1])

        dat = np.fromfile(directory+'output_init.dat')
        self.dens0 = dat[:s].reshape(ny,nx)
        self.vy0 = dat[s:2*s].reshape(ny,nx)
        self.vx0= dat[2*s:3*s].reshape(ny,nx)
        return
class Fargo():

    def __init__(self,i,directory='',nx=384,ny=128):

        self.nx = nx
        self.ny = ny
        self.dens = np.fromfile('gasdens%d.dat'%i).reshape(ny,nx)
        self.vy = np.fromfile('gasvy%d.dat'%i).reshape(ny,nx)
        self.vx = np.fromfile('gasvx%d.dat'%i).reshape(ny,nx)

        self.ymin = np.loadtxt('domain_y.dat')[3:-3]
        self.xmin = np.loadtxt('domain_x.dat')
        self.xmed = .5*(self.xmin[1:]+self.xmin[:-1])
        self.ymed = .5*(self.ymin[1:]+self.ymin[:-1])

        return
class Torque():

    def __init__(self,directory='temp_files/',nx=384,ny=128):

        dat = np.loadtxt('param_file.txt')
        self.alpha = dat[3]
        self.mp = dat[4]
        self.q = dat[5]
        self.omf = dat[6]
        self.flaringindex = dat[7]
        self.soft = dat[9]
        self.nx = nx
        self.ny = ny
        dat = np.fromfile(directory+'torque.dat')
        self.y = dat[:ny]
        dat=dat[ny:]
        self.Lt = dat[:ny]
        self.Ld = dat[ny:ny*2]
        self.Lw = dat[ny*2:ny*3]
        self.drFt = dat[ny*3:ny*4]
        self.drFd = dat[ny*4:ny*5]
        self.drFw = dat[ny*5:ny*6]
        self.Lamex = dat[ny*6:ny*7]
        self.Lamdep = dat[ny*7:ny*8]
        self.dtLt = dat[ny*8:ny*9]
        self.dtLd = dat[ny*9:ny*10]
        self.dtLw = dat[ny*10:ny*11]

        return

def compare(dat,fld,i=-1,relerror=False):
    fig,axes=plt.subplots(3,1,sharex=True)
    if type(fld.dens) != list:
        dens = fld.dens
        vx = fld.vx
        vy = fld.vy
    else:
        dens = fld.dens[i]
        vx = fld.vx[i]
        vy = fld.vy[i]

    try:
        errd = (dat.dens[3:-3,:] - dens)
        errvy = (dat.vy[3:-3,:] - vy)
        errvx = (dat.vx[3:-3,:] - vx)
    except ValueError:
        dens = dens[3:-3,:]
        vx = vx[3:-3,:]
        vy = vy[3:-3,:]
        errd = (dat.dens[3:-3,:] - dens)
        errvy = (dat.vy[3:-3,:] - vy)
        errvx = (dat.vx[3:-3,:] - vx)

    if relerror:
        errd /= dens
        errvy /= vy
        errvx /= vx


    ims = [None]*len(axes)
    ims[0]=axes[0].imshow(np.log10(abs(errd)+1e-16),origin='lower')
    ims[1]=axes[1].imshow(np.log10(abs(errvy)+1e-16),origin='lower')
    ims[2]=axes[2].imshow(np.log10(abs(errvx)+1e-16),origin='lower')
    axes[0].set_ylabel('$\\Delta\\Sigma$',fontsize=15,rotation=0)
    axes[1].set_ylabel('$\\Delta v_y$',fontsize=15,rotation=0)
    axes[2].set_ylabel('$\Delta v_x$',fontsize=15,rotation=0)
    dividers = [make_axes_locatable(ax) for ax in axes]
    caxes = [d.append_axes("top", size="20%", pad=0.05) for d in dividers]
    cbars = [plt.colorbar(im, cax=cax,orientation='horizontal') for im,cax in zip(ims,caxes)]
    for cb,ax in zip(cbars,axes):
        ax.minorticks_on()
        cb.ax.xaxis.set_ticks_position('top')
        ax.yaxis.set_label_position('right')
        ax.yaxis.labelpad = 20


