class Mesh():
    """
    Mesh class, for keeping all the mesh data.
    Input: directory [string] -> place where the domain files are.
    """
    def __init__(self, directory=""):
        if len(directory) > 1:
            if directory[-1] != '/':
                directory += '/'
        try:
            domain_x = np.loadtxt(directory+"domain_x.dat")
        except IOError:
            print "IOError with domain_x.dat"
        try:
            #We have to avoid ghost cells!
            domain_y = np.loadtxt(directory+"domain_y.dat")[3:-3]
        except IOError:
            print "IOError with domain_y.dat"
        self.xm = domain_x #X-Edge
        self.ym = domain_y #Y-Edge

        self.xmed = 0.5*(domain_x[:-1] + domain_x[1:]) #X-Center
        self.ymed = 0.5*(domain_y[:-1] + domain_y[1:]) #Y-Center

        #(Surfaces taken from the edges)
        #First we make 2D arrays for x & y, that are (theta,r)
        T,R = meshgrid(self.xm, self.ym)
        R2  = R*R
        self.surf = 0.5*(T[:-1,1:]-T[:-1,:-1])*(R2[1:,:-1]-R2[:-1,:-1])

class Parameters():
    """
    Class for reading the simulation parameters.
    input: string -> name of the parfile, normally variables.par
    """
    def __init__(self, directory=''):
        if len(directory) > 1:
            if directory[-1] != '/':
                directory += '/'
        try:
            params = open(directory+"variables.par",'r') #Opening the parfile
        except IOError:                  # Error checker.
            print  paramfile + " not found."
            return
        lines = params.readlines()     # Reading the parfile
        params.close()                 # Closing the parfile
        par = {}                       # Allocating a dictionary
        for line in lines:             #Iterating over the parfile
            name, value = line.split() #Spliting the name and the value (first blank)
            try:
                float(value)           # First trying with float
            except ValueError:         # If it is not float
                try:
                    int(value)         #                   we try with integer
                except ValueError:     # If it is not integer, we know it is string
                    value = '"' + value + '"'
            par[name] = value          # Filling the dictory
        self._params = par             # A control atribute, actually not used, good for debbuging
        for name in par:               # Iterating over the dictionary
            exec("self."+name.lower()+"="+par[name]) #Making the atributes at runtime


class Field(Mesh, Parameters):
    """
    Field class, it stores all the mesh, parameters and scalar data
    for a scalar field.
    Input: field [string] -> filename of the field
           staggered='c' [string] -> staggered direction of the field.
                                      Possible values: 'x', 'y', 'xy', 'yx'
           directory='' [string] -> where filename is
           dtype='float64' (numpy dtype) -> 'float64', 'float32',
                                             depends if FARGO_OPT+=-DFLOAT is activated
    """
    def __init__(self, field, staggered='c', directory='', dtype='float64'):
        if len(directory) > 1:
            if directory[-1] != '/':
                directory += '/'
        Mesh.__init__(self, directory) #All the Mesh attributes inside Field!
        Parameters.__init__(self, directory) #All the Parameters attributes inside Field!

        #Now, the staggering:
        if staggered.count('x')>0:
            self.x = self.xm[:-1] #Do not dump the last element
        else:
            self.x = self.xmed
        if staggered.count('y')>0:
            self.y = self.ym[:-1]
        else:
            self.y = self.ymed

        self.data = self.__open_field(directory+field,dtype) #The scalar data is here.

    def __open_field(self, f, dtype):
        """
        Reading the data
        """
        field = fromfile(f, dtype=dtype)
        return field.reshape(self.ny, self.nx)

    def plot(self, log=False, cartesian=False, cmap='Oranges_r', **karg):
        """
        A layer to plt.imshow or pcolormesh function.
        if cartesian = True, pcolormesh is launched.
        """
        ax = gca()
        if log:
            data = np.log(self.data)
        else:
            data = self.data
        if cartesian:
            T,R = meshgrid(self.x,self.y)
            X = R*cos(T)
            Y = R*sin(T)
            pcolormesh(X,Y,data,cmap=cmap,**karg)
        else:
            ax.imshow(data, cmap = cmap, origin='lower',aspect='auto',
                      extent=[self.x[0],self.x[-1],self.y[0],self.y[-1]],
                      **karg)

    def contour(self, log=False, cartesian=False, **karg):
        if log:
            data = np.log(self.data)
        else:
            data = self.data
        ax = gca()
        T,R = meshgrid(self.x,self.y)
        if cartesian:
            X = R*cos(T)
            Y = R*sin(T)
            contour(X,Y,data,**karg)
        else:
            contour(T,R,data,**karg)


    def shift_field(Field,direction):
        import copy
        """
        Half cell shifting along the direction provided by direction
        direction can be ('x','y', 'xy', 'yx').

        After a call of this function, Field.xm/xmed has not
        sense anymore (it is not hard to improve).
        """
        F = copy.deepcopy(Field)
        if direction.count('x')>0:
            F.data = 0.5*(Field.data[:,1:]+Field.data[:,:-1])
            F.x = 0.5*(Field.x[1:]+Field.x[:-1])
        if direction.count('y')>0:
            F.data = 0.5*(F.data[1:,:]+F.data[:-1,:])
            F.y = 0.5*(F.y[1:]+F.y[:-1])

        F.nx = len(F.x)
        F.ny = len(F.y)

        return F


    def cut_field(Field, direction, side):
        """
        Cutting a field:
        Input: field --> a Field class
               axis  --> 'x', 'y' or 'xy'
               side  --> 'p' (plus), 'm' (minnus), 'pm' (plus/minnus)
        """
        import copy

        cutted_field = copy.deepcopy(Field)
        ny,nx = Field.ny, Field.nx
        mx = my = px = py = 0

        if direction.count('x')>0:
            if side.count('m')>0:
                mx = 1
            if side.count('p')>0:
                px = 1
        if direction.count('y')>0:
            if side.count('m')>0:
                my = 1
            if side.count('p')>0:
                py = 1

        cutted_field.data = Field.data[my:ny-py,mx:nx-px]
        cutted_field.x = cutted_field.x[mx:nx-px]
        cutted_field.y = cutted_field.y[my:ny-py]

        return cutted_field

    def vector_field(vx,vy, **karg):
        nsx = nsy = 3
        T,R = meshgrid(vx.x[::nsx],vx.y[::nsy])
        X = R*cos(T)
        Y = R*sin(T)
        vx = vx.data[::nsy,::nsx]
        vy = vy.data[::nsy,::nsx]
        U = vy*cos(T) - vx*sin(T)
        V = vy*sin(T) + vx*cos(T)
        ax = gca()
        ax.quiver(X,Y,U,V,scale=5,pivot='midle', **karg)

    def euler(vx, vy, x, y, reverse):
        """
        Euler integrator for computing the streamlines.
        Parameters:
        ----------

        x,y: Float.
             starter coordinates.
        reverse: Boolean.
                 If reverse is true, the integratin step is negative.
        Reverse inverts the sign of the velocity

        Output
        ------

        (dx,dy): (float,float).
                 Are the azimutal and radial increment.
                 Only works for cylindrical coordinates.
        """
        sign = 1.0
        if reverse:
            sign = -1
        vphi = get_v(vx,x,y)
        vrad = get_v(vy,x,y)
        if vphi == None or vrad == None: #Avoiding problems...
            return None,None
        l = np.min((((vx.xmax-vx.xmin)/vx.nx),((vx.ymax-vx.ymin)/vx.ny)))
        h = 0.5*l/np.sqrt((vphi**2+vrad**2))

#        return sign*h*np.array([vphi/y,vrad])
        return sign*h*np.array([vphi,vrad])


    def bilinear(x,y,f,p):
        """
        Computing bilinear interpolation.
        Parameters
        ----------
        x = (x1,x2); y = (y1,y2)
        f = (f11,f12,f21,f22)
        p = (x,y)
        where x,y is the interpolated point and
        fij is the value of the function at the
        point (xi,yj).
        Output
        ------
        f(p): Float.
              The interpolated value of the function f(p) = f(x,y)
        """

        xp  = p[0]; yp   = p[1]; x1  = x[0]; x2  = x[1]
        y1  = y[0]; y2  = y[1];  f11 = f[0]; f12 = f[1]
        f21 = f[2]; f22 = f[3]
        t = (xp-x1)/(x2-x1);    u = (yp-y1)/(y2-y1)

        return (1.0-t)*(1.0-u)*f11 + t*(1.0-u)*f12 + t*u*f22 + u*(1-t)*f21

    def get_v(v, x, y):
        """
        For a real set of coordinates (x,y), returns the bilinear
        interpolated value of a Field class.
        """

        i = int((x-v.xmin)/(v.xmax-v.xmin)*v.nx)
        j = int((y-log(v.ymin))/(log(v.ymax/v.ymin))*v.ny)

        if i<0 or j<0 or i>v.data.shape[1]-2 or j>v.data.shape[0]-2:
            return None

        f11 = v.data[j,i]
        f12 = v.data[j,i+1]
        f21 = v.data[j+1,i]
        f22 = v.data[j+1,i+1]
        try:
            x1  = v.x[i]
            x2  = v.x[i+1]
            y1  = log(v.y[j])
            y2  = log(v.y[j+1])
            return bilinear((x1,x2),(y1,y2),(f11,f12,f21,f22),(x,y))
        except IndexError:
            return None



    def get_stream(vx, vy, x0, y0, nmax=1000000, maxlength=4*np.pi, bidirectional=True, reverse=False):
        """
        Function for computing a streamline.
        Parameters:
        -----------

        x0,y0: Floats.
              Initial position for the stream
        nmax: Integer.
              Maxium number of iterations for the stream.
        maxlength: Float
                   Maxium allowed length for a stream
        bidirectional=True
                      If it's True, the stream will be forward and backward computed.
        reverse=False
                The sign of the stream. You can change it mannualy for a single stream,
                but in practice, it's recommeneded to use this function without set reverse
                and setting bidirectional = True.

        Output:
        -------

        If bidirectional is False, the function returns a single array, containing the streamline:
        The format is:

                                          np.array([[x],[y]])

        If bidirectional is True, the function returns a tuple of two arrays, each one with the same
        format as bidirectional=False.
        The format in this case is:

                                (np.array([[x],[y]]),np.array([[x],[y]]))

        This format is a little bit more complicated, and the best way to manipulate it is with iterators.
        For example, if you want to plot the streams computed with bidirectional=True, you should write
        some similar to:

        stream = get_stream(x0,y0)
        ax.plot(stream[0][0],stream[0][1]) #Forward
        ax.plot(stream[1][0],stream[1][1]) #Backward

        """

        if bidirectional:
            s0 = get_stream(vx, vy, x0, y0, reverse=False, bidirectional=False, nmax=nmax,maxlength=maxlength)
            s1 = get_stream(vx, vy, x0, y0, reverse=True,  bidirectional=False, nmax=nmax,maxlength=maxlength)
            return (s0,s1)

        l = 0
        x = [x0]
        y = [y0]

        for i in xrange(nmax):
            ds = euler(vx, vy, x0, y0, reverse=reverse)
            if ds[0] == None:
                if(len(x)==1):
                    print "There was an error getting the stream, ds was NULL (see get_stream)."
                break
            l += np.sqrt(ds[0]**2+ds[1]**2)
            dx = ds[0]
            dy = ds[1]
            if(np.sqrt(dx**2+dy**2)<1e-13):
                print "Warning (get_stream): ds is very small, maybe you're in a stagnation point."
                print "Try selecting another initial point."
                break
            if l > maxlength:
                print "maxlength reached: ", l
                break
            x0 += dx
            y0 += dy
            x.append(x0)
            y.append(y0)
        return np.array([x,y])

    def get_random_streams(vx, vy, xmin=None, xmax=None, ymin=None, ymax=None, n=30, nmax=100000):
        if xmin == None:
            xmin = vx.xmin
        if ymin == None:
            ymin = vx.ymin
        if xmax == None:
            xmax = vx.xmax
        if ymax == None:
            ymax = vx.ymax
        X = xmin + np.random.rand(n)*(xmax-xmin)
        Y = ymin + np.random.rand(n)*(ymax-ymin)
        streams = []
        counter = 0
        for x,y in zip(X,Y):
            stream = get_stream(vx, vy, x, y, nmax=nmax, bidirectional=True)
            streams.append(stream)
            counter += 1
            print "stream ",counter, "done"
        return streams

    def plot_random_streams(streams, cartesian=False, **kargs):
        ax = plt.gca()
        print np.shape(streams)
        for stream in streams:
            for sub_stream in stream:
                if cartesian:
                    ax.plot(sub_stream[1]*cos(sub_stream[0]),sub_stream[1]*sin(sub_stream[0]),**kargs)
                else:
                    ax.plot(sub_stream[0],sub_stream[1],**kargs)

from scipy.interpolate import interp1d,interp2d
class Sim():
    def __init__(self,i,directory='',p=0):
    	if directory != '':
    		if directory[-1] != '/':
    			directory += '/'
        self.directory = directory
    	self.dens = Field('gasdens{0:d}.dat'.format(i),directory=directory)
    	self.vp = Field('gasvx{0:d}.dat'.format(i),directory=directory,staggered='x')
    	self.vr = Field('gasvy{0:d}.dat'.format(i),directory=directory,staggered='y')

        # self.vp.data = 0.5*(self.vp.data[:,1:]+self.vp.data[:,:-1])
        # self.vp.x = 0.5*(self.vp.x[1:]+self.vp.x[:-1])
        # self.vr.data = 0.5*(self.vr.data[1:,:]+self.vr.data[:-1,:])
        # self.vr.y = 0.5*(self.vr.y[1:]+self.vr.y[:-1])
        #
        # self.vr.nx = len(self.vr.x)
        # self.vr.ny = len(self.vr.y)
        # self.vp.nx = len(self.vp.x)
        # self.vp.ny = len(self.vp.y)

        self.vp = self.vp.shift_field('x')
        self.vr = self.vr.shift_field('y')
        self.vp = self.vp.cut_field(direction='y',side='p')
        self.vr = self.vr.cut_field(direction='x',side='p')
        self.dens = self.dens.cut_field(direction='xy', side='p')

    	self.r  = self.dens.y
    	self.phi = self.dens.x
        self.pp, self.rr = meshgrid(self.phi,self.r)
    	self.dbar = self.dens.data.mean(axis=1)
    	self.vrbar = self.vr.data.mean(axis=1)
    	self.vpbar = self.vp.data.mean(axis=1)
    	self.mdot = -2*pi*self.r*(self.vr.data*self.dens.data).mean(axis=1)
    	_,self.px,self.py,self.pz,self.pvx,self.pvy,self.pvz,self.mp,self.t,self.omf  = loadtxt(directory+'planet{0:d}.dat'.format(p))[i,:]
    	self.a = sqrt(self.px**2  + self.py**2)
    	self.nu0 = self.dens.alphaviscosity*self.dens.aspectratio*self.dens.aspectratio
    	self.tvisc = self.dens.ymax**2/self.nu(self.dens.ymax)
        self.rh = (self.mp/3)**(1./3) * self.a
    	self.omega = zeros(self.vp.data.shape)
    	self.vpc = zeros(self.vp.data.shape)
    	for i in range(len(self.r)):
    		self.omega[i,:] = self.vp.data[i,:]/self.r[i] + self.a**(-1.5)
    		self.vpc[i,:] = self.omega[i,:]*self.r[i]

        self.dTr = (-2*pi*self.dens.data * self.dp_potential(self.rr,self.pp))
        self.dTr_mean = self.dTr.mean(axis=1)

    def nu(self,r):
    	return self.nu0 * r**(2*self.dens.flaringindex+.5)
    def scaleH(self,r):
        return self.dens.aspectratio * r**(self.dens.flaringindex + 1)

    def potential(self,r,phi):
        if self.dens.rochesmoothing:
            smoothing = self.dens.thicknesssmoothing*self.scaleH(r)
        else:
            smoothing = self.dens.thicknesssmoothing*self.rh
        smoothing *= smoothing
        rad = r**2 + smoothing + self.a**2 - 2*r*self.a*cos(phi)
        rad = sqrt(rad)
        return -self.mp/rad
    def dp_potential(self,r,phi):
        if self.dens.rochesmoothing:
            smoothing = self.dens.thicknesssmoothing*self.scaleH(r)
        else:
            smoothing = self.dens.thicknesssmoothing*self.rh
        smoothing *= smoothing
        rad = r**2 + self.a**2 - 2*r*self.a*cos(phi) + self.dens.thicknesssmoothing*self.scaleH(r)**2
        rad = rad**(1.5)
        return self.mp * r*self.a*sin(phi)/rad


    def animate_mdot(self,irange,logspacing=True):

        self.__init__(0,directory=self.directory)
        temp0,temp1,temp2,temp3,temp4,temp5,temp6  = self.calculate_total_mdot(logspacing)
        dens_list = [log10(temp0)]
        divm_list = [temp1]
        divm_mean_list = [temp2]
        drmx_mean_list = [temp3]
        dpmy_mean_list = [temp4]
        mx_mean_list = [temp5]
        my_mean_list = [temp6]
        t = [self.t]

        fig,axes = subplots(1,4,figsize=(15,10))

        axes[3].set_title('t = %.1f' % t[0],color='w')
        line0=axes[0].imshow(dens_list[0],aspect='auto',origin='lower')
        line1=axes[1].imshow(divm_list[0],aspect='auto',origin='lower');

        line2,=axes[2].plot(self.r,divm_mean_list[0],label='< $ \\nabla \\cdot \\dot{M} = \\dot{\\Sigma}$ >')
        line3,=axes[2].plot(self.r,drmx_mean_list[0],label='<$\\nabla_r (r \\Sigma v_r)$>')
        line4,=axes[2].plot(self.r,dpmy_mean_list[0],label='-<$\\nabla_\\phi(\\Sigma v_\\phi)$>')
        axes[2].legend(loc='best')
        axes[2].set_ylim(-.0001,.0001)
        line5,=axes[3].plot(self.r,mx_mean_list[0],label='<-$r \\Sigma v_r$>')
        line6,=axes[3].plot(self.r,my_mean_list[0],label='-<$\\Sigma v_\\phi$>')
        axes[3].legend(loc='best')

        for i in irange:
            if i > 0:
                self.__init__(i,directory=self.directory)
                t.append(self.t)
                temp0,temp1,temp2,temp3,temp4,temp5,temp6  = self.calculate_total_mdot(logspacing)
                dens_list.append(log10(temp0))
                divm_list.append(temp1)
                divm_mean_list.append(temp2)
                drmx_mean_list.append(temp3)
                dpmy_mean_list.append(temp4)
                mx_mean_list.append(temp5)
                my_mean_list.append(temp6)


        for i in range(len(t)):
            axes[2].set_title('t = %.1f' % t[i])
            line0.set_data(dens_list[i])
            line1.set_data(divm_list[i])
            line2.set_ydata(divm_mean_list[i])
            line3.set_ydata(drmx_mean_list[i])
            line4.set_ydata(dpmy_mean_list[i])
            line5.set_ydata(mx_mean_list[i])
            line6.set_ydata(my_mean_list[i])
            fig.canvas.draw()


    def plot_mdot(self,logspacing=True):
        dens,divm, divm_mean, drmx_mean, dpmy_mean, mx_mean,my_mean = self.calculate_total_mdot(logspacing)
        fig,axes = subplots(1,4,figsize=(15,10))

        axes[0].imshow(log10(dens),aspect='auto',origin='lower');
        axes[1].imshow(divm,aspect='auto',origin='lower');
        axes[3].plot(self.r,mx_mean,self.r,my_mean)
        axes[3].legend(['<-$r \\Sigma v_r$>', '-<$\\Sigma v_\\phi$>'],loc='best')
        axes[2].plot(self.r,divm_mean,self.r,drmx_mean,self.r,dpmy_mean)
        axes[2].legend(['< $ \\nabla \\cdot \\dot{M} = \\dot{\\Sigma}$ >', '-<$\\nabla_r (r \\Sigma v_r)$>', '-<$\\nabla_\\phi(\\Sigma v_\\phi)$>'],loc='best')


    def calculate_total_mdot(self,logspacing=True):
        pp,rr = meshgrid(self.phi,self.r)
        if logspacing:
            dr = diff(log(self.r))[0]
            norm = rr
        else:
            dr = diff(self.r)[0]
            norm = 1

        dp = diff(self.phi)[0]

        dens = self.dens.data
        vp = self.vpc
        vr = self.vr.data

        mx = -rr * dens*vr
        my = -dens*vp

        drmx,_ = gradient(mx,dlr,dp)
        _,dpmy = gradient(my,dlr,dp)

        drmx /= (rr*norm)
        dpmy /= rr
        divm = drmx + dpmy
        return dens,divm, divm.mean(axis=1),drmx.mean(axis=1), dpmy.mean(axis=1), mx.mean(axis=1), my.mean(axis=1)



    def summary(self):
    	fig,(axd,axm,axv) = subplots(3,1,sharex='col')
    	lined, = axd.plot(self.r,self.dbar,linewidth=3)
    	linev, = axv.plot(self.r,self.vrbar,linewidth=3)
    	linem, = axm.plot(self.r,self.mdot,linewidth=3)
    	axv.set_xlabel('$r$',fontsize=20)
    	axv.set_ylabel('$<v_r$>',fontsize=20)
    	axd.set_ylabel('$<\\Sigma>$',fontsize=20)
    	axm.set_ylabel('<$\\dot{M}$>',fontsize=20)
    #		axm.set_yscale('symlog',linthreshy=1e-7)
    #		axv.set_ylim(-.001,.001)
    	axd.set_title('t = %.1f = %.1f P = %.1f t_visc' % (self.t,self.t/(2*pi*self.a**(1.5)),self.t/self.tvisc))
    def streams(self,rlims=None,plims=None,ax=None,**kargs):
        if ax == None:
            fig=figure()
            ax=fig.add_subplot(111)

        if rlims == None:
            rlims = (self.dens.ymin,self.dens.ymax)
        if plims == None:
            plims = (-pi,pi)

        rinds = (self.r<=rlims[1])&(self.r>=rlims[0])
        pinds = (self.phi<=plims[1])&(self.phi>=plims[0])

        vr = self.vr.data[rinds,:][:,pinds]
        vp = self.vp.data[rinds,:][:,pinds]
        dens = self.dens.data[rinds,:][:,pinds]
        rr,pp = meshgrid(self.r[rinds],self.phi[pinds])
        line2d= ax.pcolormesh(log(rr),pp,log10(dens.transpose()))
        cbar = colorbar(line2d,ax=ax)
        cbar.set_label('$\\log_{10}{\\Sigma}$',fontsize=20)
        ax.streamplot(log(self.r[rinds]),self.phi[pinds],vr.transpose(),vp.transpose(),**kargs)
        ax.set_xlabel('$\ln(r)$',fontsize=20)
        ax.set_ylabel('$\phi$',fontsize=20)
        ax.set_ylim(plims)
        ax.set_xlim((log(rlims[0]),log(rlims[1])))
        ax.axvline(log(self.a+2*self.rh),color='k',linewidth=3)
        ax.axvline(log(self.a-2*self.rh),color='k',linewidth=3)


        rh,ph = self.calculate_circle(self.rh,rlims)
        rs,ps = self.calculate_circle(self.dens.aspectratio*self.a*self.dens.thicknesssmoothing,rlims)

        ax.plot(log(rh), ph,'-w',linewidth=3)
        ax.plot(log(rs),ps,'--w',linewidth=3)
    def calculate_circle(self,d,rlims=None):
        if rlims == None:
            rlims = (self.dens.ymin,self.dens.ymax)

        rinds = (self.r>=rlims[0])&(self.r<=rlims[1])
        newinds = (self.r[rinds]>=(self.a-d))&(self.r[rinds]<=(self.a+d))
        circle_upper = arccos((self.r[rinds][newinds]**2 + self.a**2 - d**2)/(2*self.a*self.r[rinds][newinds]))
        circle_lower=  -circle_upper

        circle_r = hstack((array([self.a-d]),self.r[rinds][newinds]))
        circle_r = hstack((circle_r,array([self.a+d])))
        circle_r = hstack((circle_r,self.r[rinds][newinds][::-1]))

        circle_phi = hstack((array([0]),circle_upper))
        circle_phi = hstack((circle_phi,array([0])))
        circle_phi = hstack((circle_phi,circle_lower[::-1]))
        return circle_r, circle_phi

    def stagnation(self,r,phi,rr,pp,vr,vp,dens):
    	ivp = interp2d(rr,pp,vp)
    	ivr = interp2d(rr,pp,vr)
    	vp_func = lambda x,y: ivp(x,y)
    	vr_func = lambda x,y: ivr(x,y)
    	root(ivp)



    def animate(self,irange):
    	fig,(axd,axm,axv) = subplots(3,1,sharex='col')
    	lined, = axd.plot(self.r,self.dbar,linewidth=3)
    	linev, = axv.plot(self.r,self.vrbar,linewidth=3)
    	linem, = axm.plot(self.r,self.mdot,linewidth=3)

        s0 = Sim(0,directory=self.directory)
        axd.plot(s0.r,s0.dbar,'--k',linewidth=3)
    	axv.plot(s0.r,s0.vrbar,'--k',linewidth=3)
    	axm.plot(s0.r,s0.mdot,'--k',linewidth=3)
        axv.set_ylim(-.05,.05)

    	axv.set_xlabel('$r$',fontsize=20)
    	axv.set_ylabel('<$v_r$>',fontsize=20)
    	axd.set_ylabel('<$\\Sigma$>',fontsize=20)
    	axm.set_ylabel('<$\\dot{M}$>',fontsize=20)

    	dbar = zeros((len(self.dbar),len(irange)))
    	vrbar = zeros((len(self.dbar),len(irange)))
    	mdot = zeros(dbar.shape)
    	t = zeros(len(irange))
    	a = zeros(len(irange))

    	t[0] = self.t
    	a[0] = self.a
    	dbar[:,0] = self.dbar
    	vrbar[:,0] = self.vrbar
    	mdot[:,0] = self.mdot
    	for i,j in enumerate(irange[1:],start=1):
            self.__init__(j,directory=self.directory)
            t[i] = self.t
            dbar[:,i] = self.dbar
            vrbar[:,i] = self.vrbar
            mdot[:,i] = self.mdot
            a[i] = self.a

#    	axv.set_ylim(-.001,.001)
#    	axm.set_ylim((mdot.min(),mdot.max()))
#    	axm.set_yscale('symlog',linthreshy=1e-8)
    	for i in range(len(t)):
    		lined.set_ydata(dbar[:,i])
    		linev.set_ydata(vrbar[:,i])
    		linem.set_ydata(mdot[:,i])
    		axd.set_title('t = %.1f = %.1f P' % (t[i],t[i]/(2*pi*a[i]**(1.5))),color='w')
    		fig.canvas.draw()
