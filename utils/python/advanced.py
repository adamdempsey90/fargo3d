
viridis = load_viridis()
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

    def plot(self, log=False, cartesian=False, cmap=viridis, **karg):
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


    def animate_mdot(self,irange,logspacing=True,cmap=viridis,**kargs):

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
        line0=axes[0].imshow(dens_list[0].transpose(),aspect='auto',origin='lower',cmap=cmap,**kargs)
        line1=axes[1].imshow(divm_list[0].transpose(),aspect='auto',origin='lower',cmap=cmap,**kargs);

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
                dens_list.append(log10(temp0).transpose())
                divm_list.append(temp1.transpose())
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

        cmap = kargs.pop('cmap',viridis)

        line2d= ax.pcolormesh(log(rr),pp,log10(dens.transpose()),cmap=cmap)
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






def load_viridis():
    from matplotlib.colors import LinearSegmentedColormap
# Used to reconstruct the colormap in viscm
    parameters = {'xp': [22.674387857633945, 11.221508276482126, -14.356589454756971, -47.18817758739222, -34.59001004812521, -6.0516291196352654],
                  'yp': [-20.102530541012214, -33.08246073298429, -42.24476439790574, -5.595549738219887, 42.5065445026178, 40.13395157135497],
                  'min_JK': 18.8671875,
                  'max_JK': 92.5}

    cm_data = [[ 0.26700401,  0.00487433,  0.32941519],
           [ 0.26851048,  0.00960483,  0.33542652],
           [ 0.26994384,  0.01462494,  0.34137895],
           [ 0.27130489,  0.01994186,  0.34726862],
           [ 0.27259384,  0.02556309,  0.35309303],
           [ 0.27380934,  0.03149748,  0.35885256],
           [ 0.27495242,  0.03775181,  0.36454323],
           [ 0.27602238,  0.04416723,  0.37016418],
           [ 0.2770184 ,  0.05034437,  0.37571452],
           [ 0.27794143,  0.05632444,  0.38119074],
           [ 0.27879067,  0.06214536,  0.38659204],
           [ 0.2795655 ,  0.06783587,  0.39191723],
           [ 0.28026658,  0.07341724,  0.39716349],
           [ 0.28089358,  0.07890703,  0.40232944],
           [ 0.28144581,  0.0843197 ,  0.40741404],
           [ 0.28192358,  0.08966622,  0.41241521],
           [ 0.28232739,  0.09495545,  0.41733086],
           [ 0.28265633,  0.10019576,  0.42216032],
           [ 0.28291049,  0.10539345,  0.42690202],
           [ 0.28309095,  0.11055307,  0.43155375],
           [ 0.28319704,  0.11567966,  0.43611482],
           [ 0.28322882,  0.12077701,  0.44058404],
           [ 0.28318684,  0.12584799,  0.44496   ],
           [ 0.283072  ,  0.13089477,  0.44924127],
           [ 0.28288389,  0.13592005,  0.45342734],
           [ 0.28262297,  0.14092556,  0.45751726],
           [ 0.28229037,  0.14591233,  0.46150995],
           [ 0.28188676,  0.15088147,  0.46540474],
           [ 0.28141228,  0.15583425,  0.46920128],
           [ 0.28086773,  0.16077132,  0.47289909],
           [ 0.28025468,  0.16569272,  0.47649762],
           [ 0.27957399,  0.17059884,  0.47999675],
           [ 0.27882618,  0.1754902 ,  0.48339654],
           [ 0.27801236,  0.18036684,  0.48669702],
           [ 0.27713437,  0.18522836,  0.48989831],
           [ 0.27619376,  0.19007447,  0.49300074],
           [ 0.27519116,  0.1949054 ,  0.49600488],
           [ 0.27412802,  0.19972086,  0.49891131],
           [ 0.27300596,  0.20452049,  0.50172076],
           [ 0.27182812,  0.20930306,  0.50443413],
           [ 0.27059473,  0.21406899,  0.50705243],
           [ 0.26930756,  0.21881782,  0.50957678],
           [ 0.26796846,  0.22354911,  0.5120084 ],
           [ 0.26657984,  0.2282621 ,  0.5143487 ],
           [ 0.2651445 ,  0.23295593,  0.5165993 ],
           [ 0.2636632 ,  0.23763078,  0.51876163],
           [ 0.26213801,  0.24228619,  0.52083736],
           [ 0.26057103,  0.2469217 ,  0.52282822],
           [ 0.25896451,  0.25153685,  0.52473609],
           [ 0.25732244,  0.2561304 ,  0.52656332],
           [ 0.25564519,  0.26070284,  0.52831152],
           [ 0.25393498,  0.26525384,  0.52998273],
           [ 0.25219404,  0.26978306,  0.53157905],
           [ 0.25042462,  0.27429024,  0.53310261],
           [ 0.24862899,  0.27877509,  0.53455561],
           [ 0.2468114 ,  0.28323662,  0.53594093],
           [ 0.24497208,  0.28767547,  0.53726018],
           [ 0.24311324,  0.29209154,  0.53851561],
           [ 0.24123708,  0.29648471,  0.53970946],
           [ 0.23934575,  0.30085494,  0.54084398],
           [ 0.23744138,  0.30520222,  0.5419214 ],
           [ 0.23552606,  0.30952657,  0.54294396],
           [ 0.23360277,  0.31382773,  0.54391424],
           [ 0.2316735 ,  0.3181058 ,  0.54483444],
           [ 0.22973926,  0.32236127,  0.54570633],
           [ 0.22780192,  0.32659432,  0.546532  ],
           [ 0.2258633 ,  0.33080515,  0.54731353],
           [ 0.22392515,  0.334994  ,  0.54805291],
           [ 0.22198915,  0.33916114,  0.54875211],
           [ 0.22005691,  0.34330688,  0.54941304],
           [ 0.21812995,  0.34743154,  0.55003755],
           [ 0.21620971,  0.35153548,  0.55062743],
           [ 0.21429757,  0.35561907,  0.5511844 ],
           [ 0.21239477,  0.35968273,  0.55171011],
           [ 0.2105031 ,  0.36372671,  0.55220646],
           [ 0.20862342,  0.36775151,  0.55267486],
           [ 0.20675628,  0.37175775,  0.55311653],
           [ 0.20490257,  0.37574589,  0.55353282],
           [ 0.20306309,  0.37971644,  0.55392505],
           [ 0.20123854,  0.38366989,  0.55429441],
           [ 0.1994295 ,  0.38760678,  0.55464205],
           [ 0.1976365 ,  0.39152762,  0.55496905],
           [ 0.19585993,  0.39543297,  0.55527637],
           [ 0.19410009,  0.39932336,  0.55556494],
           [ 0.19235719,  0.40319934,  0.55583559],
           [ 0.19063135,  0.40706148,  0.55608907],
           [ 0.18892259,  0.41091033,  0.55632606],
           [ 0.18723083,  0.41474645,  0.55654717],
           [ 0.18555593,  0.4185704 ,  0.55675292],
           [ 0.18389763,  0.42238275,  0.55694377],
           [ 0.18225561,  0.42618405,  0.5571201 ],
           [ 0.18062949,  0.42997486,  0.55728221],
           [ 0.17901879,  0.43375572,  0.55743035],
           [ 0.17742298,  0.4375272 ,  0.55756466],
           [ 0.17584148,  0.44128981,  0.55768526],
           [ 0.17427363,  0.4450441 ,  0.55779216],
           [ 0.17271876,  0.4487906 ,  0.55788532],
           [ 0.17117615,  0.4525298 ,  0.55796464],
           [ 0.16964573,  0.45626209,  0.55803034],
           [ 0.16812641,  0.45998802,  0.55808199],
           [ 0.1666171 ,  0.46370813,  0.55811913],
           [ 0.16511703,  0.4674229 ,  0.55814141],
           [ 0.16362543,  0.47113278,  0.55814842],
           [ 0.16214155,  0.47483821,  0.55813967],
           [ 0.16066467,  0.47853961,  0.55811466],
           [ 0.15919413,  0.4822374 ,  0.5580728 ],
           [ 0.15772933,  0.48593197,  0.55801347],
           [ 0.15626973,  0.4896237 ,  0.557936  ],
           [ 0.15481488,  0.49331293,  0.55783967],
           [ 0.15336445,  0.49700003,  0.55772371],
           [ 0.1519182 ,  0.50068529,  0.55758733],
           [ 0.15047605,  0.50436904,  0.55742968],
           [ 0.14903918,  0.50805136,  0.5572505 ],
           [ 0.14760731,  0.51173263,  0.55704861],
           [ 0.14618026,  0.51541316,  0.55682271],
           [ 0.14475863,  0.51909319,  0.55657181],
           [ 0.14334327,  0.52277292,  0.55629491],
           [ 0.14193527,  0.52645254,  0.55599097],
           [ 0.14053599,  0.53013219,  0.55565893],
           [ 0.13914708,  0.53381201,  0.55529773],
           [ 0.13777048,  0.53749213,  0.55490625],
           [ 0.1364085 ,  0.54117264,  0.55448339],
           [ 0.13506561,  0.54485335,  0.55402906],
           [ 0.13374299,  0.54853458,  0.55354108],
           [ 0.13244401,  0.55221637,  0.55301828],
           [ 0.13117249,  0.55589872,  0.55245948],
           [ 0.1299327 ,  0.55958162,  0.55186354],
           [ 0.12872938,  0.56326503,  0.55122927],
           [ 0.12756771,  0.56694891,  0.55055551],
           [ 0.12645338,  0.57063316,  0.5498411 ],
           [ 0.12539383,  0.57431754,  0.54908564],
           [ 0.12439474,  0.57800205,  0.5482874 ],
           [ 0.12346281,  0.58168661,  0.54744498],
           [ 0.12260562,  0.58537105,  0.54655722],
           [ 0.12183122,  0.58905521,  0.54562298],
           [ 0.12114807,  0.59273889,  0.54464114],
           [ 0.12056501,  0.59642187,  0.54361058],
           [ 0.12009154,  0.60010387,  0.54253043],
           [ 0.11973756,  0.60378459,  0.54139999],
           [ 0.11951163,  0.60746388,  0.54021751],
           [ 0.11942341,  0.61114146,  0.53898192],
           [ 0.11948255,  0.61481702,  0.53769219],
           [ 0.11969858,  0.61849025,  0.53634733],
           [ 0.12008079,  0.62216081,  0.53494633],
           [ 0.12063824,  0.62582833,  0.53348834],
           [ 0.12137972,  0.62949242,  0.53197275],
           [ 0.12231244,  0.63315277,  0.53039808],
           [ 0.12344358,  0.63680899,  0.52876343],
           [ 0.12477953,  0.64046069,  0.52706792],
           [ 0.12632581,  0.64410744,  0.52531069],
           [ 0.12808703,  0.64774881,  0.52349092],
           [ 0.13006688,  0.65138436,  0.52160791],
           [ 0.13226797,  0.65501363,  0.51966086],
           [ 0.13469183,  0.65863619,  0.5176488 ],
           [ 0.13733921,  0.66225157,  0.51557101],
           [ 0.14020991,  0.66585927,  0.5134268 ],
           [ 0.14330291,  0.66945881,  0.51121549],
           [ 0.1466164 ,  0.67304968,  0.50893644],
           [ 0.15014782,  0.67663139,  0.5065889 ],
           [ 0.15389405,  0.68020343,  0.50417217],
           [ 0.15785146,  0.68376525,  0.50168574],
           [ 0.16201598,  0.68731632,  0.49912906],
           [ 0.1663832 ,  0.69085611,  0.49650163],
           [ 0.1709484 ,  0.69438405,  0.49380294],
           [ 0.17570671,  0.6978996 ,  0.49103252],
           [ 0.18065314,  0.70140222,  0.48818938],
           [ 0.18578266,  0.70489133,  0.48527326],
           [ 0.19109018,  0.70836635,  0.48228395],
           [ 0.19657063,  0.71182668,  0.47922108],
           [ 0.20221902,  0.71527175,  0.47608431],
           [ 0.20803045,  0.71870095,  0.4728733 ],
           [ 0.21400015,  0.72211371,  0.46958774],
           [ 0.22012381,  0.72550945,  0.46622638],
           [ 0.2263969 ,  0.72888753,  0.46278934],
           [ 0.23281498,  0.73224735,  0.45927675],
           [ 0.2393739 ,  0.73558828,  0.45568838],
           [ 0.24606968,  0.73890972,  0.45202405],
           [ 0.25289851,  0.74221104,  0.44828355],
           [ 0.25985676,  0.74549162,  0.44446673],
           [ 0.26694127,  0.74875084,  0.44057284],
           [ 0.27414922,  0.75198807,  0.4366009 ],
           [ 0.28147681,  0.75520266,  0.43255207],
           [ 0.28892102,  0.75839399,  0.42842626],
           [ 0.29647899,  0.76156142,  0.42422341],
           [ 0.30414796,  0.76470433,  0.41994346],
           [ 0.31192534,  0.76782207,  0.41558638],
           [ 0.3198086 ,  0.77091403,  0.41115215],
           [ 0.3277958 ,  0.77397953,  0.40664011],
           [ 0.33588539,  0.7770179 ,  0.40204917],
           [ 0.34407411,  0.78002855,  0.39738103],
           [ 0.35235985,  0.78301086,  0.39263579],
           [ 0.36074053,  0.78596419,  0.38781353],
           [ 0.3692142 ,  0.78888793,  0.38291438],
           [ 0.37777892,  0.79178146,  0.3779385 ],
           [ 0.38643282,  0.79464415,  0.37288606],
           [ 0.39517408,  0.79747541,  0.36775726],
           [ 0.40400101,  0.80027461,  0.36255223],
           [ 0.4129135 ,  0.80304099,  0.35726893],
           [ 0.42190813,  0.80577412,  0.35191009],
           [ 0.43098317,  0.80847343,  0.34647607],
           [ 0.44013691,  0.81113836,  0.3409673 ],
           [ 0.44936763,  0.81376835,  0.33538426],
           [ 0.45867362,  0.81636288,  0.32972749],
           [ 0.46805314,  0.81892143,  0.32399761],
           [ 0.47750446,  0.82144351,  0.31819529],
           [ 0.4870258 ,  0.82392862,  0.31232133],
           [ 0.49661536,  0.82637633,  0.30637661],
           [ 0.5062713 ,  0.82878621,  0.30036211],
           [ 0.51599182,  0.83115784,  0.29427888],
           [ 0.52577622,  0.83349064,  0.2881265 ],
           [ 0.5356211 ,  0.83578452,  0.28190832],
           [ 0.5455244 ,  0.83803918,  0.27562602],
           [ 0.55548397,  0.84025437,  0.26928147],
           [ 0.5654976 ,  0.8424299 ,  0.26287683],
           [ 0.57556297,  0.84456561,  0.25641457],
           [ 0.58567772,  0.84666139,  0.24989748],
           [ 0.59583934,  0.84871722,  0.24332878],
           [ 0.60604528,  0.8507331 ,  0.23671214],
           [ 0.61629283,  0.85270912,  0.23005179],
           [ 0.62657923,  0.85464543,  0.22335258],
           [ 0.63690157,  0.85654226,  0.21662012],
           [ 0.64725685,  0.85839991,  0.20986086],
           [ 0.65764197,  0.86021878,  0.20308229],
           [ 0.66805369,  0.86199932,  0.19629307],
           [ 0.67848868,  0.86374211,  0.18950326],
           [ 0.68894351,  0.86544779,  0.18272455],
           [ 0.69941463,  0.86711711,  0.17597055],
           [ 0.70989842,  0.86875092,  0.16925712],
           [ 0.72039115,  0.87035015,  0.16260273],
           [ 0.73088902,  0.87191584,  0.15602894],
           [ 0.74138803,  0.87344918,  0.14956101],
           [ 0.75188414,  0.87495143,  0.14322828],
           [ 0.76237342,  0.87642392,  0.13706449],
           [ 0.77285183,  0.87786808,  0.13110864],
           [ 0.78331535,  0.87928545,  0.12540538],
           [ 0.79375994,  0.88067763,  0.12000532],
           [ 0.80418159,  0.88204632,  0.11496505],
           [ 0.81457634,  0.88339329,  0.11034678],
           [ 0.82494028,  0.88472036,  0.10621724],
           [ 0.83526959,  0.88602943,  0.1026459 ],
           [ 0.84556056,  0.88732243,  0.09970219],
           [ 0.8558096 ,  0.88860134,  0.09745186],
           [ 0.86601325,  0.88986815,  0.09595277],
           [ 0.87616824,  0.89112487,  0.09525046],
           [ 0.88627146,  0.89237353,  0.09537439],
           [ 0.89632002,  0.89361614,  0.09633538],
           [ 0.90631121,  0.89485467,  0.09812496],
           [ 0.91624212,  0.89609127,  0.1007168 ],
           [ 0.92610579,  0.89732977,  0.10407067],
           [ 0.93590444,  0.8985704 ,  0.10813094],
           [ 0.94563626,  0.899815  ,  0.11283773],
           [ 0.95529972,  0.90106534,  0.11812832],
           [ 0.96489353,  0.90232311,  0.12394051],
           [ 0.97441665,  0.90358991,  0.13021494],
           [ 0.98386829,  0.90486726,  0.13689671],
           [ 0.99324789,  0.90615657,  0.1439362 ]]

    return LinearSegmentedColormap.from_list('viridis', cm_data)
