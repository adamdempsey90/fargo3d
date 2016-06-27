import numpy as np
import sys
from matplotlib.mlab import griddata
import multiprocessing as mp


def write_vtk_file((n,option)):
    print 'Converting number %d'%n
    ifile_dens = "gasdens" + "{:d}".format(n) + '.dat'
    ifile_ener = "gasenergy" + "{:d}".format(n) + '.dat'
    ifile_vx   = "gasvx" + "{:d}".format(n) + '.dat'
    ifile_vy   = "gasvy" + "{:d}".format(n) + '.dat'
    ifile_vz   = "gasvz" + "{:d}".format(n) + '.dat'
    #ifile_dens = "Density" + "{:d}".format(n) + '.dat'
    #ifile_ener = "Energy" + "{:d}".format(n) + '.dat'
    #ifile_vx   = "Vx" + "{:d}".format(n) + '.dat'
    #ifile_vy   = "Vy" + "{:d}".format(n) + '.dat'
    #ifile_vz   = "Vz" + "{:d}".format(n) + '.dat'

    density = np.fromfile(ifile_dens, "d")
    energy  = np.fromfile(ifile_ener, "d")
    vx      = np.fromfile(ifile_vx  , "d")
    vy      = np.fromfile(ifile_vy  , "d")
    vz      = np.fromfile(ifile_vz  , "d")


    try:
        domain_x = open("domain_x.dat",'r')
    except:
        print 'Error!!!',domain_x.dat, 'cannot be opened.'
    try:
        domain_y = open("domain_y.dat",'r')
    except:
        print 'Error!!!',domain_y.dat, 'cannot be opened.'
    try:
        domain_z = open("domain_z.dat",'r')
    except:
        print 'Error!!!',domain_y.dat, 'cannot be opened.'
    try:
        dims = open("dimensions.dat",'r')
    except:
        print 'Error!!!',dimensions.dat, 'cannot be opened.'

    DIMS = dims.readlines() #Read all dims file
    DIMS = DIMS[1].split()

    DOMAIN_X = domain_x.readlines()
    DOMAIN_Y = domain_y.readlines()
    DOMAIN_Z = domain_z.readlines()

    domain_x.close()
    domain_y.close()
    domain_z.close()
    dims.close()

    NGHY = 3
    NGHZ = 3

    Nx = int(DIMS[6])#+1 # Add a point to avoid border discontinuity
    Ny = int(DIMS[7])
    Nz = int(DIMS[8])

    x_arr = np.ndarray([Nx], dtype = float)
    y_arr = np.ndarray([Ny], dtype = float)
    z_arr = np.ndarray([Nz], dtype = float)
    mesh  = np.ndarray([3*Nx*Ny*Nz], dtype = float)
    var1  = np.ndarray([Nx*Ny*Nz], dtype = float)

    #for i in range(Nx-1):
    for i in range(Nx):
        x_arr[i] = float(DOMAIN_X[i])
    #x_arr[Nx-1] = x_arr[0]  #Periodic boundary
    for i in range(Ny):
        y_arr[i] = float(DOMAIN_Y[i+NGHY])
    for i in range(Nz):
        z_arr[i] = float(DOMAIN_Z[i+NGHZ])


    l = 0
    for k in range(Nz):
        for j in range(Ny):
            for i in range(Nx):
                if(option == 0):
                    mesh[l] = np.cos(x_arr[i])*np.sin(z_arr[k])*y_arr[j]
                    mesh[l+1] = np.sin(x_arr[i])*np.sin(z_arr[k])*y_arr[j]
                    mesh[l+2] = y_arr[j]*np.cos(z_arr[k])
                if(option == 1):
                    mesh[l] = np.cos(x_arr[i])*y_arr[j];
                    mesh[l+1] = np.sin(x_arr[i])*y_arr[j];
                    mesh[l+2] = z_arr[k];
                l+=3;

    #data1  = np.ndarray([(Nx-1)*(Ny-1)*(Nz-1)], dtype = float)
    #l=0
    #for k in range(Nz-1):
    #    for j in range(Ny-1):
    #        for i in range(Nx-1):
    #            if(i < Nx-1):
    #                data1[l] = data[l];
    #                l+=1;
    #            else:
    #                data1[l] = data[l-(Nx-1)]
    #
    #data1 = tuple(data1)
    #var = (("Density", 1, 0, data1),) #NODAL

    #Mesh construction--------------------------------------
    connectivity = [] #Unstructured Mesh!, avoid periodic problem!
    stride = Nx*Ny
    for k in range(Nz-1):
        for j in range(Ny-1):
            for i in range(Nx):
                l = i+j*Nx+k*stride
                if(i<Nx-1):
                    lxp = l+1
                else:
                    lxp = l-(Nx-1) #WARNINR Nx-1?
                connectivity += [[12,l,lxp,lxp+Nx,l+Nx,
                                  l+stride, lxp+stride,
                                  lxp+Nx+stride, l+Nx+stride]]
    #Mesh construction--------------------------------------
    #print connectivity[0], connectivity[Nx-1]

    dims = (Nx,Ny,Nz)
    var = (("Density", 1, 1, tuple(density)),
           ("Energy" , 1, 1, tuple(energy )),
           ("Vtheta" , 1, 1, tuple(vx     )),
           ("Vrad"   , 1, 1, tuple(vy     )),
           ("Vpolar" , 1, 1, tuple(vz     ))) #ZONAL

    #vw.WriteCurvilinearMesh("vwcurv3d.vtk", 1, dims, tuple(mesh), var)
    vw.WriteUnstructuredMesh("vwucd3d{:d}.vtk".format(n), 1, tuple(mesh),connectivity, var)


if __name__ == "__main__":
    try:
        import visit_writer as vw
    except ImportError:
        sys.path.append('/home/amd616/lib/python')
        import visit_writer as vw

    nstart = int(sys.argv[1])  #input from term
    nend = int(sys.argv[2]) + 1  #input from term
    option = int(sys.argv[3]) #1 = Cilindrical, 0 = Spherical

    print 'Converting outputs %d to %d'%(nstart,nend-1)
    if len(sys.argv) > 3:
        numprocs = int(sys.argv[4])
        print 'Using %d procs'%numprocs
        p = mp.Pool(numprocs)
        args = [(i,option) for i in range(nstart,nend)]
        p.map(write_vtk_file,args)
    else:
        args = [(i,option) for i in range(nstart,nend)]
        for a in args:
            write_vtk_file(a)

