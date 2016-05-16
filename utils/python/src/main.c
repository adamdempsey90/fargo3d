#include "evolve.h"

int main(int argc, char *argv[]) {

    int n;
    int status;
    char directory[256],outputdir[256];
    char param_fname[100];

    if (argc > 4) {
        printf("Too many arguments.Only using first 3.\n");
    }
    n = atoi(argv[1]);
    time_step = atof(argv[2]);
    strcpy(directory,argv[3]);

    nb = 2;

    sprintf(outputdir,"%stemp_files/",directory);
    status = mkdir(outputdir,S_IRWXU | S_IRWXG | S_IRWXO);
    if (status == -33)  {
        printf("Status is -33 for some reason.\n");
    }
    //sprintf(param_fname,"%sparam_file.txt",directory);


    printf("Reading\n");
    read_param_file(directory);
    nx = params.nx;
    ny = params.ny;
    nz = params.nz;
    CFL = params.cfl;
    IndirectTerm = params.indirect;

    size_x = nx;
    size_y = ny+2*NGHY;
    size_z = nz;
    stride = size_x*size_y;
    pitch = size_x;
    dx = 2*M_PI/nx;

    allocate_all();
    read_domain(directory);
    read_files(n,directory);
    //read_single_file(n,1,directory);
    /*
    if (cfl_dt < dt) {
        printf("Using cfl limited timestep of %.4e instead of %.4e\n",cfl_dt,dt);
        dt = cfl_dt;
    }
    */

    printf("time step = %lg\n",time_step);
    //move_to_com();
   output_init(outputdir);
   printf("output initial.\n");
   init_rk5();
    nsteps = 0;
    double tstart = psys[0].t;
    double tend = psys[0].t + time_step;
    double time = tstart;
    printf("Starting time = %lg\nEnding time = %lg\n",tstart,tend);
    while ((time < tend) && nsteps<MAXSTEPS) {
        dt = cfl();
        if (dt <= MINDT) {
            printf("Timestep has fallen below minimum!\n");
            break;
        }
        if (time + dt > tend) {
            dt = tend-time;
        }
        printf("Time %lg dt = %lg\n",time,dt);
        set_bc();
        set_avg(0);
    
        potential();
        move_planet();
    
       source_step();
       
        set_Lamex();
       viscosity();
#ifdef ARTIFICIALVISCOSITY
        artificial_visc();
#endif
    
        temp_to_vel();
        set_bc();
        
        vel_to_temp();
        
        transport_step();
      
        set_avg(1);

        set_Lamdep();

        stockholm();

        time += dt;
        printf("Advancing by %lg to %lg\n",dt,time);
        nsteps++;
    }
//    temp_to_vel();   

    dt = tend - tstart;
    int j;
    for(j=0;j<size_y;j++) {
        dbart[j]/=dt;
        //Lt[j]/=dt;
        //Lt[j+size_y]/=dt;
        //Ld[j]/=dt;
        //Ld[j+size_y]/=dt;
        //Lw[j]/=dt;
        //Lw[j+size_y]/=dt;
        Lw[j] = Lt[j] - Ld[j];
        Lw[j+size_y] = Lt[j+size_y] - Ld[j+size_y];
        drFt[j]/=dt;
        drFd[j]/=dt;
        drFw[j]/=dt;
        Lamex[j]/=dt;
        Lamex[j+size_y]/=dt;
        Lamdep[j]/=dt;
        dtLt[j]/=dt;
        dtLd[j]/=dt;
        dtLw[j]/=dt;


    }

    output(outputdir);
    output_torque(outputdir);
    free_rk5();
    //free_all();
    return 0;
}
