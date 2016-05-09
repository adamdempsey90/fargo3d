#include "evolve.h"

int main(int argc, char *argv[]) {

    int n;
    int status;
    double cfl_dt;
    char directory[256],outputdir[256];
    char param_fname[100];

    if (argc > 4) {
        printf("Too many arguments.Only using first 3.\n");
    }
    n = atoi(argv[1]);
    nsteps = atoi(argv[2]);
    strcpy(directory,argv[3]);

    nb = 2;

    sprintf(outputdir,"%stemp_files/",directory);
    status = mkdir(outputdir,S_IRWXU | S_IRWXG | S_IRWXO);
    if (status == -33)  {
        printf("Status is -33 for some reason.\n");
    }
    sprintf(param_fname,"%sparam_file.txt",directory);


    read_param_file(param_fname);

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
    cfl_dt = cfl();
    dt = cfl_dt;
    /*
    if (cfl_dt < dt) {
        printf("Using cfl limited timestep of %.4e instead of %.4e\n",cfl_dt,dt);
        dt = cfl_dt;
    }
    */

   output_init(outputdir);
   init_rk5();
    int i;
    for(i=0; i <nsteps; i++) {
        printf("Step %d of %d\n",i+1,nsteps);
        set_bc();
        set_avg(0);
    
        potential();
        move_planet();
    
       source_step();
       
        set_Lamex();
       viscosity();
    
        temp_to_vel();
        set_bc();
        
        vel_to_temp();
        
        transport_step();
      
        set_avg(1);

        set_Lamdep();
    }
//    temp_to_vel();   

    time_avg();

    output(outputdir);
    output_torque(outputdir);
    free_rk5();
    //free_all();
    return 0;
}
