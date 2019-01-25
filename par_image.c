/**** B145208 ****/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

#include "src/arralloc.h"
#include "src/pgmio.h"
#include "src/boundary.h"
#include "src/average.h"

#define MAXITER 200
#define INTERVAL 200

#define N_DIMS 2
#define DIM_Y 1
#define DIM_X 0


static void get_image(int argc, char **argv);
char *filename;

int main (int argc, char **argv)
{
    float **old, **new, **edge, **masterbuf, **buf;
    int rank, size, left, right, up, down;
    //MPI Initiation
    MPI_Init(NULL, NULL);
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &size);
    MPI_Request request;
    MPI_Status status;        
    
    int TAG = 0;
    int num_iters = MAXITER, print_interval = INTERVAL;
    int i, j, iter, sawtooth_value;
    double THRESHOLD = 0.001;
    int M, N, MP, NP;
    float val;
    
    //time monitoirng variables
    double iteration_start_time, iteration_end_time;
    double parallel_code_start_time, parallel_code_end_time;
    
    get_image(argc, argv);
    pgmsize(filename, &M, &N);

    //stopping iteration criteria variables
    float current_delta, max_delta = 0.0;
    float max_delta_all_procs=0.0;
    int check_interval = num_iters/5; //frequency of running this check
    int stop_loop = 0; //to stop the iterations on fulfilling the criteria

    // 2-D Cartesian topology
    //cyclic on first dimension, and non cyclic on second dimension
    float avg_pixel, avg_pixel_sum = 0.0, avg_pixel_value;
    int new_comm, direction_1d, direction_2d, disp;
    int dims[N_DIMS], period[N_DIMS], coords[N_DIMS];
    int coords_x[size], coords_y[size];
    int reorder;
    int x0,x1,y0,y1;
    dims[0] = 0;
    dims[1] = 0;
    coords[0] = 0;
    coords[1] = 0;
    period[0] = DIM_Y;   // DIM_Y, Cyclic horizontal
    period[1] = DIM_X;   // non cyclic vertical
    reorder = DIM_X;
    direction_1d = 0;    // shift along first index
    direction_2d = 1;	   //shift along second index
    disp = 1;            // Shift by 1

    MPI_Dims_create(size, N_DIMS, dims);
    MPI_Cart_create(comm, N_DIMS, dims, period, reorder, &new_comm);
    MPI_Comm_rank(new_comm, &rank);
    MPI_Cart_shift(new_comm, direction_1d, disp, &left, &right);
    MPI_Cart_shift(new_comm, direction_2d, disp, &up, &down);
    MPI_Cart_coords(new_comm, rank, N_DIMS, coords);

    //Initial Coordinate and Rank values
    printf("\n");
    printf("Dimensions of Coordinates : X = %d, Y = %d\n",dims[0],dims[1]);
    printf("Main Rank = %d -- left=%d, right=%d, up=%d, down=%d\n",
           rank, left, right, up, down);
    printf("Cartesian Grid rank = %d, X = %d, Y = %d\n",
           rank, coords[0], coords[1]);
    printf("*****************************************************\n");

    //setting the array size values based on the filename and number of process
    //using static array allocation, as malloc can not handle situation where MxN not divisible by dimension
    MP = M/dims[0];
    NP = N/dims[1];

    masterbuf = (float **) arralloc(sizeof(float), 2, M , N );
    buf       = (float **) arralloc(sizeof(float), 2, MP, NP);

    new  = (float **) arralloc(sizeof(float), 2, MP+2, NP+2);
    old  = (float **) arralloc(sizeof(float), 2, MP+2, NP+2);
    edge = (float **) arralloc(sizeof(float), 2, MP+2, NP+2);

    //MPI derived data type
    //vector datatype of size MPxNP to enable communication between masterbuf & buf
    //vector datatype of size MPx1 for up down halo swaps. Note stride length includes Halos
    MPI_Datatype vectorMPxNP, vectorMPx1;
    MPI_Type_vector(MP, NP, NP*dims[1], MPI_FLOAT, &vectorMPxNP);
    MPI_Type_commit( &vectorMPxNP );
    MPI_Type_vector(MP, 1, NP+2, MPI_FLOAT, &vectorMPx1);
    MPI_Type_commit( &vectorMPx1);

    printf("*****************************************************\n");
    //reading is done on rank 0 only, else would take a lot of time for IO
    //file reading time is measured
    if (rank == 0) {
        printf("\nImage file is :  <%s>\n", filename);
        pgmread(filename, &masterbuf[0][0], M, N);
    }
    
    /*--------------start of parallel region--------------*/
    parallel_code_start_time = MPI_Wtime();

    //keeping track of all cartesian cordinates by rank 0
    //all process send their coordinates to rank 0
    MPI_Gather(&(coords[0]), 1, MPI_INT, &coords_x[0], 1, MPI_INT, 0, new_comm);
    MPI_Gather(&(coords[1]), 1, MPI_INT, &coords_y[0], 1, MPI_INT, 0, new_comm);
    for (i=0; i<size && rank==0; i++) {
        printf("coords_x = %d, coords_y = %d from rank = %d\n", coords_x[i], coords_y[i], i);
    }

    if (rank == 0) {
        //masterbuf array is copied to buf array to be used by every process
        //starting index of masterbuf is based for topology (2d decomposition coordinates)
        for (i=0; i<MP; i++) {
            for(j=0; j<NP; j++) {
                buf[i][j] = masterbuf[i][j]; //copy the initial portion of masterbuf array to buf array for rank 0
            }
        }
        //for all the other ranks, we do a synchronous send using MPI Vector Datatype (vectorMPxNP)
        for (i=1; i<size; i++) {
            MPI_Ssend(&(masterbuf[coords_x[i]*MP][coords_y[i]*NP]), 1, vectorMPxNP, i, TAG, new_comm);
        }
    }
    else {
        MPI_Recv(&buf[0][0], MP*NP, MPI_FLOAT, 0,TAG, new_comm, &status);
    }


    //setting up initial values in array
    //initializing old array to be 255 (white)
    //combining the three loops into a single loop for efficiency
    for (i=0; i<MP+2; i++){
        for (j=0; j<NP+2; j++){
            old[i][j] = 255.0;
            if ((i>=1 && i<=MP) && (j>=1 && j<=NP)) {
                edge[i][j] = buf[i-1][j-1]; //creating copy or temporary array edge
            }
        }
        /* Set fixed boundary conditions on the top and bottom edges
         * compute sawtooth value.
         * sawtooth conditions
         */
        if (i>=1 && i<=MP) {
            sawtooth_value = i + coords[0]*MP;
            val = boundaryval(sawtooth_value, M);
            if (up == MPI_PROC_NULL) {
                old[i][0]   = 255.0*val;
            }
            if (down == MPI_PROC_NULL) {
                old[i][NP+1] = 255.0*(1.0-val);
            }
        }
    }

    /*--------------start of main working iteration--------------*/
    iteration_start_time = MPI_Wtime();
 
    for (iter = 1; iter <= num_iters && stop_loop == 0; iter++){
        if(iter%print_interval == 0){
            printf("Iteration %d, rank = %d\n", iter, rank);

            //calculate average pixel on each process and then reduce it to rank 0
            avg_pixel = average_pixel(MP+2, NP+2, &old[0][0], DIM_Y);
            printf("Rank = %d, Avg pixel = %f\n", rank, avg_pixel);
            MPI_Reduce(&avg_pixel, &avg_pixel_sum, 1, MPI_FLOAT, MPI_SUM, 0, new_comm); //blocking communication
            printf("Avg pixel value = %f\n", avg_pixel_sum/size);
        }

        //halo Swaps Implementation with combination of ISsend & Recvs
        //implementing periodic boundary conditions on left and right sides
        //contiguous memory set
        MPI_Issend(&old[MP][1], NP, MPI_FLOAT, right, TAG, new_comm, &request);
        MPI_Recv(&old[0][1], NP, MPI_FLOAT, left, TAG, new_comm, &status);
        MPI_Wait(&request, &status);

        MPI_Issend(&old[1][1], NP, MPI_FLOAT, left, TAG, new_comm, &request);
        MPI_Recv(&old[MP+1][1], NP, MPI_FLOAT, right, TAG, new_comm, &status);
        MPI_Wait(&request, &status);

        //implement periodic boundary conditions on up and down sides
        //Up and Down swaps are difficult in C, as it's non-contiguous
        MPI_Issend(&old[1][NP], 1, vectorMPx1, down, TAG, new_comm, &request);
        MPI_Recv(&old[1][0], 1, vectorMPx1, up, TAG, new_comm, &status);
        MPI_Wait(&request, &status);

        MPI_Issend(&old[1][1], 1, vectorMPx1, up, TAG, new_comm, &request);
        MPI_Recv(&old[1][NP+1], 1, vectorMPx1, down, TAG, new_comm, &status);
        MPI_Wait(&request, &status);

        for (i=1; i<MP+1; i++){
            for (j=1; j<NP+1; j++){
                new[i][j] = 0.25*(old[i-1][j] + old[i+1][j] + old[i][j-1] + old[i][j+1] - edge[i][j]);
                current_delta = fabs(new[i][j] - old[i][j]);
                if (current_delta > max_delta) 
                    max_delta = current_delta;
            }
        }
        
        //copying new array to old array
        //resetting the array for next iteration
        for (i=1; i<MP+1; i++){
            for (j=1; j<NP+1; j++){
                old[i][j]=new[i][j];
            }
        }

        //stopping criteria check w.r.t THRESHOLD value
        //if stop_loop is set then all the process will exit the iteration
        if ((iter != 0) && (iter%check_interval) == 0) {
            MPI_Allreduce(&max_delta, &max_delta_all_procs, 1, MPI_FLOAT, MPI_MAX, new_comm);
            if (max_delta_all_procs <= THRESHOLD) {
                printf("Iter = %d, max_delta = %f is lower than THRESHOLD=%f. \n"
                       "Stopping further calculation.\n",
                       iter, max_delta_all_procs, THRESHOLD);
                stop_loop=1; //for ending the iteration when stop_loop is set
            }
            max_delta_all_procs = 0.0;
            max_delta = 0.0;
        }

    }
    
    iteration_end_time = MPI_Wtime();
    /*--------------end of main working iteration--------------*/


    printf("Finished %d iterations\n", iter-1);


    for (i=1; i<MP+1; i++){
        for (j=1; j<NP+1; j++){
            buf[i-1][j-1]=old[i][j];
        }
    }

    //all the data is now gathered from different processes
    //all ranks except for 0 sends their buf array to rank 0 (masterbuf) using Ssend
    //masterbuf array recives the data using Recv
    if (rank != 0) {
        MPI_Ssend(&(buf[0][0]), MP*NP, MPI_FLOAT, 0, TAG, new_comm);
    }
    else {
        for (i=0; i<MP; i++) {
            for(j=0; j<NP; j++) {
                masterbuf[i][j] = buf[i][j];
            }
        }
        for (i=1; i<size; i++) {
            MPI_Recv(&(masterbuf[(coords_x[i])*MP][(coords_y[i])*NP]), 1, vectorMPxNP, i, TAG, new_comm, &status);
        }
    }

    parallel_code_end_time = MPI_Wtime();
    /*--------------end of parallel region--------------*/
    

    //masterbuf writing into output file
    //file writing time is measured
    if (rank == 0) {
        sprintf(filename,"output/image%dx%d.pgm",M, N);
        printf("Writting masterbuf to output file\n");
        pgmwrite(filename, &masterbuf[0][0], M, N);
    }


    //time calculation and printing
    printf("Iterative func time (seconds) : %f\n", iteration_end_time-iteration_start_time);
    printf("Parallel code Run time(seconds) : %f\n", parallel_code_end_time - parallel_code_start_time);


    MPI_Finalize();

    return 0;
}


void get_image(int argc, char **argv)
{
    if (argc < 2)
    {
        fprintf(stderr, "Input & Output file name required.\n");
        exit(1);
    }
    filename = argv[1];
}
