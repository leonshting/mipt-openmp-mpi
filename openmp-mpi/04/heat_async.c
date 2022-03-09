#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define N_x 100
#define N_y 100
#define t 0.00001
#define h 0.01
#define T 0.0001

double func(double x, double y)
{
    if(x >= 0.4 && x <= 0.6 && y >= 0.4 && y <= 0.6) return 1.0;
    else return 0.0;
}

int dir_map_num(int i, int j)
{
    return i + j*N_x;
}

void iter_area(double *inp, double *outp, int size_y, int start_y)
{
    int i = 0, j, e_size = size_y-1;
    for(j = 1; j<e_size; j++)
    {

        i = 0;
        while(i<N_x)
        {
           // printf("start:%d %d %d %d\n", inp, i,j, e_size); fflush(stdout);
            if(i != 0 && i != (N_x-1))
            {
                outp[dir_map_num(i,j)] = inp[dir_map_num(i,j)]+
                                         t/pow(h,2.0)*(inp[dir_map_num(i-1,j)] + inp[dir_map_num(i,j-1)]
                                         -4*inp[dir_map_num(i,j)]+ inp[dir_map_num(i+1,j)] + inp[dir_map_num(i,j+1)]);
            }
            else
            {
                outp[dir_map_num(i,j)] = inp[dir_map_num(i,j)];
            }
            i++;
          //  fflush(stderr);
        }
    }
    if((size_y + start_y) == N_y  || start_y == 0)
    {
        i = 0;
        j=(start_y!=0)?(size_y-1):0;
        while(i<N_x)
        {
            outp[dir_map_num(i,j)] = inp[dir_map_num(i,j)];
            i++;
        }
    }

}

int main(int argc, char * argv[])
{
    int N_tot = N_x*N_y, i, j, active = 0, inactive = 1;
    int size, rank, size_y, start_y, end_y, r;
    char str[100];
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    sprintf(str,"file%d.txt", rank);
    FILE * f = fopen(str,"w");
    double *dir_map[2], t_tot=0.0;
    r = atoi(argv[1]);
    size_y = (rank == size-1)?(N_y/size+N_y%size):N_y/size;
    start_y = (rank == 0)?0:(N_y/size)*rank-r;
    size_y = (rank == size-1 || rank == 0)? size_y+r: size_y+2*r;
    end_y = (rank == size-1)?(N_y):(start_y + size_y);
   // printf("MYID:%d ,size:%d, start:%d, end:%d\n", rank, size_y, start_y, end_y);
    for(i = 0; i < 2; i++)
        dir_map[i] = malloc(sizeof(double)*N_x*size_y);
    for(i = start_y; i < end_y; i++)
    {
        for(j = 0; j < N_x; j++)
            dir_map[active][dir_map_num(j, i-start_y)] = func((double)j/(double)N_x, (double)i/(double)N_y);
    }
    MPI_Request req_l, req_r;
    MPI_Request req_l1, req_r1;
    double * tmp[2];
    double * tmp2[2];
    int incr=0;
    for(incr = 0; incr<2; incr++)
    {
        tmp[incr] = malloc(sizeof(double) * (3*r)*N_x);
        tmp2[incr] = malloc(sizeof(double) * (3*r)*N_x);
    }

    while(t_tot < T)
    {
        int old_active = active;
        int old_inactive = inactive;
        if(rank == 0)
        {
            memcpy(tmp[active], &(dir_map[old_active][size_y*N_x-3*r*N_x]), N_x*3*r*sizeof(double));
            memcpy(tmp[inactive], &(dir_map[old_inactive][size_y*N_x-3*r*N_x]), N_x*3*r*sizeof(double));
        }
        else if(rank == size-1)
        {
            memcpy(tmp2[active], &(dir_map[old_active][0]), N_x*3*r*sizeof(double));
            memcpy(tmp2[inactive], &(dir_map[old_inactive][0]), N_x*3*r*sizeof(double));
        }
        else
        {
            memcpy(tmp[active], &(dir_map[old_active][size_y*N_x-3*r*N_x]), N_x*3*r*sizeof(double));
            memcpy(tmp2[active], &(dir_map[old_active][0]), N_x*3*r*sizeof(double));
            memcpy(tmp[inactive], &(dir_map[old_inactive][size_y*N_x-3*r*N_x]), N_x*3*r*sizeof(double));
            memcpy(tmp2[inactive], &(dir_map[old_inactive][0]), N_x*3*r*sizeof(double));
        }
        for(i = 0; i < r ; i++)
        {
            if(rank==0)
            {
                iter_area(&tmp[active][i*N_x], &tmp[inactive][i*N_x],
                        3*r-2*i,start_y+size_y-3*r+i);
            }
            else if(rank == size-1)
            {
                iter_area(&tmp2[active][i*N_x], &tmp2[inactive][i*N_x],
                        3*r-2*i,start_y+i);
            }
            else
            {
                iter_area(&tmp[active][i*N_x], &tmp[inactive][i*N_x],
                        3*r-2*i,start_y+size_y-3*r+i);
                iter_area(&tmp2[active][i*N_x], &tmp2[inactive][i*N_x],
                        3*r-2*i,start_y+i);
            }
            int sw = active;
            active = inactive;
            inactive = sw;
        }
        if(rank==0)
        {
            MPI_Isend(&tmp[active][r*N_x], r*N_x, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &req_r);
            MPI_Irecv(&dir_map[old_active][(size_y-r)*N_x], r*N_x, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &req_l);
        }
        else if(rank==size-1)
        {
            MPI_Isend(&tmp2[active][r*N_x], r*N_x, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &req_l);
            MPI_Irecv(dir_map[old_active], r*N_x, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &req_r);
        }
        else
        {
            MPI_Isend(&tmp[active][r*N_x], r*N_x, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &req_r);
            MPI_Irecv(&dir_map[old_active][(size_y-r)*N_x], r*N_x, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD,&req_l);
            MPI_Isend(&tmp2[active][r*N_x], r*N_x, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &req_l1);
            MPI_Irecv(dir_map[old_active], r*N_x, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &req_r1);
        }

        //async part end
        int sw1 = active, sw2 = inactive;
        active = old_active;
        inactive = old_inactive;
        for(i = 0; i < r ; i++)
        {
            iter_area(&(dir_map[active][(rank == 0)?0:((i+r)*N_x)]), &(dir_map[inactive][(rank == 0)?0:((i+r)*N_x)]),
                      (rank == size-1 || rank == 0)?(size_y-i-r):(size_y-2*i-2*r),(rank == 0)?0:start_y+i+r);
            int sw = active;
            active = inactive;
            inactive = sw;
            t_tot+=t;
        }

        memcpy(&(dir_map[active][size_y*N_x-2*r*N_x]), &tmp[sw1][r*N_x], N_x*r*sizeof(double));
        memcpy(&(dir_map[active][r*N_x]), &tmp2[sw1][r*N_x], N_x*r*sizeof(double));

        MPI_Wait(&req_l, MPI_STATUS_IGNORE);
        MPI_Wait(&req_r, MPI_STATUS_IGNORE);
        if(rank!=0 && rank!= size-1)
        {
            MPI_Wait(&req_l1, MPI_STATUS_IGNORE);
            MPI_Wait(&req_r1, MPI_STATUS_IGNORE);
        }

        printf("t=%lf\n",t_tot);
    }
    for(i=start_y;i<end_y;i++)
    {
        for(j=0;j<N_x;j++)
        {
            fprintf(f,"%lf\t%lf\t%lf\n",(double)j/(double)N_x, (double)(i)/(double)N_y, dir_map[active][dir_map_num(j,i-start_y)]);
            fflush(f);
        }
    }
    MPI_Finalize();
    return 0;
}
