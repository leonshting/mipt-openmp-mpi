#include <stdio.h>
#include <math.h>

#include "mpi.h"

#define BEGIN 0.0
#define END 1.0
int pow2(int k)
{
    int retval;
    retval = 1 << k;
    return retval;
}
int main(int argc, char **argv)
{
    int mask = 0, done = 0, n, myid, numprocs, i, d;
    double h, sum=0.0;
    int namelen;
    char processor_name[MPI_MAX_PROCESSOR_NAME];
    MPI_Status stat;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Get_processor_name(processor_name,&namelen);
    sscanf(argv[1],"%d", &n);
    h = (END - BEGIN)/n;
    d = int(trunc(log2(double(numprocs))));
    if(myid < pow2(d))
    {
        for(i = myid+1;i<=n;i+=pow2(d))
        {
            double x = BEGIN + i*h;
            sum += h * sqrt(4-x*x);
        }
    }
    for(i=0; i < d; i++)
    {
        int partner;
        if((myid & mask) == 0)
        {
            partner = myid ^ pow2(i);
            if(myid & pow2(i) != 0)
            {
                MPI_Send(&sum, 1, MPI_DOUBLE, partner, 0, MPI_COMM_WORLD);
            }
            else
            {
                double tmp;
                MPI_Recv(&tmp, 1, MPI_DOUBLE, partner , 0, MPI_COMM_WORLD, &stat);
                sum += tmp;
            }
        }
        mask = mask ^ pow2(i);
    }
    if(myid==0)
    {
        printf("Result: %lf\n", sum);
    }
  /*  if(myid != 0)
    {
        MPI_Send(&sum, 1, MPI_DOUBLE, 0, myid, MPI_COMM_WORLD);
    }
    else
    {
        for(i=1; i<numprocs;i++)
        {
            double tmp;
            MPI_Recv(&tmp, 1, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &stat);
            sum +=tmp;
            printf("Proc %d's result:%lf\n", i, tmp);
        }
        printf("Result: %lf\n", sum);
    }*/
    MPI_Finalize();
    return 0;
}
