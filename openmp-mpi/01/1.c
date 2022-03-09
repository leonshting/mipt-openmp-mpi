#include <stdio.h>
#include <math.h>

#include "mpi.h"

#define BEGIN 0.0
#define END 1.0
int main(int argc, char **argv)
{
    int done = 0, n, myid, numprocs, i;
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
    for(i = myid+1;i<=n;i+=numprocs)
    {
        double x = BEGIN + i*h;
        sum += h * sqrt(4-x*x);
    }
    if(myid != 0)
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
    }
    MPI_Finalize();
    return 0;
}
