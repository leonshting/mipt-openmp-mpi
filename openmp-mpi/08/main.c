#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

double SENDRECV(int myid, int size, int count)
{
    int i;
    char *buffer;
    buffer = malloc(size * sizeof(char));
    if(myid == 0)
        memset(buffer, '-', size);
    MPI_Barrier(MPI_COMM_WORLD);
    int start = clock();
    for(i=0; i < count; i++)
    {
        if(myid == 0)
        {
            MPI_Send(buffer, size, MPI_CHAR, 1, size, MPI_COMM_WORLD);
        }
        if(myid == 1)
        {
            MPI_Recv(buffer, size, MPI_CHAR, 0, size, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    int end = clock();
    return (double) (end - start)/(1.0 * count);
}


double ISENDRECV(int myid, int size, int count)
{
    int i;
    char **buffer;
    MPI_Request* requests;
    buffer = malloc(count * sizeof(char *));
    requests = malloc(count * sizeof(MPI_Request));
    for (i = 0; i < count; i++) {
        buffer[i] = malloc(size * sizeof(char));
        memset(buffer[i], (unsigned char) i % 255, size);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    int start = clock();
    for(i=0; i < count; i++)
    {
        if(myid == 0)
        {
            MPI_Isend(buffer[i], size, MPI_CHAR, 1, i, MPI_COMM_WORLD, &requests[i]);
        }
        if(myid == 1)
        {
            MPI_Irecv(buffer[i], size, MPI_CHAR, 0, i, MPI_COMM_WORLD, &requests[i]);
        }
    }
    if(myid == 0)
        MPI_Waitall(count, &requests[0], MPI_STATUSES_IGNORE);

    int end = clock();
    free(buffer);
    return (double) (end - start)/count;
}

double BSENDRECV(int myid, int size, int count)
{
    int i;
    char *buffer, *buf;
    buffer = malloc(size * sizeof(char));
    int tmp_size;
    MPI_Pack_size(size, MPI_CHAR, MPI_COMM_WORLD, &tmp_size);
    buf = malloc(count * tmp_size + count * MPI_BSEND_OVERHEAD);
    if(myid == 0) {
        memset(buffer, '-', size);
    }
    MPI_Buffer_attach(buf, count * tmp_size + count * MPI_BSEND_OVERHEAD);
    MPI_Barrier(MPI_COMM_WORLD);
    int start = clock();
    for(i=0; i < count; i++)
    {
        if(myid == 0)
        {
            MPI_Bsend(buffer, size, MPI_CHAR, 1, size, MPI_COMM_WORLD);
        }
        if(myid == 1)
        {
            MPI_Recv(buffer, size, MPI_CHAR, 0, size, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    int end = clock();
    int some;
    MPI_Buffer_detach(buf, &some);
    free(buf);
    free(buffer);
    return (double) (end - start)/(1.0 * count);
}

double SSENDRECV(int myid, int size, int count)
{
    int i;
    char *buffer;
    buffer = malloc(size * sizeof(char));
    if(myid == 0)
        memset(buffer, '-', size);
    MPI_Barrier(MPI_COMM_WORLD);
    int start = clock();
    for(i=0; i < count; i++)
    {
        if(myid == 0)
        {
            MPI_Ssend(buffer, size, MPI_CHAR, 1, size, MPI_COMM_WORLD);
        }
        if(myid == 1)
        {
            MPI_Recv(buffer, size, MPI_CHAR, 0, size, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    int end = clock();
    free(buffer);
    return (double) (end - start)/(1.0 * count);
}


int main(int argc, char *argv[]) {
    int numprocs, size, myid, mode;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    int count = 1000;
    for(size = 1; size<20000000; size*=2)
    {
        double avg_time[4] = {0.0};

        avg_time[0] = SENDRECV(myid, size, count);
        //avg_time[1] = ISENDRECV(myid, size, count);
        avg_time[2] = SSENDRECV(myid, size, count);
        avg_time[3] = BSENDRECV(myid, size, count);
        if(myid == 0)
        {
            printf("%d\t%lf\t%lf\t%lf\n",size, size/avg_time[0], size/avg_time[2], size/avg_time[3]);
        }
    }
    MPI_Finalize();
    return 0;
}
