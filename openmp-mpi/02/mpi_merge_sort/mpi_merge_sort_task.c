/*
 * Николай Хохлов, k_h@inbox.ru, 2012.
 * Штанько Леонид, 2015.
 * Реализация алгоритма сортировки слиянием.
 * Изначально массив распределен между процессами, после работы сбор у 0.
 * Сбор всех данных и обработка у одного процесса.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

int pow2(int i)
{
	return 1 << i;
}

int my_log2(int i)
{
	int ii=0;
	while(i!=1) 
	{
		i = i >> 1;
		ii++;
	}
	if(i != pow2(ii))
	 return ii+1;
	else
	 return ii;
}

void merge(int *a, int *b, int *c, int na, int nb);
void merge_sort(int *a, int na);
void print(int *a, int na);
int check_sort(int *a, int n);
double timer();
int n_count(int n, int rank, int size);

int main(int argc, char *argv[])
{
	int d, i, n, nlocal, *a, *b, *c;
	int size, rank, mask = 0;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double t;
	if (argc < 2) {
		if (rank == 0) {
			printf("Usage: %s num_elements.\n", argv[0]);
		}
		return 1;
	}
	n = atoi(argv[1]);
	nlocal = n_count(n, rank, size);
	a = (int*)malloc(sizeof(int) * nlocal);
	
	srand(time(NULL));
	for (i = 0; i < nlocal; i++) a[i] = rand() % 100;
	
	t = timer();
	merge_sort(a, nlocal);
	d = my_log2(size);
	/* Начало сбора. */
	int myid = rank;
    	for(i=0; i < d; i++)
    	{
           if((myid & mask) == 0)
       	   {
        	int npartner = (myid ^ pow2(i));
            	if((myid & pow2(i)) != 0)
            	{
			MPI_Send(&nlocal, 1, MPI_INT, npartner, 1, MPI_COMM_WORLD);
			MPI_Send(a, nlocal, MPI_INT, npartner, 0, MPI_COMM_WORLD);
		}
            	else if(npartner < size)
            	{
			int ni;
			MPI_Recv(&ni, 1, MPI_INT, npartner, 1 , MPI_COMM_WORLD, MPI_STATUS_IGNORE);
			b = (int*)malloc(sizeof(int) * ni);

			MPI_Recv(b, ni, MPI_INT, npartner, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            		//printf("recv from %d to %d\n", npartner, myid);
			c = (int*)malloc(sizeof(int) * (ni + nlocal));
			merge(a, b, c, nlocal, ni);
			free(a);
			a = c;
			nlocal += ni;
			free(b);
            	}
       	    }
       	    mask = (mask ^ pow2(i));
    	}
	/* Конец сбора. */

	if (rank == 0) {
		t = timer() - t;
		if (n < 100) print(a, n);
		printf("Time: %f sec, sorted: %d\n", t, check_sort(a, n));
	}
	free(a);
	MPI_Finalize();
	return 0;
}

int n_count(int n, int rank, int size)
{
	return (n / size) + (n % size > rank);
}

double timer()
{
	struct timeval ts;
	gettimeofday(&ts, NULL);
	return (double)ts.tv_sec + 1e-6 * (double)ts.tv_usec;
}

int check_sort(int *a, int n)
{
	int i;
	for (i = 0; i < n - 1; i++) {
		if (a[i] > a[i+1]) {
			return 0;
		}
	}
	return 1;
}

void print(int *a, int na)
{
	int i;
	for (i = 0; i < na; i++) printf("%d ", a[i]);
	printf("\n");
}

/*
 * Процедура слияния массивов a и b в массив c.
 */
void merge(int *a, int *b, int *c, int na, int nb)
{
	int i = 0, j = 0;
	while (i < na && j < nb) {
		if (a[i] <= b[j]) {
			c[i + j] = a[i];
			i++;
		} else {
			c[i + j] = b[j];
			j++;
		}
	}
	if (i < na) {
		memcpy(c + i + j, a + i, (na - i) * sizeof(int));
	} else {
		memcpy(c + i + j, b + j, (nb - j) * sizeof(int));
	}
}

/*
 * Процедура сортировки слиянием.
 */
void merge_sort(int *a, int na)
{
	if(na < 2) return;
	if(na == 2) {
		if(a[0] > a[1]) { int t = a[0]; a[0]=a[1]; a[1]=t; }
		return;
	}
	merge_sort(a, na / 2);
	merge_sort(a + na / 2, na - na / 2);

	int *b = (int*)malloc(sizeof(int) * na);
	
	merge(a, a + na / 2, b, na / 2, na - na / 2);
	
	memcpy(a, b, sizeof(int) * na);
	
	free(b);
}
