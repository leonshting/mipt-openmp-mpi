/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, 2016
 */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "omp.h"


#define ind(i, j) (((i + l->nx) % l->nx) + ((j + l->ny) % l->ny) * (l->nx))

typedef struct {
	int nx, ny;
	int *u0;
	int *u1;
	int steps;
	int save_steps;
} life_t;

int *start_x;
int *start_y;
int *len_x;
int *len_y;
int OMP_NUM_THREADS;

void life_init(const char *path, life_t *l);
void life_free(life_t *l);
void life_step(life_t *l);
void life_save_vtk(const char *path, life_t *l);
int primeFactors(int n, int * ret)
{
    int count = 0;
    while (n%2 == 0)
    {
        n = n/2;
        ret[count++] = 2;
    }
    int i;
    for (i = 3; i <= sqrt(n); i = i+2)
    {
        while (n%i == 0) {
            n = n / i;
            ret[count++] = i;
        }
    }
    return count;
}

int main(int argc, char **argv)
{
	if (argc != 2) {
		printf("Usage: %s input file.\n", argv[0]);
		return 0;
	}
    omp_set_num_threads(16);
	life_t l;
	life_init(argv[1], &l);
    #pragma omp parallel
    {
        OMP_NUM_THREADS = omp_get_num_threads();
    }
    start_x = malloc(sizeof(int) * OMP_NUM_THREADS);
    start_y = malloc(sizeof(int) * OMP_NUM_THREADS);
    len_x = malloc(sizeof(int) * OMP_NUM_THREADS);
    len_y = malloc(sizeof(int) * OMP_NUM_THREADS);

    int *factorized = malloc(sizeof(int) * OMP_NUM_THREADS);
    int fac_num = primeFactors(OMP_NUM_THREADS, factorized);
    int div_x=1, div_y=1, i;
    for(i = 0; i<fac_num; i++)
    {
        if(i%2 == 0)
            div_x *= factorized[i];
        else
            div_y *= factorized[i];
     //   printf("%d\t", factorized[i]);
    }
    //printf("%d\n", div_x);
    //printf("%d\n", OMP_NUM_THREADS);
    int x = 0;
    int y = 0;
    int rank = 0;
    int tmpx = 0;
    for(x = 0; x < div_x; x++)
    {
        int tmpy = 0;
        for(y = 0; y < div_y; y++)
        {
            len_x[rank] = l.nx/div_x + ((x < (l.nx % div_x))?1:0);
            len_y[rank] = l.ny/div_y + ((y < (l.ny % div_y))?1:0);
            start_y[rank] = tmpy;
            start_x[rank] = tmpx;
            tmpy+=len_y[rank];
           // printf("%d %d %d\n", div_x,  start_x[rank], start_y[rank]);
            rank++;

        }
        tmpx += l.nx/div_x + ((x<(l.nx % div_x))?1:0);
    }


	char buf[100];
    #pragma omp parallel private (i)
    {
        int rank = omp_get_thread_num();
        for (i = 0; i < l.steps; i++) {
		/*if (i % l.save_steps == 0) {
			sprintf(buf, "life_%06d.vtk", i);
			printf("Saving step %d to '%s'.\n", i, buf);
			life_save_vtk(buf, &l);
		}*/
		    life_step(&l);
	    }
        #pragma omp barrier
    }
	life_free(&l);
	return 0;
}

/**
 * Загрузить входную конфигурацию.
 * Формат файла, число шагов, как часто сохранять, размер поля, затем идут координаты заполненых клеток:
 * steps
 * save_steps
 * nx ny
 * i1 j2
 * i2 j2
 */
void life_init(const char *path, life_t *l)
{
	FILE *fd = fopen(path, "r");
	assert(fd);
	assert(fscanf(fd, "%d\n", &l->steps));
	assert(fscanf(fd, "%d\n", &l->save_steps));
//	printf("Steps %d, save every %d step.\n", l->steps, l->save_steps);
	assert(fscanf(fd, "%d %d\n", &l->nx, &l->ny));
//	printf("Field size: %dx%d\n", l->nx, l->ny);

	l->u0 = (int*)calloc(l->nx * l->ny, sizeof(int));
	l->u1 = (int*)calloc(l->nx * l->ny, sizeof(int));
	
	int i, j, r, cnt;
	cnt = 0;
	while ((r = fscanf(fd, "%d %d\n", &i, &j)) != EOF) {
		l->u0[ind(i, j)] = 1;
		cnt++;
	}
//	printf("Loaded %d life cells.\n", cnt);
	fclose(fd);
}

void life_free(life_t *l)
{
	free(l->u0);
	free(l->u1);
	l->nx = l->ny = 0;
}

void life_save_vtk(const char *path, life_t *l)
{
	FILE *f;
	int i1, i2, j;
	f = fopen(path, "w");
	assert(f);
	fprintf(f, "# vtk DataFile Version 3.0\n");
	fprintf(f, "Created by write_to_vtk2d\n");
	fprintf(f, "ASCII\n");
	fprintf(f, "DATASET STRUCTURED_POINTS\n");
	fprintf(f, "DIMENSIONS %d %d 1\n", l->nx+1, l->ny+1);
	fprintf(f, "SPACING %d %d 0.0\n", 1, 1);
	fprintf(f, "ORIGIN %d %d 0.0\n", 0, 0);
	fprintf(f, "CELL_DATA %d\n", l->nx * l->ny);
	
	fprintf(f, "SCALARS life int 1\n");
	fprintf(f, "LOOKUP_TABLE life_table\n");
	for (i2 = 0; i2 < l->ny; i2++) {
		for (i1 = 0; i1 < l->nx; i1++) {
			fprintf(f, "%d\n", l->u0[ind(i1, i2)]);
		}
	}
	fclose(f);
}

void life_step(life_t *l)
{
        int i, j;
        int rank = omp_get_thread_num();
    //    printf("%d\t", rank);
        fflush(stdout);
        for (j = start_y[rank]; j < start_y[rank] + len_y[rank]; j++) {
            for (i = start_x[rank]; i < start_x[rank] + len_x[rank]; i++) {
                int n = 0;
                n += l->u0[ind(i + 1, j)];
                n += l->u0[ind(i + 1, j + 1)];
                n += l->u0[ind(i, j + 1)];
                n += l->u0[ind(i - 1, j)];
                n += l->u0[ind(i - 1, j - 1)];
                n += l->u0[ind(i, j - 1)];
                n += l->u0[ind(i - 1, j + 1)];
                n += l->u0[ind(i + 1, j - 1)];
                l->u1[ind(i, j)] = 0;
                if (n == 3 && l->u0[ind(i, j)] == 0) {
                    l->u1[ind(i, j)] = 1;
                }
                if ((n == 3 || n == 2) && l->u0[ind(i, j)] == 1) {
                    l->u1[ind(i, j)] = 1;
                }
            }
        }
    #pragma omp master
        {int *tmp;
	tmp = l->u0;
	l->u0 = l->u1;
	l->u1 = tmp;}
}


