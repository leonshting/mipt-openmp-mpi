/*
 * Author: Nikolay Khokhlov <k_h@inbox.ru>, 2016
 */

#include <stdio.h>
#include <assert.h>
#include <stdlib.h>
#include <stddef.h>
#include "mpi.h"


#define BSM 5
#define USM 6
#define LSM 7
#define RSM 8
#define URCM 9
#define ULCM 10
#define BRCM 11
#define BLCM 12




typedef struct {
	int nx, ny;
    int beg_x, beg_y;
	int *u0;
	int *u1;
	int steps;
	int save_steps;
} life_t;

typedef struct {
    int nx, ny;
    int steps;
    int save_steps;
    int beg_x, beg_y;
} life_t_reduced;

int ind(int i, int j, life_t *l)
{
    return (((i + l->nx) % l->nx) + ((j + l->ny) % l->ny) * (l->nx));
}

int ind2(int i, int j, life_t_reduced *l)
{
    return (((i + l->nx) % l->nx) + ((j + l->ny) % l->ny) * (l->nx));
}

void life_init(const char *path, life_t *l);
void life_free(life_t *l);
void life_step(life_t *l);
void life_save_vtk(const char *path, life_t *l);
void life_save_vtk2(const char *path, int * storage, life_t_reduced * lives, int numlives, int * blocklens, int * blockstrides, life_t * l)
{
    FILE *f;
    int i, i1, i2, j;
    f = fopen(path, "w");
    assert(f);
    fprintf(f, "# vtk DataFile Version 3.0\n");
    fprintf(f, "Created by write_to_vtk2d\n");
    fprintf(f, "ASSII\n");
    fprintf(f, "DATASET STRUCTURED_POINTS\n");
    fprintf(f, "DIMENSIONS %d %d 1\n", l->nx+1, l->ny+1);
    fprintf(f, "SPACING %d %d 0.0\n", 1, 1);
    fprintf(f, "ORIGIN %d %d 0.0\n", 0, 0);
    fprintf(f, "CELL_DATA %d\n", l->nx * l->ny);

    fprintf(f, "SCALARS life int 1\n");
    fprintf(f, "LOOKUP_TABLE life_table\n");
    for(i = 0; i < numlives ; i++) {
        for (j = 0; j < (lives[i].nx-2) * (lives[i].ny-2); j++) {
            l->u0[ind(lives[i].beg_x+1 + j % (lives[i].nx-2), lives[i].beg_y+1 + j / (lives[i].nx-2), l)] =
                    storage[blockstrides[i] + ind2(lives[i].beg_x + 1 + j % (lives[i].nx-2), lives[i].beg_y + 1 + j / (lives[i].nx-2),&lives[i])];
        }
    }
    for (i2 = 0; i2 < l->ny; i2++) {
        for (i1 = 0; i1 < l->nx; i1++) {
            /*if(l->u0[ind(i1, i2, l)] == 1)*/ fprintf(f, "%d\n",l->u0[ind(i1, i2, l)]);
        }
    }
    fclose(f);
}
void make_part(life_t *old, life_t *new, int coords[], int proc_sizes[]);
void part_of_life(life_t *old, life_t *new, int beg[]);
void make_reduced_copy_of_life(life_t *some, life_t_reduced * some2);
void init_life_from_reduced(life_t *some, life_t_reduced * some2);
int * ls(life_t * some)
{
    int * tmp = (int *) calloc(some->ny-2, sizeof(int));
    int x = 1;
    int i = 1;
    for(i = 1; i < some->ny-1; i++)
        tmp[i-1] = some->u0[ind(x,i, some)];
    return tmp;
}
int * rs(life_t * some)
{
    int * tmp = (int *) calloc(some->ny-2, sizeof(int));
    int x = some->nx-2;
    int i = 1;
    for(i = 1; i < some->ny-1; i++)
        tmp[i-1] = some->u0[ind(x,i, some)];
    return tmp;
}
int * us(life_t * some)
{
    int * tmp = (int *) calloc(some->nx-2, sizeof(int));
    int y = 1;
    int i = 1;
    for(i = 1; i < some->nx-1; i++)
        tmp[i-1] = some->u0[ind(i,y, some)];
    return tmp;
}
int * bs(life_t * some)
{
    int * tmp = (int *) calloc(some->nx-2, sizeof(int));
    int y = some->ny-2;
    int i = 1;
    for(i = 1; i < some->nx-1; i++)
        tmp[i-1] = some->u0[ind(i,y,some)];
    return tmp;
}
int ul(life_t * some)
{
    return some->u0[ind(1,1,some)];
}
int ur(life_t * some)
{
    return some->u0[ind(some->nx-1,1, some)];
}
int bl(life_t * some)
{
    return some->u0[ind(1,some->ny-1, some)];
}
int br(life_t * some)
{
    return some->u0[ind(some->nx-1,some->ny-1, some)];
}

void add_to_cell(int ulc, int urc, int blc, int brc, int *rs, int *bs, int *us, int* ls, life_t * some)
{
    int i;
    for(i = 1; i < some->nx-1; i++)
    {
        some->u0[ind(i, 0, some)] = us[i-1];
        some->u0[ind(i, some->ny-1, some)] = bs[i-1];
    }
    for(i = 1; i < some->ny-1; i++)
    {
        some->u0[ind(0, i, some)] = us[i-1];
        some->u0[ind(some->nx-1, i, some)] = bs[i-1];
    }
    some->u0[ind(0,0, some)] = ulc;
    some->u0[ind(some->nx-1,0, some)] = urc;
    some->u0[ind(some->nx-1,some->ny-1, some)] = brc;
    some->u0[ind(0,some->ny-1, some)] = blc;
}


int main(int argc, char **argv)
{
    int myid, numprocs, dim_size[2] = {0};
    int period_tor[2] = {1, 1};
    int period_cylinder[2] = {1,0};
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);
    MPI_Comm new_comm;
    MPI_Request *req;
    req = malloc(sizeof(MPI_Request)*numprocs);


    MPI_Datatype lifetype;

	if (myid == 0 && argc != 2) {
		printf("Usage: %s input file.\n", argv[0]);
		return 0;
	}

    int blocklens[] = {1, 1, 1, 1, 1, 1};
    MPI_Aint indices[] = {(MPI_Aint)offsetof(life_t_reduced, nx), (MPI_Aint)offsetof(life_t_reduced, ny), (MPI_Aint)offsetof(life_t_reduced, steps),
                          (MPI_Aint)offsetof(life_t_reduced, save_steps), (MPI_Aint)offsetof(life_t_reduced, beg_x),(MPI_Aint)offsetof(life_t_reduced, beg_y)};
    MPI_Datatype old_types[] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_INT};
    MPI_Type_create_struct(6, blocklens, indices, old_types, &lifetype);
    MPI_Type_commit(&lifetype);

	life_t l;
    life_t_reduced reduced;
    life_t some;
    life_t_reduced dimreduced[numprocs];
    int *blocklens_forroot;
    int *strides_forroot;
    int sumsize = 0;
    if(myid == 0) {
        life_init(argv[1], &l);
    }

    MPI_Dims_create(numprocs, 2, dim_size);
    MPI_Cart_create(MPI_COMM_WORLD, 2, dim_size, period_tor, 0, &new_comm);
    if(myid == 0)
    {
        int i;
        blocklens_forroot = calloc(numprocs, sizeof(int));
        strides_forroot = calloc(numprocs, sizeof(int));
        for(i = 0; i < numprocs; i++)
        {
            life_t *tmp; life_t_reduced to_send;
            tmp = malloc(sizeof(life_t));
            int coords[2];
            MPI_Cart_coords(new_comm, i, 2, coords);
            make_part(&l, tmp, coords, dim_size);
            make_reduced_copy_of_life(tmp, &to_send);
            strides_forroot[i] = ((i-1)<0)?0:(strides_forroot[i-1] + tmp->nx * tmp->ny);
            blocklens_forroot[i] = tmp->nx * tmp->ny;
            sumsize += blocklens_forroot[i];
            printf("%d %d\n", strides_forroot[i], blocklens_forroot[i]);
            dimreduced[i] = to_send;
            if(i!=0)
            {
                MPI_Send(&to_send, 1, lifetype, i, 0, new_comm);
                MPI_Send(tmp->u0, tmp->nx * tmp->ny, MPI_INT, i, 1, new_comm);
            }
            else
            {
                some = *tmp;
            }
        }

    }
    else
    {
        MPI_Recv(&reduced, 1, lifetype, 0, 0, new_comm, MPI_STATUS_IGNORE);

        init_life_from_reduced(&some, &reduced);
        MPI_Recv(some.u0, some.ny * some.nx, MPI_INT, 0, 1, new_comm, MPI_STATUS_IGNORE);
    }
    MPI_Barrier(new_comm);
	int i, mc[2];
    int *bsp, *usp, *lsp, *rsp;
    int ulc, urc, brc, blc;
    int bsc[2], usc[2], lsc[2], rsc[2];
    int ulcc[2], urcc[2], brcc[2], blcc[2];
    int bsr, usr, lsr, rsr;
    int ulcr, urcr, brcr, blcr;


    int *bsp_r, *usp_r, *lsp_r, *rsp_r;
    int ulc_r, blc_r, urc_r, brc_r;
    bsp_r = (int *) calloc(some.nx-2, sizeof(int));
    usp_r = (int *) calloc(some.nx-2, sizeof(int));
    lsp_r = (int *) calloc(some.ny-2, sizeof(int));
    rsp_r = (int *) calloc(some.ny-2, sizeof(int));

    MPI_Cart_coords(new_comm, myid, 2, mc);

    bsc[0]  = mc[0];   bsc[1]  = mc[1]+1;
    usc[0]  = mc[0];   usc[1]  = mc[1]-1;
    lsc[0]  = mc[0]-1; lsc[1]  = mc[1];
    rsc[0]  = mc[0]+1; rsc[1]  = mc[1];
    ulcc[0] = mc[0]-1; ulcc[1] = mc[1]-1;
    urcc[0] = mc[0]+1; urcc[1] = mc[1]-1;
    blcc[0] = mc[0]-1; blcc[1] = mc[1]+1;
    brcc[0] = mc[0]+1; brcc[1] = mc[1]+1;

    MPI_Cart_rank(new_comm, bsc, &bsr);
    MPI_Cart_rank(new_comm, usc, &usr);
    MPI_Cart_rank(new_comm, lsc, &lsr);
    MPI_Cart_rank(new_comm, rsc, &rsr);
    MPI_Cart_rank(new_comm, ulcc, &ulcr);
    MPI_Cart_rank(new_comm, urcc, &urcr);
    MPI_Cart_rank(new_comm, brcc, &brcr);
    MPI_Cart_rank(new_comm, blcc, &blcr);

    MPI_Request send[8];
    MPI_Request sendi[8];
    int *storage = calloc(sumsize, sizeof(int));

	char buf[100];
	for (i = 0; i < some.steps; i++) {
        bsp = bs(&some); usp = us(&some); lsp = ls(&some); rsp = rs(&some);
        ulc = ul(&some); urc = ur(&some); brc = br(&some); blc = bl(&some);
        if (i % some.save_steps == 0) {

            MPI_Gatherv(some.u0, some.nx*some.ny, MPI_INT, storage, blocklens_forroot, strides_forroot, MPI_INT, 0, new_comm);
            if(myid == 0)
            {   
                sprintf(buf, "life_%06d.vtk", i);
               printf("Saving step %d to '%s'.\n", i, buf);
               life_save_vtk2(buf, storage, dimreduced, numprocs, blocklens_forroot, strides_forroot, &l);
            }
        }

        MPI_Isend(bsp, some.nx-2, MPI_INT, bsr, BSM, new_comm, &send[BSM-5]);
        MPI_Isend(lsp, some.ny-2, MPI_INT, lsr, LSM, new_comm, &send[LSM-5]);
        MPI_Isend(usp, some.nx-2, MPI_INT, usr, USM, new_comm, &send[USM-5]);
        MPI_Isend(rsp, some.ny-2, MPI_INT, rsr, RSM, new_comm, &send[RSM-5]);

        MPI_Isend(&ulc, 1, MPI_INT, ulcr, ULCM, new_comm, &send[ULCM-5]);
        MPI_Isend(&urc, 1, MPI_INT, urcr, URCM, new_comm, &send[URCM-5]);
        MPI_Isend(&brc, 1, MPI_INT, brcr, BRCM, new_comm, &send[BRCM-5]);
        MPI_Isend(&blc, 1, MPI_INT, blcr, BLCM, new_comm, &send[BLCM-5]);

        MPI_Irecv(bsp_r, some.nx-2, MPI_INT, bsr, USM, new_comm, &sendi[USM-5]);
        MPI_Irecv(lsp_r, some.ny-2, MPI_INT, lsr, RSM, new_comm, &sendi[RSM-5]);
        MPI_Irecv(usp_r, some.nx-2, MPI_INT, usr, BSM, new_comm, &sendi[BSM-5]);
        MPI_Irecv(rsp_r, some.ny-2, MPI_INT, rsr, LSM, new_comm, &sendi[LSM-5]);

        MPI_Irecv(&ulc_r, 1, MPI_INT, brcr, ULCM, new_comm, &sendi[ULCM-5]);
        MPI_Irecv(&urc_r, 1, MPI_INT, blcr, URCM, new_comm, &sendi[URCM-5]);
        MPI_Irecv(&brc_r, 1, MPI_INT, ulcr, BRCM, new_comm, &sendi[BRCM-5]);
        MPI_Irecv(&blc_r, 1, MPI_INT, urcr, BLCM, new_comm, &sendi[BLCM-5]);

        MPI_Waitall(8, sendi, MPI_STATUSES_IGNORE);
        add_to_cell(ulc_r, urc_r, blc_r, brc_r, rsp_r, bsp_r, usp_r, lsp_r, &some);
		life_step(&some);
	}
    free(bsp); free(bsp_r);
    free(lsp); free(lsp_r);
    free(usp); free(usp_r);
    free(rsp); free(rsp_r);
    free(storage);
    MPI_Finalize();
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
void init_life_from_reduced(life_t *some, life_t_reduced * some2)
{
    some->beg_x = some2->beg_x;
    some->beg_y = some2->beg_y;
    some->nx = some2->nx;
    some->steps = some2->steps;
    some->save_steps = some2->save_steps;
    some->ny = some2 -> ny;
    some->u0 = (int*)calloc(some->nx * some->ny, sizeof(int));
    some->u1 = (int*)calloc(some->nx * some->ny, sizeof(int));
}

void make_part(life_t *old, life_t *new, int coords[], int proc_sizes[])
{
    int beg[2], end[2];
    beg[0] = (int)(old->nx * coords[0])/proc_sizes[0]-1;
    beg[1] = (int)(old->ny * coords[1])/proc_sizes[1]-1;
    end[0] = (int)(old->nx * (coords[0]+1))/proc_sizes[0]+1;
    end[1] = (int)(old->ny * (coords[1]+1))/proc_sizes[1]+1;
    new -> nx = end[0] - beg[0];
    new -> ny = end[1] - beg[1];
    printf("%d %d\n", new->nx, new->ny);
    new->beg_x = beg[0];
    new->beg_y = beg[1];
    new -> save_steps = old -> save_steps;
    new -> steps = old -> steps;
    int i,j;
    new->u0 = (int*)calloc(new->nx * new->ny, sizeof(int));
    new->u1 = (int*)calloc(new->nx * new->ny, sizeof(int));
    for(i = 1; i < new->nx-1; i++)
    {
        for(j = 1; j < new->ny-1; j++)
        {
            new->u0[ind(i,j, new)] = old->u0[ind(beg[0] + i, beg[1] +j, old)];
        }
    }
}
void make_reduced_copy_of_life(life_t *some, life_t_reduced * some2)
{
    some2->nx = some->nx;
    some2->steps = some->steps;
    some2->save_steps = some->save_steps;
    some2->ny = some -> ny;
    some2->beg_x = some->beg_x;
    some2->beg_y = some->beg_y;
}
void life_init(const char *path, life_t *l)
{
	FILE *fd = fopen(path, "r");
	assert(fd);
	assert(fscanf(fd, "%d\n", &l->steps));
	assert(fscanf(fd, "%d\n", &l->save_steps));
	printf("Steps %d, save every %d step.\n", l->steps, l->save_steps);
	assert(fscanf(fd, "%d %d\n", &l->nx, &l->ny));
	printf("Field size: %dx%d\n", l->nx, l->ny);
    l->beg_x = 0;
    l->beg_y = 0;
	l->u0 = (int*)calloc(l->nx * l->ny, sizeof(int));
	l->u1 = (int*)calloc(l->nx * l->ny, sizeof(int));
	
	int i, j, r, cnt;
	cnt = 0;
	while ((r = fscanf(fd, "%d %d\n", &i, &j)) != EOF) {
		l->u0[ind(i, j, l)] = 1;
		cnt++;
	}
	printf("Loaded %d life cells.\n", cnt);
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
			fprintf(f, "%d\t%d\t%d\n", i1+l->beg_x, i2 + l->beg_y, l->u0[ind(i1, i2, l)]);
		}
	}
	fclose(f);
}

void life_step(life_t *l)
{
	int i, j;
	for (j = 0; j < l->ny; j++) {
		for (i = 0; i < l->nx; i++) {
			int n = 0;
			n += l->u0[ind(i+1, j, l)];
			n += l->u0[ind(i+1, j+1, l)];
			n += l->u0[ind(i,   j+1, l)];
			n += l->u0[ind(i-1, j, l)];
			n += l->u0[ind(i-1, j-1, l)];
			n += l->u0[ind(i,   j-1, l)];
			n += l->u0[ind(i-1, j+1, l)];
			n += l->u0[ind(i+1, j-1, l)];
			l->u1[ind(i,j, l)] = 0;
			if (n == 3 && l->u0[ind(i,j, l)] == 0) {
				l->u1[ind(i,j, l)] = 1;
			}
			if ((n == 3 || n == 2) && l->u0[ind(i,j, l)] == 1) {
				l->u1[ind(i,j, l)] = 1;
			}
		}
	}
	int *tmp;
	tmp = l->u0;
	l->u0 = l->u1;
	l->u1 = tmp;
}



