/*
 * Николай Игоревич Хохлов, k_h@inbox.ru, 2011-2014.
 * Реализация алгоритма сортировки слиянием.
 * gcc -O2 merge_sort.c
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

void merge(int *a, int *b, int *c, int na, int nb);
void merge_sort(int *a, int na);
void print(int *a, int na);
int check_sort(int *a, int n);
double timer();

int main (int argc, char *argv[])
{
	int i, n, *a;
	double t;
	if (argc < 2) {
		printf("Usage: %s num_elements.\n", argv[0]);
		return 1;
	}
	n = atoi(argv[1]);
	a = (int*)malloc(sizeof(int) * n);
	srand(time(NULL));

	for (i = 0; i < n; i++) a[i] = rand() % 1000;

	if (n < 101) print(a, n);

	t = timer();
	merge_sort(a, n);
	t = timer() - t;

	if (n < 101) print(a, n);
	printf("Time: %f sec, sorted: %d\n", t, check_sort(a, n));
	free(a);
	return 0;
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
	/* Реализовать процедуру слияния. */
}

/*
 * Процедура сортировки слиянием.
 */
void merge_sort(int *a, int na)
{
	/* Реализовать процедуру сортировки, используя процедуру слияния. */
}
