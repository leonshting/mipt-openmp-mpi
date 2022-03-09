#!./roundup

gcc -O2 -fopenmp -o merge_sort merge_sort_omp.c -lm

it_check_1() {
	test `OMP_NUM_THREADS=1 ./merge_sort 1 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_2() {
	test `OMP_NUM_THREADS=1 ./merge_sort 2 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_3() {
	test `OMP_NUM_THREADS=1 ./merge_sort 3 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_10() {
	test `OMP_NUM_THREADS=1 ./merge_sort 10 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_100() {
	test `OMP_NUM_THREADS=1 ./merge_sort 100 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_100() {
	test `OMP_NUM_THREADS=1 ./merge_sort 1000 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_10_2() {
	test `OMP_NUM_THREADS=2 ./merge_sort 10 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_10_7() {
	test `OMP_NUM_THREADS=7 ./merge_sort 10 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_10_8() {
	test `OMP_NUM_THREADS=8 ./merge_sort 10 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_10_10() {
	test `OMP_NUM_THREADS=10 ./merge_sort 10 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_100_2() {
	test `OMP_NUM_THREADS=2 ./merge_sort 100 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_100_7() {
	test `OMP_NUM_THREADS=7 ./merge_sort 100 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_100_8() {
	test `OMP_NUM_THREADS=8 ./merge_sort 100 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_100_10() {
	test `OMP_NUM_THREADS=10 ./merge_sort 100 | tail -n 1 | awk '{print $5}'` "=" 1
}
