#!./roundup

mpicc -O2 -o merge_sort mpi_merge_sort_task.c -lm

it_check_1() {
	test `mpirun -np 1 ./merge_sort 1 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_2() {
	test `mpirun -np 1 ./merge_sort 2 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_3() {
	test `mpirun -np 1 ./merge_sort 3 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_10() {
	test `mpirun -np 1 ./merge_sort 10 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_100() {
	test `mpirun -np 1 ./merge_sort 100 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_1000() {
	test `mpirun -np 1 ./merge_sort 1000 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_100_3() {
	test `mpirun -np 3 ./merge_sort 100 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_100_5() {
	test `mpirun -np 5 ./merge_sort 100 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_100_8() {
	test `mpirun -np 8 ./merge_sort 100 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_100_13() {
	test `mpirun -np 13 ./merge_sort 100 | tail -n 1 | awk '{print $5}'` "=" 1
}
