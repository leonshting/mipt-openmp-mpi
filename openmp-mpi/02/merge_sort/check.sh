#!./roundup

gcc -O2 -o merge_sort merge_sort.c

it_check_1() {
	test `./merge_sort 1 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_2() {
	test `./merge_sort 2 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_3() {
	test `./merge_sort 3 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_10() {
	test `./merge_sort 10 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_100() {
	test `./merge_sort 100 | tail -n 1 | awk '{print $5}'` "=" 1
}
it_check_100() {
	test `./merge_sort 1000 | tail -n 1 | awk '{print $5}'` "=" 1
}
