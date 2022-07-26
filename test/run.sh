#!/bin/bash
for ((i=0;i<3;i++));
	do
		rm -f test.in
		./gen $1 $2
		./a.out `cat test.in` > test1.out
		./origin_a.out `cat test.in` > test2.out
		echo "---------------------------------------"
		cat test1.out
		echo "---------------------------------------"
		cat test2.out
		echo "---------------------------------------"
		diff test1.out test2.out
		sleep 0.001
	done
