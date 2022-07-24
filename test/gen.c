#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<time.h>
int main(int argc, char *argv[]) {
	srand(time(0));
	int M = atoi(argv[1]);
	int N = atoi(argv[2]);
	FILE *f = fopen("test.in", "w");
	for(int i = 0; i < N * M + N + M; i++)
		fprintf(f, "%f ", (float)rand());
	fprintf(f, "\n");
	return 0;
}
