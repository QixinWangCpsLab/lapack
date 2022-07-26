#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<time.h>
#include<sys/time.h>

int main(int argc, char *argv[]) {
	struct timeval tv;
	gettimeofday(&tv, NULL);
	srand(tv.tv_usec);
	int M = atoi(argv[1]);
	int N = atoi(argv[2]);
	FILE *f = fopen("test.in", "w");
	fprintf(f, "%d ", M);
	fprintf(f, "%d ", N);
	for(int i = 0; i < N * M + N + M; i++)
		fprintf(f, "%f ", (float)rand());
	fprintf(f, "\n");
	return 0;
}
