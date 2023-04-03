#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/wait.h>

#define A_FILE "data/dgesvd_a_200x200.data"
#define EXEC_PATH "./dgesvd_user_refer"

int main(int argc, char *argv[]) {
  
  if(argc != 2) {
    fprintf(stderr, "Usage: %s <num> ...\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  int num_clients = atoi(argv[1]);
  for(int i = 0; i < num_clients; i++) {
    if(fork() == 0) {
      fprintf(stderr, "Subtask %d running on core %d.\n", getpid(), 0);
      execlp(EXEC_PATH, EXEC_PATH, "0", "-1", A_FILE, NULL);
      fprintf(stderr, "errno: %s\n", strerror(errno));
    }
  }

  fprintf(stderr, "Waiting %d to finish.\n", num_clients);
  int alive = num_clients;
  while(alive > 0) {
    int status;
    wait(&status);
    if(WIFEXITED(status) == true)
      alive--;
  }
  fprintf(stderr, "Terminated!\n");

  return 0;  
}

