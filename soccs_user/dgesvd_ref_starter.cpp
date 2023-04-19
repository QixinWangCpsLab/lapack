#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/wait.h>

#define A_FILE "data/dgesvd_a_200x200.data"
#define EXEC_PATH "./dgesvd_user_refer"

int main(int argc, char *argv[]) {
  
  if(argc < 2) {
    fprintf(stderr, "Usage: %s <pr_1> <pr_2> ...\n", argv[0]);
    exit(EXIT_FAILURE);
  }

  int num_clients = argc - 1;
  for(int i = 1; i <= num_clients; i++) {
    if(fork() == 0) {
      fprintf(stderr, "Subtask %d running on core 0 - P=%s.\n", getpid(), argv[i]);
      execlp(EXEC_PATH, EXEC_PATH, "0", "-1", argv[i], NULL);
      fprintf(stderr, "fail: errno: %s\n", strerror(errno));
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

