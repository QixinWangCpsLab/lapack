#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <sys/wait.h>

#define AP_FILE "data/dppsv_ap_200x200.data"
#define BX_FILE "data/dppsv_bx_200x200.data"
#define EXEC_PATH "./dppsv_user_refer"

int main(int argc, char *argv[]) {
  
  if(argc < 2) {
    fprintf(stderr, "Usage: ./dppsv_ref_starter <index_1> <index_2> ...\n");
    exit(EXIT_FAILURE);
  }

  int num_clients = argc - 1;
  for(int i = 0; i < num_clients; i++) {
    int index = atoi(argv[i]);
    if(fork() == 0) {
      fprintf(stderr, "Subtask %d running on core %d.\n", getpid(), index);
      execlp(EXEC_PATH, EXEC_PATH, argv[i], AP_FILE, BX_FILE, NULL);
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

