#include "soccs_sce_tools.h"

uint32_t pthread_t_to_uint32_t(pthread_t pthread_t_id) {
  uint8_t *p = (uint8_t*) &pthread_t_id;
  uint32_t result = 0;
  uint8_t *pp = (uint8_t*) &result;
  for (uint64_t i = 0; i < sizeof(uint32_t); i++) {
    *pp = *p;
    pp = pp + 1;
    p = p + 1;
  }
  return result;
}

void calc_timespec_difference(
    const struct timespec *p_a, 
    const struct timespec *p_b, 
    struct timespec *p_result) {
  if (p_a->tv_sec == p_b->tv_sec) {
    p_result->tv_sec = 0;
    if (p_a->tv_nsec >= p_b->tv_nsec) {
      p_result->tv_nsec = 
        p_a->tv_nsec - p_b->tv_nsec;
    } else {
      perror("p_a->tv_sec == p_b->tv_sec && p_a->tv_nsec < p_b->tv_nsec");
      exit(EXIT_FAILURE);
    }
  } else if (p_a->tv_sec > p_b->tv_sec) {
    p_result->tv_sec = p_a->tv_sec - p_b->tv_sec;
    p_result->tv_nsec = p_a->tv_nsec - p_b->tv_nsec;
    if (p_result->tv_nsec < 0) {
      p_result->tv_sec -= 1;
      p_result->tv_nsec += ONE_BILLION;
    }
  } else {
    perror("p_a->tv_sec < p_b->tv_sec");
    exit(EXIT_FAILURE);
  }
}

void calc_timespec_sum(
    const struct timespec *p_a, 
    const struct timespec *p_b, 
    struct timespec *p_result) {
  p_result->tv_sec = p_a->tv_sec + p_b->tv_sec;
  p_result->tv_nsec = p_a->tv_nsec + p_b->tv_nsec;
  while (p_result->tv_nsec >= ONE_BILLION) {
    p_result->tv_sec += 1;
    p_result->tv_nsec -= ONE_BILLION;
  }
}

void calc_timespec_sum (
    const struct timespec *p_a, 
    struct timespec *p_result) {
  p_result->tv_sec += p_a->tv_sec;
  p_result->tv_nsec += p_a->tv_nsec ;
  while (p_result->tv_nsec >= ONE_BILLION) {
    p_result->tv_sec += 1;
    p_result->tv_nsec -= ONE_BILLION;
  }
}

double PYTHAG(double a, double b) {
  double at = fabs(a), bt = fabs(b), ct, result;
	if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
	else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
	else result = 0.0;
	return(result);
}

