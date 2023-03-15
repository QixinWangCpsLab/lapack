#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <pthread.h>
#include <time.h>
#include <math.h>

#define max(a, b) \
  ({ __typeof__ (a) _a = (a); \
   __typeof__ (b) _b = (b); \
   _a > _b ? _a : _b; })

#define min(a, b) \
  ({ __typeof__ (a) _a = (a); \
   __typeof__ (b) _b = (b); \
   _a > _b ? _b : _a; })
   
#define SIGN(a, b) ((b) >= 0.0 ? \
   fabs(a): -fabs(a))

#define EPSILON 0.02
#define ONE_BILLION 1000000000
#define ONE_HUNDRED_THOUSAND 100000
#define TEN_MILLION 10000000

const struct timespec wait_timeout = {.tv_sec = 0, .tv_nsec = TEN_MILLION};

uint32_t pthread_t_to_uint32_t(pthread_t pthread_t_id);

/**
 * @assumption *p_a >= *p_b.
 * or the program exits with EXIT_FAILURE.
 */
void calc_timespec_difference(
    const struct timespec *p_a, 
    const struct timespec *p_b, 
    struct timespec *p_result);

void calc_timespec_sum(
    const struct timespec *p_a, 
    const struct timespec *p_b, 
    struct timespec *p_result);

void calc_timespec_sum(
    const struct timespec *p_a,
    struct timespec *p_result);

inline void copy_timespec(
    const struct timespec *src,
    struct timespec *des) {
  if (src == NULL || des == NULL) return;
  des->tv_sec = src->tv_sec;
  des->tv_nsec = src->tv_nsec;
}

double PYTHAG(double a, double b);
