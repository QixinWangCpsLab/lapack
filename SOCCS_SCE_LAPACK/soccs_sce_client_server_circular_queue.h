#ifndef _CIRCULAR_QUEUE_H_
#define _CIRCULAR_QUEUE_H_

#include <stdint.h>
#include <pthread.h>
#include <assert.h>

/**
 * Unit: bytes
 *
 * The size of a page frame (aka physical memory page).
 */
#define PAGE_FRAME_SIZE (4096)

/**
 * No more than one page frame (aka physical memory page),
 * unit: byte. Must be big enough to hold both the c2s_queue
 * and the s2c_queue.
 * 101 for c2s and 101 for s2c
 */
#define SHARED_MEM_SIZE (PAGE_FRAME_SIZE * 202) 

/**
 * Unit: bytes
 *
 * @see Comments on the struct circular_queue. 
 * The maximum circular queue packet size includes the 
 * circular_queue_packet_common_header size.
 */
#define MAX_CIRCULAR_QUEUE_PACKET_SIZE (PAGE_FRAME_SIZE * 50)

/**
 * The preamble must not contain zero valued byte(s), to differentiate it from a 
 * (zero valued) padding byte.
 */
#define PREAMBLE (0xAAAAAAAA)

/**
 * This is the common header of all circular queue packets.
 * 
 * @see Comments on struct circular_queue.
 */
struct circular_queue_packet_common_header {
  uint32_t preamble = PREAMBLE;
  uint32_t packet_size;
} __attribute__((packed));

/**
 * This is the common meta info stored at the begining of
 * a circular queue.
 *
 * @see Comments on struct circular_queue.
 */
struct circular_queue_meta_info {
  uint64_t head = 0;
  uint64_t tail = 0;
  pthread_mutex_t mutex;
  pthread_cond_t can_produce;
  pthread_cond_t can_consume;
} __attribute__((packed));

/**
 * Unit: bytes 
 *
 * The size of a circular_queue's data block.
 *
 * Must ensure the c2s_circular_queue size  + s2c_circular_queue size does 
 * not exceed the SHARED_MEM_SIZE.
 * Must ensure QUEUE_SIZE > 2 x MAX_CIRCULAR_QUEUE_PACKET_SIZE, where 
 * MAX_CIRCULAR_QUEUE_PACKET_SIZE is the maximum circular queue packet size 
 * (including the circular_queue_packet_common_header). 
 * This is to ensure ability to add a MAX_CIRCULAR_QUEUE_PACKET_SIZE
 * circular queue packet straightly, no matter where the head and tail is 
 * pointing at (when the queue is empty).
 */
const int QUEUE_SIZE = ((SHARED_MEM_SIZE / 2) - sizeof(struct circular_queue_meta_info));

/**
 * This circular queue can only hold packets (referred to as ``circular queue packet'') 
 * of the following format:
 * 
 *     struct circular_queue_packet_common_header {
 *         uint32_t preamble that must be the nonzero value of PREAMBLE
 *         uint32_t packet_size (in bytes) including the preamble
 *     }
 *     ...     cutomized remainder of the packet header
 *     uint8_t payload[packet_size - packet header size]
 *
 * The nonzero preamble is to differentiate the packet from padding bytes.
 * Correspondingly, padding bytes must be zero, and the byte pointed to by
 * meta_info.tail must also be zero. The byte pointed to by meta_info.tail 
 * can never store any content.
 */
struct circular_queue {
  struct circular_queue_meta_info meta_info;
  uint8_t data[QUEUE_SIZE];
} __attribute__((packed));

/**
 * As the mutex is not yet initialized, we cannot lock before calling this function.
 * Make sure there is only one thread calling this function.
 *
 * Every byte in the data block is zeroed. This is to ensure all 
 * non-circular-queue-packet contents are regarded as padding, 
 * and to ensure meta_info.tail points to a zero valued byte.
 */
inline void init_circular_queue(struct circular_queue *p_circular_queue) {

  for (int i = 0; i < QUEUE_SIZE; i++)
    p_circular_queue->data[i] = 0;

  p_circular_queue->meta_info.head = 0;
  p_circular_queue->meta_info.tail = 0;

  pthread_mutexattr_t mutex_attr;
  pthread_mutexattr_init(&mutex_attr);
  pthread_mutexattr_setpshared(&mutex_attr, PTHREAD_PROCESS_SHARED);
  pthread_mutex_init(&(p_circular_queue->meta_info.mutex), &mutex_attr);

  pthread_condattr_t cond_can_produce_attr;
  pthread_condattr_init(&cond_can_produce_attr);
  pthread_condattr_setpshared(&cond_can_produce_attr, PTHREAD_PROCESS_SHARED);
  pthread_cond_init(&(p_circular_queue->meta_info.can_produce), &cond_can_produce_attr);

  pthread_condattr_t cond_can_consume_attr;
  pthread_condattr_init(&cond_can_consume_attr);
  pthread_condattr_setpshared(&cond_can_consume_attr, PTHREAD_PROCESS_SHARED);
  pthread_cond_init(&(p_circular_queue->meta_info.can_consume), &cond_can_consume_attr);
}


/**
 * Prerequisite: must have locked the mutex before calling this function! I.e. call 
 *     pthread_mutex_lock(&p_circular_queue->meta_info.mutex)
 * before calling this function.
 *
 * @return the data size (in bytes) of the queue, i.e. bytes occupied by data in 
 *         the data block of the queue.
 */
inline int len(struct circular_queue *p_circular_queue) {
  return (p_circular_queue->meta_info.tail 
      + QUEUE_SIZE - p_circular_queue->meta_info.head) % QUEUE_SIZE;
}

/**
 * Prerequisite: must have locked the mutex before calling this function! I.e. call
 *     pthread_mutex_lock(&p_circular_queue->meta_info.mutex)
 * before calling this function.
 *
 * @return the number of free bytes in the data block of this queue.
 */
inline int num_free_bytes(struct circular_queue *p_circular_queue) {
  return QUEUE_SIZE - 1 - len(p_circular_queue); //the 1 extra byte is because
                                                 //the meta_info.tail must be pointing to an unused byte of value 0.
}

/**
 * Prerequisite: must have locked the mutex before calling this function! I.e. call
 *     pthread_mutex_lock(&p_circular_queue->meta_info.mutex)
 * before calling this function.
 *
 * @return maximum number of bytes can be enqueued without causing circling (using 
 * conventional circular queue enqueueing method).
 */
inline int num_free_bytes_till_end_of_circular_queue_data_block(
    struct circular_queue *p_circular_queue) {
  if (p_circular_queue->meta_info.tail >= p_circular_queue->meta_info.head)
    return QUEUE_SIZE - p_circular_queue->meta_info.tail - 1;
  else //p_circular_queue->meta_info.tail < p_circular_queue->meta_info.head
    return p_circular_queue->meta_info.head - p_circular_queue->meta_info.tail - 1;
}


/**
 * Prerequisite: must have locked the mutex before calling this function! I.e. call
 *     pthread_mutex_lock(&this->mutex) 
 * before calling this function.
 *
 * @return true if there is enough space to add count_of_bytes (must be >=
 *     0) bytes of data to the circular queue (circling allowed). Return false 
 *     if otherwise.
 * @param count_of_bytes must be >= 0; or an exception happens.
 */
inline bool can_enqueue_with_circling_allowed(struct circular_queue *p_circular_queue, 
    int count_of_bytes) {
  if (count_of_bytes < 0) {
    throw "data size must be non-negetive";
    // perror("data size must be non-negetive\n");
    // exit(-1);
  }
  return len(p_circular_queue) + count_of_bytes <= QUEUE_SIZE - 1;
}

/**
 * Prerequisite: must have locked the mutex before calling this function! I.e. call
 *     pthread_mutex_lock(&p_circular_queue->meta_info.mutex) 
 * before calling this function.
 *
 * Assumption: a meaningful circular queue packet's initial byte (i.e. preamble) 
 * must be nonzero.
 *
 * Assumption: paddings are zero valued uint_8 bytes to differentiate from 
 * a circular queue packet preamble.
 *
 * Assumption: the byte pointed to by the tail is always zero valued.
 *
 * Note if the queue is found to be empty, this function will set 
 * p_circular_queue->meta_info.tail = p_circular_queue->meta_info.head = 0.
 *
 * @return true if there is enough space to malloc count_of_bytes (must be >=
 *     0) bytes of data from the circular queue without circling these bytes. Return 
 *     false if otherwise.
 * @param count_of_bytes must be >= 0; or an exception happens.
 */
extern bool can_malloc_straight(struct circular_queue *p_circular_queue, 
    int count_of_bytes);

/**
 * Prerequisite: must have locked the mutex before calling this function! I.e. call
 *     pthread_mutex_lock(&p_circular_queue->meta_info.mutex) 
 * before calling this function.
 *
 * Assumption: a meaningful circular queue packet's initial byte (i.e. preamble) 
 * must be nonzero.
 *
 * Assumption: paddings are zero valued uint_8 bytes to differentiate from 
 * a circular queue packet preamble.
 *
 * Assumption: the byte pointed to by the tail is always zero valued.
 *
 * Note if the queue is found to be empty, this function will set 
 * p_circular_queue->meta_info.tail = p_circular_queue->meta_info.head = 0, and then
 * do the malloc.
 *
 * @return The pointer to the allocated space; or NULL otherwise.
 */
extern uint8_t *malloc_straight(struct circular_queue *p_circular_queue, 
    int count_of_bytes /* number of bytes to malloc */);

/**
 * Prerequisite: must have locked the mutex before calling this function! I.e. call
 *     pthread_mutex_lock(&p_circular_queue->meta_info.mutex) 
 * before calling this function.
 *
 * Assumption: a meaningful circular queue packet's initial byte (i.e. preamble) 
 * must be nonzero.
 *
 * Assumption: paddings are zero valued uint_8 bytes to differentiate from 
 * a circular queue packet preamble.
 *
 * Assumption: the byte pointed to by the tail is always zero valued.
 *
 * Note if the queue is found to be empty, this function will set 
 * p_circular_queue->meta_info.tail = p_circular_queue->meta_info.head = 0.
 *
 * @return true if there is a straight circular queue packet to shallow dequeue 
 *     (shallow dequeue means the bytes are not actually copied out of the circular 
 *     queue data block); false otherwise.
 */
extern bool can_shallow_dequeue_straight(struct circular_queue *p_circular_queue);

/**
 * Prerequisite: must have locked the mutex before calling this function! I.e. call
 *     pthread_mutex_lock(&p_circular_queue->meta_info.mutex) 
 * before calling this function.
 *
 * Assumption: a meaningful circular queue packet's initial byte (i.e. preamble) 
 * must be nonzero.
 *
 * Assumption: paddings are zero valued uint_8 bytes to differentiate from 
 * a circular queue packet preamble.
 *
 * Assumption: the byte pointed to by the tail is always zero valued.
 *
 * Note1: if the queue is found to be empty, this function will set 
 * p_circular_queue->meta_info.tail = p_circular_queue->meta_info.head = 0.
 *
 * Note2: after returned from this function, the packet is no longer considered part
 * of the enqueued data, hence can be overwritten by other enqueuing threads. It is
 * up to the user to carry out necessary mutual exclusion. If you do not want any 
 * damaging writing, use @see peek_shallow_dequeue_straight instead.
 * 
 * @param (*size) is an output: number of bytes of the shallow dequeued (i.e. 
 *     the bytes are not actually copied out of the circular queue data block) 
 *     straight circular queue packet from this circular queue, 
 *     0 if no such packet dequeued.
 * @return pointer to the preamble of the shallow dequeued straight 
 *     circular queue packet (i.e. the bytes are not actually copied out of the 
 *     circular queue data block). Note this pointer is only valid when the 
 *     outputed (*size) > 0. NULL otherwise. Denote this returned pointer as p.
 *     If (p != NULL and (*size) > 0) then 
 *     p->circular_queue_packet_common_header.packet_size == (*size).
 */
extern uint8_t * shallow_dequeue_straight(struct circular_queue *p_circular_queue, 
    int* size);

/**
 * Prerequisite: must have locked the mutex before calling this function! I.e. call
 *     pthread_mutex_lock(&p_circular_queue->meta_info.mutex) 
 * before calling this function.
 * 
 * Assumption: a meaningful circular queue packet's initial byte (i.e. preamble) 
 * must be nonzero.
 *
 * Assumption: paddings are zero valued uint_8 bytes to differentiate from 
 * a circular queue packet preamble.
 *
 * Assumption: the byte pointed to by the tail is always zero valued.
 *
 * @see Same as shallow_dequeue_straight, except that p_circular_queue->meta_info.head
 *     is not actually updated (except the Note1 below), i.e. there is no actual dequeue.
 * 
 * Note1: if the queue is found to be empty, this function will set 
 * p_circular_queue->meta_info.tail = p_circular_queue->meta_info.head = 0.
 *
 * Note2: unlike Note2 of shallow_dequeue_straight, after returned from this function, 
 * the packet is still considered part of the enqueued data, hence cannot be 
 * overwritten by other enqueuing threads. 
 * 
 */
extern uint8_t * peek_shallow_dequeue_straight(struct circular_queue *p_circular_queue, 
    int* size);

#endif

// code review reached here on 22/9/2022
