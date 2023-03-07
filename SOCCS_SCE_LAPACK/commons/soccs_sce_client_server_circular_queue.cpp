#include "soccs_sce_client_server_circular_queue.h"

bool can_malloc_straight(
    struct circular_queue *p_circular_queue,
    int count_of_bytes)
{
  if (count_of_bytes < 0)
  {
    throw "data size must be non-negative";
  }

  if (!can_enqueue_with_circling_allowed(p_circular_queue, count_of_bytes))
    return false;
  if (p_circular_queue->meta_info.tail > p_circular_queue->meta_info.head)
  {
    if (QUEUE_SIZE - p_circular_queue->meta_info.tail - 1 >= (uint64_t)count_of_bytes)
      return true;
    else if (p_circular_queue->meta_info.head - 1 >= (uint64_t)count_of_bytes)
      return true;
    else
      return false;
  }
  else if (p_circular_queue->meta_info.tail == p_circular_queue->meta_info.head)
  {
    p_circular_queue->meta_info.tail = p_circular_queue->meta_info.head = 0;
    p_circular_queue->data[p_circular_queue->meta_info.tail] = 0;
    return (QUEUE_SIZE - 1 >= count_of_bytes);
  }
  else
    return (p_circular_queue->meta_info.head >
        p_circular_queue->meta_info.tail + count_of_bytes);
}

uint8_t *malloc_straight(
    struct circular_queue *p_circular_queue,
    int count_of_bytes /* number of bytes to malloc */)
{
  if (count_of_bytes < 0)
    throw "data size must be non-negative";

  if (!can_enqueue_with_circling_allowed(p_circular_queue, (uint64_t)count_of_bytes))
    return NULL;

  // by now, enough free bytes are available, if circling is allowed.

  if (p_circular_queue->meta_info.tail > p_circular_queue->meta_info.head)
  {
    // e.g. a 10 byte circular queue:
    //| | | | | |h|x|t| | |
    if (QUEUE_SIZE - p_circular_queue->meta_info.tail - 1 >= (uint64_t)count_of_bytes)
    {
      uint8_t *result = ((uint8_t *)(p_circular_queue->data)) + p_circular_queue->meta_info.tail;
      p_circular_queue->meta_info.tail += count_of_bytes;
      // the byte pointed to by tail should be zero valued.
      p_circular_queue->data[p_circular_queue->meta_info.tail] = 0;
      return result;
    }
    else if (p_circular_queue->meta_info.head - 1 >= (uint64_t)count_of_bytes)
    {
      // add padding till the end of the circular queue block.
      for (int i = p_circular_queue->meta_info.tail; i < QUEUE_SIZE; i++)
        p_circular_queue->data[i] = 0;
      uint8_t *result = (uint8_t *)(p_circular_queue->data);
      p_circular_queue->meta_info.tail = count_of_bytes;
      // the byte pointed to by tail should be zero valued.
      p_circular_queue->data[p_circular_queue->meta_info.tail] = 0;
      return result;
    }
    else
      return NULL;
  }
  else if (p_circular_queue->meta_info.tail == p_circular_queue->meta_info.head)
  {
    // queue empty, reset
    p_circular_queue->meta_info.tail = p_circular_queue->meta_info.head = 0;
    p_circular_queue->data[p_circular_queue->meta_info.tail] = 0;
    if (QUEUE_SIZE - 1 >= count_of_bytes)
    {
      p_circular_queue->meta_info.tail = count_of_bytes;
      p_circular_queue->data[p_circular_queue->meta_info.tail] = 0;
      return (uint8_t *)(p_circular_queue->data);
    }
    else
      return NULL;
  }
  else
  {
    if (p_circular_queue->meta_info.head > p_circular_queue->meta_info.tail +
        (uint64_t)count_of_bytes)
    {
      uint8_t *result = ((uint8_t *)(p_circular_queue->data)) + p_circular_queue->meta_info.tail;
      p_circular_queue->meta_info.tail += count_of_bytes;
      p_circular_queue->data[p_circular_queue->meta_info.tail] = 0;
      return result;
    }
    else
      return NULL;
  }
}

bool can_shallow_dequeue_straight(struct circular_queue *p_circular_queue)
{
  if (p_circular_queue->meta_info.tail == p_circular_queue->meta_info.head)
  {
    p_circular_queue->meta_info.tail = p_circular_queue->meta_info.head = 0;
    p_circular_queue->data[p_circular_queue->meta_info.tail] = 0;
    return false;
  }
  while (p_circular_queue->data[p_circular_queue->meta_info.head] == 0)
  {
    p_circular_queue->meta_info.head = (p_circular_queue->meta_info.head + 1) % QUEUE_SIZE;
    if (p_circular_queue->meta_info.head == p_circular_queue->meta_info.tail)
    {
      p_circular_queue->meta_info.tail = p_circular_queue->meta_info.head = 0;
      p_circular_queue->data[p_circular_queue->meta_info.tail] = 0;
      return false;
    }
  }
  uint8_t *p_packet = ((uint8_t *)(p_circular_queue->data)) + p_circular_queue->meta_info.head;
  struct circular_queue_packet_common_header *p_common_header =
    (struct circular_queue_packet_common_header *)p_packet;
  assert(p_common_header->preamble == PREAMBLE);
  uint32_t packet_size = p_common_header->packet_size;

#ifdef DEBUGGING
  assert(p_circular_queue->meta_info.head + packet_size <= QUEUE_SIZE - 1);
  assert(packet_size > 0);
#endif

  return true;
}

uint8_t *shallow_dequeue_straight(struct circular_queue *p_circular_queue, int *size)
{
  if (p_circular_queue->meta_info.tail == p_circular_queue->meta_info.head)
  {
    p_circular_queue->meta_info.tail = p_circular_queue->meta_info.head = 0;
    p_circular_queue->data[p_circular_queue->meta_info.tail] = 0;
    (*size) = 0;
    return NULL;
  }
  while (p_circular_queue->data[p_circular_queue->meta_info.head] == 0)
  {
    p_circular_queue->meta_info.head = (p_circular_queue->meta_info.head + 1) % QUEUE_SIZE;
    if (p_circular_queue->meta_info.head == p_circular_queue->meta_info.tail)
    {
      p_circular_queue->meta_info.tail = p_circular_queue->meta_info.head = 0;
      p_circular_queue->data[p_circular_queue->meta_info.tail] = 0;
      (*size) = 0;
      return NULL;
    }
  }
  uint8_t *p_packet = ((uint8_t *)(p_circular_queue->data)) + p_circular_queue->meta_info.head;
  struct circular_queue_packet_common_header *p_common_header =
    (struct circular_queue_packet_common_header *)p_packet;
  assert(p_common_header->preamble == PREAMBLE);
  uint32_t packet_size = p_common_header->packet_size;

#ifdef DEBUGGING
  assert(p_circular_queue->meta_info.head + packet_size <= QUEUE_SIZE - 1);
  assert(packet_size > 0);
#endif

  p_circular_queue->meta_info.head = (p_circular_queue->meta_info.head + packet_size) % QUEUE_SIZE;
  if (p_circular_queue->meta_info.head == p_circular_queue->meta_info.tail)
  {
    p_circular_queue->meta_info.tail = p_circular_queue->meta_info.head = 0;
    p_circular_queue->data[p_circular_queue->meta_info.tail] = 0;
  }
  (*size) = packet_size;
  return p_packet;
}

uint8_t *peek_shallow_dequeue_straight(struct circular_queue *p_circular_queue, int *size)
{
  if (p_circular_queue->meta_info.tail == p_circular_queue->meta_info.head)
  {
    p_circular_queue->meta_info.tail = p_circular_queue->meta_info.head = 0;
    // the byte pointed to by tail should be zero valued.
    p_circular_queue->data[p_circular_queue->meta_info.tail] = 0;
    (*size) = 0;
    return NULL;
  }
  // queue is not empty

  uint64_t old_tail = p_circular_queue->meta_info.tail;
  uint64_t old_head = p_circular_queue->meta_info.head;

  while (p_circular_queue->data[p_circular_queue->meta_info.head] == 0)
  {
    p_circular_queue->meta_info.head = (p_circular_queue->meta_info.head + 1) % QUEUE_SIZE;
    if (p_circular_queue->meta_info.head == p_circular_queue->meta_info.tail)
    {
      p_circular_queue->meta_info.tail = p_circular_queue->meta_info.head = 0;
      p_circular_queue->data[p_circular_queue->meta_info.tail] = 0;
      (*size) = 0;
      return NULL;
    }
  }
  uint8_t *p_packet = ((uint8_t *)(p_circular_queue->data)) + p_circular_queue->meta_info.head;
  struct circular_queue_packet_common_header *p_common_header =
    (struct circular_queue_packet_common_header *)p_packet;
  assert(p_common_header->preamble == PREAMBLE);
  uint32_t packet_size = p_common_header->packet_size;

#ifdef DEBUGGING
  assert(p_circular_queue->meta_info.head + packet_size <= QUEUE_SIZE - 1);
  assert(packet_size > 0);
#endif

  (*size) = packet_size;

  p_circular_queue->meta_info.head = old_head; // necesary, see the above while loop.
  p_circular_queue->meta_info.tail = old_tail; // unnecesary, put here just for extra
  return p_packet;
}

