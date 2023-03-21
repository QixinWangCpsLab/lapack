#include "soccs_sce_client_server_circular_queue.h"
#include "soccs_sce_client_server_protocol.h"
#include "soccs_sce_tools.h"
#include "soccs_sce_server_lapack_dppsv.h"

using namespace std;
using namespace std::chrono;

static int fd_shm;
static void *p_shm = nullptr;
static struct circular_queue *request_queue = nullptr;
static struct circular_queue *reply_queue = nullptr;
static pthread_t server_pthread;
static bool terminated = false;

static void sigint_handler(int sig)
{
  if (sig != SIGINT)
    return;
  terminated = true;
  fprintf(stderr, "SIGINT received!\n");
}

static void *server_soccs(void *p_arg)
{
  prctl(PR_SET_NAME, "soccs_server");
  /* server soccs_sce_lapack_dppsv */

  struct sched_param server_pthread_sched_param = {
    .sched_priority = sched_get_priority_max(SCHED_RR)
  };
  pthread_setschedparam(pthread_self(), SCHED_RR,
      &server_pthread_sched_param);

  struct circular_queue *request_queue =
    ((struct soccs_sce_server_pthread_argument *)p_arg)->p_c2s_queue;
  struct circular_queue *reply_queue =
    ((struct soccs_sce_server_pthread_argument *)p_arg)->p_s2c_queue;
  auto start = system_clock::now(), end = system_clock::now();
  auto duration = duration_cast<nanoseconds> (end - start);

  try
  {
    while (terminated == false)
    {
      /* work loop */

      // TODO timer

      /* acquire request queue lock */
      pthread_mutex_lock(&request_queue->meta_info.mutex);

#ifdef DEBUGGING
      fprintf(stderr, "server: acquire request queue.\n");
#endif

      while (!can_shallow_dequeue_straight(request_queue))
      {
        /* wait until any request packets or timeout */
        struct timespec abs_timeout = {
          .tv_sec = 0, .tv_nsec = 0};
        if (clock_gettime(CLOCK_REALTIME, &abs_timeout) < 0)
          throw "clock_gettime failed.\n";
        calc_timespec_sum(&wait_timeout, &abs_timeout);
        pthread_cond_timedwait(&request_queue->meta_info.can_consume,
            &(request_queue->meta_info.mutex), &abs_timeout);

        if (terminated == true)
          break;
      }

      if (terminated == true)
      {
        pthread_mutex_unlock(&(request_queue->meta_info.mutex));
        break;
      }

#ifdef DEBUGGING
      fprintf(stderr, "server: can shallow dequeue straightly now.\n");
#endif

      /* receive request packet header */
      int request_packet_size = 0;
      struct request_packet_header *request_header =
        (struct request_packet_header *)
        shallow_dequeue_straight(request_queue, &request_packet_size);

      if (request_packet_size <= 0 || request_header == nullptr)
        throw "request_packet_size <= 0 || request_header == nullptr";

      /* work */
      switch (request_header->function_id)
      {
        case LAPACK_DPPSV:
          {

#ifdef DEBUGGING
            fprintf(stderr, "received LAPACK_DPPSV from client %u\n",
                request_header->client_pthread_id);
#endif
            struct request_lapack_dppsv_fixed *req_pkt_fixed =
              reinterpret_cast<struct request_lapack_dppsv_fixed *>
              (request_header + 1);
            const char uplo         = req_pkt_fixed->uplo;
            const lapack_int n      = req_pkt_fixed->n;
            const lapack_int nrhs   = req_pkt_fixed->nrhs;
            const lapack_int ldb    = req_pkt_fixed->ldb;
            const lapack_int len_AP = n * (n + 1) / 2;
            const lapack_int len_B  = ldb * nrhs;
            lapack_int info = 0;
            double *ap = (double *) malloc (len_AP * sizeof(double));
            double *b = (double *) malloc (len_B * sizeof(double));
            double *req_pkt_flexible =
              reinterpret_cast<double *>(req_pkt_fixed + 1);
            double *req_pkt_flexible_pos = req_pkt_flexible;

            memcpy(ap, req_pkt_flexible_pos, len_AP * sizeof(double));
            req_pkt_flexible_pos += len_AP;
            memcpy(b, req_pkt_flexible_pos, len_B * sizeof(double));
            /*
#ifdef DEBUGGING
fprintf(stderr, "ap[%d]={", len_AP);
for (int i = 0; i < len_AP; i++)
fprintf(stderr, " %lf", ap[i]);
fprintf(stderr, " }\n");
fprintf(stderr, "b[%d]={", len_B);
for (int i = 0; i < len_B; i++)
fprintf(stderr, " %lf", b[i]);
fprintf(stderr, " }\n");
#endif
             */
            /* sanity checking */
            if (uplo != 'U' && uplo != 'u' && uplo != 'L' && uplo != 'l')
              info = -1;
            else if (n < 0)
              info = -2;
            else if (nrhs < 0)
              info = -3;
            else if (ldb < max(1, n))
              info = -6;
            if (info != 0)
              goto reply_dppsv; // no need to calculate

            /* calculation starts from here */
            // start            /* 190: dpptrf( uplo, n, ap) */
            start = system_clock::now();
            if (uplo == 'L' || uplo == 'l')
            {
              ap[0] = sqrt(ap[0]);
              for (int i = 1; i < n; i++)
                ap[i] = ap[i] / ap[0];
              double sum = 0;
              for (int j = 2; j <= n; j++)
              {
                sum = 0;
                for (int k = 1; k < j; k++)
                  sum += ap[j+(k-1)*(2*n-k)/2-1] * ap[j+(k-1)*(n*2-k)/2-1];
                ap[j+(j-1)*(2*n-j)/2-1] = sqrt(ap[j+(j-1)*(2*n-j)/2-1] - sum);
                for (int i = j + 1; i <= n; i++)
                {
                  sum = 0;
                  for (int k = 1; k <= j - 1; k++)
                  {
                    sum += ap[i+(k-1)*(2*n-k)/2-1] * ap[j+(k-1)*(2*n-k)/2-1];
                  }
                  ap[i+(j-1)*(2*n-j)/2-1] = 
                    (ap[i+(j-1)*(2*n-j)/2-1] - sum) / ap[j+(j-1)*(2*n-j)/2-1];
                }
              }
              // solve LLT=b
              b[0] = b[0] / ap[0];
              for (int i = 2; i <= n; i++)
              {
                sum = 0;
                for (int j = 1; j <= i - 1; j++)
                  sum += ap[i+(j-1)*(2*n-j)/2-1] * b[j-1];
                b[i-1] = (b[i-1] - sum) / ap[i+(i-1)*(2*n-i)/2-1];
              }
              b[n-1] = b[n-1] / ap[n+n*(n-1)/2-1];
              for (int i = n - 1; i > 0; i--)
              {
                sum = 0;
                for (int j = i + 1; j <= n; j++)
                  sum += ap[j+(i-1)*(2*n-i)/2-1] * b[j-1];
                b[i-1] = (b[i-1] - sum) / ap[i+(i-1)*(2*n-i)/2-1];
              }
            }
            end = system_clock::now();
            duration = duration_cast<nanoseconds> (end - start);

#ifdef TIMER
            fprintf(stderr, "\033[032mTotal Cost\n%ld (s) : %-9ld (ns)\033[0m\n",
                duration.count() / (long) 1e9, duration.count() % (long) 1e9);
#endif

            // finish            // end of dpptrs()--------------------------------------

reply_dppsv:
            /* prepare reply packets */
#ifdef DEBUGGING
            fprintf(stderr, "server: ready to reply\n");
#endif
            int rpl_packet_payload_size = sizeof(reply_lapack_dppsv_fixed) +
              len_AP * sizeof(double) + len_B * sizeof(double);
            int rpl_packet_total_size = sizeof(struct reply_packet_header) +
              rpl_packet_payload_size;

            pthread_mutex_lock(&(reply_queue->meta_info.mutex));

#ifdef DEBUGGING
            fprintf(stderr, "server: acquire reply queue.\n");
#endif

            while (!can_malloc_straight(reply_queue, rpl_packet_total_size))
            {
              /* wait until there is space in reply queue or timeout */
              struct timespec abs_timeout = {
                .tv_sec = 0, .tv_nsec = 0};
              if (clock_gettime(CLOCK_REALTIME, &abs_timeout) < 0)
                throw "clock_gettime failed.\n";
              calc_timespec_sum(&wait_timeout, &abs_timeout);
              pthread_cond_timedwait(&(reply_queue->meta_info.can_produce),
                  &(reply_queue->meta_info.mutex), &abs_timeout);

              if (terminated == true)
                break;
            }

            if (terminated == true)
            {
              pthread_mutex_unlock(&(reply_queue->meta_info.mutex));
              pthread_mutex_unlock(&(request_queue->meta_info.mutex));
              break; // break switch
            }

#ifdef DEBUGGING
            fprintf(stderr, "server: can enqueue now.\n");
#endif

            uint8_t *raw_reply_packet =
              malloc_straight(reply_queue, rpl_packet_total_size);
            if (raw_reply_packet == nullptr)
              throw "raw_reply_packet == nullptr";

            /* load packet header -> reply queue */
            struct reply_packet_header *reply_header =
              (struct reply_packet_header *)raw_reply_packet;
            reply_header->common_header.preamble = (uint32_t)(PREAMBLE);
            reply_header->common_header.packet_size = rpl_packet_total_size;
            reply_header->client_pthread_id = request_header->client_pthread_id;
            reply_header->transaction_id = request_header->transaction_id;
            reply_header->function_id = request_header->function_id;
            reply_header->payload_len = rpl_packet_payload_size;

            /* load fixed part -> reply queue */
            struct reply_lapack_dppsv_fixed *rpl_pkt_fixed =
              reinterpret_cast<reply_lapack_dppsv_fixed *>(reply_header + 1);
            rpl_pkt_fixed->info = info;

            /* load flexible part -> reply queue */
            double *rpl_pkt_flexible =
              reinterpret_cast<double *>(rpl_pkt_fixed + 1);
            double *rpl_pkt_flexible_pos = rpl_pkt_flexible;

            memcpy(rpl_pkt_flexible_pos, ap, len_AP * sizeof(double));
            rpl_pkt_flexible_pos += len_AP;
            memcpy(rpl_pkt_flexible_pos, b, len_B * sizeof(double));

            /* broadcast to clients and release lock */
            pthread_cond_broadcast(&(reply_queue->meta_info.can_consume));
            pthread_mutex_unlock(&(reply_queue->meta_info.mutex));

            pthread_cond_broadcast(&(request_queue->meta_info.can_produce));
            pthread_mutex_unlock(&(request_queue->meta_info.mutex));

#ifdef DEBUGGING
            fprintf(stderr, "finish 1 DPPSV.\n");
#endif
            /* free resources */
            free(ap);
            free(b);
            break;
          }

        case LAPACK_DGESVD:
          {

#ifdef DEBUGGING
            fprintf(stderr, "received LAPACK_DGESVD from client %u\n",
                request_header->client_pthread_id);
#endif
            struct request_lapack_dgesvd_fixed *req_pkt_fixed =
              reinterpret_cast<struct request_lapack_dgesvd_fixed *>
              (request_header + 1);
            const char jobu         = req_pkt_fixed->jobu;
            const char jobvt        = req_pkt_fixed->jobvt;
            const lapack_int n      = req_pkt_fixed->n;
            const lapack_int m      = (req_pkt_fixed->m) > n ? (req_pkt_fixed->m) : n;
            const lapack_int lda    = m;
            const lapack_int ldu    = m;
            const lapack_int ldvt   = req_pkt_fixed->ldvt;
            const lapack_int lwork  = req_pkt_fixed->lwork;
            const lapack_int len_A  = lda * n;
            const lapack_int len_S  = m <= n ? m : n;
            const lapack_int len_U  = m * n;
            const lapack_int len_VT = n * n;
            const lapack_int len_W  = lwork >= 1 ? lwork : 1;
            lapack_int info = 0;
            double *a     = (double *) malloc (len_A * sizeof(double));
            double *s     = (double *) malloc (len_S * sizeof(double));
            double *u     = (double *) malloc (len_U * sizeof(double));
            double *vt    = (double *) malloc (len_VT * sizeof(double));
            double *work  = (double *) malloc (len_W * sizeof(double));

            double *req_pkt_flexible =
              reinterpret_cast<double *>(req_pkt_fixed + 1);

            memcpy(a, req_pkt_flexible, len_A * sizeof(double));

#ifdef DEBUGGING
            fprintf(stderr, "a[%d]={", len_A);
            for (int i = 0; i < len_A; i++)
              fprintf(stderr, " %lf", a[i]);
            fprintf(stderr, " }\n");
            fprintf(stderr, "s[%d]={", len_S);
            for (int i = 0; i < len_S; i++)
              fprintf(stderr, " %lf", s[i]);
            fprintf(stderr, " }\n");
            fprintf(stderr, "u[%d]={", len_U);
            for (int i = 0; i < len_U; i++)
              fprintf(stderr, " %lf", u[i]);
            fprintf(stderr, " }\n");
            fprintf(stderr, "vt[%d]={", len_VT);
            for (int i = 0; i < len_VT; i++)
              fprintf(stderr, " %lf", vt[i]);
            fprintf(stderr, " }\n");
            fprintf(stderr, "work[%d]={", len_W);
            for (int i = 0; i < len_W; i++)
              fprintf(stderr, " %lf", work[i]);
            fprintf(stderr, " }\n");
#endif

            /* calculation starts from here */
            start = system_clock::now();
            //TODO calculation

            int flag, i, its, j, jj, k, l, nm;
	          double c, f, h, ss, x, y, z;
          	double anorm = 0.0, g = 0.0, scale = 0.0;
           	double *rv1;
            rv1 = (double *) malloc ((unsigned int) n  *sizeof(double));
            /* Householder reduction to bidiagonal form */
	          for (i = 0; i < n; i++)
            {
		        /* left-hand reduction */
				      l = i + 1;
				      rv1[i] = scale * g;
				      g = ss = scale = 0.0;
				      if (i < m)
              {
					      for (k = i; k < m; k++)
						      scale += fabs(a[k * n + i]);
					      if (scale)
                {
					      	for (k = i; k < m; k++)
                  {
					      		a[k * n + i] = a[k * n  + i] / scale;
					      		ss += a[k * n + i] * a[k * n + i];
					      	}
					      	f = a[i * n + i];
					      	g = -SIGN(sqrt(ss), f);
					      	h = f * g - ss;
					      	a[i * n + i] = f - g;

                  if (i != n - 1)
                  {
							      for (j = l; j < n; j++)
                    {
						      		for (ss = 0.0, k = i; k < m; k++)
						      			ss += a[k * n + i] * a[k * n + j];
						      		f = ss / h;
						      		for (k = i; k < m; k++)
						      			a[k * n + j] += f * a[k * n  + i];
						      	}						      
                  }
						      for (k = i; k < m; k++)
						      	a[k * n + i] = a[k * n + i] * scale;
					      }
				      }
				      s[i] = scale * g;

              /* right-hand reduction */
		          g = ss = scale = 0.0;
		          if (i < m && i != n - 1)
		          {
		          	for (k = l; k < n; k++)
		          		scale += fabs(a[i * n + k]);
		          	if (scale)
		          	{
		          		for (k = l; k < n; k++)
		          		{
		          			a[i * n + k] = a[i * n + k] / scale;
		          			ss += a[i * n + k] * a[i * n + k];
		          		}
		          		f = a[i * n + l];
		          		g = -SIGN(sqrt(ss), f);
		          		h = f * g - ss;
			          	a[i * n + l] = f - g;
			          	for (k = l; k < n; k++)
			          		rv1[k] = a[i * n + k] / h;
			          	if (i != m - 1)
		          		{
			          		for (j = l; j < m; j++)
			          		{
			          			for (ss = 0.0, k = l; k < n; k++)
			          				ss += a[j * n + k] * a[i * n + k];
			          			for (k = l; k < n; k++)
			          				a[j * n + k] += (ss * rv1[k]);
			          		}
		          		}
	          			for (k = l; k < n; k++)
		          			a[i * n + k] = a[i * n + k] * scale;
		          	}
		          }
		          anorm = max(anorm, (fabs(s[i]) + fabs(rv1[i]))); 
            }

            /* accumulate the right-hand transformation */
	          for (i = n - 1; i >= 0; i--)
	          {
	          	if (i < n - 1)
	          	{
		          	if (g)
		          	{
		          		for (j = l; j < n; j++)
		          			vt[j * n + i] = (a[i * n + j] / a[i * n + l]) / g;
		          		/* double division to avoid underflow */
		          		for (j = l; j < n; j++)
		          		{
		          			for (ss = 0.0, k = l; k < n; k++)
		          				ss += a[i * n + k] * vt[k * n + j];
		          			for (k = l; k < n; k++)
			          			vt[k * n + j] += ss * vt[k * n + i];
		          		}
	          		}
		          	for (j = l; j < n; j++)
			          	vt[i * n + j] = vt[j * n + i] = 0.0;
		          }
		          vt[i * n + i] = 1.0;
		          g = rv1[i];
	          	l = i;
	          }

	          /* accumulate the left-hand transformation */
	          for (i = n - 1; i >= 0; i--)
	          {
	          	l = i + 1;
		          g = s[i];
		          if (i < n - 1)
		          	for (j = l; j < n; j++)
		          		a[i * n + j] = 0.0;
		          if (g)
		          {
		          	g = 1.0 / g;
		          	if (i != n - 1)
		          	{
		          		for (j = l; j < n; j++)
		          		{
		          			for (ss = 0.0, k = l; k < m; k++)
		          				ss += a[k * n + i] * a[k * n + j];
			          		f = (ss / a[i * n + i]) * g;
			          		for (k = i; k < m; k++)
		          				a[k * n + j] += f * a[k * n + i];
		          		}
	          		}
	          		for (j = i; j < m; j++)
		          		a[j * n + i] = a[j * n + i] * g;
	          	}
		          else
	          	{
	          		for (j = i; j < m; j++)
		          		a[j * n + i] = 0.0;
		          }
		          ++a[i * n + i];
	          }

            /* diagonalize the bidiagonal form */
	          for (k = n - 1; k >= 0; k--)
	          {                             /* loop over singular values */
	          	for (its = 0; its < 30; its++)
	          	{                         /* loop over allowed iterations */
	          		flag = 1;
	          		for (l = k; l >= 0; l--)
	          		{                     /* test for splitting */
		          		nm = l - 1;
		          		if (fabs(rv1[l]) + anorm == anorm)
		          		{
		          			flag = 0;
		          			break;
		          		}
		          		if (fabs(s[nm]) + anorm == anorm)
		          			break;
		          	}
		          	if (flag)
		          	{
		          		c = 0.0;
		          		ss = 1.0;
		      				for (i = l; i <= k; i++)
		      				{
			      				f = ss * rv1[i];
			      				if (fabs(f) + anorm != anorm)
		      					{
		      						g = s[i];
		      						h = PYTHAG(f, g);
		      						s[i] = h;
		      						h = 1.0 / h;
			      					c = g * h;
			      					ss = (-f * h);
			      					for (j = 0; j < m; j++)
			      					{
			      						y = a[j * n + nm];
			      						z = a[j * n + i];
			      						a[j * n + nm] = y * c + z * ss;
			      						a[j * n + i] = z * c - y * ss;
			      					}
		      					}
		      				}
		      			}
		      			z = s[k];
		      			if (l == k)
		      			{                  /* convergence */
		      				if (z < 0.0)
		      				{              /* make singular value nonnegative */
				      			s[k] = -z;
			      				for (j = 0; j < n; j++)
			      					vt[j * n + k] = (-vt[j * n + k]);
			      			}
			      			break;
		      			}
		      			if (its >= 30) {
		      				free((void*)rv1);
		      				fprintf(stderr, "No convergence after 30,000! iterations \n");
	      					return(0);
		      			}
      		
			      		/* shift from bottom 2 x 2 minor */
		      			x = s[l];
		      			nm = k - 1;
		      			y = s[nm];
	      				g = rv1[nm];
	      				h = rv1[k];
	      				f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
	      				g = PYTHAG(f, 1.0);
	      				f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

		      			/* next QR transformation */
		      			c = ss = 1.0;
		      			for (j = l; j <= nm; j++)
		      			{
		      				i = j + 1;
	      					g = rv1[i];
	      					y = s[i];
	      					h = ss * g;
	      					g = c * g;
		      				z = PYTHAG(f, h);
	      					rv1[j] = z;
		      				c = f / z;
		      				ss = h / z;
		      				f = x * c + g * ss;
		      				g = g * c - x * ss;
		      				h = y * ss;
		      				y = y * c;
		      				for (jj = 0; jj < n; jj++)
			      			{
		      					x = vt[jj * n + j];
		      					z = vt[jj * n + i];
		      					vt[jj * n + j] = (x * c + z * ss);
		      					vt[jj * n + i] = (z * c - x * ss);
		      				}
	      					z = PYTHAG(f, h);
		      				s[j] = z;
	      					if (z)
	      					{
	      						z = 1.0 / z;
		      					c = f * z;
			      				ss = h * z;
		      				}
		      				f = (c * g) + (ss * y);
		      				x = (c * y) - (ss * g);
	      					for (jj = 0; jj < m; jj++)
				      		{
			      				y = a[jj * n + j];
			      				z = a[jj * n + i];
			      				a[jj * n + j] = (y * c + z * ss);
		      					a[jj * n + i] = (z * c - y * ss);
		      				}
		      			}
		      			rv1[l] = 0.0;
	      				rv1[k] = f;
		      			s[k] = x;
	      			}
	      		}


fprintf(stderr, "------------------------\n");
for (int i = 0; i < len_A; i++)
  fprintf(stderr, "%lf ", a[i]);
fprintf(stderr, "\n------------------------\n");

fprintf(stderr, "s:\n");
for (int i = 0; i < len_S; i++)
  fprintf(stderr, "%lf ", s[i]);
fprintf(stderr, "\nu:\n");
for (int i = 0; i < len_U; i++)
  u[i] = a[i], fprintf(stderr, "%lf ", u[i]);
fprintf(stderr, "\nvt:\n");
for (int i = 0; i < len_VT; i++)
  fprintf(stderr, "%lf ", vt[i]);
fprintf(stderr, "\n");

            /* calculation ends here */
            end = system_clock::now();
            duration = duration_cast<nanoseconds> (end - start);

#ifdef TIMER
            fprintf(stderr, "\033[032mTotal Cost\n%ld (s) : %-9ld (ns)\033[0m\n",
                duration.count() / (long) 1e9, duration.count() % (long) 1e9);
#endif

            // finish            // end of dpptrs()--------------------------------------

reply_dgesvd:
            /* prepare reply packets */
#ifdef DEBUGGING
            fprintf(stderr, "server: ready to reply\n");
#endif
            int rpl_packet_payload_size = sizeof(reply_lapack_dppsv_fixed) +
              (len_A + len_S + len_U + len_VT + len_W) * sizeof(double);
            int rpl_packet_total_size = sizeof(struct reply_packet_header) +
              rpl_packet_payload_size;

            pthread_mutex_lock(&(reply_queue->meta_info.mutex));
            while (!can_malloc_straight(reply_queue, rpl_packet_total_size))
            {
              /* wait until there is space in reply queue or timeout */
              struct timespec abs_timeout = {
                .tv_sec = 0, .tv_nsec = 0};
              if (clock_gettime(CLOCK_REALTIME, &abs_timeout) < 0)
                throw "clock_gettime failed.\n";
              calc_timespec_sum(&wait_timeout, &abs_timeout);
              pthread_cond_timedwait(&(reply_queue->meta_info.can_produce),
                  &(reply_queue->meta_info.mutex), &abs_timeout);
              if (terminated == true)
                break;
            }

            if (terminated == true)
            {
              pthread_mutex_unlock(&(reply_queue->meta_info.mutex));
              pthread_mutex_unlock(&(request_queue->meta_info.mutex));
              break; // break switch
            }

#ifdef DEBUGGING
            fprintf(stderr, "server: can enqueue now.\n");
#endif

            uint8_t *raw_reply_packet =
              malloc_straight(reply_queue, rpl_packet_total_size);
            if (raw_reply_packet == nullptr)
              throw "raw_reply_packet == nullptr";

            /* load packet header -> reply queue */
            struct reply_packet_header *reply_header =
              (struct reply_packet_header *)raw_reply_packet;
            reply_header->common_header.preamble = (uint32_t)(PREAMBLE);
            reply_header->common_header.packet_size = rpl_packet_total_size;
            reply_header->client_pthread_id = request_header->client_pthread_id;
            reply_header->transaction_id = request_header->transaction_id;
            reply_header->function_id = request_header->function_id;
            reply_header->payload_len = rpl_packet_payload_size;

            /* load fixed part -> reply queue */
            struct reply_lapack_dppsv_fixed *rpl_pkt_fixed =
              reinterpret_cast<reply_lapack_dppsv_fixed *>(reply_header + 1);
            rpl_pkt_fixed->info = info;

            /* load flexible part -> reply queue */
            double *rpl_pkt_flexible =
              reinterpret_cast<double *>(rpl_pkt_fixed + 1);
            double *rpl_pkt_flexible_pos = rpl_pkt_flexible;

            memcpy(rpl_pkt_flexible_pos, a, len_A * sizeof(double));
            rpl_pkt_flexible_pos += len_A;
            memcpy(rpl_pkt_flexible_pos, s, len_S * sizeof(double));
            rpl_pkt_flexible_pos += len_S;
            memcpy(rpl_pkt_flexible_pos, u, len_U * sizeof(double));
            rpl_pkt_flexible_pos += len_U;
            memcpy(rpl_pkt_flexible_pos, vt, len_VT * sizeof(double));
            rpl_pkt_flexible_pos += len_VT;
            memcpy(rpl_pkt_flexible_pos, work, len_W * sizeof(double));

            /* broadcast to clients and release lock */
            pthread_cond_broadcast(&(reply_queue->meta_info.can_consume));
            pthread_mutex_unlock(&(reply_queue->meta_info.mutex));

            pthread_cond_broadcast(&(request_queue->meta_info.can_produce));
            pthread_mutex_unlock(&(request_queue->meta_info.mutex));

#ifdef DEBUGGING
            fprintf(stderr, "finish 1 DGESVD.\n");
#endif
            /* free resources */
            free(a);
            free(s);
            free(u);
            free(vt);
            free(work);
            break;
          }

      } // end of switch

      if (terminated == true)
      {
        pthread_mutex_unlock(&(reply_queue->meta_info.mutex));
        pthread_mutex_unlock(&(reply_queue->meta_info.mutex));
        break; // break while
      }

    } // end of while
  }
  catch (const char *msg)
  {
    fprintf(stderr, "worker thread exception: %s", msg);
    return nullptr;
  }

  return nullptr;
}

int main(int argc, char *argv[], char *envp[])
{

  /* server */
  if (argc != 2)
  {
    fprintf(stderr, "Usage: ./server <server_core_index>\n");
    exit(EXIT_FAILURE);
  }
  signal(SIGINT, sigint_handler);

  try
  {
    /* bind cpu core to server thread */
    int server_core_index = atoi(argv[1]);
    cpu_set_t *cpuset = CPU_ALLOC(NUM_CPUS);
    size_t cpuset_size = CPU_ALLOC_SIZE(NUM_CPUS);
    CPU_ZERO_S(cpuset_size, cpuset);
    int ci = server_core_index == -1 ? 0 : server_core_index;
    CPU_SET_S(ci, cpuset_size, cpuset);
    pthread_attr_t server_thread_attr;
    pthread_attr_init(&server_thread_attr);
    pthread_attr_setaffinity_np(&server_thread_attr, cpuset_size, cpuset);

    char shm_name[32];
    if(server_core_index == -1)
      strcpy(shm_name, SHARED_MEM_NAME_LAPACK);
    else
      sprintf(shm_name, "%s_%d", SHARED_MEM_NAME_LAPACK, server_core_index);
    if ((fd_shm = shm_open(shm_name, O_RDWR | O_CREAT, 0666)) == -1)
      throw "shm_open failed.\n";
    if (ftruncate(fd_shm, SHARED_MEM_SIZE) == -1)
      throw "ftrancate failed.\n";
    if ((p_shm = mmap(NULL, SHARED_MEM_SIZE, PROT_READ | PROT_WRITE,
            MAP_SHARED, fd_shm, 0)) == MAP_FAILED)
      throw "mmap failed.\n";

    request_queue = (struct circular_queue *)p_shm;
    init_circular_queue(request_queue);
    reply_queue = (struct circular_queue *)
      (((uint8_t *)p_shm) + sizeof(struct circular_queue));
    init_circular_queue(reply_queue);

    struct soccs_sce_server_pthread_argument pthread_arg = {
      .p_c2s_queue = request_queue,
      .p_s2c_queue = reply_queue
    };

#ifdef DEBUGGING
    fprintf(stderr, "Shared memory '%s'of %u bytes "
        "created and mmaped at virtual address %p.\n",
        SHARED_MEM_NAME_LAPACK, SHARED_MEM_SIZE, p_shm);
    fprintf(stderr, "c2s_queue of %lu bytes "
        "allocated and initialized at virtual address %p.\n",
        sizeof(struct circular_queue), request_queue);
    fprintf(stderr, "s2c_queue of %lu bytes "
        "allocated and initialized at virtual address %p.\n",
        sizeof(struct circular_queue), reply_queue);
#endif

    /* create server thread */
    pthread_create(&server_pthread, &server_thread_attr, server_soccs,
        (void *)&pthread_arg);

    /* wait for the server thread */
    pthread_join(server_pthread, NULL);
    munmap(p_shm, SHARED_MEM_SIZE);
    close(fd_shm);
    shm_unlink(SHARED_MEM_NAME_LAPACK);
    fprintf(stderr, "Main server thread: exited.\n");
  }
  catch (const char *msg)
  {
    fprintf(stderr, "\n Exception: %s", msg);
    fprintf(stderr, "\n errno: %s\n", strerror(errno));
  }

  fprintf(stderr, "\nOK\n");

  return 0;
}

