#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

#include <sys/time.h>

#define MM 0
#define GapM 1
#define MX 2
#define XX 3
#define MY 4
#define YY 5

//#define DEBUG

/*
	q: read quality
	i: insertion penalty
	d: deletion penalty
	c: gap continuation penalty
*/

typedef struct 
{
	int rslen, haplen, *q, *i, *d, *c;
	char *hap, *rs;
} testcase;

int normalize(char c)
{
	return ((int) (c - 33));
}

int read_testcase(testcase *tc)
{
	char *q, *i, *d, *c, *line = NULL;
	int _q, _i, _d, _c;
	int x, size = 0;
	ssize_t read;

	read = getline(&line, (size_t *) &size, stdin);
	if (read == -1)
		return -1;

	tc->hap = (char *) malloc(size);
	tc->rs = (char *) malloc(size);
	q = (char *) malloc(size);
	i = (char *) malloc(size);
	d = (char *) malloc(size);
	c = (char *) malloc(size);

	if (sscanf(line, "%s %s %s %s %s %s\n", tc->hap, tc->rs, q, i, d, c) != 6)
		return -1;

	tc->haplen = strlen(tc->hap);
	tc->rslen = strlen(tc->rs);
	tc->q = (int *) malloc(sizeof(int) * tc->rslen);
	tc->i = (int *) malloc(sizeof(int) * tc->rslen);
	tc->d = (int *) malloc(sizeof(int) * tc->rslen);
	tc->c = (int *) malloc(sizeof(int) * tc->rslen);

	for (x = 0; x < tc->rslen; x++)
	{
		_q = normalize(q[x]);
		_i = normalize(i[x]);
		_d = normalize(d[x]);
		_c = normalize(c[x]);
		tc->q[x] = (_q < 6) ? 6 : _q;
		tc->i[x] = _i;
		tc->d[x] = _d;
		tc->c[x] = _c;
	}

	free(q);
	free(i);
	free(d);
	free(c);
	free(line);

	return 0;
}

template<class T>
struct Context{};

template<>
struct Context<double>
{
	Context()
	{
		for (int x = 0; x < 128; x++)
			ph2pr[x] = pow(10.0, -((double)x) / 10.0);

		INITIAL_CONSTANT = ldexp(1.0, 1020.0);
		LOG10_INITIAL_CONSTANT = log10(INITIAL_CONSTANT);
		RESULT_THRESHOLD = 0.0;
	}

	double LOG10(double v){ return log10(v); }

	static double _(double n){ return n; }
	static double _(float n){ return ((double) n); }
	double ph2pr[128];
	double INITIAL_CONSTANT;
	double LOG10_INITIAL_CONSTANT;
	double RESULT_THRESHOLD;
};

template<>
struct Context<float>
{
	Context()
	{
		for (int x = 0; x < 128; x++)
			ph2pr[x] = powf(10.f, -((float)x) / 10.f);

		INITIAL_CONSTANT = ldexpf(1.f, 120.f);
		LOG10_INITIAL_CONSTANT = log10f(INITIAL_CONSTANT);
		RESULT_THRESHOLD = ldexpf(1.f, -110.f);
	}

	float LOG10(float v){ return log10f(v); }

	static float _(double n){ return ((float) n); }
	static float _(float n){ return n; }
	float ph2pr[128];
	float INITIAL_CONSTANT;
	float LOG10_INITIAL_CONSTANT;
	float RESULT_THRESHOLD;
};

template<class NUMBER>
NUMBER compute_full_prob(testcase *tc, NUMBER *before_last_log = NULL)
{
	int r, c;
	int ROWS = tc->rslen + 1;
	int COLS = tc->haplen + 1;

	Context<NUMBER> ctx;

	NUMBER M[ROWS][COLS];
	NUMBER X[ROWS][COLS];
	NUMBER Y[ROWS][COLS];
	NUMBER p[ROWS][6];

	p[0][MM] = ctx._(0.0);
	p[0][GapM] = ctx._(0.0);
	p[0][MX] = ctx._(0.0);
	p[0][XX] = ctx._(0.0);
	p[0][MY] = ctx._(0.0);
	p[0][YY] = ctx._(0.0);
	for (r = 1; r < ROWS; r++)
	{
		int _i = tc->i[r-1] & 127;
		int _d = tc->d[r-1] & 127;
		int _c = tc->c[r-1] & 127;
		p[r][MM] = ctx._(1.0) - ctx.ph2pr[(_i + _d) & 127];
		p[r][GapM] = ctx._(1.0) - ctx.ph2pr[_c];
		p[r][MX] = ctx.ph2pr[_i];
		p[r][XX] = ctx.ph2pr[_c];
		p[r][MY] = (r == ROWS - 1) ? ctx._(1.0) : ctx.ph2pr[_d];
		p[r][YY] = (r == ROWS - 1) ? ctx._(1.0) : ctx.ph2pr[_c];
	}

	for (c = 0; c < COLS; c++)
	{
		M[0][c] = ctx._(0.0);
		X[0][c] = ctx._(0.0);
		Y[0][c] = ctx.INITIAL_CONSTANT / (tc->haplen);
	}

	for (r = 1; r < ROWS; r++)
	{
		M[r][0] = ctx._(0.0);
		X[r][0] = X[r-1][0] * p[r][XX];
		Y[r][0] = ctx._(0.0);
	}

	 NUMBER result = ctx._(0.0);

	for (r = 1; r < ROWS; r++)
		for (c = 1; c < COLS; c++)
		{
			char _rs = tc->rs[r-1];
			char _hap = tc->hap[c-1];
			int _q = tc->q[r-1] & 127;
			NUMBER distm = ctx.ph2pr[_q];
			if (_rs == _hap || _rs == 'N' || _hap == 'N')
				distm = ctx._(1.0) - distm;
			M[r][c] = distm * (M[r-1][c-1] * p[r][MM] + X[r-1][c-1] * p[r][GapM] + Y[r-1][c-1] * p[r][GapM]);
			X[r][c] = M[r-1][c] * p[r][MX] + X[r-1][c] * p[r][XX];
			Y[r][c] = M[r][c-1] * p[r][MY] + Y[r][c-1] * p[r][YY];
		}

	for (c = 0; c < COLS; c++)
	{
		result += M[ROWS-1][c] + X[ROWS-1][c];
	}

	if (before_last_log != NULL)
		*before_last_log = result;

	return result; //ctx.LOG10(result) - ctx.LOG10_INITIAL_CONSTANT;
}

#define BATCH_SIZE  10000
#define RUN_HYBRID

int main()
{
        testcase tc[BATCH_SIZE];
        float result[BATCH_SIZE], result_avxf;
        double result_avxd;
        struct timeval start, end;
        long long aggregateTimeRead = 0L;
        long long aggregateTimeCompute = 0L;
        long long aggregateTimeWrite = 0L;

        bool noMoreData = false;
        int count =0;
        while (!noMoreData)
        {
                int read_count = BATCH_SIZE;
                gettimeofday(&start, NULL);
                for (int b=0;b<BATCH_SIZE;b++)
                        if (read_testcase(&tc[b])==-1)
                        {
                                read_count = b;
                                noMoreData = true;
                                break;
                        }
                gettimeofday(&end, NULL);
                aggregateTimeRead += ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));

                gettimeofday(&start, NULL);
                for (int b=0;b<read_count;b++)
                {
                        result_avxf = compute_full_prob<float>(&tc[b]);

                        #ifdef RUN_HYBRID
                                #define MIN_ACCEPTED 1e-28f
                                if (result_avxf < MIN_ACCEPTED) {
                                      count++;
                                      result_avxd = compute_full_prob<double>(&tc[b]);
                                      result[b] = log10(result_avxd) - log10(ldexp(1.0, 1020.f));
                                }
                                else
                                      result[b] = log10f(result_avxf) - log10f(ldexpf(1.f, 120.f));
                        #endif

                        #ifndef RUN_HYBRID
                                result[b] = log10f(result_avxf) - log10f(ldexpf(1.f, 120.f));
                        #endif

                }
                gettimeofday(&end, NULL);
                aggregateTimeCompute += ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));

                gettimeofday(&start, NULL);
                for (int b=0;b<read_count;b++)
                        printf("%E\n", result[b]);
                gettimeofday(&end, NULL);
                aggregateTimeWrite += ((end.tv_sec * 1000000 + end.tv_usec) - (start.tv_sec * 1000000 + start.tv_usec));

        }

        printf("AVX Read Time: %ld\n", aggregateTimeRead);
        printf("AVX Compute Time: %ld\n", aggregateTimeCompute);
        printf("AVX Write Time: %ld\n", aggregateTimeWrite);
        printf("AVX Total Time: %ld\n", aggregateTimeRead + aggregateTimeCompute + aggregateTimeWrite);
        printf("# Double called: %d\n", count);

        return 0;
}

