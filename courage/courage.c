#include <errno.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <complex.h>
#include <math.h>
#include <time.h>

#ifndef M_PI
#define M_PI (3.14159265358979323846)
#endif

#define MIN(a,b) (((a)<(b))?(a):(b))

struct STATE
{
  int32_t size;
  double code[100];
  double smax, s;
  double x, y;
  double vmin1, vmin2, dmin2, smin2;
};

const int32_t w = 65537;

int32_t mrand(int32_t x)
{
  return (75 * x) % w;
}

void init_coordinates(double *x, double *y)
{
  int32_t a = 31;
  int32_t i;
  double complex db[8];

  for(i = 0; i < 8; ++i)
  {
    a = mrand(a);
    db[i] = (8.0 * a / w - 4.0) / 3.0;
    a = mrand(a);
    db[i] += I * (2.0 * a / w - 1.0);
    x[i] = creal(db[i]);
    y[i] = cimag(db[i]);
  }
}

void score(struct STATE *state)
{
  double b = 0.085;
  int32_t a = 31;
  double s = 0.0, smax = 0.0;
  double complex current = CMPLX(0.0, 0.0), z = CMPLX(0.0, 0.0);
  double complex fragment, previous, t, db[8];
  double v = 0.0;
  int32_t counter = 0;
  int32_t i, k;

  for(i = 0; i < 8; ++i)
  {
    a = mrand(a);
    db[i] = (8.0 * a / w - 4.0) / 3.0;
    a = mrand(a);
    db[i] += I * (2.0 * a / w - 1.0);
    db[i] = cabs(db[i]) + I * carg(db[i]) * M_PI / carg(-1);
  }

  for(k = 0; k < state->size; ++k)
  {
    v = state->code[k];
    t = z;
    a = mrand(a);
    z += v + I * (0.04 + (double) a / w) / (9.0 * fabs(creal(z)) + 1.0);
    current = creal(z) * cexp(I * cimag(z));
    previous = creal(t) * cexp(I * cimag(t));
    s -= 500.0 * cabs(current - previous);
    if(counter == 8 && creal(z) * creal(t) <= 0.0)
    {
      ++counter;
      s += 500.0 * (3.0 - MIN(1.0, cabs(current)) - MIN(1.0, fabs(creal(z))));
    }
    for(i = 0; i < 8; ++i)
    {
      fragment = creal(db[i]) * cexp(I * cimag(db[i]));
      if(M_PI >= cimag(db[i]) && b >= cabs(current - fragment))
      {
        db[i] = creal(db[i]) + I * (2.0 * M_PI + cimag(db[i]));
        ++counter;
        s += 500.0 * (3.0 - MIN(1.0, cabs(current - fragment)) - MIN(1.0, fabs(creal(z) - creal(db[i]))));
      }
    }
    if(smax < s)
    {
      smax = s;
    }
  }
  state->smax = smax;
  state->s = s;
  state->x = creal(current);
  state->y = cimag(current);
}

void random_walk(struct STATE *state, double step)
{
  int32_t i, k;
  double v;
  for(i = 0; i < 3; ++i)
  {
    k = rand() % state->size;
    v = (rand() % 3 - 1) * step;
    state->code[k] += v;
  }

  score(state);
}

void local_search(struct STATE *state, double start, double stop, double step)
{
  int32_t i, k, kmax;
  double backup, cmax, s, smax, v;

  smax = state->smax;

  for(i = 0; i < 10; ++i)
  {
    kmax = -1;
    cmax = -1;
    for(k = 0; k < state->size; ++k)
    {
      backup = state->code[k];
      for(v = start; v < stop; v += step)
      {
        state->code[k] = backup + v;
        score(state);
        s = state->smax;
        if(smax < s)
        {
          smax = s;
          kmax = k;
          cmax = backup + v;
        }
      }
      state->code[k] = backup;
    }
    if(kmax < 0) break;
    state->code[kmax] = cmax;
    score(state);
  }
}

double round_value(double value, int32_t digits)
{
  return floor(pow(10, digits) * value + 0.5) / pow(10, digits);
}

void round_code(struct STATE *state, int32_t digits)
{
  int32_t k;

  for(k = 0; k < state->size; ++k)
  {
    state->code[k] = round_value(state->code[k], digits);
  }

  score(state);
}

void make_move(struct STATE *state, double v)
{
  int32_t k;

  k = state->size;
  state->code[k] = v;
  state->size = k + 1;

  score(state);
}

void analyze_moves(struct STATE *state, double x1, double y1, double x2, double y2, double start, double stop, double step, int32_t digits)
{
  double a1, a2, x3, y3, d, v, s;
  double dmin1 = 10.0, dmin2 = 10.0, vmin1 = 0.0, vmin2 = 0.0, smin2 = 0.0;

  a1 = atan2(y2 - y1, x2 - x1);

  for(v = start; v < stop; v += step)
  {
    v = round_value(v, digits);
    make_move(state, v);
    --state->size;
    s = state->s;
    x3 = state->x;
    y3 = state->y;
    a2 = atan2(round_value(y3 - y1, digits), round_value(x3 - x1, digits));
    d = fabs(a1 - a2);
    if(dmin1 > d)
    {
      dmin1 = d;
      vmin1 = v;
    }
    d = fabs(hypot(x3 - x2, y3 - y2));
    if(dmin2 > d)
    {
      dmin2 = d;
      vmin2 = v;
      smin2 = s;
    }
  }
  state->vmin1 = vmin1;
  state->vmin2 = vmin2;
  state->dmin2 = dmin2;
  state->smin2 = smin2;
}

void optimize(struct STATE *state, double start, double stop, double step, int32_t size)
{
  int32_t i;
  double smax = 0.0;
  struct STATE next;

  for(i = 0; i < 20000; ++i)
  {
    next = *state;
    random_walk(&next, step);
    local_search(&next, start, stop, step);
    round_code(&next, size);
    if(smax < next.smax)
    {
      smax = next.smax;
    }
    if(next.smax > state->smax)
    {
      *state = next;
    }
  }
}

void print_state(struct STATE *state, int32_t digits)
{
  int32_t k;
  double v;

  printf("%.4f %d ", state->smax, state->size);
  for(k = 0; k < state->size; ++k)
  {
    v = state->code[k];
    printf("%g,", round_value(v, digits));
  }
  printf("\n");
}

int main(int argc, char *argv[])
{
  int32_t i, size;
  unsigned seed;
  double x[8], y[8];
  double x1, y1, x2, y2, s1, s2;
  struct STATE current, next;
  const int32_t route[9] = {1, 7, 2, 0, 3, 6, 5, 4, 8};

  seed = time(NULL);
  printf("seed: %d\n", seed);

  init_coordinates(x, y);

  memset(&current, 0, sizeof(struct STATE));
  current.size = 9;
  current.dmin2 = 10.0;

  // departure coordinates
  x1 = current.x;
  y1 = current.y;

  // loop over the first seven destinations
  for(i = 0; i < 7; ++i)
  {
    // destination coordinates
    x2 = x[route[i]];
    y2 = y[route[i]];
    // fly to the destination
    for(;;)
    {
      next = current;
      analyze_moves(&next, x1, y1, x2, y2, -1.000, 1.001, 0.001, 3);
      if(current.dmin2 <= 0.085)
      {
        // if both current and next moves reach the destination,
        // then choose the one with the highest score
        if(next.dmin2 <= 0.085 && current.smin2 < next.smin2)
        {
          make_move(&next, next.vmin2);
        }
        else
        {
          --next.size;
          make_move(&next, current.vmin2);
        }

        print_state(&next, 3);
        optimize(&next, -0.002, 0.003, 0.001, 3);
        print_state(&next, 3);

        next.dmin2 = 10.0;

        // next departure coordinates
        x1 = next.x;
        y1 = next.y;

        current = next;
        break;
      }
      else
      {
        make_move(&next, next.vmin1);
        current = next;
      }
    }
  }

  // last two destinations
  make_move(&next, -0.4);
  make_move(&next, -0.4);

  optimize(&next, -0.002, 0.003, 0.001, 3);

  print_state(&next, 3);

  current = next;
  size = current.size;
  s1 = 0.0;
  for(i = 0; i <= size; ++i)
  {
    current.size = i;
    score(&current);
    s2 = current.s;
    if(s2 > s1 && i != 82)
    {
      next.size = i;
      score(&next);
      print_state(&next, 4);
      if(i == 83)
      {
        next.code[80] = 0.05;
        optimize(&next, -0.002, 0.003, 0.001, 4);
      }
      optimize(&next, -0.0002, 0.0003, 0.0001, 4);
      print_state(&next, 4);
    }
    s1 = s2;
  }

  return EXIT_SUCCESS;
}
