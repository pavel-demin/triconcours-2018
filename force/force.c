#include <errno.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

int32_t sla[252] = {352947, 563271111, 789655004, 282722171, 403967413, 210760300, 1047836600, 550097380, 569565325, 276060772, 114984064, 875623171, 966952333, 485362058, 263484080, 265652446, 202157809, 347166309, 111086464, 866911205, 573819363, 906725491, 791368125, 84109180, 578756937, 140669476, 926129608, 498856603, 506777136, 483189153, 250691744, 13225994, 24723693, 1066772783, 982945305, 843786587, 612913148, 825965600, 669553395, 369607888, 403333072, 286809348, 377674548, 688209527, 390605950, 56701174, 567642377, 169328441, 487284015, 353560674, 206996669, 69507743, 1059274456, 586356765, 95701675, 1066550349, 465714134, 753298646, 200510837, 583795316, 1068834692, 203969692, 734721235, 428834798, 471339057, 802107443, 191200458, 870570097, 865541671, 104648122, 31887964, 143849425, 684441276, 408378980, 267781048, 544091447, 559589206, 111161609, 1056131396, 374296697, 824984416, 284912856, 715588209, 982605463, 500771211, 473330684, 989377500, 501763779, 1049293707, 343632788, 852749536, 929341701, 797404076, 598606228, 887738028, 570399063, 330060036, 368773254, 334283139, 479505024, 608561404, 700648190, 67559409, 525879623, 489882936, 10202754, 752736541, 416888513, 479848057, 1005199514, 134380849, 459874322, 314092704, 441280205, 261638924, 392628435, 762610041, 1004553791, 19132623, 513189534, 882179989, 571984638, 135273474, 429837168, 138302471, 872332446, 420570256, 81876895, 641698361, 361486496, 272309312, 414963430, 337232781, 662005834, 219276931, 305441065, 1062063125, 210868011, 710643327, 552094098, 842416626, 136548511, 384547557, 238763256, 320849172, 180604590, 1026959322, 779262318, 632765800, 555780656, 525362765, 392985118, 314921877, 418547111, 438608828, 460962272, 345648951, 378662725, 116639594, 782828470, 439540348, 10883156, 377092897, 576077582, 205898174, 934463504, 980209620, 1036035606, 853028569, 250340214, 547520628, 211775804, 935543152, 872242557, 1057281077, 369624037, 674749315, 724339221, 958242731, 131969240, 730159556, 1070098734, 1047620536, 139935202, 396354952, 43413673, 582912358, 187603292, 542534848, 167625174, 851496377, 268593034, 232504847, 362476216, 800397596, 444839567, 1028035765, 617428820, 485195436, 684285193, 1006309066, 529010313, 494040414, 91725376, 808877709, 162436233, 614052694, 650962399, 369608755, 417896337, 238476849, 875915841, 517139502, 697290705, 546617256, 61179726, 676729792, 724228035, 163281452, 861013125, 229036565, 53119111, 493443248, 792574649, 1034805225, 575109749, 45656296, 693704893, 409417088, 535393028, 394140666, 408864290, 908201857, 908589546, 982036032, 594020805, 56203135, 787035492, 279990347, 661099855, 24881810, 503019732, 683259308, 944137450, 361462870, 948968954, 1033905432, 484674327, 515940818, 952103363, 26834724, 39640398};

int32_t st[5];
int32_t ll[252];
int32_t ls[30];

void pot(int32_t k, int32_t v)
{
  int32_t l, ol;
  st[0] += v - ls[k];
  for(l = 0; l < 252; ++l)
  {
    if(sla[l] & (1 << k))
    {
      ol = ll[l];
      ll[l] += v - ls[k];
      st[1] += ((93 <= ll[l] && ll[l] <= 186) && ((ol < 93) || (ol > 186))) - (((ll[l] < 93) || (ll[l] > 186)) && (93 <= ol && ol <= 186));
      st[2] += (ll[l] > 93 ? ll[l] - 93 : 0) - (ol > 93 ? ol - 93 : 0);
      st[3] += ll[l]*(ll[l] < 93) - ol*(ol < 93);
      st[4] += ((ll[l] > 186) && (ol <= 186)) - ((ll[l] <= 186) && (ol > 186));
    }
  }
  ls[k] = v;
}

double score()
{
  return st[1] - (st[0] + st[4] / 2.0 + st[2] / 5.0 + st[3] / 10.0) / 93.0;
}

struct STATE
{
  int8_t code[30];
  double score;
};

void init(struct STATE *state)
{
  int32_t k, v;

  for(k = 0; k < 30; ++k)
  {
    v = 8 + k % 2;
    state->code[k] = v;
    pot(k, v);
  }

  state->score = score();
}

void random_walk(struct STATE *state)
{
  int32_t k, v;

  for(k = 0; k < 30; ++k)
  {
    v = state->code[k] - 1 + rand() % 3;
    if(v < 0) v = 0;
    if(v > 93) v = 93;
    state->code[k] = v;
    pot(k, v);
  }

  state->score = score();
}

int32_t local_search(struct STATE *state)
{
  int32_t i, k, l, v, kmax, vmax;
  double s, smax;

  smax = state->score;

  for(i = 0; i < 10; ++i)
  {
    kmax = -1;
    vmax = -1;
    for(k = 0; k < 30; ++k)
    {
      for(l = -3; l < 4; ++l)
      {
        v = state->code[k] + l;
        if(l == 0 || v < 0 || v > 93) continue;
        pot(k, v);
        s = score();
        if(smax < s)
        {
          smax = s;
          kmax = k;
          vmax = v;
        }
      }
      pot(k, state->code[k]);
    }
    if(kmax < 0) break;
    state->code[kmax] = vmax;
    pot(kmax, vmax);
    state->score = smax;
  }

  return kmax < 0;
}

void print(struct STATE *state)
{
  int32_t k;

  printf("%.3f ", state->score);
  for(k = 0; k < 30; ++k)
  {
    printf("%c", state->code[k] + 33);
  }
  printf("\n");
}

int main(int argc, char *argv[])
{
  struct STATE current, next;
  double delta, result, t;
  unsigned seed;

  seed = 1537611153;
  srand(seed);
  printf("seed: %d\n", seed);

  memset(st, 0, 5*4);
  memset(ll, 0, 252*4);
  memset(ls, 0, 30*4);

  init(&current);

  t = 1.0;
  result = 0.0;
  while(t > 0.0)
  {
    t *= 0.999999;

    next = current;

    random_walk(&next);

    if(!local_search(&next)) continue;

    if(result < next.score)
    {
      result = next.score;
      print(&next);
    }

    delta = next.score - current.score;
    if(delta >= 0 || rand() < RAND_MAX * exp(delta / t))
    {
      current = next;
    }
  }

  return EXIT_SUCCESS;
}
