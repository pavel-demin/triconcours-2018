mrandmax=2147483647
mrand=0
def mfloor(f):
  return round(f)-(round(f)>f)
def mceil(f):
  return round(f)+(round(f)<f)
def mseed(s):
  global mrand
  mrand=s%mrandmax
def mrandom():
  mseed(mrand*16807)
  return mrand/mrandmax
def muniform(mini,maxi):
  return mrandom()*(maxi-mini)+mini
def mrandint(mini,maxi):
  return round(muniform(mceil(mini),mfloor(maxi)))
mfmax=93
def mf2f(n):
  return float(n/mfmax)
def smf(l):
  s=""
  for k in range(len(l)):
    s=s+chr(33+l[k])
  return s

mseed(42)
st = [0,0,0,0,0]
sla = [0] * 252
ll = [0] * 252
ls = [0] * 30

for k in range(252):
  sla[k] = mrandint(0, 2**30 - 1)

def pot(k, v=93):
  st[0] += v - ls[k]
  for l in range(252):
    if sla[l] & (1 << k):
      ol = ll[l]
      ll[l] += v - ls[k]
      st[1] += ((93 <= ll[l] <= 186) and ((ol < 93) or (ol > 186))) - (((ll[l] < 93) or (ll[l] > 186)) and (93 <= ol <= 186))
      st[2] += max(ll[l] - 93, 0) - max(ol - 93, 0)
      st[3] += ll[l] * (ll[l] < 93) - ol * (ol < 93)
      st[4] += ((ll[l] > 186) and (ol <= 186)) - ((ll[l] <= 186) and (ol > 186))
  ls[k] = v

def score():
  return float(st[1] - mf2f(st[0] + st[4] / 2 + st[2] / 5 + st[3] / 10))

for j in range(30):
  pot(j, 8 + j % 2)

print(score())

while True:
  smax = score()
  kmax = -1
  vmax = -1
  for k in range(30):
    backup = ls[k]
    for v in range(94):
      pot(k, v)
      s = score()
      if smax < s:
        smax = s
        kmax = k
        vmax = v
    pot(k, backup)
  if kmax >= 0:
    pot(kmax, vmax)
    print(smax)
  else:
    for k in range(30):
      v = ls[k] + mrandint(-1, 1)
      if v < 0: v = 0
      if v > 93: v = 93
      pot(k, v)
