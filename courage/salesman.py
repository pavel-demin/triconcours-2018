import math
import cmath
import random, copy
import numpy as np
import matplotlib.pyplot as plt

w = 65537
p = 75

def mrand(x):
  return (p * x) % w

fragments = []
a = 31
db = [0j] * 8
for i in range(8):
  a = mrand(a)
  db[i] = (8.0 * a / w - 4.0) / 3.0
  a = mrand(a)
  db[i] += 1j * (2.0 * a / w - 1.0)
  x = db[i].real
  y = db[i].imag
  fragments.append((x, y))
  db[i] = abs(db[i]) + 1j * cmath.phase(db[i])

fragments.append((0.0, 0.0))

# the following code is copied from
# https://ericphanson.com/blog/2016/the-traveling-salesman-and-10-lines-of-python

cities = copy.copy(fragments)

xmin = min(pair[0] for pair in cities)
xmax = max(pair[0] for pair in cities)

ymin = min(pair[1] for pair in cities)
ymax = max(pair[1] for pair in cities)

def transform(pair):
  x = pair[0]
  y = pair[1]
  return [(x-xmin)*50, (y-ymin)*50]

rescaled_cities = [transform(b) for b in cities]

cities = rescaled_cities

route = random.sample(range(9), 9)

for temperature in np.logspace(5, 0, num = 100000):
  [i, j] = sorted(random.sample(range(9), 2))
  newroute = route[:i] + route[j:j+1] +  route[i+1:j] + route[i:i+1] + route[j+1:]
  oldDistance =  sum([math.sqrt(sum([(cities[route[(k+1) % 9]][d] - cities[route[k % 9]][d])**2 for d in [0, 1] ])) for k in [j, j-1, i, i-1]])
  newDistance = sum([ math.sqrt(sum([(cities[newroute[(k+1) % 9]][d] - cities[newroute[k % 9]][d])**2 for d in [0, 1] ])) for k in [j, j-1, i, i-1]])
  if math.exp((oldDistance - newDistance) / temperature) > random.random():
    distance = sum([math.sqrt(sum([(fragments[newroute[(k+1) % 9]][d] - fragments[newroute[k % 9]][d])**2 for d in [0, 1] ])) for k in [j, j-1, i, i-1]])
    route = copy.copy(newroute)

distance = sum([math.sqrt(sum([(fragments[route[(k+1)%9]][d] - fragments[route[k%9]][d])**2 for d in [0, 1] ])) for k in range(9)])

print('distance:', distance)
print('route:', route)

fig, ax = plt.subplots()
x = [fragments[route[i % 9]][0] for i in range(10)]
y = [fragments[route[i % 9]][1] for i in range(10)]
ax.plot(x, y, 'ob-', lw = 2, ms = 20, mfc = 'w', mec = 'b', mew = 2)
ax.set_xlim(-1.5, 1.0)
ax.set_ylim(-0.9, 1.2)

for i in range(9):
  x = fragments[i][0]
  y = fragments[i][1]
  ax.annotate(str(i), xy = (x, y - 0.01), va = 'center', ha='center', fontsize = 14)

plt.show()
