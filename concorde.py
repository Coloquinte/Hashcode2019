########## tsp.py ##########

import localsolver
import sys

import os

def parse_file(filename):
    with open(filename) as f:
        lines = f.readlines()
        nb_images = int(lines[0])
        tags = []
        horizontal = []
        for l in lines[1:]:
          vals = l.split()
          assert len(vals) == int(vals[1]) + 2
          this_tags = set(vals[2:])
          this_horizontal = vals[0] == "H"
          tags.append(this_tags)
          horizontal.append(this_horizontal)
        assert len(tags) == nb_images
        assert len(horizontal) == nb_images
        return (tags, horizontal)

def write_solution(sol, filename):
    with open(filename, "w") as f:
        print(len(sol), file=f)
        for s in sol:
          if len(s) == 1:
            print(s[0], file=f)
          else:
            print(s[0], sol[1], file=f)


tags, horizontal = parse_file(sys.argv[1])

appar = [(i,) for i in range(len(tags))]

appar = appar[:200]

nb_cities = len(appar)

if len(sys.argv) >= 3:
  write_solution(sol, sys.argv[2])

m = max([len(tags[i]) for i in range(nb_cities)])

def get_tags_from_slide(tup):
    if len(tup) == 1:
        return tags[tup[0]]
    if len(i) == 2:
        return tags[tup[0]] | tags[tup[1]]

def cost(i, j):
    ti = get_tags_from_slide(i)
    tj = get_tags_from_slide(j)
    return min( len(ti&tj), len(ti-tj), len(tj-ti) )



os.remove("xx.xx") 

with open('xx.xx','a') as file:
  file.write('NAME:  test\n')
  file.write('TYPE: TSP\n')
  file.write('DIMENSION:  '+str(nb_cities) +'\n')
  file.write('EDGE_WEIGHT_TYPE: EXPLICIT\n')
  file.write('EDGE_WEIGHT_FORMAT: FULL_MATRIX\n')
  file.write('EDGE_WEIGHT_SECTION\n')
  for i in range(nb_cities):
    file.write( ' | '.join([ str(m-cost(appar[i],appar[j])) for j in range(nb_cities) ]) )
    file.write('\n')
  file.write('EOF\n')

os.system('./concorde xx.xx')

tsp = []

with open('xx.sol') as sol:
  next(sol)
  for line in sol:
    line = line[:-2]
    tsp = tsp + [ int(i) for i in line.split(' ') ]

print(tsp)
