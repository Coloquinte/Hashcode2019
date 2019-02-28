#!/usr/bin/python3

import sys
import networkx
import localsolver

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
            print(s[0], s[1], file=f)

def slide_tags(tags, slide):
    if len(slide) == 1:
      return tags[slide[0]]
    else:
      return tags[slide[0]] | tags[slide[1]]

def cost_between_tags(ti, tj):
	return min( len(ti&tj), len(ti-tj), len(tj-ti) )

def sol_cost(tags, sol):
    sol_tags = [slide_tags(tags, s) for s in sol]
    cost = 0
    for i in range(len(sol) - 1):
      cost += cost_between_tags(sol_tags[i], sol_tags[i+1])
    return cost

def make_paired_solution(tags, horizontal):
    vimgs = [i for i, h in enumerate(horizontal) if not h]
    himgs = [i for i, h in enumerate(horizontal) if h]
    return [(i,) for i in himgs] + [(vimgs[2*j], vimgs[2*j+1]) for j in range(int(len(vimgs)/2))]

def make_paired_solution_matching(tags, horizontal):
    vimgs = [i for i, h in enumerate(horizontal) if not h]
    himgs = [i for i, h in enumerate(horizontal) if h]
    vslides = [(vimgs[2*j], vimgs[2*j+1]) for j in range(int(len(vimgs)/2))]
    sum_common = sum(len(tags[v[0]] & tags[v[1]]) for v in vslides)
    #print("Common", sum_common)
    return [(i,) for i in himgs] + vslides

def optimize(tags, horizontal, sol, start, end):
    subtour = sol[start:end]
    nb_cities = len(subtour)
    sol_tags = [slide_tags(tags, s) for s in subtour]
    distance_weight = [[cost_between_tags(sol_tags[i], sol_tags[j]) for i in range(nb_cities)] for j in range(nb_cities)]
    with localsolver.LocalSolver() as ls:
        model = ls.model
        cities = model.list(end - start)
        model.constraint(model.count(cities) == end-start)
        distance_array = model.array(distance_weight)
        dist_selector = model.function(lambda i: model.at(distance_array, cities[i-1], cities[i]))
        obj = model.sum(model.range(1, nb_cities), dist_selector)
        model.maximize(obj)
        model.close()
        ls.param.time_limit = 10
        for c in range(nb_cities):
          cities.value.add(c)
        ls.solve()
        new_subtour = [subtour[c] for c in cities.value]
        sol[start:end] = new_subtour

if len(sys.argv) < 2:
  print("solver.py input [output]")
  sys.exit(1)

tags, horizontal = parse_file(sys.argv[1])

sol = [(i,) for i in range(len(tags))]
sol = make_paired_solution(tags, horizontal)
#sol = make_paired_solution_matching(tags, horizontal)

print (sol_cost(tags, sol))

subsize = 500
for i in range(10):
    if (i+1) * subsize > len(sol):
      break
    optimize(tags, horizontal, sol, i * subsize, min( (i+1) * subsize, len(sol)))
    print (sol_cost(tags, sol))

if len(sys.argv) >= 3:
  write_solution(sol, sys.argv[2])





