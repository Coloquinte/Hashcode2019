#!/usr/bin/python3

import sys
import localsolver
import random

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

def parse_solution(filename):
    with open(filename) as f:
        lines = f.readlines()
        nb_slides = int(lines[0])
        slides = []
        for l in lines[1:]:
          slides.append(tuple([int(i) for i in l.split()]))
        assert len(slides) == nb_slides
        return slides

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

def slide_dist(tags, slide1, slide2):
    t1 = slide_tags(tags, slide1)
    t2 = slide_tags(tags, slide2)
    return cost_between_tags(t1, t2)

def sol_cost(tags, sol):
    sol_tags = [slide_tags(tags, s) for s in sol]
    cost = 0
    for i in range(len(sol) - 1):
      loc = cost_between_tags(sol_tags[i], sol_tags[i+1])
      cost += loc
    return cost

def make_paired_solution(tags, horizontal):
    vimgs = [i for i, h in enumerate(horizontal) if not h]
    himgs = [i for i, h in enumerate(horizontal) if h]
    return [(i,) for i in himgs] + [(vimgs[2*j], vimgs[2*j+1]) for j in range(int(len(vimgs)/2))]

def optimize_tsp(tags, horizontal, sol, start, end):
    print("TSP from ", start, " to ", end, " out of ", len(sol))
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
        ls.param.verbosity = 0
        model.close()
        ls.param.time_limit = 15
        for c in range(nb_cities):
          cities.value.add(c)
        ls.solve()
        new_subtour = [subtour[c] for c in cities.value]
        sol[start:end] = new_subtour
    print (sol_cost(tags, sol))

def shuffle_vertical_images(sol):
    for i in range(len(sol)):
      if len(sol[i]) == 2 and random.random() < 0.5:
        sol[i] = (sol[i][1], sol[i][0])

def optimize_alloc(tags, horizontal, sol, start, end):
    subtour = sol[start:end]
    last_vertical = -2
    verticals = []
    for i, s in enumerate(subtour):
        if len(s) == 2 and i > last_vertical + 1 and random.random() < 0.8:
          last_vertical = i
          verticals.append(i)
    if len(verticals) < 2:
      return
    print("Alloc from ", start, " to ", end, " out of ", len(sol))
    shuffle_vertical_images(subtour)

    nb_cities = len(verticals)
    placement_weight = []
    for v in verticals:
      weights = []
      orig = subtour[v][1]
      for pl in verticals:
        other = subtour[pl][0]
        res = (other, orig)
        cost = 0
        if pl >= 1:
            cost += slide_dist(tags, res, subtour[pl-1])
        if pl + 1 < len(subtour):
            cost += slide_dist(tags, res, subtour[pl+1])
        weights.append(cost)
      placement_weight.append(weights)
        
    with localsolver.LocalSolver() as ls:
        model = ls.model
        cities = model.list(len(verticals))
        model.constraint(model.count(cities) == len(verticals))
        placement_array = model.array(placement_weight)
        dist_selector = model.function(lambda i: model.at(placement_array, i, cities[i]))
        obj = model.sum(model.range(0, nb_cities), dist_selector)
        model.maximize(obj)
        ls.param.verbosity = 0
        model.close()
        ls.param.time_limit = 15
        for c in range(nb_cities):
          cities.value.add(c)
        ls.solve()
        new_subtour = list(subtour)
        allocs = [verticals[c] for c in cities.value]
        for v, alloc in zip(verticals, allocs):
          new_subtour[alloc] = (subtour[alloc][0], subtour[v][1])
        sol[start:end] = new_subtour
    print (sol_cost(tags, sol))

def save_solution(sol):
  if len(sys.argv) >= 3:
    write_solution(sol, sys.argv[2])

def reoptimize_alloc(tags, horizontal, sol, subsize, start):
    for i in range(200):
        if i * subsize + start >= len(sol):
          break
        optimize_alloc(tags, horizontal, sol, i * subsize + start, min( (i+1) * subsize + start, len(sol)))
        save_solution(sol)

def reoptimize_tsp(tags, horizontal, sol, subsize, start):
    for i in range(200):
        if i * subsize + start >= len(sol):
          break
        optimize_tsp(tags, horizontal, sol, i * subsize + start, min( (i+1) * subsize + start, len(sol)))
        save_solution(sol)

def reoptimize(tags, horizontal, sol, subsize, start):
    lst = list(range(200))
    random.shuffle(lst)
    for i in lst:
        if i * subsize + start >= len(sol):
          continue
        optimize_alloc(tags, horizontal, sol, i * subsize + start, min( (i+1) * subsize + start, len(sol)))
        optimize_tsp(tags, horizontal, sol, i * subsize + start, min( (i+1) * subsize + start, len(sol)))
        save_solution(sol)

def shuffle_tsp(sol):
  rng = 1000
  n = len(sol)
  solc = [ sol[rng*s:rng*(s+1)] for s in range(int(len(sol)/rng)) ]
  random.shuffle(solc)
  for s in range(int(len(sol)/rng)):
    sol[rng*s:rng*(s+1)] = solc[s] 

if len(sys.argv) < 2:
  print("solver.py input [output]")
  sys.exit(1)

tags, horizontal = parse_file(sys.argv[1])

sol = [(i,) for i in range(len(tags))]
sol = make_paired_solution(tags, horizontal)
if len(sys.argv) >= 3:
  sol = parse_solution(sys.argv[2])

print ("Initial solution at ", sol_cost(tags, sol))
shuffle_tsp(sol)
print ("Shuffled solution at ", sol_cost(tags, sol))


#reoptimize_alloc(tags, horizontal, sol, 1000, 0)
#reoptimize_alloc(tags, horizontal, sol, 1000, 500)
#reoptimize_tsp(tags, horizontal, sol, 300, 0)
#reoptimize_tsp(tags, horizontal, sol, 300, 150)

reoptimize(tags, horizontal, sol, 1000, 0)
reoptimize(tags, horizontal, sol, 1000, 500)
print ("Final solution at ", sol_cost(tags, sol))





