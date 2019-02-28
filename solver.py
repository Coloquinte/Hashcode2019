#!/usr/bin/python3

import sys

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

def slide_tags(tags, slide):
    if len(slide) == 1:
      return tags[slide[0]]
    else:
      return tags[slide[0]] & tags[slide[1]]

def cost_between_tags(ti, tj):
	return min( len(ti&tj), len(ti-tj), len(tj-ti) )

def sol_cost(tags, sol):
    sol_tags = [slide_tags(tags, s) for s in sol]
    cost = 0
    for i in range(len(sol) - 1):
      cost += cost_between_tags(sol_tags[i], sol_tags[i+1])
    return cost

if len(sys.argv) < 2:
  print("solver.py input [output]")
  sys.exit(1)

tags, horizontal = parse_file(sys.argv[1])

sol = [(i,) for i in range(len(tags))]

print (sol_cost(tags, sol))

if len(sys.argv) >= 3:
  write_solution(sol, sys.argv[2])





