########## tsp.py ##########

import localsolver
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


tags, horizontal = parse_file(sys.argv[1])

appar = [(i,) for i in range(len(tags))]

appar = appar[:2000]

nb_cities = len(appar)

if len(sys.argv) >= 3:
  write_solution(sol, sys.argv[2])


def get_tags_from_slide(tup):
    if len(tup) == 1:
        return tags[tup[0]]
    if len(i) == 2:
        return tags[tup[0]] | tags[tup[1]]

def cost(i, j):
    ti = get_tags_from_slide(i)
    tj = get_tags_from_slide(j)
    return min( len(ti&tj), len(ti-tj), len(tj-ti) )



if len(sys.argv) < 2:
    print ("Usage: python tsp.py inputFile [outputFile] [timeLimit]")
    sys.exit(1)


def read_elem(filename):
    with open(filename) as f:
        return [str(elem) for elem in f.read().split()]


with localsolver.LocalSolver() as ls:

    #
    # Reads instance data 
    #

    file_it = iter(read_elem(sys.argv[1]))

    # The input files follow the TSPLib "explicit" format.
    for pch in file_it:
        if (pch == "DIMENSION:"):
            nb_cities = int(next(file_it))
        if (pch == "EDGE_WEIGHT_SECTION"):
            break

    # Distance from i to j
    distance_weight = [[cost(appar[i], appar[j]) for i in range(nb_cities)] for j in range(nb_cities)] 

    #
    # Declares the optimization model
    #
    model = ls.model

    # A list variable: cities[i] is the index of the ith city in the tour
    cities = model.list(nb_cities) 

    # All cities must be visited
    model.constraint(model.count(cities) == nb_cities)

    # Create a LocalSolver array for the distance matrix in order to be able to 
    # access it with "at" operators.
    distance_array = model.array(distance_weight)

    # Minimize the total distance
    dist_selector = model.function(lambda i: model.at(distance_array, cities[i-1], cities[i]))
    obj = (model.sum(model.range(1, nb_cities), dist_selector)
            + model.at(distance_array, cities[nb_cities - 1], cities[0]));
    model.maximize(obj)

    model.close()

    #
    # Parameterizes the solver
    #
    if len(sys.argv) >= 4: ls.param.time_limit = int(sys.argv[3])
    else: ls.param.time_limit = 50

    ls.solve()

    #
    # Writes the solution in a file
    #
    if len(sys.argv) >= 3:
        # Writes the solution in a file 
        with open(sys.argv[2], 'w') as f:
            f.write("%d\n" % obj.value)
            for c in cities.value:
                f.write("%d " % c)
            f.write("\n")