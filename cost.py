def cost(i, j):
	ti = tags[i]
	tj = tags[j]
	return min( len(ti&tj), len(ti-tj), len(tj-ti) )