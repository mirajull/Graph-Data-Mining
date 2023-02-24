import networkx as nx
import sys
import operator

def calculateWeight(c):
    return float(2 * nx.number_of_edges(c) / nx.number_of_nodes(c))

def orderVertices(g):
    d = nx.pagerank(g)
    sorted_v = list(map(lambda x: x[0],sorted(d.items(), key = operator.itemgetter(1), reverse=True)))
    return sorted_v

def linkAggregate(G):
	clusters = []
	vertex = orderVertices(G)

	for v in vertex:
		add = False

		for j in range(len(clusters)):
			U = clusters[j] + [v]
			UW = calculateWeight(G.subgraph(U))
			W = calculateWeight(G.subgraph(clusters[j]))
			if UW > W:
				clusters[j] += [v]
				add = True

		if add == False:
			clusters.append([v])

	return clusters


def improvedIterativeScan(cluster,G):
	cur = G.subgraph(cluster)
	W = calculateWeight(cur)
	increasing = True

	while increasing:
		N = cur.nodes()
		for vertex in N:
			adj = G.neighbors(vertex)
			N = list(set(N).union(set(adj)))

		for vertex in N:
			original_vertex = cur.nodes()
			if vertex in original_vertex:
				original_vertex.remove(vertex)
			else:
				original_vertex.append(vertex)
			if not original_vertex:
				new_cur_w=0
			else:
				new_cur = G.subgraph(original_vertex)
				new_cur_w = calculateWeight(new_cur)
			cur_w = calculateWeight(cur)
			if new_cur_w > cur_w:
				cur = new_cur.copy()
		new_W = calculateWeight(cur)

		if new_W == W:
			increasing = False
		else:
			W = new_W

	return cur.nodes()

def main():
    file = sys.argv[1]
    g = nx.Graph()
    with open(file) as f:
        next(f)
        for line in f:
            line = line.split()
            g.add_edge(int(line[0]),int(line[1]))
 
    clusters = linkAggregate(g)
    final_clusters = []
 
    for cluster in clusters:
        final_clusters.append(improvedIterativeScan(cluster, g))
 
    final = []
    for fc in final_clusters:
        fc = sorted (fc)
        if fc not in final:
            final.append(fc)
 
    with open("./result.txt", 'w') as f:
        for fwd in final:
            line = " ".join(map(str, fwd))
            f.write(line + '\n')


if __name__ == "__main__":
    main()