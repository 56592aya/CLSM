#Pkg.add("Graphs")
using Graphs
using PyPlot
g = simple_graph(3)
add_edge!(g, 1, 2)
add_edge!(g, 3, 2)
add_edge!(g, 3, 1)
plot(g)
num_vertices(g)
num_edges(g)
for v in vertices(g)
  for e in out_edges(v,g)
    t=target(e,g)
    println("$v is connected to $t")
  end
end

adjacency_matrix(g)
adjacency_matrix_sparse(g)

pairs = [(1,2), (1,3), (2,3), (2,4), (3,5), (4,5), (2,5)]
eds = Edge{Int}[Edge(i,p[1],p[2]) for (i,p) in enumerate(pairs)]
println(eds)
for e in eds
  println(e.source," is connected to ",e.target)
end
typeof(eds)
g_eds_u = simple_edgelist(5, eds; is_directed=false)

for v in vertices(g_eds_u)
  for e in edges(g_eds_u)
    println(e.source,"  is connected to   ", e.target)
  end
end
