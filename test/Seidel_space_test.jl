G = GKMproj_space(2)
w = gens(G.M)
cMax = 6

println("S1")
S1 = Seidel_element(G, w[1], cMax)
println("S2")
S2 = Seidel_element(G, w[2], cMax)
println("S3")
S3 = Seidel_element(G, w[1]+w[2], cMax)
println("S4")
S4 = Seidel_element(G, 2*w[1], cMax)

for b in keys(S4)
  println("$b: ")
  for f in S4[b]
    if f == 0
      println("\t0")
    else
      println("\t$(factor(numerator(f))) // $(factor(denominator(f)))")
    end
  end
end

SG4 = Seidel_space(G, 2*w[1])
P = ev(1, pointClass(4, SG4))
integrateGKM(SG4, edgeCurveClass(SG4, Edge(1, 4)), 1, P; show_bar=false)