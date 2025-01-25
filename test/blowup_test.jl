G = GKMproj_space(3)
con = build_GKM_connection(G)

S = GKMsubgraph_from_vertices(G, [1, 2])
(blowupSub, blowupCon) = blowupGKM(S, con)

Spoint = GKMsubgraph_from_vertices(G, [1])
(blowupPt, blowupConPt) = blowupGKM(Spoint, con)