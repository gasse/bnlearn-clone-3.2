
pc.global = function(data, test, alpha, B, debug=FALSE) {
  
  nodes = names(data)
  nb.nodes = length(nodes)
  
  # start with a completely connected graph
  edges = diag(nrow=nb.nodes, ncol=nb.nodes)
  edges = edges + 1 - 2 * edges
  
  if (nb.nodes >1) {
    
    # 0-degree neighbourhood filtering
    for (i in 1:(nb.nodes-1)) {
      for (j in (i+1):nb.nodes) {
        
        a = conditional.test(nodes[i], nodes[j], NULL, data = data, test = test, B = B,
                             alpha = alpha, debug = debug)
        if (a > alpha) {
          edges[i, j] = edges[j, i] = 0
        }
        else {
          edges[i, j] = edges[j, i] = min(edges[i, j], 1 - a/alpha)
        }
      }
    }
    
    # k-degree neighbourhood filtering
    for (n in 1:(nb.nodes-1)) {
      
      # process edges from weakest to strongest
      edg.ids = 1:length(edges)
      edg.ids = edg.ids[(edg.ids-1) %% nb.nodes < (edg.ids-1) %/% nb.nodes]
      edg.ids = edg.ids[order(edges[edg.ids])]
      edg.ids = edg.ids[edges[edg.ids] > 0]
      
      for (edge in  edg.ids) {
        
        i = (edge-1) %% nb.nodes + 1
        j = (edge-1) %/% nb.nodes + 1
        
        i.nbr = (edges[i, ] > 0)
        j.nbr = (edges[j, ] > 0)
        
        combs = NULL
        if(n == 1) {
          seps = which(i.nbr & j.nbr)
          seps = seps[order(edges[i, seps] * edges[j, seps])]
          combs = matrix(seps, nrow=1)
        }
        else {
          if (sum(i.nbr[-j] > 0) >= n) {
            # try to break edges with strongest to weakest neighbours
            seps = which(i.nbr & 1:nb.nodes != j)
            seps = seps[order(edges[i, seps])]
            combs = cbind(combs, combn(seps, n))
          }
          if (sum(j.nbr[-i] > 0) >= n) {
            # try to break edges with strongest to weakest neighbours
            seps = which(j.nbr & 1:nb.nodes != i)
            seps = seps[order(edges[j, seps])]
            combs = cbind(combs, combn(seps, n))
          }
        }
        
        if (length(combs) > 0) {
          for (k in 1:ncol(combs)) {
            
            a = conditional.test(nodes[i], nodes[j], nodes[combs[, k]],
                                 data = data, test = test, B = B,
                                 alpha = alpha, debug = debug)
            if (a > alpha) {
              edges[i, j] = edges[j, i] = 0
              break
            }
            else {
              edges[i, j] = edges[j, i] = min(edges[i, j], 1 - a/alpha)
            }
          }
        }
      }
      
      if (max(apply(edges > 0, 1, sum)) - 1 > n)
        next
      
      break
    }
  }
  
  mb = lapply(1:ncol(edges), function(i) {
    list(nbr = nodes[edges[i, ] > 0], mb = nodes[edges[i, ] > 0])
  })
  names(mb) = nodes
  
  return(mb)
  
}#PC.GLOBAL
