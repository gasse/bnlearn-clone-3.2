
rpc.global = function(targets, data, levels, test, alpha, B, debug=FALSE) {
  
#   edg.ids = 1:length(edges)
#   plot((edg.ids-1) %/% nb.nodes + 1, (edg.ids-1) %% nb.nodes + 1, ylim=c(nb.nodes, 1), xlim=c(1, nb.nodes))
#   edg.ids = edg.ids[edges[edg.ids] > 0]
#   points((edg.ids-1) %/% nb.nodes + 1, (edg.ids-1) %% nb.nodes + 1, pch="+")
#   
#   edges[-edg.ids] = 0
#   graphviz.plot(mat2pdag(edges, nodes), shape="ellipse")
  
  pval.prec = 10000
  nodes = names(data)
  nb.nodes = length(nodes)
  
  targ.ids = which(nodes %in% targets)
  nb.targets = length(targets)
  
  # start with a completely unconnected graph
  edges = matrix(as.integer(0), nrow=nb.nodes, ncol=nb.nodes)
  spouses = matrix(FALSE, nrow=nb.nodes, ncol=nb.nodes) # just X spouse of Y, we don't care which nodes are common children
  
  if (nb.nodes > 1) {
    
    # Step 1 : 0-degree target's neighbourhood filtering
    
    done = rep(FALSE, nb.nodes)
    for (i in targ.ids) {
      
      done[i] = TRUE
      for (j in (1:nb.nodes)[!done]) {
        
        a = conditional.test(nodes[i], nodes[j], NULL, data = data,
                             test = test, B = B, alpha = alpha, debug = debug)
        if (a <= alpha) {
          edges[i, j] = edges[j, i] = as.integer(pval.prec * (1 - a/alpha))
        }
      }
    }
    
    # Step 2 : 1-degree target's neighbourhood filtering
    
    # try to separate weakest neighbours first
    edg.ids = seq(length(edges))
#     edg.ids = edg.ids[(edg.ids-1) %% nb.nodes < (edg.ids-1) %/% nb.nodes]
    edg.ids = edg.ids[((edg.ids-1) %% nb.nodes + 1) %in% targ.ids]
    edg.ids = edg.ids[order(edges[edg.ids])]
    edg.ids = edg.ids[edges[edg.ids] > 0]
    
    for (edge in  edg.ids) {
      
      i = (edge-1) %% nb.nodes + 1
      j = (edge-1) %/% nb.nodes + 1
      
      i.seps = c(edges[1:i, i] > 0, edges[i, i:nb.nodes] > 0)[-i]
      j.seps = c(edges[1:j, j] > 0, edges[j, j:nb.nodes] > 0)[-j]
      
      # try to separate with common neighbours, strongest first
      seps = which(i.seps & j.seps)
      seps = seps[order(edges[i, seps] * edges[j, seps], decreasing=TRUE)]
      
      if (length(seps) > 0) {
        for (k in 1:length(seps)) {
          
          a = conditional.test(nodes[i], nodes[j], nodes[seps[k]],
                               data = data, test = test, B = B,
                               alpha = alpha, debug = debug)
          if (a > alpha) {
            edges[i, j] = edges[j, i] = - seps[k] # trick: save separator id there
            break # TODO: try to separate with all neighbours and keep them ?
          }
          else {
            edges[i, j] = edges[j, i] = min(edges[i, j], as.integer(pval.prec * (1 - a/alpha)))
          }
        }
      }
    }
    
    targ.nbr = apply(edges[targ.ids, ] > 0, 2, any)
    
    # Step 3 : 0-degree target's neighbourhood's neighbourhood (potential spouses) filtering
    
    done = rep(FALSE, nb.nodes)
    done[targ.ids] = TRUE
    for (i in which(targ.nbr & !done)) {
      
      done[i] = TRUE
      for (j in (1:nb.nodes)[!done]) {
        
        a = conditional.test(nodes[i], nodes[j], NULL, data = data,
                             test = test, B = B, alpha = alpha, debug = debug)
        if (a <= alpha) {
          edges[i, j] = edges[j, i] = as.integer(pval.prec * (1 - a/alpha))
        }
      }
    }
    
    # Step 4 : 1-degree target's neighbourhood's neighbourhood (potential spouses) filtering
    
    done = rep(FALSE, nb.nodes)
    done[targ.ids] = TRUE
    
    # try to separate weakest neighbours first
    edg.ids = seq(length(edges))
    #     edg.ids = edg.ids[(edg.ids-1) %% nb.nodes < (edg.ids-1) %/% nb.nodes]
    edg.ids = edg.ids[((edg.ids-1) %% nb.nodes + 1) %in% which(targ.nbr & !done)]
    edg.ids = edg.ids[order(edges[edg.ids])]
    edg.ids = edg.ids[edges[edg.ids] > 0]
    
    for (edge in  edg.ids) {
      
      i = (edge-1) %% nb.nodes + 1
      j = (edge-1) %/% nb.nodes + 1
      
      i.seps = c(edges[1:i, i] > 0, edges[i, i:nb.nodes] > 0)[-i]
      i.seps[j] = FALSE
      i.seps = which(i.seps)
      
      # try to separate with strongest neighbours first
      i.seps = i.seps[order(edges[i, i.seps], decreasing=TRUE)]
      
      if (length(i.seps) > 0) {
        for (k in 1:length(i.seps)) {
          
          a = conditional.test(nodes[i], nodes[j], nodes[i.seps[k]],
                               data = data, test = test, B = B,
                               alpha = alpha, debug = debug)
          if (a > alpha) {
            edges[i, j] = edges[j, i] = - i.seps[k] # trick: save separator id there
            break # TODO: try to separate with all neighbours and keep them ?
          }
          else {
            edges[i, j] = edges[j, i] = min(edges[i, j], as.integer(pval.prec * (1 - a/alpha)))
          }
        }
      }
    }
    
    # Step 5 : 1-2-degree target's spouses
    
    # parse targets
    # parse their neighbours
    # parse potential spouses (neighbours' neighbours' \ already or impossible spouses)
    done = rep(FALSE, nb.nodes)
    for (i in targ.ids) {
      
      done[i] = TRUE
      
      # my neighbours'
      i.nbr = which(edges[i, ] > 0)
      i.nbr = i.nbr[order(edges[i, i.nbr], decreasing=TRUE)]
      
      for (j in i.nbr) {
        
        # my neighbour's neighbours', not my neighbours, not separated by
        # my neighbour, not already my spouse, and not already processed
        i.j.sp = which((edges[j, ] > 0) & !(edges[i, ] > 0) &
                         !(edges[j, ] == -i) & !(spouses[i, ]) & !done)
        
        for (k in i.j.sp) {
          
          sep = -edges[j, k]
          
          a = conditional.test(nodes[i], nodes[k], nodes[c(j, if(sep > 0) sep else NULL)], data = data,
                               test = test, B = B, alpha = alpha, debug = debug)
          if (a <= alpha) {
            spouses[i, k] = spouses[k, i] = TRUE
          }
        }
      }
    }
  }
    
  # Step 6 : PC filtering from target's MB
  
  # parse targets' edges
  # choose smallest MB ? choose both ? (then OR)
  # perform reverse-MB
  # perform PC filtering
  
  neighbours = matrix(FALSE, ncol=nb.nodes, nrow=nb.nodes)
  for (i in targ.ids) {
    
    t = nodes[i]
    pcs = nodes[edges[i, ] > 0]
    rsps = nodes[spouses[i, ]]
    
    if(length(c(pcs, rsps)) < 3) {
      neighbours[i, which(edges[i, ] > 0)] = TRUE
      neighbours[which(edges[i, ] > 0), i] = TRUE
    }
    
    # my neighbours in my MB
    # TODO: start with my two first neighbours
    
    # 1. [Forward Phase (I)]
    pct = maxmin.pc.forward.phase(x=t, data = data, nodes = c(t, pcs, rsps),
                                  alpha = alpha, B = B, whitelist = NULL, blacklist = NULL,
                                  test = test, debug = debug)
    
    # 2. [Backward Phase (II)]
    pct = hybrid.pc.filter(t=t, pcs=pct, rsps=NULL, data=data, alpha=alpha, B=B,
                           whitelist=NULL, blacklist=NULL, backtracking=NULL,
                           test=test, debug=debug)
    
    for (n in pct) {
      
      # those neighbours should also see my in my MB
      # TODO: start with me (and the other first neighbour if any)
      
      # 1. [Forward Phase (I)]
      pcn = maxmin.pc.forward.phase(x = n, data = data, nodes = c(t, pcs, rsps),
                                    alpha = alpha, B = B, whitelist = NULL, blacklist = NULL,
                                    test = test, debug = debug)
      
      # 2. [Backward Phase (II)]
      pct = hybrid.pc.filter(t=n, pcs=pcn, rsps=NULL, data=data, alpha=alpha, B=B,
                             whitelist=NULL, blacklist=NULL, backtracking=NULL,
                             test=test, debug=debug)
      
      # Logical AND : add nodes only if we can both see each other as a neighbour (MMPC correctness)
      if(t %in% pcn) {
        
        j = which(nodes == n)
        neighbours[i, j] = neighbours[j, i] = TRUE
        
      }#THEN
      
    }#FOR
  }
  
  mb = lapply(1:nb.nodes, function(i) {
    list(
      nbr = nodes[neighbours[i, ]],
      mb = character(0))
  })
  names(mb) = nodes
  
  return(mb)
  
}#RPC.GLOBAL
