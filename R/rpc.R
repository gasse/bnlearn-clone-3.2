
rpc.global = function(targets, data, level, test, alpha, B, debug=FALSE, pc.method, nbr.join) {
  
#   ed.to.sep = 1:length(edges)
#   plot((ed.to.sep-1) %/% nb.nodes + 1, (ed.to.sep-1) %% nb.nodes + 1, ylim=c(nb.nodes, 1), xlim=c(1, nb.nodes))
#   ed.to.sep = ed.to.sep[edges[ed.to.sep] > 0]
#   points((ed.to.sep-1) %/% nb.nodes + 1, (ed.to.sep-1) %% nb.nodes + 1, pch="+")
#   
#   edges[-ed.to.sep] = 0
#   graphviz.plot(mat2pdag(edges, nodes), shape="ellipse")
  
  nodes = names(data)
  nb.nodes = length(nodes)
  
  targ.ids = which(nodes %in% targets)
  nb.targets = length(targets)
  
  if (nb.nodes < 2) {
    mb = list()
    mb[[targ.ids]] = list(nbr = character(0), mb = character(0))
    return(mb)
  }
  
  
  # STEP 0 : start with a disconnected graph
  edges = matrix(as.integer(0), nrow=nb.nodes, ncol=nb.nodes)
  edges.str = matrix(0, nrow=nb.nodes, ncol=nb.nodes)
  spouses = matrix(FALSE, nrow=nb.nodes, ncol=nb.nodes) # we don't care who are the common children
  
  todo = rep(FALSE, nb.nodes)
  todo[targ.ids] = TRUE
  
  
  # STEP 1 : 0-degree targets' (n-levels) neighbourhood discovery
  done = rep(FALSE, nb.nodes)
  
  # go 1 level deeper
  for (l in 1:(level + 1)) {
    for (i in (1:nb.nodes)[todo & !done]) {
      
      done[i] = TRUE
      for (j in (1:nb.nodes)[!done]) {
        
        a = conditional.test(nodes[i], nodes[j], NULL, data = data,
                             test = test, B = B, alpha = alpha, debug = debug)
        if (a <= alpha) {
          edges[i, j] = edges[j, i] = as.integer(1)
          edges.str[i, j] = edges.str[j, i] = 1 - a/alpha
          todo[j] = TRUE
        }#THEN
      }#FOR
    }#FOR
  }#FOR
  
  
  # STEP 2 : 1-degree targets' (n-levels) neighbourhood filtering
  ed.to.sep = which(edges > 0)
  ed.to.sep = ed.to.sep[(ed.to.sep-1) %% nb.nodes + 1 > (ed.to.sep-1) %/% nb.nodes + 1] # avoid double-checks
  ed.to.sep = ed.to.sep[order(edges.str[ed.to.sep])] # try to remove weakest edges first
  
  for (edge in ed.to.sep) {
    
    i = (edge-1) %% nb.nodes + 1
    j = (edge-1) %/% nb.nodes + 1
    
    i.seps = edges[i, ] > 0
    j.seps = edges[j, ] > 0
    
    # 1-degree separation is possible with common neighbours only
    edge.seps = which(i.seps & j.seps)
    # try to separate by strongest neighbours first
    edge.seps = edge.seps[order(edges.str[i, edge.seps] * edges.str[j, edge.seps], decreasing=TRUE)]
    
    for (k in edge.seps) {
      
      a = conditional.test(nodes[i], nodes[j], nodes[k],
                           data = data, test = test, B = B,
                           alpha = alpha, debug = debug)
      if (a > alpha) {
        edges[i, j] = edges[j, i] = as.integer(-k) # trick: save separator node id there
        edges.str[i, j] = edges.str[j, i] = 0
        break # TODO: try to separate with all neighbours and record them ?
      }#THEN
      else {
        edges.str[i, j] = edges.str[j, i] = min(edges.str[i, j], 1 - a/alpha)
      }#ELSE
    }#FOR
  }#FOR
  
  targets.nbr = apply(edges > 0, 2, any)
  
  # STEP 3 : 0-degree target's (n-levels) neighbourhood's neighbourhood (potential spouses) discovery
  border.nbr = todo & !done & apply(edges > 0, 1, any)
  
  for (i in which(border.nbr)) {
    
    done[i] = TRUE
    for (j in (1:nb.nodes)[!done]) {
      
      a = conditional.test(nodes[i], nodes[j], NULL, data = data,
                           test = test, B = B, alpha = alpha, debug = debug)
      if (a <= alpha) {
        edges[i, j] = edges[j, i] = as.integer(1)
        edges.str[i, j] = edges.str[j, i] = 1 - a/alpha
      }#THEN
    }#FOR
  }#FOR
  
  # STEP 4 : 1-degree target's (n-levels) neighbourhood's neighbourhood (potential spouses) filtering
  ed.to.sep = which(edges > 0 & rep(border.nbr, times=nb.nodes) & rep(border.nbr, each=nb.nodes))
  ed.to.sep = ed.to.sep[(ed.to.sep-1) %% nb.nodes + 1 > (ed.to.sep-1) %/% nb.nodes + 1] # avoid double-checks
  ed.to.sep = ed.to.sep[order(edges.str[ed.to.sep])] # try to remove weakest edges first
  
  for (edge in ed.to.sep) {
    
    i = (edge-1) %% nb.nodes + 1
    j = (edge-1) %/% nb.nodes + 1
    
    i.seps = edges[i, ] > 0
    j.seps = edges[j, ] > 0
    
    # 1-degree separation is possible with common neighbours only
    edge.seps = which(i.seps & j.seps)
    # try to separate by strongest neighbours first
    edge.seps = edge.seps[order(edges.str[i, edge.seps] * edges.str[j, edge.seps], decreasing=TRUE)]
    
    for (k in edge.seps) {
      
      a = conditional.test(nodes[i], nodes[j], nodes[k],
                           data = data, test = test, B = B,
                           alpha = alpha, debug = debug)
      if (a > alpha) {
        edges[i, j] = edges[j, i] = as.integer(-k) # trick: save separator node id there
        edges.str[i, j] = edges.str[j, i] = 0
        break # TODO: try to separate with all neighbours and record them ?
      }#THEN
      else {
        edges.str[i, j] = edges.str[j, i] = min(edges.str[i, j], 1 - a/alpha)
      }#ELSE
    }#FOR
  }#FOR
  
  # Step 5 : 1-2-degree (targets + targets' (n-levels) neighbourhoods)' spouses
  done = rep(FALSE, nb.nodes)
  for (i in which(targets.nbr)) {
    
    done[i] = TRUE
    
    # node neighbourhood
    i.nbr = which(edges[i, ] > 0)
    i.nbr = i.nbr[order(edges[i, i.nbr], decreasing=TRUE)]
    
    for (j in i.nbr) {
      
      # my neighbour's neighbours', not my neighbours, not separated by
      # my neighbour, not already my spouse, and not already processed
      i.j.sp = which((edges[j, ] > 0) & !(edges[i, ] > 0) &
                       !(edges[j, ] == -i) & !(spouses[i, ]) & !done)
      
      # parse potential spouses (neighbours' neighbours' \ already or impossible spouses)
      for (k in i.j.sp) {
        
        sep = -edges[j, k]
        
        a = conditional.test(nodes[i], nodes[k], nodes[c(j, if(sep > 0) sep else NULL)], data = data,
                             test = test, B = B, alpha = alpha, debug = debug)
        if (a <= alpha) {
          spouses[i, k] = spouses[k, i] = TRUE
        }#THEN
      }#FOR
    }#FOR
  }#FOR
  
  # Step 6 : PC filtering from target's MB
  
  # parse targets' edges
  # choose smallest MB ? choose both ? (then OR)
  # perform reverse-MB
  # perform PC filtering
  
  neighbours = matrix(FALSE, ncol=nb.nodes, nrow=nb.nodes)
  for (i in which(targets.nbr)) {
    
    t = nodes[i]
    pcs = nodes[edges[i, ] > 0]
    rsps = nodes[spouses[i, ]]
    
    # early stop : avoid PC search from complete OR-neighbourhood
    if (all(neighbours[i, edges[i, ] > 0]))
      next
    
    if(length(c(pcs, rsps)) < 3) {
      neighbours[i, edges[i, ] > 0] = TRUE
      neighbours[edges[i, ] > 0, i] = TRUE
      next
    }
    
    # my neighbours in my MB
    # TODO: start with my two first neighbours
    #   start = which(edges[i, ] > 0)
    #   start = start[order(edges[i, start], decreasing=TRUE)]
    #   start = start[1:min(2, length(start))]
    #   start = nodes[start]
    
    if (pc.method == "mmpc") {
      
      # 1. [Forward Phase (I)]
      pct = maxmin.pc.forward.phase(x=t, data=data, nodes=c(t, pcs, rsps),
                                    alpha=alpha, B=B, whitelist=NULL, blacklist=NULL,
                                    test=test, debug=debug)
      
      # 2. [Backward Phase (II)]
      pct = hybrid.pc.filter(t=t, pcs=pct, rsps=NULL, data=data, alpha=alpha, B=B,
                             whitelist=NULL, blacklist=NULL, backtracking=NULL,
                             test=test, debug=debug)
      
      # Logical AND : add nodes only if we can both see each other as a neighbour (MMPC correctness)
      for (n in pct) {
        
        j = which(nodes == n)
        
        # early stop : avoid PC search from OR-neighbourhood
        if (neighbours[i, j])
          next
        
        # those neighbours should also see my in my MB
        # TODO: start with me (and the other first neighbour if any)
        
        # 1. [Forward Phase (I)]
        pcn = maxmin.pc.forward.phase(x = n, data = data, nodes = c(t, pcs, rsps),
                                      alpha = alpha, B = B, whitelist = NULL, blacklist = NULL,
                                      test = test, debug = debug)
        
        # early stop : avoid useless filtering
        if (!t %in% pcn)
          next
        
        # 2. [Backward Phase (II)]
        pct = hybrid.pc.filter(t=n, pcs=pcn, rsps=NULL, data=data, alpha=alpha, B=B,
                               whitelist=NULL, blacklist=NULL, backtracking=NULL,
                               test=test, debug=debug)
        
        if (!t %in% pct)
          next
        
        neighbours[i, j] = TRUE
        if (nbr.join == "OR")
          neighbours[j, i] = TRUE
        
      }#FOR
    }
    else if (pc.method == "fdr.iapc") {
      
      #optimisation : avoid the 2 first rounds of the PC algorithm
      start = which(edges[i, ] > 0)
      start = start[order(edges[i, start], decreasing=TRUE)]
      start = start[1:min(2, length(start))]
      start = nodes[start]
      
      # 3. [PC] Get the Parents and Children from nodes within PCS and RSPS
      pc = hybrid.pc.nbr.search(t = t, data = data, nodes = c(t, pcs, rsps),
                                test = test, alpha = alpha, B = B, whitelist = NULL,
                                blacklist = NULL, start = start, backtracking = NULL,
                                debug = debug, method = pc.method)
      
      j = which(nodes %in% pc)
      neighbours[i, j] = TRUE
      if (nbr.join == "OR")
        neighbours[j, i] = TRUE
      
      # 4. [Neighbourhood OR] Detect and add false-negatives to PC, by checking if
      #     the target is present in potential neighbours' neighbours
      for (n in pcs[!pcs %in% pc]) {
        
        j = which(nodes == n)
        
        # early stop (avoid PC search from OR-neighbourhood)
        if (neighbours[i, j])
          next
        
        pcn = hybrid.pc.nbr.search(t = n, data = data, nodes = c(t, pcs, rsps),
                                   test = test, alpha = alpha, B = B, whitelist = NULL,
                                   blacklist = NULL, backtracking = NULL,
                                   debug = debug, method = pc.method, looking.for = t)
        
        # Logical OR : add the nodes from my PCS which I don't see
        # in my PC but which see me in theirs
        if(t %in% pcn) {
          
          neighbours[i, j] = TRUE
          if (nbr.join == "OR")
            neighbours[j, i] = TRUE
          
        }#THEN
        
      }#FOR
    }#THEN
  }#FOR
  
  # neighbourhood filtering
  neighbours = neighbours & t(neighbours)
  
  mb = list()
  todo = rep(FALSE, nb.nodes)
  done = rep(FALSE, nb.nodes)
  todo[targ.ids] = TRUE
  
  for (l in 1:(level+1)) {
    for (i in which(todo & !done)) {
      done[i] = TRUE
      todo[neighbours[i, ]] = TRUE
      mb[[nodes[i]]] = neighbours[i, ]
    }#FOR
  }#FOR
  
  mb = lapply(mb, function(x){
    list(
      nbr = nodes[x & done],
      mb = character(0))
  })
  
  return(mb)
  
}#RPC.GLOBAL
