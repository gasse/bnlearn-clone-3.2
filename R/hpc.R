
hybrid.pc.global = function(
  data, whitelist, blacklist, test, alpha, B, strict, pc.method, nbr.join,
  debug=FALSE) {
  
  nodes = names(data)
  
  mb = lapply(as.list(nodes), hybrid.pc, data = data, nodes = nodes,
              alpha = alpha, B = B, whitelist = whitelist, blacklist = blacklist,
              test = test, debug = debug, pc.method = pc.method)
  names(mb) = nodes
  
  # check neighbourhood sets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, debug = debug,
                   filter = nbr.join)
  
  return(mb)
  
}#HYBRID.PC.GLOBAL

hybrid.pc.global.optimized = function(
  data, whitelist, blacklist, test, alpha, B, strict, pc.method, nbr.join,
  debug=FALSE) {
  
  nodes = names(data)
  mb = list()
  
  for(node in nodes) {
    
    backtracking = unlist(sapply(mb, function(x){ node %in% x$nbr  }))
    
    # depending on the neighbourhood consistency filter used, a full
    # backtracking is prohibited.
    #   AND filter : known good backtracking is forbidden
    #   OR filter : known bad backtracking is forbidden
    if (!is.null(backtracking)) {
      if (nbr.join == "AND")
        backtracking = backtracking[!backtracking]
      if (nbr.join == "OR")
        backtracking = backtracking[backtracking]
      if (length(backtracking) == 0)
        backtracking = NULL
    }#THEN
    
    mb[[node]] = hybrid.pc(node, data = data, nodes = nodes, alpha = alpha,
                           B = B, whitelist = whitelist, blacklist = blacklist,
                           backtracking = backtracking, test = test,
                           debug = debug, pc.method = pc.method)
    
  }#FOR
  
  # check neighbourhood sets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, debug = debug,
                   filter = nbr.join)
  
  return(mb)
  
}#HYBRID.PC.GLOBAL.OPTIMIZED

hybrid.pc.global.cluster = function(
  data, cluster, whitelist, blacklist, test, alpha, B, strict, pc.method,
  nbr.join, debug=FALSE) {
  
  nodes = names(data)
  
  mb = clusterApplyLB(cluster, as.list(nodes), hybrid.pc, data = data,
                      nodes = nodes, alpha = alpha, B = B, whitelist = whitelist,
                      blacklist = blacklist, test = test, debug = debug, pc.method = pc.method)
  names(mb) = nodes
  
  # check neighbourhood sets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, debug = debug,
                   filter = nbr.join)
  
  return(mb)
  
}#HYBRID.PC.GLOBAL.CLUSTER

hybrid.pc = function(t, data, nodes, whitelist, blacklist, test, alpha,
                     B, pc.method, backtracking = NULL, debug = FALSE) {
  
  # 1. [PCS] Search parents and children superset
  tmp = hybrid.pc.de.pcs(t, data, nodes, alpha, B, whitelist, blacklist,
                         backtracking, test, debug)
  pcs = tmp$pcs
  pvals = tmp$pvals
  dsep = tmp$dsep
  
  #optimisation : 0 or 1 node in PCS --> PC == PCS
  if(length(pcs) < 2)
    return(list(nbr = pcs, mb = NULL))
  
  # 2. [RSPS] Search remaining spouses superset, those not already in PCS
  rsps = hybrid.pc.de.rsps(t, data, nodes, pcs, dsep, alpha, B, test, debug)
  
  #optimisation : 2 nodes in PC and no SP --> PC == PCS
  if(length(c(pcs, rsps)) < 3)
    return(list(nbr = pcs, mb = pcs))
  
  #optimisation : avoid the 2 first rounds of the PC algorithm
  start = pcs[order(pvals)][1:min(length(pcs), 2)]
  
  # 3. [PC] Get the Parents and Children from nodes within PCS and RSPS
  pc = hybrid.pc.nbr.search(t = t, data = data, nodes = c(t, pcs, rsps),
                            test = test, alpha = alpha, B = B, whitelist = whitelist,
                            blacklist = blacklist, start = start, backtracking = backtracking,
                            debug = debug, method = pc.method)
  
  # 4. [Neighbourhood OR] Detect and add false-negatives to PC, by checking if
  #     the target is present in potential neighbours' neighbours
  for (node in pcs[!pcs %in% pc]) {
    
    pcn = hybrid.pc.nbr.search(t = node, data = data, nodes = c(t, pcs, rsps),
                               test = test, alpha = alpha, B = B, whitelist = whitelist,
                               blacklist = blacklist, backtracking = NULL,
                               debug = debug, method = pc.method, looking.for = t)
    
    # Logical OR : add the nodes from my PCS which I don't see
    # in my PC but which see me in theirs
    if(t %in% pcn) {
      
      pc = c(pc, node)
      mb = c(mb, node)
      
      if (debug)
        cat("  @", node, "added to the parents and children. (HPC's OR)\n")
      
    }#THEN
    
  }#FOR
  
  if (debug)
    cat(" * PC of", t, "=", pc, "\n")
  
  res = list(nbr = pc, mb = c(pcs, rsps))
  
  return(res)
  
}#HYBRID.PC

hybrid.pc.nbr.search = function(t, data, nodes, test, alpha, B,
                                whitelist = NULL, blacklist = NULL, backtracking = NULL, debug = FALSE,
                                start = character(0), method = "inter.iapc", looking.for = NULL) {
  
  pc = c()
  
  if (method == "mmpc") {
    
    # incorrect : MMPC should be re-run on each neighbour and combine the results
    mb = maxmin.pc.forward.phase(x=t, data = data, nodes = nodes, alpha = alpha,
                                 B = B, whitelist = whitelist, blacklist = blacklist,
                                 test = test, debug = debug)
    
  }#THEN
  else if (method == "si.hiton.pc") {
    
    mb = si.hiton.pc.heuristic(x=t, data = data, nodes = nodes, alpha = alpha,
                                B = B, whitelist = whitelist, blacklist = blacklist, test = test,
                                optimized = optimized, debug = debug)
    
  }#THEN
  else if (method == "gs.pc") {
    
    mb = gs.markov.blanket(x = t, data = data, nodes = nodes,
                           alpha = alpha, B = B, whitelist = whitelist, blacklist = NULL,
                           backtracking = NULL, test = test, debug = debug)
    
  }#THEN
  else if (method == "inter.iapc") {
    
    mb = inter.ia.markov.blanket(x = t, data = data, nodes = nodes,
                                 alpha = alpha, B = B, whitelist = whitelist, blacklist = blacklist,
                                 start = start, backtracking = backtracking, test = test, debug = debug)
    
  }#THEN
  else if (method == "iapc") {
    
    mb = ia.markov.blanket(x = t, data = data, nodes = nodes,
                           alpha = alpha, B = B, whitelist = whitelist, blacklist = blacklist,
                           start = start, backtracking = backtracking, test = test, debug = debug)
    
  }#THEN
  else if (method == "fast.iapc") {
    
    mb = fast.ia.markov.blanket(x = t, data = data, nodes = nodes,
                                alpha = alpha, B = B, whitelist = whitelist, blacklist = blacklist,
                                start = start, backtracking = backtracking, test = test, debug = debug)
    
  }#THEN
  else if (method == "fdr.iapc") {
    
    mb = iambfdr(x = t, data = data, nodes = nodes,
                 alpha = alpha, B = B, whitelist = whitelist, blacklist = blacklist,
                 start = start, backtracking = backtracking, test = test, debug = debug)
    
  }#THEN
  else {
    
    stop(paste("the neighbourhood method ", method, " is not supported yet by hpc.", sep=""))
    
  }#ELSE
  
  #optimisation : no need to filter if the node we are looking for is not in MB
  if (is.null(looking.for) || looking.for %in% mb)
    pc = hybrid.pc.filter(t, pcs = mb, rsps = NULL, data, alpha, B, whitelist,
                          blacklist, backtracking, test, debug)
  
  return(pc)
  
}#HYBRID.PC.NBR.SEARCH

hybrid.pc.de.pcs = function(t, data, nodes, alpha, B, whitelist, blacklist,
                            backtracking = NULL, test, debug = FALSE) {
  
  pcs = vector()
  pvals = vector()
  dsep = list()
  
  whitelisted = nodes[vapply(nodes,
                             function(x) { is.whitelisted(whitelist, c(t, x), either = TRUE) }, logical(1))]
  blacklisted = nodes[vapply(nodes,
                             function(x) { is.blacklisted(blacklist, c(t, x), both = TRUE) }, logical(1))]
  known.good = c()
  
  if (debug) {
    cat("----------------------------------------------------------------\n")
    cat("* learning the PCS of", t, ".\n")
  }#THEN
  
  # use backtracking for a further screening of the nodes to be checked.
  # don't use known bad nodes, we want to check them anyway because of
  # HPC's OR filter on neighbours
  if (!is.null(backtracking)) {
    
    known.good = names(backtracking[backtracking])
    
    if(debug) {
      cat(" * known good (backtracking): '", known.good, "'.\n")
    }#THEN
    
  }#THEN
  
  if (debug) {
    cat(" * nodes to be tested for inclusion: '",
        nodes[!(nodes %in% c(t, known.good))], "'.\n")
  }#THEN
  
  nodes.to.check = setdiff(nodes, c(t, known.good, whitelisted, blacklisted))
  
  # Phase (I): remove X if Ind(T,X) (0-degree d-separated nodes)
  for (x in nodes.to.check) {
    
    a = conditional.test(t, x, c(), data = data, test = test, B = B,
                         alpha = alpha, debug = debug)
    
    if (a <= alpha) {
      
      if (debug)
        cat("  @ node", x, "added to the parents and children superset\n")
      
      # keep pcs with p-values in order to order nodes in the next phase
      pcs = c(pcs, x)
      pvals = c(pvals, a)
      
    }#THEN
    
  }#FOR
  
  # heuristic 1 : sort the PC candidates in decreasing p-value order
  # this way we are more prone to remove less correlated nodes first
  ord = order(pvals, decreasing = TRUE)
  pcs = pcs[ord]
  pvals = pvals[ord]
  
  nodes.to.check = pcs
  
  if (debug) 
    cat(" * nodes to be tested for exclusion: '", nodes.to.check, "'.\n")
  
  pcs = union(known.good, pcs)
  pcs = union(whitelisted, pcs)
  pvals = c(rep(0, length(pcs) - length(pvals)), pvals)
  
  # heuristic 2 : sort the d-separating canditates in increasing p-value order
  # this way we are more prone to remove with highly correlated nodes first
  nodes.to.check.against = pcs[order(pvals)]
  
  # Phase (II): remove X if Ind(T,X|Y) (1-degree d-separated nodes)
  for (x in nodes.to.check) {
    
    # heuristic 3 : sort the d-separating canditates in increasing p-value order
    # this way we are more prone to remove with highly correlated nodes first
    nodes.to.check.against = setdiff(pcs[order(pvals)], x)
    
    for (y in nodes.to.check.against) {
      
      #     for (y in setdiff(nodes.to.check.against, x)) {
      
      a = conditional.test(t, x, y, data = data, test = test, B = B,
                           alpha = alpha, debug = debug)
      
      x.ind = which(pcs == x)
      if (a > alpha) {
        
        pcs = pcs[-x.ind]
        pvals = pvals[-x.ind]
        dsep[[x]] = c(y)
        
        if (debug) {
          cat("  @ node", x, "removed from the parents and children superset\n")
        }#THEN
        
        break
        
      }#THEN
      else {
        
        if(a > pvals[x.ind])
          pvals[x.ind] = a
        
      }#ELSE
    }#FOR
  }#FOR
  
  if (debug) {
    cat(" * PCS of", t, "= '", pcs, "'.\n")
  }#THEN
  
  ret = list(pcs=pcs, dsep=dsep, pvals = pvals)
  return(ret)
  
}#HYBRID.PC.DE.PCS

hybrid.pc.de.rsps = function(t, data, nodes, pcs, dsep, alpha, B, test,
                             debug = FALSE) {
  
  rsps = vector()
  
  if (debug) {
    
    cat("----------------------------------------------------------------\n")
    cat("* learning the RSPS of", t, ".\n")
    cat(" * PCS =", pcs, "\n")
    cat(" * nodes still to be tested for inclusion: ", nodes[!(nodes %in% c(t, pcs))], "\n")
    
  }#THEN
  
  for (x in pcs) {
    
    rspsx = vector()
    pvalsx = vector()
    
    # Phase (I): search spouses Y in the form T->X<-Y from the
    # remaining nodes (not in pcs)
    for (y in setdiff(nodes, c(t, pcs))) {
      
      # optimisation : avoid irrelevant tests
      if (x %in% dsep[[y]])
        next
      
      a = conditional.test(t, y, c(dsep[[y]], x), data = data, test = test,
                           B = B, alpha = alpha, debug = debug)
      
      if (a <= alpha) {
        
        rspsx = c(rspsx, y)
        pvalsx = c(pvalsx, a)
        
        if (debug) {
          
          cat("  @ node", y, "added to the spouses superset by", x, "\n")
          
        }#THEN
        
      }#THEN
    }#FOR
    
    # heuristic : sort the candidates in decreasing p-value order
    # this way we are more prone to remove less correlated nodes first
    rspsx = rspsx[order(pvalsx, decreasing = TRUE)]
    
    # Phase (II): remove false positive spouses Y in the form T->X<-Z<-...<-Y
    # (descendants or ancestors of Z)
    for (y in rspsx) {
      for (z in rspsx[rspsx != y]) {
        
        a = conditional.test(t, y, c(x, z), data = data, test = test, B = B,
                             alpha = alpha, debug = debug)
        
        if (a > alpha) {
          
          rspsx = rspsx[rspsx != y]
          
          if (debug) {
            
            cat("  @ node", y, "removed from the spouses superset by", z, "\n")
            
          }#THEN
          
          break
          
        }#THEN
      }#FOR
    }#FOR
    
    rsps = c(rsps, rspsx)
    
  }#FOR
  
  rsps = unique(rsps)
  
  if(debug)
    cat(" * RSPS =", rsps, "\n")
  
  return(rsps)
  
}#HYBRID.PC.DE.RSPS

# Build the parents and children (PC) set of a node from it's parents and
# children superset (PCS) and it's remaining spouses superset (RSPS).
# we assume intersection(t, pcs, rsps) is empty
hybrid.pc.filter = function(t, pcs, rsps, data, alpha, B = NULL,
                            whitelist, blacklist, backtracking = NULL, test, debug = FALSE) {
  
  nodes = names(data)
  
  known.good = c()
  blacklisted = nodes[vapply(nodes, function(x) {
    is.blacklisted(blacklist, c(t, x), both = TRUE) }, logical(1))]
  whitelisted = nodes[vapply(nodes, function(x) {
    is.whitelisted(whitelist, c(t, x), either = TRUE) }, logical(1))]
  
  # markov boundary superset (markov blanket)
  mbs = c(pcs, rsps)
  
  # use backtracking for a further screening of the nodes to be checked.
  # known bad nodes are not used, we want to keep it possible to add them,
  # because of HPC's OR filter on neighbours
  if (!is.null(backtracking)) {
    
    # neighbourhood is symmetric, so X \in N(Y) <=> Y \in N(X)
    known.good = names(backtracking[backtracking])
    
  }#THEN
  
  if (debug) {
    
    cat("----------------------------------------------------------------\n")
    cat("* filtering PC of", t, ".\n")
    cat(" * blacklisted nodes: '", blacklisted, "'\n")
    cat(" * whitelisted nodes: '", whitelisted, "'\n")
    cat(" * starting with neighbourhood superset: '", pcs, "'\n")
    cat(" * with remaining spouses superset: '", rsps, "'\n")
    
    if (!is.null(backtracking))
      cat(" * known good (backtracking): '", known.good, "'.\n")
    
  }#THEN
  
  # blacklisted nodes are removed
  pcs = pcs[!(pcs %in% blacklisted)]
  
  # whitelisted nodes are added
  pcs = union(pcs, whitelisted)
  mbs = union(mbs, whitelisted)
  
  pc = vector()
  if (length(mbs) > 0) {
    
    nbr = function(x) {
      
      # excluding the node to be tested for exclusion from the conditioning set
      dsep.set = setdiff(mbs, x)
      
      # try to d-separate with all the possible conditioning subsets
      for (k in 0:length(dsep.set)) {
        
        # create all possible subsets of the conditioning set of size k.
        dsep.subsets = subsets(length(dsep.set), k, dsep.set)
        
        for (s in 1:nrow(dsep.subsets)) {
          
          a = conditional.test(t, x, dsep.subsets[s,], data = data, test = test,
                               B = B, alpha = alpha, debug = debug)
          
          if (a > alpha) {
            
            if (debug) {
              cat("  @ node", x, "is not a neighbour of", t, " any more\n")
            }
            
            return(FALSE)
            
          }#THEN
          
        }#FOR
        
      }#REPEAT
      
      # all tests failed, the nodes could not be d-separated
      return(TRUE)
      
    }#NBR
    
    # apply the PC filtering
    # do not even try to remove whitelisted and backtracked (good) nodes.
    canpc = pcs[!(pcs %in% c(whitelisted, known.good))]
    pc = canpc[vapply(canpc, nbr, logical(1))]
    
    # whitelisted nodes are included (arc orientation is irrelevant),
    # same for known good nodes
    pc = unique(c(pc, whitelisted, known.good))
    
  }#THEN
  
  return(pc)
  
}#HYBRID.PC.FILTER
