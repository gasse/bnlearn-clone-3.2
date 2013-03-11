
hybrid.pc.cached.global = function(
  data, whitelist, blacklist, test, alpha, B, strict, pc.method, nbr.join,
  debug=FALSE) {
  
  nodes = names(data)
  
  cache = new.env()
  mb = lapply(as.list(nodes), hybrid.pc.cached, data = data, nodes = nodes,
              alpha = alpha, B = B, whitelist = whitelist, blacklist = blacklist,
              test = test, pc.method = pc.method, cache = cache, debug = debug)
  rm(cache)
  names(mb) = nodes
  
  # check neighbourhood sets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, debug = debug,
                   filter = nbr.join)
  
  return(mb)
  
}#HYBRID.PC.CACHED.GLOBAL

hybrid.pc.cached = function(t, data, nodes, whitelist, blacklist, test, alpha, B,
  pc.method, cache, debug = FALSE) {

  # Initialise cache if needed
  if(is.null(cache$init)) {
    cache$init = TRUE
    cache$pcs.checked = rep(FALSE, length(nodes))
    names(cache$pcs.checked) = nodes
    cache$pcs0 = list()
    cache$pcs = list()
    cache$pvals = list()
    cache$dsep = list()
    for (node in nodes) {
      cache$pcs0[[node]] = vector()
      cache$pcs[[node]] = vector()
      cache$pvals[[node]] = vector()
      cache$dsep[[node]] = list()
    }#FOR
    if (debug)
      cat("- cache initialised.\n")
  }#THEN
  
  # 1. [PCS] Search parents and children superset
  if (!cache$pcs.checked[t])
    hybrid.pc.de.pcs.cached(t, data, nodes, alpha, B, whitelist, blacklist,
                            test, cache, debug)

  pcs = cache$pcs[[t]]
  dsep = cache$dsep[[t]]
  
  if (debug) {
    cat("- PCS of", t, "= '", pcs, "'.\n")
  }#THEN

  #optimisation : 0 or 1 node in PCS --> PC == PCS
  if(length(pcs) < 2)
    return(list(nbr = pcs, mb = NULL))
  
  # Search my neighbours' PCS too, for the RSPS search
  for (node in setdiff(cache$pcs[[t]], names(cache$pcs.checked[cache$pcs.checked])))
    hybrid.pc.de.pcs.cached(node, data, nodes, alpha, B, whitelist, blacklist,
                            test, cache, debug)
  
  # 2. [RSPS] Search remaining spouses superset, those not already in PCS
  rsps = hybrid.pc.de.rsps.cached(t, data, nodes, alpha, B, test, cache, debug)
  
  #optimisation : 2 nodes in PC and no SP --> PC == PCS
  if(length(c(pcs, rsps)) < 3)
    return(list(nbr = pcs, mb = pcs))
  
  # 3. [PC] Get the Parents and Children from nodes within PCS and RSPS
  pc = hybrid.pc.nbr.search(t = t, data = data, nodes = c(t, pcs, rsps),
                            test = test, alpha = alpha, B = B,
                            whitelist = whitelist, blacklist = blacklist,
                            backtracking = NULL, debug = debug,
                            method = pc.method)
  
  # 4. [Neighbourhood OR] Detect and add false-negatives to PC, by checking if
  #     the target is present in potential neighbours' neighbours
  for (node in pcs[!pcs %in% pc]) {
    
    pcn = hybrid.pc.nbr.search(t = node, data = data, nodes = c(t, pcs, rsps),
                               test = test, alpha = alpha, B = B,
                               whitelist = whitelist, blacklist = blacklist,
                               backtracking = NULL, debug = debug,
                               method = pc.method)
    
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
    cat("  * PC =", pc, "\n")
  
  res = list(nbr = pc, mb = c(pcs, rsps))
  
  return(res)
  
}#HYBRID.PC.CACHED

hybrid.pc.de.pcs.cached = function(t, data, nodes, alpha, B, whitelist, blacklist,
                                   test, cache, debug = FALSE) {
  
  pcs0 = cache$pcs0[[t]]
  pcs = cache$pcs[[t]]
  pvals = cache$pvals[[t]]
  dsep = cache$dsep[[t]]
  
  cache$pcs.checked[t] = TRUE
  
  whitelisted = nodes[vapply(nodes,
                             function(x) { is.whitelisted(whitelist, c(t, x), either = TRUE) }, logical(1))]
  blacklisted = nodes[vapply(nodes,
                             function(x) { is.blacklisted(blacklist, c(t, x), both = TRUE) }, logical(1))]
  
  if (debug) {
    cat("----------------------------------------------------------------\n")
    cat("* learning the PCS of", t, ".\n")
    cat(" * nodes already in PCS0 (cache): '", cache$pcs0[[t]], "'.\n")
    cat(" * nodes already not in PCS0 (cache): '", setdiff(names(cache$pcs.checked[cache$pcs.checked]), c(t, cache$pcs0[[t]])), "'.\n")
  }#THEN
  
  # cache optimisation, don't 0-degree check nodes already processed
  nodes.to.check = names(cache$pcs.checked[!cache$pcs.checked])
  
  if (debug) {
    cat(" * nodes to be tested for inclusion: '", nodes.to.check, "'.\n")
  }#THEN
  
  nodes.to.check = setdiff(nodes.to.check, c(whitelisted, blacklisted))
  
  # Phase (I): add X if Dep(T,X) (0-degree d-separated nodes)
  for (x in nodes.to.check) {
    
    a = conditional.test(t, x, c(), data = data, test = test, B = B,
                         alpha = alpha, debug = debug)
    
    if (a <= alpha) {
      
      if (debug)
        cat("  @ node", x, "added to the parents and children superset\n")
      
      # keep pcs with p-values in order to order nodes in the next phase
      pcs0 = c(pcs0, x)
      pcs = c(pcs, x)
      pvals = c(pvals, a)
      
      # cache update
      cache$pcs0[[x]] = c(cache$pcs0[[x]], t)
      cache$pcs[[x]] = c(cache$pcs[[x]], t)
      cache$pvals[[x]] = c(cache$pvals[[x]], a)
      
    }#THEN
    
  }#FOR
  
  # cache update
  cache$pcs0[[t]] = pcs0
  cache$pcs[[t]] = pcs
  cache$pvals[[t]] = pvals
  
  # heuristic 1 : sort the PC candidates in decreasing p-value order
  # this way we are more prone to remove less correlated nodes first
  ord = order(pvals, decreasing = TRUE)
  pcs = pcs[ord]
  pvals = pvals[ord]
  
  nodes.to.check = pcs
  
  if (debug) {
    cat(" * nodes to be tested for exclusion: '", nodes.to.check, "'.\n")
  }#THEN
  
  pcs = union(whitelisted, pcs)
  pvals = c(rep(0, length(pcs) - length(pvals)), pvals)
  
  # heuristic 2 : sort the d-separating canditates in increasing p-value order
  # this way we are more prone to remove with highly correlated nodes first
  nodes.to.check.against = pcs[order(pvals)]
  
  # Phase (II): remove X if Ind(T,X|Y) (1-degree d-separated nodes)
  for (x in nodes.to.check) {
    
    for (y in setdiff(nodes.to.check.against, x)) {
      
      a = conditional.test(t, x, y, data = data, test = test, B = B,
                           alpha = alpha, debug = debug)
      
      x.ind = which(pcs == x)
      
      if (a > alpha) {
        
        nodes.to.check.against = setdiff(nodes.to.check.against, x)
        
        pcs = pcs[-x.ind]
        pvals = pvals[-x.ind]
        dsep[[x]] = y
        
        if (debug) {
          cat("  @ node", x, "removed from the parents and children superset\n")
        }#THEN
        
        break
        
      }#THEN
      else {
        
        if (a > pvals[x.ind]) {
          
          pvals[x.ind] = a
          
        }#THEN
        
      }#ELSE
    }#FOR
  }#FOR
  
  # cache update
  cache$pcs[[t]] = pcs
  cache$pvals[[t]] = pvals
  cache$dsep[[t]] = dsep
  
  return(NULL)
  
}#HYBRID.PC.DE.PCS.CACHED

hybrid.pc.de.rsps.cached = function(t, data, nodes, alpha, B, test, cache,
                                    debug = FALSE) {
  
  pcs = cache$pcs[[t]]
  dsep = cache$dsep[[t]]
  rsps = character(0)
  
  if (debug) {
    cat("----------------------------------------------------------------\n")
    cat("* learning the RSPS of", t, ".\n")
    cat(" * PCS =", pcs, "\n")
  }#THEN
  
  for (x in pcs) {
    # optimisation : don't re-check nodes already in RSPS
    for (y in setdiff(cache$pcs[[x]], c(t, pcs, rsps))) {
        
      # optimisation : avoid irrelevant tests
      if (x %in% dsep[[y]])
        next
      
      a = conditional.test(t, y, c(dsep[[y]], x), data = data, test = test,
                           B = B, alpha = alpha, debug = debug)
      
      if (a <= alpha) {
        
        rsps = c(rsps, y)
        
        if (debug) {
          cat("  @ node", y, "added to the spouses superset\n")
        }#THEN
      
      }#THEN
      
    }#FOR
    
  }#FOR
  
  return(rsps)
  
}#HYBRID.PC.DE.RSPS.CACHED
