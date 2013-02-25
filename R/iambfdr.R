
iambfdr.global = function(x, whitelist, blacklist, test, alpha,
                          B, strict, debug = FALSE) {
  
  nodes = names(x)
  
  # 1. [Compute Markov Blankets]
  mb = lapply(as.list(nodes), iambfdr, data = x, nodes = nodes,
              alpha = alpha, B = B, whitelist = whitelist, blacklist = blacklist,
              test = test, debug = debug)
  names(mb) = nodes
  
  # check markov blankets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, mb = TRUE, debug = debug)
  
  # 2. [Compute Graph Structure]
  mb = lapply(as.list(nodes), neighbour, mb = mb, data = x, alpha = alpha,
              B = B, whitelist = whitelist, blacklist = blacklist, test = test,
              debug = debug)
  names(mb) = nodes
  
  # check neighbourhood sets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, debug = debug)
  
  return(mb)
  
}#IAMBFDR.GLOBAL

iambfdr.global.optimized = function(x, whitelist, blacklist,
                                    test, alpha, B, strict, debug = FALSE) {
  
  nodes = names(x)
  mb2 = mb = list()
  
  # 1. [Compute Markov Blankets]
  for (node in nodes) {
    
    backtracking = unlist(sapply(mb, function(x){ node %in% x  }))
    
    mb[[node]] = iambfdr(
      node, data = x, nodes = nodes, alpha = alpha, B = B,
      whitelist = whitelist, blacklist = blacklist,
      backtracking = backtracking, test = test, debug = debug)
    
  }#FOR
  
  # check markov blankets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, mb = TRUE, debug = debug)
  
  # 2. [Compute Graph Structure]
  for (node in nodes) {
    
    backtracking = unlist(sapply(mb2, function(x){ node %in% x[["nbr"]]  }))
    
    # save results in a copy of mb.
    mb2[[node]] = neighbour(
      node, mb = mb, data = x, alpha = alpha, B = B,
      whitelist = whitelist, blacklist = blacklist,
      backtracking = backtracking, test = test, debug = debug)
    
  }#FOR
  
  # update mb with the results of neighbour().
  mb = mb2
  
  # check neighbourhood sets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, debug = debug)
  
  return(mb)
  
}#IAMBFDR.GLOBAL.OPTIMIZED

iambfdr.global.cluster = function(x, cluster, whitelist,
                                  blacklist, test, alpha, B, strict, debug = FALSE) {
  
  nodes = names(x)
  
  # 1. [Compute Markov Blankets]
  mb = parLapply(cluster, as.list(nodes), iambfdr, data = x,
                 nodes = nodes, alpha = alpha, B = B, whitelist = whitelist,
                 blacklist = blacklist, test = test, debug = debug)
  names(mb) = nodes
  
  # check markov blankets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, mb = TRUE, debug = debug)
  
  # 2. [Compute Graph Structure]
  mb = parLapply(cluster, as.list(nodes), neighbour, mb = mb, data = x,
                 alpha = alpha, B = B, whitelist = whitelist,
                 blacklist = blacklist, test = test, debug = debug)
  names(mb) = nodes
  
  # check neighbourhood sets for consistency.
  mb = bn.recovery(mb, nodes = nodes, strict = strict, debug = debug)
  
  return(mb)
  
}#IAMBFDR.GLOBAL.CLUSTER

iambfdr = function(x, data, nodes, alpha, B, whitelist, blacklist,
                   start = character(0), backtracking = NULL, test, debug = FALSE) {
  
  whitelisted = nodes[sapply(
    nodes, function(y) { is.whitelisted(whitelist, c(x, y), either = TRUE) })]
  known.good = known.bad = character(0)
  mb = start
  
  if (debug) {
    
    cat("----------------------------------------------------------------\n")
    cat("* learning the markov blanket of", x, ".\n")
    
    if (length(start) > 0)
      cat("* initial set includes '", mb, "'.\n")
    
  }#THEN
  
  # whitelisted nodes are included by default (if there's a direct arc
  # between them of course they are in each other's markov blanket).
  # arc direction is irrelevant here.
  mb = union(mb, whitelisted)
  
  # blacklist is not checked, not all nodes in a markov blanket must be
  # neighbours.
  
  # use backtracking for a further screening of the nodes to be checked.
  if (!is.null(backtracking)) {
    
    # nodes whose markov blanket includes this node are included, because
    # X \in MB(Y) <=> Y \in MB(X)
    known.good = names(backtracking[backtracking])
    
    # and vice versa X \not\in MB(Y) <=> Y \not\in MB(X)
    known.bad = names(backtracking[!backtracking])
    
    if (debug) {
      
      cat("    * known good (backtracking): '", known.good, "'.\n")
      cat("    * known bad (backtracking): '", known.bad, "'.\n")
      cat("    * nodes still to be tested for inclusion: '", nodes, "'.\n")
      
    }#THEN
    
  }#THEN
  
  # known good nodes are included by default
  mb = union(mb, known.good)
  
  q = length(nodes)
  nodes = setdiff(nodes, x)
  just.added.node = NULL
  just.removed.node = NULL
  
  fdr.factor = (q-1) * sum(1/(1:(q-1)))
  
  culprit = character(0)
  loop.counter = 1
  state = vector(5 * length(nodes), mode = "list")
  
  repeat {
    
    state[[loop.counter]] = mb
    next.round = FALSE
    
    # get an association measure for each of the available nodes.
    association = sapply(nodes, function(node) {
      conditional.test(x, node, sx = setdiff(mb, node), test = test,
                       data = data, B = B, alpha = alpha)})
    
    pvalues.order = order(association)
    association = association[pvalues.order]
    nodes = nodes[pvalues.order]
    
    # 1 - Check nodes for exclusion
    for (i in (q-1):1) {
      
      pvalue = association[i]
      node = nodes[i]
      
      # candidate nodes for removal from the Markov Blanket
      # whitelisted and known.good nodes cannot be removed
      # a node added in previous round cannot be removed neither
      if (node %in% mb && !(node %in% c(whitelisted, known.good, just.added.node))) {
        if (pvalue * fdr.factor / i > alpha) {
          
          if (debug)
            cat("  @", node, "removed from the markov blanket\n")
          
          mb = setdiff(mb, node)
          
          next.round = TRUE
          just.removed.node = node
          just.added.node = NULL
          break
          
        }#THEN
      }#THEN
    }#FOR
    
    # 2 - Check nodes for inclusion
    if (!next.round) {
      for (i in 1:(q-1)) {
        
        pvalue = association[i]
        node = names(pvalue)
        
        if (!(node %in% c(mb, culprit))) {
          if (pvalue * fdr.factor / i <= alpha) {
            
            if (debug)
              cat("  @", node, "added to the markov blanket\n")
            
            mb = c(mb, node)
            
            next.round = TRUE
            just.removed.node = NULL
            just.added.node = node
            break
            
          }#THEN
        }#THEN
      }#FOR
    }#THEN
    
    # 3 - Repeat the process
    if (next.round) {
      
      # Prevent infinite loops
      for(prev.mb in state[1:loop.counter]) {
        if(setequal(mb, prev.mb)) {
          mb = setdiff(mb, c(just.added.node, just.removed.node))
          culprit = c(culprit, just.added.node, just.removed.node)
          loop.counter = loop.counter - 1
        }#THEN
      }#THEN
      
      loop.counter = loop.counter + 1
      
      next
      
    }#THEN
    
    break
    
  }#REPEAT
  
  return(mb)
  
}#IAMBFDR
