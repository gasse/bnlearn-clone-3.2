
# load suggested packages and initialize global variables.
.onLoad = function(lib, pkg) {

  # set the test/score counters.
  assign(".test.counter", 0, envir = .GlobalEnv)
  assign(".test.counter.permut", 0, envir = .GlobalEnv)
  assign(".score.counter", 0, envir = .GlobalEnv)

  # load the shared library.
  library.dynam("bnlearn", package = pkg, lib.loc = lib)

}#.ONLOAD

