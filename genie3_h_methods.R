
methods <- list(
  incgraph = list(
    name = "IncGraph",
    start = function(state) {},
    delta = function(state, ei, ej, start) {
      delta.calc <- incgraph::calculate.delta(state, ei, ej)
      delta.calc$rem - delta.calc$add
    }
  ), 
  orca = list(
    name = "Orca",
    start = function(state) incgraph::calculate.orbit.counts(state),
    delta = function(state, ei, ej, start) {
      new.gr <- incgraph::orca.halfdelta(state, ei, ej)
      new.gr - start
    }
  )
)