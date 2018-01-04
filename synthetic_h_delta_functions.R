delta.functions <- list(
  list(
    name = "incgraph",
    init=function(amnt.nodes, net.df) {
      incgraph::new.incgraph.network(amnt.nodes, net.df)
    },
    calc.delta=function(state, edge, change.edge) {
      tf <- edge[,1]
      tg <- edge[,2]
      
      incgraph::flip(state, tf, tg)
      delta <- incgraph::calculate.delta(state, tf, tg)
      if (!change.edge) {
        incgraph::flip(state, tf, tg)
      }
      
      list(state=state, delta=delta)
    }
  ),
  list(
    name = "orca",
    init=function(amnt.nodes, net.df) {
      incgraph::new.incgraph.network(amnt.nodes, net.df)
    },
    calc.delta=function(state, edge, change.edge) {
      tf <- edge[,1]
      tg <- edge[,2]
      
      incgraph::flip(state, tf, tg)
      delta <- incgraph::calculate.orbit.counts(state)
      if (!change.edge) {
        incgraph::flip(state, tf, tg)
      }
      
      list(state=state, delta=delta)
    }
  )
)