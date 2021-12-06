# This code is a modified version of a function from the castor package originally written by Stilianos Louca (https://cran.r-project.org/web/packages/castor/)
# Maximum parsimony ancestral state reconstruction for discrete traits.
# Modification of Sankoff algorithm for reconstructing discrete ancestral states (Weighted Small Parsimony Problem)
# Sankoff's algorithm allows the inclusion of a cost matrix:
#  	transition_costs[i,j] is the cost of transitioning i-->j (ignoring edge length)
# 	If transition_costs is "all_equal", then all transitions are penalized equally (same as if transition_costs[i,j] = 1-delta_{ij})
# 	If transition_costs is "sequential", then only single-step transitions (i-->i+1) are allowed, and all are penalized equally
# 	If transition_costs is "proportional", then all transition are allowed but they are penalized proportionally to the number of steps.
#  The modification of this function is that optionally, edge lengths can be used to weight the transition costs:
#  	Longer edges imply smaller transition costs between states
#  	Specifically, the cost of transitioning is transition_cost[i,j]/(edge_length^edge_exponent)
# 	where edge_exponent can be e.g. 0 (Sankoff's original algorithm), 1 (linear weighting) or 0.5 (square-root weighting, corresponding to a Brownian motion)
#  Requirements:
# 	Tree can be multifurcating, and can also include nodes with a single child
#	Tree must be rooted.
# 	If (edge_exponent>0) then: All edges must have length > 0
#  For a description of the original Sankoff algorithm, see: 
# 	http://telliott99.blogspot.ca/2010/03/fitch-and-sankoff-algorithms-for.html
# 	(page 11) https://cs.brown.edu/courses/csci1950-z/slides/CSCI1950ZFall09_Lecture2.pdf
#  The function returns ancestral state probabilities as a (non-flattened) NumericMatrix of size Nnodes x Nstates.
asr_max_parsimony = function(	tree, 
                              tip_states, 			# integer vector of size Ntips
                              Nstates				= NULL, 
                              transition_costs	= "all_equal", 
                              edge_exponent		= 0,
                              weight_by_scenarios	= TRUE,
                              check_input			= TRUE){
  Ntips  = length(tree$tip.label)
  Nedges = nrow(tree$edge)
  
  # basic error checking
  if(length(tip_states)!=Ntips) stop(sprintf("ERROR: Length of tip_states (%d) is not the same as the number of tips in the tree (%d)",length(tip_states),Ntips));
  if(!is.numeric(tip_states)) stop(sprintf("ERROR: tip_states must be integers"))	
  if(is.null(Nstates)) Nstates = max(tip_states);
  if(check_input){
    min_tip_state = min(tip_states)
    max_tip_state = max(tip_states)
    if((min_tip_state<1) || (max_tip_state>Nstates)) stop(sprintf("ERROR: tip_states must be integers between 1 and %d, but found values between %d and %d",Nstates,min_tip_state,max_tip_state))
    if((!is.null(names(tip_states))) && any(names(tip_states)!=tree$tip.label)) stop("ERROR: Names in tip_states and tip labels in tree don't match (must be in the same order).")
  }
  
  # construct transition_costs matrix if needed
  if(is.character(transition_costs)){
    if(transition_costs=="all_equal"){
      # all transitions penalized equally
      transition_costs = matrix(1, nrow=Nstates, ncol=Nstates)
    }else if(transition_costs=="sequential"){
      # only single-step transitions are allowed, and all are penalized equally
      transition_costs = matrix(Inf, nrow=Nstates, ncol=Nstates)
      for(i in 1:Nstates){
        if(i<Nstates) transition_costs[i,i+1] = 1
        if(i>1) transition_costs[i,i-1] 	  = 1
      }
    }else if(transition_costs=="proportional"){
      # all transitions are allowed, but penalized proportional to the number of steps
      transition_costs = matrix(0, nrow=Nstates, ncol=Nstates)
      for(i in 1:Nstates){
        transition_costs[i,] = sapply(1:Nstates, function(j) abs(j-i))
      }
    }else if(transition_costs=="exponential"){
      # all transitions are allowed, but penalized exponentially to the number of steps
      transition_costs = matrix(0, nrow=Nstates, ncol=Nstates)
      for(i in 1:Nstates){
        transition_costs[i,] = sapply(1:Nstates, function(j) exp(abs(j-i)))
      }
    }else{
      stop(sprintf("ERROR: Uknown transition_costs '%s'",transition_costs));
    }
    diag(transition_costs) = 0; # no cost for staying in the same state
  }else{
    if(nrow(transition_costs)!=Nstates || ncol(transition_costs)!=Nstates) stop(sprintf("ERROR: Transition costs has wrong size (%d x %d), expected %d x %d",nrow(transition_costs),ncol(transition_costs),Nstates,Nstates));
    if(check_input){
      if(any(transition_costs<0)) stop(sprintf("ERROR: Some transition costs are negative (found value %g)",min(transition_costs)))
      if(any(diag(transition_costs)!=0)) stop(sprintf("ERROR: The diagonal of transition_costs includes non-zero values, which makes no sense in a maximum-parsimony model"))
      if((edge_exponent!=0) && (!is.null(tree$edge.length)) && (any(tree$edge.length==0))) stop(sprintf("ERROR: edge_exponent is non-zero, but some edges in the tree have zero length"))
    }
  }
  
  results = WMPR_ASR_CPP(	Ntips			 						= Ntips,
                          Nnodes			 						= tree$Nnode,
                          Nedges			 						= Nedges,
                          Nstates			 						= Nstates,
                          tree_edge 		 						= as.vector(t(tree$edge)) - 1, # flatten in row-major format and adjust clade indices to 0-based
                          edge_length		 						= (if(is.null(tree$edge.length)) numeric() else tree$edge.length),
                          tip_states		 						= tip_states-1,
                          transition_costs 						= as.vector(t(transition_costs)),
                          branch_length_exponent 					= edge_exponent,
                          weight_posteriors_by_scenario_counts	= weight_by_scenarios,	# (INPUT) if true, then the posterior_probability of a state (in a specific node) is proportional to the number of scenarios in which the trait is at that state
                          verbose									= FALSE,
                          verbose_prefix							= "");
  
  scenario_counts <- matrix(results$scenario_counts, ncol=Nstates, byrow=TRUE)
  transition_counts <- matrix(results$transition_counts, byrow = TRUE, ncol = 2)
  
  isabovetip <- tree$edge[,1][match(1:Ntips, tree$edge[,2])]
  
  changes <- matrix(0, ncol = 1, nrow = Ntips + tree$Nnode)
  changes[(Ntips+1):nrow(changes),] <- rowSums(transition_counts)/rowSums(scenario_counts)
  changes[1:Ntips] <- (results$posterior_probabilities[isabovetip-Ntips,1] * (tip_states==2)
  ) + (results$posterior_probabilities[isabovetip-Ntips,2] * (tip_states==1))
  
  
  return(list(success 				= TRUE,
              ancestral_likelihoods 	= results$posterior_probabilities,
              scenario_counts			= scenario_counts,
              total_cost 				= results$best_root_cost,
              changes = changes));
}
