// This code was originally written by Stilianos Louca has been modified from https://cran.r-project.org/web/packages/castor/

#include <new>
#include <limits>
#include <cmath>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <complex>
#include <algorithm>
#include <Rcpp.h>
#include <time.h>
#include <ctype.h>
#include <sys/time.h>

#ifndef INFTY_D
#define INFTY_D numeric_limits<double>::infinity()
#endif

#ifndef NAN_D
#define NAN_D std::numeric_limits<double>::quiet_NaN()
#endif

#ifndef RELATIVE_EPSILON
#define RELATIVE_EPSILON 1e-10
#endif

#ifndef STRANDOM_EPSILON 
#define STRANDOM_EPSILON 1e-30
#endif


typedef std::complex<double> cdouble;
typedef std::vector<double> dvector;
typedef std::vector<long> lvector;

using namespace Rcpp;
using namespace std;


inline double string2Double(const string &number){
  return strtod(number.c_str(), NULL);
}

inline double string2Long(const string &number){
  return strtol(number.c_str(), NULL, 0);
}

template<class TYPE> 
string makeString(const TYPE &data){
  ostringstream stream;
  stream << data;
  return stream.str();
}


// Formated string creation
string vstringprintf(const char *format, va_list args){
  va_list temp;
  va_copy(temp, args);
  char *buffer = new char[vsnprintf(NULL, 0, format, temp) + 1];
  va_end(temp);
  
  vsprintf(buffer, format, args);
  string s(buffer);
  delete [] buffer;
  return s;
}

// Formated string creation
string stringprintf(const char *format, ...){
  string s;
  va_list args;
  va_start(args, format);
  s = vstringprintf(format, args);
  va_end(args);
  return s;
}



string trim_whitespace(const std::string &haystack){
  long right = haystack.length()-1;
  long left = 0;
  while(((haystack[right]==' ') || (haystack[right]=='\t') || (haystack[right]=='\n')) && (right>=0)){
    --right;
  }
  while(((haystack[left]==' ') || (haystack[left]=='\t') || (haystack[left]=='\n')) && (left<right)){
    ++left;
  }
  return haystack.substr(left,right-left+1);
}


template<class ARRAY_TYPE>
void multiply_array_with_scalar(ARRAY_TYPE &X, const double scalar){
  for(long i=0; i<X.size(); ++i) X[i] *= scalar;
}




template<class ARRAY_TYPE>
void reverse_array(ARRAY_TYPE &X){
  const long N = X.size();
  long scratch;
  for(long n=0; n<(N/2.0); ++n){
    scratch = X[n];
    X[n] = X[N-n-1];
    X[N-n-1] = scratch;
  }
}



template<class ARRAY_TYPE>
double array_min(const ARRAY_TYPE &X){
  const long N = X.size();
  if(N==0) return NAN_D;
  double minX = X[0];
  for(long n=0; n<N; ++n){
    if(X[n]<minX) minX = X[n];
  }
  return minX;
}


template<class ARRAY_TYPE>
double array_min(const ARRAY_TYPE &X, const long start, const long end){
  if(end<start) return NAN_D;
  double minX = X[start];
  for(long n=start; n<=end; ++n){
    if(X[n]<minX) minX = X[n];
  }
  return minX;
}


template<class ARRAY_TYPE>
inline double array_nonzero_min(const ARRAY_TYPE &X){
  const long N = X.size();
  double minX = NAN_D;
  for(long n=0; n<N; ++n){
    if((X[n]!=0) && (std::isnan(minX) || (X[n]<minX))) minX = X[n];
  }
  return minX;
}


template<class ARRAY_TYPE>
inline double array_max(const ARRAY_TYPE &X){
  const long N = X.size();
  if(N==0) return NAN_D;
  double maxX = X[0];
  for(long n=0; n<N; ++n){
    if(X[n]>maxX) maxX = X[n];
  }
  return maxX;
}



double array_min(const std::vector<double> &X, long start_index, long end_index){
  if(end_index<start_index) return NAN_D;
  double minX = X[start_index];
  for(long n=start_index; n<=end_index; ++n){
    if(X[n]<minX) minX = X[n];
  }
  return minX;
}


template<class ARRAY_TYPE>
inline long array_argmax(const ARRAY_TYPE &X){
  const long N = X.size();
  if(N==0) return -1;
  long max_n = 0;
  for(long n=1; n<N; ++n){
    if(X[n]>X[max_n]) max_n = n;
  }
  return max_n;
}


template<class TYPE>
inline long array_argmax(const long N, const TYPE X[]){
  if(N==0) return -1;
  long max_n = 0;
  for(long n=1; n<N; ++n){
    if(X[n]>X[max_n]) max_n = n;
  }
  return max_n;
}


template<class ARRAY_TYPE>
inline bool arrays_are_equal(const ARRAY_TYPE &A, const ARRAY_TYPE &B){
  if(A.size()!=B.size()) return false;
  for(long i=0; i<A.size(); ++i){
    if(A[i]!=B[i]) return false;
  }
  return true;
}


template<class TYPE>
inline TYPE vector_sum(const std::vector<TYPE> &values){
  TYPE S = 0;
  for(long i=0; i<values.size(); ++i) S += values[i];
  return S;
}


template<class TYPE>
inline TYPE vector_sum(const std::vector<TYPE> &values, const lvector &indices){
  TYPE S = 0;
  for(long i=0; i<indices.size(); ++i) S += values[indices[i]];
  return S;
}


inline long vector_sum(const std::vector<char> &values){
  long S = 0;
  for(long i=0; i<values.size(); ++i) S += values[i];
  return S;
}


template<class TYPE>
inline TYPE scalar_product(const std::vector<TYPE> &A, const std::vector<TYPE> &B){
  TYPE S = 0;
  for(long i=0; i<min(A.size(),B.size()); ++i) S += A[i] * B[i];
  return S;
}


double smallest_nonzero_step(const std::vector<double> &times){
  double S = INFTY_D;
  for(long i=0; i<times.size(); ++i){
    if(times[i+1]>times[i]){
      S = min(S, times[i+1]-times[i]);
    }
  }
  return S;
}


inline double vector_mean(const double values[], const long N){
  double S = 0;
  for(long i=0; i<N; ++i) S += values[i];
  return (S/N);
}


inline double vector_mean(const std::vector<double> &values){
  double S = 0;
  for(long i=0; i<values.size(); ++i) S += values[i];
  return (S/values.size());
}


inline double vector_abs_mean(const std::vector<double> &values){
  double S = 0;
  for(long i=0; i<values.size(); ++i) S += abs(values[i]);
  return (S/values.size());
}



inline double vector_mean(const std::vector<double> &values, const long first, const long last){
  double S = 0;
  for(long i=first; i<=last; ++i) S += values[i];
  return (S/(last-first+1.0));
}


template<class TYPE>
inline TYPE vector_min(const std::vector<TYPE> &values){
  TYPE A = values[0];
  for(long i=0; i<values.size(); ++i) A = (A>values[i] ? values[i] : A);
  return A;
}

template<class TYPE>
inline TYPE vector_max(const std::vector<TYPE> &values){
  TYPE A = values[0];
  for(long i=0; i<values.size(); ++i) A = (A<values[i] ? values[i] : A);
  return A;
}

inline double vector_max_abs(const std::vector<double> &values){
  double A = 0;
  for(long i=0; i<values.size(); ++i) A = max(A, abs(values[i]));
  return A;
}

// determine the maximum modulus of any non-diagonal entries in a 2D matrix
// The matrix should be of size NR x NC, in row-major format
inline double matrix_nondiagonal_max_abs(const std::vector<double> &matrix, const long NR, const long NC){
  double A = 0;
  for(long r=0; r<NR; ++r){
    for(long c=0; c<NC; ++c){
      if(r==c) continue;
      A = max(A, abs(matrix[r*NC+c]));
    }
  }
  return A;
}


// calculate sum of a single row in a 2D matrix of size NR x NC
// matrix must be provided in row-major format, i.e. matrix[r*NC+c] is the entry in row r & column c
inline double row_sum(const std::vector<double> &matrix, const long NC, const long row){
  double S = 0;
  for(long c=0; c<NC; ++c){
    S += matrix[row*NC + c];
  }
  return S;
}


// make sure no entry is negative
void make_vector_positive(std::vector<double> &values){
  for(long i=0; i<values.size(); ++i) values[i] = max(0.0, values[i]);
}

// replace negative entries
void replace_negatives(std::vector<double> &values, const double replacement){
  for(long i=0; i<values.size(); ++i){
    if(values[i]<0) values[i] = replacement;
  }
}

// replace non-strictly positive entries
void replace_non_positives(std::vector<double> &values, const double replacement){
  for(long i=0; i<values.size(); ++i){
    if(values[i]<=0) values[i] = replacement;
  }
}


// make sure entries in a vector are within the specified limits [min_value:max_value]
void cap_values(const double 		min_value,
                const double 		max_value,
                std::vector<double> &values){ // (INPUT/OUTPUT) the vector to be modified in-situ
  for(long i=0; i<values.size(); ++i){
    values[i] = max(min_value, min(max_value, values[i]));
  }
}


// extract a specific row from a 2D matrix in row-major format
template<class TYPE>
void extract_row(const std::vector<TYPE> &matrix, const long NC, const long row, std::vector<TYPE> &extracted_row){
  extracted_row.resize(NC);
  for(long c=0; c<NC; ++c){
    extracted_row[c] = matrix[row*NC+c];
  }
}


// calculate the best (lowest) cost of any transition parent-->child, given a particular parent state and a particular child cost table (extracted from master_cost_table)
// this function assumes that the cost table for the child has already been calculated (hence, you should move tips-->root)
double aux_get_cost_of_parent_state_transitioning_to_one_child(	const long					Nstates,
                                                                const long					parent_state, 
                                                                const long 					edge,
                                                                const double				edge_weight,
                                                                const long 					child,
                                                                const std::vector<double>			transition_costs, 		// (INPUT) 2D array of size Nstates x Nstates (in row-major format)
                                                                const std::vector<double>	&master_cost_table,		// (INPUT) 2D array of size (Ntips+Nnodes) x Nstates (in row-major format)
                                                                std::vector<double>			&scratch_space,			// temporary space for intermediate operations
                                                                std::vector<long>			&master_transitions,				// (INPUT/OUTPUT) 1D array (preferably reserved up to size Nnodes*Nstates*Nstates)
                                                                std::vector<long>			&edge_and_state2first_transition,	// (INPUT/OUTPUT) 1D array of size Nedges*Nstates.
                                                                std::vector<long>			&edge_and_state2last_transition){	// (INPUT/OUTPUT) 1D array of size Nedges*Nstates.
  std::vector<double> &choice_costs = scratch_space;
  choice_costs.resize(Nstates);
  for(long state=0; state<Nstates; ++state){
    choice_costs[state] = transition_costs[parent_state*Nstates + state]*edge_weight + master_cost_table[child*Nstates + state];
  }
  const double best_cost = array_min(choice_costs);
  edge_and_state2first_transition[edge*Nstates + parent_state] = master_transitions.size();
  for(long transition=0; transition<Nstates; ++transition){
    if(abs(choice_costs[transition]-best_cost)<=RELATIVE_EPSILON*best_cost){
      master_transitions.push_back(transition);
    }
  }
  edge_and_state2last_transition[edge*Nstates + parent_state]	= master_transitions.size()-1;	
  return best_cost;
}




// calculate the cost of a particular state in a particular node, best on the best (lowest) costs of transitions to each of the children
// this function assumes that the cost table of each child has already been calculated (hence, you should move tips-->root)
double aux_get_cost_of_parent_state_transitioning_to_all_children(	const long					Nstates,
                                                                   const long 					node, 							// (INPUT) integer in 0:(Nnodes-1)
                                                                   const long 					parent_state,					// (INPUT) integer in 0:(Nstates-1)
                                                                   const double				branch_length_exponent,			// (INPUT) non-negative number
                                                                   const std::vector<double>			transition_costs,	 			// (INPUT) 2D array of size Nstates x Nstates (in row-major format)
                                                                   const std::vector<double>	&master_cost_table,				// (INPUT) 2D array of size (Ntips+Nnodes) x Nstates (in row-major format)
                                                                   const std::vector<long>			&tree_edge,						// (INPUT) 2D array of size Nedges x 2 (in row-major format), in similar format as tree$edge in R "phylo" trees.
                                                                   const std::vector<double>			&edge_length,				// (INPUT) 1D array of size Nedges
                                                                   const std::vector<long>		&traversal_node2first_edge,		// (INPUT) 1D array of size Nnodes, with values in 0:(Nedges-1). Generated using the function get_tree_traversal_root_to_tips().
                                                                   const std::vector<long>		&traversal_node2last_edge,		// (INPUT) 1D array of size Nnodes, with values in 0:(Nedges-1). Generated using the function get_tree_traversal_root_to_tips().
                                                                   const std::vector<long>		&traversal_edges,				// (INPUT) 1D array of size Nedges, with values in 0:(Nedges-1). Generated using the function get_tree_traversal_root_to_tips().
                                                                   std::vector<double>			&scratch_space,						// temporary space for intermediate operations
                                                                   std::vector<long>			&master_transitions,				// (INPUT/OUTPUT) 1D array (preferably reserved up to size Nnodes*Nstates*Nstates)
                                                                   std::vector<long>			&edge_and_state2first_transition,	// (INPUT/OUTPUT) 1D array of size Nedges*Nstates.
                                                                   std::vector<long>			&edge_and_state2last_transition){	// (INPUT/OUTPUT) 1D array of size Nedges*Nstates.
  double S = 0;
  double edge_weight;
  long edge, child;
  for(long ei=traversal_node2first_edge[node]; ei<=traversal_node2last_edge[node]; ++ei){
    edge  = traversal_edges[ei];
    child = tree_edge[edge*2+1];
    edge_weight = (branch_length_exponent==0 ? 1.0 : 1/pow((edge_length.size()==0 ? 1.0 : edge_length[edge]),branch_length_exponent));
    S += aux_get_cost_of_parent_state_transitioning_to_one_child(	Nstates,
                                                                  parent_state, 
                                                                  edge, 
                                                                  edge_weight, 
                                                                  child,
                                                                  transition_costs,
                                                                  master_cost_table,
                                                                  scratch_space,
                                                                  master_transitions,
                                                                  edge_and_state2first_transition,
                                                                  edge_and_state2last_transition);
  }
  return S;
}




// determine root of a tree
// Assuming that the tree is connected (when edge directions are ignored), this function will return -1 if the tree is not properly rooted
// Hence, this function can also be used to check if the tree is properly rooted (provided that it is connected)
template<class ARRAY_TYPE>
long get_root_clade(const long			Ntips,
                    const long 			Nnodes,
                    const long			Nedges,
                    const ARRAY_TYPE	&tree_edge){			// (INPUT) 2D array (in row-major format) of size Nedges x 2
  const long Nclades = Ntips+Nnodes;
  std::vector<long> Nparents_per_clade(Nclades,0);
  for(long edge=0; edge<Nedges; ++edge){
    Nparents_per_clade[tree_edge[edge*2+1]] += 1;
  }
  long root = -1;
  for(long c=0; c<Nclades; ++c){
    if(Nparents_per_clade[c]>1) return -1; // found a clade with multiple parents, which cannot be in a rooted tree
    if(Nparents_per_clade[c]==0){
      // found clade with no parents, so this may be root
      if(root>=0) return -1; // multiple roots found, which cannot be
      root = c;
    }
  }
  return root;
}


// Calculate lookup tables mapping nodes to their outgoing (children) edges
// Requirements:
//    The tree can be multifurcating, and can also include nodes with a single child
//    The tree can be rooted or unrooted (all information on edge direction is taken from the tree_edge[] table)
// Returned values:
//	  node2first_edge[p] will be an index pointing node p (p=0:(Nnodes-1)) to edges[]
//	  node2last_edge[p] will be an index pointing node p (p=0:(Nnodes-1)) to edges[]
// 	  edges[] will be a list of edge indices (i.e. in 0:(Nedges-1)), such that edges[node2first_edge[p]],...,edges[node2last_edge[p]] is the set of edges leaving node p
template<class ARRAY_TYPE>
void get_node2edge_mappings(const long			Ntips,
                            const long 			Nnodes,
                            const long			Nedges,
                            const ARRAY_TYPE	&tree_edge, 		// (INPUT) 2D array (in row-major format) of size Nedges x 2
                            std::vector<long>	&node2first_edge,	// (OUTPUT) 1D array of size Nnodes, with values in 0:(Nedges-1), mapping nodes to their first outgoing edge.
                            std::vector<long>	&node2last_edge,	// (OUTPUT) 1D array of size Nnodes, with values in 0:(Nedges-1), mapping nodes to their last outgoing edge.
                            std::vector<long>	&edges){			// (OUTPUT) 1D array of size Nedges, with values in 0:(Nedges-1). Maps internal edge indices (i.e. as listed in node2first_edge[] and node2last_edge[]) to original edge indices.
  // Terminology in this function:
  // 	'node' runs from 0:(Nnodes-1)
  // 	'tip' runs from 0:(Ntips-1)
  // 	'parent' and 'child' runs from 0:(Ntips+Nnodes-1)
  // 	'edge' runs from 0:(Nedges-1)
  // Recall that:
  // 	tree_edge[] is of size Nedge x 2 (flattened in row-major-format), with entries in 0:(Ntips+Nnodes-1)
  
  edges.resize(Nedges);
  node2first_edge.resize(Nnodes);
  node2last_edge.resize(Nnodes);
  
  // determine number of children/edges per parent
  // child_count_per_node[n] will be the number of direct children of node n (n=0:(Nnodes-1))
  std::vector<long> child_count_per_node(Nnodes, 0);
  for(long e=0; e<Nedges; ++e){
    child_count_per_node[tree_edge[e*2+0] - Ntips] += 1;
  }
  // collect children per parent
  node2first_edge[0] = 0;
  node2last_edge[0]  = node2first_edge[0]+child_count_per_node[0] - 1;
  if(Nnodes>1){
    for(long n=1; n<Nnodes; ++n){
      node2first_edge[n] = node2last_edge[n-1]+1;
      node2last_edge[n]  = node2first_edge[n]+child_count_per_node[n] - 1;
    }
  }
  for(long e=0, node; e<Nedges; ++e){
    node = tree_edge[e*2+0] - Ntips;
    edges[node2first_edge[node]+child_count_per_node[node]-1] = e;
    child_count_per_node[node] -= 1;
  }
}



// Returns a list of all nodes (and optionally tips) of a tree, such that each node appears prior to its children, and such that nodes closer to the root (in terms of branching counts) appear first
// Also returns a list mapping nodes to their outgoing (children) edges (e.g. as listed in tree_edge)
// Nodes and tips are explored in a breadth-first-search order, from root to tips.
// the tree can be multifurcating, and can also include nodes with a single child
// root must specify the root index in the tree (typically root = Ntips)
// Returned values:
//	queue: A 1D array of integers in 0:(Ntips+Nnodes-1) if include_tips==true, or (Ntips):(Ntips+Nnodes-1) if include_tips==false
//	node2first_edge[p] will be an index pointing node p (p=0:(Nnodes-1)) to edges[]
//	node2last_edge[p] will be an index pointing node p (p=0:(Nnodes-1)) to edges[]
// 	edges[] will be a list of edge indices (i.e. in 0:(Nedges-1)), such that edges[node2first_edge[p]],...,edges[node2last_edge[p]] is the set of edges leaving node p
template<class ARRAY_TYPE>
void get_tree_traversal_root_to_tips(	const long			Ntips,
                                      const long 			Nnodes,
                                      const long			Nedges,
                                      const long 			root, 							// (INPUT) index of root node, i.e. an integer in 0:(Ntips+Nnodes-1)
                                      const ARRAY_TYPE	&tree_edge, 					// (INPUT) 2D array (in row-major format) of size Nedges x 2
                                      const bool			include_tips,					// (INPUT) if true, then tips are included in the returned queue[]. This does not affect the returned arrays node2first_edge[], node2last_edge[], edges[].
                                      const bool			precalculated_edge_mappings,	// (INPUT) if true, then the edge mapping tables node2first_edge[], node2last_edge[] and edges[] are taken as is. Otherwise, they are calculated from scratch.
                                      std::vector<long>	&queue,							// (OUTPUT) 1D array of size Nnodes if include_tips==false, or size (Ntips+Nnodes) if include_tips==true.
                                      std::vector<long>	&node2first_edge,				// (INPUT/OUTPUT) 1D array of size Nnodes, with values in 0:(Nedges-1), mapping nodes to their first outgoing edge. Either pre-calculated, or to be calculated by this function.
                                      std::vector<long>	&node2last_edge,				// (INPUT/OUTPUT) 1D array of size Nnodes, with values in 0:(Nedges-1), mapping nodes to their last outgoing edge. Either pre-calculated, or to be calculated by this function.
                                      std::vector<long>	&edges,							// (INPUT/OUTPUT) 1D array of size Nedges, with values in 0:(Nedges-1). Maps internal edge indices (i.e. as listed in node2first_edge[] and node2last_edge[]) to original edge indices. Either pre-calculated, or to be calculated by this function.
                                      const bool			verbose,
                                      const string		&verbose_prefix){
  // get node-->edge mappings if needed
  if(!precalculated_edge_mappings){
    get_node2edge_mappings(	Ntips,
                            Nnodes,
                            Nedges,
                            tree_edge,
                            node2first_edge,
                            node2last_edge,
                            edges);
  }
  
  // fill queue from root to tips
  long child,node;
  queue.clear();
  queue.reserve(include_tips ? Ntips+Nnodes : Nnodes);
  queue.push_back(root);
  long queue_pointer = 0;
  while(queue_pointer<queue.size()){
    node = queue[queue_pointer] - Ntips;
    queue_pointer += 1;
    if(node<0) continue; // queue[queue_pointer] was actually a tip, not an internal node
    if(node2first_edge[node]>node2last_edge[node]){
      // this should not happen (every node needs to have at least one child)
      if(verbose) Rcpp::Rcout << verbose_prefix << "WARNING: Node " << node << " has no children\n";
      continue;
    }
    for(long ei=node2first_edge[node]; ei<=node2last_edge[node]; ++ei){
      child = tree_edge[edges[ei]*2+1];
      if((!include_tips) && (child<Ntips)) continue; // this child is a tip, so skip as requested
      // append child to queue
      queue.push_back(child);
    }
  }
}







// Weighted maximum parsimony ansestral state reconstruction for discrete traits.
// Modification of Sankoff algorithm for reconstructing discrete ancestral states (Weighted Small Parsimony Problem)
// Sankoff's algorithm allows the inclusion of a cost matrix: transition_costs[i,j] is the cost of transitioning i-->j (ignoring edge length)
// The modification of this function is that optionally, edge lengths can be used to weight the transition costs:
// 	Longer edges imply smaller transition costs between states
// 	Specifically, the cost of transitioning is transition_cost[i,j]/(edge_length^branch_length_exponent)
//	where branch_length_exponent can be e.g. 0 (Sankoff's original algorithm), 1 (linear weighting) or 0.5 (square-root weighting, corresponding to a Brownian motion)
// Requirements:
//	Tree can be multifurcating, and can also include nodes with a single child
//	If (branch_length_exponent!=0) then: All branches must have length > 0
// For a description of the original Sankoff algorithm, see: 
//	http://telliott99.blogspot.ca/2010/03/fitch-and-sankoff-algorithms-for.html
//	(page 11) https://cs.brown.edu/courses/csci1950-z/slides/CSCI1950ZFall09_Lecture2.pdf
// The function returns a (non-flattened) NumericMatrix of size Nnodes x Nstates.
//
// Attention: Be carefull to use the C++ style indexing (0-based) when passing index-variables or index arrays to this function.
// For example, root must be a 0-based index, and tree_edge[] must have values in 0:(Ntips+Nnodes-1) instead of 1:(Ntips+Nnodes)
// [[Rcpp::export]]
Rcpp::List WMPR_ASR_CPP(const long					Ntips,
                        const long 					Nnodes,
                        const long					Nedges,
                        const long					Nstates, 				// (INPUT) number of possible states
                        const std::vector<long>		&tree_edge, 			// (INPUT) 2D array of size Nedges x 2 (in row-major format), in similar format as tree$edge in R "phylo" trees. This array holds the topology of the tree (apart from branch lengths).
                        const std::vector<double>	&edge_length,			// (INPUT) 1D array of size Nedges, synchronized with the rows of tree_edge[,], i.e. with edge_length[e] being the length of edge e. Can also be an empty vector (all edges have length 1.0).
                        const std::vector<long>		&tip_states, 			// (INPUT) 1D array of size Ntips, with values being in 0:(Nstates-1)
                        const std::vector<double>	&transition_costs,	 	// (INPUT) 2D array of size Nstates x Nstates (in row-major format), with transition_costs[i,j] being the cost of transition i-->j. Normally transition_cost[i,i]=0 for all i. Some transitions may be vorbitten, in which case the transition cost should be set to infinity (INFTY_D).
                        const double 				branch_length_exponent, // (INPUT) exponent for weighting transition costs by branch length. To ignore branch lengths (i.e. to obtain the non-weighted MPR algorithm), set this to 0.
                        bool						weight_posteriors_by_scenario_counts,	// (INPUT) if true, then the posterior_probability of a state (in a specific node) is proportional to the number of scenarios in which the trait is at that state
                        bool						verbose,
                        const std::string			&verbose_prefix){
  // Terminology in this function:
  // 	'node' runs from 0:(Nnodes-1)
  // 	'tip' runs from 0:(Ntips-1)
  // 	'parent' and 'child' runs from 0:(Ntips+Nnodes-1)
  // 	'edge' runs from 0:(Nedges-1)
  // 	'state' runs from 0:(Nstates-1)
  const long Nclades = Ntips+Nnodes;
  long node, state, parent, transition, child, edge;
  std::vector<double> scratch_space;
  
  // determine root
  const long root = get_root_clade(Ntips, Nnodes, Nedges, tree_edge);
  
  // create tree-access structures and determine order in which to traverse tree
  std::vector<long> traversal_queue_root2tips, traversal_node2first_edge, traversal_node2last_edge, traversal_edges;
  get_tree_traversal_root_to_tips(	Ntips,
                                   Nnodes,
                                   Nedges,
                                   root,
                                   tree_edge,
                                   false,
                                   false,
                                   traversal_queue_root2tips,
                                   traversal_node2first_edge,
                                   traversal_node2last_edge,
                                   traversal_edges,
                                   verbose,
                                   verbose_prefix);
  
  // get traversal route tips --> root
  std::vector<long> traversal_queue_tips2root = traversal_queue_root2tips;
  reverse_array(traversal_queue_tips2root); 
  
  
  // master_cost_table[,] should be a 2D numeric array of size (Ntips+Nnodes) x Nstates (in row-major-format)
  // the row master_cost_table[r,] is the cost table of tip or node r
  std::vector<double> master_cost_table(Nclades*Nstates, 0.0);
  
  // fill costs for tips (this is easy, since states are known)
  for(long tip=0; tip<Ntips; ++tip){
    for(state=0; state<Nstates; ++state){
      master_cost_table[tip*Nstates + state] = INFTY_D;
    }
    master_cost_table[tip*Nstates + tip_states[tip]] = 0.0;
  }
  
  // edge_and_state2first_transition[,] and edge_and_state2last_transition[,] contain indices mapping to master_transitions[]
  // for any edge e connecting parent node n to some child, and assuming n is at state s, the integers
  //	master_transitions[edge_and_state2first_transition[e,s]:edge_and_state2last_transition[e,s]]
  // are within in 1:Nstates, and are the optimal states to which node n switched during edge e.
  // Typically there will be only one "best transition" (i.e. edge_and_state2first_transition[e,s]==edge_and_state2last_transition[e,s]), 
  // but in case of multiple MPR solutions some edges may have multiple alternative best transitions (i.e. edge_and_state2first_transition[e,s] > edge_and_state2last_transition[e,s]
  
  // pre-allocate space at upper bound of possible need
  std::vector<long> master_transitions;
  master_transitions.reserve(Nnodes*Nstates*Nstates);
  std::vector<long> edge_and_state2first_transition(Nedges*Nstates);
  std::vector<long> edge_and_state2last_transition(Nedges*Nstates);
  
  
  // traverse tree (tips-->root) and build master_cost_table
  for(long parent_i=0; parent_i<Nnodes; ++parent_i){
    parent = traversal_queue_tips2root[parent_i];
    // calculate the cost associated with any state in this particular node
    for(state=0; state<Nstates; ++state){
      master_cost_table[parent*Nstates + state] = aux_get_cost_of_parent_state_transitioning_to_all_children(	Nstates,
                                                                                                              (parent - Ntips),
                                                                                                              state,
                                                                                                              branch_length_exponent,
                                                                                                              transition_costs,
                                                                                                              master_cost_table,
                                                                                                              tree_edge,
                                                                                                              edge_length,
                                                                                                              traversal_node2first_edge,
                                                                                                              traversal_node2last_edge,
                                                                                                              traversal_edges,
                                                                                                              scratch_space,
                                                                                                              master_transitions,
                                                                                                              edge_and_state2first_transition,
                                                                                                              edge_and_state2last_transition);
    }
    if(parent_i % 100 == 0) Rcpp::checkUserInterrupt();
  }
  
  
  // count number of scenarios (MPR solutions) implying each state in each node (based on lowest cost in the root, and then the associated transitions to the children)
  // scenario_count_per_node_and_state[n,s] will be the number of MPR solutions ("scenarios") in which node n is at state s
  // scenario_count_per_node_and_state[,] will be filled in the order root-->tips
  // See pages 18-19 in: https://cs.brown.edu/courses/csci1950-z/slides/CSCI1950ZFall09_Lecture2.pdf
  // Note: This should be floating point, not int, because in the latter case you risk integer overflow and thus the spontaneous introduction of negative values! This cost me 1 day of bug-hunting!
  std::vector<double> scenario_count_per_node_and_state(Nnodes*Nstates, 0.0); // scenario_count_per_node_and_state[Nnodes x Nstates] in row-major format
  std::vector<double> transition_count_per_node_and_state(Nnodes*Nstates, 0.0);
  
  const double best_root_cost = array_min(master_cost_table, root*Nstates, (root*Nstates+Nstates-1));
  for(state=0; state<Nstates; ++state){
    if(abs(master_cost_table[root*Nstates+state]-best_root_cost)<=RELATIVE_EPSILON*best_root_cost){
      scenario_count_per_node_and_state[(root-Ntips)*Nstates + state] = 1;
    }
  }
  
  for(long q=0; q<traversal_queue_root2tips.size(); ++q){
    parent 	= traversal_queue_root2tips[q];
    node	= parent-Ntips;
    for(long ei=traversal_node2first_edge[node]; ei<=traversal_node2last_edge[node]; ++ei){
      edge  = traversal_edges[ei];
      child = tree_edge[edge*2+1];
      if(child<Ntips) continue;
      for(state=0; state<Nstates; ++state){
        if(scenario_count_per_node_and_state[node*Nstates+state]>0){
          // examine all optimal transitions parent --> child, when parent is at this particular state
          double diffm = 1.0/edge_and_state2last_transition[edge*Nstates+state] - edge_and_state2first_transition[edge*Nstates+state] + 1;
          for(long transition_i=edge_and_state2first_transition[edge*Nstates+state]; transition_i<=edge_and_state2last_transition[edge*Nstates+state]; ++transition_i){
            transition = master_transitions[transition_i];
            // increment scenario_count for the corresponding state in this particular child
            scenario_count_per_node_and_state[(child-Ntips)*Nstates + transition] += scenario_count_per_node_and_state[node*Nstates + state];
            if (transition!=state){
              transition_count_per_node_and_state[(child-Ntips)*Nstates + transition] += scenario_count_per_node_and_state[node*Nstates + state];
            }
          }
        }
      }
    }
    if(q % 100 == 0) Rcpp::checkUserInterrupt();
  }
  
  
  // For a given tree, there may be multiple alternative scenarios (MPR solutions) for the ancestral states
  // based on the scenario count per node and per state, define posterior_probabilities for nodes
  double mass;
  NumericMatrix posterior_probabilities(Nnodes,Nstates); // this will be a 2D array of size Nnodes x Nstates. Note that in this case we're not flattening, for convenience, because we're returning this to R and there we like to have a non-flattened 2D matrix. 
  for(node=0; node<Nnodes; ++node){
    mass = 0;
    if(weight_posteriors_by_scenario_counts){
      // weight states proportional to the number of scenarios
      for(state=0, mass=0; state<Nstates; ++state) mass += scenario_count_per_node_and_state[node*Nstates + state];
      if(mass==0){
        //if(verbose) Rcout << verbose_prefix << "WARNING: Node " << node << " (clade " << (node+Ntips) << ") has max-parsimony mass = 0 (i.e. no predicted state). It's posterior probabilities will all be set to NaN\n";
        for(state=0; state<Nstates; ++state) posterior_probabilities(node,state) = NAN_D;
      }else{
        for(state=0; state<Nstates; ++state) posterior_probabilities(node,state) = scenario_count_per_node_and_state[node*Nstates + state]/mass;
      }
    }else{
      // all states with non-zero scenario count are weighted equally
      for(state=0, mass=0; state<Nstates; ++state) mass += (scenario_count_per_node_and_state[node*Nstates + state]>0 ? 1.0 : 0.0);
      for(state=0; state<Nstates; ++state) posterior_probabilities(node,state) = (scenario_count_per_node_and_state[node*Nstates + state]>0 ? 1.0 : 0.0)/mass;			
    }
    if(node % 100 == 0) Rcpp::checkUserInterrupt();
  }
  
  return Rcpp::List::create(	Rcpp::Named("posterior_probabilities") 	= Rcpp::wrap(posterior_probabilities),
                             Rcpp::Named("scenario_counts") 			= Rcpp::wrap(scenario_count_per_node_and_state),
                             Rcpp::Named("best_root_cost")			= best_root_cost,
                             Rcpp::Named("transition_counts")			= Rcpp::wrap(transition_count_per_node_and_state));
}
