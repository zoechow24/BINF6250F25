# BINF6250 - Project 06
# Introduction
Description of the project

# Pseudocode
Put pseudocode in this box:

```
For building the distance matrix:

- We already did SW implementation. However we probably don't need the traceback
    - Just need to return the max score for each pair of alignments
    - From there, we need to take our "similarity" scores, normalize them, and convert to probs
        - 1 - (Alignment of A/B) / (max((A,A),(B,B))) <--- A,A and B,B are essentially just match score * length(A) or length(B)
        - Compute this for all pairs. Maybe what we can do is do SW on all the alignments, then use np to vectorize above calculation

- Finding minimum non-diagonal is pretty simple, we can choose our pick of the following options:
    https://stackoverflow.com/questions/29394377/minimum-of-numpy-array-ignoring-diagonal

Data Structure discussion:

sequences = list of 'n' Dictionaries (name seq) #from fasta file

Dmatrix = n x n numpy array of normalized distances between the sequences.
D is the Dmatrix (distance matrix) n x n [ixj] where i and j are sequences
    Calculated from smith-waterman "score" going from one sequence i to sequence j
Distance from S(i) to S(j) = D(i,j)
If k is an internal node on our graph:
  Distance from S(k,j) (aka S(j,k)) = D(i,j) - D(i,k) 
  
 T:
    will need n leaves (a node with only one edge to rest of graph)
    will need internal_nodes a node that has more than 1 edge
    will need edges associated with all nodes 
    edges have lengths
    
    Tree:
    class TreeNode
      name = from sequences[].name 
      children: List of nodes
      edges:  list of Edges?
    
     class Edge 
      length: int # this gets recalculated every time.
    
    will need methods to create_leaves, merge_nodes, prune 

  Q is the Qmatrix n x n  numpy array.
  taking/measuring how central the cluster is; (penalizes outliers)
    where Q[i,j] = eqn1 (see above)
  
Other things to know:
  # note:  D(k,j) = D(j,k) for all values of j and k; symetry.
  # also: There will always be a distance from i to j; though they may not be directly connected.
  # also:  D(k,k) = 0 for all values of k.
  # also:  nodes will have an idea of 'degrees' as well as to which other nodes they are connected, and which may be good if something is leaf... 
  # also:  merge nodes -- we will know who the common ancestor is 
    
**Input** (begin pseudo code)

Fasta file has been imported, and returns sequences a list of Dictionarys 
                with name as key, seq as value

  # note: n = 20 as their are 20 seq in our fasta file. 
  
  Next step will be:

Initialize the Distance Matrix:

Micheal has used his Smith-Waterman code as well as creating a normalization process.
  Dmatrix is created.
  
To decrease the insanity of the next section, I have pulled out the equations....

eqn1 = $Q(i,j) = (n-2)d(i,j) - \sum_{k=1}^n d(i,k) - \sum_{k=1}^n d(j,k)$
eqn2 = Branch_length =  $d_i = \frac{d(i,j) + (r_i - r_j)/(n-2)}{2}$
   $d_j = d(i,j) - d_i$
   where $r_i = \sum_{k=1}^n d(i,k)$
eqn3 = $d(k,x) = \frac{d(i,x) + d(j,x) - d(i,j)}{2}$

**Input** \## While working with Nikaela (project 4), I implemented a read_fasta that was a hacked up ##version of Marcus's read fasta from project 3. \## \## I will put that file in our project06 directory as my_read_fasta.py \## feel free to mess with it, make it more error tolerant, etc, as I haven't checked ##. for errors on openning or ??? ##.I am adding my_read_fasta.py and the appropriate import

  for t in set from 0 to n
    T.leaf = create_leaf( S[name] ); tuple from dictionary associated with name
  
    
**Iteration**: 
q = n (number of sequences)

While q > 2:
1. Calculate Q-matrix where:
  for i in ... do
    for j in ... do
      Q[i,j] = "eqn1" (given S(i) and S(j))

2. Find pair (i,j) with minimum Q(i,j) -- 
Note that we need to be aware that the diagonal are 0s and we only need to look 
at half the array 
  find S(i) and S(j) that represent the entry in Q

3. Calculate branch lengths:    
  Now that we have i and j; 

  from "eqn2" given S(i), S(j):
  
  calculate ?2? branch lengths ( this is distance from i to j );

4. Create new node k
     "Add branches from k to i and j with branch lengths D(i) and D(j)"

5. Update distances to remaining nodes:
   using eqn3 above, re-calculate Dmatrix????

6. Remove nodes i and j from active nodes list
  we have now connected them with node k, and their edges have lengths.
 (should this happen before item 5?) ?remove column i and column j ; row i and row j

7. Add internal node k to active nodes list

8. q = q - 1

**Termination**:
When q = 2 with remaining nodes i and j:     # last node to add
- Add final branch between i and j with length d(i,j)
- Return unrooted tree T
```

# Successes
Description of the team's learning points

# Struggles
Description of the stumbling blocks the team experienced

# Personal Reflections
## Group Leader
Group leader's reflection on the project

## Other member
Other members' reflections on the project

# Generative AI Appendix
As per the syllabus
