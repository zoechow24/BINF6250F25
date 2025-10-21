# BINF6250 - Project 06

## Zoe Chow, Michael Bambha, Jacque Caldwell

# Introduction

We implemented a distance-based phylogenetic tree construction method using the Neighbor-Joining (NJ) algorithm and Smith-Waterman local alignment scores. Unlike UPGMA (Unweighted Pair Group Method using Arithmetic averages), NJ produces unrooted trees and does not assume a constant evolutionary rate across lineages, making it more biologically realistic for analyzing sequence relationships.

We will use HIV-1 reverse transcriptase sequences to construct our phylogenetic tree. The Smith-Waterman algorithm will be used to generate pairwise local alignment scores, which will then be converted into distances for tree construction. This approach is particularly suitable for HIV sequence analysis as it: 1. Handles sequence variations effectively through local alignment 2. Accounts for potential rate heterogeneity across different viral strains 3. Does not assume a molecular clock

The ultimate output will be an unrooted phylogenetic tree representing the evolutionary relationships between the HIV-1 sequences, visualized using the ete3 library. This method provides insights into viral diversity and evolutionary patterns while avoiding the assumptions of simpler hierarchical clustering approaches.

The key innovations of this implementation are: - Use of Smith-Waterman for sensitive local alignment scoring - Implementation of Neighbor-Joining for unrooted tree construction - More biologically realistic evolutionary model

# Pseudocode

```         
For building the distance matrix:

- We already did SW implementation. However we probably don't need the traceback
    - Just need to return the max score for each pair of alignments
    - From there, we need to take our "similarity" scores, normalize them w/ cosine similarity.
```

In geometry, cosine similarity between two vectors **A** and **B** is defined as:

$ \text{cosine similarity}(A, B) = \frac{A \cdot B}{\sqrt{A^\top A} \, \sqrt{B^\top B}} = \frac{\sum_{i=1}^{n} A_i B_i}{\sqrt{\sum_{i=1}^{n} A_i^2} \, \sqrt{\sum_{i=1}^{n} B_i^2}} $

We can apply this concept to our Smith-Waterman alignment as well!

\[ S\_{\mathrm{norm}}(A, B) = \frac{S(A, B)}{\sqrt{S(A, A)\,S(B, B)}} \]

Since SW really is a measure of 'similarity' not distance, we can convert these normalized scores to distances:

\[ D(A,B) ;=; 1 - S\_{\mathrm{norm}}(A,B) ;=; 1 - \frac{S(A,B)}{\sqrt{S(A,A)\,S(B,B)}} \]​

```         
- Finding minimum non-diagonal is pretty simple, since we only use Q-matrix for finding min, we can fill diag with inf
    - Probably not an issue modifying it in-place since we need to recalculate it again from DM anyway.

- Q-matrix:
    1. Find number of current rows or columns (# of species/taxa) = n
    2. The summations in the Q-formula are just the row or column sums ("net divergence")
    3. Given that i,j is just a point in our matrix, this is actually pretty straightforward.
```

-   Branch length:

$d_i = \frac{d(i,j) + (r_i - r_j)/(n-2)}{2}$ $d_j = d(i,j) - d_i$ where $r_i = \sum_{k=1}^n d(i,k)$

```         
This is straightforward, just some linear algebra and numpy involved. Find the row sums, find the column sums,
find d_ij and apply the formula.

- Updating distances
    - Add a new row/col k (perhaps to the end to avoid indexing issues)
    - Compute distances using the formula provided to our new node K - get back an np array
        - Probably doesn't make sense to compute i,k and j,k. So k will be smaller and can only be added
        once we delete i,j
    - Delete row/cols associated w/ i, j (np.delete)
    - Add in our computed k-distance array


- Neighbor joining:
    - Initialize all of the labels to be leaf nodes (just create a list of `Node` instantiations)
    - Find n_taxa, which is either the # of rows or cols of our DM
    - Build q-matrix
    - Find min distance
    - Find branch lengths
    - Create k 
    - Remove nodes i, j from our originally instantiated list and add k instead.
    - Update DM and decrement n_taxa by 1
    - Continue everything until we have only 2x2 matrix, add the final distances.
    - Convert our final class representation to a newick string and plot!
```

# Successes

Huge dive into `numpy` methods, and we were able to find efficient ways to manipulate our matrices using linear algebra / `numpy` instead of using `for` loops. Successfully able to create our tree and implement NJ.

Small dive into ete3 and why it doesn't work with current versions of python, and many other options available (and some not so available) for visualizing Newick trees.

# Struggles

Description of the stumbling blocks the team experienced

# Personal Reflections

## Group Leader

Group leader's reflection on the project

## Other member

## Michael

This was pretty tough. Honestly, a large portion of this was simply just linear algebra and I actually think I learned quite a bit about some matrix operations and ways to implement them in Python. One of them being the self-score calculation trick using `np.outer(np.diag, np.diag)` which I thought was really cool. The hardest part was definitely the OOP and the branch lengths for me - while figuring out how to efficiently implement the formulas in some of the steps was somewhat time-consuming, it wasn't necessarily conceptually difficult or confusing. That part was kind of plug-and-chug, and then optimize. But figuring out how to properly implement our `Node` class was tough, and then figuring out how to actually use this class once we made it to map our branch lengths to our new nodes and update the node list as we continued through the algorithm was conceptually more complicated. The recursion for the Newick strings was also pretty time consuming as well.

Overall, it felt like most of the individual steps of the algorithm were actually pretty easy, but actually putting them together was much harder.

## Jacque

I was away for a large chunk of time during this project, I did spend a bunch of time early on to understand the algorithm, and feel that I have a pretty good handle on what it is doing.  A lot of unexpected things came up during this that were difficult to forsee during the pseudocoding part of this project.  In particular planning the data structure that was needed for the tree, and planning how that was going to be used was much harder up front, but made itself clearer as the project continued.  I did struggle a bit to try to get ete3 going, but the fact that the CGI libraries that it depends on are no longer part of the base python code made it much harder to get going.  While we did have some early success working with biopython and matplotlib, they stopped working for us towards the time for wrapping up this project, so I used the Newick string that we generated to create a tree on the on-line Ete3 web site.  Asynchronous programming continues to be difficult, especially when using a single file for programming between three people in three different time zones.  Rstudio also continues to make moving between github and Explorer difficult at times, with many errors popping up just because spaces changed to tabs or vice versa when moving files between filesystems.  

# Generative AI Appendix

ChatGPT was used for the LaTeX formatting on the cosine similarity / Smith-Waterman normalization formulas in our Pseudocode section above. (This is just so that our reviewers can see this in one place and don't necessarily need to look it up).