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
