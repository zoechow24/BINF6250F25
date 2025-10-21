"""
File: smith_waterman.py
Description: Small Python module for finding
the maximum score of a Smith-Waterman alignment
"""
from enum import Enum
from itertools import filterfalse
from functools import lru_cache
import math
from typing import Tuple, Dict, Callable
import warnings
import numpy as np

try:
    import pandas as pd
    _HAS_PANDAS = True
except ImportError:
    _HAS_PANDAS = False
    

class Direction(Enum):
    END = 0
    DIAG = 1
    UP = 2
    LEFT = 3

def cal_score(
    matrix: np.ndarray,
    seq1: str,
    seq2: str,
    i: int,
    j: int,
    match_score: int,
    mismatch_score: int,
    gap_score: int,
):
    """Calculate score for position (i,j) in scoring matrix, also record move to trace back

    Args:
        matrix (numpy array): scoring matrix
        seq1 (str): sequence 1
        seq2 (str): sequence 2
        i (int): current row number
        j (int): current column number

    Returns:
        score in position (i,j)
        move to trace back: 0-END, 1-DIAG, 2-UP, 3-LEFT

    Pseudocode:
        Calculate scores based on upper-left, up, and left neighbors:
            diag_score = upper-left + (match or mismatch)
            up_score = up + gap
            left_score = left + gap
        score = max(0, diag_score, up_score, left_score)
        traceback = maximum direction or end

    """
    up_neighbor = (i - 1, j)
    diag_neighbor = (i - 1, j - 1)
    left_neighbor = (i, j - 1)

    up_score = matrix[up_neighbor] + gap_score
    left_score = matrix[left_neighbor] + gap_score
    diag_score = matrix[diag_neighbor] + (
        match_score if seq1[i - 1] == seq2[j - 1] else mismatch_score
    )

    score_dict = {
        Direction.DIAG: diag_score,
        Direction.UP: up_score,
        Direction.LEFT: left_score,
        Direction.END: 0,
    }

    max_direction, max_score = max(score_dict.items(), key=lambda score: score[1])
    return max_direction, max_score

def traceback(
    seq1: str, seq2: str, traceback_matrix: np.ndarray, maximum_position: int
):
    """Find the optimal path through scoring marix

        diagonal: match/mismatch
        up: gap in seq1
        left: gap in seq2

    Args:
        seq1 (str) : First sequence being aligned
        seq2 (str) : Second sequence being aligned
        traceback_matrix (numpy array): traceback matrix
        maximum_position (tuple): starting position to trace back from

    Returns:
        aligned_seq1 (str): e.g. GTTGAC
        aligned_seq2 (str): e.g. GTT-AC

    Pseudocode:

        while current_move != END:
            current_move = traceback_matrix[current_row][current_col]
            if current_move == DIAG:
               ...
            elif current_move == UP:
                ...
            elif current_move == LEFT:
                ...

               base1 =  seq1[current_row] | base2 = seq2[current_col]

        aligned_seq1 = aligned_seq1[::-1]
    """
    aligned_seq1 = ""
    aligned_seq2 = ""
    current_move = traceback_matrix[maximum_position]
    current_row, current_col = maximum_position

    while current_move != Direction.END:

        if current_move == Direction.DIAG:
            aligned_seq1 += seq1[current_row - 1]
            aligned_seq2 += seq2[current_col - 1]
            current_row -= 1
            current_col -= 1

        elif current_move == Direction.LEFT:
            aligned_seq1 += seq1[current_row - 1]
            aligned_seq2 += "-"
            current_col -= 1

        elif current_move == Direction.UP:
            aligned_seq1 += "-"
            aligned_seq2 += seq2[current_col - 1]
            current_row -= 1

        current_move = traceback_matrix[current_row, current_col]

    return aligned_seq1[::-1], aligned_seq2[::-1]


def smith_waterman(
    seq1: str,
    seq2: str,
    match_score: int = 1,
    mismatch_score: int = -1,
    gap_score: int = -1,
):
    """Smith-Waterman algorithm for local alignment

    Args:
        seq1 (str): input seq 1
        seq2 (str): input seq 2
        match: default = +1
        mismatch: default = -1
        gap: default = -1

    Returns:
        aligned_seq1 (str)
        aligned_seq2 (str)
        score_matrix (numpy array): scoring matrix
    """
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    score_matrix = np.zeros((rows, cols))
    # using an object array isn't necessarily efficient, but is readable
    # since we aren't really using np vectorization methods, shouldn't matter
    # if performance was a concern, we could use IntEnum instead
    traceback_matrix = np.full((rows, cols), Direction.END, dtype=object)

    max_score = 0
    max_pos = (0, 0)

    for i in range(1, rows):
        for j in range(1, cols):
            direction, score = cal_score(
                score_matrix, seq1, seq2, i, j, match_score, mismatch_score, gap_score
            )
            score_matrix[i, j] = score
            traceback_matrix[i, j] = direction

            if score > max_score:
                max_score = score
                max_pos = (i, j)

    aligned_seq1, aligned_seq2 = traceback(
        seq1, seq2, traceback_matrix, maximum_position=max_pos
    )

    return aligned_seq1, aligned_seq2, score_matrix


def smith_waterman_max_score(
    seq1: str,
    seq2: str,
    match_score: int = 1,
    mismatch_score: int = -1,
    gap_score: int = -1,
) -> float | int:
    """Find the maximum score of a local alignment between
    two sequences using Smith-Waterman. Does not do
    traceback or find the actual alignment -- just returns
    the maximum score.

    Args:
        seq1 (str): input seq 1
        seq2 (str): input seq 2
        match (float | int): default = +1
        mismatch (float | int): default = -1
        gap (float | int): default = -1

    Returns:
        int: Maximum score for SW alignment.
    """
    rows = len(seq1) + 1
    cols = len(seq2) + 1

    score_matrix = np.zeros((rows, cols))
    max_score = 0

    for i in range(1, rows):
        for j in range(1, cols):
            _, score = cal_score(
                score_matrix, seq1, seq2, i, j, match_score, mismatch_score, gap_score
            )
            score_matrix[i, j] = score
            max_score = max(max_score, score)

    return max_score


def calculate_p_distance(seq1: str, seq2: str) -> float:
    """Calculate the p-distance between two aligned sequences.

    Args:
        seq1 (str): First aligned sequence 1
        seq2 (str): Second aligned sequence

    Returns:
        float: p distance between 0 and 1
    """
    valid_pairs = ((a, b) for a, b in zip(seq1, seq2) if "-" not in (a, b))
    mismatches = total = 0
    for a, b in valid_pairs:
        total += 1
        mismatches += a != b
    return round(mismatches / total, 6) if total else 0.0


def jukes_cantor(aligned_seq1: str, aligned_seq2: str) -> float:
    """Corrects p-distances of two sequences based on Jukes-Cantor evolutionary model.

    Args:
        aligned_seq1 (str): Aligned sequence 1
        aligned_seq2 (str): Aligned sequence 2

    Returns:
        float: Corrected p-distances based on JC
    """
    p_dist = calculate_p_distance(aligned_seq1, aligned_seq2)
    if p_dist >= 0.75:
        return float('inf')
    dist = (-3/4) * math.log(1 - (4 * p_dist / 3))
    return round(dist, 6)


def calculate_transition_proportions(seq1: str, seq2: str) -> Tuple[float, float]:
    """Find the transition and transversion proportions of two sequences.

    Args:
        seq1 (str): Aligned sequence 1
        seq2 (str): Aligned sequence 2

    Returns:
        Tuple[float, float]: Tuple of (transition, transversion) proportions.
    """
    purines = {'A', 'G'}
    pyrimidines = {'C', 'T'}

    pairs = filterfalse(lambda ab: '-' in ab, zip(seq1.upper(), seq2.upper()))

    transitions = transversions = valid_sites = 0

    for a, b in pairs:
        valid_sites += 1
        if a == b:
            continue
        if (a in purines and b in purines) or (a in pyrimidines and b in pyrimidines):
            transitions += 1
        else:
            transversions += 1

    if valid_sites == 0:
        return 0.0, 0.0

    return transitions / valid_sites, transversions / valid_sites


def kimura_two_parameter(aligned_seq1: str, aligned_seq2: str) -> float:
    """Find the K2P distance based on proportion of
    transitions and transversions from two sequences.

    Args:
        aligned_seq1 (str): Aligned sequence 1
        aligned_seq2 (str): Aligned sequence 2

    Returns:
        float: K2P corrected p-distances
    """
    P, Q = calculate_transition_proportions(aligned_seq1, aligned_seq2)
    if (1 - 2*P - Q) <= 0 or (1 - 2*Q) <= 0:
        return float('inf')
    dist = -0.5 * math.log(1 - 2*P - Q) - 0.25 * math.log(1 - 2*Q)
    return round(dist, 6)


@lru_cache(maxsize=1)
def load_blosum(blosum_f: str = "BLOSUM62") -> pd.DataFrame:
    return pd.read_csv(blosum_f, sep=r"\s+", index_col=0, header=0)


def blosum_dist(aligned_seq1: str = None, aligned_seq2: str = None) -> float:
    """Calculate normalized distances from a BLOSUM62 matrix between two aligned sequences.

    Args:
        aligned_seq1 (str, optional): First aligned sequence. Defaults to None.
        aligned_seq2 (str, optional): Second aligned sequence. Defaults to None.

    Returns:
        float: Score based on BLOSUM62 matrix, normalized by the maximum score of aligned_seq1 to itself
    """
    if not _HAS_PANDAS:
        warnings.warn("Pandas not found -- BLOSUM distances will be skipped.")
        return float("nan")
    blosum_df = load_blosum() # cached
    # to vectorize, we can use .to_numpy() 
    blosum_matrix = blosum_df.to_numpy()
    residues = blosum_df.index.tolist()
    res_to_idx = {res: i for i, res in enumerate(residues)}

    # ignore pairs which have a gap on either seq
    pairs = [(res_to_idx[r1], res_to_idx[r2]) 
            for r1, r2 in zip(aligned_seq1, aligned_seq2) 
            if '-' not in (r1, r2)]

    if not pairs:
        return 1.0

    # transpose -- pairs will give a column vector and we need rows
    idx1, idx2 = np.array(pairs).T
    # below is where the actual speed gains from vectorized lookups come up
    scores = blosum_matrix[idx1, idx2]
    max_scores = blosum_matrix[idx1, idx1]

    return round(1 - (scores.sum() / max_scores.sum()), 6)

DISTANCE_FUNCTIONS: Dict[str, Callable[[str, str], float]] = {
    "p": calculate_p_distance,
    "jc": jukes_cantor,
    "jukes_cantor": jukes_cantor,
    "k2p": kimura_two_parameter,
    "kimura": kimura_two_parameter,
    "blosum": blosum_dist,
}