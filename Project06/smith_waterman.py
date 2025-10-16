"""
File: smith_waterman.py
Description: Small Python module for finding
the maximum score of a Smith-Waterman alignment
"""
from enum import Enum
import numpy as np

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


def smith_waterman_max_score(
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
    # traceback_matrix = np.full((rows, cols), Direction.END, dtype=object)

    max_score = 0
    # max_pos = (0, 0)

    for i in range(1, rows):
        for j in range(1, cols):
            direction, score = cal_score(
                score_matrix, seq1, seq2, i, j, match_score, mismatch_score, gap_score
            )
            score_matrix[i, j] = score
            # traceback_matrix[i, j] = direction

            max_score = max(max_score, score)
                # max_pos = (i, j)

    return max_score
