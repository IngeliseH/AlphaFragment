"""
Contains functions to split a protein sequence into fragments.
"""

def recursive_fragmentation(protein, fragment_start, min_len, max_len,
                            overlap, min_overlap, max_overlap, cutpoints=None):
    """
    Recursively splits a protein sequence into overlapping fragments, avoiding breaking domains.

    Parameters:
    - protein (Protein): The protein object to fragment.
    - fragment_start (int): The starting position for fragmentation.
    - min_len (int): Minimum allowed fragment length.
    - max_len (int): Maximum allowed fragment length. (May be increased in the
      fragmentation process)
    - overlap (int): Ideal overlap between fragments.
    - min_overlap (int): Minimum allowable overlap between fragments.
    - max_overlap (int): Maximum allowable overlap between fragments.
    - cutpoints (list of tuples, optional): Accumulator for storing fragment cutpoints.

    Returns:
    - list of tuples or None: The list of fragment cutpoints if successful; otherwise, None.
    """
    def find_next_start(res):
        # Use ideal overlap if possible
        if check_valid_cutpoint(res - overlap, protein.domain_list, protein.last_res):
            return res - overlap
        # Force None if moving fragment end would allow better overlap
        for forwards_res in range(max_overlap, min_overlap - 1, -1):
            if check_valid_cutpoint(res + forwards_res, protein.domain_list, protein.last_res):
                return None
        # Attempt to find a valid cutpoint by first increasing, then decreasing overlap
        for adjusted_overlap in (list(range(overlap + 1, max_overlap + 1)) +
                                 list(range(overlap - 1, min_overlap - 1, -1))):
            if check_valid_cutpoint(res - adjusted_overlap, protein.domain_list, protein.last_res):
                return res - adjusted_overlap
        # If no valid cutpoint is found within overlap boundaries, return None
        return None
    
    # Base case: if previous fragment end is the end of the protein, we've
    # reached the end and return the cutpoints
    if cutpoints:
        if cutpoints[-1][1] >= protein.last_res:
            return cutpoints

    if cutpoints is None:
        cutpoints = []

    for res in range(fragment_start + min_len, min(fragment_start + max_len, protein.last_res) + 1):
        if check_valid_cutpoint(res, protein.domain_list, protein.last_res):
            next_start = find_next_start(res)
            if next_start:
                cutpoints.append((fragment_start, res))

                # Recursively process the next segment
                result = recursive_fragmentation(protein, next_start, min_len, max_len,
                                                 overlap, max_overlap, min_overlap,
                                                 cutpoints)

                # If a valid fragmentation pattern is found, return the result
                if result is not None:
                    return result
                # If the recursive call did not find a valid pattern, remove the
                # last added cutpoints
                cutpoints.pop()

    # If no valid cut is found in the loop, return None to indicate failure
    return None

def check_valid_cutpoint(res, domains, sequence_end):
    """
    Checks if a residue position is a valid cutpoint.

    Parameters:
    - res (int): The residue position to check.
    - domains (list of Domain): The domains within the protein.
    - sequence_end (int): The last residue position in the protein sequence.

    Returns:
    - bool: True if the residue position is a valid cutpoint; False otherwise.
    """
    # Check if res is beyond the end or before the start of the sequence
    if res > sequence_end:
        return False
    if res < 0:
        return False

    # If res is at end of sequence, doesn't matter if it is in a domain
    if res == (sequence_end):
        return True

    for domain in domains:
        # Check if current and previous res are within the same domain
        # if both are in same domain, cutting at res will split the domain -
        # if they are in different domains, this is fine
        if domain.start <= res <= domain.end and domain.start <= res-1 <= domain.end:
            return False

    return True
