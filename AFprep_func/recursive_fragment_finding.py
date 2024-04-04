"""
Contains functions to split a protein sequence into fragments.
Leverages the Protein class for accessing protein attributes like
domain list and sequence endpoints.
"""

from AFprep_func.classes import Protein

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
    # Base case: if previous fragment end is the end of the protein, we've
    # reached the end and return the cutpoints
    if cutpoints:
        if cutpoints[-1][1] >= protein.last_res:
            return cutpoints

    if cutpoints is None:
        cutpoints = []

    for res in range(fragment_start + min_len, min(fragment_start + max_len, protein.last_res) + 1):
        if check_valid_cutpoint(res, protein.domain_list, protein.last_res):
            next_start = find_next_start(res, overlap, min_overlap, max_overlap,
                                         protein.domain_list, protein.last_res)
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

def find_next_start(res, overlap, min_overlap, max_overlap, domains, sequence_end):
    """
    Determines the next start position for a fragment based on overlap criteria.
    (Start position will be before the end of the previous fragment)

    Parameters:
    - res (int): The current residue position.
    - overlap (int): Ideal overlap between fragments.
    - min_overlap (int): Minimum allowable overlap between fragments.
    - max_overlap (int): Maximum allowable overlap between fragments.
    - domains (list of Domain): The domains within the protein.
    - sequence_end (int): The last residue position in the protein sequence.

    Returns:
    - int or None: The next start position if an allowed start position is
      found; otherwise, None.

    Note:
    - None is returned if moving previous fragment end would allow better
      overlap - this is to prevent defaulting to very short overlaps when a
      better solution is available
    """
    if check_valid_cutpoint(res - overlap, domains, sequence_end):
        return res - overlap
    #force none if moving fragment end would allow better overlap
    for forwards_res in range(max_overlap, min_overlap, -1):
        if check_valid_cutpoint(res + forwards_res, domains, sequence_end):
            return None
    for increased_overlap in range(overlap + 1, max_overlap):
        if check_valid_cutpoint(res - increased_overlap, domains, sequence_end):
            return res - increased_overlap
    for decreased_overlap in range(overlap - 1, min_overlap, -1):
        if check_valid_cutpoint(res - decreased_overlap, domains, sequence_end):
            return res - decreased_overlap
    return None

def fragment_protein(protein, min_len, max_len, overlap, min_overlap, max_overlap):
    """
    Fragments a protein using recursive fragmentation.

    Parameters:
    - protein (Protein): The protein object to fragment.
    - min_len (int): Minimum fragment length.
    - max_len (int): Maximum fragment length. (May be increased in the fragmentation process)
    - overlap (int): Ideal overlap between fragments.
    - min_overlap (int): Minimum allowable overlap between fragments.
    - max_overlap (int): Maximum allowable overlap between fragments.

    Returns:
    - list of tuples: The list of fragment start and end positions.

    Note:
    - Does not adjust well for long domains - ideally used after long domains
      have been identifed and handled.
    """
    if not isinstance(protein, Protein):
        raise ValueError("protein must be an instance of Protein.")
    fragments = None

    #recursive fragmentation
    while fragments is None:
        #deal with short proteins by classifying whole protein as one fragment
        #included in while loop so as max length increases this is still happening
        if len(protein.sequence) <= max_len:
            fragments = [(protein.first_res, protein.last_res)]
            continue

        fragments = recursive_fragmentation(protein, protein.first_res, min_len,
                                            max_len, overlap, min_overlap, max_overlap)
        if fragments is None:
            #min_len = max (min_len-10, max_overlap + 1)
            #print(min_len)
            max_len = min(max_len+10, len(protein.sequence))
            print(max_len)
    return fragments
