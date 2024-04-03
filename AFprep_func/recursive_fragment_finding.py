from AFprep_func.classes import Domain

def recursive_fragmentation(domains, fragment_start, first_res, last_res, min_len, max_len, overlap, max_overlap, min_overlap, cutpoints=None, overlaps=None):
    # Base case: if previous fragment end is at or beyond the last_res, we've reached the end and return the cutpoints
    if cutpoints:
        if cutpoints[-1][1] >= last_res:
            return cutpoints
    
    if cutpoints is None:
        cutpoints = []
    if overlaps is None:
        overlaps = []

    valid_cut_found = False

    for res in range(fragment_start + min_len, min(fragment_start + max_len, last_res) + 1):
        if check_valid_cutpoint(res, domains, last_res):
            next_start = find_next_start(res, overlap, max_overlap, min_overlap, domains, last_res)
            if next_start:
                valid_cut_found = True
                cutpoints.append((fragment_start, res))
                overlaps.append(res - next_start)
                
                # Recursively process the next segment
                result = recursive_fragmentation(domains, next_start, first_res, last_res, min_len, max_len, overlap, max_overlap, min_overlap, cutpoints)
                
                # If a valid fragmentation pattern is found, return the result
                if result is not None:
                    return result
                else:
                    # If the recursive call did not find a valid pattern, remove the last added cutpoints and overlaps
                    cutpoints.pop()
                    overlaps.pop()

    # If no valid cut is found in the loop, check if this is the first call and return None to indicate failure
    if not valid_cut_found and fragment_start == first_res:
        return None
    # If this is a recursive call, return None to backtrack
    elif not valid_cut_found:
        return None

def check_valid_cutpoint(res, domains, sequence_end):
    """
    Checks if a cutpoint is valid based on the specified conditions.
    
    Parameters:
    - res (int): The position to check for a valid cutpoint.
    - domains (list of Domain): The list of domain objects.
    - sequence_end (int): The end position of the protein sequence.
    
    Returns:
    - bool: True if the cutpoint is valid, False otherwise.
    """
    # Check if res is beyond the end of the sequence
    if res > sequence_end + 1:
        return False
    
    for domain in domains:
        # Check if res and res-1 are within the same domain
        if domain.start <= res <= domain.end and domain.start <= res-1 <= domain.end:
            return False
            
    return True


def find_next_start(res, overlap, max_overlap, min_overlap, domains, sequence_end):
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
