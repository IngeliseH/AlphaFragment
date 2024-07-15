""""
Internal utility functions for protein fragmention process.

Functions:
    - validate_fragmentation_parameters: Validates the parameters used for
      protein fragmentation.
    - merge_overlapping_domains: Merges overlapping domains within a list of
      domains.
    - check_valid_cutpoint: Helper function to validate potential fragment
      boundaries.
    - recursive_fragmentation: Main function for recursively generating fragments.
    
Dependencies:
    - Domain: A class representing a domain within a protein sequence.
    - Protein: A class representing a protein sequence.
"""
from .classes import Domain, Protein

def validate_fragmentation_parameters(protein, length, overlap):
    """
    Validates the parameters used for protein fragmentation.
    
    Parameters:
        - protein (Protein): The protein object to be fragmented.
        - length (dict): Dictionary containing the ideal, minimum, and maximum
          length values, in the format:
          {'min': min_len, 'ideal': ideal_len, 'max': max_len}
          where min_len, ideal_len, and max_len are all integers,
          with min_len <= ideal_len <= max_len.
        - overlap (dict): Dictionary containing the ideal, minimum, and maximum
          overlap values, in the format:
          {'min': min_overlap, 'ideal': ideal_overlap, 'max': max_overlap}
          where min_overlap, ideal_overlap, and max_overlap are all integers, with
          min_overlap <= ideal_overlap <= max_overlap.

    Returns:
        - None

    Raises:
        - ValueError: If any of the parameter validations fail.
        - TypeError: If the protein input is not an instance of the Protein class.
    """
    # Check that the protein input is an instance of the Protein class
    if not isinstance(protein, Protein):
        raise TypeError("Input protein must be an instance of the Protein class.")

    # Check that length is a dictionary
    if not isinstance(length, dict):
        raise TypeError("Length must be a dictionary.")

    # Check that the length dictionary contains the required keys
    if not all(key in length for key in ['min', 'ideal', 'max']):
        raise ValueError("Length dictionary must contain keys 'min', 'ideal', and 'max'.")

    # Check that the length values are integers
    if not all(isinstance(length[key], int) for key in ['min', 'ideal', 'max']):
        raise TypeError("Length values must be integers.")

    # Check that the min length is less than the max length
    if length['min'] > length['max']:
        raise ValueError(f"Minimum fragment length ({length['min']}) must be less than "
                         f"maximum fragment length ({length['max']}).")

    # Check that the ideal length is within the min and max length bounds
    if not length['min'] <= length['ideal'] <= length['max']:
        raise ValueError("Ideal length must be within the min and max length bounds.")

    # Check that the min_len is greater than 0
    if length['min'] <= 0:
        raise ValueError("Minimum fragment length must be greater than 0.")

    # Check that overlap is a dictionary
    if not isinstance(overlap, dict):
        raise TypeError("Overlap must be a dictionary.")

    # Check that the overlap dictionary contains the required keys
    if not all(key in overlap for key in ['min', 'ideal', 'max']):
        raise ValueError("Overlap dictionary must contain keys 'min', 'ideal', and 'max'.")

    # Check that the overlap values are integers
    if not all(isinstance(overlap[key], int) for key in ['min', 'ideal', 'max']):
        raise TypeError("Overlap values must be integers.")

    # Check that the minimum overlap is less than or equal to the maximum overlap
    if overlap['min'] > overlap['max']:
        raise ValueError(f"Minimum overlap ({overlap['min']}) must be less than "
                         f"or equal to maximum overlap ({overlap['max']}).")

    # Check that the ideal overlap is within the min and max overlap bounds
    if not overlap['min'] <= overlap['ideal'] <= overlap['max']:
        raise ValueError("Ideal overlap must be within the min and max overlap bounds.")

    # Check that the maximum overlap is less than the minimum fragment length
    if overlap['max'] >= length['min']:
        raise ValueError(f"Maximum overlap ({overlap['max']}) must be less than the minimum "
                         f"fragment length ({length['min']}) to avoid overlap-length conflicts.")
def merge_overlapping_domains(domains):
    """
    Merges overlapping domains within a list of domains.

    Parameters:
        - domains (list of Domain): List of domain objects.

    Returns:
        - list of Domain: A list of domains where overlapping domains have been
          merged into single entries.
    """
    # Sort domains by their start positions
    sorted_domains = sorted(domains, key=lambda x: x.start)
    combined_domains = []

    for domain in sorted_domains:
        if not combined_domains:
            combined_domains.append(domain)
        else:
            last = combined_domains[-1]
            # Check if the current domain overlaps with the last one in the list
            if domain.start <= last.end:
                # Merge the two domains by updating the end of the last domain
                combined_domains[-1] = Domain(last.id, last.start,
                                              max(last.end, domain.end), last.type)
            else:
                combined_domains.append(domain)

    return combined_domains

def check_valid_cutpoint(res, domains, sequence_end):
    """
    Checks if a slicing index is a valid cutpoint.

    Parameters:
        - res (int): The residue position to check (will be sliced before this residue).
        - domains (list of Domain): The domains within the protein.
        - sequence_end (int): The last residue position in the protein sequence.

    Returns:
        - bool: True if the residue position is a valid cutpoint; False otherwise.
    """
    # Check if res is beyond the end or before the start of the sequence
    if res > sequence_end + 1:
        return False
    if res < 0:
        return False

    # If slicing index will cut at end of sequence, this is always valid
    if res == (sequence_end+1):
        return True

    for domain in domains:
        # Check if current and previous res are within the same domain
        # (if both are in same domain, cutting before res would split the domain)
        if domain.start <= res <= domain.end and domain.start <= res-1 <= domain.end:
            return False

    return True

def recursive_fragmentation(protein, domains, fragment_start, length, overlap, cutpoints=None):
    """
    Recursively splits a protein sequence into overlapping fragments, avoiding
    breaking domains.

    Parameters:
        - protein (Protein): The protein object to fragment.
        - domains (list of Domain): The list of domains within the protein -
          doesn't use protein.domain_list as overlapping domains should be merged.
        - fragment_start (int): The starting position for fragmentation.
        - length (dict): Dictionary containing the ideal, minimum, and maximum
          length values, in the format:
          {'min': min_len, 'ideal': ideal_len, 'max': max_len}
          where min_len, ideal_len, and max_len are all integers,
          with min_len <= ideal_len <= max_len.
        - overlap (dict): Dictionary containing the ideal, minimum, and maximum
          overlap values, in the format:
          {'min': min_overlap, 'ideal': ideal_overlap, 'max': max_overlap}
          where min_overlap, ideal_overlap, and max_overlap are all integers,
          with min_overlap <= ideal_overlap <= max_overlap.
        - cutpoints (list of tuples, optional): Accumulator for storing fragment
          cutpoints.

    Returns:
        - list of tuples or None: The list of fragment cutpoints if successful;
          otherwise, None.
    """
    validate_fragmentation_parameters(protein, length, overlap)

    def find_next_start(res):
        # Use ideal overlap if possible
        if check_valid_cutpoint(res - overlap['ideal'], domains, protein.last_res):
            return res - overlap['ideal']
        # Force None if moving current fragment end would allow better overlap with new fragment
        for forwards_res in range(overlap['max'], overlap['min'] - 1, -1):
            if check_valid_cutpoint(res + forwards_res, domains, protein.last_res):
                return None
        # Attempt to find a valid cutpoint by first increasing, then decreasing overlap
        for adjusted_overlap in (list(range(overlap['ideal'] + 1, overlap['max'] + 1)) +
                                 list(range(overlap['ideal'] - 1, overlap['min'] - 1, -1))):
            if check_valid_cutpoint(res - adjusted_overlap, domains, protein.last_res):
                return res - adjusted_overlap
        # If no valid cutpoint is found within overlap boundaries, return None
        return None

    if cutpoints is None:
        cutpoints = []

    # Iterate over possible fragment end cutpoints, iterating from ideal length to max, then ideal to min
    for l in (list(range(length['ideal'], length['max'] + 1)) +
              list(range(length['ideal'] - 1, length['min'] - 1, -1))):

        res = fragment_start + l

        if check_valid_cutpoint(res, domains, protein.last_res):
            # If the current fragment end is at the end of the protein, finalize here.
            if res == protein.last_res + 1:
                cutpoints.append((fragment_start, res))
                return cutpoints
            # If a valid cutpoint to start the next fragment is found, add it to the list and continue
            next_start = find_next_start(res)
            if next_start:
                cutpoints.append((fragment_start, res))

                # Recursively process the next segment
                result = recursive_fragmentation(protein, domains, next_start,
                                                 length, overlap, cutpoints)

                # If a valid fragmentation pattern is found, return the result
                if result is not None:
                    return result
                # If the recursive call did not find a valid pattern, remove the
                # last added cutpoints
                cutpoints.pop()

    # If no valid cut is found in the loop, return None to indicate failure
    return None
