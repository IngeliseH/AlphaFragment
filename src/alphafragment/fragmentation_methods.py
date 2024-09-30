""""
Internal utility functions for protein fragmention process.

Functions:
    - validate_fragmentation_parameters: Validates the parameters used for
      protein fragmentation.
    - merge_overlapping_domains: Merges overlapping domains within a list of
      domains.
    - check_valid_cutpoint: Helper function to validate potential fragment
      boundaries.
    - find_next_start: Finds the next valid fragment start position.
    - recursive_fragmentation: Main function for recursively generating fragments.
    - break_in_half: Splits a given Protein or ProteinSubsection object into two, ensuring
      no domains are broken and that the new subsections overlap.
    
Dependencies:
    - time: Used for setting a time limit for recursive fragmentation.
    - Domain: A class representing a domain within a protein sequence.
    - Protein: A class representing a protein sequence.
    - ProteinSubsection: A class representing a subsection of a protein sequence.
"""
import time
from .classes import Domain, Protein, ProteinSubsection

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
    
    # Check that the min overlap is not negative
    if overlap['min'] < 0:
        raise ValueError("Minimum overlap cannot be negative.")

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

def find_next_start(res, protein, domains, overlap):
        """
        Finds the next valid fragment start position.

        Parameters:
            - res (int): The current residue position.
            - protein (Protein): The protein object to be fragmented.
            - domains (list of Domain): The domains within the protein.
            - overlap (dict): Dictionary containing the ideal, minimum, and maximum
              overlap values, in the format:
              {'min': min_overlap, 'ideal': ideal_overlap, 'max': max_overlap}
              where min_overlap, ideal_overlap, and max_overlap are all integers,
              with min_overlap <= ideal_overlap <= max_overlap.
        
        Returns:
            - int or None: The next valid fragment start position if found; otherwise, None.
        """
        # Attempt to find a valid cutpoint by increasing overlap from ideal to max
        for overlap_adjusted in range(overlap['ideal'], overlap['max'] + 1):
            if check_valid_cutpoint(res - overlap_adjusted, domains, protein.last_res):
                return res - overlap_adjusted
        # Check for largest valid cutpoint
        possible_current_overlap = None
        for overlap_adjusted in range(overlap['ideal'] - 1, overlap['min'] - 1, -1):
            if check_valid_cutpoint(res - overlap_adjusted, domains, protein.last_res):
                possible_current_overlap = res - overlap_adjusted
                break
        if possible_current_overlap:
            # Force None if moving current fragment end would allow better overlap with new fragment
            #for forwards_res in range((overlap['ideal'] - possible_current_overlap), (overlap['max'] - possible_current_overlap + 1)):
            #    if check_valid_cutpoint(res + forwards_res, domains, protein.last_res):
            #        return None
            return res - overlap_adjusted
        # If no valid cutpoint is found within overlap boundaries, return None
        return None

def recursive_fragmentation(protein, domains, fragment_start, length, overlap, cutpoints=None, time_limit=None, start_time=None):
    """
    Recursively splits a protein sequence into overlapping fragments, avoiding
    breaking domains. If the process exceeds a specified time limit, it returns
    None, signaling that the protein should be split further.

    Parameters:
        - protein (Protein): The protein object to fragment.
        - domains (list of Domain): The list of domains within the protein -
          doesn't use protein.domain_list as overlapping domains should be merged.
        - fragment_start (int): The starting position for fragmentation.
        - length (dict): Dictionary containing the ideal, minimum, and maximum
          length values.
        - overlap (dict): Dictionary containing the ideal, minimum, and maximum
          overlap values.
        - cutpoints (list of tuples, optional): Accumulator for storing fragment
          cutpoints.
        - time_limit (float, optional): The maximum time allowed for the operation, in seconds.
        - start_time (float, optional): The start time of the operation.

    Returns:
        - list of tuples or "TIME_LIMIT_EXCEEDED": The list of fragment cutpoints if successful;
          otherwise, "TIME_LIMIT_EXCEEDED" if the time limit is exceeded, or None if no valid fragmentation pattern is found.
    """
    validate_fragmentation_parameters(protein, length, overlap)

    if cutpoints is None:
        cutpoints = []

    # Initialize the start time if not provided
    if start_time is None:
        start_time = time.time()

    # Check if the time limit has been reached
    if time_limit and (time.time() - start_time) > time_limit:
        return "TIME_LIMIT_EXCEEDED"  # Signal that the time limit has been exceeded

    # Iterate over possible fragment end cutpoints
    for l in (list(range(length['ideal'], length['max'] + 1)) +
              list(range(length['ideal'] - 1, length['min'] - 1, -1))):

        res = fragment_start + l

        if check_valid_cutpoint(res, domains, protein.last_res):
            # If the current fragment end is at the end of the protein, finalize here.
            if res == protein.last_res + 1:
                cutpoints.append((fragment_start, res))
                return cutpoints

            next_start = find_next_start(res, protein, domains, overlap)
            if next_start:
                cutpoints.append((fragment_start, res))

                # Recursively process the next segment
                result = recursive_fragmentation(protein, domains, next_start,
                                                 length, overlap, cutpoints, time_limit, start_time)

                # If a valid fragmentation pattern is found, return the result
                if result == "TIME_LIMIT_EXCEEDED":
                    return "TIME_LIMIT_EXCEEDED"  # Propagate the time limit signal

                if result is not None:
                    return result
                # If the recursive call did not find a valid pattern, remove the
                # last added cutpoints
                cutpoints.pop()

    return None  # Return None to indicate that no valid fragmentation was found

def break_in_half(protein, length, overlap):
    """
    Splits a given Protein or ProteinSubsection object into two subsections, ensuring no domains are broken
    and that the subsections overlap. The split is as close to the center as possible.

    Parameters:
        - protein (Protein or ProteinSubsection): The protein or protein subsection object to be split.
        - length (dict): Dictionary containing the ideal, minimum, and maximum length values.
        - overlap (dict): Dictionary containing the ideal, minimum, and maximum overlap values.

    Returns:
        - tuple: A tuple containing two new ProteinSubsection objects if valid; otherwise, None.
    """
    # Error handling for invalid input types
    if not isinstance(protein, (Protein, ProteinSubsection)):
        raise ValueError("Input must be either a Protein or ProteinSubsection object.")

    # Determine the parent protein and relevant residue indices
    parent_protein = (protein.parent_protein 
                      if isinstance(protein, ProteinSubsection) 
                      else protein)
    
    first_res = protein.first_res
    last_res = protein.last_res
    domains = protein.domain_list

    # Calculate the initial midpoint
    midpoint = (first_res + last_res) // 2

    # Initialize variables for the cutpoints
    valid_cutpoint = None
    start_second_subsection = None

    # Function to check a given shift
    def check_shift(shift):
        nonlocal valid_cutpoint, start_second_subsection

        # Check if the shifted position is a valid cutpoint
        for direction in [1, -1]:  # 1 for forward, -1 for backward
            current_res = midpoint + (shift * direction)
            if check_valid_cutpoint(current_res, domains, last_res):
                start_second_subsection = find_next_start(current_res, protein, domains, overlap)
                if start_second_subsection is not None:
                    valid_cutpoint = current_res
                    return True
            # Check length constraint
            if direction == 1 and (last_res - current_res < length['min']):
                return False
            elif direction == -1 and (current_res - first_res < length['min']):
                return False
        return True

    # Search for a valid cutpoint near the midpoint
    shift = 0
    while check_shift(shift):
        shift += 1
        if valid_cutpoint is not None:
            break

    # If a valid cutpoint is found, create the two subsections
    if valid_cutpoint is not None:
        first_half = ProteinSubsection(
            parent_protein=parent_protein,
            start=first_res,
            end=valid_cutpoint
        )

        second_half = ProteinSubsection(
            parent_protein=parent_protein,
            start=start_second_subsection,
            end=last_res
        )

        return first_half, second_half

    return None