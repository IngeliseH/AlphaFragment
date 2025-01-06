"""
This module is designed to fragment a given protein into smaller, manageable
sections or fragments, ideal for use in AlphaFold predictions. It employs a
strategy that initially identifies and handles long domains within the protein
and then fragments the remaining protein sections. The module ensures that each
fragment is within specified size limits where possible, overlapping fragments
within specified boundaries.

Functions:
  - fragment_protein: Fragments a given protein into smaller sections based on
    specified parameters.

Dependencies:
  - time: Used for setting a time limit for recursive fragmentation.
  - .fragmentation_methods.validate_fragmentation_parameters: Validates
    the parameters used for protein fragmentation.
  - .fragmentation_methods.recursive_fragmentation: Used for recursively fragmenting
    protein sections that are not classified as long domains.
  - .fragmentation_methods.merge_overlapping_domains: Merges overlapping
    domains within a list of domains, and outputs as a new list
  - .long_domains.handle_long_domains: Handles the fragmentation of long domains
    within the protein.
  - .fragmentation_methods.break_in_half: Splits a protein or subsection in half
"""

def fragment_protein(protein, length=None, overlap=None, len_increase=10, time_limit=0.1):
    """
    Fragments a given protein into smaller, manageable sections. Initially, it
    identifies long domains within the protein and organizes fragments around
    these. Subsequently, it fragments the remaining protein sections. This
    function is designed to ensure each fragment falls within the specified size
    limits (where possible - long domains will always have fragments above the
    max length, and max length will be increased if a solution cannot be found).
    Protein fragments overlap within specified boundaries. Uses a time limit for
    recursive fragmentation to prevent infinite loops.

    Parameters:
      - protein (Protein): The protein object to be fragmented, containing domain
        and sequence information.
      - length (dict, optional): Dictionary containing the ideal, minimum, and
        maximum length values, in the format:
        {'min': min_len, 'ideal': ideal_len, 'max': max_len}
        where min_len, ideal_len, and max_len are all integers,
        with min_len <= ideal_len <= max_len. Default is None, in which case
        the default values are used
      - overlap (dict, optional): Dictionary containing the ideal, minimum, and
        maximum overlap values, in the format:
        {'min': min_overlap, 'ideal': ideal_overlap, 'max': max_overlap}
        where min_overlap, ideal_overlap and max_overlap are all integers, with
        min_overlap <= ideal_overlap <= max_overlap. Default is None, in which case
        the default values are used
      - len_increase (int, optional): The amount by which to incrementally increase
        the maximum fragment length if a solution cannot be found. Default is 10.
      - time_limit (float, optional): The maximum time allowed for recursive
        fragmentation, in seconds.

    Returns:
      - list of tuples: A sorted list of tuples, where each tuple represents a
        fragment with its start and end positions within the protein sequence.
    """
    import time
    from .long_domains import handle_long_domains
    from .fragmentation_methods import validate_fragmentation_parameters, recursive_fragmentation, merge_overlapping_domains, break_in_half

    if not length:
        length = {'min': 200, 'ideal': 384, 'max': 400}
    if not overlap:
        overlap = {'min': 20, 'ideal': 30, 'max': 40}

    # Validate the input parameters
    validate_fragmentation_parameters(protein, length, overlap)

    subsections, fragments = handle_long_domains(protein, length, overlap)

    while subsections:  # Continue processing until all subsections are fragmented
        next_subsections = []
        for subsection in subsections:
            merged_domains = merge_overlapping_domains(subsection.domain_list)
            subsection_fragments = None
            max_len = length['max']
            original_max_len = length['max']

            while subsection_fragments is None:
                if len(subsection.sequence) <= max_len:
                    subsection_fragments = [(subsection.first_res, subsection.last_res + 1)]
                    break

                # Start recursive fragmentation with a time limit
                start_time = time.time()
                subsection_fragments = recursive_fragmentation(subsection, merged_domains,
                                                               subsection.first_res,
                                                               length, overlap, original_max_len,
                                                               time_limit=time_limit,
                                                               start_time=start_time)

                if subsection_fragments == "TIME_LIMIT_EXCEEDED":
                    # If the time limit is exceeded, split protein + add new subsections to the list
                    subsections_to_process = break_in_half(subsection, length, overlap)
                    if subsections_to_process is None:
                        subsection_fragments = [(subsection.first_res, subsection.last_res + 1)]
                        print("Recommended to increase time limit. Fragmentation ineffective.")
                        break
                    next_subsections.extend(subsections_to_process)
                    break  # Exit the while loop to process the new subsections
                if subsection_fragments is None:
                    # If no valid fragmentation pattern was found, increase max length and try again
                    max_len = min(max_len + len_increase, len(subsection.sequence))
                    length['max'] = max_len
                else:
                    # Successfully found fragments, exit the loop
                    break

            # If the subsection was successfully fragmented, add the fragments to the list
            if subsection_fragments and subsection_fragments != "TIME_LIMIT_EXCEEDED":
                fragments.extend(subsection_fragments)
        if max_len > original_max_len:
            print(f"Max length increased to {max_len} for section of {protein.name}")

        # Update the list of subsections for the next round of processing
        subsections = next_subsections

    fragments.sort()
    print(protein.name, " is ", len(protein.sequence), " residues long and has ",
          len(fragments), "fragments", ":", fragments)

    return fragments
