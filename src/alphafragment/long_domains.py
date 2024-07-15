"""
This module provides functionality for handling long domains within proteins in
the context of protein fragmentation. It identifies domains exceeding a
specified maximum length and generates appropriate fragments from these long
domains. The module also identifies and outputs unfragmented regions of the
protein around these long domains for further processing.

Functions:
  - handle_long_domains: Main function to handle long domains in proteins, creating
    fragments and identifying adjacent unfragmented regions which are output as a
    list of ProteinSubsection objects.

Dependencies:
  - .classes.ProteinSubsection: Used to represent sections of a protein,
    including domains and fragments.
  - .fragmentation_methods.validate_fragmentation_parameters: Used to
    validate the input parameters for protein fragmentation.
  - .fragmentation_methods.check_valid_cutpoint: Utilized to ensure
    proposed fragmentation points are valid based on domain boundaries and protein.
  - A.fragmentation_methods.merge_overlapping_domains: Used to merge
    overlapping domains within a protein.
"""

from .classes import ProteinSubsection
from .fragmentation_methods import validate_fragmentation_parameters, check_valid_cutpoint, merge_overlapping_domains

def handle_long_domains(protein, length, overlap):
    """
    Identifies long domains within a protein and generates fragments that
    include these domains, along with the unfragmented regions of the protein
    surrounding these domains. It takes into consideration addition of overlap.
    This function is designed to work as part of a protein fragmentation workflow.

    Parameters:
      - protein (Protein): The protein object containing domain and sequence
        information.
      - length (dict): Dictionary containing the ideal, minimum, and maximum
        length values, in the format:
        {'min': min_len, 'ideal': ideal_len, 'max': max_len}
        where min_len, ideal_len, and max_len are all integers,
        with min_len <= ideal_len <= max_len.
      - overlap (dict): Dictionary containing the ideal, minimum, and maximum
        overlap values, in the format:
        {'min': min_overlap, 'ideal': ideal_overlap, 'max': max_overlap}
        where min_overlap, ideal_overlap and max_overlap are all integers, with
        min_overlap <= ideal_overlap <= max_overlap.

    Returns:
      - tuple: A tuple containing two lists: (1) unfragmented subsections of the
        protein surrounding the long domains, and (2) fragments that include the
        long domains with appropriate overlaps.

    Function Logic:
      - Merge overlapping domains into a single domain to simplify processing, and add
        these to a new list.
      - Iterate over the list of domains in the protein to identify those exceeding
        the specified maximum length (`length['max']`), categorizing them as long domains.
      - For each identified long domain:
          - Determine if it should be merged with sequence between this and adjacent
            long domains, or protein ends
              - If distance to these points is less than the minimum fragment length
                (`length['min']`), merge to avoid leaving very short fragments. If a region
                between two long domains is < length['min'], merge this with the shorter
                of the two long domains.
          - For the region before the current long domain, create a subsection if it
            hasn't been included in a previous fragment or subsection and if its
            length is viable.
      - Adjust the start and end points of each long domain fragment to include
        overlaps:
          - Attempt to add an overlap to both the start and end points within the
            bounds of 'overlap['max']' and 0, ensuring the selected points are valid
            cutpoints.
          - If an overlap cannot be added, the original start and end points are
            used.
          - Ensure the adjusted start and end points do not extend beyond the protein's
            boundaries.
      - Add the long domain, neighbouring merged regions, and overlap, as a fragment.
      - After processing all long domains, check for any remaining unfragmented
        protein sequence at the end and create a subsection if necessary.
      - Return the created subsections and fragments.

    Notes:
      - If two long domains are adjacent and the distance between them is less than
        'overlap['min']', the domains will still be created as separate fragments with
        as much overlap as is allowed by the space between them.
    """
    # Validate the input parameters
    validate_fragmentation_parameters(protein, length, overlap)

    # Nested function for adjusting adding overlap around long domains
    def add_overlap(position, direction):
        for adjusted_overlap in (list(range(overlap['ideal'], overlap['max'] + 1)) +
                                 list(range(overlap['ideal'] - 1, - 1, -1))):
            adjusted_position = position + (direction * adjusted_overlap)
            if check_valid_cutpoint(adjusted_position, protein.domain_list, protein.last_res):
                return adjusted_position
        return position  # Return original if no adjustment is valid (should never be used as overlap will be reduced to 0 which will return original anyway)

    long_domain_list = []
    subsections = []
    fragments = []

    # Merge overlapping domains to simplify processing
    combined_domains = merge_overlapping_domains(protein.domain_list)

    prev_end = 0
    for domain in combined_domains:
        domain_len = domain.end - domain.start + 1
        if domain_len >= length['max']:
            long_domain_list.append(domain)

    # Process each long domain
    for long_domain in long_domain_list:
        # Check distances to protein ends
        distance_to_start = long_domain.start
        distance_to_end = protein.last_res - long_domain.end

        # Determine if the long domain should be extended to the start or end of the protein
        should_merge_start = distance_to_start < length['min']
        should_merge_end = distance_to_end < length['min']

        start = protein.first_res if should_merge_start else long_domain.start
        end = protein.last_res if should_merge_end else long_domain.end

        # Check for distance to other long domains
        long_domain_len = long_domain.end - long_domain.start + 1
        for other_domain in long_domain_list:
            other_domain_len = other_domain.end - other_domain.start + 1
            if other_domain == long_domain:
                continue
            if (long_domain.end < other_domain.start) and (other_domain.start - long_domain.end - 1 < length['min']):
                if long_domain_len <= other_domain_len:
                    end = max(end, other_domain.start)
            if (long_domain.start > other_domain.end) and (long_domain.start - other_domain.end - 1 < length['min']):
                if long_domain_len < other_domain_len:
                    start = min(start, other_domain.end)

        # Create subsection for sequence between end of previous long domain and
        # start of current one
        # Subsections are inclusive of both start and end
        if prev_end < (start - 1):
            subsection_sequence_start = prev_end + 1
            subsection_sequence_end = long_domain.start - 1
            subsections.append(ProteinSubsection(protein,
                                                 subsection_sequence_start,
                                                 subsection_sequence_end))

        prev_end = end

        # Adjust to add overlap if possible
        adjusted_start = add_overlap(start, -1)
        adjusted_end = add_overlap(end, 1)

        # Make sure that start and end are within boundaries of protein
        adjusted_start = max(0, adjusted_start)
        adjusted_end = min(protein.last_res, adjusted_end)

        # Add the domain with adjusted overlap as a fragment - converting to slicing indexing
        fragments.append((adjusted_start, adjusted_end + 1))

    # Handle the tail end of the protein sequence, if any, after the last long domain
    # If no long domains were found, return whole protein
    if prev_end == 0:
        # No long domains were found, so the entire protein is a subsection
        subsections.append(protein)
    elif prev_end < protein.last_res:
        subsections.append(ProteinSubsection(protein, prev_end + 1, protein.last_res))

    return subsections, fragments
