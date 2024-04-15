"""
This module provides functionality for handling long domains within proteins in
the context of protein fragmentation. It identifies domains exceeding a
specified maximum length and generates appropriate fragments from these long
domains. The module also identifies and outputs unfragmented regions of the
protein around these long domains for further processing.

Functions:
  - handle_long_domains(protein, min_len, max_len, overlap, min_overlap,
    max_overlap): Main function to handle long domains in proteins, creating
    fragments and identifying adjacent unfragmented regions which are output as a
    list of ProteinSubsection objects.

Dependencies:
  - AFprep_func.classes.ProteinSubsection: Used to represent sections of a protein,
    including domains and fragments.
  - AFprep_func.classes.Domain: Used to represent protein domains.
  - AFprep_func.fragmentation_methods.check_valid_cutpoint: Utilized to ensure
    proposed fragmentation points are valid based on domain boundaries and protein.
"""

from AFprep_func.classes import ProteinSubsection, Domain
from AFprep_func.fragmentation_methods import check_valid_cutpoint

def handle_long_domains(protein, min_len, max_len, overlap, min_overlap, max_overlap):
    """
    Identifies long domains within a protein and generates fragments that
    include these domains, along with the unfragmented regions of the protein
    surrounding these domains. It takes into consideration addition of overlap.
    This function is designed to work as part of a protein fragmentation workflow.

    Parameters:
      - protein (Protein): The protein object containing domain and sequence
        information.
      - min_len (int): The minimum acceptable length for a protein fragment.
      - max_len (int): The maximum acceptable standard length for a protein
        fragment, beyond which a domain is considered long.
      - overlap (int): The initial, 'ideal' overlap size to try when adjusting the
        start and end of a fragment around a long domain.
      - min_overlap (int): The minimum allowed overlap size.
      - max_overlap (int): The maximum allowed overlap size.

    Returns:
      - tuple: A tuple containing two lists: (1) unfragmented subsections of the
        protein surrounding the long domains, and (2) fragments that include the
        long domains with appropriate overlaps.
    
    Function Logic:
      - Merge overlapping domains into a single domain to simplify processing, and add
        these to a new list.
      - Iterate over the list of domains in the protein to identify those exceeding
        the specified maximum length (`max_len`), categorizing them as long domains.
      - For each identified long domain:
          - Determine if it should be merged with sequence between this and adjacent
            long domains, or protein ends
              - If distance to these points is less than the minimum fragment length
                (`min_len`), merge to avoid leaving very short fragments. If a region
                between two long domains is < min_len, merge this with the shorter
                of the two long domains.
          - For the region before the current long domain, create a subsection if it
            hasn't been included in a previous fragment or subsection and if its
            length is viable.
      - Adjust the start and end points of each long domain fragment to include
        overlaps:
          - Attempt to add an overlap to both the start and end points within the
            bounds of `min_overlap` and `max_overlap`, ensuring the selected points
            are valid cutpoints.
          - Ensure the adjusted start and end points do not extend beyond the protein's
            boundaries.
      - Add the long domain, neighbouring merged regions, and overlap, as a fragment.
      - After processing all long domains, check for any remaining unfragmented
        protein sequence at the end and create a subsection if necessary.
      - Return the created subsections and fragments.
    """
    # Nested function for adjusting adding overlap around long domains
    def add_overlap(position, direction):
        for adjusted_overlap in (list(range(overlap, max_overlap + 1)) +
                                 list(range(overlap - 1, min_overlap - 1, -1))):
            adjusted_position = position + (direction * adjusted_overlap)
            if check_valid_cutpoint(adjusted_position, protein.domain_list, protein.last_res):
                return adjusted_position
        return position  # Return original if no adjustment is valid

    long_domain_list = []
    subsections = []
    fragments = []

    # Merge overlapping domains to simplify processing
    sorted_domains = sorted(protein.domain_list, key=lambda x: x.start)
    combined_domains = []
    for domain in sorted_domains:
        if not combined_domains:
            combined_domains.append(domain)
        else:
            last = combined_domains[-1]
            # Check if the current domain overlaps with the last one in the list
            if domain.start <= last.end:
                # Merge the two domains by updating the end of the last domain
                combined_domains[-1] = Domain(last.num, last.start, max(last.end, domain.end), last.type)
            else:
                combined_domains.append(domain)

    prev_end = 0
    for domain in protein.domain_list:
        domain_len = domain.end - domain.start + 1
        if domain_len >= max_len:
            long_domain_list.append(domain)

    # Process each long domain
    for long_domain in long_domain_list:
        # Check distances to protein ends and nearby long domains
        distance_to_start = long_domain.start
        distance_to_end = protein.last_res - long_domain.end

        # Determine if the long domain should be extended to include the gap to
        # other long domains, or to the start of end of the protein
        should_merge_start = distance_to_start < min_len
        should_merge_end = distance_to_end < min_len

        start = protein.first_res if should_merge_start else long_domain.start
        end = protein.last_res if should_merge_end else long_domain.end

        # Check for distance to other long domains
        long_domain_len = long_domain.end - long_domain.start + 1
        for other_domain in long_domain_list:
            other_domain_len = other_domain.end - other_domain.start + 1
            if other_domain == long_domain:
                continue
            if abs(long_domain.end - other_domain.start) < min_len:
                if long_domain_len <= other_domain_len:
                    end = max(end, other_domain.start)
            if abs(long_domain.start - other_domain.end) < min_len:
                if long_domain_len < other_domain_len:
                    start = min(start, other_domain.end)

        # Create subsection for sequence between end of previous long domain and
        # start of current one
        if prev_end < (start - 1):
            subsection_sequence_start = prev_end
            subsection_sequence_end = long_domain.start - 1
            if subsection_sequence_end - subsection_sequence_start > 0:
                subsections.append(ProteinSubsection(protein,
                                                     subsection_sequence_start,
                                                     subsection_sequence_end))

        prev_end = end + 1

        # Adjust to add overlap if possible
        adjusted_start = add_overlap(start, -1)
        adjusted_end = add_overlap(end, 1)

        # Make sure that start and end are within boundaries of protein
        adjusted_start = max(0, adjusted_start)
        adjusted_end = min(protein.last_res, adjusted_end)

        # Add the domain with adjusted overlap as a fragment
        fragments.append((adjusted_start, adjusted_end))

    # Handle the tail end of the protein sequence, if any, after the last long domain
    if prev_end < len(protein.sequence):
        subsections.append(ProteinSubsection(protein, prev_end, len(protein.sequence)))

    return subsections, fragments
