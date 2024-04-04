from AFprep_func.classes import ProteinSubsection
from AFprep_func.recursive_fragment_finding import check_valid_cutpoint
from AFprep_func.fragment_protein import fragment_protein

def adjust_long_domain_overlap(position, overlap, min_overlap, max_overlap, domains, sequence_length, is_start):
    """
    Adjusts the overlap for a fragment's start or end position, ensuring the position does not split any domains.
    """
    valid_position = position - overlap if is_start else position + overlap
    if check_valid_cutpoint(valid_position, domains, sequence_length):
        return valid_position

    # If the desired overlap isn't valid, attempt to find a valid position within the allowed overlap range
    if is_start:  # Adjusting start position
        for adjusted_overlap in range(overlap + 1, max_overlap + 1):
            if check_valid_cutpoint(position - adjusted_overlap, domains, sequence_length):
                return position - adjusted_overlap
        for adjusted_overlap in range(overlap - 1, min_overlap - 1, -1):
            if check_valid_cutpoint(position - adjusted_overlap, domains, sequence_length):
                return position - adjusted_overlap
    else:  # Adjusting end position
        for adjusted_overlap in range(overlap + 1, max_overlap + 1):
            if check_valid_cutpoint(position + adjusted_overlap, domains, sequence_length):
                return position + adjusted_overlap
        for adjusted_overlap in range(overlap - 1, min_overlap - 1, -1):
            if check_valid_cutpoint(position + adjusted_overlap, domains, sequence_length):
                return position + adjusted_overlap

    # If no valid adjustment is found, return the original position
    return position


def handle_long_domains(protein, min_len, max_len, overlap, min_overlap, max_overlap):
    long_domain_list = []
    subsections = []
    fragments = []

    prev_end = 0
    for domain in protein.domain_list:
        domain_len = domain.end - domain.start + 1
        if domain_len >= max_len:
            long_domain_list.append(domain)

    # Process each long domain
    for long_domain in long_domain_list:
        # Check distances to protein ends and nearby long domains
        distance_to_start = long_domain.start
        distance_to_end = len(protein.sequence) - long_domain.end

        # Determine if the long domain should be merged with adjacent regions or other domains
        should_merge_start = distance_to_start < min_len
        should_merge_end = distance_to_end < min_len

        start = protein.first_res if should_merge_start else long_domain.start
        end = protein.last_res if should_merge_end else long_domain.end

        # Check for adjacency with other long domains
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
        
        # Create subsection for the sequence between the end of the previous long domain and the start of the current one
        if prev_end < (start - 1):
            subsection_sequence_start = prev_end
            subsection_sequence_end = long_domain.start - 1  # Adjust for the start of the current long domain
            if subsection_sequence_end - subsection_sequence_start > 0:
                subsections.append(ProteinSubsection(protein, subsection_sequence_start, subsection_sequence_end))
        
        prev_end = end + 1

        # Now adjust for overlap, ensuring we don't split any domains
        adjusted_start = adjust_long_domain_overlap(start, overlap, min_overlap, max_overlap, protein.domain_list, len(protein.sequence), True)
        adjusted_end = adjust_long_domain_overlap(end, overlap, min_overlap, max_overlap, protein.domain_list, len(protein.sequence), False)

        adjusted_start = max(0, adjusted_start)
        adjusted_end = min(protein.last_res, adjusted_end)

        # Add the domain with adjusted overlap as a fragment
        fragments.append((adjusted_start, adjusted_end))

    # Handle the tail end of the protein sequence, if any, after the last long domain
    if prev_end < len(protein.sequence):
        subsections.append(ProteinSubsection(protein, prev_end, len(protein.sequence)))

    return subsections, fragments

def produce_and_compile_fragments(protein, min_len = 150, max_len = 250, overlap = 10, min_overlap = 0, max_overlap = 30):
    subsections, fragments = handle_long_domains(protein, min_len, max_len, overlap, min_overlap, max_overlap)
    
    for subsection in subsections:
        print(subsection.first_res)
        subsection_fragments = fragment_protein(subsection, min_len, max_len, overlap, min_overlap, max_overlap)
        fragments.extend(subsection_fragments)
    
    fragments.sort()

    return fragments
