from AFprep_func.classes import ProteinSubsection
from AFprep_func.fragmentation_methods import check_valid_cutpoint

def handle_long_domains(protein, min_len, max_len, overlap, min_overlap, max_overlap):
    # Nested function for adjusting adding overlap around long domains
    def add_overlap(position, direction, overlap, min_overlap, max_overlap):
        for adjusted_overlap in list(range(overlap, max_overlap + 1)) + list(range(overlap - 1, min_overlap - 1, -1)):
            adjusted_position = position + (direction * adjusted_overlap)
            if check_valid_cutpoint(adjusted_position, protein.domains, protein.last_res):
                return adjusted_position
        return position  # Return original if no adjustment is valid
    
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

        # Create subsection for the sequence between the end of the previous long domain and the start of the current one
        if prev_end < (start - 1):
            subsection_sequence_start = prev_end
            subsection_sequence_end = long_domain.start - 1  # Adjust for the start of the current long domain
            if subsection_sequence_end - subsection_sequence_start > 0:
                subsections.append(ProteinSubsection(protein, subsection_sequence_start, subsection_sequence_end))

        prev_end = end + 1

        # Adjust to add overlap if possible
        adjusted_start = add_overlap(start, overlap, min_overlap, max_overlap, protein.domain_list, len(protein.sequence), True)
        adjusted_end = add_overlap(end, overlap, min_overlap, max_overlap, protein.domain_list, len(protein.sequence), False)

        # Make sure that start and end are within boundaries of protein
        adjusted_start = max(0, adjusted_start)
        adjusted_end = min(protein.last_res, adjusted_end)

        # Add the domain with adjusted overlap as a fragment
        fragments.append((adjusted_start, adjusted_end))

    # Handle the tail end of the protein sequence, if any, after the last long domain
    if prev_end < len(protein.sequence):
        subsections.append(ProteinSubsection(protein, prev_end, len(protein.sequence)))

    return subsections, fragments

