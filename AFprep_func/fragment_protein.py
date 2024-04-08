from fragmentation_methods import recursive_fragmentation
from long_domains import handle_long_domains

def fragment_protein(protein, min_len = 150, max_len = 250, overlap = 10,
                     min_overlap = 0, max_overlap = 30):
    subsections, fragments = handle_long_domains(protein, min_len, max_len,
                                                 overlap, min_overlap, max_overlap)

    for subsection in subsections:
        subsection_fragments = None
        while subsection_fragments is None:
            #deal with short proteins/sections by classifying as one fragment
            #included in while loop so as max length increases this is still happening
            if len(subsection.sequence) <= max_len:
                subsection_fragments = [(subsection.first_res, subsection.last_res)]
                continue

            subsection_fragments = recursive_fragmentation(subsection, subsection.first_res,
                                                           min_len, max_len, overlap,
                                                           min_overlap, max_overlap)
            if subsection_fragments is None:
                max_len = min(max_len+10, len(subsection.sequence))
                print(max_len)
        fragments.extend(subsection_fragments)

    fragments.sort()

    return fragments
