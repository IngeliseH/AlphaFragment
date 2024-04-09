"""
This module is designed to fragment a given protein into smaller, manageable
sections or fragments, ideal for use in AlphaFold predictions. It employs a
strategy that initially identifies and handles long domains within the protein
and then fragments the remaining protein sections. The module ensures that each
fragment is within specified size limits where possible, overlapping fragments
within specified boundaries.

Functions:
  - fragment_protein(protein, min_len, max_len, overlap, min_overlap, max_overlap):
    Fragments a given protein into smaller sections based on specified parameters.

Dependencies:
  - AFprep_func.fragmentation_methods.recursive_fragmentation: Used for recursively fragmenting
    protein sections that are not classified as long domains.
  - AFprep_func.long_domains.handle_long_domains: Handles the fragmentation of long domains
    within the protein.
"""

from AFprep_func.fragmentation_methods import recursive_fragmentation
from AFprep_func.long_domains import handle_long_domains

def fragment_protein(protein, min_len = 150, max_len = 250, overlap = 10,
                     min_overlap = 0, max_overlap = 30):
    """
    Fragments a given protein into smaller, manageable sections. Initially, it
    identifies long domains within the protein and organises fragments around
    these. Subsequently, it fragments the remaining protein sections. This
    function is designed to ensure each fragment falls within the specified size
    limits (where possible - long domains will always have fragments above the
    max length, and max length will be increased if a solution cannot be found).
    Protein fragments overlap within specified boundaries.

    Parameters:
      - protein (Protein): The protein object to be fragmented, containing domain
        and sequence information.
      - min_len (int): The minimum acceptable length for a protein fragment.
      - max_len (int): The initial maximum acceptable length for a protein
        fragment, adjusted dynamically for certain subsections.
      - overlap (int): The ideal size of the overlap between fragments.
      - min_overlap (int): The minimum allowed overlap size.
      - max_overlap (int): The maximum allowed overlap size.

    Returns:
      - list of tuples: A sorted list of tuples, where each tuple represents a
        fragment with its start and end positions within the protein sequence.
    """
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
        fragments.extend(subsection_fragments)

    fragments.sort()

    return fragments
