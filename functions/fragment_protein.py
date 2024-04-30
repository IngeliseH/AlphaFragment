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
  - .fragmentation_methods.validate_fragmentation_parameters: Validates
    the parameters used for protein fragmentation.
  - .fragmentation_methods.recursive_fragmentation: Used for recursively fragmenting
    protein sections that are not classified as long domains.
  - .fragmentation_methods.merge_overlapping_domains: Merges overlapping
    domains within a list of domains, and outputs as a new list
  - .long_domains.handle_long_domains: Handles the fragmentation of long domains
    within the protein.
"""

from .fragmentation_methods import validate_fragmentation_parameters, recursive_fragmentation, merge_overlapping_domains
from .long_domains import handle_long_domains

def fragment_protein(protein, min_len = 150, max_len = 250, overlap = None, len_increase = 10):
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
      - overlap (dict, optional): Dictionary containing the ideal, minimum, and
        maximum overlap values, in the format:
        {'min':min_overlap, 'ideal':ideal_overlap, 'max':max_overlap}
        where min_overlap, ideal_overlap and max_overlap are all integers, with
        min_overlap<=ideal_overlap<=max_overlap. Default is None, in which case
        the default values are used: {'min':0, 'ideal':10, 'max':30}
      - len_increase (int, optional): The amount by which to incrementally increase
        the maximum fragment length if a solution cannot be found. Default is 10.

    Returns:
      - list of tuples: A sorted list of tuples, where each tuple represents a
        fragment with its start and end positions within the protein sequence.
    """
    if not overlap:
        overlap = {'min':0, 'ideal':10, 'max':30}

    # Validate the input parameters
    validate_fragmentation_parameters(protein, min_len, max_len, overlap)

    subsections, fragments = handle_long_domains(protein, min_len, max_len, overlap)

    for subsection in subsections:
        merged_domains = merge_overlapping_domains(subsection.domain_list)
        subsection_fragments = None
        while subsection_fragments is None:
            # Deal with short proteins/sections by classifying as one fragment
            # included in while loop so as max length increases this is still happening
            if len(subsection.sequence) <= max_len:
                subsection_fragments = [(subsection.first_res, subsection.last_res+1)]
                continue

            subsection_fragments = recursive_fragmentation(subsection, merged_domains,
                                                           subsection.first_res,
                                                           min_len, max_len, overlap)
            if subsection_fragments is None:
                max_len = min(max_len + len_increase, len(subsection.sequence))
        fragments.extend(subsection_fragments)

    fragments.sort()
    print(protein.name, " is ", len(protein.sequence), " residues long and has ",
          len(fragments), "fragments", ":", fragments)

    return fragments
