from AFprep_func.classes import Protein
from AFprep_func.recursive_fragment_finding import recursive_fragmentation
from AFprep_func.plot_fragments import plot_fragments

def fragment_protein(protein, min_len = 150, max_len = 250, overlap = 10, max_overlap = 30, min_overlap = 0):
    if not isinstance(protein, Protein):
        raise ValueError("domain must be an instance of Domain.")
    fragments = None
    first_res = 0
    last_res = len(protein.sequence)
    #deal with short proteins by classifying whole protein as one fragment
    if len(protein.sequence) <= max_len:
        fragments = [(first_res, last_res)]

    #recursive fragmentation
    while fragments == None:
        fragments = recursive_fragmentation(protein.domain_list, first_res, first_res, last_res, min_len, max_len, overlap, max_overlap, min_overlap)
        if fragments == None:
            #min_len = max (min_len-10, max_overlap + 1)
            #print(min_len)
            max_len = min(max_len+10, len(protein.sequence))
            print(max_len)
    if fragments is not None:
        print("Cutpoints found for protein ", protein.name, ":", fragments)
        plot_fragments(protein, fragments)
        return fragments
