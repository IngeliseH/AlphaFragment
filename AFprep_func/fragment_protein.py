from AFprep_func.classes import Protein
from AFprep_func.recursive_fragment_finding import recursive_fragmentation
from AFprep_func.plot_fragments import plot_fragments

def fragment_protein(protein, min_len, max_len, overlap, min_overlap, max_overlap):
    if not isinstance(protein, Protein):
        raise ValueError("protein must be an instance of Protein.")
    fragments = None
    
    #recursive fragmentation
    while fragments == None:
        #deal with short proteins by classifying whole protein as one fragment
        #included in while loop so as max length increases this is still happening
        if len(protein.sequence) <= max_len:
            fragments = [(protein.first_res, protein.last_res)]
            continue

        fragments = recursive_fragmentation(protein.domain_list, protein.first_res, protein.first_res, protein.last_res, min_len, max_len, overlap, min_overlap, max_overlap)
        if fragments == None:
            #min_len = max (min_len-10, max_overlap + 1)
            #print(min_len)
            max_len = min(max_len+10, len(protein.sequence))
            print(max_len)
    if fragments is not None:
        return fragments
