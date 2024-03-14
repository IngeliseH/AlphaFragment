class Domain:
    def __init__(self, start, end):
        """
        Initializes a Domain instance.

        Parameters:
        - start (int): The starting position of the domain in the protein sequence. Must be > 0.
        - end (int): The ending position of the domain in the protein sequence. Must be > 0.
        """
        if start <= 0 or end <= 0:
            raise ValueError("Start and end must be greater than 0.")
        self.start = start
        self.end = end

    def __str__(self):
        # Format the string representation of the Domain instance
        return f"({self.start}, {self.end})"
    

def recursive_fragmentation(fragment_start, first_res, last_res, min_len, max_len, overlap, max_overlap, min_overlap, cutpoints=None, overlaps=None):
    # Base case: if previous fragment end is at or beyond the last_res, we've reached the end and return the cutpoints
    if cutpoints:
        if cutpoints[-1][1] >= last_res:
            return cutpoints
    
    if cutpoints is None:
        cutpoints = []
    if overlaps is None:
        overlaps = []

    valid_cut_found = False

    for res in range(fragment_start + min_len, min(fragment_start + max_len, last_res) + 1):
        if check_valid_cutpoint(res, domains, last_res):
            next_start = find_next_start(res, overlap, max_overlap, min_overlap, domains, last_res)
            if next_start:
                valid_cut_found = True
                cutpoints.append((fragment_start, res))
                overlaps.append(res - next_start)
                
                # Recursively process the next segment
                result = recursive_fragmentation(next_start, first_res, last_res, min_len, max_len, overlap, max_overlap, min_overlap, cutpoints)
                
                # If a valid fragmentation pattern is found, return the result
                if result is not None:
                    return result
                else:
                    # If the recursive call did not find a valid pattern, remove the last added cutpoints and overlaps
                    cutpoints.pop()
                    overlaps.pop()

    # If no valid cut is found in the loop, check if this is the first call and return None to indicate failure
    if not valid_cut_found and fragment_start == first_res:
        return None
    # If this is a recursive call, return None to backtrack
    elif not valid_cut_found:
        return None

def check_valid_cutpoint(res, domains, sequence_end):
    """
    Checks if a cutpoint is valid based on the specified conditions.
    
    Parameters:
    - res (int): The position to check for a valid cutpoint.
    - domains (list of Domain): The list of domain objects.
    - sequence_end (int): The end position of the protein sequence.
    
    Returns:
    - bool: True if the cutpoint is valid, False otherwise.
    """
    # Check if res is beyond the end of the sequence
    if res > sequence_end + 1:
        return False
    
    for domain in domains:
        # Check if res and res-1 are within the same domain
        if domain.start <= res <= domain.end and domain.start <= res-1 <= domain.end:
            return False
            
    return True


def find_next_start(res, overlap, max_overlap, min_overlap, domains, last_res):
    if check_valid_cutpoint(res - overlap, domains, last_res):
        return res - overlap
    for increased_overlap in range(overlap + 1, max_overlap):
        if check_valid_cutpoint(res - increased_overlap, domains, last_res):
            return res - increased_overlap
    for decreased_overlap in range(overlap - 1, min_overlap, -1):
        if check_valid_cutpoint(res - decreased_overlap, domains, last_res):
            return res - decreased_overlap
    return None

# Example call (values for first_res, last_res, min_len, max_len need to be defined)
first_res = 1
last_res = 1000  # Example end of sequence
domains = [
    Domain(1, 190),
    Domain(400, 430),
    Domain(435, 500),
    Domain(300, 350),
    Domain(565, 700)
]
min_len = 100
max_len = 200
overlap = 10
max_overlap = 20
min_overlap = 0
result = recursive_fragmentation(first_res, first_res, last_res, min_len, max_len, overlap, max_overlap, min_overlap)

if result is not None:
    print("Cutpoints found:", result)
else:
    print("No valid fragmentation pattern found.")