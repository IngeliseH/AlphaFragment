"""
Defines classes for representing proteins and their structural components. This
includes handling domains, fragments, and subsections within protein sequences,
intended for use as part of the AFprep workflow.

Classes:
    - Domain: Represents a domain within a protein sequence.
    - Protein: Models a protein, including sequence, domains, and fragments.
    - ProteinSubsection: Represents a subsection of a protein sequence.
"""

class Domain:
    """
    Represents a domain in a protein sequence, defined by start/end positions and a type.

    Attributes:
        - num (str): Identifier for the domain.
        - start (int): Start position of the domain in the sequence. Must be >= 0.
        - end (int): End position of the domain in the sequence. Must be >= 0.
        - domain_type (str): Type of the domain.

    Raises:
        - ValueError: If `start` or `end` is less than 0, or 'start' > 'end'.
    """
    def __init__(self, num, start, end, domain_type):
        """
        Initializes a new instance of the Domain class.
        """
        if start <0 or end < 0:
            raise ValueError("Domain start and end must be greater than 0.")
        if start > end:
            raise ValueError("Domain start cannot be after end.")
        if start == end:
            raise ValueError("Domain must be longer than 1 residue.")
        self.num = num
        self.start = start
        self.end = end
        self.type = domain_type

    def __str__(self):
        """
        Returns a string representation of the Domain instance, formatting it
        as 'domain_type + num' followed by the range '(start, end)'.
        """
        return f"{self.type}{self.num} ({self.start}, {self.end})"

    def __eq__(self, other):
        """
        Checks if this Domain instance is equal to another by comparing their attributes.

        Parameters:
            - other (Domain): Another Domain instance to compare against.

        Returns:
            - bool: True if the domains have the same attributes, False otherwise.
        """
        if not isinstance(other, Domain):
            return NotImplemented
        return (self.num == other.num and
                self.start == other.start and
                self.end == other.end and
                self.type == other.type)

    def __repr__(self):
        """
        Returns a formal string representation of the Domain instance.
        """
        return f"Domain(num={self.num}, start={self.start}, end={self.end}, domain_type='{self.type}')"

class Protein:
    """
    Models a protein, including sequence, domains and fragments. 
    Domains are significant structural or functional units within the protein,
    while fragments refer to specific subsequences of the protein sequence
    created as part of the AFprep workflow.

    Attributes:
        - name (str): Name of the protein.
        - accession_id (str): UniProt accession ID.
        - sequence (str): Amino acid sequence of the protein.
        - first_res (int): Index of the first residue.
        - last_res (int): Index of the last residue, defaults to the sequence length.
        - domain_list (list of Domain instances, optional): Domains within the protein.
        - fragment_list (list of tuples, optional): Fragments identified in the
          protein sequence, represented as a tuple in the form (start_pos, end_pos)
    """
    def __init__(self, name, accession_id, sequence, first_res=0, last_res=None,
                domain_list=None, fragment_list=None):
        """
        Initializes a new instance of the Protein class.
        """
        self.name = name
        self.accession_id = accession_id
        self.sequence = sequence
        self.first_res = first_res
        self.last_res = last_res if last_res is not None else len(sequence)
        self.domain_list = domain_list if domain_list is not None else []
        self.fragment_list = fragment_list if fragment_list is not None else []

    def add_domain(self, domain):
        """
        Adds a Domain instance to the protein's domain list.

        Parameters:
            - domain (Domain): The Domain instance to be added.

        Raises:
            - ValueError: If the input is not an instance of the Domain class.
        """
        if not isinstance(domain, Domain):
            raise ValueError("domain must be an instance of Domain.")
        self.domain_list.append(domain)

    def add_fragment(self, start, end):
        """
        Adds a fragment to the protein's fragment list. Ensures that the
        fragment's start and end positions are positive integers and that the
        start is less than the end. Also checks for continuity with
        the last fragment added.

        Parameters:
            - start (int): The start position of the fragment in the protein sequence.
              Must be > 0.
            - end (int): The end position of the fragment in the protein sequence.
              Must be > start.

        Raises:
            - ValueError: If start or end are not positive integers or if start
              is not less than end. Also raises ValueError if the fragment does not
              follow sequentially after the last added fragment.
        """
        if not (isinstance(start, int) and isinstance(end, int) and 0 <= start < end):
            raise ValueError("Start and end must be positive integers, and start "
                             "must be less than end.")

        if self.fragment_list and start <= self.fragment_list[-1][0]:
            raise ValueError("Start of the new fragment must be greater than the "
                             "start of the previous fragment.")

        self.fragment_list.append((start, end))

    def __str__(self):
        """
        Returns a string representation of the Protein instance, listing its
        name, accession ID, and the number of associated domains and fragments.
        """
        return (f"Protein Name: {self.name}, Accession ID: {self.accession_id}, "
                f"Domains: {len(self.domain_list)}, Fragments: {len(self.fragment_list)}")

class ProteinSubsection(Protein):
    """
    Represents a specific subsection of a protein sequence, inheriting from the
    Protein class. Initialized based on a parent protein and specified start/end
    positions within the parent's sequence. Takes sequence within the specified
    region and all domains and fragments from parent.

    Attributes:
        - parent_protein (Protein): The original protein from which the subsection is derived.
        - start (int): Start position of the subsection in the parent protein sequence.
        - end (int): End position of the subsection in the parent protein sequence.
    """
    def __init__(self, parent_protein, start, end):
        """
        Initializes a new ProteinSubsection instance, including all parent domains and fragments.
        Validates that the start and end indices are within the parent protein's sequence boundaries.
        """
        if start < 0 or end > len(parent_protein.sequence) or start >= end:
            raise ValueError(f"Invalid start ({start}) or end ({end}) for the parent protein sequence length {len(parent_protein.sequence)}. Start must be less than end.")

        super().__init__(parent_protein.name,
                         parent_protein.accession_id,
                         parent_protein.sequence[start:end],
                         start,
                         end,
                         parent_protein.domain_list,
                         parent_protein.fragment_list)
        self.parent_protein = parent_protein
