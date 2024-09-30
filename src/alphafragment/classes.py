"""
Defines classes for representing proteins and their structural components. This
includes handling domains, fragments, and subsections within protein sequences,
intended for use as part of the AlphaFragment workflow.

Classes:
    - Domain: Represents a domain within a protein sequence.
    - Protein: Models a protein, including sequence, domains, and fragments.
    - ProteinSubsection: Represents a subsection of a protein sequence.
"""

class Domain:
    """
    Represents a domain in a protein sequence, defined by start/end positions and a type.

    Attributes:
        - id (str): Identifier for the domain.
        - start (int): Start position of the domain in the sequence. Must be >= 0.
        - end (int): End position of the domain in the sequence. Must be >= 0.
        - domain_type (str): Type of the domain.

    Raises:
        - ValueError: If `start` or `end` is less than 0, or 'start' > 'end'.
    """
    def __init__(self, identifier, start, end, domain_type):
        """
        Initializes a new instance of the Domain class.
        """
        if start <0 or end < 0:
            raise ValueError("Domain start and end cannot be less than 0.")
        if start > end:
            raise ValueError("Domain start cannot be after end.")
        self.id = identifier
        self.start = start
        self.end = end
        self.type = domain_type

    def __str__(self):
        """
        Returns a string representation of the Domain instance, formatting it
        as 'domain_type + id' followed by the range '(start, end)'.
        """
        return f"{self.type}{self.id} ({self.start}, {self.end})"

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
        return (self.id == other.id and
                self.start == other.start and
                self.end == other.end and
                self.type == other.type)

    def __repr__(self):
        """
        Returns a formal string representation of the Domain instance.
        """
        return f"Domain(id={self.id}, start={self.start}, end={self.end}, domain_type='{self.type}')"

class Protein:
    """
    Models a protein, including sequence, domains and fragments. 
    Domains are significant structural or functional units within the protein,
    while fragments refer to specific subsequences of the protein sequence
    created as part of the AlphaFragment workflow.

    Attributes:
        - name (str): Name of the protein.
        - accession_id (str): UniProt accession ID.
        - sequence (str): Amino acid sequence of the protein.
        - first_res (int): Index of the first residue.
        - last_res (int): Index of the last residue, defaults to the sequence length.
        - domain_list (list of Domain instances, optional): Domains within the protein.
        - fragment_list (list of tuples, optional): Fragments identified in the
          protein sequence, represented as a tuple in the form (start_pos, end_pos),
          using pythonic slice notation (ie inclusive of start, exclusive of end).
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
        if last_res is not None:
            self.last_res = last_res
        elif last_res is None and sequence:
            self.last_res = len(sequence) - 1
        else:
            self.last_res = None
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

    def add_fragment(self, *args):
        """
        Adds a fragment to the protein's fragment list. Allows input as either 
        separate start and end integers or a tuple (start, end).
        Ensures that the fragment's start and end positions are positive integers
        and that the start is less than the end. Also checks for continuity with
        the last fragment added.

        Parameters:
            - args: Can be two integers (start, end) or a single tuple (start, end).

        Raises:
            - ValueError: If start or end are not positive integers or if start
              is not less than end. Also if the fragment does not
              follow sequentially after the last added fragment, or does not fit within
              sequence bounds
        """
        if len(args) == 1 and isinstance(args[0], tuple):
            start, end = args[0]
        elif len(args) == 2:
            start, end = args
        else:
            raise ValueError("Invalid arguments. Provide (start, end) as either two arguments or a tuple.")

        if not (isinstance(start, int) and isinstance(end, int) and 0 <= start < end):
            raise ValueError("Start and end must be positive integers, and start must be less than end.")

        if end > self.last_res + 1:
            raise ValueError("End of the new fragment must be within the protein sequence bounds.")
        if start < self.first_res:
            raise ValueError("Start of the new fragment must be within the protein sequence bounds.")

        if self.fragment_list and start < self.fragment_list[-1][0]:
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

    def __repr__(self):
        """
        Returns a formal string representation of the Protein instance.
        """
        domain_reprs = [repr(d) for d in self.domain_list]
        fragment_reprs = [repr(f) for f in self.fragment_list]
        sequence_repr = self.sequence[:10] + '...' if len(self.sequence) > 10 else self.sequence
        return (f"Protein(name={repr(self.name)}, accession_id={repr(self.accession_id)}, "
                f"sequence={repr(sequence_repr)}, "
                f"first_res={self.first_res}, last_res={self.last_res}, "
                f"domain_list={domain_reprs}, fragment_list={fragment_reprs})")


    def __eq__(self, other):
        """
        Checks if this Protein instance is equal to another by comparing their attributes.

        Parameters:
            - other (Protein): Another Protein instance to compare against.

        Returns:
            - bool: True if the proteins have the same attributes, False otherwise.
        """
        if not isinstance(other, Protein):
            return NotImplemented
        return (self.name == other.name and self.accession_id == other.accession_id and
                self.sequence == other.sequence and self.first_res == other.first_res and
                self.last_res == other.last_res and self.domain_list == other.domain_list and
                self.fragment_list == other.fragment_list)

class ProteinSubsection(Protein):
    """
    Represents a specific subsection of a protein sequence, inheriting from the
    Protein class. Initialized based on a parent protein and specified first_res/last_res
    positions within the parent's sequence. Takes sequence within the specified
    region and inherits all domains and fragments from the parent.

    Attributes:
        - parent_protein (Protein): The original protein from which the subsection is derived.
        - first_res (int): Start position of the subsection in the parent protein sequence.
        - last_res (int): End position of the subsection in the parent protein sequence.

    Note:
        - first_res and last_res are expected to be in 0-based indexing, and inclusive of the start and end.
    """
    def __init__(self, parent_protein, first_res, last_res):
        """
        Initializes a new ProteinSubsection instance, including all parent domains and fragments.
        Validates that the first_res and last_res indices are within the parent protein's sequence boundaries.
        """
        # Validate the boundaries
        if first_res < 0 or last_res > len(parent_protein.sequence) - 1 or first_res >= last_res:
            raise ValueError(
                f"Invalid first_res ({first_res}) or last_res ({last_res}) for the parent "
                f"protein sequence length {len(parent_protein.sequence)}."
                f"first_res must be less than last_res.")

        # Call the Protein constructor to initialize the subsection
        super().__init__(parent_protein.name,
                         parent_protein.accession_id,
                         parent_protein.sequence[first_res:last_res+1],  # Subset of the sequence
                         first_res,  # first_res of the subsection
                         last_res,   # last_res of the subsection
                         parent_protein.domain_list,  # Keep the same domain list
                         parent_protein.fragment_list)  # Keep the same fragment list

        # Store reference to the parent protein
        self.parent_protein = parent_protein
