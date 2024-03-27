"""
This module defines classes for representing and manipulating protein sequences,
focusing on the identification and organization of domains and fragments within
these sequences. It provides a structured approach to model proteins, their
distinct domains, and specific subsequences or fragments.

Classes:
- Domain: Represents a specific domain within a protein sequence, characterized
  by a unique identifier, start and end positions, and a domain type.

- Protein: Models a protein with its amino acid sequence, including information
  about its name, UniProt accession ID, and optionally, lists of its domains and
  fragments. This class allows for the addition of domains and fragments to
  a protein.

Usage:
This module is intended for use in bioinformatics applications, and
particularly as part of the workflow for domain aware fragmentation of protein
sequences for AlphaFold input. It enables the organization of protein data into
structured objects, making it easier to work with complex protein information,
such as domain architecture and sequence fragments.

Dependencies:
The classes within this module do not rely on external libraries for their
basic functionality.

Note:
This module does not include functionality for automatic domain prediction or
sequence retrieval from databases. This module is designed primarily for use as
part of the AFprep workflow.
"""

class Domain:
    """
    Represents a domain within a protein sequence, identified by a unique string,
    with specified start and end positions, and the domain's type.

    A domain is a distinct part of a protein sequence. Each domain in this
    context is defined by its starting and ending positions in the protein
    sequence, along with a type or name that categorizes the domain.

    Parameters:
    - num (str): A string identifying the domain.
    - start (int): The start position of the domain in the protein sequence. Must
      be >= 0.
    - end (int): The end position of the domain in the protein sequence. Must be
      >= 0.
    - domain_type (str): The type of the domain.

    Raises:
    - ValueError: If `start` or `end` is less than 0.
    """
    def __init__(self, num, start, end, domain_type):
        """
        Initializes a new instance of the Domain class.
        """
        if start <0 or end < 0:
            raise ValueError("Start and end must be greater than 0.")
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

class Protein:
    """
    Represents a protein with a specified amino acid sequence, including
    information on its name, uniprot accession ID, and optionally, domains and
    fragments within the protein sequence.

    Domains are significant structural or functional units within the protein,
    while fragments refer to specific subsequences of the protein sequence.

    Parameters:
    - name (str): The name of the protein.
    - accession_id (str): The uniprot accession ID of the protein.
    - sequence (str): The amino acid sequence of the protein.
    - domain_list (list of Domain instances, optional): A list representing
      the domains identified within the protein sequence. Defaults to an empty
      list.
    - fragment_list (list of tuples, optional): A list of tuples, each
      representing a fragment's start and end positions within the protein
      sequence. Defaults to an empty list.
    """

    def __init__(self, name, accession_id, sequence, domain_list=None, fragment_list=None):
        """
        Initializes a new instance of the Protein class.
        """
        self.name = name
        self.accession_id = accession_id
        self.sequence = sequence
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
        if not (isinstance(start, int) and isinstance(end, int) and 0 < start < end):
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
