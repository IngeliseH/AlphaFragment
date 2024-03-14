class Domain:
    def __init__(self, num, start, end, domain_type):
        """
        Initializes a Domain instance.

        Parameters:
        - start (int): The starting position of the domain in the protein sequence. Must be > 0.
        - num (string): A string identifying the domain
        - end (int): The ending position of the domain in the protein sequence. Must be > 0.
        - domain_type (str): The type or name of the domain.
        """
        if start <= 0 or end <= 0:
            raise ValueError("Start and end must be greater than 0.")
        self.start = start
        self.num = num
        self.end = end
        self.type = domain_type

    def __str__(self):
        # Format the string representation of the Domain instance
        return f"{self.type}{self.num} ({self.start}, {self.end})"

class Protein:
    def __init__(self, name, accession_id, sequence, domain_list=None, fragment_list=None):
        """
        Initializes a Protein instance.

        Parameters:
        - name (str): The name of the protein.
        - accession_id (str): The accession ID of the protein.
        - sequence (str): The amino acid sequence of the protein.
        - domain_list (list of Domain instances, optional): A list of domains within the protein. Default is None.
        - fragment_list (list of tuples, optional): A list of fragment tuples, where each tuple contains integers representing a fragment's start and end positions in the sequence. Default is None.
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
        """
        if not isinstance(domain, Domain):
            raise ValueError("domain must be an instance of Domain.")
        self.domain_list.append(domain)

    def add_fragment(self, start, end):
        """
        Adds a fragment to the protein's fragment list. Ensures that the fragment's start and end positions are positive integers and that the start is less than or equal to the end. Also checks for continuity with the last fragment added.

        Parameters:
        - start (int): The start position of the fragment in the protein sequence. Must be > 0.
        - end (int): The end position of the fragment in the protein sequence. Must be >= start.
        """
        if not (isinstance(start, int) and isinstance(end, int) and start > 0 and end >= start):
            raise ValueError("Start and end must be positive integers, and start must be less than or equal to end.")

        if self.fragment_list and start <= self.fragment_list[-1][-1]:
            raise ValueError("Start of the new fragment must be greater than the last value of the previous fragment.")

        self.fragment_list.append((start, end))

    def __str__(self):
        """
        Returns a string representation of the Protein instance, listing its name, accession ID, the count of domains, and the count of fragments.
        """
        return f"Protein Name: {self.name}, Accession ID: {self.accession_id}, Domains: {len(self.domain_list)}, Fragments: {len(self.fragment_list)}"