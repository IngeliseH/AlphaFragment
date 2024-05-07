=========================
Setting up Custom Classes
=========================

This section of the documentation outlines the custom classes used in the AlphaFragment library, specifically designed for managing proteins and their structural components such as domains and fragments. These classes are integral to the AFprep workflow.

.. automodule:: alphafragment.classes
   :members:
   :undoc-members:
   :show-inheritance:

Introduction
------------
The AlphaFragment library includes several classes that model various aspects of protein structures:

- **Domain**: Represents a domain within a protein sequence.
- **Protein**: Models the entire protein, including its sequence, domains, and fragments.
- **ProteinSubsection**: Represents a specific subsection of a protein sequence.

Classes and Their Descriptions
------------------------------

.. _class-domain:

Domain
^^^^^^
The ``Domain`` class models a specific domain within a protein sequence, characterized by start and end positions along with a type designation.

**Attributes:**
- ``id`` (str): Identifier for the domain.
- ``start`` (int): Start position of the domain in the sequence. Must be non-negative.
- ``end`` (int): End position of the domain in the sequence. Must be non-negative and greater than the start.
- ``domain_type`` (str): Type of the domain.

**Initialization:**
.. code-block:: python

    domain = Domain(identifier='D01', start=100, end=200, domain_type='helix')

**Exceptions:**
- Raises ``ValueError`` if start or end is out of valid range or if 'start' > 'end'.

.. _class-protein:

Protein
^^^^^^^
The ``Protein`` class models a complete protein and includes methods for managing its domains and fragments.

**Attributes:**
- ``name`` (str): Name of the protein.
- ``accession_id`` (str): UniProt accession ID.
- ``sequence`` (str): Amino acid sequence of the protein.
- ``first_res`` (int): Index of the first residue (default is 0).
- ``last_res`` (int): Index of the last residue (defaults to the sequence length).
- ``domain_list`` (list): List of ``Domain`` instances.
- ``fragment_list`` (list): List of fragments identified in the sequence.

**Methods:**
- ``add_domain(domain)``: Adds a domain to the protein.
- ``add_fragment(start, end)``: Adds a fragment to the protein's fragment list.

**Example:**
.. code-block:: python

    protein = Protein(name='Example', accession_id='P00001', sequence='MKT...')
    protein.add_domain(domain)
    protein.add_fragment(1, 100)

.. _class-protein-subsection:

ProteinSubsection
^^^^^^^^^^^^^^^^^
The ``ProteinSubsection`` class represents a subsection of a protein sequence, derived from a parent protein.

**Attributes:**
- ``parent_protein`` (Protein): The original protein from which the subsection is derived.
- ``start`` (int): Start position of the subsection in the parent protein's sequence.
- ``end`` (int): End position of the subsection in the parent protein's sequence.

**Initialization:**
.. code-block:: python

    subsection = ProteinSubsection(parent_protein=protein, start=100, end=200)

**Note:**
- Start and end positions should follow 0-based indexing and include the start and end positions.

Conclusion
----------
These classes provide a structured way to manage and analyze protein data within the AlphaFragment library, facilitating the AFprep workflow.