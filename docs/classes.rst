=========================
Setting up Custom Classes
=========================

This section explains the classes used by AlphaFragment.

Introduction
------------
AlphaFragment uses custom classes to represent proteins, domains, and protein
subsections.

Detailed Class Information
--------------------------

Domain Class
^^^^^^^^^^^^

The ``Domain`` class represents a specific domain within a protein sequence.

.. autoclass:: alphafragment.classes.Domain
   :members:
   :undoc-members:
   :show-inheritance:

   .. automethod:: __init__

Protein Class
^^^^^^^^^^^^^

The ``Protein`` class models a complete protein including its sequence, domains, and fragments.

.. autoclass:: alphafragment.classes.Protein
   :members:
   :undoc-members:
   :show-inheritance:

   Methods
   +++++++

   .. automethod:: add_domain
       :noindex:

   .. automethod:: add_fragment
       :noindex:

ProteinSubsection Class
^^^^^^^^^^^^^^^^^^^^^^^

The ``ProteinSubsection`` class represents a specific subsection of a protein sequence, derived from a parent protein.

.. autoclass:: alphafragment.classes.ProteinSubsection
   :members:
   :undoc-members:
   :show-inheritance:

   .. automethod:: __init__

Examples
--------

**Domain Initialization Example:**

.. code-block:: python

    # Initializing a domain within a protein sequence
    domain = Domain(identifier='D01', start=100, end=200, domain_type='helix')

**Protein and Subsection Usage:**

.. code-block:: python

    # Creating a protein and adding a domain
    protein = Protein(name='Example', accession_id='P00001', sequence='MKT...')
    protein.add_domain(domain)
    protein.add_fragment(1, 100)

**Note:** Most of these classes are used internally by AlphaFragment and are
not intended to be used directly by the user, past basic initialization.