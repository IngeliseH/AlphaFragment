============
Introduction
============

AlphaFragment: Protein Fragmentation and Pair Generation for AlphaFold and AlphaPulldown
--------------------------------------------------------------------------------------

**AlphaFragment** is a Python package created to improve the preparation of
protein data for structural analysis using tools like AlphaFold and
AlphaPulldown. The package equips users with efficient methods to import
protein data, accurately identify protein domains using information from
AlphaFoldDB and UniProt, and perform protein fragmentation while preserving
domain integrity. It also facilitates the generation of output files that are
immediately ready for structure prediction applications.

Features
--------

AlphaFragment includes:

- **Domain Compilation**: Efficiently identifies protein domains using a
  combination of predictions from AlphaFoldDB, information from UniProt and
  manual input, ensuring accurate domain recognition.

- **Protein Fragmentation**: Provides tools to fragment proteins according to
  specified lengths, carefully avoiding disruption of identified domains and
  regions of interest, and adding overlap between fragments to maximise chances
  of identifying PPIs.

- **Visualization**: Includes functionality to generate detailed graphical
  representations of protein domains and their fragmentation patterns, for easy
  visual inspection and validation, as well as for use when analysing predictions.

- **Output Preparation**: Prepares and exports data in formats compatible with
  AlphaFold and AlphaPulldown.

Getting Started
---------------

This documentation is designed to guide you through the installation process
and getting started with using AlphaFragment. It also provides detailed
information on each class and function included the package.