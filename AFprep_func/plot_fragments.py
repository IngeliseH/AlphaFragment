"""
This module provides functionality to plot the AFprep fragmentation output of a protein.
It includes functions to plot domains and fragments of a protein on a matplotlib figure.
"""

import itertools
import os
import matplotlib.pyplot as plt
from matplotlib import patches

def plot_domain(ax, domain, base_y_position, domain_height, domain_colors):
    """
    Plots a rectangle representing a single protein domain, on a given
    matplotlib axis. The color of the domain is taken from a cycle of predefined
    colors.

    Parameters:
    - ax (matplotlib.axes.Axes): The matplotlib axis on which to plot the domain.
    - domain (Domain): Expected to be a domain type object
    - base_y_position (float): The base Y-axis position for the domain rectangle.
    - domain_height (float): The height of the domain rectangle.
    - domain_colors (iterator): An iterator that cycles through a list of colors.

    Note:
    - This function adds to the provided axis but does not show it. The display
      is managed by the caller.
    """
    color = next(domain_colors)
    rect = patches.Rectangle((domain.start, base_y_position - domain_height / 2),
                             domain.end - domain.start,
                             domain_height,
                             edgecolor='none',
                             facecolor=color,
                             alpha=0.5)
    ax.add_patch(rect)

def plot_fragment(ax, fragment, index, base_y_position, fragment_height, offset):
    """
    Plots a rectangle representing a single protein fragment, on a given
    matplotlib axis. Plots each fragment with a specified vertical offset to
    prevent overlap between consecutive fragments.

    Parameters:
    - ax (matplotlib.axes.Axes): The matplotlib axis on which to plot the
      fragment.
    - fragment (tuple): A tuple containing fragment start and end positions.
    - index (int): The index of the fragment, used to alternate vertical offset.
    - base_y_position (float): The base Y-axis position for fragment rectangles.
    - fragment_height (float): The height of the fragment rectangle.
    - offset (float): The vertical offset for each fragment (to distinguish
      between overlapping fragments).

    Note:
    - This function adds to the provided axis but does not show it. The display
      is managed by the caller.
    """
    start, end = fragment
    vertical_position = (base_y_position - (fragment_height / 2) +
                     (offset if index % 2 == 0 else -offset))
    fragment_rect = patches.Rectangle((start, vertical_position),
                                      end - start,
                                      fragment_height,
                                      edgecolor='black',
                                      facecolor='red',
                                      linewidth=1)
    ax.add_patch(fragment_rect)

def plot_fragmentation_output(protein, fragments, save_location=None, figsize=(12, 4)):
    """
    Creates and optionally saves a visualization of protein domains and fragments.

    This function generates a matplotlib figure depicting the domains and
    fragments of a protein.
    Domains are plotted as colored rectangles, and fragments are plotted as red
    rectangles with a vertical offset so they can be distinguished. 
    The resulting plot can be saved to a file by specifiying a save location.

    Parameters:
    - protein (Protein): The protein object, expected to have
      'domain_list' and 'last_res' attributes.
    - fragments (list of tuples): A list of tuples, each containing the start
      and end positions of a fragment.
    - save_location (str, optional): The directory path where the plot image
      will be saved. If not provided, the image is not saved.
    - figsize (tuple, optional): The size of the output figure.

    Returns:
    - matplotlib.figure.Figure: The figure object containing the plot.

    Note:
    - If `save_location` is provided and the directory does not exist, it will
      be created.
    - The function sets the figure's x-axis to represent the protein sequence
      position and removes the y-axis labels for clarity.
    """
    fig, ax = plt.subplots(figsize=figsize)
    base_y_position = 0.55
    domain_height = 0.4
    fragment_height = 0.05
    offset = 0.02  # Vertical offset for fragments
    domain_colors = itertools.cycle(['skyblue', 'pink', 'cyan', 'gold',
                                     'purple', 'silver', 'tan'])

    for domain in protein.domain_list:
        plot_domain(ax, domain, base_y_position, domain_height, domain_colors)

    for index, fragment in enumerate(fragments):
        plot_fragment(ax, fragment, index, base_y_position, fragment_height, offset)

    ax.set_xlim(0, protein.last_res)
    ax.set_ylim(0, 0.8)
    ax.set_xlabel('Protein Sequence Position')
    ax.set_yticks([])
    ax.set_title(f'Protein Domains and Fragments in {protein.name}')

    if save_location:
        if not os.path.exists(save_location):
            os.makedirs(save_location)
        fig.savefig(f"{save_location}/{protein.name}fragments.png", bbox_inches='tight')

    return fig
