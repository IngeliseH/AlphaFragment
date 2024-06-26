"""
This module provides functionality to plot the AFprep fragmentation output of a protein.

Functions:
  - plot_domain: Plots a rectangle representing a single protein domain
  - plot_fragment: Plots a rectangle representing a single protein fragment
  - draw_label: Adds a label to the plot
  - calculate_tick_freq: Calculates the frequency of x axis ticks based on sequence length
  - plot_fragmentation_output: Creates and optionally saves a visualization of
    protein domains and fragments

Dependencies:
  - itertools: Used for cycling through colors for the protein domains.
  - os: Used for creating directories if they do not exist.
  - matplotlib: Used for plotting the protein domains and fragments.
"""

import itertools
import os
import matplotlib.pyplot as plt
from matplotlib import patches

def plot_domain(ax, domain, base_y_position, domain_height, domain_color_cycle, color_mode='type'):
    """
    Plots a rectangle representing a single protein domain, on a given
    matplotlib axis.  The color of the domain can be either cycled through a set
    of predefined colors or set based on domain type.

    Parameters:
      - ax (matplotlib.axes.Axes): The matplotlib axis on which to plot the domain.
      - domain (Domain): Expected to be a domain type object
      - base_y_position (float): The base Y-axis position for the domain rectangle.
      - domain_height (float): The height of the domain rectangle.
      - domain_color_cycle (itertools.cycle): An itertools.cycle object that
        cycles through a set of colors for the protein domains.
      - color_mode (str): 'cycle' for cycling colors or 'type' for color based on domain type.
      

    Note:
      - This function adds to the provided axis but does not show it. The display
        is managed by the caller.
    """
    if color_mode not in ['type', 'cycle']:
        raise ValueError("color_mode must be 'type' or 'cycle'")
    type_colors = {'AF': 'orange', 'UniProt': 'green', 'manually_defined': 'blue'}
    default_color = 'gray'
    if color_mode == 'type':
        if domain.type in type_colors:
            color = type_colors[domain.type]
        else:
            print(f"Domain type {domain.type} not in type_colours. Used default color {default_color} instead.")
            color = default_color
    elif color_mode == 'cycle':
        color = next(domain_color_cycle)

    rect = patches.Rectangle((domain.start+0.5, base_y_position - domain_height / 2),
                             domain.end-domain.start+1,
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
      - Adds 0.5 to either end of the fragment to center it on the residue position.
    """
    start, end = fragment
    vertical_position = (base_y_position - (fragment_height / 2) +
                     (offset if index % 2 == 0 else -offset))
    fragment_rect = patches.Rectangle((start + 0.5, vertical_position),
                                      end - start,
                                      fragment_height,
                                      edgecolor='black',
                                      facecolor='red',
                                      linewidth=1)
    ax.add_patch(fragment_rect)

def draw_label(label, x_left, x_right, y, bracket_height, ax):
    """
    Adds a bracket to the plot labelling a specified region with the given text

    Parameters:
      - label (str): The text to be displayed on the bracket.
      - x_left (float): The left x-coordinate of the bracket.
      - x_right (float): The right x-coordinate of the bracket.
      - y (float): The y-coordinate of the bracket.
      - bracket_height (float): The height of the bracket.
      - ax (matplotlib.axes.Axes): The matplotlib axis on which to plot the bracket.
  
    Note:
      - This function adds to the provided axis but does not show it. The display
        is managed by the caller.
    """
    mid_x = (x_left + x_right) / 2
    ax.plot([x_left, mid_x], [y, y - bracket_height], color='black', lw=1)
    ax.plot([x_right, mid_x], [y, y - bracket_height], color='black', lw=1)

    text_y_position = y - bracket_height
    ax.text((x_left + x_right) / 2, text_y_position, label, ha='center',
            va='top', fontsize=8, rotation=45)

def calculate_tick_freq(n):
    """
    Calculates the frequency of x axis ticks for a given number of residues.

    Parameters:
      - n (int): The number of residues in the protein sequence.

    Returns:
      - int: The frequency of x-axis ticks.
    """
    # Check if the number is a positive integer
    if not isinstance(n, int) or n <= 0:
        raise ValueError("Input must be a positive integer")

    # Convert integer to string to easily access and count digits
    n_str = str(n)
    digit_count = len(n_str)  # Count the number of digits
    leading_digit = int(n_str[0])  # Get the leading digit

    # Special case for single digit numbers
    if digit_count == 1:
        return 1

    # Apply conditions based on the leading digit and number length
    if leading_digit == 1:
        return 10 ** (digit_count - 2)
    if leading_digit in {5, 6, 7, 8, 9}:
        return 10 ** (digit_count - 1)
    if leading_digit in {2, 3, 4}:
        return 5*10 ** (digit_count - 2)
    raise ValueError("Unexpected leading digit")

def plot_fragmentation_output(protein, fragments, save_location=None,
                              figsize=(12, 4), color_mode='type', label=None):
    """
    Creates and optionally saves a visualization of protein domains and fragments.
    Domains are plotted as colored rectangles, and fragments as red
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
      - color_mode ('type' or 'cycle): Whether to color domains by type (AlphaFold,
        UniProt, manually defined), or using a series of colours to distinguish
        nearby domains
      - label (list, optional): List of sources for which domains should be labelled.
        Can include 'UniProt', 'AF' and 'manually_defined'. Defaults to None.


    Returns:
      - matplotlib.figure.Figure: The figure object containing the plot.

    Note:
      - If `save_location` is provided and the directory does not exist, it will
        be created.
      - x-axis represents the protein sequence position, with 1-based indexing.
    """

    fig, ax = plt.subplots(figsize=figsize)
    base_y_position = 0.55
    domain_height = 0.4
    fragment_height = 0.05
    offset = 0.025  # Vertical offset for fragments
    domain_color_cycle = itertools.cycle(['skyblue', 'pink', 'cyan', 'gold',
                                     'purple', 'silver', 'tan'])

    for domain in protein.domain_list:
        plot_domain(ax, domain, base_y_position, domain_height, domain_color_cycle, color_mode)

        if label and domain.type in label:
            # Draw label
            bracket_height = 0.05
            draw_label(domain.id, domain.start+0.5, domain.end+1.5,
                       base_y_position - domain_height / 2,
                       bracket_height, ax)

    for index, fragment in enumerate(fragments):
        plot_fragment(ax, fragment, index, base_y_position, fragment_height, offset)

    ax.set_xlim(0, protein.last_res + 2)
    ax.set_ylim(0, 0.8)
    ax.set_xlabel('Protein Sequence Position')
    # Setting the X-axis
    ax.set_xlim(1, protein.last_res + 1)  # +1 to shift end limit to 1-based index

    # Setting ticks
    tick_freq = calculate_tick_freq(protein.last_res)
    ticks = list(range(0, protein.last_res + 1, tick_freq))  # Shift by 1 for 1-based index
    ticks.remove(0)
    ticks.append(1)
    ticks.append(protein.last_res + 1)
    ax.set_xticks(ticks)
    ax.set_xticklabels([str(tick) for tick in ticks])

    ax.set_yticks([])
    ax.set_title(f'Protein Domains and Fragments in {protein.name}')

    if save_location:
        if not os.path.exists(save_location):
            os.makedirs(save_location)
        fig.savefig(f"{save_location}/{protein.name}fragments.png", bbox_inches='tight')

    return fig
