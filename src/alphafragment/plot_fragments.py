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
import textwrap
import re

def plot_domain(ax, domain, base_y_position, domain_height, domain_color_cycle, color_mode='origin', type_colors={}):
    """
    Plots a rectangle representing a single protein domain, on a given
    matplotlib axis. The color of the domain can be either cycled through a set
    of predefined colors or set based on domain origin or type.

    Parameters:
      - ax (matplotlib.axes.Axes): The matplotlib axis on which to plot the domain.
      - domain (Domain): Expected to be a domain type object.
      - base_y_position (float): The base Y-axis position for the domain rectangle.
      - domain_height (float): The height of the domain rectangle.
      - domain_color_cycle (itertools.cycle): An itertools.cycle object that
        cycles through a set of colors for the protein domains.
      - color_mode (str): 'cycle' for cycling colors, 'origin' for color based on domain
        origin, or 'type' for color based on domain identifier.
      - type_colors (dict): A dictionary to store colors for 'type' color mode.

    Returns:
      - dict: The updated type_colors dictionary, only used if color_mode is 'type'.

    Note:
      - Although this function does have a return, its primary purpose is to add
        a rectangle to the provided axis. The return is used to maintain color
        consistency when plotting multiple domains with the 'type' color mode.
    """

    if color_mode not in ['origin', 'cycle', 'type']:
        raise ValueError("color_mode must be 'origin', 'cycle', or 'type'")
    
    origin_colors = {'AF': 'skyblue', 'UniProt': 'pink', 'manually_defined': '#9FD18A'}
    default_color = 'gray'
    
    if color_mode == 'origin':
        color = origin_colors.get(domain.type, default_color)
    elif color_mode == 'type':
        # Extract the text part of the identifier, ignoring any numbers
        identifier_text = re.sub(r'\d+', '', domain.id)
        if identifier_text not in type_colors:
            type_colors[identifier_text] = next(domain_color_cycle)
        color = type_colors[identifier_text]
    elif color_mode == 'cycle':
        color = next(domain_color_cycle)

    rect = patches.Rectangle((domain.start + 0.5, base_y_position - domain_height / 2),
                             domain.end - domain.start + 1,
                             domain_height,
                             edgecolor='none',
                             facecolor=color,
                             alpha=0.5)
    ax.add_patch(rect)
    return type_colors

def plot_fragment(ax, fragment, index, base_y_position, fragment_height, offset):
    """
    Plots a rectangle representing a single protein fragment, on a given
    matplotlib axis. Plots each fragment with a specified vertical offset to
    prevent overlap between consecutive fragments.

    Parameters:
      - ax (matplotlib.axes.Axes): The matplotlib axis on which to plot the fragment.
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

def draw_label(ax, label, x, x_left, x_right, y, bracket_height, max_x):
    """
    Adds a bracket to the plot labeling a specified region with the given text.

    Parameters:
      - label (str): The text to be displayed on the bracket.
      - x (float): The x-coordinate of the label.
      - x_left (float): The left x-coordinate of the bracket.
      - x_right (float): The right x-coordinate of the bracket.
      - y (float): The y-coordinate of the bracket.
      - bracket_height (float): The height of the bracket.
      - max_x (float): The maximum x-coordinate allowed for the label to prevent it from spilling out of the plot axes.
      - ax (matplotlib.axes.Axes): The matplotlib axis on which to plot the bracket.
  
    Note:
      - This function adds to the provided axis but does not show it. The display
        is managed by the caller.
    """
    ax.plot([x_left, x], [y, y - bracket_height], color='black', lw=1)
    ax.plot([x_right, x], [y, y - bracket_height], color='black', lw=1)

    text_y_position = y - bracket_height
    wrapped_label = "\n".join(textwrap.wrap(label, width=20))
    x = min(x, max_x)  # Ensure the label does not spill out of the plot axes
    ax.text(x, text_y_position, wrapped_label, ha='center',
            va='top', fontsize=8, rotation=90, linespacing=0.8)

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
        return 5 * 10 ** (digit_count - 2)
    raise ValueError("Unexpected leading digit")

def plot_fragmentation_output(protein, fragments, save_location=None,
                              figsize=(12, 4), color_mode='type', label=['UniProt', 'manually_defined']):
    """
    Creates and optionally saves a visualization of protein domains and fragments.
    Domains are plotted as colored rectangles, and fragments as red rectangles 
    with a vertical offset so they can be distinguished. The resulting plot can 
    be saved to a file by specifying a save location.

    Parameters:
      - protein (Protein): The protein object, expected to have 'domain_list' and
        'last_res' attributes.
      - fragments (list of tuples): A list of tuples, each containing the start and
        end positions of a fragment.
      - save_location (str, optional): The directory path where the plot image will
        be saved. If not provided, the image is not saved.
      - figsize (tuple, optional): The size of the output figure.
      - color_mode ('origin', 'cycle', or 'type'): Whether to color domains by origin,
        cycle colors, or type. Defaults to 'type'.
      - label (list, optional): List of sources for which domains should be labeled.
        Can include 'UniProt', 'AF', and 'manually_defined'. Defaults to labelling
        UniProt and manually defined domains.

    Returns:
      - matplotlib.figure.Figure: The figure object containing the plot.

    Note:
      - If `save_location` is provided and the directory does not exist, it will be
        created.
      - x-axis represents the protein sequence position, with 1-based indexing.
    """

    fig, ax = plt.subplots(figsize=figsize)
    base_y_position = 0.55
    fragment_height = 0.05
    offset = 0.025
    domain_color_cycle = itertools.cycle(['#57BBEA', '#ED7AB0', '#8DC640', '#F68B1F',
                                          '#9F83BC', '#FFCC71', '#6DC8BF', '#A74399',
                                          '#A6ADD3', '#E64425', '#C2C1B1', '#00A45D',
                                          '#BA836E', '#3E4291'])

    type_colors = {}  # Dictionary to store colors for 'type' color mode

    for domain in protein.domain_list:
        if domain.type == 'AF':
            domain_height =0.4
        else:
            domain_height = 0.36
        # retain the updated type_colors dictionary to maintain color consistency if 'type' color mode is used, and for use in legend
        type_colors = plot_domain(ax, domain, base_y_position, domain_height, domain_color_cycle, color_mode, type_colors)

    y = base_y_position - domain_height / 2
    tentative_label_positions = []
    max_x = protein.last_res  # Maximum x-coordinate for the label to prevent spillover
    intertext_distance = protein.last_res / 70 # Minimum distance between labels

    if label:
        # Check label is a list, if not ignore and print error
        if not isinstance(label, list):
            print("Label must be a list of sources for which domains should be labeled." +
                  "Source options are 'UniProt', 'AF', and 'manually_defined'. Switching to default.")
            label = ['UniProt', 'manually_defined']
        for domain in protein.domain_list:
            if domain.type in label:
                bracket_height = 0.05
                x_left, x_right = domain.start + 0.5, domain.end + 1.5
                
                # Adjust horizontal position to avoid overlap
                x = (x_left + x_right) / 2
                # Check if any value in label_positions is within set distance of x
                while any(abs(x - pos) < intertext_distance for pos in tentative_label_positions):
                    x += intertext_distance / 4
                tentative_label_positions.append(x)

        # Adjust label positions to prevent overlap at the right-hand side
        adjusted_positions = [max_x]
        for pos in reversed(tentative_label_positions):
            while any(abs(x - pos) < intertext_distance for x in adjusted_positions):
                    pos -= intertext_distance / 4
            pos = min(pos, max_x-intertext_distance)  # Ensure the label does not spill out of the plot axes
            adjusted_positions.append(pos)
        adjusted_positions = list(reversed(adjusted_positions))
        i = 0
        for _, domain in enumerate(protein.domain_list):
            if domain.type in label:
                bracket_height = 0.05
                x_left, x_right = domain.start + 0.5, domain.end + 1.5
                x = adjusted_positions[i]
                i += 1
                draw_label(ax, domain.id, x, x_left, x_right, y, bracket_height, max_x)

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
    ticks.sort()

    if ticks[-1] - ticks[-2] < tick_freq / 3:
        ticks.pop(-2)

    ax.set_xticks(ticks)
    ax.set_xticklabels([str(tick) for tick in ticks])

    ax.set_yticks([])
    ax.set_title(f'Protein Domains and Fragments in {protein.name}')

    if color_mode in ['origin', 'type']:
        if color_mode == 'origin':
            colors = {'AF': 'skyblue', 'UniProt': 'pink', 'manually_defined': '#9FD18A'}
            legend_labels = {'AF': 'AlphaFold', 'UniProt': 'UniProt', 'manually_defined': 'Manually Defined'}
        elif color_mode == 'type':
            colors = type_colors
            legend_labels = {key: 'AlphaFold predicted' if key == 'AF_D' else key for key in colors.keys()}

        domain_handles = [patches.Patch(color=color, label=legend_labels[key]) for key, color in colors.items()]
        fragment_handle = patches.Patch(facecolor='red', edgecolor='black', linewidth=1, label='Fragment')

        # Combine legends
        handles = [fragment_handle] + domain_handles
        labels = [handle.get_label() for handle in handles]

        ax.legend(handles, labels, loc='upper left', bbox_to_anchor=(1, 1), title='Domains/Fragments')

    if save_location:
        if not os.path.exists(save_location):
            os.makedirs(save_location)
        fig.savefig(f"{save_location}/{protein.name}fragments.png", bbox_inches='tight')

    return fig
