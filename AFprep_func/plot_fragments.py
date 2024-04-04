import matplotlib.pyplot as plt
import matplotlib.patches as patches
import itertools
import os

def plot_fragments(protein, fragments, save_location = None):
    fig, ax = plt.subplots(figsize=(12, 4))
    base_y_position = 0.55
    domain_height = 0.4
    fragment_height = 0.05
    domain_alpha = 0.5  # Opacity for 'DOMAIN' and 'REGION' types
    domain_colors = itertools.cycle(['skyblue', 'pink', 'cyan', 'gold', 'purple', 'silver', 'tan'])

    # Plot domains         
    for domain in protein.domain_list:
        start_pos = domain.start
        end_pos = domain.end
        color = next(domain_colors)
        height = domain_height
        alpha = domain_alpha
        rect = patches.Rectangle((start_pos, base_y_position - height / 2), end_pos - start_pos, height, edgecolor='none', facecolor=color, alpha=alpha)
        ax.add_patch(rect)

     # Calculate and plot fragment positions with offset
    sequence_length = len(protein.sequence)
    fragment_positions = []
    for fragment in range(len(fragments)):
        start = fragments[fragment][0]
        end = fragments[fragment][1]
        fragment_positions.append([start, end])

    offset = 0.02  # Vertical offset for fragments
    for index, (start, end) in enumerate(fragment_positions):
        vertical_position = base_y_position - fragment_height / 2 + (offset if index % 2 == 0 else -offset)
        fragment_rect = patches.Rectangle((start, vertical_position), end - start, fragment_height, edgecolor='black', facecolor='red', linewidth=1)
        ax.add_patch(fragment_rect)

    # Formatting the plot
    ax.set_xlim(0, sequence_length)
    ax.set_ylim(0, 0.8)
    ax.set_xlabel('Protein Sequence Position')
    ax.set_yticks([])
    ax.set_title(f'Protein Domains and Fragments in {protein.name}')

    if save_location:
        # Check if the directory exists, and create it if it doesn't
        if not os.path.exists(save_location):
            os.makedirs(save_location)
    
    # Save the figure
    fig.savefig(f"{save_location}/{protein.name}fragments.png")


    return fig