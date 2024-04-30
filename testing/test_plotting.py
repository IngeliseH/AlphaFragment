"""
Test file for the plot_fragments module
"""
from unittest.mock import patch
import itertools
import pytest
from matplotlib import pyplot as plt
from matplotlib import colors
from functions.classes import Protein, Domain
from functions.plot_fragments import plot_domain, plot_fragment, draw_label, calculate_tick_freq, plot_fragmentation_output

@pytest.mark.parametrize("color_mode, domain_type, expected_color", [
    # color_mode 'type' with domain type 'AF'
    ('type', 'AF', 'orange'),
    # color_mode 'type' with domain type 'UniProt'
    ('type', 'UniProt', 'green'),
    # color_mode 'type' with domain type 'manually_defined'
    ('type', 'manually_defined', 'blue'),
    # color_mode 'type' with domain type 'not_in_dict'
    ('type', 'not_in_dict', 'gray'),
    # color_mode 'cycle'
    ('cycle', 'any', 'skyblue')
])
def test_plot_domain_color_modes(color_mode, domain_type, expected_color):
    """
    Test that the function correctly assigns colors to the domain rectangle based on the color mode and domain type
    """
    ax = plt.figure().add_subplot(111)
    domain = Domain("1", 1, 100, domain_type)
    base_y_position = 0.5
    domain_height = 0.1
    domain_color_cycle = itertools.cycle(['skyblue', 'pink'])

    plot_domain(ax, domain, base_y_position, domain_height, domain_color_cycle, color_mode=color_mode)
    rect = ax.patches[0]

    expected_rgba = tuple(list(colors.to_rgba(expected_color))[:-1] + [0.5])
    assert expected_rgba == rect.get_facecolor(), "Color does not match expected for color mode"

@pytest.mark.parametrize("color_mode, domain_type", [
    # invalid color mode
    ('invalid', 'AF'),
    # completely wrong type for color_mode
    (123, 'AF')
])
def test_plot_domain_invalid_color_mode(color_mode, domain_type):
    """
    Test that the function raises a ValueError when an invalid color mode is provided
    """
    ax = plt.figure().add_subplot(111)
    domain = Domain("1", 1, 100, domain_type)
    base_y_position = 0.5
    domain_height = 0.1
    domain_color_cycle = itertools.cycle(['skyblue', 'pink', 'cyan', 'gold', 'purple', 'silver', 'tan'])

    with pytest.raises(ValueError):
        plot_domain(ax, domain, base_y_position, domain_height, domain_color_cycle, color_mode=color_mode), "Function did not raise error for invalid color mode"

@pytest.mark.parametrize("index, expected_offset", [
    # even index, expect positive offset
    (0, 0.02),
    # odd index, expect negative offset
    (1, -0.02)
])
def test_plot_fragment_offsets(index, expected_offset):
    """
    Test that the fragment rectangle has the correct vertical position based on the index
    """
    ax = plt.figure().add_subplot(111)
    fragment = (50, 150)
    base_y_position = 0.5
    fragment_height = 0.05
    offset = 0.02

    plot_fragment(ax, fragment, index, base_y_position, fragment_height, offset)
    fragment_rect = ax.patches[0]

    expected_y_position = base_y_position - (fragment_height / 2) + expected_offset
    assert expected_y_position == fragment_rect.get_y(), f"Vertical position does not match expected for fragment index, expected {expected_y_position} but got {fragment_rect.get_y()}"

def test_plot_fragment_visual_properties():
    """
    Test that the fragment rectangle has the correct visual properties
    """
    ax = plt.figure().add_subplot(111)
    fragment = (200, 250)
    base_y_position = 0.5
    fragment_height = 0.05
    index = 0
    offset = 0.02

    plot_fragment(ax, fragment, index, base_y_position, fragment_height, offset)
    fragment_rect = ax.patches[0]

    assert fragment_rect.get_width() == 50, f"Fragment width incorrect, expected 50 but got {fragment_rect.get_width()}"
    assert fragment_rect.get_facecolor() == (1, 0, 0, 1), f"Fragment color incorrect, expected red but got {fragment_rect.get_facecolor()}"
    assert fragment_rect.get_linewidth() == 1, f"Fragment border width incorrect, expected 1 but got {fragment_rect.get_linewidth()}"

def test_draw_label():
    """
    Test that the labelling function draws a two line bracket on the plot
    """
    fig, ax = plt.subplots()
    # Coordinates for the bracket
    x_left = 100
    x_right = 300
    y = 0.5
    height = 0.1
    label_text = "TestLabel"
    draw_label(label_text, x_left, x_right, y, height, ax)

    # Check two lines are drawn for a bracket
    assert len(ax.lines) == 2, f"Expected 2 lines to be drawn for the curly bracket, but got {len(ax.lines)} lines"

    # Check line positions and properties
    for line in ax.lines:
        xdata = line.get_xydata()[0]
        assert xdata[0] == x_left or xdata[0] == x_right, f"Line x-coordinates do not match expected start or end of the bracket, got {xdata[0]}"
        assert line.get_color() == 'black', f"Line color does not match expected. Expected black, got {line.get_color()}"
        assert line.get_linewidth() == 1, f"Linewidth does not match expected. Expected linewidth 1, got {line.get_linewidth()}"

    # Check label placement
    assert len(ax.texts) == 1, f"Expected one label to be placed on the bracket, but got {len(ax.texts)} labels"
    text = ax.texts[0]
    assert text.get_text() == label_text, f"Label text does not match expected, expected {label_text} but got {text.get_text()}"
    assert text.get_position() == ((x_left + x_right) / 2, y - height), f"Label position is incorrect, expected {(x_left + x_right) / 2, y - height} but got {text.get_position()}"
    assert text.get_ha() == 'center' and text.get_va() == 'top', f"Text alignment is not as expected, expected center and top but got {text.get_ha()} and {text.get_va()}"
    assert text.get_rotation() == 45, f"Label rotation is not as expected, expected 45 degrees but got {text.get_rotation()}"
    assert text.get_fontsize() == 8, f"Font size of the label is not as expected, expected 8 but got {text.get_fontsize()}"


@pytest.mark.parametrize("num, expected", [
    # Basic tests with range of inputs
    (5, 1),
    (10, 1),
    (199, 10),
    (3000, 500),
    (7123456, 1000000),
    # input = 1
    (1, 1),
    # very large input
    (10000000, 1000000)
])
def test_tick_frequency(num, expected):
    """
    Test that the tick frequency function correctly calculates the tick frequency based on the number of residues
    """
    tick_freq = calculate_tick_freq(num)
    assert tick_freq == expected, f"Expected tick frequency of {expected} for input {num}, got {tick_freq}"

@pytest.mark.parametrize("num", [
    # String input
    "100",
    # Float input
    100.0,
    # Negative input
    -100,
    # Zero input
    0
])
def test_tick_frequency_invalid_input(num):
    """
    Test that the tick frequency function raises a ValueError when invalid input is provided
    """
    with pytest.raises(ValueError):
        calculate_tick_freq(num), "Invalid input should raise ValueError"

@patch("os.path.exists", return_value=False)
@patch("os.makedirs")
@patch("matplotlib.figure.Figure.savefig")
def test_plot_fragmentation_output_saving(mock_savefig, mock_makedirs, mock_exists):
    """
    Test that the function creates the save location directory and saves the figure when a save location is provided
    """
    fragments = [(50, 75), (180, 200)]
    save_location = "/fake/path"

    domain1 = Domain("1", 1, 100, 'AF')
    domain2 = Domain("2", 120, 200, 'UniProt')
    protein = Protein("TestProtein", "accession", "sequence", first_res=1, last_res=300, domain_list=[domain1, domain2], fragment_list=fragments)
    plot_fragmentation_output(protein, protein.fragment_list, save_location=save_location)
    mock_makedirs.assert_called_once_with(save_location), "Directory was not created when save location was provided"
    mock_savefig.assert_called_once(), "Figure was not saved when save location was provided"

def test_plot_fragmentation_output_no_save():
    """
    Test that the function does not save the figure when no save location is provided
    """
    fragments = [(50, 150)]
    domain1 = Domain("1", 1, 100, 'AF')
    domain2 = Domain("2", 120, 200, 'UniProt')
    protein = Protein("TestProtein", "accession", "sequence", first_res=1, last_res=300, domain_list=[domain1, domain2], fragment_list=fragments)

    with patch("matplotlib.figure.Figure.savefig") as mock_savefig:
        plot_fragmentation_output(protein, protein.fragment_list)
        mock_savefig.assert_not_called(), "Figure was saved when no save location was provided"

def test_plot_fragmentation_output_integration():
    """
    Test that the function plots correct number of domains and fragments on the figure
    """
    fragments = [(50, 75), (80, 150), (155, 175)]
    domain1 = Domain("1", 1, 100, 'AF')
    domain2 = Domain("2", 120, 200, 'UniProt')
    protein = Protein("TestProtein", "accession", "sequence", first_res=1, last_res=300, domain_list=[domain1, domain2], fragment_list=fragments)

    fig = plot_fragmentation_output(protein, protein.fragment_list)
    ax = fig.axes[0]

    assert len(ax.patches) == 5, f"Should have 5 elements (3 fragments and 2 domains) plotted, but got {len(ax.patches)} elements"

def test_large_protein_many_domains():
    """
    Test that the function can handle a very long protein with a large number of domains and fragments
    """
    # Generating a large number of domains
    domain_list = [Domain(identifier=str(i), start=i * 100, end=i * 100 + 49, domain_type='AF') for i in range(100)]
    fragments = [(i * 100 + 25, i * 100 + 75) for i in range(100)]
    protein = Protein("VeryLargeProtein", "accession", "sequence", first_res=1, last_res=10000, domain_list=domain_list, fragment_list=fragments)
    fig = plot_fragmentation_output(protein, protein.fragment_list)
    assert len(fig.axes[0].patches) == 200, f"Expected 200 patches on the plot (100 domains + 100 fragments), but got {len(fig.axes[0].patches)} patches"

def test_overlapping_domains():
    """
    Test that the function can handle overlapping domains
    """
    fragments = [(200, 350)]
    domains = [Domain("1", 100, 300, "AF"), Domain("2", 250, 450, "UniProt")]
    protein = Protein('OverlappingProtein', 'accession', 'sequence', 1, 500, domain_list=domains, fragment_list=fragments)

    fig = plot_fragmentation_output(protein, fragments)
    assert len(fig.axes[0].patches) == 3, "Expected 3 patches on the plot (2 domains + 1 fragment)"


def test_visual_output():
    """
    Test the visual output of the plot_fragmentation_output function
    """
    # Set up large protein
    domains = [Domain("AF1", 10, 30, "AF"), Domain("UP2", 50, 90, "UniProt"), Domain("M3", 250, 500, "manually_defined"), Domain("O4", 700, 900, "other"), Domain("O5", 1300, 1900, "other")]
    fragments = [(0, 110), (100, 250), (240, 550), (530, 700), (690, 1000), (990, 1250), (1240, 1952)]
    protein = Protein("LargeProtein", "accession", "sequence", first_res=1, last_res=1952, domain_list=domains, fragment_list=fragments)

    #close any existing plots
    plt.close('all')

    plot_fragmentation_output(protein, fragments)
    plot_fragmentation_output(protein, fragments, label=['AF', 'UniProt'])
    plot_fragmentation_output(protein, fragments, color_mode='cycle')

    # Set up small protein for checking boundaries
    domains = [Domain("AF1", 1, 5, "AF"), Domain("UP2", 10, 15, "UniProt")]
    fragments = [(0, 7), (6, 10), (10, 17)]
    protein = Protein("SmallProtein", "accession", "sequence", first_res=0, last_res=16, domain_list=domains, fragment_list=fragments)

    plot_fragmentation_output(protein, fragments, label=['AF', 'UniProt', 'manually_defined'])

    plt.show()  # This will display the plot window for manual verification
