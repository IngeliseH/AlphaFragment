"""
Test file for the fragmentation_methods module.
"""
import pytest
from unittest import mock
from alphafragment.classes import Domain, Protein
from alphafragment.fragmentation_methods import validate_fragmentation_parameters, merge_overlapping_domains, check_valid_cutpoint, recursive_fragmentation, break_in_half

@pytest.mark.parametrize(
    "domains, expected",
    [
        # No Overlap
        ([Domain(1, 1, 5, 'TYPE'), Domain(2, 7, 10, 'TYPE')],
         [Domain(1, 1, 5, 'TYPE'), Domain(2, 7, 10, 'TYPE')]),
        # Partial Overlap
        ([Domain(1, 1, 5, 'TYPE'), Domain(2, 4, 10, 'TYPE')],
         [Domain(1, 1, 10, 'TYPE')]),
        # Complete Overlap
        ([Domain(1, 1, 10, 'TYPE'), Domain(2, 2, 9, 'TYPE')],
         [Domain(1, 1, 10, 'TYPE')]),
        # Multiple Merging
        ([Domain(1, 1, 5, 'TYPE'), Domain(2, 4, 8, 'TYPE'), Domain(3, 7, 10, 'TYPE')],
         [Domain(1, 1, 10, 'TYPE')]),
        # Touching but not overlapping
        ([Domain(1, 1, 5, 'TYPE'), Domain(2, 6, 10, 'TYPE')],
         [Domain(1, 1, 5, 'TYPE'), Domain(2, 6, 10, 'TYPE')]),
        # Order Independence
        ([Domain(2, 4, 8, 'TYPE'), Domain(1, 1, 5, 'TYPE'), Domain(3, 7, 10, 'TYPE')],
         [Domain(1, 1, 10, 'TYPE')]),
    ]
)
def test_merge_overlapping_domains(domains, expected):
    """
    Tests the merge_overlapping_domains function with various arrangements of domain
    overlaps to ensure domains are merged correctly or left separate as appropriate.
    """
    result = merge_overlapping_domains(domains)
    assert len(result) == len(expected), "Number of resulting domains does not match expected"
    for res_domain, exp_domain in zip(result, expected):
        assert res_domain.id == exp_domain.id, f"Domain id mismatch: expected {exp_domain.id}, got {res_domain.id}"
        assert res_domain.start == exp_domain.start, f"Start mismatch: expected {exp_domain.start}, got {res_domain.start}"
        assert res_domain.end == exp_domain.end, f"End mismatch: expected {exp_domain.end}, got {res_domain.end}"
        assert res_domain.type == exp_domain.type, f"Type mismatch: expected {exp_domain.type}, got {res_domain.type}"

@pytest.mark.parametrize(
    "res, domains, sequence_end, expected",
    [
        # Within Domain - not allowed
        (7, [Domain(1, 5, 10, 'TYPE')], 15, False),
        # At domain end, allowed
        (11, [Domain(1, 5, 10, 'TYPE')], 15, True),
        # At domain start, allowed
        (5, [Domain(1, 5, 10, 'TYPE')], 15, True),
        # Between Domains
        (11, [Domain(1, 1, 10, 'TYPE'), Domain(2, 11, 20, 'TYPE')], 20, True),
        # Cut before sequence start
        (-1, [], 10, False),
        # Cut at sequence start
        (0, [], 10, True),
        # Cut beyond sequence end
        (12, [], 10, False),
        # Cut at sequence end
        (11, [], 10, True),
        # End of one domain is start of another (= 1 res overlap, so not allowed)
        (10, [Domain(1, 5, 10, 'TYPE'), Domain(2, 10, 15, 'TYPE')], 15, False),
    ]
)
def test_check_valid_cutpoint(res, domains, sequence_end, expected):
    """
    Tests the check_valid_cutpoint function with various residue positions and
    domain configurations, verifying that the function correctly determines if a
    cut is valid.
    """
    result = check_valid_cutpoint(res, domains, sequence_end)
    assert result == expected, f"Expected {expected} but got {result} for residue {res} with sequence end {sequence_end}"

@pytest.mark.parametrize(
    "domains, length, overlap",
    [
        # Basic Functionality: No domains, simple split
        ([], {'min': 10, 'ideal': 15, 'max': 20}, {'min': 1, 'ideal': 2, 'max': 3}),
        # Domain Handling: No split through domain
        ([Domain(1, 5, 10, 'TYPE')], {'min': 5, 'ideal': 7, 'max': 10}, {'min': 1, 'ideal': 2, 'max': 3}),
        # Overlap Rules: Ideal overlap not possible
        ([], {'min': 5, 'ideal': 7, 'max': 10}, {'min': 1, 'ideal': 3, 'max': 4}),
    ]
)
def test_recursive_fragmentation_lengths_and_overlaps(domains, length, overlap):
    """
    Tests the recursive_fragmentation function to ensure that all fragments meet
    the length and overlap requirements specified by the parameters.
    """
    result = recursive_fragmentation(Protein("Protein1", "example_acc_id", 'A'*20),
                                     domains, 0, length, overlap, original_max_len=length['max'])

    # Check if any fragments were generated
    for i in range(len(result)):
        # Check fragment length
        start, end = result[i]
        fragment_length = end - start
        assert length['min'] <= fragment_length <= length['max'], f"Fragment length {fragment_length} out of bounds ({length['min']}, {length['max']}) at index {i}"

        # Check overlap with next fragment if not the last one
        if i < len(result) - 1:
            next_start = result[i + 1][0]
            actual_overlap = start + fragment_length - next_start
            assert overlap['min'] <= actual_overlap <= overlap['max'], f"Overlap {actual_overlap} out of bounds ({overlap['min']}, {overlap['max']}) between fragments {i} and {i+1}"

@pytest.mark.parametrize("length, overlap, expected_error",
    [
        # max_len < min_len
        ({'min': 60, 'ideal': 55, 'max': 50}, {'min': 1, 'ideal': 2, 'max': 5}, ValueError),
        # max overlap < min_overlap
        (None, {'min': 2, 'ideal': 2, 'max': 0}, ValueError),
        # max overlap = min_len
        ({'min': 10, 'ideal': 15, 'max': 20}, {'min': 1, 'ideal': 2, 'max': 10}, ValueError),
        # ideal overlap not between min and max overlap
        (None, {'min': 1, 'ideal': 7, 'max': 5}, ValueError),
        # length values not integers
        ({'min': 10.5, 'ideal': 15, 'max': 20}, None, TypeError),
        # overlap values not integers
        (None, {'min': 1.5, 'ideal': 2, 'max': 5}, TypeError),
        # length values string
        ({'min': '10', 'ideal': '15', 'max': '20'}, None, TypeError),
        # overlap values string
        (None, {'min': '1', 'ideal': '2', 'max': '5'}, TypeError),
        # missing keys in length
        ({'min': 10, 'max': 20}, None, ValueError),
        # missing keys in overlap
        (None, {'min': 1, 'ideal': 2}, ValueError),
        # length not a dictionary
        (20, None, TypeError),
        # overlap not a dictionary
        (None, 5, TypeError),
        # min_len < 0
        ({'min': -10, 'ideal': 15, 'max': 20}, None, ValueError),
        # min_len == 0
        ({'min': 0, 'ideal': 15, 'max': 20}, None, ValueError),
        # min_overlap < 0
        (None, {'min': -1, 'ideal': 2, 'max': 5}, ValueError)
    ]
)
def test_validation_error_handling(length, overlap, expected_error):
    """
    Tests that the validate_fragmentation_parameters function raises correct
    errors for invalid length and overlap parameters.
    """
    if length is None:
        length = {'min': 10, 'ideal': 15, 'max': 20}
    if overlap is None:
        overlap = {'min': 1, 'ideal': 2, 'max': 5}
    
    try:
        validate_fragmentation_parameters(Protein("Protein1", "example_acc_id", 'A'*20),
                                          length, overlap)
    except expected_error:
        pass
    else:
        pytest.fail(f"Expected {expected_error} but no error was raised for input {length}, {overlap}")

@pytest.mark.parametrize("length, overlap",
    [
        # ideal_len == min_len
        ({'min': 10, 'ideal': 10, 'max': 20}, {'min': 1, 'ideal': 2, 'max': 5}),
        # ideal_len == max_len
        ({'min': 10, 'ideal': 20, 'max': 20}, {'min': 1, 'ideal': 2, 'max': 5}),
        # ideal_overlap == min_overlap
        ({'min': 10, 'ideal': 15, 'max': 20}, {'min': 1, 'ideal': 1, 'max': 5}),
        # ideal_overlap == max_overlap
        ({'min': 10, 'ideal': 15, 'max': 20}, {'min': 1, 'ideal': 5, 'max': 5}),
    ]
)
def test_valid_parameters_no_error(length, overlap):
    """
    Tests that the validate_fragmentation_parameters function does not raise
    an error for valid length and overlap parameters where the ideal values
    are equal to the min or max values.
    """
    try:
        validate_fragmentation_parameters(Protein("Protein1", "example_acc_id", 'A'*20),
                                          length, overlap)
    except Exception as e:
        pytest.fail(f"Unexpected error {e} raised for input {length}, {overlap}")

def test_break_in_half_does_not_split_through_central_domain():
    """
    Tests that break_in_half does not split a protein through a central domain.
    Specifically checks that the two halves created are larger than the minimum length
    and no domain is split by the function.
    """
    # Create a protein with a central domain
    sequence = 'A' * 2000
    domains = [Domain(identifier="1", start=950, end=1050, domain_type="TYPE")]  # Central domain (100 residues)

    protein = Protein("Protein1", "example_acc_id", sequence)

    # Fragmentation parameters
    length = {'min': 100, 'ideal': 150, 'max': 200}
    overlap = {'min': 10, 'ideal': 20, 'max': 30}

    # Break the protein in half
    first_half, second_half = break_in_half(protein, length, overlap)

    # Ensure the halves were created and their lengths are above the minimum
    assert first_half is not None and second_half is not None, "Failed to split the protein into halves"
    assert (first_half.last_res - first_half.first_res + 1) >= length['min'], "First half is smaller than the minimum fragment length"
    assert (second_half.last_res - second_half.first_res + 1) >= length['min'], "Second half is smaller than the minimum fragment length"

    # Ensure the central domain is not split
    assert not (first_half.first_res < domains[0].start < first_half.last_res and
                first_half.first_res < domains[0].end < first_half.last_res), "Central domain was split by the first half"
    assert not (second_half.first_res < domains[0].start < second_half.last_res and
                second_half.first_res < domains[0].end < second_half.last_res), "Central domain was split by the second half"


def test_recursive_fragmentation_ana1():
    """
    Tests recursive fragmentation on the Ana1 protein (1729 residues long) with multiple domains
    including AF and UniProt domains. Ensures that the break_in_half function does not split through
    a domain and that both halves are above the minimum length.
    """
    # Create the Ana1 protein sequence
    ana1_sequence = (
        "MALQLTVNGKFGANLRRTIGPLDDPEQVRKDLESERRITRLKQVREKSNHLARQIREDVAAEKRRQIHKLEQVKQRELNAWREHVLVKKHQDYRSAIFQVGAAHRAAKEENERIE"
        "QQRQERIDKIRRCRKQALKRSAKASVELRTTNTMNLNAEGRASAGTQTPIVEDKENRMESGCNKPCCKGKRKRKTSCGCASTDAEDEKEDDASSDDDSSILIIESPSSNRRLQKT"
        "TPVILDVEIEETISETDKNPEGMDINDRFMQTNRKFSHVVRPSDDEPQRRPRFTQISDLVRKTETTVRAGPSQRREADPTPPVSAPPSPTKSPPRSPRKCVEREEPVPAPAPRK"
        "SPTKTAKGSPKKTVTAPRNQPSVKRTGLKLNPAPAKVIDAGIHRGQKAKAPATLPKVTAIQKPAVEPPRNLPPIPEQAMPAVPNPSMTHCYPMPQNQPYMHPYAQLPMQPYAMPY"
        "SMPFPVQAQFPQPHPMYAPQQMMPNQPAVLGPAITAPPSTATQSTVSTTTITMSSRQDPRTTQCGRVQFYDHSNKYHRTYEAPTQSVQCNEKDATQLTAMDHARIENQLRDLRE"
        "QELDKLRKISDDRGQKALEREQVRRDCAELTEKLDALTQQQPQLLPTDANIIPSHRYADAAARREQKMNDAMEEMLLRPAIITCPEVRTKPTPPNLSRSKSKASVAVNLGEPPV"
        "QIGRDNLGSSESCSSILLDYVNDQSKQLKSDLQAEQSNSLKSMRLKSLLERIEKIRVQLLEELKAGESGASKGDNAQELNNIRQERADILSERTRTLNERESDLQQKEAILEQR"
        "LRKFYKETKNAKSSESDKPKSNAKEDGPVEIIIKVRSDGTVKQYVPRTKSKAKTKPSAIDNEKDVPLGSTDSTPRDDEPKPTEQDHGKRPFLQDQRQISIDSNSTAYRSLPPVS"
        "YKNMTTAAPPSAPLHPMVVQYINRLLGMTRQNIDEMGVSSSCVTTPSTSVINSLRNVTPCPEPEAGRAERDQDENTEMNEHRMERVQTFIADNRSLINDLEESIRCQQQLQREQ"
        "QSLDAEKSMRAFDQIWNKRLAKNPEDCQQTNKEKEQTREKMPRERPTPQGDRKNAAPPQESGKIRRTQEGDRQQQTNGAVQQKREQKTVKQGSGADKSTNRFAPEKSKALSTSE"
        "ESSARNMERYAQLTENCTQRIAELTELITKVREEKQRLVEVTLTSNSDGERQSTEYLELPTGQQQTRCRTISDRSDSQTPSTSEALPLQKHKPTAASRDSGIADSRPITAQGQV"
        "CIDVEPISLGSSTQNNTARGRVRAPPATIRRYSPQLNAEDLAHELSTITEVETPGQSHIVAATPVPKPFPSFDQYAKELHLDLSRLDADQSQRLQGEFNDLILAIRQRNQGTDY"
        "REFPSINAYLHNMTSTRIHVEIGQDADQATLSPGELMRQLRVITHSIQDFPKRREYFEKLMAKQPPDQRDLIDSASLENDSTDSFNVEEELRQRKILKNSFRRGNATETTLTAQ"
        "EVASSTRRESVAPVTNDQPNESGIDPLSGSNFSSDADQRSPCWHATMHERQQQVDELCSSTTASSPERRPRKSRLRGGHEHNNRSKDDSTVQDASQIGRSLNLREFLTRELLKH"
        "RVHDDVASESSDESLRGHFLKSVLHSLSPSNNTQTPGVGVSHATGATNDRQKTSTPVGSFLSIPEKAGSMSNTGSQLFSGESRISLVNYPDGTPPIPFEQQSKLSTNQTRQSSG"
        "QIGSGVRTTHRKSPRK"
    )

    # Define the Ana1 domains
    domains = [
        Domain(identifier="AF_D1", start=23, end=136, domain_type="AF"),
        Domain(identifier="AF_D2", start=554, end=615, domain_type="AF"),
        Domain(identifier="AF_D3", start=625, end=651, domain_type="AF"),
        Domain(identifier="AF_D4", start=699, end=722, domain_type="AF"),
        Domain(identifier="AF_D5", start=724, end=752, domain_type="AF"),
        Domain(identifier="AF_D6", start=761, end=812, domain_type="AF"),
        Domain(identifier="AF_D7", start=927, end=951, domain_type="AF"),
        Domain(identifier="AF_D8", start=988, end=1050, domain_type="AF"),
        Domain(identifier="AF_D9", start=1141, end=1186, domain_type="AF"),
        Domain(identifier="AF_D10", start=1323, end=1366, domain_type="AF"),
        Domain(identifier="AF_D11", start=1373, end=1386, domain_type="AF"),
        Domain(identifier="AF_D12", start=1401, end=1412, domain_type="AF"),
        Domain(identifier="AF_D13", start=1419, end=1435, domain_type="AF"),
        Domain(identifier="AF_D14", start=1455, end=1467, domain_type="AF"),
        Domain(identifier="AF_D15", start=1586, end=1600, domain_type="AF"),
        Domain(identifier="AF_D16", start=1611, end=1627, domain_type="AF"),
        Domain(identifier="AF_D17", start=1723, end=1728, domain_type="AF"),
        Domain(identifier="Disordered", start=170, end=212, domain_type="UniProt"),
        Domain(identifier="Disordered", start=260, end=370, domain_type="UniProt"),
        Domain(identifier="Disordered", start=496, end=520, domain_type="UniProt"),
        Domain(identifier="Disordered", start=807, end=831, domain_type="UniProt"),
        Domain(identifier="Disordered", start=843, end=906, domain_type="UniProt"),
        Domain(identifier="Disordered", start=967, end=989, domain_type="UniProt"),
        Domain(identifier="Disordered", start=1046, end=1148, domain_type="UniProt"),
        Domain(identifier="Disordered", start=1184, end=1247, domain_type="UniProt"),
        Domain(identifier="Disordered", start=1477, end=1523, domain_type="UniProt"),
        Domain(identifier="Disordered", start=1537, end=1580, domain_type="UniProt"),
        Domain(identifier="Disordered", start=1624, end=1655, domain_type="UniProt"),
        Domain(identifier="Disordered", start=1683, end=1728, domain_type="UniProt"),
        Domain(identifier="COILED", start=59, end=130, domain_type="UniProt"),
        Domain(identifier="COILED", start=1156, end=1183, domain_type="UniProt"),
        Domain(identifier="Basic and acidic residues", start=267, end=286, domain_type="UniProt"),
        Domain(identifier="Pro residues", start=309, end=323, domain_type="UniProt"),
        Domain(identifier="Polar residues", start=350, end=369, domain_type="UniProt"),
        Domain(identifier="Basic and acidic residues", start=869, end=892, domain_type="UniProt"),
        Domain(identifier="Basic and acidic residues", start=1053, end=1078, domain_type="UniProt"),
        Domain(identifier="Polar residues", start=1099, end=1148, domain_type="UniProt"),
        Domain(identifier="Polar residues", start=1184, end=1230, domain_type="UniProt"),
        Domain(identifier="Polar residues", start=1696, end=1720, domain_type="UniProt"),
    ]

    # Create the Protein object for Ana1
    ana1_protein = Protein("Ana1", "ana1_acc_id", ana1_sequence)
    for domain in domains:
        ana1_protein.add_domain(domain)

    # Fragmentation parameters
    length = {'min': 100, 'ideal': 150, 'max': 200}
    overlap = {'min': 20, 'ideal': 30, 'max': 40}

    # Run recursive fragmentation with a time limit to force `break_in_half`
    result = recursive_fragmentation(ana1_protein, domains, 0, length, overlap, time_limit=0.01)

    # Check if the fragmentation timed out, and use break_in_half if necessary
    if result == "TIME_LIMIT_EXCEEDED":
        first_half, second_half = break_in_half(ana1_protein, length, overlap)

        # Ensure the halves are valid and their lengths are above the minimum length
        assert first_half is not None and second_half is not None, "Failed to split the protein into halves"
        assert (first_half.last_res - first_half.first_res + 1) >= length['min'], "First half is smaller than the minimum fragment length"
        assert (second_half.last_res - second_half.first_res + 1) >= length['min'], "Second half is smaller than the minimum fragment length"

        # Check if any domain is split between the first and second halves
        for domain in domains:
            if first_half.first_res <= domain.start < first_half.last_res and domain.end > first_half.last_res:
                # Domain starts in the first half and ends in the second half, which means it's split
                assert False, f"Domain {domain.id} was split between the first and second halves (domain = {domain.start}-{domain.end}, first = {first_half.first_res}-{first_half.last_res}, second = {second_half.first_res}-{second_half.last_res})"
    else:
        pytest.fail("Expected TIME_LIMIT_EXCEEDED but fragmentation completed without timeout")

def mock_recursive_fragmentation_timeout(*args, **kwargs):
    """
    Mock function that simulates a timeout during recursive fragmentation,
    forcing the system to invoke `break_in_half` instead of completing normally.
    """
    return "TIME_LIMIT_EXCEEDED"

def test_two_rounds_of_fragmentation():
    """
    Tests two rounds of fragmentation by mocking a timeout in the first round,
    ensuring that the start and end of the protein are covered exactly once in the final fragments
    and that there are overlaps between the fragments.
    """
    # Create a long protein sequence (~2000 residues)
    sequence = 'A' * 2000
    protein = Protein("Protein1", "example_acc_id", sequence)

    # Fragmentation parameters
    length = {'min': 100, 'ideal': 150, 'max': 200}
    overlap = {'min': 20, 'ideal': 30, 'max': 40}

    # First, we mock the `recursive_fragmentation` function to simulate a timeout
    with mock.patch('alphafragment.fragmentation_methods.recursive_fragmentation', side_effect=mock_recursive_fragmentation_timeout):
        
        # Perform the first round of break_in_half
        first_half, second_half = break_in_half(protein, length, overlap)

        # Perform the second round of break_in_half on the resulting halves
        sub_first_half_1, sub_second_half_1 = break_in_half(first_half, length, overlap)
        sub_first_half_2, sub_second_half_2 = break_in_half(second_half, length, overlap)

    # Ensure the original start and end residues are covered exactly once in the final fragments
    assert sub_first_half_1.first_res == 0, "Start residue not covered correctly."
    assert sub_second_half_2.last_res == len(sequence) - 1, "End residue not covered correctly."

    # Ensure the entire sequence is covered
    assert sub_first_half_1.first_res == protein.first_res, "First residue of the original sequence not covered."
    assert sub_second_half_2.last_res == protein.last_res, "Last residue of the original sequence not covered."

    # Check that the four resulting fragments together cover the entire sequence
    fragments = [
        (sub_first_half_1.first_res, sub_first_half_1.last_res),
        (sub_second_half_1.first_res, sub_second_half_1.last_res),
        (sub_first_half_2.first_res, sub_first_half_2.last_res),
        (sub_second_half_2.first_res, sub_second_half_2.last_res),
    ]
    
    # Sort fragments by starting position
    fragments.sort(key=lambda x: x[0])

    # Ensure overlaps between fragments
    for i in range(len(fragments) - 1):
        current_fragment_end = fragments[i][1]
        next_fragment_start = fragments[i + 1][0]

        # Ensure overlap
        assert current_fragment_end - overlap['min'] >= next_fragment_start, (
            f"Fragments {i} and {i+1} do not overlap correctly: "
            f"end of fragment {i} = {current_fragment_end}, start of fragment {i+1} = {next_fragment_start}"
        )
        assert current_fragment_end - overlap['max'] <= next_fragment_start, (
            f"Fragments {i} and {i+1} have an overlap larger than allowed: "
            f"end of fragment {i} = {current_fragment_end}, start of fragment {i+1} = {next_fragment_start}"
        )

    # Ensure the fragments cover the full sequence
    assert fragments[0][0] == protein.first_res, "Fragments do not start at the beginning of the sequence."
    assert fragments[-1][1] == protein.last_res, "Fragments do not end at the last residue of the sequence."
