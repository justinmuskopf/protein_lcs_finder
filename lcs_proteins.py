"""
    Authors: Justin Muskopf, Aaron Johnson
    Course: CSCE 4110, Fall 2018
    Assignment: 2
    REQUIRED: Sequences file that formats every sequence like so:
        >SEQUENCE_NAME
        ...SEQUENCE__DATA...
      Where every sequence is followed by a blank line
"""
from sys import argv
from os.path import isfile

# Constants
LENGTH   = 0
NAME     = 0
SEQUENCE = 1


def die(die_string = None, code = 1):
    if die_string:
        print("Error: {}".format(die_string))
    else:
        print("Error: An Undefined Error Occurred!")

    exit(code)


def get_sequences_from_file(filename):
    try:
        with open(filename, "r") as f:
            file_string = f.read()
    except IOError:
        die("Error! Could not open {}".format(filename))

    # Split by empty line and check length
    sequences = file_string.split("\n\n")
    if len(sequences) == 0:
        die("Could not find any sequences!")

    # Remove class descriptions
    if sequences[0][0] != '>':
        sequences = sequences[1:]

    # Get the sequence data for all sequences
    sequences = [get_sequence_data(i) for i in sequences if i != '']

    return sequences


# Get the LCS of two sequences; O(mn) || O(n^2) (if same len)
def lcs(sequence_1, sequence_2):
    # find the length of the strings
    sequence_1_len = len(sequence_1)
    sequence_2_len = len(sequence_2)

    # Will hold the LCS substrings
    lcs_array = [[""] * (sequence_2_len + 1) for _ in range(sequence_1_len + 1)]

    # Fill lcs_array
    for i in range(sequence_1_len + 1):
        for j in range(sequence_2_len + 1):
            if i == 0 or j == 0:
                continue

            # Get each sequence's characters for comparison
            sequence_1_char = sequence_1[i - 1]
            sequence_2_char = sequence_2[j - 1]

            # Match!
            if sequence_1_char == sequence_2_char:
                lcs_array[i][j] = lcs_array[i - 1][j - 1] + sequence_1_char
                continue
            # No Match
            else:
                # Get current string from lcs_array
                first_string = lcs_array[i - 1][j]
                first_len = len(first_string)

                # Get current string from lcs_array
                second_string = lcs_array[i][j - 1]
                second_len = len(second_string)

                # Assign longest string
                lcs_array[i][j] = first_string if first_len > second_len else second_string

    # The LCS string is stored at the last index of the array
    lcs_string = lcs_array[sequence_1_len][sequence_2_len]

    return len(lcs_string), lcs_string


def get_sequence_data(sequence):
    # Name ends at first newline
    first_newline_index = sequence.find('\n')

    # Remove > and get name
    name = sequence[1:first_newline_index]

    # Get the rest as sequence
    sequence = sequence[first_newline_index:]

    # Remove newlines from sequence
    sequence = sequence.replace("\n", "")

    return name, sequence


def get_sequence_max_lcs(sequence):
    max_lcs = 0
    for key in sequence:
        # Get sequence object and LCS length
        sequence_object = sequence[key]
        lcs_length = sequence_object[LENGTH]

        # New max
        if lcs_length > max_lcs:
            max_lcs_object = (key, sequence_object)
            max_lcs = lcs_length

    return max_lcs_object


def add_to_group(name, to_add, groups):
    if name not in groups:
        groups[name] = []
    groups[name].append(to_add)


def get_grouped_sequences(sequences):
    groups = {}

    for sequence in sequences:
        max_sequence_name, max_sequence_object = get_sequence_max_lcs(sequences[sequence])

        to_add = (sequence, max_sequence_object[LENGTH], max_sequence_object[SEQUENCE])

        add_to_group(max_sequence_name, to_add, groups)

    return groups


def get_filename(args):
    num_args = len(args)

    # Get filename
    if num_args < 2:
        filename = input("Please enter the filename containing protein sequences: ")
    else:
        filename = argv[1]

    # Die if not a file
    if not isfile(filename):
        die("Could not find file '{}'".format(filename))

    return filename


def lcs_all_sequences(sequences):
    # Create dictionary of sequences
    lcs_dict = {}
    for name, _ in sequences:
        lcs_dict[name] = {}

    num_sequences = len(sequences)

    # For all sequences
    for sequence_idx, first_sequence_data in enumerate(sequences):
        first_name = first_sequence_data[NAME]
        first_sequence = first_sequence_data[SEQUENCE]

        print("[{}/{}] {}...".format(sequence_idx + 1, num_sequences, first_name))

        # Only loop through those that have not been LCSd
        for second_name, second_sequence in sequences[sequence_idx:]:
            # Same sequence, skip
            if first_name == second_name:
                continue

            lcs_data = lcs(first_sequence, second_sequence)

            # Place LCS data under both names
            lcs_dict[first_name][second_name] = lcs_data
            lcs_dict[second_name][first_name] = lcs_data

            print("   {}... LCS Length: {}".format(second_name, lcs_data[LENGTH]))
        print()

    return lcs_dict


def print_grouped_item(item, padding=2, offset=20):
    padded = " " * padding

    name_string = "Name: {}".format(item[NAME])
    offset = max(offset - len(name_string), 0) * " "

    name_string += offset

    print("{} {} Length of LCS: {}".format(padded, name_string, item[1]))


def main(args):
    filename = get_filename(args)

    sequences = get_sequences_from_file(filename)

    lcs_sequence_dict = lcs_all_sequences(sequences)

    grouped_sequences = get_grouped_sequences(lcs_sequence_dict)

    for group in grouped_sequences:
        print("Group {}:".format(group))
        for grouped_item in grouped_sequences[group]:
            print_grouped_item(grouped_item)


if __name__ == "__main__":
    main(argv)
