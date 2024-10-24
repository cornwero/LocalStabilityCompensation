import argparse
import sys
import bpRNAStructure as ST
"""
filename: generateBulgeData.py
Author: Micheal Hathaway
Updated: Robert Cornwell-Arquitt

Description: Python script to gather bulge data for testing local
            energy compensation hypothesis using the bpRNAStructure module.

Usage: python3 generateHairpinData.py <RNA_Identifiers_File> <Input_Data_File>
        <Output_file_name>
"""


def get_rna_identifiers(rna_identifiers_file):
    """
    Function: get_rna_identifiers(filename)
    Description: Function to read in RNA Identifier data in from external file
                and store as a python dictionary.
    Parameters:
        (rna_identifiers_file) - str - filename containing the RNA Identifiers
    Return Value:
        Dictionary of RNA Identifiers
    """
    # create python dictionary to hold RNA Identifier data
    rna_identifiers = {}

    # open file containing RNA Identifiers
    with open(rna_identifiers_file, 'r') as f:
        for line in f:
            data = line.strip().split('\t')  # strip and split on tabs
            if data[0] != 'rna_name':
                # add data to the rna_identifiers dictionary
                rna_identifiers[data[0]] = {
                    'rna_name': data[0], 'rna_id': data[1], 'rna_type': data[2]}

        return rna_identifiers


def parse_arguments():
    """
    Function: parse_arguments()
    Description: Function to parse command line arguments provided by the user.
    Parameters:
            None
    Return Value:
            tuple(rna_identifiers_File, Input_Data, output_file)
    """
    parser = argparse.ArgumentParser(
        description="Script to generate Local Energy hairpin data.")
    parser.add_argument('rna_identifiers_File',
                        help="Specify RNA Identifiers File to be parsed.",
                        type=str)
    parser.add_argument(
        'Input_Data', help="Structure Type file. txt file of .st file names \
        will also work", type=str)
    parser.add_argument('Output_file_name',
                        help="Specify name of new file", type=str)

    args = parser.parse_args()
    return (args.rna_identifiers_File, args.Input_Data, args.Output_file_name)


def multipleFiles(input_name):
    """
    Function: multipleFile(input_name)
    Description: Function returns True if the input file is a txt file of
                multiple file paths to .st files and False is a single .st
                file is provided.
    Parameters:
        (input_name) -- str -- input filename provided by the user as a
                            command line argument
    Return Value:
            bool
    """
    if input_name[-3:] == '.st':
        return False
    elif input_name[-4:] == '.txt':
        return True
    else:
        print('Unrecognized input data format.')
        sys.exit()


def create_output_file(output_file):
    """
    Function: create_output_file(output_file)
    Description: Function to create the file that data will be outputted to
                and prep file with correct data headers
    Parameters:
        (output_file) -- str -- name of the file to output data to.
    Return Value:
        None
    """
    with open(output_file, 'w') as f:
        headers = 'rna_name' + '\t' + 'rna_id' + '\t' + 'rna_type' + '\t' + \
            'hairpinLabel' + '\t' + 'hairpinLen' + '\t' + 'hairpinEnergy' + '\t' + \
            'stemLabel' + '\t' + 'stemLen' + '\t' + 'stemEnergy' + '\n'

        f.write(headers)


def extractDataFromHairpin(hairpin, structureObject):
    """
    Function: extractDataFromHairpin(hairpin, structureObject)
    Description: Function to extract data from a given hairpin object
    Parameters:
        (hairpin) -- hairpin object -- hairpin of interest
        (structureObject) -- Structure object -- Structure object that the
                                bulge is a part of
    Return Value:
        tuple(bulgeLabel, bulgeLen, bulgeEnergy, energy5p, energy3p) OR None
    """
    if hairpin.canonical():  # check for valid hairpin parameter
        hairpinLabel = hairpin.label()  # store hairpin label
        hairpinEnergy = round(hairpin.energy(), 3)  # store hairpin energy
        hairpinLen = hairpin.sequenceLen()  # store length of hairpin sequence

        neighbor5p, neighbor3p = structureObject.neighbors(
            hairpinLabel, object=True)  # get the 5' and 3' neighbor of the hairpin
        # check for valid energy parameters in both neighbors
        if (neighbor5p.canonical()) and (neighbor3p.canonical()):
            # round energy values to 3 decimal places
            energy5p = round(neighbor5p.energy(), 3)
            label5p = neighbor5p.label()
            # round energy values to 3 decimal places
            energy3p = round(neighbor3p.energy(), 3)
            label3p = neighbor3p.label()

            # only return data if all values can be accuratley calculated
            return (hairpinLabel, hairpinLen, hairpinEnergy, label5p, energy5p,
                    label3p, energy3p)

    return None


def get_data_from_structure(structure, rna_identifiers):
    """
    Function: get_data_from_structure(structure, rna_identifiers)
    Description: Function to extract data from Structure object
    Parameters:
        (structure) -- Structure object -- Structure object that is
                            being examined
        rna_identifiers) -- dict -- dictionary of RNA Identifier data
    Return Value:
        tuple(rna_name, rna_id, rna_type, bulges)
    """
    rna_name = structure.name()  # get rna_name from structure
    rna_id = rna_identifiers[rna_name]['rna_id']
    rna_type = rna_identifiers[rna_name]['rna_type']

    # get a list of all the bulges in the Structure object
    hairpins = structure.hairpins()

    return (rna_name, rna_id, rna_type, hairpins)


def process_single_structure(filename, rna_identifiers, output_file):
    """
    Function: process_single_structure(filename, rna_identifiers, output_file)
    Description: Function to process a single .st file and write all hairpin data
                 to a file
    Parameters:
        (filename) -- str -- .st filename
        (rna_identifiers) -- dict -- dictionary of RNA Identifier data
        (output_file) -- str -- name of file to output data to
    Return Value:
        None
    """
    try:
        structureObject = ST.Structure(filename)  # create Structure object
    except Exception as e:
        print(
            f'An error occurred when creating a Structure object \
            for {filename}: {e}.')
        return

    rna_name, rna_id, rna_type, hairpins = get_data_from_structure(
        structureObject, rna_identifiers)  # get data from Structure Object

    with open(output_file, 'a') as f:
        for hairpin in hairpins:
            hairpinData = extractDataFromHairpin(hairpin, structureObject)
            if hairpinData:  # make sure all energy data is present
                hairpinLabel, hairpinLen, hairpinEnergy, label5p, energy5p, label3p, energy3p = hairpinData
                netE = energy5p+hairpinEnergy
                dataString = f'{rna_name}\t{rna_id}\t{rna_type}\t{hairpinLabel}\t{hairpinLen}\t{hairpinEnergy}\t{label5p}\t{energy5p}\t{netE}\n'
                f.write(dataString)


def process_multiple_files(file_list, rna_identifiers, output_file):
    """
    Function: process_multiple_files(file_list, rna_identifiers, output_file)
    Description: Function to process a list of .st file and write all hairpin
                data to a file
    Parameters:
            (file_list) -- str -- .txt file of .st filenames to collect data for.
            (rna_identifiers) -- dict -- dictionary of RNA Identifier data
            (Output_file_name) -- str -- name of file to output data to
    Return Value: None
    """
    with open(file_list, 'r') as f:
        for line in f:
            filename = 'stFiles/' + line.strip()
            processSingleStructure(filename, RNA_Identifiers, Output_Data_File)


if __name__ == '__main__':
    # get command line arguments
    rna_identifiers_file, input_file, output_file = parse_arguments()
    # parse RNA Identifier data
    rna_identifiers = get_rna_identifiers(rna_identifiers_file)

    # create file to output data and prep with appropriate data headers
    create_output_file(output_file)

    # check if multiple .st files are provided in .txt file, or just one
    multipleFiles = multipleFiles(input_file)
    if(multipleFiles):
        process_multiple_files(
            input_file, rna_identifiers, output_file)
    else:
        process_single_structure(
            input_file, rna_identifiers, output_file)
