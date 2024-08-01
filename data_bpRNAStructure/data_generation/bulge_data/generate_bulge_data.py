import argparse
import sys
"""
filename: generateBulgeData.py
Author: Michael Hathaway

Description: Python script to gather bulge data for testing local
            energy compensation hypothesis using the bpRNAStructure module.

Usage: python3 bulgeLocalEnergy.py <RNA_Identifiers_File> <Input_Data_File>
        <Output_file_name>
"""


def path_setup():
    """
    This is to get around issue with local imports. These are the paths to
    bpRNAStructue and to the Turner Parameters
    """
    sys.path.insert(0, '/nfs0/BB/Hendrix_Lab/bpRNA/LocalEnergy/bpRNAStructure')
    sys.path.insert(
        0, '/nfs0/BB/Hendrix_Lab/bpRNA/LocalEnergy/bpRNAStructure/TurnerParameters/parameters')


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
                try:
                    rna_identifiers[data[0]] = {
                        'rna_name': data[0], 'rna_id': data[1], 'rna_type': data[2]}
                except Exception as e:
                    print(f"Error for: {data}")

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
        description="Script to generate Local Energy bulge data.")
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
            'bulgeLabel' + '\t' + 'bulgeLen' + '\t' + 'bulgeEnergy' + '\t' + \
            'label5p' + '\t' + 'energy5p' + '\t' + 'label3p' + '\t' + \
            'energy3p' + '\n'

        f.write(headers)


def extractDataFromBulge(bulge, structureObject):
    """
    Function: extractDataFromBulge(bulge, structureObject)
    Description: Function to extract data from a given bulge object
    Parameters:
        bulge) -- Bulge object -- bulge of interest
        (structureObject) -- Structure object -- Structure object that the
                                bulge is a part of
    Return Value:
        tuple(bulgeLabel, bulgeLen, bulgeEnergy, energy5p, energy3p) OR None
    """
    if bulge.canonical():  # check for valid bulge parameter
        bulgeLabel = bulge.label()  # store bulge label
        bulgeEnergy = round(bulge.energy(), 3)  # store bulge energy
        bulgeLen = bulge.sequenceLen()  # store length of bulge sequence

        neighbor5p, neighbor3p = structureObject.neighbors(
            bulgeLabel, object=True)  # get the 5' and 3' neighbor of the bulge
        # check for valid energy parameters in both neighbors
        if (neighbor5p.canonical()) and (neighbor3p.canonical()):
            # round energy values to 3 decimal places
            energy5p = round(neighbor5p.energy(), 3)
            label5p = neighbor5p.label()
            # round energy values to 3 decimal places
            energy3p = round(neighbor3p.energy(), 3)
            label3p = neighbor3p.label()

            # only return data if all values can be accuratley calculated
            return (bulgeLabel, bulgeLen, bulgeEnergy, label5p, energy5p,
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
    bulges = structure.bulges()

    return (rna_name, rna_id, rna_type, bulges)


def process_single_structure(filename, rna_identifiers, output_file):
    """
    Function: process_single_structure(filename, rna_identifiers, output_file)
    Description: Function to process a single .st file and write all bulge data
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

    rna_name, rna_id, rna_type, bulges = get_data_from_structure(
        structureObject, rna_identifiers)  # get data from Structure Object

    with open(output_file, 'a') as f:
        for bulge in bulges:
            bulgeData = extractDataFromBulge(bulge, structureObject)
            if bulgeData:  # make sure all energy data is present
                bulgeLabel, bulgeLen, bulgeEnergy, label5p, energy5p, label3p, energy3p = bulgeData
                dataString = f'{rna_name}\t{rna_id}\t{rna_type}\t{bulgeLabel}\t{bulgeLen}\t{bulgeEnergy}\t{label5p}\t{energy5p}\t{label3p}\t{energy3p}\n'
                f.write(dataString)


def process_multiple_files(file_list, rna_identifiers, output_file):
    """
    Function: process_multiple_files(file_list, rna_identifiers, output_file)
    Description: Function to process a list of .st file and write all bulge
                data to a file
    Parameters:
            (file_list) -- str -- .txt file of file paths to .st files
            (rna_identifiers) -- dict -- dictionary of RNA Identifier data
            (Output_file_name) -- str -- name of file to output data to
    Return Value: None
    """
    with open(file_list, 'r') as f:
        for line in f:
            filename = line.strip()
            process_single_structure(
                filename, rna_identifiers, output_file)


if __name__ == '__main__':
    path_setup()
    import Structure as ST
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
