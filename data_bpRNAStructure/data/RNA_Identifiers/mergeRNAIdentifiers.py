"""
filename: mergeRNA_ID_type.py
Author: Michael Hathaway

Descrition: python script to merge allbpRNA_OriginURLs.txt and cleaned_typelist.txt into a single file containing bpRNA ID, rnaID, and rna type
for each RNA moleculue in the bpRNA database.

Usage: python3 mergeRNA_ID_Type.py <allbpRNA_OriginURLs.txt> <cleaned_typelist.txt> <outputfile>
"""

#Imports
import argparse


'''
Function: parseArgs()
Description: Function to handle command line arguments
parameters: None
Return Value: tuple containing the two input file names and the output file name
'''
def parseArgs():
    parser = argparse.ArgumentParser(description="Script to combine.merge allbpRNA_OriginURLs.txt and cleaned_typelist.txt into a single file containing bpRNA ID, rnaID, and rna type for each RNA moleculue in the bpRNA database.")
    parser.add_argument('Original_URL_file', help="Specify URL file to be parsed.", type=str)
    parser.add_argument('RNA_Type_file', help="Specify file containing all RNA types.", type=str)
    parser.add_argument('Output_file_name', help="Specify name of new file", type=str)

    args = parser.parse_args()
    return (args.Original_URL_file, args.RNA_Type_file, args.Output_file_name)


"""
Function: parseURLFile(filename, combinedData)
Description: Function to parse the allbpRNA_OriginURLs.txt file to get rna_ID and rna_name
parameters:
        (filename) - str - file to be parsed(allbpRNA_OriginURLs.txt)
        (combinedData) - dict - dictionary containing rna_ID, rna_name, and rna_type data
Return Value:
        None
"""
def parseURLFile(filename, combinedData):
    with open(filename, 'r') as f:
        for line in f:
            data = line.strip().split('\t')
            if data[0] != 'rna_ID': #skip column headers
                combinedData[data[0]] = {}
                combinedData[data[0]]['rna_name'] = data[1]



"""
Function: parseRNATypesFile(filename, combinedData)
Description: Function to parse the cleaned_typelist.txt file to get rna_type
parameters:
        (filename) - str - file to be parsed(allbpRNA_OriginURLs.txt)
        (combinedData) - dict - dictionary containing rna_ID, rna_name, and rna_type data
Return Value:
        None
"""
def parseRNATypesFile(filename, combinedData):
    with open(filename, 'r') as f:
        for line in f:
            data = line.strip().split('\t')
            if data[0] in combinedData:
                combinedData[data[0]]['rna_ID'] = data[0]
                combinedData[data[0]]['rna_type'] = data[1]


"""
Function: writeToCombinedFile(outputFileName, combinedData)
Description: Function to write the data from the combinedData dictionary to a file that contains rna_name, rna_ID, and rna_type data.
parameters:
        (outputFileName) - str - file to output data to
        (combinedData) - dict - dictionary containing rna_ID, rna_name, and rna_type data
Return Value:
        None
"""
def writeToCombinedFile(outputFileName, combinedData):
    with open(outputFileName, 'w') as f:
        headers = 'rna_name' + '\t' + 'rna_ID' + '\t' + 'rna_type' + '\n'
        f.write(headers)
        for rna in combinedData.keys():
            output = combinedData[rna]['rna_name'] + '\t' + combinedData[rna]['rna_ID'] + '\t' + combinedData[rna]['rna_type'] + '\n'
            f.write(output)


"""
Main Function
"""
if __name__ == '__main__':
    combinedData = {} #create dictionary to store the combined date of the two files
    args = parseArgs() #get command line arguments

    parseURLFile(args[0], combinedData) #get rna_name and rna_id
    parseRNATypesFile(args[1], combinedData) #get rna_type

    writeToCombinedFile(args[2], combinedData) #write data to output file
