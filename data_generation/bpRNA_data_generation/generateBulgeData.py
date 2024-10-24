"""
filename: generateBulgeData.py
Author: Michael Hathaway
Updated: Robert Cornwell-Arquitt

Description: Python script to gather bulge data for testing local stability compensation hypothesis using the bpRNAStructure module.

Usage: python generateBulgeData.py <RNA_Identifiers_File> <Input_Data_File> <Output_file_name>

"""

"""
Imports
"""
import argparse
import logging
import sys
import bpRNAStructure.Structure as ST

"""
Function: getRNAIdentifiers(filename)
Description: Function to read in RNA Identifier data in from external file and store as a python dictionary.
Parameters:
        (RNAIdentifiersFile) -- str -- filename containing the RNA Identifiers
Return Value:
        Dictionary of RNA Identifiers
"""
def getRNAIdentifiers(RNAIdentifiersFile):
    RNAIdentifiers = {} #create python dictionary to hold RNA Identifier data
    with open(RNAIdentifiersFile, 'r') as f: #open file containing RNA Identifiers
        for line in f:
            data = line.strip().split('\t') #strip each line in file and split on tabs
            if data[0] != 'rna_name':
                RNAIdentifiers[data[0]] = {'rna_name': data[0], 'rna_ID': data[1], 'rna_Type':data[2]} #add data to the RNAIdentifiers dictionary

        return RNAIdentifiers


"""
Function: parseArguments()
Description: Function to parse command line arguments provided by the user.
Parameters:
        None
Return Value:
        tuple(RNA_Identifiers_File, Input_Data, Output_Data_File)
"""
def parseArguments():
    parser = argparse.ArgumentParser(description="Script to generate Local Energy bulge data.")
    parser.add_argument('RNA_Identifiers_File', help="Specify RNA Identifiers File to be parsed.", type=str)
    parser.add_argument('Input_Data', help="Structure Type file. txt file of .st file names will also work", type=str)
    parser.add_argument('Output_file_name', help="Specify name of new file", type=str)

    args = parser.parse_args()
    return (args.RNA_Identifiers_File, args.Input_Data, args.Output_file_name)


"""
Function: multipleFile(inputData)
Description: Function returns True if the input file is a txt file of multiple file paths to .st files and False is a single .st file is provided.
Parameters:
        (inputData) -- str -- input filename provided by the user as a command line argument
Return Value:
        bool
"""
def multipleFiles(inputData):
    if inputData[-3:] == '.st':
        return False
    elif inputData[-4:] == '.txt':
        return True
    else:
        print('Unrecognized input data format.')
        sys.exit()


"""
Function: createOutputDataFile(Output_Data_File)
Description: Function to create the file that data will be outputted to and prep file with correct data headers
Parameters:
        (Output_Data_File) -- str -- name of the file to output data to.
Return Value:
        None
"""
def createOutputDataFile(Output_Data_File):
    with open(Output_Data_File, 'w') as f:
        headers = 'rna_name\trna_ID\trna_Type\tbulgeLabel\tbulgeLen\tbulgeEnergy\tlabel5p\tenergy5p\tlabel3p\tenergy3p\tnetE' + '\n'

        f.write(headers)


"""
Function: extractDataFromBulge(bulge, structureObject)
Description: Function to extract data from a given bulge object
Parameters:
        (bulge) -- Bulge object -- bulge of interest
        (structureObject) -- Structure object -- Structure object that the bulge is a part of
Return Value:
        tuple(bulgeLabel, bulgeLen, bulgeEnergy, energy5p, energy3p) OR None
"""
def extractDataFromBulge(bulge, structureObject):
    if bulge.canonical(): #check for valid bulge parameter
        bulgeLabel = bulge.label() #store bulge label
        bulgeEnergy = round(bulge.energy(), 3) #store bulge energy
        bulgeLen = bulge.sequenceLen() #store length of bulge sequence

        neighbor5p, neighbor3p = structureObject.neighbors(bulgeLabel, object=True) #get the 5' and 3' neighbor of the bulge
        if (neighbor5p.canonical()) and (neighbor3p.canonical()): #check for valid energy parameters in both neighbors
            energy5p = round(neighbor5p.energy(), 3) #round energy values to 3 decimal places
            label5p = neighbor5p.label()
            energy3p = round(neighbor3p.energy(), 3) #round energy values to 3 decimal places
            label3p = neighbor3p.label()

            #only return data if all values can be accuratley calculated
            return (bulgeLabel, bulgeLen, bulgeEnergy, label5p, energy5p, label3p, energy3p)

    return None


"""
Function: extractDataFromStructure(structureObject, RNA_Identifiers)
Description: Function to extract data from Structure object
Parameters:
        (structureObject) -- Structure object -- Structure object that is being examined
        (RNA_Identifiers) -- dict -- dictionary of RNA Identifier data
Return Value:
        tuple(rna_name, rna_ID, rna_Type, bulges)
"""
def extractDataFromStructure(structureObject, RNA_Identifiers):
    rna_name = structureObject.name() #get rna_name from structureObject
    rna_ID = RNA_Identifiers[rna_name]['rna_ID']
    rna_Type = RNA_Identifiers[rna_name]['rna_Type']
    bulges = structureObject.bulges() #get a list of all the bulges in the Structure object

    return (rna_name, rna_ID, rna_Type, bulges)


"""
Function: processSingleStructure(structureTypeFile, RNA_Identifiers, Output_file_name)
Description: Function to process a single .st file and write all bulge data to a file
Parameters:
        (structureTypeFile) -- str -- .st filename
        (RNA_Identifiers) -- dict -- dictionary of RNA Identifier data
        (Output_file_name) -- str -- name of file to output data to
Return Value:
        None
"""
def processSingleStructure(structureTypeFile, RNA_Identifiers, Output_file_name):
    try:
        structureObject = ST.Structure(structureTypeFile) #create Structure object
    except Exception as e:
        print(f'An error occurred when creating a Structure object for {structureTypeFile}: {e}.')
        return

    rna_name, rna_ID, rna_Type, bulges = extractDataFromStructure(structureObject, RNA_Identifiers) #get data from Structure Object

    with open(Output_file_name, 'a') as f:
        for bulge in bulges:
            bulgeData = extractDataFromBulge(bulge, structureObject)
            if bulgeData: #make sure all energy data is present
                bulgeLabel, bulgeLen, bulgeEnergy, label5p, energy5p, label3p, energy3p = bulgeData
                netE = (energy5p+energy3p)/2 + bulgeEnergy
                dataString = f'{rna_name}\t{rna_ID}\t{rna_Type}\t{bulgeLabel}\t{bulgeLen}\t{bulgeEnergy}\t{label5p}\t{energy5p}\t{label3p}\t{energy3p}\t{netE}\n'
                f.write(dataString)
           # else:
           #    f.write('Bad stem or bulge\n')

"""
Function: processMultipleFiles(fileList, RNA_Identifiers, Output_Data_File)
Description: Function to process a list of .st file and write all bulge data to a file
Parameters:
        (fileList) -- str -- .txt file of file paths to .st files
        (RNA_Identifiers) -- dict -- dictionary of RNA Identifier data
        (Output_file_name) -- str -- name of file to output data to
Return Value: None
"""
def processMultipleFiles(fileList, RNA_Identifiers, Output_Data_File):
    with open(fileList, 'r') as f:
        for line in f:
            filename = 'stFiles/' + line.strip()
            processSingleStructure(filename, RNA_Identifiers, Output_Data_File)



"""
Main
"""
if __name__ == '__main__':
    #get command line arguments
    RNA_Indentifiers_file, Input_Data_File, Output_Data_File = parseArguments()
    #parse RNA Identifier data
    RNA_Identifiers = getRNAIdentifiers(RNA_Indentifiers_file)

    #create file to output data and prep with appropriate data headers
    createOutputDataFile(Output_Data_File)

    #check if multiple .st files are provided in .txt file, or just one
    multipleFiles = multipleFiles(Input_Data_File)
    if(multipleFiles):
        processMultipleFiles(Input_Data_File, RNA_Identifiers, Output_Data_File)
    else:
        processSingleStructure(Input_Data_File, RNA_Identifiers, Output_Data_File)
