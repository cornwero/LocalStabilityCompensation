'''
Filename: Structure.py
Author: Michael Hathaway

Description: python module that defines the Structure Object.
The Structure Object provides a user friendly mechanism for working with
RNA structure type files in the python programming language.
'''

## Module Imports ##
import numpy as np
import sys
import re

## Structure Type Component Imports ##
from StructureComponents import Stem, Hairpin, Bulge, InternalLoop, ExternalLoop, MultiLoop, PseudoKnot, End, NCBP

'''
## About the structure object ##
The Structure object is a python object-oriented representation of the information contained within an RNA Structure Type file.
The object provides a mechanism to easily access and work with the data in the python programming language. In addition it includes
functionality for calculating the energy associated with certain RNA secondary structures within the RNA molecule.
'''
class Structure:
    #__init__() method for the Structure object
    def __init__(self, filename=None):
        #RNA Molecule basic info
        #all values are stored as strings
        self._name = None
        self._length = None
        self._pageNum = None

        #structural representations of RNA molecule
        #all values are stored as strings
        self._sequence = None
        self._DBN = None
        self._structureArray = None
        self._varna = None

        '''
        secondary structure information
        secondary structure information for each RNA molecule is stored as a dictionary where
        the key is the label for the secondary structure and the value is an object that contains
        the data for the secondary structure with the given label.

        For information on how each secondary structure class was implemented and how to access
        secondary structure data, see the StructureTypeComponents section below.
        '''
        self._stems = {}
        self._hairpins = {}
        self._bulges = {}
        self._internalLoops = {}
        self._multiLoops = {}
        self._externalLoops = {}
        self._pk = {}
        self._ncbp = {}
        self._ends = {}

        '''
        Component Array
        The component array is a numpy array of the same length as the molecule where each index
        contains the label for the secondary structure tha that index is a part of.
        the component array is initialized as None. When the length of the molecule is
        parsed from the .st file, a numpy array of that length is generated
        '''
        self._componentArray = None

        #load data from file if file is specified by user
        if filename != None:
            self._loadFile(filename)


    #define string representation of the molecule
    def __str__(self):
        return f'RNA: {self._name}'

    #define len function for Structure object
    def __len__(self):
        return self._length


####################################################
###### Load File and Associated Functions ##########
####################################################

    '''
    Function: _resetStructure()
    Description: Function resets all Structure member variables to None. Function is called when a file cannot be completley parsed.
    Parameters:
            None
    Return Value:
            None
    '''
    def _resetStructure(self):
        #reset basic molecule information
        self._name = None
        self._length = None
        self._pageNum = None

        #reset all structure representations of the molecule
        self._sequence = None
        self._DBN = None
        self._structureArray = None
        self._varna = None

        #reset all StructureComponent dictionaries
        self._stems.clear()
        self._hairpins.clear()
        self._bulges.clear()
        self._internalLoops.clear()
        self._multiLoops.clear()
        self._externalLoops.clear()
        self._pk.clear()
        self._ncbp.clear()
        self._ends.clear()

        #reset component array
        self._componentArray = None


    '''
    Function Name: loadFile(filename)
    Description: user accessible function that can be used to load data from a structure type file into
    the Structureobject if no file is provided at object instantiation.
    Parameters:
            (filename) - str - name of the structure type file to be loaded into the object
    Return Type:
            None
    '''
    def loadFile(self, filename):
        self._loadFile(filename)


    '''
    Function Name: _loadFile(filename)
    Description: Internal method to parse the data in an RNA structure tyoe file into a Structureobject
    Parameters:
            (filename) - str, name of the file to be parsed
    Return Type:
            Structure object

    #Note: At this point, the data for segments and pseudoknots is not parsed from the .st file. Multiloop data is parsed, but the multiloop object is incomplete.
    '''
    def _loadFile(self, filename):

        # check that file is valid structure type
        if filename[-3::] != '.st':
            print('Must provide a valid structure type file.')
            return

        #try to open the provided file + error handling
        try:
            f = open(filename, 'r')
        except OSError: #error finding or opening file
            print('An error ocurred when trying to access the file. Check to make sure that the file exists and that the correct filepath was provided.')
            return
        except: #something weird happened
            print('Something unexpected ocurred when accessing the file')
            return


        #Variables to validate all features have been read
        sequenceRead = False
        dotBracketRead = False
        structureArrayRead = False
        varnaRead = False

        #iterate through all of the lines in the file
        for line in f:

            if line[0] == '#':
                #get name of RNA molecule
                if line[0:6] == '#Name:':
                    self._name = line[6:].strip()

                #get length of the RNA sequence
                elif line[0:8] == '#Length:':
                    self._length = int(line[8:].strip().strip(','))
                    self._componentArray = np.empty(self._length, dtype=object)

                #get page number for molecule
                elif line[0:12] == '#PageNumber:':
                    self._pageNum = int(line[12:].strip())

                else:
                    continue

            #get actual RNA sequence
            elif (sequenceRead == False):
                    self._sequence = line.strip() #drop the newline character
                    sequenceRead = True

            #get Dot Bracket Notation for the molecule
            elif (dotBracketRead == False):
                    self._DBN = line.strip() #drop the newline character
                    dotBracketRead = True

            #get Annotated symbol form of the molecule
            elif (structureArrayRead == False):
                    self._structureArray = line.strip() #drop the newline character
                    structureArrayRead = True

            #varna notation for the molecule
            elif (varnaRead == False):
                    self._varna = line.strip() #drop the newline characters
                    varnaRead = True

                    if (sequenceRead and dotBracketRead and structureArrayRead and varnaRead):
                        #when all identifying data has been parsed, parse the StructureComponents
                        features = f.read() #read the rest of the file into features variable
                        break
                    else:
                        print(f'File: {filename} is not proper .st format')
                        self._resetStructure() #reset the Structure object


        features = features.split('\n')[:-1] #split rest of file contents into a list of strings and drop last line(empty string)
        i = 0 #while loop allows for indexing multiple file lines ahead of current. Used for Multiloops and Internal Loops that have multiple components
        while i < (len(features)): #iterate through the individual string

            ##stems##
            if features[i][0] == 'S' and features[i][1].isdigit():
                self._parseStemData(features[i].strip().split(' '))

            ##Hairpins##
            elif features[i][0] == 'H':
                self._parseHairpinData(features[i].split(' '))

            ##Bulges##
            elif features[i][0] == 'B':
                self._parseBulgeData(features[i].split(' '))

            ##Inner Loops##
            elif features[i][0] == 'I' and re.search('I\d{1,3}.1', features[i]):
                self._parseInternalLoopData(features[i].split(' '), features[i+1].split(' ')) #pass both inner loop components

            ##MultiLoops##
            elif features[i][0] == 'M':
                parentLabel = self._getMultiloopParentLabel(features[i]) #get parent label of the multiloop
                subcomponents = [] #array to temporarily store multiloop subcomponents
                while self._getMultiloopParentLabel(features[i]) == parentLabel: #linear probe for other multiloop subcomponents
                    subcomponents.append(features[i].split(' ')) #append each subcomponent to the subcomponents list
                    i += 1

                self._parseMultiLoopData(subcomponents) #parse the entire multiloop subcomponent list
                continue #once all components parsed, continue to next iteration without affecting counter

            ##external loops##
            elif features[i][0] == 'X':
                self._parseExternalLoopData(features[i].split(' '))

            ##NCBP##
            elif features[i][0:4] == 'NCBP':
                self._parseNCBPData(features[i].split(' '))

            ##Ends##
            elif features[i][0] == 'E':
                self._parseEndData(features[i].strip().split(' '))

            i += 1 #increment counter

        #add neighbors to structure component Objects
        self._addStructureComponentNeighbors()
        #add stem neighboring bulge boolean controls
        self._addStemBulgeNeighborBooleans()

        f.close() #close the file



    '''
    Function Name: _parseStemData()
    Description: Internal method used by _loadFile() to parse all the stem information in the
    structure type file into the Structureobject
    Parameters:
            (stemData) - list - list of data from a line in the structure type file describing a stem
    Return Type:
            None
    '''
    def _parseStemData(self, stemData):
        #get stem Label
        stemLabel = stemData[0]

        #get start of first segment
        part5p_start = ''
        for char in stemData[1]:
            if char.isnumeric():
                part5p_start += char
            else:
                break
        part5p_start = int(part5p_start)

        #get stop of first segment
        part5p_stop = ''
        for char in reversed(stemData[1]):
            if char.isnumeric():
                part5p_stop += char
            else:
                break
        part5p_stop = int(part5p_stop[::-1])

        #get sequence for first segment
        part5p_seq = ''
        for char in stemData[2]:
            if char.isalpha():
                part5p_seq += char

        #get start of second segment
        part3p_start = ''
        for char in stemData[3]:
            if char.isnumeric():
                part3p_start += char
            else:
                break
        part3p_start = int(part3p_start)

        #get stop of second segment
        part3p_stop = ''
        for char in reversed(stemData[3]):
            if char.isnumeric():
                part3p_stop += char
            else:
                break
        part3p_stop = int(part3p_stop[::-1])

        #get sequence for second segment
        part3p_seq = ''
        for char in stemData[4]:
            if char.isalpha():
                part3p_seq += char

        #add data to the stems dictionary
        newStem = Stem(stemLabel, part5p_seq, part3p_seq, (part5p_start, part5p_stop), (part3p_start, part3p_stop))
        self._addStemToComponentArray(newStem)
        self.addStem(stemLabel, newStem)



    '''
    Function Name: _parseHairpinData()
    Description: Internal method used by _loadFile() to parse all the hairpin information in the
    structure type file into the Structureobject
    Parameters:
            (hairpinData) - list - list of data from a line in the structure type file describing a hairpin
    Return Type:
            None
    '''
    def _parseHairpinData(self, hairpinData):
        #get hairpin label
        hairpinLabel = hairpinData[0]

        #get index of start of hairpin
        hairpin_start = ''
        for char in hairpinData[1]:
            if char.isnumeric():
                hairpin_start += char
            else:
                break
        hairpin_start = int(hairpin_start)

        #get index of end of hairpin
        hairpin_stop = ''
        for char in reversed(hairpinData[1]):
            if char.isnumeric():
                hairpin_stop += char
            else:
                break
        hairpin_stop = int(hairpin_stop[::-1])

        #get sequence of hairpin
        hairpin_seq = ''
        for char in hairpinData[2]:
            if char.isalpha():
                hairpin_seq += char

        #get index of 5' closing base
        close_5_prime_index = ''
        for char in hairpinData[3]:
            if char.isnumeric():
                close_5_prime_index += char
            elif char == ',':
                break
        close_5_prime_index = int(close_5_prime_index)

        #get index of 3' closing base
        close_3_prime_index = ''
        for char in reversed(hairpinData[3]):
            if char.isnumeric():
                close_3_prime_index += char
            elif char == ',':
                break
        close_3_prime_index = int(close_3_prime_index[::-1])

        close_5_prime_base = hairpinData[4][0] #get 5' closing base
        close_3_prime_base = hairpinData[4][2] #get 3' closing base

        #get pk info
        if hairpinData[5] != '':
            pk = hairpinData[5][3]
        else:
            pk = None


        newHairpin = Hairpin(hairpinLabel, hairpin_seq, (hairpin_start, hairpin_stop),
            (close_5_prime_base, close_3_prime_base), (close_5_prime_index, close_3_prime_index), pk)
        self._addHairpinToComponentArray(newHairpin)
        self.addHairpin(hairpinLabel, newHairpin)



    '''
    Function Name: _parseBulgeData()
    Description: Internal method used by _loadFile() to parse all the bulge information in the
    structure type file into the Structureobject
    Parameters:
            (bulgeData) - list - list of data from a line in the structure type file describing a bulge
    Return Type:
            None
    '''
    def _parseBulgeData(self, bulgeData):
        #get bulge label
        bulgeLabel = bulgeData[0]

        #get start index of bulge
        bulge_start = ''
        for char in bulgeData[1]:
            if char.isnumeric():
                bulge_start += char
            else:
                break
        bulge_start = int(bulge_start)

        #get stop index of bulge
        bulge_stop = ''
        for char in reversed(bulgeData[1]):
            if char.isnumeric():
                bulge_stop += char
            else:
                break
        bulge_stop = int(bulge_stop[::-1])

        #get bulge sequence
        bulge_seq = ''
        for char in bulgeData[2]:
            if char.isalpha():
                bulge_seq += char

        #get index of 5' base in preceding pair
        precedingPair5pIndex = ''
        for char in bulgeData[3]:
            if char.isnumeric():
                precedingPair5pIndex += char
            elif char == ',':
                break
        precedingPair5pIndex = int(precedingPair5pIndex)

        #get index of 3' base in preceding pair
        precedingPair3pIndex = ''
        for char in reversed(bulgeData[3]):
            if char.isnumeric():
                precedingPair3pIndex += char
            elif char == ',':
                break
        precedingPair3pIndex = int(precedingPair3pIndex[::-1])

        #get 5' base of preceding pair
        precedingPair5pBase = bulgeData[4][0]

        #get 3' base in preceding pair
        precedingPair3pBase = bulgeData[4][2]

        #get index of 5' base in trailing pair
        trailingPair5pIndex = ''
        for char in bulgeData[5]:
            if char.isnumeric():
                trailingPair5pIndex += char
            elif char == ',':
                break
        trailingPair5pIndex = int(trailingPair5pIndex)

        #get index of 3' base in trailing pair
        trailingPair3pIndex = ''
        for char in reversed(bulgeData[5]):
            if char.isnumeric():
                trailingPair3pIndex += char
            elif char == ',':
                break
        trailingPair3pIndex = int(trailingPair3pIndex[::-1])

        #get 5' and 3' bases of trailing pair
        trailingPair5pBase = bulgeData[6][0]
        trailingPair3pBase = bulgeData[6][2]

        #need to make sure trailing pair is ordered correctly
        if(abs(bulge_stop - trailingPair5pIndex) > abs(bulge_stop - trailingPair3pIndex)):
            trailingBasePair = (trailingPair3pBase, trailingPair5pBase)
            trailingBasePairIndex = (trailingPair3pIndex, trailingPair5pBase)
        else:
            trailingBasePair = (trailingPair5pBase, trailingPair3pBase)
            trailingBasePairIndex = (trailingPair5pIndex, trailingPair3pIndex)

        #get pk info
        if bulgeData[7] != '':
            pk = bulgeData[7][3]
        else:
            pk = None


        newBulge = Bulge(bulgeLabel, bulge_seq, (bulge_start, bulge_stop), (precedingPair5pBase, precedingPair3pBase),
                        (precedingPair5pIndex, precedingPair3pIndex), trailingBasePair,
                        trailingBasePairIndex, pk)
        self._addBulgeToComponentArray(newBulge)
        self.addBulge(bulgeLabel, newBulge)




    '''
    Function Name: _parseInternalLoopData()
    Description: Internal method used by _loadFile() to parse all the inner loop information in the
    structure type file into the Structureobject
    Parameters:
            (loop1) - list - list of data from a line in the structure type file describing the 5' loop of an Internal Loop
            (loop2) - list - list of data from a line in the structure type file describing the 3' loop of an Internal Loop
    Return Type:
            None
    '''
    def _parseInternalLoopData(self, loop1, loop2):
        #get Inner Loop parent Label
        parentLabel = ''
        for char in loop1[0]:
            if char == '.':
                break
            else:
                parentLabel += char

        #get inner loop subunit label
        loop1SubunitLabel = loop1[0][-1]

        loop2SubunitLabel = loop2[0][-1]

        #get start index of loop subunit
        loop1StartIndex = ''
        for char in loop1[1]:
            if char.isnumeric():
                loop1StartIndex += char
            else:
                break
        loop1StartIndex = int(loop1StartIndex)

        #get stop index of loop subunit
        loop1StopIndex = ''
        for char in reversed(loop1[1]):
            if char.isnumeric():
                loop1StopIndex += char
            else:
                break
        loop1StopIndex = int(loop1StopIndex[::-1])

        loop1Seq = ''
        for char in loop1[2]:
            if char.isalpha():
                loop1Seq += char

        #get start index of loop subunit
        loop2StartIndex = ''
        for char in loop2[1]:
            if char.isnumeric():
                loop2StartIndex += char
            else:
                break
        loop2StartIndex = int(loop2StartIndex)

        #get stop index of loop subunit
        loop2StopIndex = ''
        for char in reversed(loop2[1]):
            if char.isnumeric():
                loop2StopIndex += char
            else:
                break
        loop2StopIndex = int(loop2StopIndex[::-1])

        loop2Seq = ''
        for char in loop2[2]:
            if char.isalpha():
                loop2Seq += char

        #get index of first base in closing pair
        loop1ClosingPairStart = ''
        for char in loop1[3]:
            if char.isnumeric():
                loop1ClosingPairStart += char
            elif char == ',':
                break
        loop1ClosingPairStart = int(loop1ClosingPairStart)

        #get index of second base in closing pair
        loop1ClosingPairEnd = ''
        for char in reversed(loop1[3]):
            if char.isnumeric():
                loop1ClosingPairEnd += char
            elif char == ',':
                break
        loop1ClosingPairEnd = int(loop1ClosingPairEnd[::-1])

        #get index of first base in closing pair
        loop2ClosingPairStart = ''
        for char in loop2[3]:
            if char.isnumeric():
                loop2ClosingPairStart += char
            elif char == ',':
                break
        loop2ClosingPairStart = int(loop2ClosingPairStart)

        #get index of second base in closing pair
        loop2ClosingPairEnd = ''
        for char in reversed(loop2[3]):
            if char.isnumeric():
                loop2ClosingPairEnd += char
            elif char == ',':
                break
        loop2ClosingPairEnd = int(loop2ClosingPairEnd[::-1])

        #store closing pair as a tuple
        closingPairs = ((loop1[4][0], loop1[4][2]), (loop2[4][2], loop2[4][0]))

        newInternalLoop = InternalLoop(parentLabel, loop1SubunitLabel, loop2SubunitLabel, loop1Seq, loop2Seq, (loop1StartIndex, loop1StopIndex),
                                (loop2StartIndex, loop2StopIndex), closingPairs,((loop1ClosingPairStart, loop1ClosingPairEnd), (loop2ClosingPairEnd, loop2ClosingPairStart)))
        self._addInternalLoopToComponentArray(newInternalLoop)
        self.addInternalLoop(parentLabel, newInternalLoop)



    '''
    Function Name: _parseExternalLoopData()
    Description: Internal method used by _loadFile() to parse all the external loop information in the
    structure type file into the Structureobject
    Parameters:
            (externalLoopData) - list - list of data from a line in the structure type file describing an external loop
    Return Type:
            None
    '''
    def _parseExternalLoopData(self, externalLoopData):
        #get label for external loop
        externalLoopLabel = externalLoopData[0]

        #get start index of loop subunit
        startIndex = ''
        for char in externalLoopData[1]:
            if char.isnumeric():
                startIndex += char
            else:
                break
        startIndex = int(startIndex)

        #get stop index of loop subunit
        stopIndex = ''
        for char in reversed(externalLoopData[1]):
            if char.isnumeric():
                stopIndex += char
            else:
                break
        stopIndex = int(stopIndex[::-1])

        #get multiloop sequence
        seq = ''
        for char in externalLoopData[2]:
            if char.isalpha():
                seq += char

        #get index of first base in 5' closing pair
        closingPair5pStart = ''
        for char in externalLoopData[3]:
            if char.isnumeric():
                closingPair5pStart += char
            elif char == ',':
                break
        closingPair5pStart = int(closingPair5pStart)

        #get index of second base in 5' closing pair
        closingPair5pEnd = ''
        for char in reversed(externalLoopData[3]):
            if char.isnumeric():
                closingPair5pEnd += char
            elif char == ',':
                break
        closingPair5pEnd = int(closingPair5pEnd[::-1])

        #store closing pair as a tuple
        closingPair5p = (externalLoopData[4][0], externalLoopData[4][2])

        #get index of first base in 3' closing pair
        closingPair3pStart = ''
        for char in externalLoopData[5]:
            if char.isnumeric():
                closingPair3pStart += char
            elif char == ',':
                break
        closingPair3pStart = int(closingPair3pStart)

        #get index of second base in 3' closing pair
        closingPair3pEnd = ''
        for char in reversed(externalLoopData[5]):
            if char.isnumeric():
                closingPair3pEnd += char
            elif char == ',':
                break
        closingPair3pEnd = int(closingPair3pEnd[::-1])

        #store closing pair as a tuple
        closingPair3p = (externalLoopData[6][0], externalLoopData[6][2])

        newExternalLoop = ExternalLoop(externalLoopLabel, seq, (startIndex, stopIndex), closingPair5p,
                                     (closingPair5pStart, closingPair5pEnd), closingPair3p, (closingPair3pStart, closingPair3pEnd))
        self._addExternalLoopToComponentArray(newExternalLoop)
        self.addExternalLoop(externalLoopLabel, newExternalLoop)


    '''
    Function Name: _getMultiloopParentLabel(self, multiloopString)
    Description:
    parameters:
    Return value
    '''
    def _getMultiloopParentLabel(self, multiloopString):
        parentLabel = ''
        for char in multiloopString:
            if char == '.':
                break
            parentLabel += char

        return parentLabel


    '''
    Function Name: _parseMultiLoopData()
    Description: Internal method used by _loadFile() to parse all the multiloop information in the
    structure type file into the Structureobject
    Parameters:
            (multiloopData) - list - list of data from a line in the structure type file describing an MultiLoop
    Return Type:
            None
    '''
    def _parseMultiLoopData(self, multiloopComponents):
        #get parent label for all multiloop components
        parentLabel = self._getMultiloopParentLabel(multiloopComponents[0][0])

        subunitLabels = [] #list to store all subunit labels
        sequences = {} #dictionary to store all sequences
        spans = {} #dictionary to store all sequence spans
        closingPairs = {} #dictionary to store closing pairs for multiloop segments
        closingPairsSpan = {} #dictionary to store closing pairs spans for multiloop segments

        for multiloopData in multiloopComponents: #iterate through subcomponents

            #get inner loop subunit label and append to subunits list
            subunitLabel = multiloopData[0][-1]
            subunitLabels.append(subunitLabel)

            #get start index of loop subunit
            startIndex = ''
            for char in multiloopData[1]:
                if char.isnumeric():
                    startIndex += char
                else:
                    break
            startIndex = int(startIndex)

            #get stop index of loop subunit
            stopIndex = ''
            for char in reversed(multiloopData[1]):
                if char.isnumeric():
                    stopIndex += char
                else:
                    break
            stopIndex = int(stopIndex[::-1])

            #add start and stop as tuple to spans dictionary
            spans[subunitLabel] = (startIndex, stopIndex)

            #get multiloop sequence and add to sequences dictionary
            seq = ''
            for char in multiloopData[2]:
                if char.isalpha():
                    seq += char

            sequences[subunitLabel] = seq #add sequence to dictionary

            #get index of first base in 5' closing pair
            closingPair5pStart = ''
            for char in multiloopData[3]:
                if char.isnumeric():
                    closingPair5pStart += char
                elif char == ',':
                    break
            closingPair5pStart = int(closingPair5pStart)

            #get index of second base in 5' closing pair
            closingPair5pEnd = ''
            for char in reversed(multiloopData[3]):
                if char.isnumeric():
                    closingPair5pEnd += char
                elif char == ',':
                    break
            closingPair5pEnd = int(closingPair5pEnd[::-1])

            #get whole 5' closing pair
            closingPair5pSpan = (closingPair5pStart, closingPair5pEnd)

            #store closing pair as a tuple
            closingPair5p = (multiloopData[4][0], multiloopData[4][2])

            #get index of first base in 3' closing pair
            closingPair3pStart = ''
            for char in multiloopData[5]:
                if char.isnumeric():
                    closingPair3pStart += char
                elif char == ',':
                    break
            closingPair3pStart = int(closingPair3pStart)

            #get index of second base in 3' closing pair
            closingPair3pEnd = ''
            for char in reversed(multiloopData[5]):
                if char.isnumeric():
                    closingPair3pEnd += char
                elif char == ',':
                    break
            closingPair3pEnd = int(closingPair3pEnd[::-1])

            #store whole 3' closing pair span
            closingPair3pSpan = (closingPair3pStart, closingPair3pEnd)

            #store closing pair as a tuple
            closingPair3p = (multiloopData[6][0], multiloopData[6][2])

            #add closing pairs and closing pair spans to dictionaries
            closingPairs[subunitLabel] = (closingPair5p, closingPair3p)
            closingPairsSpan[subunitLabel] = (closingPair5pSpan, closingPair3pSpan)

        #create new multiloop, add to Structure object, and update component array
        newMultiLoop = MultiLoop(parentLabel, subunitLabels, sequences, spans, closingPairs, closingPairsSpan)
        self._addMultiLoopToComponentArray(newMultiLoop)
        self.addMultiLoop(parentLabel, newMultiLoop)


    '''
    Function Name: _parseNCBPData()
    Description: Internal method used by _loadFile() to parse all the NCBP information in the
    structure type file into the Structureobject
    Parameters:
            (ncbpData) - list - list of data from a line in the structure type file describing an NCBP
    Return Type:
            None
    '''
    def _parseNCBPData(self, ncbpData):
        #get label for ncbp
        ncbpLabel = ncbpData[0]

        #get index of 5' base
        base1Span = int(ncbpData[1])

        #get 5' base
        base1 = ncbpData[2]

        #get index 3' base
        base2Span = int(ncbpData[3])

        #get 3' base
        base2 = ncbpData[4]

        #get location of NCBP in other secondary structure
        if ncbpData[5] == '':
            loc = None
        else:
            loc = ncbpData[5]


        newNCBP = NCBP(ncbpLabel, (base1, base2), (base1Span, base2Span), loc)
        self.addNCBP(ncbpLabel, newNCBP)




    '''
    Function Name: _parseEndData()
    Description: Internal method used by _loadFile() to parse all the End information in the
    structure type file into the Structureobject
    Parameters:
            (endData) - list - list of data from a line in the structure type file describing an end
    Return Type:
            None
    '''
    def _parseEndData(self, endData):
        #get label for End
        endLabel = endData[0]

        #get start index
        startIndex = ''
        for char in endData[1]:
            if char.isnumeric():
                startIndex += char
            else:
                break
        startIndex = int(startIndex)

        #get stop index
        stopIndex = ''
        for char in reversed(endData[1]):
            if char.isnumeric():
                stopIndex += char
            else:
                break
        stopIndex = int(stopIndex[::-1])

        #get sequence
        seq = ''
        for char in endData[2]:
            if char.isalpha():
                seq += char


        newEnd = End(endLabel, seq, (startIndex, stopIndex))
        self._addEndToComponentArray(newEnd)
        self.addEnd(endLabel, newEnd)




    '''
    Function Name: _parsePseudoknotData()
    Description:
    Parameters:
    Return Type:
    '''
    def _parsePsuedoknotData(self, pkData):
        pass



    '''
    Function Name: _parseSegmentData
    Description:
    Parameters:
    Return Type:
    '''
    def _parseSegmentData(self, segData):
        pass


##############################################
###### Add StructureComponent Neighbors ######
##############################################

    '''
    Function Name: _addStructureComponentNeighbors()
    Description: Function fills in the neighboring structure information for each of the StructureComponent objects contained in the Structure object
    Parameters:
            None
    Return Value:
            None
    '''
    def _addStructureComponentNeighbors(self):
        allStructureComponentLabels = self.features() #get all labels for StructureComponents in Structure

        for feature in allStructureComponentLabels: #iterate through all labels
            try:
                structureComponent = self.component(feature) #get the StructureComponent object for a given label
                neighbors = self.neighbors(feature) #get neighbors of the feature
                structureComponent._addNeighbors(neighbors[0], neighbors[1]) #add neighbors to StructureComponent object
            except: #will skip NCBPs and MultiLoops
                continue


#############################
###### COMPONENT ARRAY ######
#############################

    '''
    Function Name: _addStemToComponentArray(stem)
    Description: Internal method used in _loadFile() that adds a given stem to the component array
    Parameters:
            (stem) - Stem object - stem to be added to the component array
    Return Type:
            None
    '''
    def _addStemToComponentArray(self, stem):
        for i in range(stem.sequence5pSpan()[0]-1, stem.sequence5pSpan()[1]):
            self._componentArray[i] = stem.label()

        for i in range(stem.sequence3pSpan()[0]-1, stem.sequence3pSpan()[1]):
            self._componentArray[i] = stem.label()

    '''
    Function Name: _addBulgeToComponentArray(bulge)
    Description:  Internal method used in _loadFile() that adds a given bulge to the component array
    Parameters:
            (bulge) - Bulge object - bulge to be added to the component array
    Return Type:
            None
    '''
    def _addBulgeToComponentArray(self, bulge):
        for i in range(bulge.span()[0]-1, bulge.span()[1]):
            self._componentArray[i] = bulge.label()

    '''
    Function Name: _addHairpinToComponentArray(hairpin)
    Description: Internal method used in _loadFile() that adds a hairpin to the component array
    Parameters:
            (hairpin) - Hairpin object - hairpin to be added to the component array
    Return Type:
            None
    '''
    def _addHairpinToComponentArray(self, hairpin):
        for i in range(hairpin.span()[0]-1, hairpin.span()[1]):
            self._componentArray[i] = hairpin.label()

    '''
    Function Name: _addEndToComponentArray(end)
    Description: Internal method used in _loadFile() that adds an end to the component array
    Parameters:
            (end) - End object - end to be added to the component array
    Return Type:
            None
    '''
    def _addEndToComponentArray(self, end):
        for i in range(end.span()[0]-1, end.span()[1],):
            self._componentArray[i] = end.label()

    '''
    Function Name: _addInternalLoopToComponentArray(InternalLoop)
    Description: Internal method used in _loadFile() that adds an inner loop to the component array
    Parameters:
            (internalLoop) - InternalLoop object - inner loop to be added to the component array
    Return Type:
            None
    '''
    def _addInternalLoopToComponentArray(self, internalLoop):
        for pair in internalLoop.span():
            for i in range(pair[0]-1, pair[1]):
                self._componentArray[i] = internalLoop.label()

    '''
    Function Name: _addExternalLoopToComponentArray(el)
    Description: Internal method used in _loadFile() that adds an external loop to the component array
    Parameters:
            (el) - ExternalLoop object - External loop to be added to the component array
    Return Type:
            None
    '''
    def _addExternalLoopToComponentArray(self, el):
        for i in range(el.span()[0]-1, el.span()[1]):
            self._componentArray[i] = el.label()

    '''
    Function Name: _addMultiLoopToComponentArray(multiloop)
    Description: Internal method used in _loadFile() that adds an multiloop to the component array
    Parameters:
            (multiloop) - MultiLoop object - multiloop to be added to the component array
    Return Type:
             None
    '''
    def _addMultiLoopToComponentArray(self, multiloop):
        for subunit in multiloop._subunitLabels: #iterate through subunit labels
            span = multiloop._spans[subunit] #get span for particular subunit
            for i in range(span[0]-1, span[1]):
                self._componentArray[i] = multiloop._parentLabel

    '''
    Function Name: componentArray()
    Description: function that returns the _componentArray for the Structureobject
    Parameters:
            None
    Return Type:
            numpy array
    '''
    def componentArray(self):
        return self._componentArray






###########################
###### SEQUENCE INFO ######
###########################



    '''
    Function Name: Name()
    Description: Function returns that name of the RNA molecule represented in the .st file
    Parameters:
            None
    Return Type:
            str
    '''
    def name(self):
        return self._name

    '''
    Function Name: Length()
    Description: Function returns the length of the RNA sequence represented in the .st file
    Parameters:
            None
    Return Type:
            int
    '''
    def length(self):
        return self._length

    '''
    Function Name: PageNum()
    Description: Function returns the page number for the RNA molecule represented in the .st file
    Parameters:
            None
    Return Type:
            int
    '''
    def pageNum(self):
        return self._pageNum


########################################
###### STRUCTURAL REPRESENTATIONS ######
########################################


    '''
    Function Name: Sequence()
    Description: Function returns the RNA sequence(A,U,G,C) for the RNA molecule represented in the .st file
    Parameters:
            None
    Return Type:
            str
    '''
    def sequence(self):
        return self._sequence

    '''
    Function Name: DotBracket()
    Description: function to get the Dot Bracket notation for the Structureobject
    Parameters:
            None
    Return Type:
            str
    '''
    def dotBracket(self):
        return self._DBN


    '''
    Function Name: StructureArray()
    Description: Function to get the structure Array for the StructureTyoe object
    Parameters:
            None
    Return Type:
            str
    '''
    def structureArray(self):
        return self._structureArray


    '''
    Function Name: VARNA()
    Description:
    Parameters:
            None
    Return Type:
             str
    '''
    def VARNA(self):
        return self._varna



###################
###### STEMS ######
###################


    '''
    Function Name: _addStemBulgeNeighborBooleans()
    Description: Function fills in boolean values for whether or not each stem object is adjacent to a bulge of length 1
    Parameters:
            None
    Return value:
            None
    '''
    def _addStemBulgeNeighborBooleans(self):
        for stem in self.stems(): #iterate through stems
            neighbor5p, neighbor3p = self.neighbors(stem.label(), object=True) #get neighbors
            bool5p, bool3p = False, False

            #check if a 5' neighbor is a length=1 bulge and if so change bool5p to True
            if(neighbor5p[0] != 'EOM' and neighbor5p[0].label()[0] == 'B'):
                if (neighbor5p[0].sequenceLen() == 1):
                    bool5p = True
            elif(neighbor5p[1] != 'EOM' and neighbor5p[1].label()[0] == 'B'):
                if (neighbor5p[1].sequenceLen() == 1):
                    bool5p = True

            #check if a 3' neighbor is a length=1 bulge and if so change bool3p to True
            if(neighbor3p[0] != 'EOM' and neighbor3p[0].label()[0] == 'B'):
                if (neighbor3p[0].sequenceLen() == 1):
                    bool3p = True
            elif(neighbor3p[1] != 'EOM' and neighbor3p[1].label()[0] == 'B'):
                if (neighbor3p[1].sequenceLen() == 1):
                    bool3p = True

            stem._addAdjacentBulgeBoolean(bool5p, bool3p)


    '''
    Function Name: addStem()
    Description: Function to add a new stem to the Structureobject
    Parameters:
            (stemLabel) - str - key for stem object in self._stems dictionary
            (newStem) - Stem object - Stem object to be stored at given key in the self._stems dictionary
    Return Value:
            None
    '''
    def addStem(self, stemLabel, newStem):
        self._stems[stemLabel] = newStem


    '''
    Function Name: stemLabels()
    Description: function to access all the stem labels for the Structureobject
    Parameters:
            None
    Return Type:
            list
    '''
    def stemLabels(self):
        return list(self._stems.keys())


    '''
    Function Name: stems(label=None)
    Description: function to get all the stem objects in the Structureobject. If label is provided the stem object
    with the matching label is returned.
    Parameters:
            (label=None) - str - the label for the stem being searched for.
    Return Type:
            list or Stem object
    '''
    def stems(self, label=None):
        if label:
            return self._getStemByLabel(label)
        else:
            return list(self._stems.values())


    '''
    Function Name: numStems()
    Description: function to get the number of stems in a Structureobject
    Parameters:
             None
    Return Value:
             int
    '''
    def numStems(self):
        return len(self._stems)


    '''
    Function Name: getStemByLabel(stemLabel)
    Description: Function to get a particular Stem object based on its label
    Parameters:
            (stemLabel) - str - label for Stem to be accessed
    Return Value:
             Stem object
    '''
    def _getStemByLabel(self, stemLabel):
        try:
            stem = self._stems[stemLabel]
            return stem
        except KeyError:
            print(f'Stem label: {stemLabel} not found.')
            return None



######################
###### HAIRPINS ######
######################



    '''
    Function Name: addHairpin(label, newHairpin)
    Description: Function to add a new hairpin to the Structureobject
    Parameters:
            (label) - str - key for Hairpin object in self._hairpins dictionary
            (newHairpin) - Hairpin object - Hairpin object to be stored at the given key in the self._haripins dicitonary
    Return Type:
            None
    '''
    def addHairpin(self, label, newHairpin):
        self._hairpins[label] = newHairpin


    '''
    Function Name: hairpinLabels()
    Description: Function to access all the hairpin labels in the Structureobject
    Parameters:
            None
    Return Type:
             list
    '''
    def hairpinLabels(self):
        return list(self._hairpins.keys())


    '''
    Function Name: hairpins()
    Description: Function to get all Hairpin objects in the Structureobject. if label is provided, the Hairpin object with the
    matching label is returned.
    Parameters:
            None
    Return Type:
             list
    '''
    def hairpins(self, label=None):
        if label:
            return self._getHairpinByLabel(label)
        else:
            return list(self._hairpins.values())


    '''
    Function Name: numHairpins()
    Description: Function to get the number of hairpins in the Structureobject
    Parameters:
            None
    Return Type:
             int
    '''
    def numHairpins(self):
        return len(self._hairpins)

    '''
    Function Name: getHairpinByLabel(hairpinLabel)
    Description: Function to access a particular Hairpin Object based on its label
    Parameters:
            (hairpinLabel) - str - label for hairpin to accessed
    Return Type:
             Hairpin Object
    '''
    def _getHairpinByLabel(self, hairpinLabel):
        try:
            hairpin = self._hairpins[hairpinLabel]
            return hairpin
        except KeyError:
            print(f'Hairpin label: {hairpinLabel} not found')



####################
###### BULGES ######
####################



    '''
    Function Name: addBulge(bulgeLabel, newBulge)
    Description: Function to add a new Bulge object to the Structureobject
    Parameters:
            (bulgeLabel) - str - the key value to be used for the new Bulge object
            (newBulge) - Bulge Object - Bulge object to be stored at the given key in the self._bulges dictionary
    Return Type:
             None
    '''
    def addBulge(self, bulgeLabel, newBulge):
        self._bulges[bulgeLabel] = newBulge


    '''
    Function Name: bulgeLabels()
    Description: Function to get all the bulge labels for the Structureobject
    Parameters:
            None
    Return Type:
             list
    '''
    def bulgeLabels(self):
        return list(self._bulges.keys())

    '''
    Function Name: bulges()
    Description: Function to get all the Bulge objects for the Structureobject
    Parameters:
            None
    Return Type:
             list
    '''
    def bulges(self, label=None):
        if label:
            return self._getBulgeByLabel(label)
        else:
            return list(self._bulges.values())

    '''
    Function Name: numBulges()
    Description: Function to get the number of bulges in a given Structureobject
    Parameters:
            None
    Return Type:
             int
    '''
    def numBulges(self):
        return len(self._bulges)

    '''
    Function Name: getBulgeByLabel(bulgeLabel)
    Description: Function to access a particular Bulge object based on its label
    Parameters:
            (bulgeLabel) - str - label for Bulge object to be accessed
    Return Type:
             Bulge object
    '''
    def _getBulgeByLabel(self, bulgeLabel):
        try:
            bulge = self._bulges[bulgeLabel]
            return bulge
        except KeyError:
            print(f'Bulge label: {bulgeLabel} not found')
            return None


###########################
###### Inner Loops ########
###########################



    '''
    Function Name: addInternalLoop(parentLabel, subunitLabel, newInternalLoop)
    Description: function to add a new InternalLoop object to the Structureobject
    Parameters:
            (parentLabel) - str - parent key value for the Inner loop to be added
            (newInternalLoop) - InternalLoop object - InternalLoop object to be stored at the the given key in the self._internalLoops dictionary
    Return Type:
            None
    '''
    def addInternalLoop(self, parentLabel, newInternalLoop):
        self._internalLoops[parentLabel] = newInternalLoop


    '''
    Function Name: InternalLoopLabels()
    Description: Function to return a list of all the inner loop labels in the Structure object
    Parameters:
            None
    Return Type:
             list
    '''
    def internalLoopLabels(self):
        return list(self._internalLoops.keys())


    '''
    Function Name: InternalLoops()
    Description: Function to return a list of the InternalLoop objects in the Structure object
    Parameters:
            None
    Return Type:
            list
    '''
    def internalLoops(self, label=None):
        if label:
            return self._getInternalLoopByLabel(label)
        else:
            return list(self._internalLoops.values())


    '''
    Function Name: numInternalLoops()
    Description: function to get the number of inner loops in a Structure object
    Parameters:
            None
    Return Type:
             int
    '''
    def numInternalLoops(self):
        return len(self._internalLoops)


    '''
    Function Name: getInternalLoopByLabel(label)
    Description: returns InternalLoop object stored at the given key value
    Parameters:
            (label) - str - the key value for the inner loop to be accessed
    Return Type:
             InternalLoop object
    '''
    def _getInternalLoopByLabel(self, label):
        try:
            internalLoop = self._internalLoops[label]
            return internalLoop
        except KeyError:
            print(f'Internal Loop: {label} not found.')
            return None


    '''
    Function Name: getInternalLoopSubunitByLabel(parentLabel, subunitLabel)
    Description:
    Parameters:
             (parentLabel) - str - parent label for inner loop object
            (subunitLabel) - str - subunit label for the inner loop object
    Return Type:
            dictionary with the following subunit information

                {
                    'label' : ,
                    'Sequence' : ,
                    'span' :
                }
    '''
    def getInternalLoopSubunit(self, parentLabel, subunitLabel):
        internalLoop = None
        try:
            internalLoop = self._internalLoops[parentLabel]
        except KeyError:
            print(f'Inner Loop: {parentLabel} not found.')


        if subunitLabel == '1':
            subunit = {
                'label' : f'{internalLoop._parentLabel}.1',
                'Sequence' : internalLoop._5pSequence,
                'span' : internalLoop._span5p,
            }
            return subunit
        elif subunitLabel == '2':
            subunit = {
                'label' : f'{internalLoop._parentLabel}.2',
                'Sequence' : internalLoop._3pSequence,
                'span' : internalLoop._span5p,
            }
            return subunit
        else:
            return None





###########################
###### MultiLoops #########
###########################


    '''
    Function Name: addMultiLoop(parentLabel, subunitLabel, newMultiLoop)
    Description: function to add a new MultiLoop object to the Structure object
    Parameters:
            (parentLabel) - str - parent multiloop label for the multiloop to be added
            (subunitLabel) - str - subunit label for the multiloop to be added
            (newMultiLoop) - MultiLoop object - MultiLoop object to be added at the given key values
    Return Type:
             None
    '''
    def addMultiLoop(self, parentLabel, newMultiLoop):
        self._multiLoops[parentLabel] = newMultiLoop

    '''
    Function Name: numMultiLoops()
    Description: Function to get the number of multiloops in a StructureTyoe object
    Parameters:
            None
    Return Type:
             int
    '''
    def numMultiLoops(self):
        return len(self._multiLoops)

    '''
    Function Name: _getMultiLoopByLabel(self, label)
    Description: Function returns the MultiLoop object identified by a given parent label.
    Parameters:
            (label) -- str -- parent label of the MultiLoop object being accessed
    Return Value:
            MultiLoop object
    '''
    def _getMultiLoopByLabel(self, label):
        try:
            multiloop = self._multiLoops[label]
            return multiloop
        except KeyError:
            print(f'MultiLoop {label} not found.')
            return None

    '''
    Function Name: multiLoops(label=None)
    Description: Function returns a list of all the multiloops stored in a Structure object. If a particular label is provided, only the MultiLoop object with that label will be returned
    Parameters:
            (label=None) -- str -- parent label of the MultiLoop object being accessed. Default value is None
    Return Value:
            if label is provided: MultiLoop object
            if label=None: list of MultiLoop objects
    '''
    def multiLoops(self, label=None):
        if(label):
            return self._getMultiLoopByLabel(label)
        else:
            return list(self._multiLoops.values())



###############################
###### External Loops #########
###############################


    '''
    Function Name: addExternalLoop(elLabel, newEL)
    Description: function to add a new External Loop object to the Structure object
    Parameters:
            (elLabel) - str - the key value to be used for the new ExternalLoop object
            (newEL) - ExternalLoop Object - External loop to be stored at given key value
    Return Type:
             None
    '''
    def addExternalLoop(self, elLabel, newEL):
        self._externalLoops[elLabel] = newEL

    '''
    Function Name: externalLoopLabels()
    Description: Function to return a list of the external loop labels for a Structure object
    Parameters:
            none
    Return Type:
             list
    '''
    def externalLoopLabels(self):
        return list(self._externalLoops.keys())

    '''
    Function Name: externalLoops()
    Description: function to return a list of all the ExternalLoop Objects in a Structure object
    Parameters:
            None
    Return Type:
             list
    '''
    def externalLoops(self, label=None):
        if label:
            return self._getExternalLoopByLabel(label)
        else:
            return list(self._externalLoops.values())

    '''
    Function Name: numExternalLoops()
    Description: function to return the number of external loops in a Structure object
    Parameters:
            None
    Return Type:
             int
    '''
    def numExternalLoops(self):
        return len(self._externalLoops)

    '''
    Function Name: getExternalLoopByLabel(elLabel)
    Description: Function to access a particular ExternalLoop Object based on its label
    Parameters:
            (elLabel) - str - label for the ExternalLoop to be accessed
    Return Type:
             ExternalLoop object
    '''
    def _getExternalLoopByLabel(self, elLabel):
        try:
            el = self._externalLoops[elLabel]
            return el
        except KeyError:
            print(f'External Loop: {elLabel} not found.')
            return None



####################
###### NCBP ########
####################



    '''
    Function Name: addNCBP(ncbpLabel, newNCBP)
    Description: Function to add a new NCBP to the Structure object
    Parameters:
            (ncbpLabel) - str - label for the new NCBP
            (newNCBP) - NCBP Object - NCBP object to be added
    Return Type:
             None
    '''
    def addNCBP(self, ncbpLabel, newNCBP):
        self._ncbp[ncbpLabel] = newNCBP

    '''
    Function Name: ncbpLabels()
    Description: function to return a list of all the ncbp labels in a given Structure object
    Parameters:
            None
    Return Type:
             list
    '''
    def ncbpLabels(self):
        return list(self._ncbp.keys())

    '''
    Function Name: NCBPs()
    Description: function to return a list of all the NCBPs in a Structure object
    Parameters:
            None
    Return Type:
             list
    '''
    def NCBPs(self, label=None):
        if label:
            return self._getNCBPByLabel(label)
        else:
            return list(self._ncbp.values())

    '''
    Function Name: numNCBPs()
    Description: Function to get the number of NCBP's in a given Structure object
    Parameters:
            None
    Return Type:
             int
    '''
    def numNCBPs(self):
        return len(self._ncbp)

    '''
    Function Name: getNCBPByLabel(ncbpLabel)
    Description: Function to get a particular NCBP object based on its label
    Parameters:
            (ncbpLabel) - str - label for NCBP object to be accessed
    Return Type:
             NCBP object
    '''
    def _getNCBPByLabel(self, ncbpLabel):
        try:
            ncbp = self._ncbp[ncbpLabel]
            return ncbp
        except KeyError:
            print(f'NCBP label: {ncbpLabel} not found.')
            return None



####################
###### ENDS ########
####################



    '''
    Function Name: addEnd(endLabel, newEnd)
    Description: Function to add a new End to the Structure object
    Parameters:
            (endLabel) - str - label for new End object
            (newEnd) - End Object - new End object to be added
    Return Type:
             None
    '''
    def addEnd(self, endLabel, newEnd):
        self._ends[endLabel] = newEnd

    '''
    Function Name: endLabels()
    Description: Function to return a list of all the end labels for the Structure object
    Parameters:
            None
    Return Type:
             list
    '''
    def endLabels(self):
        return list(self._ends.keys())

    '''
    Function Name: ends()
    Description: function to return a list of all the End objects for the Structure object
    Parameters:
            (label=None) - str - label for the end being accessed
    Return Type:
            End object or list
    '''
    def ends(self, label=None):
        if label:
            return self._getEndByLabel(label)
        else:
            return list(self._ends.values())

    '''
    Function Name: numEnds()
    Description: function to get the number of ends in Structure object
    Parameters:
            None
    Return Type:
             int
    '''
    def numEnds(self):
        return len(self._ends)

    '''
    Function Name: getEndByLabel(endLabel)
    Description: function to access a particular End Object based on its label
    Parameters:
            (endLabel) - str - the label for the End object to be accessed
    Return Type:
             End Object
    '''
    def _getEndByLabel(self, endLabel):
        try:
            end = self._ends[endLabel]
            return end
        except KeyError:
            print(f'End: {endLabel} not found.')
            return None



################################
######## OTHER FUNCTIONs #######
################################

    '''
    Function Name: component(label, subLabel=None)
    Description: More general form of get<secondarystructure>ByLabel(). Allows you to access a given
    structure type component based on its label. Useful when using the componentArray because you may not know which
    type of structure type component you will be accessing if you are(for example) looping through the entire array
    Parameters:
            (label) - str - label of the feature to be accessed
    Return Type:
            StructureComponent object
    '''
    def component(self, label):
        if label[0] == 'S':
            return self._getStemByLabel(label)
        elif label[0] == 'H':
            return self._getHairpinByLabel(label)
        elif label[0] == 'B':
            return self._getBulgeByLabel(label)
        elif label[0] == 'X':
            return self._getExternalLoopByLabel(label)
        elif label[0] == 'E':
            return self._getEndByLabel(label)
        elif label[0] == 'N':
            return self._getNCBPByLabel(label)
        elif label[0] == 'I':
            return self._getInternalLoopByLabel(label)
        elif label[0] == 'M':
            return self._getMultiLoopByLabel(label)
        else:
            #if label is not handled by any of these blocks
            print(f'Label: {label} not found in Structure object.')
            return None


    '''
    Function Name: neighbors(label, object=False)
    Description: Function to get the secondary structures adjacent to the feature of interest.
    Parameters:
            (label) - str - label for the feature of interest
            (object) - bool - optional argument that causes the function to return the actual StructureTypeComponent objects instead of just the object label
    Return Type:
            Returns a tuple containing the labels for the adjacent features in order of 5' to 3' locations
    '''
    def neighbors(self, label, object=False):
        adjacentFeatures = [] #list to store the adjacent RNA features

        if label in self._componentArray: #check if the feature is valid
            span = self.component(label).span() #get index locations of the feature
            if all(type(i) is int for i in span): #tuple only containes integer index locations(example: bulge location)
                try: #try/except block will handle ends which only have one neighbor and one out of range index
                    neighbor5p = (self._componentArray[span[0]-2] if not object else self.component(self._componentArray[span[0]-2]))
                except:
                    neighbor5p = 'EOM' # 'End of Molecule'

                try: #try/except block will handle ends which only have one neighbor and one out of range index
                    neighbor3p = (self._componentArray[span[1]] if not object else self.component(self._componentArray[span[1]]))
                except:
                    neighbor3p = 'EOM'

                adjacentFeatures = (neighbor5p, neighbor3p)

            else: #tuple containes other tuples within it(example: InternalLoop locations)
                try: #try/except block will handle ends which only have one neighbor and one out of range index
                    seq1_neighbor5p = (self._componentArray[span[0][0]-2] if not object else self.component(self._componentArray[span[0][0]-2]))
                except:
                    seq1_neighbor5p = 'EOM'

                try: #try/except block will handle ends which only have one neighbor and one out of range index
                    seq1_neighbor3p = (self._componentArray[span[0][1]] if not object else self.component(self._componentArray[span[0][1]]))
                except:
                    seq1_neighbor3p = 'EOM'

                try: #try/except block will handle ends which only have one neighbor and one out of range index
                    seq2_neighbor5p = (self._componentArray[span[1][0]-2] if not object else self.component(self._componentArray[span[1][0]-2]))
                except:
                    seq2_neighbor5p = 'EOM'

                try: #try/except block will handle ends which only have one neighbor and one out of range index
                    seq2_neighbor3p = (self._componentArray[span[1][1]] if not object else self.component(self._componentArray[span[1][1]]))
                except:
                    seq2_neighbor3p = 'EOM'

                adjacentFeatures = ((seq1_neighbor5p, seq2_neighbor3p), (seq2_neighbor5p, seq1_neighbor3p))

            return adjacentFeatures
        else: #otherwise print error and return None
            return None


    """
    Function: features()
    Description: Function to return a list of all the StructureComponent labels in a bpRNAStructure object
    Parameters:
            None
    Return Value:
            list of strings, each string is the label for a StructureComponent in the bpRNAStructure
            *Currently only return the StructureComponents: stems, bulges, hairpins, and inner loops
    """
    def features(self):
        featureLabels = []

        featureLabels.extend(self.stemLabels())
        featureLabels.extend(self.bulgeLabels())
        featureLabels.extend(self.hairpinLabels())
        featureLabels.extend(self.internalLoopLabels())

        return featureLabels
