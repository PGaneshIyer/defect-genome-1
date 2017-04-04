import os, sys
import re
import csv
import tarfile
import numpy as np
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Poscar


class InputOutput:
    """
    TODO:replace default python file io with monty.io
    """

    def poscar_cart_to_direct(self, inputFilePath, outputFilePath,
                              vasp4Compatible=False):
        noDecimals = 6
        direct = True
        crystalStruc = Structure.from_file(inputFilePath)
        poscarFile = Poscar(crystalStruc)
        poscarString = poscarFile.get_string(direct, vasp4Compatible,
                                             noDecimals)
        with open(outputFilePath, "w") as f:
            f.write(poscarString)

    def poscar_direct_to_cart(self, inputFilePath, outputFilePath,
                              vasp4Compatible=False):
        noDecimals = 6
        direct = False
        crystalStruc = Structure.from_file(inputFilePath)
        poscarFile = Poscar(crystalStruc)
        poscarString = poscarFile.get_string(direct, vasp4Compatible,
                                             noDecimals)
        with open(outputFilePath, "w") as f:
            f.write(poscarString)

    # @staticmethod
    # def get_occurence_infile_from_end(filePath, str1, regexSplit=' ', noOccurence=1):
    #     lineList = []
    #     count = 0
    #     with open(filePath, 'r') as f:
    #         for line in reversed(f.readlines()):
    #             if str1 in line and count < noOccurence:
    #                 line_val = re.split(regexSplit, line)
    #                 lineList.append(line_val)
    #                 count += 1
    #             elif str1 in line and count >= noOccurence:
    #                 break
    #     return lineList

    @staticmethod
    def get_occurence_infile(filePath, str1, regexSplit=' ', noOccurence=None,
                             reverse=False):
        lineList = []
        count = 0
        with open(filePath, 'r') as f:
            if reverse:
                for line in reversed(f.readlines()):
                    if str1 in line:
                        if (noOccurence is not None):
                            if (count < noOccurence):
                                line_val = re.split(regexSplit, line)
                                lineList.append(line_val)
                                count += 1
                            else:
                                break
                        else:
                            line_val = re.split(regexSplit, line)
                            lineList.append(line_val)
            else:
                for line in f.readlines():
                    if str1 in line:
                        if (noOccurence is not None):
                            if (count < noOccurence):
                                line_val = re.split(regexSplit, line)
                                lineList.append(line_val)
                                count += 1
                            else:
                                break
                        else:
                            line_val = re.split(regexSplit, line)
                            lineList.append(line_val)
        return lineList

    @staticmethod
    def get_occurence_inlist_after_str(list1, str1, noOccurence=1,
                                       dataType=None):
        occurenceList = []
        listValWithStr = [x for x in list1 if str1 in x]
        indexValWithStr = list1.index(listValWithStr[0])
        count = 0
        for index, val in enumerate(list1):
            if index > indexValWithStr:
                try:
                    if dataType == float:
                        val = float(val)
                    elif dataType == int:
                        val = int(float(val))
                except ValueError:
                    continue
                if (isinstance(val, dataType) and count < noOccurence):
                    occurenceList.append(val)
                    count += 1
                elif (isinstance(val, dataType) and count >= noOccurence):
                    break
        return occurenceList

    @staticmethod
    def list_to_csv(filePath, listVal, headerVal=None,
                    delimiter=','):
        with open(filePath, 'w') as f:
            writer = csv.writer(f, delimiter=delimiter)
            if (headerVal is not None):
                if (len(headerVal) != len(listVal[0])):
                    print 'WARNING: Header length smaller than list'
                writer.writerow(headerVal)
            writer.writerows(listVal)

    @staticmethod
    def csv_to_list(filePath, delimiterVal=',', skipHeader=True,
                    convertFloat=True, csvExcel=False, transpose=False,
                    convertNumpy=False):
        outputList = []
        if csvExcel:
            reader = csv.reader(open(filePath, 'rU'),
                                delimiter=delimiterVal,
                                dialect=csv.excel_tab)
            if skipHeader:
                next(reader, None)
            for row in reader:
                if convertFloat:
                    for count, val in enumerate(row):
                        try:
                            row[count] = float(val)
                        except ValueError:
                            continue
                outputList.append(row)
        else:
            with open(filePath, "r") as f:
                reader = csv.reader(f, delimiter=delimiterVal)
                if skipHeader:
                    next(reader, None)
                for row in reader:
                    if convertFloat:
                        for count, val in enumerate(row):
                            try:
                                row[count] = float(val)
                            except ValueError:
                                continue
                    outputList.append(row)
        if transpose:
            if convertNumpy:
                raise Warning('Transpose is not compatible with convertNumpy')
            else:
                outputList = list(map(list, zip(*outputList)))
        if convertNumpy:
            raise Warning('Convertion to Numpy does not create 2D array. Use with care')
            # Adapted from http://stackoverflow.com/a/32305363/1652217
            # FIXME: ASSuMEs equal number of columns
            # FIXME: ASSuMEs same datatype for all rows
            # FIXME: Only differentiates between string and float
            # FIXME: Generalize through function call for every row
            # FIXME: This screws up the creation of 2D array
            dtypeVal = ''
            for val in outputList[0]:
                if isinstance(val, basestring):
                    dtypeVal += 'S5,'
                else:
                    dtypeVal += 'f,'
            dtypeVal = dtypeVal[:-1]
            outputList = np.array([tuple(x) for x in outputList],
                                  dtype=dtypeVal)
        return outputList

    @staticmethod
    def file_to_list(filePath, splitVal=None, lineStart=0,
                     lineEnd=None, colStart=0, colEnd=None,
                     convertFloat=True):
        outputList = []
        with open(filePath) as f:
            content = f.readlines()
            if lineEnd is None:
                lineEnd = len(content)
            if colEnd is None:
                colEnd = len(content[lineStart])
            for count, val in enumerate(content):
                if (count >= lineStart and count <= lineEnd):
                    if splitVal is None:
                        valList = val.split()[colStart:colEnd]
                    else:
                        print 'WARNING: Only default split value can be used'
                    if convertFloat:
                        for count2, val2 in enumerate(valList):
                            try:
                                valList[count2] = float(val2)
                            except ValueError:
                                continue
                    outputList.append(valList)
        return outputList

    # HAS BUGS
    # def replace_line_with_string(key_to_replace, new_value, filePath,
    #                              backup='.bak'):
    #     # 'import fileinput' is required
    #     # for line in fileinput.input(filePath, inplace=True, backup=backup):
    #     keyPresent = False
    #     f = fileinput.FileInput(filePath, inplace=True, backup='.bak')
    #     for line in f:
    #         if line.find(key_to_replace) >= 0:
    #             print new_value
    #             keyPresent = True
    #         else:
    #             print line[:-1]
    #             f.close()
    #     if (keyPresent is False):
    #         with open(filePath, "a") as f:
    #             f.write(new_value)

    def make_tarfile(dirPath, inFileList, tarFileName, option=None):
        with tarfile.open(os.path.join(dirPath, tarFileName),
                          "w:gz") as tar:
            for inFile in inFileList:
                tar.add(inFile)

