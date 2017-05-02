#  Authors: Ram Balachandran.
#  Center for Nanophase Materials Sciences, Oak Ridge National Laboratory, Oak Ridge, TN.

import os,sys
import math
from operator import itemgetter
import numpy as np
import pymatgen as mg
import copy
import csv
from collections import namedtuple
from enum import Enum
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.core.periodic_table import Element
#VaspIO
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.core import Spin, OrbitalType

import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Ellipse, Polygon, Rectangle
import matplotlib.colors as colors
from matplotlib import cm
from monty.io import zopen
from pprint import pprint
import warnings
# Class inputOutput
import fileinput

# Class VaspIO
import re


class AnalyzeStructure:
    """
    Performs analysis not available in the structure class of pymatgen.
    TODO: Should be merged into analysis.chemenv  of Pymatgen.
    """
    def find_farthest_atoms(self, crystalStruc, ref_atomType, ref_atomIndex,
                           find_atomType, noAtoms, maxDistance=None,
                           verbose=False):
        ref_atom = crystalStruc[ref_atomIndex]
        find_atomList = []
        for index, atom in enumerate(crystalStruc):
            if atom.species_string in find_atomType:
                distance = crystalStruc.get_distance(ref_atomIndex, index)
                if maxDistance:
                    if (distance < maxDistance):
                        find_atomList.append([atom, index, distance])
                else:
                    find_atomList.append([atom, index, distance])
        find_atomList.sort(key=itemgetter(1), reverse=True)
        return find_atomList[0:noAtoms]

    def find_closest_atoms(self, crystalStruc, ref_atomType, ref_atomIndex,
                           find_atomType, noAtoms, maxDistance=None,
                           verbose=False, removeSelf=True):
        ref_atom = crystalStruc[ref_atomIndex]
        find_atomList = []
        for index, atom in enumerate(crystalStruc):
            if atom.species_string == find_atomType:
                distance = crystalStruc.get_distance(ref_atomIndex, index)
                if maxDistance:
                    if (distance < maxDistance):
                        # find_atomList.append([atom, index, distance])
                        if atom.species_string == ref_atomType and removeSelf == True:
                            if distance > 1.0E-03:
                                find_atomList.append([atom, index, distance])
                        else:        
                            find_atomList.append([atom, index, distance])
                else:
                    # find_atomList.append([atom, index, distance])
                    if atom.species_string == ref_atomType and removeSelf == True:
                        if distance > 1.0E-03:
                            find_atomList.append([atom, index, distance])
                    else:        
                        find_atomList.append([atom, index, distance])
                    
        find_atomList.sort(key=itemgetter(2))
        return find_atomList[0:noAtoms]

    def return_index_of_atomtype(self, crystalStruc, atomTypes,
                                 noIndices=None):
        """
        :param crystalStruc: (pymatgen structure object) contains 
         the structure information
        :param atomTypes: [list] a list of strings of atom 
         Types of interest
        :return: dictionary containing atom indices of interest.
        Issues & Comments:
            1. In case of noIndices, break the loop, efficiently
        """
        atomIndices = {}
        for i in range(0, len(atomTypes)):
            atomIndices[atomTypes[i]] = []
            indexCount = 0
            for index, atom in enumerate(crystalStruc):
                if atom.species_string in atomTypes[i]:
                    atomIndices[atomTypes[i]].append(index)
                    indexCount += 1
                    if (isinstance(noIndices, int) and
                        noIndices >= 0 and indexCount >= noIndices):
                        break
        if (isinstance(noIndices, int) and noIndices >= 0):
            return atomIndices[noIndices]
        else:
            return atomIndices

    def is_atom_in_subspace(self, crystalStruc, atomIndex, xRange, yRange,
                            zRange, coordType='frac',jimage=None):
        valRangeList = [xRange, yRange, zRange]
        inRangeList = [False, False, False]
        if coordType=='frac':
            frac_coords = crystalStruc[atomIndex].frac_coords
            latticeVec = [1.0, 1.0, 1.0]
            imageVec = [0, 0, 0]
        else:
            #latticeVec = list(crystalStruc._lattice.abc)
            raise ValueError('Provide range values in Fractional Coordinates')
        for count, coordVal in enumerate(frac_coords):
            coordValImgList = [coordVal-latticeVec[count], coordVal,
                            coordVal+latticeVec[count]]
            inRangeFlag = False
            for imgVal in coordValImgList:
                #print latticeVec, frac_coords, coordVal, coordValImgList, imgVal
                inRangeFlag = Miscellaneous.isValueIn(imgVal, valRangeList[count])
                if inRangeFlag:
                    inRangeList[count] = inRangeFlag
                    imageVec[count] = imgVal- coordVal
                    break
        if all(inRangeList):
            return [True, imageVec]
        else:
            return [False, imageVec]

    """      
    TODO: Needs to be debugged
    def return_index_in_subspace(self, crystalStruc, atomTypes, spaceVal='frac',
                                 xRange=None, yRange=None, zRange=None, noIndices=None,
                                 accuracy = 1.0E-02):
        ## COMMENT BEGINS
        :param crystalStruc: (pymatgen structure object) contains 
         the structure information
        :param atomTypes: [list] a list of strings of atom 
         Types of interest
        :return: dictionary containing atom indices of interest.  
        Issues & Comments:
            1. In case of noIndices, break the loop, efficiently
            2. REQUIRES more bug testing.
        ## COMMENT ENDS
        tmpIndices = self.return_index_of_atomtype(crystalStruc, atomTypes, noIndices)
        atomIndices = {}
        for key,val in tmpIndices.iteritems():
            atomIndexList = []
            for atomIndex in val:
                frac_coords = crystalStruc[atomIndex].frac_coords
                if xRange is None:
                    isValInXrange = True
                else:
                    isValInXrange = Miscellaneous.isValueIn(frac_coords[0],xRange,accuracy=1.0E-03)
                if yRange is None:
                    isValInYrange = True
                else:
                    isValInYrange = Miscellaneous.isValueIn(frac_coords[1],yRange,accuracy=1.0E-03)

                if zRange is None:
                    isValInZrange = True
                else:
                     isValInZrange = Miscellaneous.isValueIn(frac_coords[2],zRange,accuracy=1.0E-03)
                if (isValInXrange and isValInYrange and isValInZrange):
                   atomIndexList.append(atomIndex)
            atomIndices[key] = atomIndexList
        return atomIndices
    """
    
    def get_distance_array(self, crystalStruc, atomType1, atomType2,
                           maxDistance=None):
        """
        Get all possible distances between atomType1 and atomType2
        Args:
            atomType1 (str): string of the element name
            atomType2 (list): list of strings of element names
        Returns:
            distanceArray (list): list of distance values in ascending order
        """
        pass

    def get_neighbor_distance(self, crystalStruc, atomType, atomIndex,
                              neighborType, noNeighbors, distanceArray,
                              dr=1.0E-03):
        """
        Return the distance from an atom we could find noNeighbors 
        values of neighborType atoms
        Args:
            atomIndex (list): list of atomIndex values
            atomType (list): list of strings of element name
            neighborType (list): list of strings of neighbor type atoms
            noNeighbors (list): list of integers of no. of neighbors
            distanceArray (list): list of input distance values   
            dr (float): shell radius 
        Returns:
           a list of list with following information
           atomType, atomIndex, neighborType1, noNeighborsType1, distanceVal
        """
        pass
    
    def change_atom_type(self, crystalStruc, oldAtomType, newAtomType,
                         atomIndex=None, sort=True):
        """
        Change atomType of the index provided from oldAtomType to newAtomType
        Args:
            olAtomType: Name of the atom type to be replaced
            atomIndex (list): list of atomIndex to be replaced
            newAtomType (str): Name of the atom type that will replace
        Returns:
           a list of list with following information
           atomType, atomIndex, neighborType1, noNeighborsType1, distanceVal
        """
        failedIndex = []
        newCrystalStruc = crystalStruc.copy()
        for index in atomIndex:
            try:
                newCrystalStruc[index].species_string in oldAtomType
            except ValueError:
                print 'Index', index, 'FAILED as its not part of atom types',
                oldAtomType
                print 'Moving to next Index'
                failedIndex.append(index)
                continue
            newCrystalStruc.replace(index, newAtomType)
        if sort:
            newCrystalStruc.sort()
        return newCrystalStruc, failedIndex
            
    def bin_atom_by_distance(self, crystalStruc, atomIndices, distanceArray, accuracy=0.01):
        totNoNeighbors = {}
        for i in atomIndices:
            noNeighbors = zerolistmaker(len(distanceArray))
            for j in atomIndices:
                if (i != j):
                    distance = crystalStruc.get_distance(i,j)
                    for index, val in enumerate(distanceArray):
                        if (abs(distance - val) < accuracy):
                            noNeighbors[index] += 1
            totNoNeighbors[i] = noNeighbors
        return totNoNeighbors

    def get_neighbors_to_coords(self, crystalStruc, coords, radius,
                                neighborType, noNeighbors=1, frac_coord=True):
        #print ' get_neighbors_to_coords', coords, radius, frac_coord #TEST
        neighborSites = []
        for count, atom in enumerate(crystalStruc):
            if atom.species_string in neighborType:
                if frac_coord:
                    distance = crystalStruc.lattice.get_distance_and_image(coords,
                                                                           atom.frac_coords)
                else:
                   raise NotImplementedError('Currently can only handle Fractional Coordinates')
                #print 'get_neighbors_to_coords - 2', count, atom.frac_coords, distance  #TEST
                if (distance[0] <= radius):
                    neighborSites.append((atom, distance[0], count))
        neighborSites.sort(key=lambda val: val[1])
        #print pprint(neighborSites) #TEST
        return neighborSites[0:noNeighbors]

    def calc_displacement(self, crystalStrucInit, crystalStrucFinal,
                          maxDisplacement, indexInit=None, atomType=None, 
                          frac_coord=True, sort_reverse=True):
        indexFinal = []
        distanceList = []
        if (indexInit == None and atomType == None):
            raise ValueError('Both index and atomType cannot be None')
        elif (indexInit == None and isinstance(atomType, str)):
            indexInit = return_index_of_atomtype(crystalStrucInit, [atomType,])
        for indexVal in indexInit:
            if frac_coord:
                coordInit = crystalStrucInit[indexVal].frac_coords
            else:
                raise NotImplementedError('Currently can only handle Fractional Coordinates')
            #print 'calc_displacement', indexVal, coordInit #TEST
            atomFinal = self.get_neighbors_to_coords(crystalStrucFinal,coordInit,
                                                     maxDisplacement,
                                                     crystalStrucInit[indexVal].species_string,
                                                     noNeighbors=1,
                                                     frac_coord=frac_coord)
            if len(atomFinal) > 0:
                indexValFinal = atomFinal[0][2]
                distanceVal = atomFinal[0][1]
                indexFinal.append(indexValFinal)
                distanceList.append([indexVal, indexValFinal, distanceVal])
            else:
                print 'AtomIndex: ', indexVal
                warnings.warn('Atom Corresponding to Index '+str(indexVal) +
                              'cannot be found')
        distanceList.sort(key=lambda x:x[2], reverse=sort_reverse)
        return distanceList

    def get_all_pair_distances(self, crystalStruc, atomType, pureLatticeConst, variationVal):
        noPairs = 2
        latticeVec = np.asarray(crystalStruc._lattice.abc)
        atomIndex = self.return_index_of_atomtype(crystalStruc, atomType)
        atomIndexList = []
        for i in atomType:
            atomIndexList += atomIndex[i]
        atomIndexList.sort()
        atomCombinations = [subset for subset in combinations(atomIndexList, noPairs)]
        atomDistancesList = []
        for combVal in atomCombinations:
            frac_coord1 = crystalStruc[combVal[0]].frac_coords
            frac_coord2 = crystalStruc[combVal[1]].frac_coords
            (distance, image) = crystalStruc._lattice.get_distance_and_image(frac_coord1,
                                                                             frac_coord2)
            coord1 = crystalStruc[combVal[0]].coords
            coord2 = crystalStruc[combVal[1]].coords + np.multiply(image, latticeVec) # works only for cubic supercells
            diffVec = coord2 - coord1
            unitDiffVec = diffVec/np.linalg.norm(diffVec)
            distance2 = np.linalg.norm(diffVec) #TESTING
            if distance < pureLatticeConst + variationVal:
                atomDistancesList.append([combVal, distance, unitDiffVec])
        return atomDistancesList

    def get_all_angle_pairs(self,crystalStruc, atomDistancesList, perfectAngle, variationVal):
        # THE LOGIC IS WEAK, Ideally you need to get the direction based on common-atom and then take dot product
        atomAngleList = []
        for i in range(0, len(atomDistancesList)):
            for j in range(i+1, len(atomDistancesList)):
                atomPair1 = atomDistancesList[i][0]
                atomPair2 = atomDistancesList[j][0]
                commonAtoms = list(set(atomPair1).intersection(atomPair2))
                if len(commonAtoms) > 0 :
                    vec1 = atomDistancesList[i][2]
                    vec2 = atomDistancesList[j][2]
                    angle = math.degrees(math.acos(np.dot(vec1,vec2)))
                    if ( perfectAngle-variationVal <= angle <= perfectAngle + variationVal):
                        atomAngleList.append([atomPair1, atomPair2, commonAtoms, angle])
                    elif ((360 - perfectAngle) - variationVal <= angle <=
                          (360 - perfectAngle) + variationVa):
                        angle = 360 - angle
                        atomAngleList.append([atomPair1, atomPair2, commonAtoms, angle])
        return atomAngleList
