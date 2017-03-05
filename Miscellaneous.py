import math
import numpy as np
from pymatgen.core.periodic_table import Element

class Miscellaneous:
    # https://stackoverflow.com/questions/38987/how-can-i-merge-two-python-dictionaries-in-a-single-expression
    def merge_two_dicts(x, y):
        '''Given two dicts, merge them into a new dict as a shallow copy.'''
        z = x.copy()
        z.update(y)
        return z


    #Code adopted from https://stackoverflow.com/questions/6802577/python-rotation-of-3d-vector
    def rotation_matrix(axis, theta):
        """
        Return the rotation matrix associated with counterclockwise rotation about
        the given axis by theta radians.
        """
        axis = np.asarray(axis)
        theta = np.asarray(theta)
        axis = axis/math.sqrt(np.dot(axis, axis))
        a = math.cos(theta/2)
        b, c, d = -axis*math.sin(theta/2)
        aa, bb, cc, dd = a*a, b*b, c*c, d*d
        bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
        return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                         [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                         [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])

        """
    TODO:
    """
    # def get_common_sublists(list1, noOfSubListsReturn):
    #     c = Counter(tuple(x) for x in list1)
    #     return c.most_common(noOfSubListsReturn)

    def compare_lists(list1, list2, accuracy=1.0E-02):
        flagVal = False
        if (len(list1) == len(list2)):
            diffList = [abs(i-j) for i, j in zip(list1, list2)]
            diffList.sort(reverse=True)
            if (diffList[0] <= accuracy):
                flagVal = True
        else:
            raise ValueError('Length of two lists not the same')

        return flagVal

    def zerolistmaker(n):
    #SO-Link: https://stackoverflow.com/questions/8528178/list-of-zeros-in-python
        listofzeros = [0] * n
        return listofzeros

    @staticmethod
    def isValueIn(val, valRange, accuracy=1.0E-03):
        valueIn = False
        if (valRange[0] <= val <= valRange[1]):
            valueIn = True
        return valueIn

    def find_recursive(self, needle, haystack):
        """
        Find index of an item in list of list
        """
        # WebLink https://stackoverflow.com/questions/14124203/index-of-item-in-a-list-of-lists-of-lists-python
        for index, item in enumerate(haystack):
            if not isinstance(item, str):
                try:
                    path = self.find_recursive(needle, item)
                    if path is not None:
                        return (index, ) + path
                except TypeError:
                    pass
            if needle == item:
                return index,
        return None

    def get_nonrecvalues(arr1, dtype='float'):
        nrValList = [arr1[0], ]
        nrIndexList = [0, ]
        for i in range(1, len(arr1)):
            if arr1[i] != arr1[i-1]:
                nrValList.append(arr1[i])
                nrIndexList.append(i)
        nrValArr = np.array(nrValList, dtype=dtype)
        nrIndexArr = np.array(nrIndexList, dtype=int)
        return nrValArr, nrIndexArr

    @staticmethod
    def get_ionic_radius(element, oxidationNo):
        if element == 'Mn' and oxidationNo == 4:
            radius = 0.67
        else:
            radius = Element(element).ionic_radii[oxidationNo]
        return radius

    def get_centroid(self, arr1):
        return np.mean(arr1, axis=0)

    def calc_pervoskite_tolerance_factor(self, siteA, oxdA, siteB,
                                         oxdB, siteC, oxdC):
        ionic_radA = self.get_ionic_radius(siteA, oxdA)
        ionic_radB = self.get_ionic_radius(siteB, oxdB)
        ionic_radC = self.get_ionic_radius(siteC, oxdC)
        toleranceFactor = ((ionic_radA + ionic_radC) /
                           (math.sqrt(2) *
                            (ionic_radB + ionic_radC)))
        return toleranceFactor

    # Hacked from http://www.saltycrane.com/blog/2007/12/how-to-convert-list-of-dictionaries-to/
    def listdict_to_listlist(self, lod, keylist):
        
        """Converts a list of dictionaries to a list of lists using the
        values of the dictionaries. Assumes that each dictionary has the
        same keys. 
           lod: list of dictionaries
           keylist: list of keys, ordered as desired
           Returns: a list of lists where the inner lists are rows. 
        i.e. returns a list of rows. """
        # declarative/functional approach
        return [[row[key] for key in keylist] for row in lod]

    @staticmethod
    def get_index_not_in_range(arr1, minVal, maxVal):
        indexList = np.empty(0, dtype=int)
        if all([minVal, maxVal]):
            indexList = np.where(np.logical_not((arr1 >= minVal) &
                                                (arr1 <= maxVal)))
        else:
            if minVal is not None:
                indexList = np.where(np.logical_not((arr1 >= minVal)))
            elif maxVal is not None:
                indexList = np.where(np.logical_not((arr1 <= maxVal)))
        return indexList

    def modify_from_value(self, data, colIndexList, colValRangeList,
                          inRange=False):
        removeIndex = []
        for cnt, val in enumerate(data):
            for colCnt, colIndex in enumerate(colIndexList):
                if (val[colIndex] < colValRangeList[colCnt][0] or
                    val[colIndex] > colValRangeList[colCnt][1]):
                    removeIndex.append(cnt)
        removeIndex = sorted(list(set(removeIndex)))
        outData = [x for cnt, x in enumerate(data)
                   if cnt not in removeIndex]
        return outData
