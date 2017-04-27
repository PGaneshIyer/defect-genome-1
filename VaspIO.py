import os, sys
import re
import shutil
from enum import Enum

from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.core.periodic_table import Element
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.core import Spin, OrbitalType

from InputOutputFns import InputOutput


class VaspIO:
    def return_orbital_spin_enum(self, str1, valType='orbital'):
        dataType = None
        if (valType == 'orbital' and isinstance(str1, str)):
            if (str1 == 's'):
                dataType = OrbitalType.s
            elif (str1 == 'p'):
                dataType = OrbitalType.p
            elif (str1 == 'd'):
                dataType = OrbitalType.d
            elif (str1 == 'f'):
                dataType = OrbitalType.f
            else:
                raise ValueError('str1 is not a valid OrbitalType')
        elif (valType == 'spin' and isinstance(str1, str)):
            if (str1 == '1'):
                dataType = Spin.up
            elif (str1 == '-1'):
                dataType = Spin.down
            else:
                raise ValueError('str1  is not a valid Spin')
        else:
            raise ValueError('valType is neither OrbitalType nor Spin')
        return dataType

    def get_force_vasp(filePath, string='g(F)'):
        string2 = 'trialstep'
        forceVal = None
        stepVal = None
        with open(filePath, 'r') as f:
            for line in reversed(f.readlines()):
                if string in line:
                    line_list = re.split('=', line)
                    for cnt, val in enumerate(line_list):
                        if string in val:
                            forceList = re.split(' ', line_list[cnt+1])
                            for val2 in forceList:
                                try:
                                    forceVal = float(val2)
                                    break
                                except:
                                    continue
                        elif string2 in val:
                            stepList = re.split(' |\)', line_list[cnt+1])
                            for val2 in stepList:
                                try:
                                    stepVal = float(val2)
                                    break
                                except:
                                    continue
                    break
        print line, forceVal, stepVal
        # print forceList, forceVal #TEST
        # print stepList, stepVal   #TEST
        return forceVal, stepVal, line

    def get_potim_vasp(filePath, lineStr=None, string='POTIM'):
        with open(filePath, 'r') as f:
            for line in f:
                if string in line:
                    val = float(line.split()[2])
                    break
        # if string in lineStr:
        #     val = 0.5
        return val

    def get_partial_dos(self, filePath, atomType, orbitalType=None,
                        spinType=None, normalized=True, energyShift=True):
        """
        Returns energy and partial Density of States corresponding to an atom, orbital, spin combination

        Args:
        :param filePath (str): path to vasprun.xml file
        :param atomType (str): Atom Symbol (only one atom at at time)
        :param orbitalType : Orbital of interest (only one orbital at a time, 
                             if none provided all orbitals are summed and returned)
        :param spinType : Spin of interest (only one spin at a time,
                             if none provided all spins are summed and returned)
        :param normalized (bool): To normalize the P-DoS by the number of atoms of the particular type
        :param energyShift (bool): Returns energies array shifted by the Fermi energy.
        
        Returns:
        pDos (numpy arr): partial density of states of the specified atom, orbital, spin combination
        energies (numpy arr): energy values employed in the VASP calculation.
        """
        
        vaspObj = Vasprun(filePath)
        vaspCompleteDos = vaspObj.complete_dos
        energies = vaspCompleteDos.energies
        efermi = vaspObj.efermi
        print 'Atom: ', atomType
        print 'Fermi Energy: ', efermi
        if energyShift:
            energies -= efermi
        if orbitalType:
            print 'Orbital: ', orbitalType
            if isinstance(orbitalType, OrbitalType):
                pass
            elif isinstance(orbitalType, str):
                orbitalType = self.return_orbital_spin_enum(orbitalType,
                                                            'orbital')
            else:
                raise ValueError('orbitalType is not a valid')
            atomDos = vaspCompleteDos.get_element_spd_dos(atomType)[orbitalType]
        else:
            atomDos = vaspCompleteDos.get_element_dos()[Element(atomType)]
        if spinType:
            print 'Spin: ', spinType
            if isinstance(spinType, Spin):
                pass
            elif isinstance(spinType, str):
                spinType = self.return_orbital_spin_enum(spinType, 'spin')
            else:
                raise ValueError('spinType not valid')
            pDos = atomDos.get_densities(spin=spinType)
        else:
            pDos = atomDos.get_densities()

        if normalized:
            noAtoms = len([x for x in vaspObj.final_structure
                           if x.species_string in atomType])
            print 'No. of Atoms: ', noAtoms
            pDos = pDos/noAtoms
        return pDos, energies

    def poscar_cart_to_direct(self, inputFilePath, outputFilePath,
                              vasp4Compatible=False):
        """
        TODO: replace default python file io with monty.io
        """
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
        """
        TODO: replace default python file io with monty.io
        """
        noDecimals = 6
        direct = False
        crystalStruc = Structure.from_file(inputFilePath)
        poscarFile = Poscar(crystalStruc)
        poscarString = poscarFile.get_string(direct, vasp4Compatible,
                                             noDecimals)
        with open(outputFilePath, "w") as f:
            f.write(poscarString)

