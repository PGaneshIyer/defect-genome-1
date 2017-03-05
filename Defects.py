import os,sys
from pymatgen.core.structure import Structure
from pymatgen.io.vasp import Poscar
from pymatgen.core.periodic_table import Element
#VaspIO
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.core import Spin, OrbitalType


class CreatePointDefects:

    def create_substitutional_dopant_structure(self, origStruc, keyVal,
                                               atomIndices, dopantAtomType):
        """Substitutes dopant at specified indices in structure.

        """

        #  Create structure.
        structure = origStruc

        #  Make substitution.
        for index in atomIndices:
            structure.replace(index, dopantAtomType)

        # Get initial key values.
        lattice_system = keyVal[0]
        vacancyAtoms = keyVal[4]
        noVacancies = keyVal[5]
        noCells = keyVal[6]
        interstitialAtoms = keyVal[7]
        noInterstitials = keyVal[8]

        #  Create dopantAtoms list.
        dopantAtoms = keyVal[2]
        noDopants = []

        if dopantAtomType not in dopantAtoms:
            dopantAtoms.append(dopantAtomType)

        # Find number of Dopants.
        for dopant in dopantAtoms:
            dopant_count = 0
            for site in structure:
                species = site.specie
                if str(dopant) == str(species):
                    dopant_count += 1
            noDopants.append(dopant_count)

        # Create basic atom list.
        basic_atoms = keyVal[1]

        for site in structure:
            species = str(site.specie)
            if species not in basic_atoms and species not in dopantAtoms:
                basic_atoms.append(species)

        # Create key.
        keyVal = self.structureKey(lattice_system, basic_atoms,
                                   dopantAtoms, noDopants,
                                   vacancyAtoms, noVacancies,
                                   noCells, interstitialAtoms,
                                   noInterstitials)

        return structure, keyVal

    def create_vacancy_defect_structure(self, origStruc, keyVal,
                                        vacancyIndices):
        """Creates vacancies at specified indices.

        """

        # Verify vacancyIndices is a list.
        if isinstance(vacancyIndices, list) is False:
            print "vacancyIndices must be a list."
            sys.exit(1)

        # Create structure.
        structure = origStruc

        #  Get initial key values.
        lattice_system = keyVal[0]
        basic_atoms = keyVal[1]
        dopantAtoms = keyVal[2]
        noDopants = keyVal[3]
        vacancyAtoms = keyVal[4]
        noVacancies = keyVal[5]
        noCells = keyVal[6]
        interstitialAtoms = keyVal[7]
        noInterstitials = keyVal[8]

        #  Create vacancyAtoms list.
        for index in vacancyIndices:
            site = structure[index]
            species = str(site.specie)
            if species not in vacancyAtoms:
                vacancyAtoms.append(species)

        # Create noVacancies list.
        preexisting_vacancy_index = 0
        for atom in vacancyAtoms:
            vacancy_count = 0
            index = 0
            for site in structure:
                species = str(site.specie)
                if str(atom) == str(species) and index in vacancyIndices:
                    vacancy_count += 1
                index += 1
            # Append vacancy_count to noVacancies.  If previous vacancies exist for atom,
            # add vacancy_count to previous vacancy_count.
            try:
                if noVacancies[preexisting_vacancy_index] > 0:
                    new_vacancy_count = noVacancies[preexisting_vacancy_index] + vacancy_count
                    noVacancies[preexisting_vacancy_index] = new_vacancy_count
            except:
                noVacancies.append(vacancy_count)
            preexisting_vacancy_index += 1

        # Remove vacancy sites.
        structure.remove_sites(vacancyIndices)

        #  Create key.
        keyVal = self.structureKey(lattice_system, basic_atoms,
                                   dopantAtoms, noDopants,
                                   vacancyAtoms, noVacancies,
                                   noCells, interstitialAtoms,
                                   noInterstitials)

        return structure, keyVal

    def create_interstitial_dopant_structure(self, origStruc, keyVal,
                                             atomsToInsert,
                                             neighboringIndices,
                                             distancesFromIndices,
                                             directionsFromIndices):
        """Creates interstitial dopants near neighboring indices.  Interstitial dopants are created at
        a specified distance and direction from the index.

        """

        # Verify atoms is a list.
        if isinstance(atomsToInsert, list) is False:
            print "atomsToInsert must be a list."
            sys.exit(1)

        # Verify interstitialIndices, distanceFromIndices and directionsFromIndices have the same length.
        if (len(neighboringIndices) != len(distancesFromIndices)) or (
                    len(neighboringIndices) != len(directionsFromIndices)) or (
                    len(distancesFromIndices) != len(directionsFromIndices)):
            print("neighboringIndices, distanceFromIndices and directionsFromIndices must be the same length.")
            sys.exit(1)

        # Create structure.
        structure = origStruc

        #  Get initial key values.
        lattice_system = keyVal[0]
        basic_atoms = keyVal[1]
        dopantAtoms = keyVal[2]
        noDopants = keyVal[3]
        vacancyAtoms = keyVal[4]
        noVacancies = keyVal[5]
        noCells = keyVal[6]
        interstitialAtoms = keyVal[7]
        noInterstitials = keyVal[8]

        # Get interstitial atom species.  Currently Hardcoded to assume their is one only one interstitial atom species.
        interstitialAtom = atomsToInsert[0]

        #  Insert interstitials.
        for key, index in enumerate(neighboringIndices):
            coordinatesAtIndex = structure[index].coords
            distance = distancesFromIndices[key]
            directionVector = directionsFromIndices[key]

            #  Check that direction is normalized.
            vectorLength = (directionVector[0] ** 2 +
                            directionVector[1] ** 2 + directionVector[2] ** 2) ** .5

            #  If not normalized, normalize vector.
            if vectorLength != 1:
                directionVector[0] /= vectorLength
                directionVector[1] /= vectorLength
                directionVector[2] /= vectorLength

            # Get coordinates to place interstitial atom.
            interstitialCoordinates = []
            for key2, component in enumerate(directionVector):
                distanceTimesComponent = component * distance
                interstitialComponent = coordinatesAtIndex[key2] + distanceTimesComponent
                interstitialCoordinates.append(interstitialComponent)

            # Insert interstitial atom into structure at interstitialCoordinates.
            structure.append(interstitialAtom, interstitialCoordinates, True, True)

        # Create key.
        #  Add interstitialAtom to interstitialAtoms if not in list.
        if interstitialAtom not in interstitialAtoms:
            interstitialAtoms.append(interstitialAtom)

        # Find number of interstitials.
        for atom in interstitialAtoms:
            atom_count = 0
            for site in structure:
                species = site.specie
                if str(atom) == str(species):
                    atom_count += 1
            noInterstitials.append(atom_count)

        keyVal = self.structureKey(lattice_system, basic_atoms,
                                   dopantAtoms, noDopants,
                                   vacancyAtoms, noVacancies,
                                   noCells, interstitialAtoms,
                                   noInterstitials)

        return structure, keyVal

    def move_interstitial_dopant(self, origStruc, keyVal, atomsToMove,
                                 indicesToMove,
                                 distancesFromIndices,
                                 directionsFromIndices):
        """Move an interstitial atom to another location.

        """

        # Verify atomsToMove is a list.
        if isinstance(atomsToMove, list) is False:
            print "atomsToMove must be a list."
            sys.exit(1)

        if (len(indicesToMove) != len(distancesFromIndices)) or (
                    len(indicesToMove) != len(directionsFromIndices)) or (
                    len(distancesFromIndices) != len(directionsFromIndices)):
            print "indicesToMove, distanceFromIndices and directionsFromIndices must be the same length."
            sys.exit(1)

        # Create structure.
        structure = origStruc

        # Get interstitial atom species.  Currently Hardcoded to assume their is one only one interstitial atom species.
        interstitialAtom = atomsToMove[0]

        #  Verify that indicesToMove match atoms.
        for index in indicesToMove:
            species = structure[index].specie
            if str(interstitialAtom) != str(species):
                print "indicesToMove does not match atomsToMove"
                sys.exit(1)

        #  Get initial key values.
        lattice_system = keyVal[0]
        basic_atoms = keyVal[1]
        dopantAtoms = keyVal[2]
        noDopants = keyVal[3]
        vacancyAtoms = keyVal[4]
        noVacancies = keyVal[5]
        noCells = keyVal[6]
        interstitialAtoms = keyVal[7]
        noInterstitials = keyVal[8]

        #  Verify atomsToMove are interstitial.
        if interstitialAtom not in interstitialAtoms:
            print ("%s is not an interstitial atom." % str(interstitialAtom))
            sys.exit(1)

        #  Verify that you are not trying to move more interstitials than currently exist.
        if len(indicesToMove) > noInterstitials:
            print "IndicesToMove is greater than the number of interstitials."
            sys.exit(1)

        #  Insert interstitials.
        for key, index in enumerate(indicesToMove):
            coordinatesAtIndex = structure[index].coords
            distance = distancesFromIndices[key]
            directionVector = directionsFromIndices[key]

            #  Check that direction is normalized.
            vectorLength = (directionVector[0] ** 2 +
                            directionVector[1] ** 2 + directionVector[2] ** 2) ** .5

            #  If not normalized, normalize vector.
            if vectorLength != 1:
                directionVector[0] /= vectorLength
                directionVector[1] /= vectorLength
                directionVector[2] /= vectorLength

            # Get coordinates to move interstitial atom.
            interstitialCoordinates = []
            for key2, component in enumerate(directionVector):
                distanceTimesComponent = component * distance
                interstitialComponent = coordinatesAtIndex[key2] + distanceTimesComponent
                interstitialCoordinates.append(interstitialComponent)

            # Insert interstitial atom into structure at interstitialCoordinates.
            structure.append(interstitialAtom, interstitialCoordinates, True, True)

        structure.remove_sites(indicesToMove)

        #  Create unchanged key.
        keyVal = self.structureKey(lattice_system, basic_atoms,
                                   dopantAtoms, noDopants,
                                   vacancyAtoms, noVacancies,
                                   noCells, interstitialAtoms,
                                   noInterstitials)

        return structure, keyVal


class PointDefectEnergies:
    def get_intst_formation_energy(self, energy_OHStruc,
                                   energy_PureStruc,
                                   efermi_PureStruc,
                                   chargeVal,
                                   chemPot_H):
        # chemPot_H (from H20, PBE) = -4.880197265
        formEnergy = ((energy_OHStruc +
                       chargeVal*efermi_PureStruc) -
                      (energy_PureStruc + chemPot_H))
        return formEnergy

    def get_vacancy_formation_energy(self, energy_VacStruc,
                                     energy_PureStruc,
                                     efermi_PureStruc,
                                     chargeVal,
                                     chemPot_O):
        # chemPot_O (from O2, PBE) = -4.4633148
        formEnergy = ((energy_VacStruc + chemPot_O +
                       chargeVal*efermi_PureStruc) -
                      energy_PureStruc)
        return formEnergy

    def get_dopant_formation_energy(self, energy_DopedStruc,
                                    bsite_EnergyPerAtom,
                                    energy_PureStruc,
                                    efermi_PureStruc,
                                    chargeVal,
                                    chemPot_Dopant,
                                    chemPot_O):
        # FIXME chemPot_BSite only works for 2/4 Perovskites
        # chemPot_Dopant (Y from Y2O3, mp-2652) = -16.0880278
        chemPot_BSite = 3*bsite_EnergyPerAtom - 2*chemPot_O
        formEnergy = ((energy_DopedStruc + chemPot_BSite +
                       chargeVal*efermi_PureStruc) -
                      (energy_PureStruc + chemPot_Dopant))
        return formEnergy

    def get_pure_formation_energy(self, energy_PureStruc,
                                  asite_EnergyPerAtom,
                                  bsite_EnergyPerAtom,
                                  noCells):
        # FIXME: noAtoms will only work for 2/4 Perovskites
        noAtoms = 5
        formEnergy = ((energy_PureStruc/noCells -
                       (3*bsite_EnergyPerAtom +
                        2*asite_EnergyPerAtom))/noAtoms)
        return formEnergy

    def get_dopant_intst_interaction_energy(self, energy_PureStruc,
                                            energy_DopedStruc,
                                            energy_OHStruc,
                                            energy_Dop_OHStruc):
        intrEnergy = ((energy_Dop_OHStruc + energy_PureStruc) -
                      (energy_DopedStruc + energy_OHStruc))
        return intrEnergy
