#!/usr/bin/env python
# coding: utf-8

# Count Cation Pi Bonds Using a Custom MDAnalysis Class 

# Jacob Graham
# Northwestern University Department of Mechanical Engineering
# Sinan Keten's Computational Nanodynamics Research Laboratory


import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis.base import (AnalysisBase, AnalysisFromFunction, analysis_class)
from MDAnalysis.analysis import distances
from MDAnalysis.analysis.base import AnalysisBase, Results
from MDAnalysis.lib.distances import capped_distance
import logging
import warnings
from MDAnalysis.lib.correlations import autocorrelation, correct_intermittency
from MDAnalysis.exceptions import NoDataError
from MDAnalysis.core.groups import AtomGroup

class CationPiBondAnalysis(AnalysisBase):
    """
    Perform an analysis of cation pi bonds in a Universe.
    """

    def __init__(self, 
                 universe, 
                 pi_sel=None, 
                 cations_sel=None,
                 c_p_cutoff=5.0,
                 update_selections=True):
        """
        Set up atom selections and geometric criteria for finding cation pi
        bonds in a Universe.
        
        Cation pi bond selections with 'pi_sel' and 'cations_sel' may be achieved
        with either a *resname*, atom *name* combination, or when those are absent, with atom *type* selections.
        
        Parameters:
        -----------
        universe : Universe
            MDAnalysis Universe object
        pi_sel : str
            Selection string for the atoms participating in the quadropole (pi-system).
        cations_sel : str
            Selection string for the atom that represents the cation(s) involved in cation pi
            interactions.
        c_p_cutoff : float
            Distance cutoff between ring positional centers and cation used for finding cation pi pairs.
        """
        
        self.u = universe # Initialize self.u as the input universe. 
        self._trajectory = self.u.trajectory # Initialize self._trajectory as the input self.u.trajectory
        self.pi_sel = pi_sel # Initialize self.pi_sel as pi_sel
        self.cations_sel = cations_sel # Initialize self.cations_sel as the input cations_sel
        
        msg = ("{} is an empty selection string - no cation pi bonds will "
               "be found. This may be intended, but please check your "
               "selection."
              )
        for sel in ['pi_sel','cations_sel']:
            val = getattr(self, sel)
            if isinstance(val, str) and not val:
                warnings.warn(msg.format(sel))
        
        self.c_p_cutoff = c_p_cutoff # Initialize self.c_p_cutoff as input c_p_cutoff
        self.update_selections = update_selections # Initialize self.update_selections as input update_selections
        self.results = Results() # Initialize self.results as input Results()
        self.results.cation_pi_bonds = None # Initialize self.results.cation_pi_bonds as None.
    
    def _get_pi_system_atoms(self):
        """
        Generate a list, where each element is an atomgroup containing only atoms for a specific ring.
        Currently, this function assumes that there are exactly six atoms in each ring and that atoms
        in the same ring are selected consecutively by the select_atoms() MDAnalysis function.

        Returns
        -------
        pi_rings : python list
            List of AtomGroups, where each AtomGroup contains the 6 atoms in a particular ring.
        """
        # Generate initial atom selection of pi system atoms based on pi_sel string.
        pi_atoms = self.u.select_atoms(self.pi_sel)
        # Split pi_atoms such that there is a unique AtomGroup for each individual ring.
        # Must assume that there are exactly six carbons in each ring.

        c_per_ring = 6
        pi_rings=[] # Blank list becomes a list of atom groups each containing only atoms in an individual ring.
        # The below for loop splits pi_atoms into AtomGroups of exactly 6 atoms.
        for i in range(c_per_ring,len(pi_atoms)+c_per_ring,c_per_ring):
            pi_rings.append(pi_atoms[(i-c_per_ring):i])
        
        return pi_rings
    
    def _prepare(self):
        # Initialize emtpy lists to store results.
        self.results.cation_pi_bonds = [[],[],[],[]]
        self.results.cation_pi_bond_counts = []

        
        # Select atom groups
        self._pi_sel = self._get_pi_system_atoms()
        self._cations_sel = self.u.select_atoms(self.cations_sel,
                                                updating=self.update_selections)
            
    def _single_frame(self):
        
        box = self._ts.dimensions
        
        # find cation and pi system center (average of 6 positions) within cutoff distance of one another.
        pi_rings = self._get_pi_system_atoms()
        centers=np.zeros([len(pi_rings),3])
        # Replace each row in the array 'centers' with the arithmetic mean position of all members of one ring (centroid).
        for i in range(len(pi_rings)):
            avg_pos = np.array(pi_rings[i].positions.mean(axis=0))
            centers[i]=avg_pos
        # Find C and P within cutoff distance of one another.
        # min_cutoff = 1.0 as an atom cannot form a cation pi bond with itself.
        c_p_indices, c_p_distances = capped_distance(
            self._cations_sel.positions,
            centers,
            max_cutoff=self.c_p_cutoff,
            min_cutoff=1.0,
            box=box,
            return_distances=True
        )

        if np.size(c_p_indices) == 0:
            warnings.warn(
                "No cation pi bonds were found given c-p cutoff of "
                f"{self.c_p_cutoff} between Cation, {self.cations_sel}, and "
                f"Pi System, {self.pi_sel} at step {self._ts}."
            )
        
        # Remove C-P pairs more than c_p_cutoff away from one another
        tmp_cations=self._cations_sel[c_p_indices.T[0]]
        tmp_pi_rings=[pi_rings[x] for x in c_p_indices.T[1]]

        # Set output to include frame, atom selection containing cation atom, atom group containg the 6 members of the pi ring, distance between cation atom and the centroid of the pi ring
        self.results.cation_pi_bond_counts.append(len(tmp_cations)) # Append the total number of cation pi bonds detected in this frame.
        
        self.results.cation_pi_bonds[0].extend(np.full_like(tmp_cations, self._ts.frame)) # frame
        self.results.cation_pi_bonds[1].extend(tmp_cations) # atom selection containing cation atom
        self.results.cation_pi_bonds[2].extend(tmp_pi_rings) # atom group containg the 6 members of the pi ring
        self.results.cation_pi_bonds[3].extend(c_p_distances) # distance between cation atom and the centroid of the pi ring

    def _conclude(self):
        
        self.results.cation_pi_bonds = np.asarray(self.results.cation_pi_bonds, dtype=object).T