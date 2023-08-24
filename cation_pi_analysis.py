#!/usr/bin/env python
# coding: utf-8

# Count Cation Pi Bonds Using a Custom MDAnalysis Class 

# Jacob Graham
# Northwestern University Department of Mechanical Engineering
# Sinan Keten's Computational Nanodynamics Research Laboratory


import numpy as np
import MDAnalysis as mda
import matplotlib.pyplot as plt
from MDAnalysis.analysis.base import (AnalysisBase, AnalysisFromFunction, analysis_class)
from MDAnalysis.analysis import distances
from MDAnalysis.analysis.base import AnalysisBase, Results
from MDAnalysis.lib.distances import capped_distance, calc_angles
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
        """Set up atom selections and geometric criteria for finding cation pi
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
            Selection string for the atom that represents the cation involved in cation pi
            interactions.
        c_p_cutoff : float
            Distance cutoff used for finding cation pi pairs.
        plane_p_c_angle_cutoff : float
            Aromatic plane-barycenter-cation angle cutoff for cation pi pairs in degrees.
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
        self.update_selections = update_selections # Initialize self.update_selections as update_selections
        self.results = Results() # Initialize self.results as Results()
        self.results.cation_pi_bonds = None # Initialize self.results.cation_pi_bonds as None.
    
    def _get_pi_system_atoms(self):
        """ Finds cation pi pairs.
        
        Returns
        -------
        pi_system_atoms : AtomGroup
            AtomGroup corresponding to all atoms in pi systems.
            AtomGroups are ordered such that, if zipped, will
            produce a list of lists of pi system atoms.
        """
        pi_atoms = self.u.select_atoms(self.pi_sel)
        # Split pi_atoms into individual rings
        # Assume that there are exactly six carbons in each ring.
        c_per_ring = 6
        pi_rings=[] # Blank list becomes a list of atom groups each containing only atoms in an individual ring.
        for i in range(c_per_ring,66,c_per_ring):
            pi_rings.append(pi_atoms[(i-c_per_ring):i])
        
        return pi_rings
    
    def _prepare(self):
        self.results.cation_pi_bonds = []

        # Select atom groups
        self._pi_sel = self._get_pi_system_atoms()
        self._cations_sel = self.u.select_atoms(self.cations_sel,
                                                updating=self.update_selections)
            
    def _single_frame(self):
        
        box = self._ts.dimensions
        
        # find cation and pi barycenter within cutoff distance of one another
        pi_rings = self._get_pi_system_atoms()
        barycenters=np.zeros([len(pi_rings),3])
        for i in range(len(pi_rings)):
            avg_pos = np.array(pi_rings[i].positions.mean(axis=0))
            barycenters[i]=avg_pos
        c_p_indices, c_p_distances = capped_distance(
            self._cations_sel.positions,
            barycenters,
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
        
        # Store data on number of cation pi bonds found at this frame
        self.results.cation_pi_bonds.append(len(tmp_cations)) # This is just the number of cation pi bonds.
    
    def _conclude(self):
        
        self.results.cation_pi_bonds = np.asarray(self.results.cation_pi_bonds).T