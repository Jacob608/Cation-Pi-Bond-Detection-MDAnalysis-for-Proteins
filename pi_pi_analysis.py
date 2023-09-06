#!/usr/bin/env python
# coding: utf-8

# Count Pi Pi Bonds Using MDAnalysis Class Building

# Jacob Graham
# Northwestern University Department of Mechanical Engineering
# Sinan Keten's Computational Nanodynamics Research Laboratory

# 9/6/2023




# import numpy as np
# import math
# import os
# import glob
# import pandas as pd
# import MDAnalysis as mda
# import argparse
# import matplotlib.pyplot as plt
# import sys
# from MDAnalysis.analysis.base import (AnalysisBase, AnalysisFromFunction, analysis_class)
# from MDAnalysis.analysis import distances
# from statistics import mean



from MDAnalysis.analysis.base import AnalysisBase, Results
from MDAnalysis.lib.distances import capped_distance, calc_angles
import logging
import warnings
import numpy as np
from MDAnalysis.lib.distances import capped_distance, calc_angles
from MDAnalysis.lib.correlations import autocorrelation, correct_intermittency
from MDAnalysis.exceptions import NoDataError
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.analysis.hydrogenbonds.hbond_autocorrel import find_hydrogen_donors

class PiPiBondAnalysis(AnalysisBase):
    """
    Perform an analysis of pi pi bonds in a Universe.
    """
    
    def __init__(self, 
                 universe, 
                 pi_sel=None, 
                 p_p_cutoff=4.0,
                 update_selections=True):
        """Set up atom selections and geometric criteria for finding cation pi
        bonds in a Universe.
        
        P pi bond selection with 'pi_sel' may be achieved
        with either a *resname*, atom *name* combination, or when those are absent, with atom *type* selections.
        
        Parameters:
        -----------
        universe : Universe
            MDAnalysis Universe object
        pi_sel : str
            Selection string for the atoms participating in the quadropole (pi-system).
        p_p_cutoff : float
            Distance cutoff used for finding pi pi pairs.
        plane_p_p_angle_cutoff : float
            Aromatic plane-barycenter-cation angle cutoff for cation pi pairs in degrees.
        """
        
        self.u = universe # Initialize self.u as the input universe. 
        self._trajectory = self.u.trajectory # Initialize self._trajectory as the input self.u.trajectory
        self.pi_sel = pi_sel # Initialize self.pi_sel as pi_sel
        
        msg = ("{} is an empty selection string - no cation pi bonds will "
               "be found. This may be intended, but please check your "
               "selection."
              )

        val = getattr(self, 'pi_sel')
        if isinstance(val, str) and not val:
            warnings.warn(msg.format('pi_sel'))
        
        self.p_p_cutoff = p_p_cutoff # Initialize self.p_p_cutoff as input p_p_cutoff
        self.update_selections = update_selections # Initialize self.update_selections as update_selections
        self.results = Results() # Initialize self.results as Results()
        self.results.pi_pi_bonds = None # Initialize self.results.pi_pi_bonds as None.
    
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
        self.results.pi_pi_bonds = []

        # Select atom groups
        self._pi_sel = self._get_pi_system_atoms()
        
    def _single_frame(self):
        
        box = self._ts.dimensions
        
        # find cation and pi barycenter within cutoff distance of one another
        pi_rings = self._get_pi_system_atoms()
        barycenters=np.zeros([len(pi_rings),3])
        for i in range(len(pi_rings)):
            avg_pos = np.array(pi_rings[i].positions.mean(axis=0))
            barycenters[i]=avg_pos
        # capped_distance currently returns duplicates of every pi-pi pair.
        p_p_indices, p_p_distances = capped_distance(
            barycenters,
            barycenters,
            max_cutoff=self.p_p_cutoff,
            min_cutoff=1.0,
            box=box,
            return_distances=True
        )
        if np.size(p_p_indices) == 0:
            warnings.warn(
                "No pi pi bonds were found given p-p cutoff of "
                f"{self.p_p_cutoff} for "
                f"Pi System, {self.pi_sel} at step {self._ts}."
            )
        
        # Remove P-P pairs more than p_p_cutoff away from one another
        tmp_pi_rings=[pi_rings[x] for x in p_p_indices.T[1]]
        
        # Store data on number of cation pi bonds found at this frame
        self.results.pi_pi_bonds.append(len(tmp_pi_rings)/2) # This is just the number of cation pi bonds.
    
    def _conclude(self):
        
        self.results.pi_pi_bonds = np.asarray(self.results.pi_pi_bonds).T




