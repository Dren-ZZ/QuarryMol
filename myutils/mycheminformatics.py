import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import DataStructs, rdFMCS, rdDistGeom, rdForceFieldHelpers
from rdkit.Chem.Descriptors import rdFingerprintGenerator
from rdkit.Chem.Draw import rdDepictor, rdMolDraw2D

import myutils.mypubchem as mypcp
import myutils.mygraphics as mygraphics

morgan_generator = rdFingerprintGenerator.GetMorganGenerator(
    radius=2, countSimulation=True, includeChirality=True
)

RAW_PALETTE_DARK = {
    1: (255, 255, 255),
    6: (255, 255, 255),
    7: (163, 234, 255),
    8: (255, 145, 145),
    9: (163, 255, 214),
    15: (255, 189, 138),
    16: (237, 224, 104),
    17: (83, 189, 131),
    35: (217, 138, 78),
    53: (245, 135, 255),
}

RAW_PALETTE_LIGHT={
    1: (0, 0, 0),
    6: (0, 0, 0),
    7: (0, 132, 255),
    8: (199, 0, 0),
    9: (54, 173, 163),
    15: (237, 111, 0),
    16: (184, 171, 0),
    17: (11, 125, 61),
    35: (115, 55, 0),
    53: (77, 0, 115),            
}

def draw_compound(
        mol, w, h, palette, 
        bg=mygraphics.rgb_rescaled(100, 100, 100, 0.05),
        highlight_color=(0, 0, 0, 0),
        highlight_atoms=[],
    ):
    mol_img = rdMolDraw2D.MolDraw2DCairo(w, h)

    rdDepictor.Compute2DCoords(mol)
    rdDepictor.StraightenDepiction(mol)

    draw_options = mol_img.drawOptions()
    draw_options.setBackgroundColour(bg)
    draw_options.setHighlightColour(highlight_color)

    rdkit_palette = mygraphics.create_color_palette(mol, palette, (0, 0, 0, 1))
    draw_options.updateAtomPalette(rdkit_palette)

    mol_img.DrawMolecule(mol, highlightAtoms=highlight_atoms,)
    mol_img.FinishDrawing()
    
    return mol_img.GetDrawingText()


def draw_3D_molecule(mol):
    mol_hs = Chem.AddHs(mol)
    rdDistGeom.EmbedMolecule(mol_hs, rdDistGeom.ETKDG())
    rdForceFieldHelpers.UFFOptimizeMolecule(mol_hs)
    return Chem.MolToMolBlock(mol_hs)


def generate_tanimoto_matrix(molecules):
    morgan_fps = [morgan_generator.GetCountFingerprint(m) for m in molecules]
    n = len(molecules)
    sim_matrix = np.full((n, n), np.nan)
    for i in range(n):
        sim_matrix[i, i] = 1.0
        for j in range(i+1, n):
            sim = DataStructs.cDataStructs.TanimotoSimilarity(morgan_fps[i], morgan_fps[j])
            sim_matrix[i, j] = sim
    
    return morgan_fps, sim_matrix



def show_sims(mol_x, mol_y, highlight_color):
    Chem.Kekulize(mol_x, clearAromaticFlags=True) # Removes aromatic ambiguity
    Chem.Kekulize(mol_y, clearAromaticFlags=True)

    params = rdFMCS.MCSParameters()
    params.BondCompare = rdFMCS.BondCompare.CompareOrderExact # Compares bonds
    params.AtomCompare = rdFMCS.AtomCompare.CompareElements # Compares atoms
    params.RingMatchesRingOnly = True # Ensures only rings match rings
    params.CompleteRingsOnly = True # Avoids partial ring matching
    params.MatchValences = True # Ensures valence compatibility
    params.MatchChiralTag = True # If stereochemistry matters
    
    mcs = rdFMCS.FindMCS([mol_x, mol_y], parameters=params)
    mcs_mol = Chem.MolFromSmarts(mcs.smartsString)
    
    match_x = mol_x.GetSubstructMatch(mcs_mol)
    match_y = mol_y.GetSubstructMatch(mcs_mol)
    target_atoms_x = list(match_x)
    target_atoms_y = list(match_y)

    mol_x_img = mygraphics.convert_to_png(
        draw_compound(
            mol_x, 800, 800, RAW_PALETTE_LIGHT,
            highlight_color=highlight_color,
            highlight_atoms=target_atoms_x,
        )
    )

    mol_y_img = mygraphics.convert_to_png(
        draw_compound(
            mol_y, 800, 800, RAW_PALETTE_LIGHT,
            highlight_color=highlight_color,
            highlight_atoms=target_atoms_y,
        )
    )
    return [mol_x_img, mol_y_img]

