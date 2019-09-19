__author__ = 'mnowotka'

# ----------------------------------------------------------------------------------------------------------------------

from rdkit.Chem import Descriptors, rdChemReactions, rdmolops, CombineMols, SanitizeMol, MolToSmiles
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from chembl_beaker.beaker.utils.functional import _call, _apply
from chembl_beaker.beaker.utils.io import _parseMolData

MAX_PASSES = 10

RDK_DESC_LIST = ['qed', 'MolWt', 'ExactMolWt', 'TPSA', 'HeavyAtomCount', 'NumHAcceptors', 'NumHDonors', 'NumRotatableBonds', 'MolLogP']

# ----------------------------------------------------------------------------------------------------------------------


def apply_rxn(mol, rxn):
    mols = [mol]
    changed = False
    for n_pass in range(MAX_PASSES):
        products = {}
        for m in mols:
            for product in [x[0] for x in rxn.RunReactants((m,))]:
                try:
                    SanitizeMol(product)
                    smiles = MolToSmiles(product, isomericSmiles=True)
                except ValueError as error:
                    # assuming an unphysical molecule has been generated
                    continue
                if smiles in products:
                    # keep only new structures
                    continue
                products[smiles] = product
        if products:
            changed = True
            # update list of mols
            mols = list(products.values())
        else:
            break
    return mols[0] if changed else mol

# ----------------------------------------------------------------------------------------------------------------------

def neutralise_sulphoxide(mol):
    smirks = '[S+1:1][O-1:2]>>[S+0:1]=[O-0:2]'
    rxn = rdChemReactions.ReactionFromSmarts(smirks)
    frags = rdmolops.GetMolFrags(mol, asMols=True)
    n_frags = list(filter(lambda x: x is not None,
                          [apply_rxn(frag, rxn) for frag in frags]))
    if len(n_frags) == 1:
        n_mol = n_frags[0]
    elif len(n_frags) == 2:
        n_mol = CombineMols(*n_frags)
        SanitizeMol(n_mol)
    else:
        n_mol = CombineMols(n_frags[0], n_frags[1])
        for i in range(2, len(n_frags)):
            n_mol = CombineMols(n_mol, n_frags[i])
        SanitizeMol(n_mol)
    return n_mol

# ----------------------------------------------------------------------------------------------------------------------

def remove_isotope_info(mol):
    """
    Removing isotpe information RDKit MolWt will calc the average weight and ExactMolWt the monoisotopic weight.
    MolWt takes average weight of each atom only if no isotope information is given.
    ExactMolWt takes most abundant isotope weight for each atom only if no isotope information is given.
    """
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)
    return mol

# ----------------------------------------------------------------------------------------------------------------------

def ro3_pass(mw_freebase, hba, hbd, alogp, rtb, psa):
    """
    It is suggested that compounds that pass all these criteria are more likely to be hits in fragment screening.
    molecular weight <=300,
    number of hydrogen bond donors <=3,
    number of hydrogen bond acceptors <=3
    AlogP <=3.
    RTB <=3
    PSA<=60
    From:
    A ‘Rule of Three’ for fragment-based lead discovery? Miles Congreve, Robin Carr,
    Chris Murray and Harren Jhoti. Drug Discovery Today, 2003,8(19), 876-877
    """
    ro3pass = None
    if mw_freebase and hba and hbd and alogp and rtb and psa:
        if mw_freebase <= 300 and hba <= 3 and hbd <= 3 and alogp <= 3 and rtb <= 3 and psa <= 60:
            ro3pass = 1
        else:
            ro3pass = 0
    return ro3pass

# ----------------------------------------------------------------------------------------------------------------------

def num_ro5_violations(alogp, mw_freebase, hba, hbd):
    """
    Number of properties defined in Lipinski’s Rule of 5 (RO5) that the compound fails.
    Conditions which violate the RO5 are:
        Molecular weight>=500
        AlogP>=5
        HBD>5
        HBA>10
    From:
    Lipinski, C. A.; Lombardo, F.; Dominy, B. W.; Feeney, P. J. Experimental and Computational Approaches to Estimate
    Solubility and Permeability in Drug Discovery and Development Settings. Adv. Drug Deliv. Rev., 1997, 23, 3-25.
    """
    violations = None
    if alogp or mw_freebase or hba or hbd:
        violations = 0
    if alogp and alogp >= 5:
        violations += 1
    if mw_freebase and mw_freebase >= 500:
        violations += 1
    if hba and hba > 10:
        violations += 1
    if hbd and hbd > 5:
        violations += 1
    return violations

# ----------------------------------------------------------------------------------------------------------------------

def _desc(mol, name):
    if name and isinstance(name, str) and hasattr(Descriptors, name):
        return getattr(Descriptors, name)(mol)

# ----------------------------------------------------------------------------------------------------------------------

def _desc_list(mol):
    mol = remove_isotope_info(mol)
    mol = neutralise_sulphoxide(mol)
    descriptors = dict()
    for name, fn in Descriptors.descList:
        if name in RDK_DESC_LIST:
            descriptors[name] = fn(mol)
    if 'MolecularFormula' not in descriptors:
        descriptors['MolecularFormula'] = CalcMolFormula(mol)
    descriptors['Ro3Pass'] = ro3_pass(descriptors['MolWt'],
                                      descriptors['NumHAcceptors'],
                                      descriptors['NumHDonors'],
                                      descriptors['MolLogP'],
                                      descriptors['NumRotatableBonds'],
                                      descriptors['TPSA'])
    descriptors['NumRo5'] = num_ro5_violations(descriptors['MolLogP'],
                                               descriptors['MolWt'],
                                               descriptors['NumHAcceptors'],
                                               descriptors['NumHDonors'])
    descriptors['MonoisotopicMolWt'] = descriptors.pop('ExactMolWt')
    return descriptors

# ----------------------------------------------------------------------------------------------------------------------


def _getNumAtoms(data, sanitize=True, removeHs=True, strictParsing=True):
    return _call(_parseMolData(data, sanitize=sanitize, removeHs=removeHs,
                               strictParsing=strictParsing), "GetNumAtoms")

# ----------------------------------------------------------------------------------------------------------------------


def _getNumBonds(data, sanitize=True, removeHs=True, strictParsing=True):
    return _call(_parseMolData(data, sanitize=sanitize, removeHs=removeHs,
                               strictParsing=strictParsing), "GetNumBonds")

# ----------------------------------------------------------------------------------------------------------------------


def _getLogP(data, sanitize=True, removeHs=True, strictParsing=True):
    return _apply(_parseMolData(data, sanitize=sanitize, removeHs=removeHs,
                                strictParsing=strictParsing), _desc, 'MolLogP')

# -----------------------------------------------------------------------------------------------------------------------


def _getTPSA(data, sanitize=True, removeHs=True, strictParsing=True):
    return _apply(_parseMolData(data, sanitize=sanitize, removeHs=removeHs,
                                strictParsing=strictParsing), _desc, 'TPSA')

# ----------------------------------------------------------------------------------------------------------------------


def _getMolWt(data, sanitize=True, removeHs=True, strictParsing=True):
    return _apply(_parseMolData(data, sanitize=sanitize, removeHs=removeHs,
                                strictParsing=strictParsing), _desc, 'MolWt')

# ----------------------------------------------------------------------------------------------------------------------


def _getDescriptors(data, sanitize=True, removeHs=True, strictParsing=True):
    return _apply(_parseMolData(data, sanitize=sanitize, removeHs=removeHs,
                                strictParsing=strictParsing), _desc_list)

# ----------------------------------------------------------------------------------------------------------------------
