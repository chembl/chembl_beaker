__author__ = 'efelix'

# ----------------------------------------------------------------------------------------------------------------------

from rdkit.Chem import Descriptors, rdChemReactions, rdmolops, CombineMols, SanitizeMol, MolToSmiles
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
from beaker.utils.functional import _call, _apply
from beaker.utils.io import _parseMolData

MAX_PASSES = 10
CBL_DESC_LIST = ['qed', 'MolWt', 'ExactMolWt', 'TPSA', 'HeavyAtomCount', 'NumHAcceptors', 'NumHDonors', 
                 'NumRotatableBonds', 'MolLogP', 'NumAromaticRings']

# ----------------------------------------------------------------------------------------------------------------------


def _apply_rxn(mol, rxn):
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


def _neutralise_sulphoxide(mol):
    smirks = '[S+1:1][O-1:2]>>[S+0:1]=[O-0:2]'
    rxn = rdChemReactions.ReactionFromSmarts(smirks)
    frags = rdmolops.GetMolFrags(mol, asMols=True)
    n_frags = list(filter(lambda x: x is not None,
                          [_apply_rxn(frag, rxn) for frag in frags]))
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


def _remove_isotope_info(mol):
    """
    Remove isotpe information so ExactMolWt will use the most abundant isotope weight for each atom and hence calculate 
    the monoisotopic mass.
    """
    for atom in mol.GetAtoms():
        atom.SetIsotope(0)
    return mol

# ----------------------------------------------------------------------------------------------------------------------


def _ro3_pass(mw_freebase, hba, hbd, alogp, rtb, psa):
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

def _num_ro5_violations(alogp, mw_freebase, hba, hbd):
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

def _chembl_desc_list(mol):
    mol = _remove_isotope_info(mol)
    mol = _neutralise_sulphoxide(mol)
    descriptors = dict()
    for name, fn in Descriptors.descList:
        if name in CBL_DESC_LIST:
            descriptors[name] = fn(mol)
    if 'MolecularFormula' not in descriptors:
        descriptors['MolecularFormula'] = CalcMolFormula(mol)
    descriptors['Ro3Pass'] = _ro3_pass(descriptors['MolWt'],
                                       descriptors['NumHAcceptors'],
                                       descriptors['NumHDonors'],
                                       descriptors['MolLogP'],
                                       descriptors['NumRotatableBonds'],
                                       descriptors['TPSA'])
    descriptors['NumRo5'] = _num_ro5_violations(descriptors['MolLogP'],
                                                descriptors['MolWt'],
                                                descriptors['NumHAcceptors'],
                                                descriptors['NumHDonors'])
    descriptors['MonoisotopicMolWt'] = descriptors.pop('ExactMolWt')
    return descriptors

# ----------------------------------------------------------------------------------------------------------------------

def _desc_list(mol, names):
    descriptors = dict()
    for name, fn in Descriptors.descList:
        if not names or name in names:
            descriptors[name] = fn(mol)
    if 'MolecularFormula' not in descriptors:
        descriptors['MolecularFormula'] = CalcMolFormula(mol)
    return descriptors

# ----------------------------------------------------------------------------------------------------------------------

def _getDescriptors(data, ds, loadMol=True, useRDKitChemistry=True):
    if ds:
        ds = ds.split(',')
    return _apply(_parseMolData(data, loadMol=loadMol, useRDKitChemistry=useRDKitChemistry), _desc_list, ds)
# ----------------------------------------------------------------------------------------------------------------------

def _getChemblDescriptors(data, loadMol=True, useRDKitChemistry=True):
    return _apply(_parseMolData(data, loadMol=loadMol, useRDKitChemistry=useRDKitChemistry), _chembl_desc_list)
