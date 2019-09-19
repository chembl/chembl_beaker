__author__ = "Juan F. Mosquera"
__email__ = "jfmosquera@ebi.ac.uk"


class UnexpectedException(Exception):
    """
        Last resort exception for unexpected cases
    """
    pass


class ChemistryError(Exception):
    """
        Error that can occur while handling chemistry data
    """
    pass
