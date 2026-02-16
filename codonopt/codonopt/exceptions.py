class CodonoptError(Exception):
    """Base exception for codonopt"""
    pass


class ConstraintError(CodonoptError):
    """Raised when constraints make optimization impossible"""
    pass


class InputFormatError(CodonoptError):
    """Raised for invalid or unsupported input formats"""
    pass
