
class ClickException(Exception):
    """
    Base exception for the Click package.he
    """
    pass


class NoProductError(ClickException):
    """
    Error raised if a reaction does not give any product.
    """
    pass


class AmbiguousProductError(ClickException):
    """
    Gets raised if a reaction leads unexpectedly to one or more products.
    """
    pass
