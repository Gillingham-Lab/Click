
class ClickException(Exception):
    pass


class NoProductError(ClickException):
    pass


class AmbiguousProductError(ClickException):
    pass
