class MultiscaleError(RuntimeError):
    """Base error for Multiscale workflow failures."""


class MissingExternalProjectError(MultiscaleError):
    """Raised when a required external project is unavailable."""
