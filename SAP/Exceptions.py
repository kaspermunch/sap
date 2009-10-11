
class Error(Exception):
    """
    Base class for exceptions in this module.
    """
    pass
    

class AnalysisTerminated(Error):
    """
    Exception raised when an analysis is to be terminated. Used instead of sys.exit() so
    we can catch the special exception in the GUI and print a message without exiting.

    Attributes:
        exit valuee
        message
    """
    def __init__(self, exitValue, message):
        self.exitValue = exitValue
        self.msg = message
