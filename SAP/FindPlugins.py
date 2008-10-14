
pkg_resources_available = True
try:
    from pkg_resources import load_entry_point, get_entry_info, get_entry_map, iter_entry_points, EntryPoint
except ImportError:
    pkg_resources_available = False

# Make sure the try statement does not hide the import from py2app


class Error(Exception):
    """
    Base class for exceptions in this module.
    """
    pass

    
class PluginNotFoundError(Error):
    """
    Raised when a plugin is not found.
    """
    def __init__(self, plugin):
        self.plugin = plugin

def findPlugin(name, entry_point_name):

    if not pkg_resources_available:
        raise PluginNotFoundError, name
    
    plugin = None
    for entryp in iter_entry_points(entry_point_name):
       if entryp.name == name:
          plugin = entryp.load()
          break
    if not plugin:
       raise PluginNotFoundError, name
    return plugin
