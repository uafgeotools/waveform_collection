# -----------------------------------------------------------------------------
# WARNING CONFIG
# -----------------------------------------------------------------------------

import warnings

# Subclass UserWarning
class CollectionWarning(UserWarning):
    pass

# Make warnings more consistent
warnings.simplefilter(action='always', category=CollectionWarning)

# Make a custom warning format for CollectionWarning instances
def _formatwarning(message, category, *args, **kwargs):
    if category == CollectionWarning:
        msg = f'{category.__name__}: {message}\n'  # Much cleaner
    else:
        import warnings
        msg_form = warnings.WarningMessage(message, category, *args, **kwargs)
        msg = warnings._formatwarnmsg_impl(msg_form)  # Default
    return msg
warnings.formatwarning = _formatwarning

# Clean up
del _formatwarning
del warnings

# -----------------------------------------------------------------------------
# EXPOSE PUBLIC MODULES
# -----------------------------------------------------------------------------

from .server import gather_waveforms, gather_waveforms_bulk, INFRASOUND_CHANNELS
from . import local