import time

# ----------------------------------------------------------------
def get_start_stamp():
    """
    Returns information about start time of an event. Nominally float seconds since the epoch,
    but articulated here as being compatible with the format_elapsed function.
    """
    return time.time()

def format_elapsed(start_stamp, message: str):
    return "%s TIME %.3f" % (message, time.time() - start_stamp)
