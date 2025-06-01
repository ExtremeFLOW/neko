import os

def get_neko():
    """
    Returns the path to the turboneko executable via the
    environmental variable NEKO_BIN, or just returns "neko" 
    """

    src_dir = os.getenv("NEKO_BIN")
    if src_dir:
        return os.path.join(src_dir, "neko")
    else:
        return "neko"
