import sys
import os

def load_simpl_module():
    build_path = os.path.abspath("../build/Release")
    if build_path not in sys.path:
        sys.path.append(build_path)
    import simpl
    return simpl
