import sys
from pathlib import Path
path_to_dev = (Path(__file__).parent.parent.parent.parent / 'dev').absolute()
sys.path.append(str(path_to_dev))
import generate_headers 
__version__ = lambda : generate_headers.get_version(root_dir=str(Path(__file__).parent.parent.parent.parent))