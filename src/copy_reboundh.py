import os
import rebound
from pathlib import Path

reb_so_libpath = rebound.__libpath__
reb_h_libpath  = Path(rebound.__file__).parent / "rebound.h"

cp_h  = f'cp {reb_h_libpath} .'

os.system(cp_h)