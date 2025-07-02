**Setting Up a Hail Environment on Sherlock with Python 3.12**

This document outlines the exact lines used to successfully set up and run Hail on the Sherlock computing cluster using Python 3.12.

---

### 1. Load Required Modules

*Note*: Python 3.9.0 is not used due to a known bug in the standard `typing` module that causes Hail to fail with `TypeError: unhashable type: 'list'` during initialization.

```bash
module load python/3.12.1 java/11.0.11 gcc/10.1.0
module load zlib libjpeg-turbo freetype libpng
```

*Notes*: The second line of modules above are prerequisites for installing Pillow, one of the dependencies of Hail. These only need to be loaded prior to installing Hail.

---

### 2. Create and Activate a Virtual Environment

```bash
python -m venv ~/hail_env_py312
source ~/hail_env_py312/bin/activate
pip install --upgrade pip
```

---

### 3. Install Hail and Dependencies

```bash
pip install hail
pip install pyspark==3.5.0
```

*Notes*: Avoid installing pyspark 4.x, which is incompatible with Hail.

---

---

### 5. Optional: Test Hail

```bash
#!/bin/bash
#SBATCH --job-name=hail_test
#SBATCH --partition=mrivas
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --time=01:00:00
#SBATCH --output=hail_test.log
#SBATCH --error=hail_test.err

ml python/3.12.1
ml java/11.0.11
ml gcc/10.1.0
source ~ls

python -c "
import sys
import hail as hl

try:
    hl.init()
    print('  Hail initialized successfully.')
    print(hl.__version__)
except Exception as e:
    print('  Hail failed to initialize:', e, file=sys.stderr)
    sys.exit(1)  # Exit with error to stop the SLURM job
"
```

*Notes*: Avoid running Hail on the login node. Always use an `salloc` or SLURM script for initialization. Errors like `JAVA_GATEWAY_EXITED` are common on login nodes.

---

With this streamlined process, you can reproduce a working Hail setup on Sherlock.
