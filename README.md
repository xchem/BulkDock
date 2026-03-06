# BulkDock
💪 BulkDock: Manage batches of Fragmenstein restrained protein-ligand docking jobs

## Usage

Clone the repo on IRIS and change to directory, picking a suitable release.
Here we clone release 1.0.0 of the repository: -

```bash
git clone --depth 1 --branch 1.0.0 git@github.com:xchem/BulkDock
cd BulkDock
```

### Configure

Configure variables

```
python -m bulkdock configure
```

### Create directories

```
python -m bulkdock create-directories
```

### Setup target

First download a target zip archive from fragalysis and place in TARGETS directory.

Extract the target:

```
python -m bulkdock extract TARGET_NAME
```

Setup the HIPPO database **from within a valid SLURM job**:

```
python -m bulkdock setup TARGET_NAME
```

On IRIS this can be achieved with:

```
sb.sh --job-name "BULKDOCK_SETUP" /opt/xchem-fragalysis-2/maxwin/slurm/run_python.sh -m bulkdock setup TARGET_NAME
```

### Placements

Place a CSV in the INPUTS directory with the following structure:

|           smiles            | inspiration 1 | inspiration 2 |
|-----------------------------|---------------|---------------|
| CC(=O)Nc1ccccc1CCS(N)(=O)=O | A0487a        | A0719a        |

**The `smiles` column is required and must come first. The name of the other columns is not used, but values in those columns must refer to the observation shortcode of any inspiration hits for this compound**

Configure the `SLURM_PYTHON_SCRIPT` variable which should point to a template SLURM file that has some basic SLURM headers, sets up the necessary python environment and then executes `python @`. On IRIS this is provided at:

```
python -m bulkdock configure SLURM_PYTHON_SCRIPT /opt/xchem-fragalysis-2/maxwin/slurm/run_python.sh
```

Optionally set the variable `DIR_SLURM_LOGS` to put the log files in a place other than `SCRATCH/logs`.


Once configured, and the inputs have been placed in the INPUTS directory, the command to launch the placement jobs is:

```
python -m bulkdock place TARGET_NAME SDF_NAME
```

Once the placement jobs have finished the individual SDF outputs will be located in the OUTPUTS directory as configured. The above command will also queue a `combine` job to run after the placement jobs, and generate a `_combined.sdf` output.

### Monitoring jobs

To monitor the jobs try:

```
python -m bulkdock status
```

### Collating outputs from multiple (failed) placement jobs

If a placement jobs do not correctly write out SDF outputs you can use a "collate" job to extract any Poses from certain jobs that were registered to the database. In this case write a json file containing the job ID's to the `SCRATCH/${TARGET}_inputs` directory containing the job ids. E.g. with python:

```
from bulkdock import BulkDock
import json
job_ids = {1,2,3,4}
engine = BulkDock()
target = "FatA"
json.dump(list(job_ids), open(engine.get_scratch_subdir(f"{target}_inputs") / "FatA_Knitwork_36_active_collate_job_ids.json", "wt"))
```

And submit a collation job:

```
sbatch --job-name "BulkDock.collate:FatA:FatA_Knitwork_36_active" ../slurm/run_python.sh -m bulkdock.batch collate "FatA_Knitwork_36_active.sdf" FatA SCRATCH/FatA_inputs/FatA_Knitwork_36_active_collate_job_ids.json
```

### Fragalysis export

The SDF output can be modified for direct upload to the Fragalysis RHS with the `to-fragalysis` command. To see the options:

```
python -m bulkdock to-fragalysis --help
```

To avoid passing the same values for submitter_institution, submitter_name, submitter_email, and ref_url you can set the config variables:

```
python -m bulkdock configure FRAGALYSIS_EXPORT_SUBMITTER_NAME "Max Winokan"
python -m bulkdock configure FRAGALYSIS_EXPORT_SUBMITTER_EMAIL "max.winokan@diamond.ac.uk"
python -m bulkdock configure FRAGALYSIS_EXPORT_SUBMITTER_INSTITUTION "DLS"
python -m bulkdock configure FRAGALYSIS_EXPORT_REF_URL "https://github.com/mwinokan/BulkDock"
```

Then to check the configuration:

```
python -m bulkdock to-fragalysis TARGET SDF_FILE METHOD_NAME
```

Once happy, submit the job (or run the above from within a notebook):

```
sb.sh --job-name "BULKDOCK_EXPORT" /opt/xchem-fragalysis-2/maxwin/slurm/run_python.sh -m bulkdock to-fragalysis TARGET SDF_FILE METHOD_NAME
```
