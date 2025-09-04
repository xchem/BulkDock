import mrich
import logging
from pathlib import Path

mrich.debug("from fragmenstein import Wictor, Laboratory")
from fragmenstein import Wictor, Laboratory  # , Igor

mrich.debug("from fragmenstein.laboratory.validator import place_input_validator")
from fragmenstein.laboratory.validator import place_input_validator

from pandas import DataFrame
from .io import mols_to_sdf
from rdkit import Chem


def fragmenstein_place(
    *,
    animal: "HIPPO",
    scratch_dir: "Path",
    smiles: str,
    inchikey: str,
    reference: "Pose",
    inspirations: "list[Pose]",
    protein_path: "Path",
    writer: "SDWriter",
    n_cores: int = 8,
    n_retries: int = 3,
    timeout: int = 300,
    write_hit_mols: bool = True,
    metadata: dict | None = None,
) -> "Pose | bool":

    metadata = metadata or {}

    # set up lab
    laboratory = setup_wictor_laboratory(
        scratch_dir=scratch_dir, protein_path=protein_path
    )

    # create inputs
    queries = create_fragmenstein_queries_df(
        smiles=smiles, inchikey=inchikey, reference=reference, inspirations=inspirations
    )

    orig_smiles = smiles

    # validate inputs
    queries = place_input_validator(queries)

    name = queries.at[0, "name"]
    smiles = queries.at[0, "smiles"]
    subdir = scratch_dir / name

    mrich.h3("Fragmenstein info")
    mrich.var("name", name)
    mrich.var("orig_smiles", orig_smiles)
    mrich.var("smiles", smiles)
    mrich.var("scratch_dir", scratch_dir)
    mrich.var("subdir", subdir)
    mrich.var("protein_path", protein_path)

    for attempt in range(n_retries):

        # run the placement
        result = laboratory.place(
            queries,
            n_cores=n_cores,
            timeout=timeout,
        )

        # process outputs

        if result is None:
            mrich.error("Placement null result")
            continue

        assert len(result) == 1

        result = result.iloc[0].to_dict()

        if result["outcome"] == "crashed" and result["error"] == "TimeoutError":
            mrich.error("Placement timed out")
            continue

        mrich.h3("Placement Result")

        mrich.var("name", result.get("name", "N/A"))
        mrich.var("error", result.get("error", "N/A"))
        mrich.var("mode", result.get("mode", "N/A"))
        mrich.var("∆∆G", result.get("∆∆G", "N/A"))
        mrich.var("comRMSD", result.get("comRMSD", "N/A"))
        mrich.var("runtime", result.get("runtime", "N/A"))
        mrich.var("outcome", result.get("outcome", "N/A"))

        # write some mols to files for debugging
        if write_hit_mols and "hit_mols" in result:
            mols_to_sdf(result["hit_mols"], subdir / "hit_mols.sdf")

        break

    ## Write data into SDF

    mol_path = subdir / f"{name}.minimised.mol"

    if mol_path.exists():

        mol_bytes = result.get("min_binary")

        if not mol_bytes:
            mrich.error("no 'min_binary' in result")
            return False

        mol = Chem.Mol(mol_bytes)

        mol.SetProp("_Name", name)
        mol.SetProp("smiles", smiles)
        mol.SetProp("orig_smiles", orig_smiles)
        mol.SetProp("inchikey", inchikey)
        mol.SetProp("target_id", str(1))
        mol.SetProp("reference_id", str(reference.id))
        mol.SetProp("inspiration_ids", str([p.id for p in inspirations]))
        mol.SetProp("energy_score", str(result.get("∆∆G", "N/A")))
        mol.SetProp("distance_score", str(result.get("comRMSD", "N/A")))
        mol.SetProp("path", str(mol_path))
        mol.SetProp("scratch_subdir", str(subdir.resolve()))
        mol.SetProp("fragmenstein_runtime", str(result.get("runtime", "N/A")))
        mol.SetProp("fragmenstein_outcome", str(result.get("outcome", "N/A")))
        mol.SetProp("fragmenstein_mode", str(result.get("mode", "N/A")))
        mol.SetProp("fragmenstein_error", str(result.get("error", "N/A")))

        writer.write(mol)

        mrich.success(f"Wrote data to SDF")

        return True

    else:

        mrich.error("Placement not successful")

        return False


# from syndirella.slipper.SlipperFitter.setup_Fragmenstein
def setup_wictor_laboratory(
    *,
    scratch_dir: "Path",
    protein_path: "Path",
    monster_joining_cutoff: float = 5,  # Å
) -> "Laboratory":

    # from fragmenstein import Laboratory, Wictor, Igor
    assert Path(scratch_dir).exists(), f"{scratch_dir=} does not exist"
    assert Path(protein_path).exists(), f"{protein_path=} does not exist"

    # set up Wictor
    Wictor.work_path = scratch_dir
    Wictor.monster_throw_on_discard = True  # stop if fragment unusable
    Wictor.monster_joining_cutoff = monster_joining_cutoff
    Wictor.quick_reanimation = False  # for the impatient
    Wictor.error_to_catch = Exception  # stop the whole laboratory otherwise
    Wictor.enable_stdout(logging.CRITICAL)
    Wictor.enable_logfile(scratch_dir / f"fragmenstein.log", logging.DEBUG)

    # os.chdir(output_path)  # needed?

    Laboratory.Victor = Wictor

    # Igor.init_pyrosetta() # needed?

    with open(protein_path) as fh:
        pdbblock: str = fh.read()

    lab = Laboratory(pdbblock=pdbblock, covalent_resi=None, run_plip=False)

    return lab


def create_fragmenstein_queries_df(
    *, smiles: str, inchikey: str, reference: "Pose", inspirations: "list[Pose]"
):

    return DataFrame(
        [
            {
                "name": f"{inchikey}-{reference}",
                "smiles": smiles,
                "hits": [pose.mol for pose in inspirations],
            }
        ]
    )
