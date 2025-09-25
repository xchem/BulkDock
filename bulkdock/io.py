import mrich


def split_input_csv(in_path: "Path", split: int, out_dir: "Path") -> "list[Path]":

    mrich.h3("bulkdock.io.split_input_csv()")
    mrich.var("input", in_path)
    mrich.var("output", out_dir)
    mrich.var("batch size", split)

    from pandas import read_csv

    # read dataframe from CSV
    df = read_csv(in_path)

    mrich.var("#compounds", len(df))

    dfs = [df[i : i + split] for i in range(0, df.shape[0], split)]

    mrich.var("#batches", len(dfs))

    paths = []

    for i, df in enumerate(dfs):

        out_file = in_path.name.removesuffix(".csv") + f"_split{split}_batch{i:03}.csv"

        out_path = out_dir / out_file

        mrich.writing(out_path)
        df.to_csv(out_path, index=False)
        paths.append(out_path)

    return paths


def parse_input_csv(
    animal: "HIPPO",
    file: "Path",
    debug: bool = False,
    reference: str | None = None,
) -> list[dict]:
    """
    Parse a BulkDock input CSV to prepare for an ensemble docking run where a compound is placed against each protein conformation from its inspirations.

    :param animal: `HIPPO` object to work within
    :param file: `Path` object to the input CSV
    :param debug: Increase verbosity of CLI output
    :returns: a list of dictionaries containing HIPPO objects:

    compound: Compound
    reference: Pose
    inspirations: PoseSet

    """

    from pandas import read_csv
    from hippo.tools import inchikey_from_smiles

    # read dataframe from CSV
    df = read_csv(file)

    assert "smiles" in df.columns
    
    mrich.var("len(df)", len(df))

    hits = animal.poses(tag="hits")
    mrich.var("hits", len(hits), unit="Poses")

    mrich.debug("get_pose_alias_id_dict()")
    alias_lookup = animal.db.get_pose_alias_id_dict(hits)
    mrich.debug("get_pose_id_obj_dict()")
    pose_lookup = animal.db.get_pose_id_obj_dict(hits)

    mrich.h1("Enumerate placement tasks")
    
    data = []

    for i, row in df.iterrows():

        try:
            # extract row values
            smiles = row.smiles
            inchikey = inchikey_from_smiles(smiles)
            inspirations = row.values[1:]

            # debug output
            if debug:
                mrich.debug("i", i)
                mrich.debug("smiles", smiles)

            inspirations = [i for i in inspirations if isinstance(i, str) and i]

            inspiration_poses = list()
            for alias in inspirations:
                pose_id = alias_lookup.get(alias, None)
                if not pose_id:
                    mrich.error(f"Could not find inspiration with {alias=}")
                    continue
                pose = pose_lookup[pose_id]
                inspiration_poses.append(pose)

            if not inspiration_poses:
                mrich.error(f"No inspirations. {i=}, {smiles=}")
                continue

            if not reference:

                # one placement against each inspiration's protein conformation
                for pose in inspiration_poses:

                    # all info needed for placement
                    data.append(
                        dict(
                            smiles=smiles,
                            inchikey=inchikey,
                            reference=pose,
                            inspirations=inspiration_poses,
                        )
                    )

                    # debug output
                    if debug:
                        mrich.h3("Placement")
                        mrich.var("smiles", smiles)
                        mrich.var("inchikey", inchikey)
                        mrich.var("protein", pose.alias)
                        mrich.var("inspirations", [p.alias for p in inspiration_poses])

            else:

                reference_pose = animal.poses[reference]

                # all info needed for placement
                data.append(
                    dict(
                        smiles=smiles,
                        inchikey=inchikey,
                        reference=reference_pose,
                        inspirations=inspiration_poses,
                    )
                )

                # debug output
                if debug:
                    mrich.h3("Placement")
                    mrich.var("smiles", smiles)
                    mrich.var("inchikey", inchikey)
                    mrich.var("protein", pose.alias)
                    mrich.var("inspirations", [p.alias for p in inspiration_poses])
        except Exception as e:
            mrich.error("Could not parse row", e)
            print(row)
            continue

    return data


def mols_to_sdf(mols, out_path):

    from rdkit.Chem import Mol
    from rdkit.Chem import PandasTools
    from pandas import DataFrame

    data = []
    for i, mol in enumerate(mols):
        data.append({"_Name": f"mol{i}", "ROMol": Mol(mol)})

    df = DataFrame(data)
    PandasTools.WriteSDF(df, out_path, "ROMol", "_Name", list(df.columns))
