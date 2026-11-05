import pandas as pd
from pathlib import Path

# ============================================================
# Settings
# ============================================================

# Folder containing all sample files
input_dir = Path("map_stats_per_sample")

# Folder where combined files will be written
output_dir = input_dir / "combined_results"
output_dir.mkdir(exist_ok=True)


# ============================================================
# Process each sample
# ============================================================

# Find all pydamage result files
for csv_file in input_dir.glob("*_all_bacteria_pydamage_results.csv"):

    # Extract sample ID
    # Example:
    # MS11669_all_bacteria_pydamage_results.csv
    # becomes:
    # MS11669
    sample = csv_file.name.replace(
        "_all_bacteria_pydamage_results.csv",
        ""
    )

    print()
    print("=" * 60)
    print("Processing:", sample)
    print("=" * 60)

    # Corresponding mapping stats file
    mapping_file = input_dir / f"{sample}_mapping_stats.txt"

    if not mapping_file.exists():
        print("WARNING: Mapping stats file not found:")
        print(mapping_file)
        continue


    # ========================================================
    # Read pydamage results
    # ========================================================

    damage = pd.read_csv(csv_file)

    # Clean column names
    damage.columns = [
        str(col).strip()
        for col in damage.columns
    ]

    print("Damage columns:")
    print(damage.columns.tolist())

    # Check that reference exists
    if "reference" not in damage.columns:
        print("ERROR: 'reference' column not found in pydamage file")
        print("Columns found:", damage.columns.tolist())
        continue


    # ========================================================
    # Read mapping stats
    # ========================================================

    mapping = pd.read_csv(
        mapping_file,
        sep=r"\s+",
        header=0,
        comment=None,
        engine="python"
    )

    # Clean column names
    mapping.columns = [
        str(col).strip().lstrip("#").strip()
        for col in mapping.columns
    ]

    # Explicitly rename the first column to "reference"
    # This handles #rname, rname, or other variations.
    mapping = mapping.rename(
        columns={
            mapping.columns[0]: "reference"
        }
    )

    print("Mapping columns:")
    print(mapping.columns.tolist())


    # ========================================================
    # Check reference columns
    # ========================================================

    if "reference" not in mapping.columns:
        print("ERROR: 'reference' column not found in mapping file")
        print("Columns found:", mapping.columns.tolist())
        continue

    print("Damage reference exists:", "reference" in damage.columns)
    print("Mapping reference exists:", "reference" in mapping.columns)


    # ========================================================
    # Clean reference names
    # ========================================================

    # Remove leading/trailing whitespace from reference names
    damage["reference"] = damage["reference"].astype(str).str.strip()
    mapping["reference"] = mapping["reference"].astype(str).str.strip()


    # ========================================================
    # Merge the two tables
    # ========================================================

    combined = pd.merge(
        mapping,
        damage,
        on="reference",
        how="outer",
        suffixes=("_mapping", "_damage")
    )


    # ========================================================
    # Report matching statistics
    # ========================================================

    damage_refs = set(damage["reference"])
    mapping_refs = set(mapping["reference"])

    matched = damage_refs.intersection(mapping_refs)
    only_damage = damage_refs - mapping_refs
    only_mapping = mapping_refs - damage_refs

    print()
    print("Reference statistics:")
    print("  Pydamage references :", len(damage_refs))
    print("  Mapping references  :", len(mapping_refs))
    print("  Matched references  :", len(matched))
    print("  Only in pydamage   :", len(only_damage))
    print("  Only in mapping    :", len(only_mapping))

    if only_damage:
        print()
        print("WARNING: References only in pydamage:")
        for ref in sorted(only_damage):
            print(" ", ref)

    if only_mapping:
        print()
        print("WARNING: References only in mapping:")
        for ref in sorted(only_mapping):
            print(" ", ref)


    # ========================================================
    # Save combined table
    # ========================================================

    output_file = output_dir / f"{sample}_combined.csv"

    combined.to_csv(
        output_file,
        index=False
    )

    print()
    print("Combined table saved to:")
    print(output_file)
    print("Number of rows:", len(combined))
    print("Number of columns:", len(combined.columns))


print()
print("=" * 60)
print("Finished.")
print("=" * 60)