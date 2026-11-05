import pandas as pd
from pathlib import Path

# ============================================================
# Settings
# ============================================================

# Folder containing the combined CSV files
input_dir = Path("map_stats_per_sample/combined_results")

# Output Excel file
output_file = input_dir / "mapping_stats_bacterial_mapping.xlsx"


# ============================================================
# Find all combined CSV files
# ============================================================

csv_files = sorted(
    input_dir.glob("*_combined.csv")
)

print("Found", len(csv_files), "combined CSV files")


# ============================================================
# Create Excel workbook
# ============================================================

with pd.ExcelWriter(
    output_file,
    engine="openpyxl"
) as writer:

    for csv_file in csv_files:

        # Extract sample ID
        # Example:
        # MS11669_combined.csv
        # becomes:
        # MS11669
        sample = csv_file.stem.replace(
            "_combined",
            ""
        )

        print("Adding:", sample)

        # Read combined CSV
        df = pd.read_csv(csv_file)

        # Excel worksheet names can be maximum 31 characters
        sheet_name = sample[:31]

        # Write to Excel
        df.to_excel(
            writer,
            sheet_name=sheet_name,
            index=False
        )


# ============================================================
# Finished
# ============================================================

print()
print("=" * 60)
print("Finished!")
print("=" * 60)
print()
print("Excel file created:")
print(output_file)
