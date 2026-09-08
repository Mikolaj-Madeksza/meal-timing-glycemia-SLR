"""Export the extraction sheets without changing their values or column order."""
import csv
from pathlib import Path
from openpyxl import load_workbook

ROOT = Path(__file__).resolve().parents[1]
SHEETS = {"Studies": "studies", "Reports": "reports", "Conditions": "conditions",
          "Outcomes": "outcomes", "Derivations": "derivations", "Data limitations": "data_limitations"}

def main():
    workbook = load_workbook(ROOT / "extraction.xlsx", read_only=True, data_only=False)
    destination = ROOT / "meta_analysis" / "data"
    for name, filename in SHEETS.items():
        sheet = workbook[name]
        rows = list(sheet.values)
        if any(cell.data_type == "f" for row in sheet for cell in row):
            raise ValueError(f"{name}: export requires source values, not unevaluated formulas")
        with (destination / f"{filename}.csv").open("w", newline="", encoding="utf-8") as handle:
            csv.writer(handle, lineterminator="\n").writerows(rows)
        print(f"{name}: {len(rows)-1} records")
    workbook.close()

if __name__ == "__main__":
    main()
