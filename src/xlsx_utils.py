"""Shared helper for the Extended Data Table scripts (scripts/plots/FigSupp/ExtendedDataTableN.py)
- every one of them saves as .xlsx with a "data" sheet (the actual table) and a "legends" sheet
(one row per column: "Column name", "Explanation"), user request."""
import pandas as pd


def save_table_with_legend(table, column_explanations, output_path):
    """table: the data DataFrame, saved as-is to the "data" sheet. column_explanations: dict
    {column name: explanation} - must cover every column in `table` (checked here), written to
    the "legends" sheet in `table`'s own column order, not dict insertion order, so the two
    sheets stay in step even if a caller built the dict out of order."""
    missing = [c for c in table.columns if c not in column_explanations]
    if missing:
        raise ValueError(f"No legend explanation provided for column(s): {missing}")
    legend = pd.DataFrame({
        "Column name": table.columns,
        "Explanation": [column_explanations[c] for c in table.columns],
    })
    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        table.to_excel(writer, sheet_name="data", index=False)
        legend.to_excel(writer, sheet_name="legends", index=False)
