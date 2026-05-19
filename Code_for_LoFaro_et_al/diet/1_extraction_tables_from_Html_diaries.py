# ANDREA CARTA 15/04/2026

# ============================================
# Generic HTML Tables Extraction Script (IDXXX)
# ============================================

import pandas as pd
from bs4 import BeautifulSoup
import zipfile
import os

# --------------------------------------------
# USER INPUT
# --------------------------------------------
# Replace IDXXX.html with your actual file (e.g., ID129.html)
html_file_path = "IDXXX.html"

# Output folder
output_dir = "output_tables"

# --------------------------------------------
# SETUP
# --------------------------------------------
# Create output directory if it does not exist
os.makedirs(output_dir, exist_ok=True)

# Extract prefix from file name (e.g., IDXXX)
prefix = os.path.splitext(os.path.basename(html_file_path))[0]

# --------------------------------------------
# LOAD AND PARSE HTML
# --------------------------------------------
with open(html_file_path, "r", encoding="utf-8") as file:
    soup = BeautifulSoup(file, "html.parser")

# Extract all tables
tables = soup.find_all("table")

print(f"Found {len(tables)} tables in {html_file_path}")

# --------------------------------------------
# EXTRACT TABLES AND SAVE AS CSV
# --------------------------------------------
csv_files = []

for i, table in enumerate(tables, start=1):
    rows = table.find_all("tr")
    table_data = []

    for row in rows:
        cells = row.find_all(["td", "th"])
        cell_text = [cell.get_text(strip=True) for cell in cells]
        table_data.append(cell_text)

    df = pd.DataFrame(table_data)

    # Save CSV
    csv_filename = f"{prefix}_table_{i}.csv"
    csv_path = os.path.join(output_dir, csv_filename)
    df.to_csv(csv_path, index=False, header=False, encoding="utf-8-sig")

    csv_files.append(csv_path)

    print(f"Saved: {csv_filename}")

# --------------------------------------------
# CREATE ZIP ARCHIVE
# --------------------------------------------
zip_filename = f"{prefix}_tables.zip"
zip_path = os.path.join(output_dir, zip_filename)

with zipfile.ZipFile(zip_path, "w", zipfile.ZIP_DEFLATED) as zipf:
    for file in csv_files:
        zipf.write(file, os.path.basename(file))

print(f"\nZIP archive created: {zip_path}")

# --------------------------------------------
# DONE
# --------------------------------------------
