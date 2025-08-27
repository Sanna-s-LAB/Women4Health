## file updated on 27/08/2025
## by Andrea Carta andrea.carta88@unica.it and Fabio Chillotti fabiochillotti@cnr.it 

import pandas as pd
from bs4 import BeautifulSoup
import zipfile
import os
 
def process_html_to_zip(file_path, prefix):
    """
    Extract HTML tables from a file, save on CSV and compress in one ZIP file.
 
    Args:
        file_path (str): path of HTML file.
        prefix (str): prefix to use for CSV and ZIP.
 
    Returns:
        str: path of generated  ZIP file
    """
    # Read HTML content
    with open(file_path, "r", encoding="utf-8") as file:
        soup = BeautifulSoup(file, "html.parser")
 
    # Extract tables
    tables = soup.find_all("table")
    dataframes = []
    for table in tables:
        rows = table.find_all("tr")
        table_data = []
        for row in rows:
            cols = row.find_all(["td", "th"])
            cols = [ele.get_text(strip=True) for ele in cols]
            table_data.append(cols)
        df = pd.DataFrame(table_data)
        dataframes.append(df)
 
    # Create ZIP file
    zip_filename = f"/mnt/data/tables_{prefix}.zip"
    csv_filenames = []
 
    for i, df in enumerate(dataframes):
        csv_filename = f"/mnt/data/tables_{prefix}_table_{i+1}.csv"
        df.to_csv(csv_filename, index=False, header=False)
        csv_filenames.append(csv_filename)
 
    with zipfile.ZipFile(zip_filename, "w") as zipf:
        for csv_file in csv_filenames:
            zipf.write(csv_file, os.path.basename(csv_file))
 
    return zip_filename