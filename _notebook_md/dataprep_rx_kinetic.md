# dataprep_rx_kinetic

_Source notebook: `dataprep_rx_kinetic.ipynb`_

# Preparation of Experimantal Data

## Testing Google Spreadsheet Connection

```python
import gspread
import pandas as pd

# 1. Trigger the browser login using your downloaded JSON file
# It will save a second file ('authorized_user.json') so you don't have to log in every single time.
gc = gspread.oauth(
    credentials_filename='credentials.json',
    authorized_user_filename='authorized_user.json'
)

# 2. Open the Google Sheet using its URL
# Replace this with the actual URL of your HPLC Google Sheet
SHEET_URL = "https://docs.google.com/spreadsheets/d/13ihUtinfoi9AmLZtnqF6BX7L-dqAk552I-KZNu2HtFc/edit#gid=0"
workbook = gc.open_by_url(SHEET_URL)

# 3. Load the Standard tab into Pandas
worksheet_name = "data" 
worksheet = workbook.worksheet(worksheet_name)

# 4. Convert the data into a Pandas DataFrame
data = worksheet.get_all_records()
df_std = pd.DataFrame(data)

# 5. Display the first 5 rows to verify!
display(df_std.head())
```

## Experimental Data Prep for Modeling

```python
import pandas as pd
import numpy as np

# 1. Define the sheets we want to process and their corresponding "Perbandingan" values
sheet_mapping = {
    "1:1": 1.0, 
    "1:1.5": 1.5, 
    "1:2": 2.0
}

# List to hold the processed dataframes before we combine them
processed_dfs = []

for sheet_name, perbandingan_val in sheet_mapping.items():
    
    # 2. Open the worksheet
    # NOTE: If your tab names have the full prefix, change this to f"Data HPLC NHS-DesB30 - {sheet_name}"
    ws = workbook.worksheet(sheet_name) 
    
    # Load into Pandas
    data = ws.get_all_records()
    df = pd.DataFrame(data)
    
    # 3. Select only the 3 specific columns
    df = df[['Kode Sampel', 'Peak Identification', 'Conc']]
    
    # 4. Filter out empty "Peak Identification" rows
    # We replace any purely blank spaces with NaN, then drop the NaNs
    df['Peak Identification'] = df['Peak Identification'].replace(r'^\s*$', np.nan, regex=True)
    df = df.dropna(subset=['Peak Identification'])
    
    # 5. Insert the "Perbandingan" column at position 0 (before "Kode Sampel")
    df.insert(0, 'Perbandingan', perbandingan_val)
    
    # 6. Ensure the data types are correct
    df['Perbandingan'] = pd.to_numeric(df['Perbandingan'])
    
    # For 'Conc', sometimes Google Sheets passes numbers as strings with commas. 
    # This ensures it forces them into strict numeric formats safely.
    df['Conc'] = pd.to_numeric(df['Conc'].astype(str).str.replace(',', ''), errors='coerce') 
    
    df['Kode Sampel'] = df['Kode Sampel'].astype(str)
    df['Peak Identification'] = df['Peak Identification'].astype(str)
    
    # Add the cleaned dataframe to our list
    processed_dfs.append(df)

# 7. Combine everything into a single DataFrame
exp_data = pd.concat(processed_dfs, ignore_index=True)

# Display the result and the data types to verify!
display(exp_data.head())
print("\n--- Data Types ---")
print(exp_data.dtypes)
```

```python
import pandas as pd
import numpy as np

# 1. Clean Time column (convert "t2" -> 2.0)
exp_data['Time'] = exp_data['Kode Sampel'].str.replace('t', '', regex=False).astype('float64')

# 2. Pivot the data WITHOUT filling missing values. 
# This finds every unique peak that ever appeared in the dataset and creates a column for it.
# Intermediate missing values remain as NaN.
exp_data_grid = exp_data.pivot_table(
    index=['Perbandingan', 'Time'], 
    columns='Peak Identification', 
    values='Conc'
).reset_index()

# 3. Melt it back into a long format
final_kinetic_data = exp_data_grid.melt(
    id_vars=['Perbandingan', 'Time'],
    value_name='Conc'
)

# 4. Tag Substrate vs. Product (Addressing Point 4)
# np.where works like an Excel IF statement: if it's desb30, label it Substrate, else Product.
final_kinetic_data['Type'] = np.where(
    final_kinetic_data['Peak Identification'].str.lower() == 'desb30', 
    'Substrate', 
    'Product'
)

# 5. The Initial Value Fix (Addressing Points 1 & 2)
# Find rows where Time is 0 AND the Type is Product, then explicitly set their Conc to 0.0
mask_t0_products = (final_kinetic_data['Time'] == 0.0) & (final_kinetic_data['Type'] == 'Product')
final_kinetic_data.loc[mask_t0_products, 'Conc'] = 0.0

# 6. Clean up intermediate missing values (Addressing Point 3)
# We drop rows that are STILL NaN (e.g., the t=40 missing data point).
final_kinetic_data = final_kinetic_data.dropna(subset=['Conc'])

# 7. Reorder columns and sort for a clean view
final_kinetic_data = final_kinetic_data[['Perbandingan', 'Time', 'Type', 'Peak Identification', 'Conc']]
final_kinetic_data = final_kinetic_data.sort_values(
    by=['Perbandingan', 'Peak Identification', 'Time']
).reset_index(drop=True)

# Display the result
display(final_kinetic_data.head(15))
```

```python
display(final_kinetic_data)
```

```python
# Save the clean data to a CSV file
file_name = 'kinetic_data_processed.csv'
final_kinetic_data.to_csv(file_name, index=False)

print(f"Dataset successfully saved as {file_name}")
```
