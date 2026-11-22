from astropy.table import Table
import re

input_path = 'import path.csv'
output_path = 'export path ipac.txt'

data = Table.read(input_path, format='ascii.csv')

# Automatically clean every column name to match IPAC rules
for i, col in enumerate(data.colnames):
    # Replace spaces, commas, or symbols with underscores
    clean_name = re.sub(r'[^a-zA-Z0-9_]', '_', col)

    # Rename blank/unnamed columns (like 'col1' or empty strings)
    if not clean_name or clean_name.startswith('_'):
        clean_name = f'col_{i}_{clean_name}'
    data.rename_column(col, clean_name)

data.write(output_path, format='ascii.ipac', overwrite=True)
