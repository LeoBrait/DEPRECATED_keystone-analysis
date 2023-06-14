import pandas as pd

# Read the CSV file into a pandas DataFrame
df = pd.read_csv('phyla.relative.matrix.csv')

# Get unique habitat values
habitats = df['habitat'].unique()

# Iterate over each unique habitat
for habitat in habitats:
    # Create a subset DataFrame based on the current habitat
    subset = df[df['habitat'] == habitat]
    
    # Get unique ecosystem values in the current habitat subset
    ecosystems = subset['ecosystem'].unique()
    
    # Iterate over each unique ecosystem in the current habitat
    for ecosystem in ecosystems:
        # Create a sub-subset DataFrame based on the current ecosystem
        sub_subset = subset[subset['ecosystem'] == ecosystem]
        
        # Generate the filename based on habitat and ecosystem names
        filename = f"{habitat}_{ecosystem}.csv"
        
        # Save the sub-subset DataFrame to a new CSV file
        sub_subset.to_csv(filename, index=False)
