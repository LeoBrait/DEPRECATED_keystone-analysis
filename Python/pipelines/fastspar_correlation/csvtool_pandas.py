import sys
import pandas as pd

f = sys.argv[1]
o = sys.argv[2]

x = pd.read_csv(f, index_col=0)

# fastspar implementation requirements (BIOME tsv format)
x.index.name = "#OTU ID"
x.columns.name = "#OTU ID"

x.T.to_csv(o,sep='\t')
