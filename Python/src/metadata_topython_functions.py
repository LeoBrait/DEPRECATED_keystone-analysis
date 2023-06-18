import numpy as np
import pandas as pd

def metadata(nodes):
    df = pd.read_csv(nodes)
    Taxa = df["label"].tolist()
    matching = [i for i,s in enumerate(Taxa) if "Candidatus" in s]
    matching2 = [i for i,s in enumerate(Taxa) if "candidate" in s]
    N = len(Taxa)
    Candidatus = np.zeros(N)
    Candidatus[matching] =1
    Candidatus[matching2] =1

    Taxa = np.array(Taxa)
    Abundance = np.array(df["abundance (%)"].tolist())
    AbsoluteAbundance = df["abundance"].tolist()
    Prevalence = df["prevalence"].tolist()

    Marker = []
    for i in range(N):
        if (Candidatus[i]==0):
            if Prevalence[i] < .5:
                Marker.append('o')
            elif Prevalence[i] < .75:
                Marker.append('s')
            else:
                Marker.append('D')

        else:
            if Prevalence[i] < .5:
                Marker.append('v')
            elif Prevalence[i] < .75:
                Marker.append('>')
            else:
                Marker.append('^')
    Marker = np.array(Marker)

    Sem = df["sem"].tolist()
    Std = df["std"].tolist()
    Med = df["med"].tolist()

    return N, Marker, Taxa, Abundance, AbsoluteAbundance, Prevalence, Sem, Std, Med
