def printLabelNodes(path=None, df=None):
    with open(path, 'w') as f:
        f.write("id,label\n")
        for i in range(0, len(df)):
            f.write(str(i)+','+str(df[i])+'\n')

def printNodes(path=None, df=None):
    with open(path, 'w') as f:
        f.write("id,label,%s\n" % ",".join(df.keys()))
        count = 0
        for ind, ser in df.iterrows():
            f.write(str(count)+','+ind+','+','.join(str(s) for s in ser)+'\n')
            count += 1

def printCorr(path=None, df=None, a=None, b=None, c=None):
    with open(path, 'w') as f:
        f.write("source,target,type,weight,a,b,c,signal\n")
        for i in range(0, df.shape[0]):
            for j in range(i+1, df.shape[1]):
                f.write("%d,%d,Undirected,%.8lf,%d,%d,%d,%s\n" %
                    (i, j, df[i,j], a[i][j], b[i][j], c[i][j],'pos' if df[i,j] >= 0. else "neg"))
