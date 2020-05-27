import pandas as pd
import numpy as np


def process2genes(path):
    fp = open(path, 'r')  
    path2genes={}
    for line in fp.readlines():
        for entry in line.split("\t"):
            if ("GO:" in entry):
                goid = entry.split("(GO:")[1].split(")")[0]
                goid = ("GO:"+goid)
                path2genes[goid] = []
                continue  
            else:
                gene = entry.split(",")[0].rstrip()
                if len(gene)>0:
                    path2genes[goid].append(gene)
    return(path2genes)

def constructConnectMat(path2genes,df):
    #first get all genes
    glist = []
    relevant_paths = df["term_ID"].values
    clusters = df["cluster"].values
    for path in relevant_paths:
        glist+=path2genes[path]
    glist = list(set(glist))
    print(len(glist))

    binarr=[];paths=[]
    for path in relevant_paths:
        subarr=[]
        for gene in glist:
            if (gene in path2genes[path]):
                subarr.append(1)
            else:
                subarr.append(0)
        binarr.append(subarr)
        paths.append(path)

    dfout = pd.DataFrame(binarr,columns=glist)
    dfout["cluster"] = clusters 
    dfout.insert(0,"go annotation",paths)
    return(dfout)
     
def insert_cluster_assign(df):
    cntr=0; assigns=[0]
    vals = df['eliminated'].values
    
    for i in range(1,vals.shape[0]):
        current = vals[i-1]
        next = vals[i]
        if ((current==0) & (next==0)):
            cntr+=1
            assigns.append(cntr);continue
        if ((current==0) & (next==1)):
            assigns.append(cntr);continue
        if ((current==1) & (next==1)):
            assigns.append(cntr);continue
        if ((current==1) & (next==0)):
            cntr+=1
            assigns.append(cntr)
    df["cluster"] = assigns
    return(df)

def return_parent_terms(df):
    clusters = np.unique(df["cluster"].values)
    assoc = {};parent_terms = []
    for cluster in clusters:
        subdf = df[df["cluster"]==cluster]
        parent_term = subdf["term_ID"].values[0]
        parent_terms.append(parent_term)
    return(parent_terms)


path2genes = process2genes("references/GO_Biological_Process_2018.txt")
pathways_df = pd.read_csv("references/REVIGO_1lfc_short.csv")


pathways_df = insert_cluster_assign(pathways_df)
parent_terms = return_parent_terms(pathways_df)
mat = constructConnectMat(path2genes,pathways_df)


mat = mat.groupby("cluster").sum()
mat.insert(0,"go annotation",parent_terms)
matvals = mat[mat.columns[1:]].values
matvals[matvals>1] = 1
mat[mat.columns[1:]] = matvals
#replace go id with go description
id2desc = dict(zip(pathways_df["term_ID"].values,pathways_df["description"].values))
newcols = []

for pid in mat["go annotation"].values:
    if (id2desc[pid]+" "+pid == "response to cesium ion GO:0010164"):
        newcols.append("response to metal ion GO:0010038")
    else:
        newcols.append(id2desc[pid]+" "+pid)

mat["go annotation"] = newcols
mat.to_csv("results/revigo_graph.csv",index=False)



