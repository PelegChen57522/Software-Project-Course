import symnmfModule
import kmeans
import sklearn.metrics as sk
import sys
import pandas as pd
import numpy as np
import math



def main():
    np.random.seed(0)
    if (len(sys.argv) < 2 or len(sys.argv) > 3):
        print("An Error Has Occurred")
        return

    k = sys.argv[1]
    file_name = sys.argv[2]

    try: 
        data = pd.read_csv(file_name, header=None)
    except Exception as e:
        print("An Error Has Occurred")
        sys.exit(1)

    vectorX = [x.tolist() for index, x in data.iterrows()]
    if int(k) >= len(vectorX) or len(vectorX) == 0:
        print("An Error Has Occurred")
        sys.exit(1)


    k=int(k)
    n=len(vectorX)
    d=len(vectorX[0])

    
    kmeans_clusters =  kmeans.InitalizeKMeans(k,300,file_name)

   
    W = symnmfModule.norm(vectorX)
    mean=np.mean(W)
    inital_H=np.random.uniform(0,2*math.sqrt(mean/k),size=(n,k))
    inital_H_list = inital_H.tolist()
    H =symnmfModule.symnmf(vectorX,inital_H_list,k, n,d)
    symnmflist = symnmfModule.analysis(H, len(vectorX),k)
        
        
        
    
    k_points = np.array(kmeans_clusters)
    
    nearest_clusters = []

    for index, row in data.iterrows():

        distances = np.linalg.norm(k_points - row.values, axis=1)

        nearest = np.argmin(distances)
        nearest_clusters.append(nearest)

    kmeans_labels_flat =  nearest_clusters


    nmf = sk.silhouette_score(vectorX, symnmflist)
    kmean = sk.silhouette_score(vectorX, kmeans_labels_flat)


    print("nmf: {:.4f}".format(nmf))
    print("kmeans: {:.4f}".format(kmean))


if __name__ == "__main__":
    main()


