import sys
import math
Epsilon=0.0001
N=-1 

def InitalizeKMeans(K,numOfIter,inputFileData):
    global N 
    with open(inputFileData, "r") as inputFile:
        
        dataPoints = [[float(x) for x in line.split(",")] for line in inputFile.readlines()]
        
    
    N = len(dataPoints)
    

   

    if (K > 1 and K < N) == False:
        print("Invalid number of clusters!")
        sys.exit()

    return mainKMeans(K,numOfIter,dataPoints)

def mainKMeans(K, numOfIter, dataPoints):
    centroids = dataPoints[:K]
    for iter in range(numOfIter):
        newCentroids = computeCentroids(dataPoints, centroids)
        convergence = True
        for j in range(K):
            dist = EuclideanDistance(newCentroids[j], centroids[j])
            if dist > Epsilon:  
                centroids = newCentroids
                convergence = False
                break
        if convergence:
            break

    return newCentroids


def computeCentroids(dataPoints, centroids):
    dimension=len(centroids[0]) 
    sums = [[0 for j in range(dimension)] for i in range(len(centroids))] 
    counts = [0 for j in range(len(centroids))] 

   
    for dataPoint in dataPoints:
        min_index = 0
        min_dist = float('inf')
        for j in range(len(centroids)): 
            dist = 0
            for l in range(dimension):
                dist += (math.pow(dataPoint[l] - centroids[j][l], 2))

            dist = math.sqrt(dist)
            if (dist < min_dist):
                min_dist = dist
                min_index = j

        counts[min_index] += 1
        for m in range(len(dataPoint)):
            sums[min_index][m] += dataPoint[m]

    
    for i in range(len(centroids)):
        if counts[i] != 0:
           for l in range(dimension):
                sums[i][l] /= counts[i]

    return sums



def EuclideanDistance(v1,v2):
    assert len(v1)==len(v2)
    sum=0
    for i in range(len(v1)):
        sum+=(math.pow((v1[i]-v2[i]),2))
    res=math.sqrt(sum)
    return res



def main(argv):
    if (len(argv) < 3 or len(argv) > 4): 
        print("Inavlid Input!")
        sys.exit()

    if (len(argv) == 4): 
       
        if (argv[1].isnumeric() == False):
            print("Invalid number of clusters!")
            sys.exit()
        if (argv[2].isnumeric() == False or ((int)(argv[2]) <= 1) or  ((int)(argv[2]) >= 1000)):
            print("Invalid maximum iteration!")
            sys.exit()
        else:
            K = (int)(argv[1])
            numOfIter = (int)(argv[2])
            return InitalizeKMeans(K,  numOfIter,argv[3])

    if (len(argv) == 3):  
        
        if (argv[1].isnumeric() == False):
            print("Invalid number of clusters!")
            sys.exit()
        else:
            K = (int)(argv[1])
            numOfIter = 200 
            return InitalizeKMeans(K,numOfIter, argv[2])

if __name__ == '__main__':
    main(sys.argv)
