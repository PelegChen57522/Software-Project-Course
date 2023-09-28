import sys
import math
Epsilon=0.001
N=-1 #number of data points

def InitalizeKMeans(K,numOfIter,inputFileData):
    global N # number of datapoints that were given
    with open(inputFileData, "r") as inputFile:
        # take the date points into the lists
        dataPoints = [[float(x) for x in line.split(",")] for line in inputFile.readlines()]
        # Initialize centroids as first k datapoints: µk = xk, ∀k ∈ K
        centroids = dataPoints[:K]

    # number of datapoints that were given
    N = len(dataPoints)

    # check if 1 < K < N

    if (K > 1 and K < N) == False:
        print("Invalid number of clusters!")
        sys.exit()

    return mainKMeans(K,numOfIter,dataPoints,centroids)

def mainKMeans(K, numOfIter,dataPoints,centroids):
    for iter in range(numOfIter):
        newCentroids=computeCentroids(dataPoints,centroids)
        convergence= False
        for j in range(K):
            dist = EuclideanDistance(newCentroids[j],centroids[j])
            if (dist >= Epsilon): #Check if not convergence
                centroids = newCentroids #update centroids
                convergence = True
                break
        if (convergence == False):
            break
    #print
    for centroid in newCentroids:
        print(','.join(['%.4f' % c for c in centroid]))


def computeCentroids(dataPoints, centroids):
    dimension=len(centroids[0]) #dimension is initialized to the length of the first centroid, which is the same as the number of features in each data point.
    sums = [[0 for j in range(dimension)] for i in range(len(centroids))] #sums is a 2D list initialized with all zeros. It will be used to store the sum of all data points assigned to each centroid.
    counts = [0 for j in range(len(centroids))] #counts is a 1D list initialized with all zeros. It will be used to store the number of data points assigned to each centroid.

    #Assign every xi to the closest cluster k: argminkd(xi, µk), ∀k 1 ≤ k ≤ K
    for dataPoint in dataPoints:
        min_index = 0
        min_dist = float('inf')
        for j in range(len(centroids)): #This is a loop that iterates over all centroids.
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

    #Update the centroids: µk =1/|k|sums of xi∈k(xi)
    for i in range(len(centroids)):
        if counts[i] != 0:
           for l in range(dimension):
                sums[i][l] /= counts[i]

    return sums


# compute the Euclidean Distance between 2 vectors
def EuclideanDistance(v1,v2):
    assert len(v1)==len(v2)
    sum=0
    for i in range(len(v1)):
        sum+=(math.pow((v1[i]-v2[i]),2))
    res=math.sqrt(sum)
    return res



def main(argv):
    if (len(argv) < 3 or len(argv) > 4): #the user gives less/more inputs from what we need
        print("Inavlid Input!")
        sys.exit()

    if (len(argv) == 4):  # argument "iter" has been given by user
        # checking if the first 2 arguments (K, iter) are integers
        if (argv[1].isnumeric() == False):
            print("Invalid number of clusters!")
            sys.exit()
        if (argv[2].isnumeric() == False or ((int)(argv[2]) <= 1) or  ((int)(argv[2]) >= 1000)):#check the number of iterations were given
            print("Invalid maximum iteration!")
            sys.exit()
        else:
            K = (int)(argv[1])
            numOfIter = (int)(argv[2])
            return InitalizeKMeans(K,  numOfIter,argv[3])

    if (len(argv) == 3):  # argument "iter" has not been given by user
        # checking if  K is integer
        if (argv[1].isnumeric() == False):
            print("Invalid number of clusters!")
            sys.exit()
        else:
            K = (int)(argv[1])
            numOfIter = 200 #default
            return InitalizeKMeans(K,numOfIter, argv[2])

if __name__ == '__main__':
    main(sys.argv)
