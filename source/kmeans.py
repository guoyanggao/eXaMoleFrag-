import math
import random
import time

######################################################################
# This section contains functions for clustering a dataset
# using the k-means algorithm.
######################################################################
# input vectors are booked in dict.

#计算Euclidean
def distance(instance1, instance2):
    if instance1 == None or instance2 == None:
        return float("inf")
    sumOfSquares = 0
    for i in range(1, len(instance1)): # starts with 1 because a label in instances1[0]
        sumOfSquares += (instance1[i] - instance2[i])**2
    return sumOfSquares

def createEmptyListOfLists(numSubLists):
    myList = []
    for i in range(numSubLists):
        myList.append([])
    return myList

#找离instance最近的“质心”    遍历所有质心
def assign(instances, instance, centroids):
    minDistance = distance(instances[instance], centroids[0])
    minDistanceIndex = 0
    for i in range(1, len(centroids)):
        d = distance(instances[instance], centroids[i])
        if (d < minDistance):
            minDistance = d
            minDistanceIndex = i
    return minDistanceIndex

def assignAll(instances, centroids):
    clusters = createEmptyListOfLists(len(centroids))    #是一个空列表
    # cycle each instance and assign to a cluster
    for instance in range(1,len(instances)+1): 
        clusterIndex = assign(instances, instance, centroids)
        clusters[clusterIndex].append(instance)    #给第i个质心，append上这个节点
    return clusters




def meanInstance(name, instances, instanceList):
    numInstances = len(instanceList)
    if (numInstances == 0):
        return

    # dimmension of array
    numAttributes = len(instances[instanceList[0]])
    # init means array
    means = [name] + [0] * (numAttributes - 1)
    for instance in instanceList:
        for i in range(1, numAttributes):
            means[i] += instances[instance][i]
    for i in range(1, numAttributes):
        means[i] /= float(numInstances)
    return tuple(means)

def computeCentroids(instances, clusters):
    centroids = []
    for i in range(0,len(clusters)):
        name = "centroid" + str(i+1)
        centroid = meanInstance(name, instances, clusters[i])
        centroids.append(centroid)
    return centroids

def kmeans(instances, k, animation=False, initClusters=None):
    # input instances is booked in a dict with N verctors
    # [1: [x1, y1, z1, ...], .... ]
    # initclusters:
    # [[1,2],[3,4..], ...]

    result = {}
    # init clusters if no input
    if (initClusters == None or len(initClusters) < k):
        # randomly select k initial centroids
        random.seed(time.time())
        initClusters = random.sample(instances, k)
        clusters = [initClusters[i:i+1] for i in range(0,len(initClusters))]
    else:
        clusters = initClusters

    prevCentroids = []
    centroids = computeCentroids(instances, clusters)
             
    #if animation:
    #    delay = 1.0 # seconds
    #    canvas = prepareWindow(instances)
    #    clusters = createEmptyListOfLists(k)
    #    clusters[0] = instances
    #    paintClusters2D(canvas, clusters, centroids, "Initial centroids")
    #    time.sleep(delay)

    #print " Iter        Withiss"
    iteration = 0
    while (centroids != prevCentroids):
        iteration += 1
        clusters = assignAll(instances, centroids)
        #print "iter", iteration
        #print clusters

        #if animation:
        #    paintClusters2D(canvas, clusters, centroids, "Assign %d" % iteration)
        #    time.sleep(delay)
        prevCentroids = centroids
        centroids = computeCentroids(instances, clusters)
        withinss  = computeWithinss(instances, clusters, centroids)
        #if animation:
        #    paintClusters2D(canvas, clusters, centroids,
        #                    "Update %d, withinss %.1f" % (iteration, withinss))
        #    time.sleep(delay)
        #print " %4d    %12.4f" % (iteration,withinss)


    result["clusters"] = clusters
    result["centroids"] = centroids
    result["withinss"] = withinss
    return result

def computeWithinss(instances, clusters, centroids):
    result = 0
    for i in range(len(centroids)):
        centroid = centroids[i]
        cluster = clusters[i]
        for instance in cluster:
            result += distance(centroid, instances[instance])
    return result

# Repeats k-means clustering n times, and returns the clustering
# with the smallest withinss
def repeatedKMeans(instances, k, n):
    bestClustering = {}
    bestClustering["withinss"] = float("inf")
    for i in range(1, n+1):
        print("k-means trial %d," % i ,)
        trialClustering = kmeans(instances, k)
        print("withinss: %.1f" % trialClustering["withinss"])
        if trialClustering["withinss"] < bestClustering["withinss"]:
            bestClustering = trialClustering
            minWithinssTrial = i
    print("Trial with minimum withinss:", minWithinssTrial)
    return bestClustering


