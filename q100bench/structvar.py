import sys

def write_structural_errors(refobj, queryobj, outputdict, bmstats, args)->str:

    structstatsfile = outputdict["structdetailsfile"]
    structerrorvcf = outputdict["structvariantsvcf"]

    alignclusters = bmstats["alignclusters"]

    for refentry in alignclusters.keys():
        numclusters = len(alignclusters[refentry])
        print("Ref " + refentry + " has " + str(numclusters) + " clusters")
        clusternum = 1
        for cluster in alignclusters[refentry]:
            clusterlength = len(cluster["aligns"])
            print("\tRef " + refentry + " cluster " + str(clusternum) + " has " + str(clusterlength) + " aligns")
            clusternum = clusternum + 1

    return 0

