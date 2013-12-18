/*
 * Main developer: Nicolas Van Cleemput
 * 
 * Copyright (C) 2013 Nicolas Van Cleemput.
 * Licensed under the GNU GPL, read the file LICENSE.txt for details.
 */

/* 
 * The vertices on the lefthand side belong to the 4-partition and will all have
 * degree 4 for each finished graph.
 */

#define WITH_PROFILING

#include "generate34Bipartite.h"

// debugging methods

void depthPrint(){
    int i;
    for(i = 0; i < rightVertexCount; i++){
        fprintf(stderr, " ");
    }
}

void setPrint(bitv s, int maxSize){
    int i;
    depthPrint();
    for(i = 0; i < maxSize; i++){
        if(CONTAINS(s, i)){
            fprintf(stderr, "%d ", i);
        }
    }
    fprintf(stderr, "\n");
}

void printCurrentGraph(){
    int i, j;
    for(i = 0; i < leftVertexCount; i++){
        depthPrint();
        fprintf(stderr, "%d: ", i);
        for(j = 0; j < rightVertexCount; j++){
            if(CONTAINS(leftNeighbourhood[i], j)){
                fprintf(stderr, "%d ", leftVertexCount + j);
            }
        }
        fprintf(stderr, "\n");
    }
    for(i = 0; i < rightVertexCount; i++){
        depthPrint();
        fprintf(stderr, "%d: ", leftVertexCount + i);
        for(j = 0; j < leftVertexCount; j++){
            if(CONTAINS(rightNeighbourhood[i], j)){
                fprintf(stderr, "%d ", j);
            }
        }
        fprintf(stderr, "\n");
    }
    fprintf(stderr, "\n");
}

void printNautySparseGraph(sparsegraph *sgraph){
    int i, j;
    depthPrint();
    fprintf(stderr, "nauty sparsegraph:\n");
    for(i = 0; i < sgraph->nv; i++){
        depthPrint();
        fprintf(stderr, "%d: ", i);
        for(j = 0; j < sgraph->d[i]; j++){
            fprintf(stderr, "%d ", sgraph->e[sgraph->v[i] + j]);
        }
        fprintf(stderr, "\n");
    }
}

void printGenerators(int depth){
    int i, j;
    depthPrint();
    fprintf(stderr, "Generators:\n");
    for(i = 0; i < generatorCount[depth]; i++){
        depthPrint();
        for(j = 0; j < leftVertexCount + depth; j++){
            fprintf(stderr, "%d ", generators[depth][i][j]);
        }
        fprintf(stderr, "\n");
    }
}

void debugAddRightVertex(int n1, int n2, int n3){
    int j;
    bitv s = EMPTY_SET;
    ADD(s, n1);
    ADD(s, n2);
    ADD(s, n3);
    rightNeighbourhood[rightVertexCount] = s;
    for(j = 0; j < leftVertexCount; j++){
        if(CONTAINS(s, j)){
            ADD(leftNeighbourhood[j], rightVertexCount);
            leftVertexDegree[j]++;
        }
    }
    rightVertexCount++; 
    generatorsDetermined[rightVertexCount] = FALSE;
}

void debugInit(int degree4Count, int degree3Count){
    leftVertexCount = degree4Count;
    maximumRightVertexCount = degree3Count;
    
    initNautyRelatedVariables();
    generatorsDetermined[0] = FALSE;
}

// end debugging methods

void outputMulticode(){
    int i, j;
    
    static boolean first = TRUE;
    
    if(first){
        fprintf(outputFile, ">>multi_code<<");
        first = FALSE;
    }
    
    putc(leftVertexCount + rightVertexCount, outputFile);
    
    for(i = 0; i < leftVertexCount; i++){
        for(j = 0; j < MAXNRIGHT; j++){
            if (CONTAINS(leftNeighbourhood[i], j)){
                putc(j + leftVertexCount + 1, outputFile); 
            }
        }
        putc(0, outputFile);
    }
    for(i = 0; i < rightVertexCount - 1; i++){
        putc(0, outputFile);
    }

}

/**
 * Method which is called each time nauty finds a generator.
 */
void storeGenerators(int count, permutation perm[], nvector orbits[], int numorbits, int stabvertex, int n) {
    memcpy(generators[rightVertexCount] + generatorCount[rightVertexCount], perm, sizeof (permutation) * n);

    generatorCount[rightVertexCount]++;
}

void initNautyRelatedVariables(){
    
    options.getcanon = TRUE;
    options.defaultptn = FALSE;
    options.userautomproc = storeGenerators;

    /* Init the nauty datastructures */
    SG_INIT(sg);
    SG_ALLOC(sg, leftVertexCount + maximumRightVertexCount, 6 * maximumRightVertexCount, "Failed to allocate memory to store graph");

    //sg.v only has to be set once for the left vertices
    int i;
    for(i = 0; i < leftVertexCount; i++) {
        sg.v[i] = i * 4;
    }
    //sg.v and sg.d only have to be set once for the right vertices
    for(i = 0; i < maximumRightVertexCount; i++) {
        sg.v[leftVertexCount + i] = leftVertexCount * 4 + i * 3;
        sg.d[leftVertexCount + i] = 3;
    }

    SG_INIT(sg_canon);
    SG_ALLOC(sg_canon, leftVertexCount + maximumRightVertexCount, 6 * maximumRightVertexCount, "Failed to allocate memory to store graph");
}

inline void prepareNautyCall(){
    generatorCount[rightVertexCount] = 0;
}

/* This method translates the internal data structure to nauty's sparse graph
 * data structure, so the graph can be passed to nauty.
 */
inline void translateCurrentGraphToNautySparseGraph(){
    sg.nv = leftVertexCount + rightVertexCount;
    sg.nde = 6 * rightVertexCount; //right vertices have degree 3, and these are all edges

    int i, j;
    for(i = 0; i < leftVertexCount;i++) {
        sg.d[i] = leftVertexDegree[i];
        int k = 0;
        for(j = 0; j < leftVertexDegree[i]; j++) {
            while(!CONTAINS(leftNeighbourhood[i], k)){
                k++;
            }
            sg.e[i * 4 + j] = leftVertexCount + k;
            k++;
        }
    }
    for(i = 0; i < rightVertexCount; i++) {
        int k = 0;
        for(j = 0; j < 3; j++) {
            while(!CONTAINS(rightNeighbourhood[i], k)){
                k++;
            }
            sg.e[leftVertexCount * 4 + i*3 + j] = k;
            k++;
        }
    }
}

void callNauty(){
    int i;
    //translate graph to nauty graph
    translateCurrentGraphToNautySparseGraph();
    
    //partition the vertices
    int labellingIndex = 0;
    for(i = 0; i < rightVertexCount; i++){
        lab[labellingIndex] = leftVertexCount + i;
        ptn[labellingIndex] = 1;
        labellingIndex++;
    }
    ptn[labellingIndex - 1] = 0;
    
    for(i = 0; i < leftVertexCount; i++){
        lab[labellingIndex] = i;
        ptn[labellingIndex] = 1;
        labellingIndex++;
    }
    ptn[labellingIndex - 1] = 0;
    
    //call nauty
    prepareNautyCall();
    nauty((graph*) &sg, lab, ptn, NULL, orbits, &options, &stats, workspace, WORKSIZE, MAXM, leftVertexCount + rightVertexCount, (graph*) &sg_canon);
    
    generatorsDetermined[rightVertexCount] = TRUE;
}

boolean doesLastVertexLieInOrbitOfSmallestLabel(bitv maxColouredVertices){
    int i;
    //translate graph to nauty graph
    translateCurrentGraphToNautySparseGraph();
    
    //partition the vertices
    int labellingIndex = 0;
    for(i = 0; i < rightVertexCount; i++){
        if(CONTAINS(maxColouredVertices, i)){
            lab[labellingIndex] = leftVertexCount + i;
            ptn[labellingIndex] = 1;
            labellingIndex++;
        }
    }
    ptn[labellingIndex - 1] = 0;
    
    if(labellingIndex!=rightVertexCount){
        for(i = 0; i < rightVertexCount; i++){
            if(!CONTAINS(maxColouredVertices, i)){
                lab[labellingIndex] = leftVertexCount + i;
                ptn[labellingIndex] = 1;
                labellingIndex++;
            }
        }
        ptn[labellingIndex - 1] = 0;
    }
    
    for(i = 0; i < leftVertexCount; i++){
        lab[labellingIndex] = i;
        ptn[labellingIndex] = 1;
        labellingIndex++;
    }
    ptn[labellingIndex - 1] = 0;
    
    //call nauty
    prepareNautyCall();
    nauty((graph*) &sg, lab, ptn, NULL, orbits, &options, &stats, workspace, WORKSIZE, MAXM, leftVertexCount + rightVertexCount, (graph*) &sg_canon);

    //check whether last vertex lies in the orbit of the canonical vertex
    int reverseLabelling[leftVertexCount + rightVertexCount];
    for (i = 0; i < leftVertexCount + rightVertexCount; i++) {
        reverseLabelling[lab[i]]=i;
    }
    int smallestLabelOrbitLastVertex = reverseLabelling[leftVertexCount + rightVertexCount - 1];
    //i.e. the smallest label of a vertex in the orbit of the last vertex
    int smallestOtherRightVertexLabel = INFINITY;
    for (i = leftVertexCount; i < leftVertexCount + rightVertexCount; i++) {
        if(orbits[i]==orbits[leftVertexCount + rightVertexCount - 1]){
            if(reverseLabelling[i]<smallestLabelOrbitLastVertex) {
                smallestLabelOrbitLastVertex = reverseLabelling[i];
            }
        } else if(CONTAINS(maxColouredVertices, i - leftVertexCount)){
            if(reverseLabelling[i]<smallestOtherRightVertexLabel) {
                smallestOtherRightVertexLabel = reverseLabelling[i];
            }
        }
    }

    
    generatorsDetermined[rightVertexCount] = TRUE;
#ifdef WITH_PROFILING
    if(smallestLabelOrbitLastVertex <= smallestOtherRightVertexLabel){
        nautyAccepted++;
    } else {
        nautyRejected++;
    }
#endif
    return smallestLabelOrbitLastVertex <= smallestOtherRightVertexLabel;
}

/* v is the index of a vertex in the righthand side partition
 */
int getNumberOfVerticesAtDistance2(int v){
    int i;
    bitv set = EMPTY_SET;
    bitv_size setSize;
    
    for(i = 0; i < MAXNLEFT; i++){
        if(CONTAINS(rightNeighbourhood[v], i)){
            ADD_ALL(set, leftNeighbourhood[i]);
        } 
    }
    
    SETSIZE(set, setSize)
    return setSize;
}

int getNumberOfNeighboursWithDegree(int v, int deg){
    //returns the number of neighbours with a given degree for a right vertex
    int i, count = 0;
    
    for(i = 0; i < MAXNLEFT; i++){
        if(CONTAINS(rightNeighbourhood[v], i)){
            if(leftVertexDegree[i]==deg) count++;
        } 
    }
    return count;
}

/* Check whether vertex at position maxReducibleVertex is canonical
 */
boolean isLastVertexCanonical(){
    int i, maxCount; 
    
#ifdef WITH_PROFILING
    canonicityCalls++;
#endif
    
    //this set holds the vertices that can be the canonical vertex
    bitv maxVertices = EMPTY_SET;
    
    //////////////////////////////////////
#ifdef NO_COLOURS
    for(i = 0; i < rightVertexCount; i++){
        ADD(maxVertices, i);
    }
    
    return doesLastVertexLieInOrbitOfSmallestLabel(maxVertices);
#endif
    ///////////////////////////////////////
    
    int colour = getNumberOfNeighboursWithDegree(rightVertexCount-1, 4);
    
    maxCount = 0;
    for(i = 0; i < rightVertexCount-1; i++){
        int iColour = getNumberOfNeighboursWithDegree(i, 4);
        if(iColour < colour){
#ifdef WITH_PROFILING
            colour1Rejected++;
#endif
            return FALSE;
        } else if(iColour == colour){
            maxCount++;
            ADD(maxVertices, i);
        }
    }
    if(maxCount==0){
#ifdef WITH_PROFILING
        colour1Accepted++;
#endif
        return TRUE;
    }
    
    colour = getNumberOfVerticesAtDistance2(rightVertexCount-1);
    
    maxCount = 0;
    for(i = 0; i < rightVertexCount-1; i++){
        if(CONTAINS(maxVertices, i)){
            int iColour = getNumberOfVerticesAtDistance2(i);
            if(iColour > colour){
#ifdef WITH_PROFILING
                colour2Rejected++;
#endif
                return FALSE;
            } else if(iColour == colour){
                maxCount++;
            } else {
                REMOVE(maxVertices, i);
            }
        }
    }
    if(maxCount==0) {
#ifdef WITH_PROFILING
        colour2Accepted++;
#endif
        return TRUE;
    }
    
    colour = getNumberOfNeighboursWithDegree(rightVertexCount-1, 3);
    
    maxCount = 0;
    for(i = 0; i < rightVertexCount-1; i++){
        if(CONTAINS(maxVertices, i)){
            int iColour = getNumberOfNeighboursWithDegree(i, 3);
            if(iColour < colour){
#ifdef WITH_PROFILING
                colour3Rejected++;
#endif
                return FALSE;
            } else if(iColour == colour){
                maxCount++;
            } else {
                REMOVE(maxVertices, i);
            }
        }
    }
    if(maxCount==0) {
#ifdef WITH_PROFILING
        colour3Accepted++;
#endif
        return TRUE;
    }
    
    colour = getNumberOfNeighboursWithDegree(rightVertexCount-1, 2);
    
    maxCount = 0;
    for(i = 0; i < rightVertexCount-1; i++){
        if(CONTAINS(maxVertices, i)){
            int iColour = getNumberOfNeighboursWithDegree(i, 2);
            if(iColour < colour){
#ifdef WITH_PROFILING
                colour4Rejected++;
#endif
                return FALSE;
            } else if(iColour == colour){
                maxCount++;
            } else {
                REMOVE(maxVertices, i);
            }
        }
    }
    if(maxCount==0) {
#ifdef WITH_PROFILING
        colour4Accepted++;
#endif
        return TRUE;
    }
    
    //if we got here, then the last vertex also has the maximum colour
    ADD(maxVertices, rightVertexCount-1);
    
    return doesLastVertexLieInOrbitOfSmallestLabel(maxVertices);
    
}

void constructAllExtendible3Sets(){
    int i, j, k;
    vertex3SetCount[rightVertexCount] = 0;
    
    for(i = 2; i < leftVertexCount; i++){
        if(leftVertexDegree[i]!=4){
            for(j = 1; j < i; j++){
                if(leftVertexDegree[j]!=4){
                    for(k = 0; k < j; k++){
                        if(leftVertexDegree[k]!=4){
                            //set i,j,k is extendible
                            bitv newSet = EMPTY_SET;
                            ADD(newSet, i);
                            ADD(newSet, j);
                            ADD(newSet, k);
                            vertex3Sets[rightVertexCount][vertex3SetCount[rightVertexCount]] = newSet;
                            vertex3SetCount[rightVertexCount]++;
                        }
                    }
                }
            }
        }
    }
}

//implementation of union-find algorithm for sets of vertices
int findRootOfElement(int *forest, int element) {
    //find with path-compression
    if(element!=forest[element]){
        forest[element]=findRootOfElement(forest, forest[element]);
    }
    return forest[element];
}

//implementation of union-find algorithm for sets of vertices
void unionElements(int *forest, int *treeSizes, int *numberOfComponents, int element1, int element2){
    int root1 = findRootOfElement(forest, element1);
    int root2 = findRootOfElement(forest, element2);

    if(root1==root2) return;

    if(treeSizes[root1]<treeSizes[root2]){
        forest[root1]=root2;
        treeSizes[root2]+=treeSizes[root1];
    } else {
        forest[root2]=root1;
        treeSizes[root1]+=treeSizes[root2];
    }
    (*numberOfComponents)--;
}

/*
 * determines the orbits of a given list of vertex sets of size 3. The sets must
 * be so that the image of a set in the list under the automorphism group of the
 * graph is also in the list.
 *
 * In the end vertex3SetOrbitRepresentative[i] will contain the index of the 
 * canonical representative of the orbit to which vertex3Sets[i] belongs and 
 * vertex3SetOrbitCount will contain the number of orbits.
 * 
 * To get all orbit representatives, just look for the sets which have
 * vertex3SetOrbitRepresentative[i]==i.
 */
void determineVertex3SetsOrbits(){
    //it is assumed that all sets have 3 elements
    int i, j, k, l;
    int orbitSize[MAXNLEFT * (MAXNLEFT - 1) * (MAXNLEFT - 2)/6];

    //initialization of the variables
    for(i = 0; i < vertex3SetCount[rightVertexCount] ; i++){
        vertex3SetOrbitRepresentative[rightVertexCount][i]=i;
        orbitSize[i]=1;
    }
    vertex3SetOrbitCount[rightVertexCount] = vertex3SetCount[rightVertexCount];

    permutation *permutation;
    bitv setImage;
    
    for(i = 0; i < generatorCount[rightVertexCount]; i++) {
        //the generators were stored in the global variable generators by the method storeGenerators
        permutation = generators[rightVertexCount][i];

        for(j = 0; j < vertex3SetCount[rightVertexCount]; j++){
            //apply permutation to current vertex pair
            setImage = EMPTY_SET;
            for(l = 0; l < leftVertexCount; l++){
                if(CONTAINS(vertex3Sets[rightVertexCount][j],l)){
                    ADD(setImage, permutation[l]);
                }
            }

            //search the pair in the list
            for(k = 0; k<vertex3SetCount[rightVertexCount]; k++){
                if(setImage == vertex3Sets[rightVertexCount][k]){
                    //union j and k
                    unionElements(vertex3SetOrbitRepresentative[rightVertexCount], orbitSize,
                            vertex3SetOrbitCount +rightVertexCount, j, k);
                    //the list of sets doesn't contain any duplicates so we can stop
                    break;
                }
            }
        }
    }
}

void addNextVertex(){
    int i, j;
    
    //output graph if we reached the end
    if(rightVertexCount == maximumRightVertexCount){
        outputMulticode();
        graphCount++;
        return;
    }
    
    //check to see if automorphism group still needs to be calculated
    if(!generatorsDetermined[rightVertexCount]){
        callNauty();
    }
    
    //determine orbits of sets of 3 vertices
    constructAllExtendible3Sets();
    determineVertex3SetsOrbits();
    
    //for each orbit of 3-sets: pick a representative and add vertex adjacent to that set
    for(i = 0; i < vertex3SetCount[rightVertexCount] ; i++){
        if(vertex3SetOrbitRepresentative[rightVertexCount][i]==i){
            //add vertex adjacent to all vertices in vertex3Sets[i]
            rightNeighbourhood[rightVertexCount] = vertex3Sets[rightVertexCount][i];
            for(j = 0; j < leftVertexCount; j++){
                if(CONTAINS(vertex3Sets[rightVertexCount][i], j)){
                    ADD(leftNeighbourhood[j], rightVertexCount);
                    leftVertexDegree[j]++;
                }
            }
            rightVertexCount++;
            generatorsDetermined[rightVertexCount] = FALSE;
    
            //check canonicity: add next vertex if canonical
            if(isLastVertexCanonical()){
                addNextVertex();
            }
            
            //undo adding of vertex
            rightVertexCount--;
            for(j = 0; j < leftVertexCount; j++){
                if(CONTAINS(vertex3Sets[rightVertexCount][i], j)){
                    REMOVE(leftNeighbourhood[j], rightVertexCount);
                    leftVertexDegree[j]--;
                }
            }
        }
    }
}

void startGeneration(int degree4Count, int degree3Count){
    leftVertexCount = degree4Count;
    maximumRightVertexCount = degree3Count;
    
    initNautyRelatedVariables();
    generatorsDetermined[0] = FALSE;
    
    //start the generation
    addNextVertex();
}

int main(int argc, char** argv) {
    outputFile = stdout;

    startGeneration(9, 12);
    fprintf(stderr, "Found %lu graphs.\n", graphCount);
    
#ifdef WITH_PROFILING
    
    fprintf(stderr, "Canonicity calls: %10llu\n\n", canonicityCalls);
    fprintf(stderr, "C1      rejected: %10llu\n", colour1Rejected);
    fprintf(stderr, "        accepted: %10llu\n", colour1Accepted);
    fprintf(stderr, "C2      rejected: %10llu\n", colour2Rejected);
    fprintf(stderr, "        accepted: %10llu\n", colour2Accepted);
    fprintf(stderr, "C3      rejected: %10llu\n", colour3Rejected);
    fprintf(stderr, "        accepted: %10llu\n", colour3Accepted);
    fprintf(stderr, "C4      rejected: %10llu\n", colour4Rejected);
    fprintf(stderr, "        accepted: %10llu\n", colour4Accepted);
    fprintf(stderr, "nauty   rejected: %10llu\n", nautyRejected);
    fprintf(stderr, "        accepted: %10llu\n", nautyAccepted);

#endif
    
    return (EXIT_SUCCESS);
}

