/*
 * Main developer: Nicolas Van Cleemput
 * 
 * Copyright (C) 2013 Nicolas Van Cleemput.
 * Licensed under the GNU GPL, read the file LICENSE.txt for details.
 */

#ifndef GENERATE34BIPARTITE_H
#define	GENERATE34BIPARTITE_H

#include <stdio.h>
#include <stdlib.h>

//can handle up to 32 vertices of degree 3
#define MAXNLEFT 24
#define MAXNRIGHT 32

#define MAXN MAXNLEFT + MAXNRIGHT

#define INFINITY 2*(MAXNLEFT + MAXNRIGHT)

#include "nauty/nausparse.h"

#include "bitvectors.h"

/* Nauty worksize */
#define WORKSIZE 50 * MAXM

/** Nauty variables */
int lab[MAXN], ptn[MAXN], orbits[MAXN];
static DEFAULTOPTIONS_SPARSEGRAPH(options);
statsblk stats;
setword workspace[WORKSIZE];

sparsegraph sg; /* Sparse graph datastructure for nauty */
sparsegraph sg_canon; /* Sparse graph datastructure for nauty */

permutation generators[MAXNRIGHT][MAXN][MAXN];
int generatorCount[MAXNRIGHT];
boolean generatorsDetermined[MAXNRIGHT];

// end of nauty variables

FILE *outputFile;

int leftVertexCount; //this number stays fixed during the run of the program
int rightVertexCount; //this number increases with the addition of each vertex
int maximumRightVertexCount; //this number increases with the addition of each vertex

bitv leftNeighbourhood[MAXNLEFT];
bitv rightNeighbourhood[MAXNRIGHT];
int leftVertexDegree[MAXNLEFT];

//3-sets (used for determining possible extensions)
bitv vertex3Sets[MAXNRIGHT][MAXNLEFT * (MAXNLEFT - 1) * (MAXNLEFT - 2)/6];
int vertex3SetCount[MAXNRIGHT];
int vertex3SetOrbitRepresentative[MAXNRIGHT][MAXNLEFT * (MAXNLEFT - 1) * (MAXNLEFT - 2)/6];
int vertex3SetOrbitCount[MAXNRIGHT];

unsigned long int graphCount = 0;


void outputMulticode();

void storeGenerators(int count, permutation perm[], nvector orbits[], int numorbits, int stabvertex, int n);

void initNautyRelatedVariables();

inline void prepareNautyCall();

inline void translateCurrentGraphToNautySparseGraph();

void callNauty();

boolean doesLastVertexLieInOrbitOfSmallestLabel(bitv maxColouredVertices);

int getNumberOfVerticesAtDistance2(int v);

int getNumberOfNeighboursWithDegree(int v, int deg);

boolean isLastVertexCanonical();

void constructAllExtendible3Sets();

void determineVertex3SetsOrbits();

void addNextVertex();

void startGeneration(int degree4Count, int degree3Count);


#endif	/* GENERATE34BIPARTITE_H */

