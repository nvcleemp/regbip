/*
 * Main developer: Nicolas Van Cleemput
 * 
 * Copyright (C) 2013 Nicolas Van Cleemput.
 * Licensed under the GNU GPL, read the file LICENSE.txt for details.
 * 
 * This program requires the multicode library files from the graphtools project.
 */

#include <stdio.h>
#include <stdlib.h>
#include"shared/multicode_base.h"
#include"shared/multicode_output.h"
#include"shared/multicode_input.h"

boolean checkGraph(GRAPH graph, ADJACENCY adj){
    int i, j;
    
    if(graph[0][0]%7){
        return FALSE;
    }
    
    for(i = 1; i <= graph[0][0]; i++){
        if(adj[i]==3){
            for(j = 0; j < 3; j++){
                if(adj[graph[i][j]]!=4){
                    return FALSE;
                }
            }
        } else if(adj[i]==4){
            for(j = 0; j < 4; j++){
                if(adj[graph[i][j]]!=3){
                    return FALSE;
                }
            }
        } else {
            return FALSE;
        }
    }
    
    return TRUE;
}

int main(int argc, char** argv) {
    GRAPH graph;
    ADJACENCY adj;
    
    int readGraphs = 0;
    int filteredGraphs = 0;
    
    unsigned short code[MAXCODELENGTH];
    int length;
    while (readMultiCode(code, &length, stdin)) {
        decodeMultiCode(code, length, graph, adj);
        
        readGraphs++;
        
        if(checkGraph(graph, adj)){
            filteredGraphs++;
            writeMultiCode(graph, adj, stdout);
        }
    }
    
    fprintf(stderr, "Read %d graph%s, filtered %d graph%s.\n",
            readGraphs, readGraphs==1 ? "" : "s",
            filteredGraphs, filteredGraphs==1 ? "" : "s");

    return (EXIT_SUCCESS);
}

