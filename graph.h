
#ifndef DMCE3_GRAPH_H
#define DMCE3_GRAPH_H
#include "gadget.h"
#include <vector>
#include <string>

class Graph
{
public:
    Graph();
    Graph(int n);
    Graph(char* fname);
    ~Graph();

    void build(char *fname);
    void addEdge(int,int);     void setDeg(int,int);
    void print();
    int dest(int v, int k);
    int deg(int v);
    int V();
    int E();
    int hindex();
    VI& adjvec(int i) {return gvec[i];}
    bool hasCore(int k);

private:
    VVI gvec;
    static const int maxe = 200000000;
    int n, m;
    int *degree, *start, *elist;
};
#endif 

