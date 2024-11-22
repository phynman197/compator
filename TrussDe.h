
#ifndef DMCE3_TRUSSDE_H
#define DMCE3_TRUSSDE_H
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <ctime>
#include <string>
#include <cstring>
#include <map>
#include <algorithm>
#include <cstring>
#include <vector>
#include "graph.h"
using namespace std;

typedef struct {
    int u,v;
} TEdge;
inline bool operator<(TEdge e1, TEdge e2) {
    return e1.u<e2.u || (e1.u==e2.u && e1.v<e2.v);
}
bool compVertex(int i, int j);


typedef map<int,int> MII;
typedef vector<int> VI;
typedef MII::iterator EdgeIter;const int maxClass=1000;

class TrussDe
{
private:
    ifstream fin;
    ofstream fout;
    string infile, outfile;

    int n, m;
    VI mapto;    VI bin;    vector<TEdge> binEdge;    vector<VI> A;    vector<MII> adj, pos;
    VI vertexorder;
    int cntClass[maxClass];



    inline void orderPair(int &u, int &v) {        if (!compVertex(u,v)) swap(u,v);
    }

    inline void printClass(int u, int v, int cls) {        ++cntClass[cls];
    }

    inline void updateSupport(int u, int v, int delta) {
        adj[u][v]+=delta;
        adj[v][u]+=delta;
    }

    inline void removeEdge(int u, int v) {
        adj[u].erase(v);
        adj[v].erase(u);
    }

    void readGraph();
    void readGraph2(VI & cand);
    void reorder();
    void intersect(const VI &a, const VI &b, VI &c);
    void countTriangles();
    void binSort();
    void updateEdge(int u, int v, int minsup);
    void trussDecomp();

public:
    int runTD(string, int*);

};
#endif 

