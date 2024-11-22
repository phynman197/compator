//
// Created by Xiaofan Li on 10/3/21.
//

#ifndef DMCE3_CLIQUEHOLDER_H
#define DMCE3_CLIQUEHOLDER_H
#include <set>
#include <map>
#include "gadget.h"
#include "graph.h"

class CliqueHolder
{
public:
    CliqueHolder();
    bool insert(VI &C, char* function);
    int topk(int k);
    void init(double _ltau, bool _filter, FILE *_fout);
    int count() { return cliques.size();}
    VI largest() { return maxcli; }
    int average() { return sumsize/cliques.size(); }
    inline VI last() { return L; }
    inline double getTau() { return tau; }

private:
    FILE *fout;
    bool filter;
    long long cnt, sumsize;
    int ls;
    double tau;
    VVI cliques;
    VI maxcli, L, D;
    VI scores;
    multimap<int,int> lv;

    int updateScore(int i, VI &D);
    bool passFilter(VI &C);
    bool exactFilter(VI &C1, VI &C2, int t);
    int hook(const VI & C1, const VI & C2);
    int hookSize(const VI & C1, const VI & C2);
    int hook(const VI & C1, const VI & C2, int C1begin, int C1end, int C2begin, int C2end);
    int hookSize(const VI & C1, const VI & C2, int C1begin, int C1end, int C2begin, int C2end);
    int partialIntersection(const VI & C1, const VI & C2, int C1begin, int C1end, int C2begin, int C2end);
    int hamDis(const VI & C1, const VI & C2, int C1begin, int C1end, int C2begin, int C2end);
    int hamDisLB(const VI & A, const VI & B, int Abegin, int Aend, int Bbegin, int Bend, int maxAllow, int currentDep, int maxDep);
    int hamDisLBloose(const VI & A, const VI & B, int Abegin, int Aend, int Bbegin, int Bend);
    int binarySearchVertex(const VI & A, int Abegin, int Aend, int v);
    void printClique(const VI & C);
};
#endif //DMCE3_CLIQUEHOLDER_H

