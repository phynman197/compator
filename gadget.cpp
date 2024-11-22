
#include "gadget.h"
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <cassert>
#include <algorithm>
#include <iostream>

typedef vector<int> VI;
using namespace std;
int* tag;
extern char order1;
bool cmps(VI &a,VI &b)
{
    return a.size()>b.size();
}

void *alloc(int n) {
    return malloc(sizeof(int)*(n+1));
}

VI operator+(VI &v1, VI &v2) {
    VI res(v1.size()+v2.size());
    sort(begin(v1), end(v1));    sort(begin(v2), end(v2));    int sz=set_union(v1.begin(), v1.end(), v2.begin(), v2.end(), res.begin()) - res.begin();
    res.resize(sz);
    return res;
};
VI operator-(VI &v1, VI &v2) {
    VI res(v1.size());
    sort(begin(v1), end(v1));    sort(begin(v2), end(v2));    int sz=set_difference(v1.begin(), v1.end(), v2.begin(), v2.end(), res.begin()) - res.begin();
    res.resize(sz);
    return res;

};
VI operator*(VI &v1, VI &v2) {
    VI res(v1.size());
    sort(begin(v1), end(v1));    sort(begin(v2), end(v2));    int sz=set_intersection(v1.begin(), v1.end(), v2.begin(), v2.end(), res.begin()) - res.begin();
    res.resize(sz);
    return res;
};

void vprint(VI v, FILE *fout) {
    for (int i=0; i<v.size(); ++i) {
        fprintf(fout, "%-8d", v[i]);
    }
    fprintf(fout, "\n");
}

void list2adj(char *ifname, char *ofname) {

    vector< VI > g;
    FILE *fi = fopen(ifname, "r");
    FILE *fo = fopen(ofname, "w");
    int n,m;
    fscanf(fi, "%d%d", &n, &m);
    g.resize(n);
    for (int i=0; i<m; ++i) {
        int v1, v2;
        fscanf(fi, "%d%d", &v1, &v2);
        g[v1].push_back(v2);
        g[v2].push_back(v1);
    }

        if (order1 != 'd')
    {
        fprintf(fo, "%d %d \n", n, m);
        for (int i=0; i<n; ++i)
        {
            sort(g[i].begin(), g[i].end());
            g[i].resize(unique(g[i].begin(), g[i].end()) - g[i].begin());
            int deg=g[i].size();
                        fprintf(fo, "%d %d ", i, deg);
            for(int j=0; j<deg; ++j)
                fprintf(fo, " %d ", g[i][j]);
            fprintf(fo, "\n");
        }
    }
        

    else{

        for(int i=0;i<n;i++)
        {
            g[i].resize(unique(g[i].begin(), g[i].end()) - g[i].begin());
            g[i].push_back(i);

        }
        sort(g.begin(),g.end(),cmps);
        tag = (int*)alloc(n);


        int x;
        for(int i=0;i<n;i++)
        {
            x=g[i].size();
            tag[g[i][x-1]]=i;            g[i].pop_back();

        }
        fprintf(fo, "%d %d \n", n, m);
        int deg;
        for (int i=0; i<n; ++i)
        {
            deg=g[i].size();
            for(int j=0; j<deg; ++j) g[i][j]=tag[g[i][j]];            sort(g[i].begin(), g[i].end());

            g[i].resize(unique(g[i].begin(), g[i].end()) - g[i].begin());
            deg=g[i].size();
                                                fprintf(fo, "%d %d", i, deg);
            for(int j=0; j<deg; ++j)
                fprintf(fo, " %d", g[i][j]);
            fprintf(fo, "\n");
        }
        delete tag;
    }
    
    cout<<"list2adj finished. "<<endl;
    fclose(fi);
    fclose(fo);
}

