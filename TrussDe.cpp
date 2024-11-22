
#include "TrussDe.h"
#include "graph.h"
extern Graph g;

extern string datasetname;
VI deg_truss;
bool compVertex(int i, int j) {
    return deg_truss[i]<deg_truss[j] || (deg_truss[i]==deg_truss[j] && i<j);}

void TrussDe::readGraph()
{
    fin >> n >> m;    int vMax=0;    int u,v;    for (int i=0; i<m; ++i) {        fin >> u >> v;
        if (u==v) continue;
        vMax=max(vMax,max(u,v));
    }
    n=vMax+1;
    fin.close();
    fin.open(infile.c_str());
    int junk;
    fin>>junk>>junk;
    deg_truss.clear();
    deg_truss.resize(n,0);
    adj.resize(n);
    for (int i=0; i<n; ++i) adj[i].clear();
    for (int i=0; i<m; ++i) {
        fin >> u >> v;
        if (u==v) continue;
        if (adj[u].find(v)==adj[u].end()) {
            adj[u][v]=0;            adj[v][u]=0;
            ++deg_truss[u]; ++deg_truss[v];
        }
    }
}
void TrussDe::reorder() {    mapto.resize(n);
    for (int i=0; i<n; ++i) mapto[i]=i;
    sort(mapto.begin(), mapto.end(), compVertex);
}

void TrussDe::intersect(const VI &a, const VI &b, VI &c) {    c.clear();
    unsigned j=0;
    for (unsigned i=0; i<a.size(); ++i) {
        while (j<b.size() && b[j]>a[i]) ++j;
        if (j>=b.size()) break;
        if (b[j]==a[i]) c.push_back(a[i]);
    }
}

void TrussDe::countTriangles() {    A.resize(n);
    for (int i=0; i<n; ++i) A[i].clear();
    int nDeltas=0;
    for (int vi=n-1; vi>=0; --vi) {
        int v=mapto[vi];
        for (EdgeIter it = adj[v].begin(); it!=adj[v].end(); ++it) {
            int u = it->first;            if (!compVertex(u,v)) continue;            VI common;
            intersect(A[u], A[v], common);
            for (unsigned i=0; i<common.size(); ++i) {
                int w=mapto[common[i]];
                ++nDeltas;
                updateSupport(u,v,1);
                updateSupport(v,w,1);
                updateSupport(w,u,1);
            }
            A[u].push_back(vi);
        }
    }
}

void TrussDe::binSort() {
    bin.clear();
    bin.resize(n,0);
    int nBins=0;
    int mp=0;
    for (int u=0; u<n; ++u) {
        MII tadj = adj[u];
        for (EdgeIter it=tadj.begin(); it!=tadj.end(); ++it) {
            int v=it->first;
            if (!compVertex(u,v)) continue;
            int sup=it->second;            if (sup==0) {
                printClass(u,v,2);                removeEdge(u,v);
                                if (adj[u].size() == 0)
                    vertexorder.push_back(u);
                if (adj[v].size() == 0)
                    vertexorder.push_back(v);
                                continue;
            }
            ++mp;            ++bin[sup];            nBins=max(sup,nBins);        }
    }
    m=mp;    ++nBins;
    int count=0;
    for (int i=0; i<nBins; ++i) {
        int binsize=bin[i];
        bin[i]=count;        count+=binsize;
    }
    pos.clear();
    pos.resize(n);
    for (int i=0; i<n; ++i) pos[i].clear();
    binEdge.resize(m);
    for (int u=0; u<n; ++u)
        for (EdgeIter it=adj[u].begin(); it!=adj[u].end(); ++it) {
            int v=it->first;
            if (!compVertex(u,v)) continue;
            int sup=it->second;
            TEdge e={u,v};
            int &b=bin[sup];            binEdge[b]=e;            pos[u][v]=b++;        }
    for (int i=nBins; i>0; --i) bin[i]=bin[i-1];
    bin[0]=0;}

void TrussDe::updateEdge(int u, int v, int minsup) {
    orderPair(u,v);
    int sup=adj[u][v];
    if (sup<=minsup) return;
    int p=pos[u][v];    int posbin=bin[sup];    TEdge se=binEdge[posbin];
    TEdge e={u,v};
    if (p!=posbin) {
        pos[u][v]=posbin;
        pos[se.u][se.v]=p;
        binEdge[p]=se;
        binEdge[posbin]=e;
    }
    ++bin[sup];
    updateSupport(u,v,-1);
}

void TrussDe::trussDecomp() {
    for (int s=0; s<m; ++s) {        int u=binEdge[s].u;        int v=binEdge[s].v;
        orderPair(u,v);
        int supuv=adj[u][v];        printClass(u,v,supuv+2);
        int nfound=0;
        for (EdgeIter it=adj[u].begin(); it!=adj[u].end(); ++it) {
            if (nfound==supuv) break;
            int w=it->first;
            if (w==v) continue;
            if (adj[v].find(w)!=adj[v].end()) {
                ++nfound;
                updateEdge(u,w,supuv);
                updateEdge(v,w,supuv);
            }
        }
        removeEdge(u,v);
                if (adj[u].size() == 0)
            vertexorder.push_back(u);
        if (adj[v].size() == 0)
            vertexorder.push_back(v);
            }
}

int TrussDe::runTD(string dataname, int*lookup)
{
    clock_t start1, finish1;
    string dataset(dataname);
    infile=dataset+".txt";
    infile = datasetname;
        outfile=dataset+"-out.txt";
    fin.open(infile.c_str());
    fout.open(outfile.c_str());

    start1 = clock();
    readGraph();
    reorder();
        countTriangles();
    binSort();
        trussDecomp();

    for (int i=0;i<maxClass; ++i)
        if (cntClass[i]>0);

    fin.close();
    fout.close();

    finish1 = clock();
    double duration1 = (double)(finish1 - start1) / CLOCKS_PER_SEC;
    cout<<"duration: "<<duration1<<endl;
        int *test = new int[n];
    int testcnt = 0;
    for (int i = 0; i < n; i++)
    {
        test[i] = 0;
    }
    for (int i = 0; i < vertexorder.size(); i++)
    {
        test[vertexorder[i]] = 1;
        lookup[vertexorder[i]] = i;
    }
    for (int i = 0; i < n; i++)
    {
        if (test[i] == 0)
        {
            testcnt++;
        }
    }


    cout<<"Totally "<<testcnt<<" vertices are not sorted. "<<endl;

        return 0;
}

