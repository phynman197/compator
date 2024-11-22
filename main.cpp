#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <cstring>
#include <queue>
#include <algorithm>
#include <cmath>
#include <ctime>
#include "graph.h"
#include "gadget.h"
#include "cliqueholder.h"
#include "TemData.h"
#include "TrussDe.h"
#include<mach/mach.h>
int getmemory()
{
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;

    if (KERN_SUCCESS != task_info(mach_task_self(),
                                  TASK_BASIC_INFO, (task_info_t)&t_info,
                                  &t_info_count))
    {
        return -1;
    }
    return t_info.resident_size;
        }

using namespace std;
string datasetname;
void cliSearch2(VI &C, VI &cand, VI &prev);
void cliSearch3(VI &C, VI &cand, VI &prev, double pprod, char* bound, char* function);
int realfinalprunefunc();
int lessthan2func();
int summarysize(void);
int visi_lastdeletefunc();
int lsnumberfunc();
typedef vector<int> VI;
typedef vector<VI> VVI;
int stimes = 0;
int times0 = 0;
int times1 = 0;
int times2 = 0;
int times3 = 0;
int times4 = 0;
int dtimes = 0;
long noCmpCnt = 0, hsCmpCnt = 0, preCmpCnt = 0, posCmpCnt = 0, sufCmpCnt = 0, lfIntCnt = 0;
int MINSIZE = 3;
clock_t star,endt,frt;
double firstduration;


vector<string> split(const string& str, const string& delim) {
    vector<string> res;
    if("" == str) return res;
        char * strs = new char[str.length() + 1] ;     strcpy(strs, str.c_str());

    char * d = new char[delim.length() + 1];
    strcpy(d, delim.c_str());

    char *p = strtok(strs, d);
    while(p) {
        string s = p;         res.push_back(s);         p = strtok(NULL, d);
    }

    return res;
}

Graph g;
double tau;
bool randomized, filter, dotopk;
char est;
char order1;
int strategy;
int comparator;int maxDep = 0;

CliqueHolder holder;VI tmp, unit(1);
int *lookup;
int* sortedV;
int* staticlookup;
int ntofilter;VI L, D, E, X, hg;
FILE *fcliques;

int *vtable;
int *vCore;



double visi = 0;int ndelete = 0;

struct doubleint{
    int x;
    int y;
};

void outputarr(int* arr, int size)
{
    fstream fout;
    fout.open("/Users/xiaofanli/Desktop/outarray.txt", ios::app);
    for (int i = 0; i < size; i++)
    {
        fout<<arr[i]<<endl;
    }
    fout<<endl;
    fout.close();
}

bool LessSort (doubleint a,doubleint b) {
    return (a.y>b.y); }

int degeneracyOrdering(Graph & g, int par, char* order) {

    int n=g.V();    int *bin = (int*)alloc(n);
    int *pv = (int*)alloc(n);
    int *pos = (int*)alloc(n);
    int *tdeg = (int*)alloc(n);
    int* templookup = (int*)alloc(n);
    int*tempdegree = (int*)alloc(n);
    memset(bin, 0, sizeof(int)*n);     memset(tdeg, 0, sizeof(int)*n);

    int maxCoreNumber = 0;


    for (int i=0; i<n; ++i) tdeg[i] = g.deg(i);    for (int i=0; i<n; ++i) ++bin[tdeg[i]];
    int c=0;
    for (int i=0; i<n; ++i) {
        int d=bin[i];
        bin[i] = c;
        c += d;
    }    for (int i=0; i<n; ++i) {
        pos[i]=bin[tdeg[i]];        tempdegree[i] = pos[i];
        pv[pos[i]]=i;        ++bin[tdeg[i]];    }
    for (int i=n-1; i>0; --i)
        bin[i]=bin[i-1];
    bin[0]=0;
    for (int i=0; i<n; ++i) {
        int v=pv[i], dv=tdeg[v];        if(dv>maxCoreNumber){
            maxCoreNumber = dv;
                    }
        if(par==0){
            if (order[0] == 'I')
                lookup[v]=i;            if (order[0] == 'R')
                lookup[v]=v;
            if (order[0] == 'D')
                lookup[v]=n-1-i;
        }
        for (int j=0; j<g.deg(v); ++j) {
            int w=(g.adjvec(v))[j], dw=tdeg[w];
            if (dw>dv) {
                int pw=pos[w];
                int pt=bin[dw];
                int t=pv[pt];
                if (t!=w) {
                    pos[w]=pt, pos[t]=pw;
                    pv[pw]=t, pv[pt]=w;
                }
                ++bin[dw];
                --tdeg[w];
            }
        }
    }




    for(int i = 0; i < n; i++)
    {
        sortedV[lookup[i]] = i;
    }

    free(bin);
    free(pos);
    free(pv);
    if(par==0){
        for(int i=0;i<g.V();i++){
            vCore[i]=tdeg[i];
        }
    }

    free(tdeg);
    free(templookup);
    free(tempdegree);


    if(par==0){

        if (order[0] == 'T')
        {
            TrussDe TrussDecomposition;
            TrussDecomposition.runTD("no_use", lookup);
            for(int i = 0; i < n; i++)
            {
                sortedV[lookup[i]] = i;
            }
        }
    }

    return maxCoreNumber;

}


int hLocalDegreeDobund(VI & cand){

    int i=0;
    for(int v: cand){
        vtable[v]=i;
        ++i;
    }

        Graph localG(cand.size());
    int maxDeg = 0;
    for(int v: cand){
        VI nvs = g.adjvec(v);         VI nvsincand = cand*nvs;
        int vInCand = vtable[v];
        localG.setDeg(vInCand,nvsincand.size());
        if(localG.deg(vInCand)>maxDeg){
            maxDeg = localG.deg(vInCand);
        }
    }

    int *count =(int*)alloc(maxDeg);
    
    
    for(int v: cand){
        int vInCand = vtable[v];
        count[localG.deg(vInCand)]+=1;
    }


    int result = maxDeg;
    int acc=0;
    if(maxDeg>0){
        for(int j=maxDeg-1; j>=0;--j){
            acc+=count[j];
            if(acc>j+1){
                result =j+1;
                break;
            }
        }

    }
    free(count);
    return result+1;


}





int MaxDegree(VI & cand){

    int i=0;
    for(int v: cand){
        vtable[v]=i;
        ++i;
    }

        Graph localG(cand.size());
    int maxDeg = 0;
    for(int v: cand){
        VI nvs = g.adjvec(v);         VI nvsincand = cand*nvs;
        int vInCand = vtable[v];
        localG.setDeg(vInCand,nvsincand.size());
        if(localG.deg(vInCand)>maxDeg){
            maxDeg = localG.deg(vInCand);
        }
    }

    return maxDeg+1;

}



TemData coreBound(VI &cand){
    TemData t;
    int bound = 0;
                        int i=0;
        for(int v: cand){
                vtable[v]=i;
        ++i;
    }
    Graph localG(cand.size());
            int maxInt =0;
    int maxV = -1;
    for(int v: cand){
        VI nvs = g.adjvec(v);         VI nvsincand = cand*nvs;
        int vInCand = vtable[v];
        localG.setDeg(vInCand,nvsincand.size());
        if(maxInt<nvsincand.size()){
            maxInt = nvsincand.size();
            maxV = v;
        }
        for(int j:nvsincand){

            int jInCand = vtable[j];
            localG.addEdge(vInCand,jInCand);
                    }
    }

    char* a = new char[5];
    a[0] = 'R';
    t.bound = degeneracyOrdering(localG,1,a)+1;
    t.maxInters = maxInt;
    t.maxV = maxV;
            return t;
}

int  hCoreBound(VI & cand){

    int maxdep=0;
    for(int v: cand){
        if(maxdep<vCore[v]){
            maxdep = vCore[v];
        }
    }

    int *count =(int*)alloc(maxdep);
    
    
    for(int v: cand){
        count[vCore[v]]+=1;
    }


    int result = maxdep;

    if(maxdep>0){
        for(int j=maxdep-1; j>=0;--j){
            if(count[j]>j){
                result =j;
                break;
            }
        }

    }
    free(count);
    return result;
}


inline double NewprobSample(double r) {
    return r<tau? (tau-r)/(1-r) :0;
}

inline double OldprobSample(double r) {
    return (1.0-r)*(2.0-tau)/(2.0-r-tau);
}

inline double probKeep(double r, double pprod, int d, char* function) {    if (function[0]=='N')
        return pow(NewprobSample(r), 1.0/d);
    else
        return pow(OldprobSample(r), 1.0/d);



}

inline bool keepBranch(double pr) {
    return double(rand())/RAND_MAX < pr;
}








void cliSearch(VI &C, VI &cand, VI &prev, double pprod, char* bound, char* function, int d) {

        if (cand.empty() && prev.empty()) {
        ++ntofilter;
        holder.insert(C, function);
        return;
    }

    if (cand.empty() && !prev.empty()) {
        return;
    }

    
        int d00 = d;    int nc=cand.size();
    int maxdep;    L=holder.last();
    X=C*L;    E= cand-L;
    int cx=X.size(), ce=E.size(), cc=C.size();


    
        if (strategy != 0)    {
        stimes++;
        maxdep = d00;        
        if(strategy & 1){            if(((double)cx)/(cc+d00-1)>=tau)
            {
                                times1++;
                return;
            }
        }


        if(strategy & 2){            if (((double)(L.size()))/(cc+1)<tau)
            {
                times2++;
                cliSearch3(C, cand, prev, 1, bound, function);
                                return;
            }
        }

        if(strategy & 4){            if (((double)(cx+d00-1))/(cc+d00-1)<tau)
            {
                times3++;
                                cliSearch3(C, cand, prev, 1, bound, function);
                                return;
            }
        }
        if (strategy & 8)        {
            int ct=(L*cand).size();
            if (((double)(cx+ct))/(cc+1)<tau)            {
                times4++;
                cliSearch3(C, cand, prev, 1, bound, function);
                                return;
            }
        }



    }

    
        if (bound[0]=='L')
        maxdep = nc;


    if (bound[0]=='C')
    {
        TemData t = coreBound(cand);
        maxdep = t.bound;
    }



    if (bound[0]=='H')
    {
        maxdep = hLocalDegreeDobund(cand);

    }

    if (bound[0]=='M')
    {
        maxdep = MaxDegree(cand);
    }


    dtimes++;    
        double minr=1.0;
    for (int t=1; t<=maxdep; ++t) {
        double tr= double(max(t-ce,0)+cx)/(cc+t);
        minr = min(tr, minr);
    }
    double pr=1.0;
    if (!randomized) {
        if (minr >= holder.getTau())
            return;
        else
        {
                        int maxdeg=-1, pivot=-1;
            tmp=cand+prev;
            for (int i=0; i<tmp.size(); ++i) {
                int v=tmp[i];
                int dv=(g.adjvec(v)*cand).size();
                if (dv>maxdeg)
                    maxdeg=dv, pivot = v;
            }
            
            
            VI doing(cand-g.adjvec(pivot));
                                    for (int i=0; i<doing.size(); ++i) {
                int v = doing[i];
                unit[0]=v;                C= C+unit;
                cand= cand-unit;
                VI candt(cand*g.adjvec(v)), prevt(prev*g.adjvec(v));                cliSearch(C, candt, prevt, pprod*pr, bound, function, maxdep);
                unit[0]=v;
                C= C-unit;
                prev= prev+unit;
            }
        }
    }
    else{
        int maxl = maxdep+C.size();
        pr = probKeep(minr, pprod, maxl, function);

        if (!keepBranch(pr)){
            return;
        }

        else{
                        int maxdeg=-1, pivot=-1;
            tmp=cand+prev;
            for (int i=0; i<tmp.size(); ++i) {
                int v=tmp[i];
                int dv=(g.adjvec(v)*cand).size();
                if (dv>maxdeg)
                    maxdeg=dv, pivot = v;
            }
                        VI doing(cand-g.adjvec(pivot));
            for (int i=0; i<doing.size(); ++i) {
                int v = doing[i];
                unit[0]=v;                C= C+unit;
                cand= cand-unit;
                VI candt(cand*g.adjvec(v)), prevt(prev*g.adjvec(v));                cliSearch(C, candt, prevt, pprod*pr, bound, function, maxdep);
                unit[0]=v;
                C= C-unit;
                prev= prev+unit;
            }
        }
    }
}





void cliSearch3(VI &C, VI &cand, VI &prev, double pprod, char* bound, char* function) {

        if (cand.empty() && prev.empty()) {
        ++ntofilter;
        holder.insert(C, function);
        return;
    }
    if (cand.empty() && !prev.empty()) return;
    int maxdeg=-1, pivot=-1;
    tmp=cand+prev;
    for (int i=0; i<tmp.size(); ++i) {
        int v=tmp[i];
        int dv=(g.adjvec(v)*cand).size();
        if (dv>maxdeg)
            maxdeg=dv, pivot = v;
    }
    
    VI doing(cand-g.adjvec(pivot));
    for (int i=0; i<doing.size(); ++i) {
        int v = doing[i];
        unit[0]=v;        C= C+unit;
        cand= cand-unit;
        VI candt(cand*g.adjvec(v)), prevt(prev*g.adjvec(v));        cliSearch3(C, candt, prevt, 1, bound, function);
        unit[0]=v;
        C= C-unit;
        prev= prev+unit;
    }

}














int pritag(void);
int main(int argc, char** argv)
{



    srand( time(NULL) );


        
    if (argc<6) {
        printf("Usage:\t    1) graph_file\n\t2) tau\n\t3) R(andomized/D(eterministic)\n\t4) G(lobal)/L(ocal)\n\t5) output_file\n");
                exit(1);
    }
    char ses[100]="/Users/xiaofanli/GraphDataset/adjmatrix.txt";
        string temppath = argv[1];
    string path = "/";
    string orders(1,argv[6][0]);
    vector<string> vs;
    vs = split(temppath, "/");
    for (int i = 0; i < vs.size() - 1; i++)
    {
        path += vs[i];
        path += "/";
    }
    path = path + orders + "_";
    path += vs[vs.size() - 1];

            char pathc[100];
    strcpy(pathc, path.c_str());
    tau = atof(argv[2]);
    randomized = (argv[3][0]=='R');    filter = (argv[4][0]=='G');    fcliques = fopen(argv[5], "w");
    order1 = argv[6][0];
    strategy = (int) atof(argv[9]);    comparator = atof(argv[12]);
    maxDep = atof(argv[13]);
    cout<<"Dataset: "<<path<<endl;
    cout<<"Bound: "<<argv[7]<<endl;
    cout<<"Comparator: "<<argv[12]<<endl;
    cout<<"Max Depth: "<<argv[13]<<endl;
    cout<<"Tau: "<<tau<<endl;


    list2adj(pathc,ses);          datasetname = argv[1];

    g.build(ses);
    star = clock();
    vCore =  (int*)alloc(g.V());
    vtable = new int[g.V()];
    

    lookup = (int*)alloc(g.V());    sortedV = (int*)alloc(g.V());
    staticlookup = (int*)alloc(g.V());
    for (int i = 0; i < g.V(); i++)
    {
        lookup[i] = i;
        sortedV[i] = i;
        staticlookup[i] = i;
    }

        VI oc, odoing, ocand, oprev;    for (int i=0; i<g.V(); ++i) {
        odoing.push_back(sortedV[i]);
    }



        holder.init(tau, filter, fcliques);

    for (int ci=0; ci<odoing.size(); ++ci)
    {
                int i=odoing[ci];

        ocand.clear();
        oprev.clear();

        vector<doubleint> vint;





        for (int j=0; j<g.deg(i); ++j) {            int v= g.dest(i,j);

            if (lookup[i]<lookup[v]){
                doubleint t;
                t.x = v;
                t.y = lookup[v];
                vint.push_back(t);
                            }
            else
                oprev.push_back(v);
        }

        if (!vint.empty())
            sort(vint.begin(), vint.end(), LessSort);
                while (!vint.empty())
        {

            doubleint tt = *(vint.end()-1);
            ocand.push_back(tt.x);
            vint.pop_back();
        }



        oc.clear();
        oc.push_back(i);
        cliSearch(oc, ocand, oprev, 1.0, argv[7], argv[8], 99999);
    }
    double  duration;
    endt = clock();
    duration = (double)(endt - star) / CLOCKS_PER_SEC;









        cout<<"Summary Size: "<<summarysize()<<endl;


    cout<<"Duration: "<<duration<<endl;

    cout<<endl<<"Total Intersections: "<<noCmpCnt<<endl;
    cout<<"HookSize: "<<hsCmpCnt<<endl;
    cout<<"Prefix: "<<preCmpCnt<<endl;
    cout<<"Position: "<<posCmpCnt<<endl;
    cout<<"Sufix: "<<sufCmpCnt<<endl;
    cout<<"Left Intersection: "<<lfIntCnt<<endl<<endl;

    double visiout = double(visi+visi_lastdeletefunc()+summarysize())/(ndelete
                                                                       +realfinalprunefunc()+summarysize());



    string boundstr = argv[7];
    string functionstr = argv[8];
    string orderstr = argv[6];
    ofstream oFile;
    ofstream oFile3;
    ofstream oFileH;
    ofstream oFilePr;
    ofstream oFilePo;
    ofstream oFileSu;
    ofstream oFileI;

        if (orderstr == "d")
    {
        oFile.open("/Users/xiaofanli/Desktop/DMCE3x/output/st_"  + boundstr + "dd" + functionstr + argv[12]+ argv[13] + ".csv", ios::app );
        oFile3.open("/Users/xiaofanli/Desktop/DMCE3x/output/mt_"  + boundstr + "dd" + functionstr + argv[12]+ argv[13] + ".csv", ios::app );
        oFileH.open("/Users/xiaofanli/Desktop/DMCE3x/output/ht_"  + boundstr + "dd" + functionstr + argv[12]+ argv[13] + ".csv", ios::app );
        oFilePr.open("/Users/xiaofanli/Desktop/DMCE3x/output/prt_"  + boundstr + "dd" + functionstr + argv[12]+ argv[13] + ".csv", ios::app );
        oFilePo.open("/Users/xiaofanli/Desktop/DMCE3x/output/pot_"  + boundstr + "dd" + functionstr + argv[12]+ argv[13] + ".csv", ios::app );
        oFileSu.open("/Users/xiaofanli/Desktop/DMCE3x/output/sut_"  + boundstr + "dd" + functionstr + argv[12]+ argv[13] + ".csv", ios::app );
        oFileI.open("/Users/xiaofanli/Desktop/DMCE3x/output/it_"  + boundstr + "dd" + functionstr + argv[12]+ argv[13] + ".csv", ios::app );
    }
    else{
        oFile.open("/Users/xiaofanli/Desktop/DMCE3x/output/st_" + boundstr + orderstr + functionstr + argv[12]+ argv[13] + ".csv", ios::app );
        oFile3.open("/Users/xiaofanli/Desktop/DMCE3x/output/mt_" + boundstr + orderstr + functionstr + argv[12]+ argv[13] + ".csv", ios::app );
        oFileH.open("/Users/xiaofanli/Desktop/DMCE3x/output/ht_"  + boundstr + orderstr + functionstr + argv[12]+ argv[13] + ".csv", ios::app );
        oFilePr.open("/Users/xiaofanli/Desktop/DMCE3x/output/prt_"  + boundstr + orderstr + functionstr + argv[12]+ argv[13] + ".csv", ios::app );
        oFilePo.open("/Users/xiaofanli/Desktop/DMCE3x/output/pot_"  + boundstr + orderstr + functionstr + argv[12]+ argv[13] + ".csv", ios::app );
        oFileSu.open("/Users/xiaofanli/Desktop/DMCE3x/output/sut_"  + boundstr + orderstr + functionstr + argv[12]+ argv[13] + ".csv", ios::app );
        oFileI.open("/Users/xiaofanli/Desktop/DMCE3x/output/it_"  + boundstr + orderstr + functionstr + argv[12]+ argv[13] + ".csv", ios::app );
    }


            if (atof(argv[2])==0.1)
    {
        oFile<<argv[1]<<","<<argv[3]<<","<<argv[6]<<","<<argv[7]<<","<<argv[8]<<","<<argv[12]<<endl;
        oFile3<<argv[1]<<","<<argv[3]<<","<<argv[6]<<","<<argv[7]<<","<<argv[8]<<","<<argv[12]<<endl;
        oFileH<<argv[1]<<","<<argv[3]<<","<<argv[6]<<","<<argv[7]<<","<<argv[8]<<","<<argv[12]<<endl;
        oFilePr<<argv[1]<<","<<argv[3]<<","<<argv[6]<<","<<argv[7]<<","<<argv[8]<<","<<argv[12]<<endl;
        oFilePo<<argv[1]<<","<<argv[3]<<","<<argv[6]<<","<<argv[7]<<","<<argv[8]<<","<<argv[12]<<endl;
        oFileSu<<argv[1]<<","<<argv[3]<<","<<argv[6]<<","<<argv[7]<<","<<argv[8]<<","<<argv[12]<<endl;
        oFileI<<argv[1]<<","<<argv[3]<<","<<argv[6]<<","<<argv[7]<<","<<argv[8]<<","<<argv[12]<<endl;

        oFile<<"0.1"<<","<<"0.1"<<","<<"0.3"<<","<<"0.3"<<","<<"0.5"<<","<<"0.5"<<","<<"0.7"<<","<<"0.7"<<","<<"0.9"<<","<<"0.9"<<endl;
        oFile3<<"0.1"<<","<<"0.1"<<","<<"0.3"<<","<<"0.3"<<","<<"0.5"<<","<<"0.5"<<","<<"0.7"<<","<<"0.7"<<","<<"0.9"<<","<<"0.9"<<endl;
        oFileH<<"0.1"<<","<<"0.1"<<","<<"0.3"<<","<<"0.3"<<","<<"0.5"<<","<<"0.5"<<","<<"0.7"<<","<<"0.7"<<","<<"0.9"<<","<<"0.9"<<endl;
        oFilePr<<"0.1"<<","<<"0.1"<<","<<"0.3"<<","<<"0.3"<<","<<"0.5"<<","<<"0.5"<<","<<"0.7"<<","<<"0.7"<<","<<"0.9"<<","<<"0.9"<<endl;
        oFilePo<<"0.1"<<","<<"0.1"<<","<<"0.3"<<","<<"0.3"<<","<<"0.5"<<","<<"0.5"<<","<<"0.7"<<","<<"0.7"<<","<<"0.9"<<","<<"0.9"<<endl;
        oFileSu<<"0.1"<<","<<"0.1"<<","<<"0.3"<<","<<"0.3"<<","<<"0.5"<<","<<"0.5"<<","<<"0.7"<<","<<"0.7"<<","<<"0.9"<<","<<"0.9"<<endl;
        oFileI<<"0.1"<<","<<"0.1"<<","<<"0.3"<<","<<"0.3"<<","<<"0.5"<<","<<"0.5"<<","<<"0.7"<<","<<"0.7"<<","<<"0.9"<<","<<"0.9"<<endl;
    }

    oFile <<duration<<","<<summarysize();
    oFileH <<hsCmpCnt<<","<<summarysize();
    oFilePr <<preCmpCnt<<","<<summarysize();
    oFilePo <<posCmpCnt<<","<<summarysize();
    oFileSu <<sufCmpCnt<<","<<summarysize();
    oFileI <<lfIntCnt<<","<<summarysize();
    oFile3 <<firstduration<<","<<(double)getmemory()/1024/1024;    if (atof(argv[2])==0.9)
    {
        oFile<<endl<<endl;
        oFileH<<endl<<endl;
        oFilePr<<endl<<endl;
        oFilePo<<endl<<endl;
        oFileSu<<endl<<endl;
        oFileI<<endl<<endl;
        oFile3<<endl<<endl;
    }
    else{
        oFile<<",";
        oFileH<<",";
        oFilePr<<",";
        oFilePo<<",";
        oFileSu<<",";
        oFileI<<",";
        oFile3<<",";
    }
    oFile.close();
    oFile3.close();
    oFileH.close();
    oFilePr.close();
    oFilePo.close();
    oFileSu.close();
    oFileI.close();



    

    





    fclose(fcliques);
    free(lookup);
    free(sortedV);
    free(staticlookup);


    return 0;
}

