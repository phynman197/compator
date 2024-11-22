
#include "cliqueholder.h"
#include "gadget.h"
#include <cstdio>
#include <algorithm>
#include <bitset>
#include <cassert>
#include <cmath>
#include <iostream>
using namespace std;

extern int maxDep;
extern int comparator;
extern long noCmpCnt;
extern long hsCmpCnt;
extern long preCmpCnt;
extern long posCmpCnt;
extern long sufCmpCnt;
extern long lfIntCnt;
extern int MINSIZE;
extern VI lookup;
extern Graph g;
extern clock_t frt;
extern clock_t star;
extern double firstduration;
double visi_lastdelete=0;

void CliqueHolder::printClique(const VI & C){
    for(vector<int>::const_iterator it = C.begin(); it != C.end(); ++it){
        cout<<(*it)<<", ";
    }
    cout<<endl;
}

int CliqueHolder::hook(const VI & C1, const VI & C2, int C1begin, int C1end, int C2begin, int C2end){
    if ((C1end - C1begin <= 0)||(C2end - C2begin <= 0)) return 0;
    if ((C1end<1)||(C2end<1)||(C1begin<0)||(C2begin<0)) return 0;
    return max(0, min(C1[C1end-1], C2[C2end-1])-max(C1[C1begin], C2[C2begin]));
}
int CliqueHolder::hookSize(const VI & C1, const VI & C2, int C1begin, int C1end, int C2begin, int C2end){
    return min(hook(C1, C2, C1begin, C1end, C2begin, C2end), (int) min(C1end-C1begin,C2end-C2begin));
}
int CliqueHolder::partialIntersection(const VI & C1, const VI & C2, int C1begin, int C1end, int C2begin, int C2end){
    int result = 0;
    int i = C1begin, j = C2begin;
    while(i != C1end && j != C2end){
        if(C1[i]<C2[j]) ++i;
        else if(C1[i]>C2[j]) ++j;
        else if(C1[i]==C2[j]){
            ++result;
            ++i; ++j;
        }
    }
    return result;
}
int CliqueHolder::hamDisLB(const VI & A, const VI & B, int Abegin, int Aend, int Bbegin, int Bend, int maxAllow, int currentDep, int maxDep){
    if(Aend-Abegin <= 0 || Bend-Bbegin <=0) return max(0,max(Aend-Abegin, Bend-Bbegin));
    int mid = (Aend - Abegin)/2;
    int v = A[mid + Abegin];
    if (Bbegin+mid-maxAllow >= Bbegin && Bbegin+mid+maxAllow <= Bend){
        if (B[Bbegin+mid-maxAllow] >= v || B[Bbegin+mid+maxAllow] <= v)
            return maxAllow+100;
    }
    int w = binarySearchVertex(B, Bbegin, Bend, v);
    int Dl = hamDisLBloose(A, B, Abegin, Abegin + mid, Bbegin, w);
    int Dr = hamDisLBloose(A, B, Abegin + mid + 1, Aend, w + 1, Bend);
    if (Dl + Dr > maxAllow) return maxAllow+100;
    if (currentDep++ >= maxDep) return Dl + Dr;
    Dl = hamDisLB(A, B, Abegin, Abegin + mid, Bbegin, w, maxAllow-Dr, currentDep, maxDep);
    if (Dl+Dr > maxAllow) return maxAllow+100;
    Dr = hamDisLB(A, B, Abegin + mid + 1, Aend, w + 1, Bend, maxAllow-Dl, currentDep, maxDep);
    if (Dl+Dr > maxAllow) return maxAllow+100;
    return Dl+Dr;
}
int CliqueHolder::hamDisLBloose(const VI & A, const VI & B, int Abegin, int Aend, int Bbegin, int Bend){
    return (Aend - Abegin) + (Bend - Bbegin) - 2*hookSize(A, B, Abegin, Aend, Bbegin, Bend);
}
int CliqueHolder::binarySearchVertex(const vector<int> & A, int Abegin, int Aend, int v){
    if (A[Aend-1] <= v) return Aend-1;
    if (A[Abegin] >= v) return Abegin;
    int mid = (Abegin + Aend - 1)/2;
    while(Abegin != Aend){
        if (A[mid]<v) {
            Abegin = mid + 1;
            mid = (Abegin + Aend - 1)/2;
        }
        if (A[mid] > v){
            Aend = mid;
            mid = (Abegin + Aend - 1)/2;
        }
        if(A[mid] == v) return mid;
    }
    return Abegin;
}

inline bool compv(int a, int b) {
    return lookup[a]<lookup[b];
}

CliqueHolder::CliqueHolder() {
}

void CliqueHolder::init(double _tau, bool _filter, FILE *_fout) {

    fout = _fout;
    cnt = sumsize = 0;
    filter = _filter;
    tau = _tau;
}

bool CliqueHolder::exactFilter(VI &C_, VI &C, int t) {    if(t == 0) return false;
    int accComVer = 0;
    int Cpend = C.size()-t+1;    int C_pend = C_.size()-t+1;    int iLastCom = 0;
    int jLastCom = 0;

    if(comparator == 0){        ++noCmpCnt;
        ++lfIntCnt;
        return  (C_*C).size()<=t;
                    }
    if(comparator & 1){        if(hookSize(C_,C, 0, C_.size(), 0, C.size())<t) {
            ++hsCmpCnt;
            return true;
        }
    }
    if(comparator & 2){
        if(hookSize(C, C_, 0, Cpend, 0, C_pend) == 0) {
            ++hsCmpCnt;
            return true;
        }
        int i=0, j=0;        while(i != Cpend && j != C_pend){
            if(C[i]<C_[j]) ++i;
            else if(C[i]>C_[j]) ++j;
            else if(C[i]==C_[j]){
                ++accComVer;
                if(accComVer+hookSize(C, C_, i+1, Cpend, j+1, C_pend)<t) {
                    ++posCmpCnt;
                    return true;
                }
                iLastCom = i++;
                jLastCom = j++;
            }
        }
        if(accComVer==0) {
            ++preCmpCnt;
            return true;
        }
    }
    if(comparator & 4){
                if(comparator & 2){
            int dMax = max(0, (int)((C.size() - 1 - iLastCom) + (C_.size() - 1 - jLastCom) + 2*(accComVer-t)));
            if (hamDisLB(C, C_, iLastCom+1, C.size(), jLastCom+1, C_.size(), dMax, 0, maxDep) > dMax) {
                ++sufCmpCnt;
                return true;
            }
        }

    }

    if(accComVer != 0) {
        ++lfIntCnt;
        return (accComVer+partialIntersection(C, C_, Cpend, C.size(), C_pend, C_.size()))<=t;
    }
    else {
        ++lfIntCnt;
        return (C_*C).size()<=t;
    }
}
int taggg=0,tagge=0;
bool CliqueHolder::passFilter(VI &C) {
    int upovl = (int)floor(C.size()*getTau());
    multimap<int,int>::iterator itlow, itup, it;    itlow= lv.begin();
        itup= lv.upper_bound(g.deg(D.front())-1);
    
    vector<multimap<int,int>::iterator> toerase;
    for (it=itlow; it!=itup; ++it)
        toerase.push_back(it);
    for (int i=0; i<toerase.size(); ++i)
        lv.erase(toerase[i]);
    toerase.clear();

    bool pass=true;
    int ncomp=0;
        for (it=itup; it!= lv.end(); ++it) {
            int i= (*it).second;
        ++ncomp;
        VI &cli=cliques[i];
        int j=cli.size()-1;

        
                
        VI P=cli;
        sort(P.begin(), P.end());
        taggg+=P.size()+C.size();
        if (!exactFilter(P, C, upovl)) {
            pass=false;
            break;
        }
    }

    for (int i=0; i<toerase.size(); ++i)
        lv.erase(toerase[i]);

    return pass;
}


double visibility=0.0;
int maxtuan;
int lsnumber=0;
int lessthan2 = 0;
int realfinalprune = 0;

bool CliqueHolder::insert(VI &C, char* function) {

        if (C.size()<=MINSIZE)
    {lessthan2++;
        lsnumber++;
        return false;}

        if (function[0]=='N')
    {
        int upovl = (int)floor(C.size()*getTau());
        if ((L*C).size()>=upovl){            realfinalprune++;
            lsnumber++;
            return false;}
    }
    if (function[0] == 'O')
    {
        double r,pr;
        r = ((double)(L*C).size()/C.size());
        pr = (double)(1-r)*(2-getTau())/(2-r-getTau());
        pr = pow(pr,1.0/C.size());

        if (( double(rand())/RAND_MAX) < 1-pr) return false;

    }
    taggg+=L.size()+C.size();
    D=C;
    sort(D.begin(), D.end());
    if (filter) {        if (!passFilter(C))
            return false;
    }
    tagge++;
    if (tagge == 1)    {
        frt = clock();
        firstduration = (double)(frt - star) / CLOCKS_PER_SEC;
        cout<<"First Result Time: "<<firstduration<<endl;
    }
    int id=cliques.size();
    cliques.push_back(D);
        lv.insert(make_pair(g.deg(D.back()),id));
    


    sumsize += C.size();
    L = C;
    if (C.size()> maxcli.size())
        maxcli = C;

    return true;
}

int CliqueHolder::updateScore(int i, VI &D) {

        VI &C = cliques[i];
    C= C-D;
    return C.size();
}

int CliqueHolder::topk(int k) {

        int quality=0;
    VI best;
    int nc=cliques.size();
    vector<bool> inans(nc, 0);
    scores.resize(nc, -1);
    for (int i=0; i<min(k,nc); ++i) {
        int scorebest=-1, ibest=-1;
        for (int j=0; j<nc; ++j)
            if (!inans[j]) {
                int score = updateScore(j, best);                 if (score>scorebest) {
                    scorebest = score;
                    ibest = j;
                }
            }
        quality+= scorebest;
        assert(scorebest>=0 && ibest>=0);
        best = cliques[ibest];
        inans[ibest]=1;
        vprint(best);
        printf("\n");
    }

    return quality;
}



int lsnumberfunc()
{
    return lsnumber;
}





int summarysize(void)
{
    return tagge;
}
int realfinalprunefunc()
{
    return realfinalprune;
}
int lessthan2func()
{
    return lessthan2;
}
int visi_lastdeletefunc()
{
    return visi_lastdelete;
}


