#include "oyd.dynamic_graph.h"
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <sys/time.h>
#include <tuple>
#include <vector>

#if defined(__i386__) || defined(__x86_64__)
#include <mmintrin.h>
#include <xmmintrin.h>
#endif

using namespace std;

namespace INDEX {
#define MAX_K 100
#define INT_MAX 999999999
#define RESERVE_TIME 1

// graph
clock_t ct;
int cnt, tree_width = 0;
const int INF = 999999999;
struct Graph {
    int n, m;
    vector<int> V;
    vector<map<int, int>> E;
    vector<vector<pair<int, int>>> Edge;
    vector<int> D;
    Graph() {
        n = m = 0;
        V.clear();
        E.clear();
    }
    Graph(int nN,  vector<tuple<int, int, int>> el) {
        Graph();
        n = nN;

        for (int i = 0; i <= n; i++) {
            map<int, int> v;
            v.clear();
            E.push_back(v);
        }

        for (int i = 0; i < el.size(); i++) {
            int x = std::get<0>(el[i]);
            int y = std::get<1>(el[i]);
            int z = std::get<2>(el[i]);

            if (x > n || y > n)
                continue;
            if (E[x].find(y) != E[x].end()) {
                if (E[x][y] > z) {
                    E[x][y] = z;
                    E[y][x] = z;
                }
            } else {
                E[x].insert(make_pair(y, z));
                E[y].insert(make_pair(x, z));
            }
        }
        D.clear();
        D.push_back(0);
        for (int j = 1; j <= n; j++)
            D.push_back(E[j].size());
    }

    void EdgeInitialize() {
        Edge.clear();
        for (int i = 0; i <= n; i++) {
            vector<pair<int, int>> Ed;
            Ed.clear();
            for (map<int, int>::iterator it = E[i].begin(); it != E[i].end();
                 it++) {
                Ed.push_back(*it);
            }
            Edge.push_back(Ed);
        }
    }
    bool isEdgeExist(int u, int v) {
        if (E[u].find(v) == E[u].end())
            return false;
        else
            return true;
    }
    void insertEdge(int u, int v, int k) {
        if (E[u].find(v) != E[u].end())
            return;
        E[u].insert(make_pair(v, k));
        E[v].insert(make_pair(u, k));
        D[u]++;
        D[v]++;
    }
    void deleteEdge(int u, int v) {
        if (E[u].find(v) == E[u].end())
            return;
        E[u].erase(E[u].find(v));
        E[v].erase(E[v].find(u));
        D[u]--;
        D[v]--;
    }
};

// tree decomposition
int *DD, *DD2, *NUM;
int *_DD, *_DD2;
bool *changed;
struct SelEle {
    int x;
    SelEle();
    SelEle(int _x) { x = _x; }
    bool operator<(const SelEle se) const {
        if (DD[x] != DD[se.x])
            return DD[x] < DD[se.x];
        if (DD2[x] != DD2[se.x])
            return DD2[x] < DD2[se.x];
        return x < se.x;
    }
};
struct Node {
    vector<int> vert, VL, pos, pos2, dis;
    vector<int> ch;
    int height;
    int pa;
    int uniqueVertex;
    Node() {
        vert.clear();
        VL.clear();
        pos.clear();
        dis.clear();
        ch.clear();
        pa = -1;
        uniqueVertex = -1;
        height = 0;
    }
};

class Tree_Decomposition {
  public:
    FILE *fout, *fin;
    Graph G, H;
    set<SelEle> deg;
    int maxSize;
    Tree_Decomposition() {}
    vector<vector<int>> neighbor, length;
    vector<int> ord;
    int heightMax;
    void reduce() {
        deg.clear();
        neighbor.clear();
        length.clear();
        vector<int> vectmp;
        vectmp.clear();
        for (int i = 0; i <= G.n; i++) {
            neighbor.push_back(vectmp);
            length.push_back(vectmp);
        }
        DD = (int *)malloc(sizeof(int) * (G.n + 1));
        DD2 = (int *)malloc(sizeof(int) * (G.n + 1));
        _DD = (int *)malloc(sizeof(int) * (G.n + 1));
        _DD2 = (int *)malloc(sizeof(int) * (G.n + 1));
        NUM = (int *)malloc(sizeof(int) * (G.n + 1));
        changed = (bool *)malloc(sizeof(bool) * (G.n + 1));
        for (int i = 0; i <= G.n; i++)
            NUM[i] = i;
        for (int i = 1; i <= G.n; i++) {
            int j = rand() % G.n + 1;
            int x = NUM[i];
            NUM[i] = NUM[j];
            NUM[j] = x;
        }
        for (int i = 1; i <= G.n; i++) {
            DD[i] = G.D[i];
            DD2[i] = G.D[i];
            _DD[i] = G.D[i];
            _DD2[i] = G.D[i];
            changed[i] = false;
            deg.insert(SelEle(i));
        }
        bool *exist;
        exist = (bool *)malloc(sizeof(bool) * (G.n + 1));
        for (int i = 1; i <= G.n; i++)
            exist[i] = true;
        ord.clear();
        int cnt = 0;
        while (!deg.empty()) {
            cnt++;
            int x = (*deg.begin()).x;
            while (true) {
                if (changed[x]) {
                    deg.erase(SelEle(x));
                    DD[x] = _DD[x];
                    DD2[x] = _DD2[x];
                    deg.insert(SelEle(x));
                    changed[x] = false;
                    x = (*deg.begin()).x;
                } else
                    break;
            }
            ord.push_back(x);
            deg.erase(deg.begin());
            exist[x] = false;
            vector<int> neigh, leng;
            neigh.clear();
            leng.clear();
            for (map<int, int>::iterator it = G.E[x].begin();
                 it != G.E[x].end(); it++) {
                int y = (*it).first;
                if (exist[y]) {
                    neigh.push_back(y);
                    leng.push_back((*it).second);
                }
            }
            int k = -1;
            for (int i = 0; i < neigh.size(); i++) {
                int y = neigh[i];
                G.deleteEdge(x, y);
                _DD[y] = G.D[y];
                changed[y] = true;
            }
            for (int pu = 0; pu < neigh.size(); pu++) {
                for (int pv = 0; pv < neigh.size(); pv++)
                    if (pu != pv) {
                        int u = neigh[pu], v = neigh[pv];
                        if (G.isEdgeExist(u, v)) {
                            if (G.E[u][v] > leng[pu] + leng[pv])
                                G.E[u][v] = leng[pu] + leng[pv];
                            if (G.E[v][u] > leng[pu] + leng[pv])
                                G.E[v][u] = leng[pu] + leng[pv];

                        } else {
                            G.insertEdge(u, v, leng[pu] + leng[pv]);
                            _DD[u] = G.D[u];
                            _DD[v] = G.D[v];
                            ++_DD2[u];
                            ++_DD2[v];
                            changed[u] = true;
                            changed[v] = true;
                        }
                    }
            }
            if (neigh.size() > tree_width)
                tree_width = neigh.size();
            neighbor[x] = neigh;
            length[x] = leng;
        }
        free(DD);
        free(DD2);
        free(exist);
    }
    int match(int x, vector<int> &neigh) {
        int nearest = neigh[0];
        for (int i = 1; i < neigh.size(); i++)
            if (rank[neigh[i]] > rank[nearest])
                nearest = neigh[i];
        int p = belong[nearest];
        vector<int> a = Tree[p].vert;
        if (Tree[p].uniqueVertex >= 0) {
            a.push_back(Tree[p].uniqueVertex);
        }
        sort(a.begin(), a.end());
        int i, j = 0;
        for (; (i < neigh.size()) && (j < a.size());) {
            if (neigh[i] == a[j]) {
                i++;
                j++;
            } else if (neigh[i] < a[j])
                break;
            else
                j++;
        }
        if (i >= neigh.size()) {
            return p;
        }
        printf("no match!\n");
    }
    vector<Node> Tree;
    int root;
    int *belong, *rank;
    void makeTree() {
        belong = (int *)malloc(sizeof(int) * (H.n + 1));
        rank = (int *)malloc(sizeof(int) * (H.n + 1));
        int len = ord.size() - 1;
        Node rootn;
        Tree.clear();
        heightMax = 0;

        int x = ord[len];
        rootn.vert = neighbor[x];
        rootn.VL = length[x];
        rootn.uniqueVertex = x;
        rootn.pa = -1;
        rootn.height = 1;
        rank[x] = 0;
        belong[x] = 0;
        Tree.push_back(rootn);
        len--;

        for (; len >= 0; len--) {
            int x = ord[len];
            Node nod;
            nod.vert = neighbor[x];
            nod.VL = length[x];
            nod.uniqueVertex = x;
            int pa = match(x, neighbor[x]);
            Tree[pa].ch.push_back(Tree.size());
            nod.pa = pa;
            nod.height = Tree[pa].height + 1;
            if (nod.height > heightMax)
                heightMax = nod.height;
            rank[x] = Tree.size();
            belong[x] = Tree.size();
            Tree.push_back(nod);
        }

        root = 0;
    }
    struct PT {
        int dis;
        int x;
        PT() {}
        PT(int _dis, int _x) {
            dis = _dis;
            x = _x;
        }

        bool operator<(const PT _pt) const {
            if (dis == _pt.dis)
                return x > _pt.x;
            return dis > _pt.dis;
        }
    };
    void calculateDistanceAll(int start, int *a) {
        for (int i = 0; i <= G.n; i++)
            a[i] = INF;
        a[start] = 0;
        priority_queue<PT> Q;
        //	heap Q;
        while (!Q.empty())
            Q.pop();
        Q.push(PT(0, start));
        while (!Q.empty()) {
            int distance = Q.top().dis;
            int x = Q.top().x;
            Q.pop();
            if (distance > a[x])
                continue;
            for (int i = 0; i < H.Edge[x].size(); i++) {
                int y = H.Edge[x][i].first, z = H.Edge[x][i].second;
                if (a[y] > distance + z) {
                    a[y] = distance + z;
                    Q.push(PT(a[y], y));
                }
            }
        }
    }

    int *toRMQ;
    vector<int> EulerSeq;
    vector<vector<int>> RMQIndex;
    void makeRMQDFS(int p, int height) {
        toRMQ[p] = EulerSeq.size();
        EulerSeq.push_back(p);
        for (int i = 0; i < Tree[p].ch.size(); i++) {
            makeRMQDFS(Tree[p].ch[i], height + 1);
            EulerSeq.push_back(p);
        }
    }
    void makeRMQ() {
        EulerSeq.clear();
        toRMQ = (int *)malloc(sizeof(int) * (G.n + 1));
        makeRMQDFS(root, 1);
        RMQIndex.clear();
        RMQIndex.push_back(EulerSeq);
        int m = EulerSeq.size();
        for (int i = 2, k = 1; i < m; i = i * 2, k++) {
            vector<int> tmp;
            tmp.clear();
            tmp.resize(EulerSeq.size());
            for (int j = 0; j < m - i; j++) {
                int x = RMQIndex[k - 1][j], y = RMQIndex[k - 1][j + i / 2];
                if (Tree[x].height < Tree[y].height)
                    tmp[j] = x;
                else
                    tmp[j] = y;
            }
            RMQIndex.push_back(tmp);
        }
    }
    int LCAQuery(int _p, int _q) {
        int p = toRMQ[_p], q = toRMQ[_q];
        if (p > q) {
            int x = p;
            p = q;
            q = x;
        }
        int len = q - p + 1;
        int i = 1, k = 0;
        while (i * 2 < len) {
            i *= 2;
            k++;
        }
        q = q - i + 1;
        if (Tree[RMQIndex[k][p]].height < Tree[RMQIndex[k][q]].height)
            return RMQIndex[k][p];
        else
            return RMQIndex[k][q];
    }
    int distanceQueryAncestorToPosterity(int p, int q) {
        if (p == q)
            return 0;
        int x = belong[p], y = belong[q];
        return Tree[y].dis[Tree[x].pos[Tree[x].pos.size() - 1]];
    }
    void calculateIndexSizeDFS(int p, int pre, int tmp, long long &result) {
        for (int i = 0; i < Tree[p].ch.size(); i++) {
            calculateIndexSizeDFS(Tree[p].ch[i], pre + 1, (pre + 1) + tmp,
                                  result);
        }
        if (tmp + (pre + 1) > result)
            result = tmp + (pre + 1);
        //		result += pre;
    }
    long long calculateIndexSize() {
        long long res = Tree[root].vert.size();
        for (int i = 0; i < Tree[root].ch.size(); i++) {
            calculateIndexSizeDFS(Tree[root].ch[i], Tree[root].vert.size(),
                                  Tree[root].vert.size(), res);
        }
        return res;
    }
    void makeIndexDFS(int p, vector<int> &list, int *toList) {
        //	cnt++;
        //	cout << endl << cnt << endl;
        Tree[p].pos.resize(Tree[p].vert.size() + 1);
        Tree[p].dis.resize(list.size());
        //	printf("step1");
        for (int i = 0; i < Tree[p].vert.size(); i++) {
            int j;
            for (j = 0; j < list.size(); j++)
                if (list[j] == Tree[p].vert[i])
                    break;
            Tree[p].pos[i] = j;
        }
        Tree[p].pos[Tree[p].vert.size()] = list.size();
        for (int i = 0; i < list.size(); i++) {
            Tree[p].dis[i] = INF;
        }
        priority_queue<PT> Q;
        while (!Q.empty())
            Q.pop();
        int x = Tree[p].uniqueVertex;

        for (int i = 0; i < Tree[p].vert.size(); i++) {
            if (Tree[p].dis[toList[Tree[p].vert[i]]] > Tree[p].VL[i])
                Tree[p].dis[toList[Tree[p].vert[i]]] = Tree[p].VL[i];
            int x = Tree[p].vert[i];
            int k;
            for (k = 0; k < list.size(); k++)
                if (list[k] == x)
                    break;
            for (int j = 0; j < list.size(); j++) {
                int y = list[j];
                int z;
                if (k < j)
                    z = distanceQueryAncestorToPosterity(x, y);
                else
                    z = distanceQueryAncestorToPosterity(y, x);
                if (Tree[p].dis[toList[y]] >
                    z + Tree[p].dis[toList[Tree[p].vert[i]]])
                    Tree[p].dis[toList[y]] =
                        z + Tree[p].dis[toList[Tree[p].vert[i]]];
            }
        }
        //	ct += clock() - t;
        //	printf("step4");
        toList[Tree[p].uniqueVertex] = list.size();
        list.push_back(Tree[p].uniqueVertex);
        for (int i = 0; i < Tree[p].ch.size(); i++) {
            makeIndexDFS(Tree[p].ch[i], list, toList);
        }
        list.pop_back();
        //	printf("step5");
        //---
        Tree[p].pos2 = Tree[p].pos;
        for (int i = Tree[p].vert.size() - 1; i >= 0; i--) {
            if (Tree[p].VL[i] > Tree[p].dis[Tree[p].pos2[i]]) {
                Tree[p].pos2.erase(Tree[p].pos2.begin() + i);
            }
        }
        //---

        sort(Tree[p].pos.begin(), Tree[p].pos.end());
        sort(Tree[p].pos2.begin(), Tree[p].pos2.end());

        fwrite(&p, SIZEOFINT, 1, fout);
        //
        Tree[p].vert.push_back(Tree[p].uniqueVertex);
        printIntVector(Tree[p].vert);
        Tree[p].vert.pop_back();
        printIntVector(Tree[p].dis);

        vector<int> v1, v2;
        v1.clear();
        v2.clear();
        //	Tree[p].pos.swap(v1);
        Tree[p].dis.swap(v2);
    }
    void makeIndex() {
        makeRMQ();

        H.EdgeInitialize();
    }
    void makeIndex2() {
        vector<int> list;
        list.clear();
        int *toList;
        toList = (int *)malloc(sizeof(int) * (G.n + 1));
        Tree[root].pos.clear();

        // dijstra ---- map -> vector
        toList[Tree[root].uniqueVertex] = 0;
        list.push_back(Tree[root].uniqueVertex);
        Tree[root].pos.push_back(0);

        fwrite(&root, SIZEOFINT, 1, fout);
        Tree[root].vert.push_back(Tree[root].uniqueVertex);
        printIntVector(Tree[root].vert);
        Tree[root].vert.pop_back();
        printIntVector(Tree[root].dis);
        cnt = 0;
        for (int i = 0; i < Tree[root].ch.size(); i++)
            makeIndexDFS(Tree[root].ch[i], list, toList);
        free(toList);
    }
    void reducePosDFS(int p) {
        //----
        if (Tree[p].ch.size() == 2) {
            int t = Tree[p].ch[0];
            if (Tree[Tree[p].ch[0]].pos.size() > Tree[Tree[p].ch[1]].pos.size())
                t = Tree[p].ch[1];
            Tree[p].pos = Tree[t].pos;
            Tree[p].pos.erase(Tree[p].pos.begin() + (Tree[p].pos.size() - 1));
        }
        //----
        for (int i = 0; i < Tree[p].ch.size(); i++)
            reducePosDFS(Tree[p].ch[i]);
        fwrite(&p, SIZEOFINT, 1, fout);
        printIntVector(Tree[p].pos);
    }
    void reducePos() { reducePosDFS(root); }

    void cntSize() {
        long long s_nonroot = 0;
        long long s_size = 0;

        long long s_dis = 0;
        for (int i = 0; i < Tree.size(); ++i) {
            s_nonroot += Tree[i].height - 1;
            s_size += Tree[i].vert.size();
            s_dis += Tree[i].height;
        }
        long long s_root = (long long)Tree[0].vert.size() * (long long)G.n;
        // printf("tree width: %d\n", tree_width);
        // printf("nonroot idx size = %0.3lfGB, avg node size=%0.3lf, avg label
        // "
        //        "size=%0.3lf\n",
        //        s_nonroot * 4.0 / 1000000000.0, s_size * 1.0 / G.n,
        //        s_dis * 1.0 / G.n);
    }
    static const int SIZEOFINT = 4;
    void printIntArray(int *a, int n) { fwrite(a, SIZEOFINT, n, fout); }
    void printIntVector(vector<int> &a) {
        if (a.size() == 0) {
            int x = 0;
            fwrite(&x, SIZEOFINT, 1, fout);
            return;
        }
        int x = a.size();
        fwrite(&x, SIZEOFINT, 1, fout);
        for (int i = 0; i < a.size(); i++) {
            fwrite(&a[i], SIZEOFINT, 1, fout);
        }
    }
    void print() {
        // G.n
        fwrite(&G.n, SIZEOFINT, 1, fout);
        // printf("G.n %d\n", G.n);
        //  Tree.size() Tree.height
        int x = Tree.size();
        fwrite(&x, SIZEOFINT, 1, fout);
        // printf("Tree.size(): %d\n", x);
        for (int i = 0; i < Tree.size(); i++) {
            fwrite(&Tree[i].height, SIZEOFINT, 1, fout);
        }
        for (int i = 0; i < Tree.size(); i++) {
            fwrite(&Tree[i].pa, SIZEOFINT, 1, fout);
        }
        for (int i = 0; i < Tree.size(); i++) {
            fwrite(&Tree[i].uniqueVertex, SIZEOFINT, 1, fout);
        }
        // belong
        printIntArray(belong, H.n + 1);
        // LCA - toRMQ - RMQIndex
        printIntArray(toRMQ, H.n + 1);
        x = RMQIndex.size();
        fwrite(&x, SIZEOFINT, 1, fout);
        // printf("RMWIndex.size(): %d\n", RMQIndex.size());
        x = EulerSeq.size();
        fwrite(&x, SIZEOFINT, 1, fout);
        // printf("EulerSeq.size(): %d\n", EulerSeq.size());
        for (int i = 0; i < RMQIndex.size(); i++)
            printIntVector(RMQIndex[i]);
        // rootDistance
        fwrite(&root, SIZEOFINT, 1, fout);
        // printf("root: %d\n", root);
        //
        for (int i = 0; i < Tree.size(); i++) {
            int t = Tree[i].ch.size();
            fwrite(&t, SIZEOFINT, 1, fout);
            for (int j = 0; j < t; j++)
                fwrite(&Tree[i].ch[j], SIZEOFINT, 1, fout);
        }
        //
    }
};
} // namespace INDEX

namespace QUERY {
#define MAX_K 100
#define INT_MAX 999999999
#define RESERVE_TIME 1

int reset_times = 0;

const int infinity = 999999999;
const int SIZEOFINT = 4;

double cnt_pre_query_time = 0;
char *insert_type, *query_type;

int *toRMQ, *height, *pa, *uniqueVertex, **RMQIndex;
int *belong;
int root, TreeSize;
int **rootToRoot, *rootSite;
int **dis, **pos, **pos2;
int *posSize, *pos2Size;
int *chSize;
int **ch;
int *LOG2, *LOGD;
int rootSize;
int *DFSList, *toDFS;
int ***BS;

struct HASH_NODE {
    int a, b, c;
    int next;
    HASH_NODE() {}
    HASH_NODE(int _a, int _b, int _c) {
        a = _a;
        b = _b;
        c = _c;
    }
};
class HASH {
  public:
    static const int P = 100000007;

    vector<HASH_NODE> nodes;
    vector<int> start;
    HASH() {
        nodes.clear();
        nodes.push_back(HASH_NODE());
        start.resize(P);
        for (int i = 0; i < P; i++)
            start[i] = 0;
    }
    inline int get_id(int a, int b) { return ((long long)(a)*b) % P; }
    inline bool is_exist(int a, int b) {
        int t = get_id(a, b);
        int p = start[t];
        while (p != 0) {
            if (nodes[p].a == a && nodes[p].b == b)
                return true;
            p = nodes[p].next;
        }
        return false;
    }
    inline int get_value(int a, int b) {

        int t = get_id(a, b);

        int p = start[t];

        while (p != 0) {
            if (nodes[p].a == a && nodes[p].b == b)
                return nodes[p].c;
            p = nodes[p].next;
        }
        return -1;
    }
    inline bool insert_node(int a, int b, int c) {

        if (is_exist(a, b) == true)
            return false;
        int t = get_id(a, b);
        HASH_NODE hn = HASH_NODE(a, b, c);
        hn.next = start[t];
        nodes.push_back(hn);
        start[t] = nodes.size() - 1;
        return true;
    }
    inline bool delete_node(int a, int b) {
        int t = get_id(a, b);
        int p = start[t];
        int pre = -1;
        while (p != 0) {
            if (nodes[p].a == a && nodes[p].b == b) {
                if (pre < 0)
                    start[t] = nodes[p].next;
                else
                    nodes[pre].next = nodes[p].next;
                return true;
            }
            pre = p;
            p = nodes[p].next;
        }
        return false;
    }
};

inline int LCAQuery(int _p, int _q) {
    int p = toRMQ[_p], q = toRMQ[_q];
    //	printf("p q : %d %d\n", p, q);
    if (p > q) {
        int x = p;
        p = q;
        q = x;
    }
    int len = q - p + 1;
    //	printf("len: %d\n", len);
    int i = LOGD[len], k = LOG2[len];
    //	printf("i, k: %d %d\n", i, k);
    q = q - i + 1;
    if (height[RMQIndex[k][p]] < height[RMQIndex[k][q]])
        return RMQIndex[k][p];
    else
        return RMQIndex[k][q];
}

long long queryCnt;
// long long aCnt;

inline int distanceQuery(int p, int q) {
    if (p == q)
        return 0;
    int x = belong[p], y = belong[q];

    if (toRMQ[x] > toRMQ[y]) {
        int k = x;
        x = y;
        y = k;
    }
    return dis[y][height[x] - 1];
}

int *Degree;
int **Neighbor, **Weight;
void readGraph(char *filename) {
    FILE *file = fopen(filename, "r");
    int n, m;
    fscanf(file, "%d %d", &n, &m);
    Degree = (int *)malloc(sizeof(int) * (n + 1));
    vector<vector<pair<int, int>>> nb;
    vector<pair<int, int>> v;
    v.clear();
    for (int i = 0; i <= n; i++) {
        //	Degree[i] = 0;
        nb.push_back(v);
    }
    for (int i = 0; i < m; i++) {
        int x, y, z;
        fscanf(file, "%d %d %d", &x, &y, &z);
        nb[x].push_back(make_pair(y, z));
    }
    Neighbor = (int **)malloc(sizeof(int *) * (n + 1));
    Weight = (int **)malloc(sizeof(int *) * (n + 1));
    for (int i = 1; i <= n; i++) {
        Degree[i] = nb[i].size();
        Neighbor[i] = (int *)malloc(sizeof(int) * nb[i].size());
        Weight[i] = (int *)malloc(sizeof(int) * nb[i].size());
        for (int j = 0; j < nb[i].size(); j++) {
            Neighbor[i][j] = nb[i][j].first;
            Weight[i][j] = nb[i][j].second;
        }
    }
}
inline int shortestPathQuery(int p, int q) {
    int res = 0;
    while (p != q) {
        res++;
        int pq = distanceQuery(p, q);
        for (int i = 0; i < Degree[p]; i++) {
            int x = Neighbor[p][i];
            //	int y = Weight[p][i];
            int xq = distanceQuery(x, q);
            if (xq + Weight[p][i] == pq) {
                p = x;
                break;
            }
        }
    }
    return res;
}

int test(char *file) {
    //	cout << "test: " << file << " BEGIN" << endl;
    FILE *fin = fopen(file, "r");
    int x, y, dis, res = 0;
    vector<int> X, Y, DIS;
    X.clear();
    Y.clear();
    DIS.clear();
    while (fscanf(fin, "%d %d %d", &x, &y, &dis) != EOF) {
        if (x <= 0 || y <= 0)
            break;
        X.push_back(x);
        Y.push_back(y);
        DIS.push_back(dis);
        res++;
    }
    cout << X.size() << endl;

    //
    queryCnt = 0;
    //	aCnt = 0;
    //
    double t = 0;
    int *ANS;
    ANS = (int *)malloc(sizeof(int) * (X.size() + 1));
    double start_time = GetTime();
    for (int i = 0; i < X.size(); i++) {
        ANS[i] = distanceQuery(X[i], Y[i]);
        //	if (ANS[i] != DIS[i]){
        //		cout << X[i] << " " << Y[i] << " " << DIS[i] << " " <<
        // ANS[i] << endl; 		while(1);
        //	}
    }
    double end_time = GetTime();
    for (int i = 0; i < X.size(); i++) {
        t += ANS[i];
    }
    cout << "Check Count: " << double(queryCnt) / res << endl;
    //	cout << "aCnt: " << double(aCnt) / res << endl;
    printf("Distance Query Time : %lf usec\n",
           (end_time - start_time) * 1e6 / res);
    printf("average distance: %.6lf\n", t / res);
    return res;
}

FILE *fin;
string TT = "";
void scanIntArray(int *a, int n) { fread(a, SIZEOFINT, n, fin); }
int *scanIntVector(int *a) {
    int _n;
    fread(&_n, SIZEOFINT, 1, fin);
    a = (int *)malloc(sizeof(int) * _n);
    scanIntArray(a, _n);
    return a;
}

int n;
int *EulerSeq;
void readIndex(const char *file) {
    double _time = GetTime();
    int tree_height = 0, tree_width = 0, most_sp = 0;
    fin = fopen(file, "rb");
    fread(&n, SIZEOFINT, 1, fin);
    //	printf("n: %d\n", n);
    int ts;
    fread(&ts, SIZEOFINT, 1, fin);
    //	printf("ts: %d\n", ts);
    TreeSize = ts;
    height = (int *)malloc(sizeof(int) * (ts + 1));
    pa = (int *)malloc(sizeof(int) * (ts + 1));
    uniqueVertex = (int *)malloc(sizeof(int) * (ts + 1));
    for (int i = 0; i < ts; i++) {
        fread(&height[i], SIZEOFINT, 1, fin);
    }
    for (int i = 0; i < ts; i++) {
        fread(&pa[i], SIZEOFINT, 1, fin);
    }
    for (int i = 0; i < ts; i++) {
        fread(&uniqueVertex[i], SIZEOFINT, 1, fin);
    }
    belong = (int *)malloc(sizeof(int) * (n + 1));
    fread(belong, SIZEOFINT, n + 1, fin);
    toRMQ = (int *)malloc(sizeof(int) * (n + 1));
    fread(toRMQ, SIZEOFINT, n + 1, fin);
    int ris;
    fread(&ris, SIZEOFINT, 1, fin);
    //	printf("ris: %d\n", ris);
    fread(&ts, SIZEOFINT, 1, fin);
    //	printf("ts: %d\n", ts);
    EulerSeq = (int *)malloc(sizeof(int) * (ts + 1));
    RMQIndex = (int **)malloc(sizeof(int *) * (ris + 1));
    for (int i = 0; i < ris; i++) {
        RMQIndex[i] = scanIntVector(RMQIndex[i]);
    }
    fread(&root, SIZEOFINT, 1, fin);
    //	cout << "root: " << root << endl;

    posSize = (int *)malloc(sizeof(int) * (n + 1));
    pos2Size = (int *)malloc(sizeof(int) * (n + 1));
    pos = (int **)malloc(sizeof(int *) * (TreeSize));
    pos2 = (int **)malloc(sizeof(int *) * (TreeSize));
    dis = (int **)malloc(sizeof(int *) * (TreeSize));
    chSize = (int *)malloc(sizeof(int) * (TreeSize));
    ch = (int **)malloc(sizeof(int *) * (TreeSize));

    for (int i = 0; i < TreeSize; i++) {
        fread(&chSize[i], SIZEOFINT, 1, fin);
        ch[i] = (int *)malloc(sizeof(int) * chSize[i]);
        for (int j = 0; j < chSize[i]; j++) {
            int x;
            fread(&x, SIZEOFINT, 1, fin);
            ch[i][j] = x;
        }
    }
    for (int i = 0; i < TreeSize; i++) {
        int x;
        fread(&x, SIZEOFINT, 1, fin);
        fread(&posSize[x], SIZEOFINT, 1, fin);
        pos[x] = (int *)malloc(sizeof(int) * (posSize[x] + 1));
        fread(pos[x], SIZEOFINT, posSize[x], fin);
        if (posSize[x] > tree_width)
            tree_width = posSize[x];
        int _n;
        fread(&_n, SIZEOFINT, 1, fin);
        dis[x] = (int *)malloc(sizeof(int) * _n);
        fread(dis[x], SIZEOFINT, _n, fin);
        if (_n > tree_height)
            tree_height = _n;
    }
    //	printf("dis read finished!\n");
    for (int i = 0; i < TreeSize; i++) {
        int x;
        fread(&x, SIZEOFINT, 1, fin);
        fread(&pos2Size[x], SIZEOFINT, 1, fin);
        pos2[x] = (int *)malloc(sizeof(int) * (pos2Size[x] + 1));
        fread(pos2[x], SIZEOFINT, pos2Size[x], fin);
        if (pos2Size[x] > most_sp)
            most_sp = pos2Size[x];
    }

    fclose(fin);
}
int cnt;
void getDFSListDFS(int p) {
    toDFS[p] = cnt;
    DFSList[cnt++] = p;
    for (int i = 0; i < chSize[p]; i++) {
        getDFSListDFS(ch[p][i]);
    }
    BS[p] = (int **)malloc(sizeof(int *) * chSize[p]);
    for (int i = 0; i < chSize[p]; i++) {
        BS[p][i] = (int *)malloc(sizeof(int) * chSize[p]);
        for (int j = 0; j < chSize[p]; j++) {
            if (posSize[ch[p][i]] < posSize[ch[p][j]])
                BS[p][i][j] = ch[p][i];
            else
                BS[p][i][j] = ch[p][j];
        }
    }
}
void getDFSList() {
    DFSList = (int *)malloc(sizeof(int) * (TreeSize + 1));
    toDFS = (int *)malloc(sizeof(int) * (TreeSize + 1));
    BS = (int ***)malloc(sizeof(int **) * (TreeSize + 1));
    cnt = 0;
    getDFSListDFS(root);
}

int *tmp_dis, *higher, *cloest_higher;
bool anc_compare(int p, int q) {
    if (tmp_dis[p] < tmp_dis[q])
        return true;
    return false;
}

int *is_current_object, *current_distance, *group_height, *current_state,
    current_stamp;

class kNN {
  public:
    struct list_node {
        int previous, next, key, dist;
        list_node() {
            previous = -1;
            next = -1;
            key = -1;  // vertex
            dist = -1; // distance
        }
    };

    struct object_saveing_structure {

        vector<list_node> a;
        vector<int> trash_can;
        int current, size_num;
        object_saveing_structure() {
            a.clear();
            list_node _a;
            _a.key = -1;
            _a.dist = -1;
            _a.previous = 0;
            _a.next = 0;
            a.push_back(_a);
            trash_can.clear();
            current = 0;
            size_num = 0;
        }
    };
    HASH hash;
    vector<double> times_period[10];
    int period = 8;
    vector<object_saveing_structure> OSS;
    vector<int> object_number;
    vector<vector<pair<int, int>>> path_from_root;
    //	set<int> object_set;
    void create_kNN_index() {
        // /	object_set.clear();
        vector<pair<int, int>> _v;
        _v.clear();
        for (int i = 0; i <= TreeSize; i++)
            path_from_root.push_back(_v);
        OSS.clear();
        for (int i = 0; i < TreeSize; i++) {
            object_saveing_structure oss;
            OSS.push_back(oss);
        }
        // printf("TreeSize: %d\n", TreeSize);
    }
    void delete_element(int p, int x) {
        OSS[p].size_num--;
        int pre = OSS[p].a[x].previous;
        int ne = OSS[p].a[x].next;
        OSS[p].a[ne].previous = pre;
        OSS[p].a[pre].next = ne;
    }
    inline void insert(int p, int x) {
        cout << "(--" << p << ", " << x << "--)";

        int y = uniqueVertex[p];
        int disxy = distanceQuery(y, x);
        if (OSS[p].size_num >= int(MAX_K))
            if (OSS[p].a[OSS[p].a[0].previous].dist <= disxy)
                return;

        int i = 0;

        while (OSS[p].a[i].previous != 0 &&
               OSS[p].a[OSS[p].a[i].previous].dist > disxy)
            i = OSS[p].a[i].previous;
        list_node _a;
        _a.next = i;
        _a.previous = OSS[p].a[i].previous;
        _a.key = x;
        _a.dist = disxy;

        OSS[p].a.push_back(_a);
        //
        hash.insert_node(p, x, OSS[p].a.size() - 1);
        //
        OSS[p].a[_a.previous].next = OSS[p].a.size() - 1;
        OSS[p].a[_a.next].previous = OSS[p].a.size() - 1;
        OSS[p].size_num++;
        if (OSS[p].size_num > MAX_K * RESERVE_TIME)
            delete_element(p, OSS[p].a[0].previous);
    }
    inline void insert(int x) {

        if (is_current_object[x] == 1)
            return;
        is_current_object[x] = 1;

        int p = belong[x];

        while (p >= 0) {
            //			printf("insert p, x: %d %d\n", p, x);
            object_number[p]++;
            insert(p, x);
            //		list.push_back(make_pair(p, OSS[p].a.size() - 1));
            p = pa[p];
        }
    }
    void get_all_object(int p, vector<int> &a) {
        int x = uniqueVertex[p];
        if (is_current_object[x]) {
            a.push_back(x);
        }
        for (int i = 0; i < chSize[p]; i++) {
            get_all_object(ch[p][i], a);
        }
    }
    void delete_object_sort(int p, int x) {}
    void delete_object(int p, int x) {
        //   cout << "delete object" << endl;
        int y = uniqueVertex[p];
        int disxy = distanceQuery(y, x);

        hash.delete_node(p, x);
        //
        if (OSS[p].a[OSS[p].a[0].previous].dist < disxy) {
            //           cout << "larger than top_k" << endl;
            return;
        }
        int i = 0;
        while (OSS[p].a[i].next != 0 &&
               OSS[p].a[OSS[p].a[i].next].dist <= disxy) {
            i = OSS[p].a[i].next;
            if (OSS[p].a[i].key == x)
                break;
        }

        if (OSS[p].a[i].key == x) {
            delete_element(p, i);

        } else {
            return;
        }
        if (OSS[p].size_num >= MAX_K) {
            return;
        }

        if (object_number[p] <= OSS[p].size_num) {
            if (object_number[p] < OSS[p].size_num)
                while (1)
                    ;
            return;
        }

        reset_times++;
        for (int i = OSS[p].a[0].next; i != 0; i = OSS[p].a[i].next)
            delete_element(p, i);
        vector<int> a;
        a.clear();
        if (is_current_object[x] == 1)
            OSS_push_front(p, x, 0);
        //    cout << insert_type << endl;
        if (strcmp(insert_type, "sort") == 0) {
            //      cout << "hehe" << endl;
            vector<int> a;
            a.clear();
            get_all_object(p, a);
            for (int i = 0; i < a.size(); i++)
                current_distance[a[i]] = distanceQuery(y, a[i]);
            int _MAX = int(MAX_K * RESERVE_TIME);
            int current_k;
            if (_MAX < a.size()) {
                nth_element(a.begin(), a.begin() + _MAX, a.end(),
                            object_compare);
                current_k = _MAX;
            } else {
                current_k = a.size();
            }
            sort(a.begin(), a.end(), object_compare);
            //  for (int j = OSS[p].a[0].next; j != 0; j = OSS[p].a[j].next )
            //    delete_element(p, j);
            for (int i = 0; i < current_k; i++) {
                OSS_push_back(p, a[i], current_distance[a[i]]);
            }

        } else {
            get_subtree(p, uniqueVertex[p], a);
            join_subtree(p, uniqueVertex[p], a);
        }
    }
    inline void delete_object(int x) {

        if (is_current_object[x] == 0)
            return;
        is_current_object[x] = 0;

        int p = belong[x];
        //	vector<pair<int, int> > list;
        //	list.clear();
        while (p >= 0) {
            object_number[p]--;
            delete_object(p, x);

            //		list.push_back(make_pair(p, OSS[p].a.size() - 1));
            p = pa[p];
        }
    }

    static bool object_compare(int a, int b) {
        if (current_distance[a] < current_distance[b])
            return true;
        else
            return false;
    }

    void dfs_merge(int p, vector<int> &a) {}
    void OSS_push_back(int p, int key, int dist) {
        //	printf("p, key, dist: %d %d %d\n", p, key, dist);
        list_node _a;
        _a.previous = OSS[p].a[0].previous;
        _a.next = 0;
        _a.key = key;
        _a.dist = dist;
        OSS[p].a.push_back(_a);
        OSS[p].a[_a.previous].next = OSS[p].a.size() - 1;
        OSS[p].a[_a.next].previous = OSS[p].a.size() - 1;
        OSS[p].size_num++;
    }
    void OSS_push_front(int p, int key, int dist) {
        list_node _a;
        _a.previous = 0;
        _a.next = OSS[p].a[0].next;
        _a.key = key;
        _a.dist = dist;
        OSS[p].a.push_back(_a);
        OSS[p].a[_a.previous].next = OSS[p].a.size() - 1;
        OSS[p].a[_a.next].previous = OSS[p].a.size() - 1;
        OSS[p].size_num++;
    }
    void get_subtree(int p, int key, vector<int> &a) {
        //	printf("p, key, a.size(): %d %d %d\n", p, key, a.size());
        bool is_containing = false;
        int x = uniqueVertex[p];
        if (x == key)
            is_containing = true;
        else {
            //		printf("posSize[p]: %d\n", posSize[p]);
            for (int i = 0; i < posSize[p]; i++) {
                //			printf("pos[p][i]: %d\n", pos[p][i]);
                if (pos[p][i] == key) {
                    is_containing = true;
                    break;
                }
            }
        }
        //	printf("is_containing: %d\n", is_containing);
        if (is_containing) {
            a.push_back(p);

            //		printf("a.size(), chSize[p]: %d %d\n", a.size(),
            // chSize[p]);
            for (int i = 0; i < chSize[p]; i++)
                get_subtree(ch[p][i], key, a);
        }
    }
    void dfs_sort(int p, vector<int> &a) {
        int _MAX = int(MAX_K * RESERVE_TIME);
        vector<int> b;
        a.clear();
        int x = uniqueVertex[p];
        if (is_current_object[x] == 1)
            a.push_back(x);
        for (int i = 0; i < chSize[p]; i++) {
            dfs_sort(ch[p][i], b);
            a.insert(a.end(), b.begin(), b.end());
        }
        int current_k;
        for (int i = 0; i < a.size(); i++)
            current_distance[a[i]] = distanceQuery(x, a[i]);
        if (_MAX < a.size()) {
            nth_element(a.begin(), a.begin() + _MAX, a.end(), object_compare);
            current_k = _MAX;
        } else {
            current_k = a.size();
        }
        //	printf("a.size(), current_k: %d %d\n", a.size(), current_k);
        // sort(a.begin(), a.begin() + current_k, object_compare);
        sort(a.begin(), a.end(), object_compare);

        //	if (OSS[p].size_num > MAX_K)
        //		delete_element(p, OSS[p].a[0].previous);
        if (p == root) {
            //		for (int i = 0; i < a.size(); i++)
            //			printf("%d ", a[i]);
            //		printf("\n");
        }
        for (int i = 0; i < current_k; i++) {
            OSS_push_back(p, a[i], current_distance[a[i]]);
        }
        if (!double_objects(p)) {
            print(p);
            stop();
        }
        //	printf("p, a.size()):%d %d\n", p, a.size());
    }
    void stop() {
        while (1)
            ;
    }
    void join_subtree(int p, int x, vector<int> &a) {

        priority_queue<PT> Q;
        while (!Q.empty())
            Q.pop();
        vector<int> iter;
        iter.clear();
        iter.push_back(0);
        current_stamp++;
        for (int i = 1; i < a.size(); i++) {
            iter.push_back(OSS[a[i]].a[0].next);
            int t;
            while (iter[i] != 0) {
                t = OSS[a[i]].a[iter[i]].key;
                //	cout << "a[i]: " << a[i] << ", iter[i]: " << iter[i] <<
                // endl; 	cout << "t: " << t << " "  << current_state[t]
                // << " " << current_stamp << " " << current_distance[t] << " "
                // << distanceQuery(t, x) << endl;;
                if (current_state[t] != current_stamp) {
                    current_state[t] = current_stamp;
                    current_distance[t] = distanceQuery(t, x);
                    break;
                } else
                    iter[i] = OSS[a[i]].a[iter[i]].next;
            }
            if (iter[i] != 0) {
                Q.push(PT(current_distance[t], i));
                //		cout << "push: " << p << ", " << a[i] << ", " <<
                // current_distance[t] << ", " << OSS[a[i]].a[iter[i]].dist <<
                // endl;
            }
        }
        vector<int> b;
        b.clear();
        for (int i = 0; i < int(RESERVE_TIME * MAX_K); i++) {
            if (!Q.empty()) {
                PT pt = Q.top();
                Q.pop();
                //	OSS_push_back(p, OSS[a[pt.x]].a[iter[pt.x]].key,
                // pt.dis);
                b.push_back(OSS[a[pt.x]].a[iter[pt.x]].key);
                //		cout << "pop: " << p << ", " <<
                // OSS[a[pt.x]].a[iter[pt.x]].key << ", " << pt.dis << endl;
                int t;
                while (iter[pt.x] != 0) {
                    t = OSS[a[pt.x]].a[iter[pt.x]].key;

                    if (current_state[t] != current_stamp) {
                        current_state[t] = current_stamp;
                        current_distance[t] = distanceQuery(t, x);
                        break;
                    } else
                        iter[pt.x] = OSS[a[pt.x]].a[iter[pt.x]].next;
                }
                if (iter[pt.x] != 0) {
                    Q.push(PT(current_distance[t], pt.x));
                }
                //	else cout << pt.x << " " << a[pt.x] << " end " << endl;
            } else
                break;
        }
        //	cout << endl;
        sort(b.begin(), b.end(), object_compare);
        for (int i = 0; i < b.size(); i++) {
            OSS_push_back(p, b[i], current_distance[b[i]]);
            hash.insert_node(p, b[i], OSS[p].a.size() - 1);
        }
    }
    void dfs_neighbor(int p) {
        //    printf("dfs_neighbor: %d\n", p);
        vector<int> a;
        a.clear();
        int x = uniqueVertex[p];
        if (is_current_object[x] == 1) {
            //           printf("%d is object\n", x);
            OSS_push_front(p, x, 0);
        }
        for (int i = 0; i < chSize[p]; i++) {
            dfs_neighbor(ch[p][i]);
        }
        get_subtree(p, x, a);

        join_subtree(p, x, a);

        if (!double_objects(p)) {
            print(p);
            stop();
        }
    }
    void traversal(int p) {

        int x = uniqueVertex[p];
        if (is_current_object[x] == 1)
            object_number[p]++;
        for (int i = 0; i < chSize[p]; i++) {
            traversal(ch[p][i]);
            object_number[p] += object_number[ch[p][i]];
        }
    }
    void compute_object_number() {
        object_number.resize(n + 1);
        for (int i = 0; i <= n; i++)
            object_number[i] = 0;
        traversal(root);
    }
    void initialize_object() {
        compute_object_number();
        // printf("start initialize:\n");
        // cout << "insert_type: " << insert_type << endl;
        //	traversal(root);
        //  type 1 sort
        vector<int> a;
        for (int i = 0; i < n; i++)
            for (int j = OSS[i].a[0].next; j != 0; j = OSS[i].a[j].next)
                delete_element(i, j);
        //
        //	dfs_merge(root, a);
        if (strcmp(insert_type, "sort") == 0)
            dfs_sort(root, a);
        else
            dfs_neighbor(root);
        // cout << "initialize object finished" << endl;
    }
    int *query_mark;
    int query_mark_stamp;
    vector<double> time_save;
    vector<pair<int, int>> query(int x, int top_k) {
        time_save.clear();
        double pre_time = GetTime();
        double start_time = pre_time;
        query_mark_stamp++;
        int cnt1 = 0, cnt2 = 0, cnt3 = 0;
        //	printf("path_from_root[1].size(): %d\n",
        // path_from_root[1].size());
        vector<pair<int, int>> result;
        result.clear();
        int p = belong[x];

        vector<int> anc;
        anc.clear();

        int MAX_DIS = INT_MAX;

        vector<int> a;
        a.clear();
        a.push_back(p);
        for (int i = 0; i < posSize[p]; i++)
            a.push_back(belong[pos[p][i]]);

        for (int i = 0; i < a.size(); i++) {
            int q = a[i];
            if (OSS[q].size_num <= top_k)
                continue;
            tmp_dis[q] = distanceQuery(uniqueVertex[q], x);
            if (tmp_dis[q] >= MAX_DIS)
                continue;
            if (OSS[q].a[0].next == 0)
                continue;
            int j = OSS[q].a[0].next;
            int _cnt = 1;
            while (_cnt < top_k) {
                _cnt++;
                j = OSS[q].a[j].next;
                cnt2++;
            }
            if (MAX_DIS > tmp_dis[q] + OSS[q].a[j].dist)
                MAX_DIS = tmp_dis[q] + OSS[q].a[j].dist;
        }

        int max_height = height[p];
        while (p >= 0 && height[p] >= max_height) {
            cnt3++;
            int t = OSS[p].a[0].next;
            tmp_dis[p] = distanceQuery(uniqueVertex[p], x);
            if (t > 0 && MAX_DIS >= tmp_dis[p] + OSS[p].a[t].dist) {
                anc.push_back(p);
                OSS[p].current = OSS[p].a[0].next;
                if (OSS[p].size_num >= top_k) {
                    int v = OSS[p].a[0].previous;
                    int tmp = tmp_dis[p] + OSS[p].a[v].dist;
                    if (tmp < MAX_DIS)
                        MAX_DIS = tmp;
                }
            }
            if (tmp_dis[p] + cloest_higher[p] < MAX_DIS)
                if (max_height > higher[p]) {
                    max_height = higher[p];
                }
            p = pa[p];
        }
        sort(anc.begin(), anc.end(), anc_compare);
        //	printf("anc.size(): %d\n", anc.size());
        p = belong[x];
        cnt_pre_query_time += GetTime() - pre_time;
        int _cnt = 0;
        for (int i = 0; i < top_k; i++) {
            //	printf("i result.size(): %d %d\n", i, result.size());
            //--------------find minimum-------------
            int k = -1, dist_k = -1;
            for (int j = 0; j < anc.size(); j++) {
                _cnt++;
                int q = anc[j];
                //	printf("q: %d\n", q);
                while (OSS[q].current != 0 &&
                       query_mark[OSS[q].a[OSS[q].current].key] ==
                           query_mark_stamp)
                    OSS[q].current = OSS[q].a[OSS[q].current].next;
                if (OSS[q].current == 0)
                    continue;
                int _dist = tmp_dis[q] + OSS[q].a[OSS[q].current].dist;
                if (k < 0 || (dist_k > _dist)) {
                    k = q;
                    dist_k = _dist;
                }
            }
            if (k < 0)
                break;
            result.push_back(make_pair(OSS[k].a[OSS[k].current].key, dist_k));

            //--------------delete and update ------------------

            int y = OSS[k].a[OSS[k].current].key;
            query_mark[y] = query_mark_stamp;
            double current_time = GetTime();
            time_save.push_back(current_time - pre_time);
            if ((i + 1) % 5 == 0)
                times_period[((i + 1) / 5) - 1].push_back(current_time -
                                                          start_time);
            //   printf("---%.6lf\n", current_time - pre_time);
            pre_time = current_time;
        }
        //   cout << "_cnt: " << _cnt << endl;

        return result;
    }

    vector<pair<int, int>> query_naive(int x, int top_k) {
        time_save.clear();
        double pre_time = GetTime();
        query_mark_stamp++;
        int cnt1 = 0, cnt2 = 0, cnt3 = 0;
        vector<pair<int, int>> result;
        result.clear();
        int p = belong[x];

        vector<int> anc;
        anc.clear();
        int MAX_DIS = INT_MAX;
        vector<int> a;
        a.clear();
        a.push_back(p);
        for (int i = 0; i < posSize[p]; i++)
            a.push_back(belong[pos[p][i]]);
        for (int i = 0; i < a.size(); i++) {
            int q = a[i];
            if (OSS[q].size_num <= top_k)
                continue;
            tmp_dis[q] = distanceQuery(uniqueVertex[q], x);
            if (tmp_dis[q] >= MAX_DIS)
                continue;
            if (OSS[q].a[0].next == 0)
                continue;
            int j = OSS[q].a[0].next;
            int _cnt = 1;
            while (_cnt < top_k) {
                _cnt++;
                j = OSS[q].a[j].next;
                cnt2++;
            }
            if (MAX_DIS > tmp_dis[q] + OSS[q].a[j].dist)
                MAX_DIS = tmp_dis[q] + OSS[q].a[j].dist;
        }

        int max_height = height[p];
        while (p >= 0 && height[p] >= max_height) {
            cnt3++;
            int t = OSS[p].a[0].next;
            tmp_dis[p] = distanceQuery(uniqueVertex[p], x);
            if (t > 0 && MAX_DIS >= tmp_dis[p] + OSS[p].a[t].dist) {
                anc.push_back(p);
                OSS[p].current = OSS[p].a[0].next;
                if (OSS[p].size_num >= top_k) {
                    int v = OSS[p].a[0].previous;
                    int tmp = tmp_dis[p] + OSS[p].a[v].dist;
                    if (tmp < MAX_DIS)
                        MAX_DIS = tmp;
                }
            }
            if (tmp_dis[p] + cloest_higher[p] < MAX_DIS)
                if (max_height > higher[p]) {
                    max_height = higher[p];
                }
            p = pa[p];
        }
        sort(anc.begin(), anc.end(), anc_compare);
        p = belong[x];
        cnt_pre_query_time += GetTime() - pre_time;

        //--------------find minimum-------------
        //    int k = -1, dist_k = -1;

        a.clear();
        for (int j = 0; j < anc.size(); j++) {
            cnt1++;
            int q = anc[j];
            //    printf("q: %d\n", q);
            int i = 0;
            OSS[q].current = OSS[q].a[0].next;
            while (OSS[q].current != 0 && i < top_k) {
                i++;
                int t_t = tmp_dis[q] + OSS[q].a[OSS[q].current].dist;

                int v = OSS[q].a[OSS[q].current].key;

                if (query_mark[OSS[q].a[OSS[q].current].key] ==
                    query_mark_stamp) {
                    if (t_t < tmp_dis[v])
                        tmp_dis[v] = t_t;
                } else {
                    query_mark[v] = query_mark_stamp;
                    a.push_back(v);
                    tmp_dis[v] = t_t;
                }
                OSS[q].current = OSS[q].a[OSS[q].current].next;
            }
        }

        sort(a.begin(), a.end(), anc_compare);
        for (int i = 0; i < a.size(); i++) {
            if (i >= top_k)
                break;
            result.push_back(make_pair(a[i], tmp_dis[a[i]]));
        }

        double current_time = GetTime();
        time_save.push_back(current_time - pre_time);

        return result;
    }
    vector<pair<int, int>> query_delay(int x, int top_k) {
        time_save.clear();
        double pre_time = GetTime();
        double start_time = pre_time;
        vector<pair<int, int>> result;
        result.clear();
        int p = belong[x];

        vector<int> anc;
        anc.clear();

        int MAX_DIS = INT_MAX;
        vector<int> a;
        a.clear();
        a.push_back(p);

        //   cout << "step2"<< endl;
        for (int i = 0; i < posSize[p]; i++)
            a.push_back(belong[pos[p][i]]);
        for (int i = 0; i < a.size(); i++) {
            int q = a[i];
            if (OSS[q].size_num <= top_k)
                continue;
            tmp_dis[q] = distanceQuery(uniqueVertex[q], x);
            if (tmp_dis[q] >= MAX_DIS)
                continue;
            if (OSS[q].a[0].next == 0)
                continue;
            int j = OSS[q].a[0].next;
            int _cnt = 1;
            while (j != 0 && _cnt < top_k) {
                j = OSS[q].a[j].next;
            }
            if (MAX_DIS > tmp_dis[q] + OSS[q].a[j].dist)
                MAX_DIS = tmp_dis[q] + OSS[q].a[j].dist;
        }

        //   cout << "step3"<< endl;
        int max_height = height[p];
        while (p >= 0 && height[p] >= max_height) {
            int t = OSS[p].a[0].next;
            tmp_dis[p] = distanceQuery(uniqueVertex[p], x);
            if (t > 0 && MAX_DIS >= tmp_dis[p] + OSS[p].a[t].dist) {
                anc.push_back(p);
                OSS[p].current = OSS[p].a[0].next;
                if (OSS[p].size_num >= top_k) {
                    int v = OSS[p].a[0].previous;
                    int tmp = tmp_dis[p] + OSS[p].a[v].dist;
                    if (tmp < MAX_DIS)
                        MAX_DIS = tmp;
                }
            }
            if (tmp_dis[p] + cloest_higher[p] < MAX_DIS)
                if (max_height > higher[p]) {
                    max_height = higher[p];
                }
            p = pa[p];
        }
        //     printf("anc.size(): %d\n", anc.size());
        //    printf("OSS[7745].current: %d\n", OSS[7745].current);
        sort(anc.begin(), anc.end(), anc_compare);

        //      cout << "step4"<< endl;
        p = belong[x];
        int _cnt = 0;
        for (int i = 0; i < top_k; i++) {

            //--------------find minimum-------------
            int k = -1, dist_k = -1;
            for (int j = 0; j < anc.size(); j++) {
                _cnt++;
                //      cnt1++;
                int q = anc[j];
                //        printf("q: %d\n", q);
                if (OSS[q].current == 0)
                    continue;
                int _dist = distanceQuery(uniqueVertex[q], x) +
                            OSS[q].a[OSS[q].current].dist;
                if (k < 0 || (dist_k > _dist)) {
                    k = q;
                    dist_k = _dist;
                }
            }
            //    printf("k dist_k: %d %d\n", k, dist_k);
            if (k < 0)
                break;
            int y = OSS[k].a[OSS[k].current].key;
            //   cout << y << " " << dist_k << " " << k << endl;
            result.push_back(make_pair(y, dist_k));

            //--------------delete and update ------------------

            //-----------------------------------

            for (int j = 0; j < anc.size(); j++) {
                _cnt++;
                int t = anc[j];

                if (OSS[t].a[OSS[t].current].key == y)
                    OSS[t].current = OSS[t].a[OSS[t].current].next;
                int pos = hash.get_value(t, y);
                //         if (t == 1476 && y == 30137)
                //            cout << t << " " << pos << endl;

                if (pos >= 0) {
                    OSS[t].a[OSS[t].a[pos].previous].next = OSS[t].a[pos].next;
                    OSS[t].a[OSS[t].a[pos].next].previous =
                        OSS[t].a[pos].previous;
                    OSS[t].trash_can.push_back(pos);
                }
            }

            //    while (1);
            double current_time = GetTime();
            time_save.push_back(current_time - pre_time);
            if ((i + 1) % 5 == 0) {
                times_period[((i + 1) / 5) - 1].push_back(current_time -
                                                          start_time);
                //      printf("%.6lf\n", current_time - start_time);
            }
            //    printf("%.6lf\n", current_time - pre_time);
            pre_time = current_time;
        }
        //    cout << "_cnt: " << _cnt << endl;
        //   get_min_double(time_save);

        for (int i = 0; i < anc.size(); i++) {
            int q = anc[i];
            for (int j = OSS[q].trash_can.size() - 1; j >= 0; j--) {
                int pos = OSS[q].trash_can[j];
                OSS[q].a[OSS[q].a[pos].next].previous = pos;
                OSS[q].a[OSS[q].a[pos].previous].next = pos;
            }
            OSS[q].trash_can.clear();
        }

        return result;
    }
    void object_setting(int n) {
        is_current_object = (int *)malloc(sizeof(int) * (n + 1));
        current_distance = (int *)malloc(sizeof(int) * (n + 1));
        group_height = (int *)malloc(sizeof(int) * (n + 1));
        current_state = (int *)malloc(sizeof(int) * (n + 1));
        current_stamp = 0;
        for (int i = 0; i <= n; i++)
            current_state[i] = 0;
    }
    void print(int root) {
        int cnt = 0;
        for (int i = OSS[root].a[0].next; i != 0; i = OSS[root].a[i].next) {
            printf("(%d, %d)", OSS[root].a[i].key, OSS[root].a[i].dist);
            cnt++;
        }
        printf("\n");

        if (cnt != OSS[root].size_num) {
            cout << "cnt: " << cnt << "     "
                 << "OSS[root].size_num: " << OSS[root].size_num << endl;
            while (1)
                ;
        }
    }
    bool double_objects(int p) {
        for (int i = OSS[p].a[0].next; i != 0; i = OSS[p].a[i].next) {
            if (OSS[p].a[i].key == OSS[p].a[OSS[p].a[i].previous].key)
                return false;
            if ((OSS[p].a[i].next != 0) &&
                (OSS[p].a[OSS[p].a[i].next].dist < OSS[p].a[i].dist))
                return false;
        }
        return true;
    }
    bool check_everyone() {
        for (int i = 0; i < n; i++)
            if (double_objects(i) == false)
                return false;
        return true;
    }
};

double get_mean_double(vector<double> &times) {
    double mean = 0.0;
    for (double val : times) {
        mean += val;
    }
    return mean / times.size();
}
double get_var_double(vector<double> &times) {
    double mean = get_mean_double(times);
    double var = 0.0;
    for (double val : times) {
        var += (val - mean) * (val - mean);
    }
    return var / times.size();
}
double get_max_double(vector<double> &times) {
    double max = 0.0;
    for (double val : times) {
        if (max < val)
            max = val;
    }
    return max;
}
} // namespace QUERY