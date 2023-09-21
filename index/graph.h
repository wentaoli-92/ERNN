#include "index.h"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <queue>
#include <random>
#include <set>
#include <thread>
#include <vector>
#include <chrono>

//#define DEBUG

using namespace std;

namespace ERNN {
const int INF = 1e9;
typedef pair<int, int> pii;

enum Prune { NoPrune = 0, Prune1 = 1, Prune2 = 2, Prune3 = 4 };

class Vertex {
  public:
    vector<pair<int, int>> adj;
    int id;
    int isPOI;

    int dist;
    int closest_node;

    bool operator<(const Vertex &rhs) const { return dist > rhs.dist; }
};

class Graph {
  public:
    vector<Vertex> graph;
    int nVertices;
    string path_n;
    int b_n;
    int f_n; 
  public:
    Graph(string path);
    ~Graph();

  public:
    void dijkstra();
    void dijkstra_index();
    int updateWeight(int u, int v, int delta);
    void dynamicDijkstra(int u, int v, int delta);
    void recomputeDijkstra(int u, int v, int delta);

  public:
    int greedy(int targetID, int budget, int pruneOP);
    int greedyIndex(int targetID, int budget, int pruneOP);

  public:
    vector<int> boundList;
    void commputeBound(int targetID);

  public:
    void setPOI(vector<int> poi);
};

// Implementation
Graph::Graph(string filename) {
    path_n = filename;	
    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "Could not open file " << filename << endl;
        return;
    }

    int n, m;
    infile >> n >> m;

    nVertices = n;
    graph.resize(n + 1);
    int u, v, w;
    set<pair<int, int>> uvPairs;

    for (int i = 0; i < m; i++) {
        infile >> u >> v >> w; // read edge
        auto [it, inserted] = uvPairs.insert(make_pair(u, v));
        if (inserted) {
            graph[u].adj.push_back(make_pair(v, w));
        }
        auto [it2, inserted2] = uvPairs.insert(make_pair(v, u));
        if (inserted2) {
            graph[v].adj.push_back(make_pair(u, w));
        }
    }
    infile.close();

    for (int i = 0; i <= n; i++) {
        graph[i].id = i;
        graph[i].isPOI = 0;
        graph[i].dist = INF;
        graph[i].closest_node = -1;
    }
    // cout << n << " " << m << endl;
}

Graph::~Graph() {}

void Graph::dijkstra_index() {
    int n = nVertices;
    vector<tuple<int, int, int>> el;
    for (int u = 0; u <= n; ++u) {
        for (const auto &edge : graph[u].adj) {
            int v = edge.first;
            int w = edge.second;
            if (u < v)
                el.push_back(make_tuple(u, v, w));
        }
    }

    string fileout = "index/tmp.index";
    size_t pos = path_n.rfind("/");
    string name;
    if (pos != string::npos) {
        name = path_n.substr(pos + 1);
        pos = name.rfind(".");
        name = name.substr(0, pos);
    }

    fileout = fileout + name + to_string(f_n);
	
    INDEX::Tree_Decomposition td;
    td.fout = fopen(fileout.c_str(), "wb");
    td.G = INDEX::Graph(nVertices, el);
    td.H = td.G;
    td.reduce();
    td.makeTree();
    td.makeIndex();
    td.print();
    td.makeIndex2();
    td.reducePos();
    td.cntSize();
    fclose(td.fout);

    // query
    // char *indexfile = "tmp.index";
    QUERY::readIndex(fileout.c_str());

    QUERY::higher = (int *)malloc(sizeof(int) * (n + 1));
    QUERY::cloest_higher = (int *)malloc(sizeof(int) * (n + 1));
    for (int i = 0; i < n; i++) {
        QUERY::higher[i] = n;
        QUERY::cloest_higher[i] = INT_MAX;
        for (int j = 0; j < QUERY::posSize[i]; j++) {
            if (QUERY::height[QUERY::belong[QUERY::pos[i][j]]] <
                QUERY::higher[i])
                QUERY::higher[i] =
                    QUERY::height[QUERY::belong[QUERY::pos[i][j]]];
            int _d =
                QUERY::distanceQuery(QUERY::uniqueVertex[i], QUERY::pos[i][j]);
            if (QUERY::cloest_higher[i] > _d)
                QUERY::cloest_higher[i] = _d;
        }
    }

    //----------------------prepare-------------
    QUERY::LOG2 = (int *)malloc(sizeof(int) * (n * 2 + 10));
    QUERY::LOGD = (int *)malloc(sizeof(int) * (n * 2 + 10));
    int k = 0, j = 1;
    for (int i = 0; i < n * 2 + 10; i++) {
        if (i > j * 2) {
            j *= 2;
            k++;
        }
        QUERY::LOG2[i] = k;
        QUERY::LOGD[i] = j;
    }

    QUERY::insert_type = "neighbor";
    QUERY::query_type = "optimal";

    //-------------test------------------
    QUERY::kNN knn;
    knn.create_kNN_index();
    int number_object = 0;
    for (int i = 0; i <= n; ++i) {
        if (graph[i].isPOI == 1) {
            number_object++;
            // cout << i << endl;
        }
    }
    knn.object_setting(n);

    for (int i = 0; i <= n; ++i) {
        if (graph[i].isPOI == 1)
            QUERY::is_current_object[i] = 1;
    }


    int cnt_object = 0;
    for (int i = 1; i <= n; i++)
        if (QUERY::is_current_object[i] == 1)
            cnt_object++;
    knn.initialize_object();

    if (knn.check_everyone()) {
        //printf("right\n");
    }

    for (int i = 0; i < knn.period; i++)
        knn.times_period[i].clear();

    int q_n = 1;
    knn.query_mark = (int *)malloc(sizeof(int) * (n + 1));
    knn.query_mark_stamp = 0;
    for (int i = 0; i <= n; i++)
        knn.query_mark[i] = 0;
   

    QUERY::tmp_dis = (int *)malloc(sizeof(int) * (n + 1));
    vector<pair<int, int>> res;

    for (int i = 1; i <= n; ++i) {
        //cout << i << endl;
        int x = i;
        int knum = 1;
     
        res = knn.query(x, knum);
        //cout << res.size() << endl;
        //cout << res.size() << endl;
        if (res.size() == 0) {
            res = knn.query_naive(x, 100);
            //cout << res.size() << endl;
        }
        if (res.size() == 0){
            continue;
        }
            
        //cout << i << " " << res[0].first << " " << res[0].second << endl;
        graph[i].closest_node = res[0].first;
        graph[i].dist = res[0].second;
    }
    //cout << endl;
}

void Graph::dijkstra() {
    for (int i = 0; i <= nVertices; i++) {
        graph[i].dist = INF;
        graph[i].closest_node = -1;
    }

    vector<int> P;
    for (int i = 0; i <= nVertices; i++) {
        if (graph[i].isPOI == 1)
            P.push_back(i);
    }

    priority_queue<Vertex> Q;
    vector<bool> visited(nVertices + 1, false);

    for (int u : P) {
        graph[u].dist = 0;
        graph[u].closest_node = u;
        Q.push(graph[u]);
    }

    while (!Q.empty()) {
        Vertex current = Q.top();
        Q.pop();
        int v = current.id;
        int d = current.dist;
        int u = current.closest_node;

        if (visited[v])
            continue;
        visited[v] = true;

        for (const auto &edge : graph[v].adj) {
            int w = edge.first;
            int dist_v_w = edge.second;

            if (!visited[w] && d + dist_v_w < graph[w].dist) {
                graph[w].dist = d + dist_v_w;
                graph[w].closest_node = u;
                Q.push(graph[w]);
            }
        }
    }
}

void Graph::setPOI(vector<int> poi) {
    for (int i : poi)
        graph[i].isPOI = 1;
}

int Graph::updateWeight(int u, int v, int delta) {
    int ans = 0;
    for (auto &edge : graph[u].adj) {
        if (edge.first == v) {
            if (delta == -1)
                edge.second = 0;
            else
                edge.second = max(edge.second - delta, 0);
            ans = edge.second;
            break;
        }
    }
    for (auto &edge : graph[v].adj) {
        if (edge.first == u) {
            if (delta == -1)
                edge.second = 0;
            else
                edge.second = max(edge.second - delta, 0);
            break;
        }
    }
    return ans;
}

void Graph::dynamicDijkstra(int u, int v, int delta) {
    int w = updateWeight(u, v, delta);

    priority_queue<Vertex> Q;
    vector<bool> visited(nVertices + 1, false);

    if (graph[u].dist + w < graph[v].dist) {
        graph[v].dist = graph[u].dist + w;
        graph[v].closest_node = graph[u].closest_node;
        Q.push(graph[v]);
    }

    if (graph[v].dist + w < graph[u].dist) {
        graph[u].dist = graph[v].dist + w;
        graph[u].closest_node = graph[v].closest_node;
        Q.push(graph[u]);
    }

    while (!Q.empty()) {
        Vertex current = Q.top();
        Q.pop();
        int v = current.id;
        int d = current.dist;
        int u = current.closest_node;

        if (visited[v])
            continue;
        visited[v] = true;

        for (const auto &edge : graph[v].adj) {
            int w = edge.first;
            int dist_v_w = edge.second;

            if (!visited[w] && d + dist_v_w < graph[w].dist) {
                graph[w].dist = d + dist_v_w;
                graph[w].closest_node = u;
                Q.push(graph[w]);
            }
        }
    }
}

void Graph::recomputeDijkstra(int u, int v, int delta) {
    int w = updateWeight(u, v, delta);
    dijkstra();
}

void Graph::commputeBound(int targetID) {
    boundList.clear();
    vector<int> a;
    a.clear();
    for (int i = 0; i <= nVertices; ++i) {
	if(graph[i].dist == INF) continue;    
        if (graph[i].closest_node != targetID and graph[i].dist != 0)
            a.push_back(graph[i].dist);
    }
    // for(auto e : a) cout << e << endl;

    vector<int> b;
    map<int, int> c;
    map<int, int> d;
    map<int, int> e;

    // Step 1: Generate vector b
    for (int i : a) {
        b.push_back(i - 1);
    }

    // Step 2: Generate map c
    for (int i : b) {
        c[i]++;
    }

    // Step 3: Generate map d (accumulate values of c in reverse order)
    int accumulated = 0;
    for (auto it = c.rbegin(); it != c.rend(); ++it) {
        accumulated += it->second;
        d[it->first] = accumulated;
    }

    // Step 4: Generate map e (fill in gaps)
    int lastKey = d.rbegin()->first;
    int lastValue = d[lastKey];
    for (int i = lastKey; i >= 0; --i) {
        if (d.find(i) != d.end()) {
            e[i] = d[i];
            lastValue = d[i];
        } else {
            e[i] = lastValue;
        }
    }

    boundList.clear();
    for (int i = 0; i < e.size(); ++i) {
        boundList.push_back(e[i]);
    }
}

int Graph::greedy(int targetID, int budget, int pruneOP) {
    int nBefore = 0;
    int nAfter = 0;
    int gain = -INF;   // each round
    int totalGain = 0; // whole round
    int cost = 0;

    graph[0].closest_node = 0;
    graph[0].dist = 0;

    for (int i = 0; i <= nVertices; ++i) {
        if (graph[i].closest_node == targetID)
            nBefore++;
    }

#ifdef DEBUG
    cout << "Before = " << nBefore << endl;
#endif
    commputeBound(targetID);

    int roundID = 0;
    while (budget > 0) {
        budget--;
        tuple<int, int, int> tempAnsEpoh(0, 0, 0);
        vector<Vertex> gtemp = graph;
        vector<int> dist(nVertices + 1, numeric_limits<int>::max());
        vector<bool> visited(nVertices + 1, false);
        vector<int> max_edge_on_path(nVertices + 1, 0);

        priority_queue<pii, vector<pii>, greater<pii>> pq;
        pq.push({0, targetID});
        dist[targetID] = 0;
	set<pair<int, int>> uvPairs;
        uvPairs.clear();

        while (!pq.empty()) {
            int u = pq.top().second;
            pq.pop();

            if (visited[u])
                continue;
            visited[u] = true;

            if (pruneOP & Prune::Prune3) {
                if (boundList[dist[u]] <= gain)
                    break;
            }

            for (const auto &edge : graph[u].adj) {
                int v = edge.first;
                int w = edge.second;

                if (pruneOP & Prune::Prune1) {
                    // prune if both u and v not contain a targetID
                    if (graph[u].closest_node != targetID and
                        graph[v].closest_node != targetID)
                        continue;
                    // prune if v not contains a targetID and is far
                    if (graph[v].closest_node != targetID and
                        graph[v].dist < graph[u].dist)
                        continue;
                    if (graph[u].closest_node != targetID and
                        graph[u].dist < graph[v].dist)
                        continue;
                }

                
                auto [it, inserted] = uvPairs.insert(make_pair(u, v));
                auto [it2, inserted2] = uvPairs.insert(make_pair(v, u));
                if (inserted == false or inserted2 == false)
                    continue;

                if (pruneOP & Prune::Prune2) {
                    // prune if the weight is smaller
                    if (w <= max_edge_on_path[u])
                        continue;
                }

                cost++;
                dynamicDijkstra(u, v, w);
#ifdef DEBUG
                cout << "(" << u << ", " << v << "), ";
#endif
                int cnt = 0;
                for (int j = 0; j <= nVertices; ++j) {
                    if (graph[j].closest_node == targetID)
                        cnt++;
                }
                int gainAfter = cnt - nBefore; // gap
                // cout << u<< " " << v << " " << gainAfter << endl;
#ifdef DEBUG
                cout << "After = " << cnt << ", Before = " << nBefore
                     << ", gain = " << gainAfter << endl;
#endif
                if (gainAfter > gain) {
                    gain = gainAfter;
                    tempAnsEpoh = make_tuple(u, v, w);
                    nAfter = cnt;
                }
                graph = gtemp;
            }

            for (const auto &edge : graph[u].adj) {
                int v = edge.first;
                int w = edge.second;

                if (dist[u] + w < dist[v]) {
                    dist[v] = dist[u] + w;
                    max_edge_on_path[v] = max(max_edge_on_path[u], w);
                    pq.push({dist[v], v});
                } else if (dist[u] + w == dist[v]) {
                    max_edge_on_path[v] =
                        max(max_edge_on_path[v], max(max_edge_on_path[u], w));
                }
            }
        }

        int u = get<0>(tempAnsEpoh);
        int v = get<1>(tempAnsEpoh);
        int w = get<2>(tempAnsEpoh);
        dynamicDijkstra(u, v, w);
        nBefore = nAfter;
        totalGain += gain;
#ifdef DEBUG
        cout << "round " << ++roundID << " choose (" << u << " , " << v << ")"
             << endl;
#endif
        gain = -INF;
        commputeBound(targetID);
    }
#ifdef DEBUG
    cout << "The final gain = " << totalGain << endl;
    cout << "cost = " << cost << endl;
#endif

    cout << "cost = " << cost << endl;

    return totalGain;
}

int Graph::greedyIndex(int targetID, int budget, int pruneOP) {
    int nBefore = 0;
    int nAfter = 0;
    int gain = -INF;   // each round
    int totalGain = 0; // whole round
    int cost = 0;
    b_n = budget;

    graph[0].closest_node = 0;
    graph[0].dist = 0;

    for (int i = 0; i <= nVertices; ++i) {
        if (graph[i].closest_node == targetID)
            nBefore++;
    }

#ifdef DEBUG
    cout << "Before = " << nBefore << endl;
#endif
    commputeBound(targetID);

    int roundID = 0;
    while (budget > 0) {
        budget--;
        tuple<int, int, int> tempAnsEpoh(0, 0, 0);
        vector<Vertex> gtemp = graph;
        vector<int> dist(nVertices + 1, numeric_limits<int>::max());
        vector<bool> visited(nVertices + 1, false);
        vector<int> max_edge_on_path(nVertices + 1, 0);

        priority_queue<pii, vector<pii>, greater<pii>> pq;
        pq.push({0, targetID});
        dist[targetID] = 0;
	set<pair<int, int>> uvPairs;
        uvPairs.clear();

        while (!pq.empty()) {
            int u = pq.top().second;
            pq.pop();

            if (visited[u])
                continue;
            visited[u] = true;

            if (pruneOP & Prune::Prune3) {
                if (boundList[dist[u]] <= gain)
                    break;
            }

            for (const auto &edge : graph[u].adj) {
                int v = edge.first;
                int w = edge.second;

                //cout << u << " " << v << " " << w << endl;

                //cout << "here 1" << endl;
                if (pruneOP & Prune::Prune1) {
                    // prune if both u and v not contain a targetID
                    if (graph[u].closest_node != targetID and
                        graph[v].closest_node != targetID)
                        continue;
                    // prune if v not contains a targetID and is far
                    if (graph[v].closest_node != targetID and
                        graph[v].dist < graph[u].dist)
                        continue;
                    if (graph[u].closest_node != targetID and
                        graph[u].dist < graph[v].dist)
                        continue;
                }

                auto [it, inserted] = uvPairs.insert(make_pair(u, v));
                auto [it2, inserted2] = uvPairs.insert(make_pair(v, u));
                if (inserted == false or inserted2 == false)
                    continue;

                // cout << "here 4" << endl;
                if (pruneOP & Prune::Prune2) {
                    // prune if the weight is smaller
                    if (w <= max_edge_on_path[u])
                        continue;
                }

                //cout << "do1" << endl;
                cost++;
                updateWeight(u, v, w);
                //cout << "do2" << endl;
                dijkstra_index();
                //cout << "here" << endl;
                graph[0].closest_node = 0;
                graph[0].dist = 0;
                //cout << "do3" << endl;
                //  dynamicDijkstra(u, v, w);
#ifdef DEBUG
                cout << "(" << u << ", " << v << "), ";
#endif
                int cnt = 0;
                for (int j = 0; j <= nVertices; ++j) {
                    if (graph[j].closest_node == targetID)
                        cnt++;
                }
                int gainAfter = cnt - nBefore; // gap
                // cout << u<< " " << v << " " << gainAfter << endl;
#ifdef DEBUG
                cout << "After = " << cnt << ", Before = " << nBefore
                     << ", gain = " << gainAfter << endl;
#endif
                if (gainAfter > gain) {
                    gain = gainAfter;
                    tempAnsEpoh = make_tuple(u, v, w);
                    nAfter = cnt;
                }
                graph = gtemp;
            }

            for (const auto &edge : graph[u].adj) {
                int v = edge.first;
                int w = edge.second;

                if (dist[u] + w < dist[v]) {
                    dist[v] = dist[u] + w;
                    max_edge_on_path[v] = max(max_edge_on_path[u], w);
                    pq.push({dist[v], v});
                } else if (dist[u] + w == dist[v]) {
                    max_edge_on_path[v] =
                        max(max_edge_on_path[v], max(max_edge_on_path[u], w));
                }
            }
        }

        int u = get<0>(tempAnsEpoh);
        int v = get<1>(tempAnsEpoh);
        int w = get<2>(tempAnsEpoh);
        updateWeight(u, v, w);
        dijkstra_index();
        graph[0].closest_node = 0;
        graph[0].dist = 0;
        // dynamicDijkstra(u, v, w);
        nBefore = nAfter;
        totalGain += gain;
#ifdef DEBUG
        cout << "round " << ++roundID << " choose (" << u << " , " << v << ")"
             << endl;
#endif
        gain = -INF;
        commputeBound(targetID);
    }
#ifdef DEBUG
    cout << "The final gain = " << totalGain << endl;
    cout << "cost = " << cost << endl;
#endif

    cout << "cost = " << cost << endl;

    return totalGain;
}
} // namespace ERNN
