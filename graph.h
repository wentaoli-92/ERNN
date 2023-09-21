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

// #define DEBUG

using namespace std;

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

  public:
    Graph(string path);
    ~Graph();

  public:
    void dijkstra();
    int updateWeight(int u, int v, int delta);
    void dynamicDijkstra(int u, int v, int delta);
    void recomputeDijkstra(int u, int v, int delta);

  public:
    int exact(int targetID, int budget);
    int greedy(int targetID, int budget, int pruneOP);

  public:
    vector<int> boundList;
    void commputeBound(int targetID);

  public:
    int maxWeight(int targetID, int budget);
    int chooseNbr(int targetID, int budget);

  public:
    void setPOI(vector<int> poi);
};

// Implementation
Graph::Graph(string filename) {
    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "Could not open file " << filename << endl;
        return;
    }

    int n, m;
    infile >> n >> m;

    nVertices = n;
    graph.resize(n);
    int u, v, w;
    set<pair<int, int>> uvPairs;

    for (int i = 0; i < m; i++) {
        infile >> u >> v >> w; // read edge
        u--, v--;              // starting from one
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

    for (int i = 0; i < n; i++) {
        graph[i].id = i;
        graph[i].isPOI = 0;
        graph[i].dist = INF;
        graph[i].closest_node = -1;
    }
    // cout << n << " " << m << endl;
}

Graph::~Graph() {}

void Graph::dijkstra() {
    for (int i = 0; i < nVertices; i++) {
        graph[i].dist = INF;
        graph[i].closest_node = -1;
    }

    vector<int> P;
    for (int i = 0; i < nVertices; i++) {
        if (graph[i].isPOI == 1)
            P.push_back(i);
    }

    priority_queue<Vertex> Q;
    vector<bool> visited(nVertices, false);

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
    vector<bool> visited(nVertices, false);

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

int Graph::exact(int targetID, int budget) {
    int nBefore = 0;
    int gain = -INF;
    int cost = 0;
    vector<tuple<int, int, int>> el;
    for (int i = 0; i < nVertices; ++i) {
        for (auto j : graph[i].adj) {
            if (i < j.first) {
                el.push_back(make_tuple(i, j.first, j.second));
            }
        }
        if (graph[i].closest_node == targetID)
            nBefore++;
    }
#ifdef DEBUG
    cout << "Before = " << nBefore << endl;
#endif

    vector<pair<int, int>> tempAns;
    int m = el.size();
    vector<int> a(m);
    iota(a.begin(), a.end(), 0);
    queue<pair<vector<int>, int>> subsetsQueue;
    subsetsQueue.push({{}, 0});

    while (!subsetsQueue.empty()) {
        auto currentSubset = subsetsQueue.front().first;
        int startIndex = subsetsQueue.front().second;
        subsetsQueue.pop();

        vector<Vertex> gtemp = graph;
#ifdef DEBUG
        cout << "Choose ";
#endif
        if (currentSubset.size() <= budget) {
            for (int element : currentSubset) {
                int u = get<0>(el[element]);
                int v = get<1>(el[element]);
                int w = get<2>(el[element]);
                dynamicDijkstra(u, v, w);
#ifdef DEBUG
                cout << "(" << u << ", " << v << "), ";
#endif
            }
        }

        int nafter = 0;
        for (int i = 0; i < nVertices; ++i) {
            if (graph[i].closest_node == targetID)
                nafter++;
        }
#ifdef DEBUG
        cout << "After = " << nafter << ", gain = " << nafter - nBefore << endl;
#endif
        int gainAfter = nafter - nBefore;
        if (gain < gainAfter) {
            gain = gainAfter;
            tempAns.clear();
            for (int element : currentSubset) {
                int u = get<0>(el[element]);
                int v = get<1>(el[element]);
                tempAns.push_back(make_pair(u, v));
            }
        }
        graph = gtemp;
        cost++;

        if (currentSubset.size() == budget) {
            continue;
        }

        for (int i = startIndex; i < m; i++) {
            vector<int> newSubset = currentSubset;
            newSubset.push_back(a[i]);
            subsetsQueue.push({newSubset, i + 1});
        }
    }
#ifdef DEBUG
    cout << "maximum gain is " << gain << endl;
    cout << "select: ";
    for (auto ans : tempAns) {
        cout << "(" << ans.first << ", " << ans.second << "), ";
    }
    cout << endl;
    cout << "cost = " << cost << endl;
#endif
    return gain;
}

void Graph::commputeBound(int targetID) {
    boundList.clear();
    vector<int> a;
    a.clear();
    for (int i = 0; i < nVertices; ++i) {
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

    // Step 4: Generate map e (fill in gaps between distances)
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
    // Step 5: Generate vector boundList (only preserve the values in each
    // position)
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

    for (int i = 0; i < nVertices; ++i) {
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
        vector<int> dist(nVertices, numeric_limits<int>::max());
        vector<bool> visited(nVertices, false);
        vector<int> max_edge_on_path(nVertices, 0);

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

                auto [it, inserted] = uvPairs.insert(make_pair(u, v));
                auto [it2, inserted2] = uvPairs.insert(make_pair(v, u));
                if (inserted == false or inserted2 == false)
                    continue;

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
                for (int j = 0; j < nVertices; ++j) {
                    if (graph[j].closest_node == targetID)
                        cnt++;
                }
                int gainAfter = cnt - nBefore; // gap
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

    //cout << "cost = " << cost << endl;

    if(totalGain <= 0) totalGain = 0;
    return totalGain;
}

struct EdgeWeightComparator {
    bool operator()(const tuple<int, int, int> &e1,
                    const tuple<int, int, int> &e2) const {
        return get<2>(e1) > get<2>(e2);
    }
};

int Graph::maxWeight(int targetID, int budget) {
    int nBefore = 0;
    vector<tuple<int, int, int>> el;
    for (int i = 0; i < nVertices; ++i) {
        for (auto j : graph[i].adj) {
            if (i < j.first) {
                el.push_back(make_tuple(i, j.first, j.second));
            }
        }
        if (graph[i].closest_node == targetID)
            nBefore++;
    }
#ifdef DEBUG
    cout << "Before = " << nBefore << endl;
#endif

    priority_queue<tuple<int, int, int>, vector<tuple<int, int, int>>,
                   EdgeWeightComparator>
        min_heap;

    for (const auto &edge : el) {
        if (min_heap.size() < budget) {
            min_heap.push(edge);
        } else if (get<2>(edge) > get<2>(min_heap.top())) {
            min_heap.pop();
            min_heap.push(edge);
        }
    }

    vector<tuple<int, int, int>> largest_edges;
    while (!min_heap.empty()) {
        tuple<int, int, int> largest_edge = min_heap.top();
        min_heap.pop();

        int u = get<0>(largest_edge);
        int v = get<1>(largest_edge);
        int w = get<2>(largest_edge);
        dynamicDijkstra(u, v, w);
#ifdef DEBUG
        cout << "(" << u << ", " << v << "), ";
#endif
    }
#ifdef DEBUG
    cout << endl;
#endif

    int nafter = 0;
    for (int i = 0; i < nVertices; ++i) {
        if (graph[i].closest_node == targetID)
            nafter++;
    }
#ifdef DEBUG
    cout << "After = " << nafter << ", gain = " << nafter - nBefore << endl;
#endif
    int gain = nafter - nBefore;
#ifdef DEBUG
    cout << "maximum gain is " << gain << endl;
#endif

    return gain;
}

int Graph::chooseNbr(int targetID, int budget) {
    int nBefore = 0;
    for (int i = 0; i < nVertices; ++i) {
        if (graph[i].closest_node == targetID)
            nBefore++;
    }
#ifdef DEBUG
    cout << "Before = " << nBefore << endl;
#endif

    vector<bool> visited(nVertices, false);
    queue<int> q;
    set<tuple<int, int, int>> el;
    visited[targetID] = true;
    q.push(targetID);

    while (!q.empty() && el.size() < budget) {
        int u = q.front();
        q.pop();

        for (const auto &edge : graph[u].adj) {
            int v = edge.first;
            int w = edge.second;
            if (!visited[v]) {
                visited[v] = true;
                q.push(v);
            }

            tuple<int, int, int> e1 = make_tuple(u, v, w);
            tuple<int, int, int> e2 = make_tuple(v, u, w);
            if (el.find(e1) == el.end() && el.find(e2) == el.end()) {
                el.insert(e1);
                if (el.size() >= budget) {
                    break;
                }
            }
        }
    }

    for (const auto &edge : el) {
        int u = get<0>(edge);
        int v = get<1>(edge);
        int w = get<2>(edge);
        dynamicDijkstra(u, v, w);
#ifdef DEBUG
        cout << "(" << u << ", " << v << "), ";
#endif
    }
#ifdef DEBUG
    cout << endl;
#endif

    int nafter = 0;
    for (int i = 0; i < nVertices; ++i) {
        if (graph[i].closest_node == targetID)
            nafter++;
    }
#ifdef DEBUG
    cout << "After = " << nafter << ", gain = " << nafter - nBefore << endl;
#endif
    int gain = nafter - nBefore;
#ifdef DEBUG
    cout << "maximum gain is " << gain << endl;
#endif
    return gain;
}
