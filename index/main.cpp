#include "graph.h"


int returnNumOfVertices(string filename) {
    ERNN::Graph g(filename.c_str());
    return g.nVertices;
}

void testDynamic(string filename, vector<int> poi) {
    ERNN::Graph g(filename.c_str());
    g.setPOI(poi);
    g.dijkstra();
}

void DEBISTAR(string filename, vector<int> poi, int targetID, int rnd) {
    ERNN::Graph g(filename.c_str());
    g.setPOI(poi);
    g.dijkstra();
    g.f_n = rnd;

    auto start = chrono::high_resolution_clock::now();
    int gain = g.greedy(targetID, 1,
                        ERNN::Prune::Prune1 | ERNN::Prune::Prune2 | ERNN::Prune::Prune3);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "DEBISTAR: Gain = " << gain
         << ", Time = " << duration.count() * 1.0 / 1000 << endl;
}

void UsingIndex(string filename, vector<int> poi, int targetID, int rnd) {
    ERNN::Graph g(filename.c_str());
    g.setPOI(poi);
    g.f_n = rnd;
    g.dijkstra_index();
    auto start = chrono::high_resolution_clock::now();
    int gain = g.greedyIndex(targetID, 1,
                        ERNN::Prune::Prune1 | ERNN::Prune::Prune2 | ERNN::Prune::Prune3);
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::milliseconds>(end - start);
    cout << "Index: Gain = " << gain
         << ", Time = " << duration.count() * 1.0 / 1000 << endl;
}

int main(int argc, char *argv[]) {
    string strArg = argv[1];
    string filename = "../data/" + strArg + ".tmp";
    int n = returnNumOfVertices(filename);

    int numOfPOI = std::atoi(argv[2]);
    int numberOfTest = std::atoi(argv[3]);
    int budget = 1;
    int rnd = std::atoi(argv[4]);

    vector<int> poi;
    poi.clear();
    srand(1234);
    for (int i = 0; i < numOfPOI; i++) {
        int randomNum = rand() % n;
        poi.push_back(randomNum);
    }
    vector<int> vid(numOfPOI);
    iota(vid.begin(), vid.end(), 0);

    default_random_engine engine(rnd);
    shuffle(vid.begin(), vid.end(), engine);
    for (int i = 0; i < numberOfTest; ++i) {
        int id = vid[i];
        int targetID = poi[id];
        cout << "Round : " << i << ", targetID = " << targetID << endl;
        //cout << targetID << endl;
        DEBISTAR(filename, poi, targetID, rnd);
        UsingIndex(filename, poi, targetID, rnd);
    }
    return 0;
}
