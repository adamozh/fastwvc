#include <algorithm>
#include <chrono>
#include <iostream>
#include <random>
#include <unordered_set>
#include <utility>
#include <vector>

using namespace std;
using namespace chrono;

typedef pair<int, int> pii;
typedef vector<bool> vb;
typedef vector<int> vi;

struct edge {
    int u;
    int v;
};

bool debug = false;
int INF = 1e9;

int n, e;
vi node_weights;
vi node_degrees;
vector<edge> edge_list;
vector<vector<pii>> AL; // pair is (neighbour, edge_index)

int best_weight = INF;
vb best_cover;

int step_cover_weight;
vb step_cover;

int step = 0;

vb node_conf_change;
vi node_dscores;
vi edge_weights;
vb tabu_list;
vi remove_candidates;
vi timestamp;
unordered_set<int> uncovered_edges;

// these variables are used for an extension which is not implemented for now
int avg_weight = 1;
int delta_total_weight = 0;
double p_scale = 0.3;
int threshold;

/* DEBUGGING METHODS */
void print_vi(vi &v, string name) {
    cout << name << ": ";
    for (auto &i : v) {
        cout << i << " ";
    }
    cout << endl;
}

void print_vb(vb &v, string name) {
    cout << name << ": ";
    for (const auto &b : v) {
        cout << b; // Automatically converts to 1 or 0
    }
    cout << endl;
}

void print_set(unordered_set<int> &s, string name) {
    cout << name << ": ";
    for (const auto &i : s) {
        cout << i << ' ';
    }
    cout << endl;
}

// generates a good initial candidate using ConstructMWVC
void init_mwvc() {
    int best_weight_init = INF;
    vb best_cover_init(n, false);

    // 1. generate a first "draft" greedily by iterating through edges and
    // picking the node with better weight to degree ratio
    int curr_weight = 0;
    vb curr_cover(n, false);
    for (edge e : edge_list) {
        if (!curr_cover[e.u] && !curr_cover[e.v]) {
            double degreeU = (double)node_degrees[e.u] / (double)node_weights[e.u];
            double degreeV = (double)node_degrees[e.v] / (double)node_weights[e.v];
            // choose the endpoint with larger degree/weight to be added
            int better_node = degreeU > degreeV ? e.u : e.v;
            curr_weight += node_weights[better_node];
            curr_cover[better_node] = true;
        }
    }
    best_weight_init = curr_weight;
    best_cover_init = curr_cover;

    // 2. repeated extending phase. first split the edges into blocks
    vi blocks(e / 1024 + 1);
    for (int i = 0; i < e / 1024 + 1; i++) {
        blocks[i] = i;
    }
    int tries = 50;
    while (--tries >= 0) {
        curr_weight = 0;
        curr_cover = vb(n, false);

        // shuffle order of blocks
        shuffle(blocks.begin(), blocks.end(), default_random_engine(0));

        for (int block_id : blocks) {
            int begin = block_id * 1024;
            int end = block_id == e / 1024 ? e - 1 : begin + 1024;
            int block_size = end - begin + 1;
            vi idx(block_size); // holds edge indices within the current block
            for (int i = 0; i < block_size; i++) {
                idx[i] = begin + i;
            }
            while (block_size > 0) {
                int i = rand() % block_size;
                edge edge = edge_list[idx[i]];
                int u = edge.u, v = edge.v;
                block_size--;
                swap(idx[i], idx[block_size]);
                // if both endpoints are not in cover, add the better one
                if (!curr_cover[u] && !curr_cover[v]) {
                    double degreeU = (double)node_degrees[u] / (double)node_weights[u];
                    double degreeV = (double)node_degrees[v] / (double)node_weights[v];
                    int better_node = degreeU > degreeV ? u : v;
                    curr_weight += node_weights[better_node];
                    curr_cover[better_node] = true;
                }
            }
        }

        if (curr_weight < best_weight_init) {
            best_weight_init = curr_weight;
            best_cover_init = curr_cover;
        }
    }

    // 3. shrinking phase
    // harm_values[i] = if node i is removed from the cover, how many edges
    // would become uncovered aka number of edges that are solely covered by it
    vi harm_values(n, 0);
    for (int i = 0; i < e; i++) {
        edge edge = edge_list[i];
        if (best_cover_init[edge.u] && !best_cover_init[edge.v]) {
            harm_values[edge.u]++;
        }
        if (best_cover_init[edge.v] && !best_cover_init[edge.u]) {
            harm_values[edge.v]++;
        }
    }

    // remove nodes from cover if harm value is 0
    for (int i = 0; i < n; i++) {
        if (best_cover_init[i] && harm_values[i] == 0) {
            best_cover_init[i] = false;
            best_weight_init -= node_weights[i];
            for (pii p : AL[i]) {
                int neighbour = p.first;
                harm_values[neighbour]++; // neighbour is now responsible for
                                          // covering one more edge
            }
        }
    }

    // set this result as global optimal (for now)
    best_weight = best_weight_init;
    best_cover = best_cover_init;
}

// calculate the initial dscore values based on init candidate, for calculating
// gain / loss initially, all are loss values since all edges are initially
// covered, hence dscores are neg
void generate_dscores_init() {
    for (int i = 0; i < e; i++) {
        edge edge = edge_list[i];
        // edges that are covered by only one node constribute to the node's
        // dscore if the edge becomes uncovered, cost will increase by the edge
        // weight
        if (best_cover[edge.u] && !best_cover[edge.v]) {
            node_dscores[edge.u] -= edge_weights[i];
        } else if (!best_cover[edge.u] && best_cover[edge.v]) {
            node_dscores[edge.v] -= edge_weights[i];
        }
    }
}

// generate an initial list of nodes that are part of the cover, to be
// maintained throughout the algorithm for ease of iterating through nodes that
// are part of the cover
void generate_remove_candidates() {
    for (int i = 0; i < n; i++) {
        if (best_cover[i]) {
            remove_candidates.push_back(i);
        }
    }
}

// in reference implementation, removal is O(1) by tracking indices, but O(n) is chosen for
// simplicity and because it is theoretically not crucial, as complexity of other steps outweighs
// this
void remove_node(int node, int &step_cover_weight, vb &step_cover) {
    if (debug) cout << "remove node " << node << endl;

    // remove from the current step cover
    step_cover_weight -= node_weights[node];
    step_cover[node] = false;

    remove_candidates.erase(remove(remove_candidates.begin(), remove_candidates.end(), node),
                            remove_candidates.end());
    node_dscores[node] = -node_dscores[node]; // now adding it will cover X amount of edges
    node_conf_change[node] = false;

    // update node_dscores and node_conf_change of neighbours
    for (pii p : AL[node]) {
        int neighbour = p.first, edge_idx = p.second;
        if (debug)
            cout << "initial dscore of " << neighbour << " = " << node_dscores[neighbour] << endl;
        node_conf_change[neighbour] = true;

        if (step_cover[neighbour]) {
            node_dscores[neighbour] -= edge_weights[edge_idx]; // neighbour has more responsibility
        } else {
            node_dscores[neighbour] += edge_weights[edge_idx];
            node_conf_change[neighbour] = true; // neighbour eligible to be added into cover
            uncovered_edges.insert(edge_idx);
        }
        if (debug)
            cout << "new dscore of " << neighbour << " = " << node_dscores[neighbour] << endl;
    }
}

// basically the opposite of remove_node
void add_node(int node, int &step_cover_weight, vb &step_cover) {
    step_cover_weight += node_weights[node];
    step_cover[node] = true;
    remove_candidates.push_back(node);
    node_dscores[node] = -node_dscores[node];
    // don't need to set it's own conf_change, only need to set when removing

    for (pii p : AL[node]) {
        int neighbour = p.first, edge_idx = p.second;
        if (step_cover[neighbour]) {
            node_dscores[neighbour] += edge_weights[edge_idx];
        } else {
            node_dscores[neighbour] -= edge_weights[edge_idx];
            node_conf_change[neighbour] = true;
            uncovered_edges.erase(edge_idx);
        }
    }
}

// main step funciton for local search. at every step, remove 2 nodes and try to reconstruct
void step_function() {
    // remove the first node by greedily choosing the one with smallest loss
    if (debug) print_vi(remove_candidates, "remove_candidates");
    double smallest_loss = INF;
    int node_to_remove_1 = remove_candidates[0];
    for (int i = 0; i < remove_candidates.size(); i++) {
        int node = remove_candidates[i];
        double loss = (double)abs(node_dscores[node]) / (double)node_weights[node];
        if (loss < smallest_loss) {
            smallest_loss = loss;
            node_to_remove_1 = node;
        }
    }
    remove_node(node_to_remove_1, step_cover_weight, step_cover);

    // choose second node to remove using BMS strategy (biased memory saving)
    int tries = 100;
    // use same variable node_to_remove from before
    int idx = rand() % remove_candidates.size();
    int node_to_remove_2 = remove_candidates[idx];
    double current_loss =
        (double)abs(node_dscores[node_to_remove_2]) / (double)node_weights[node_to_remove_2];
    while (--tries >= 0) {
        int node = remove_candidates[rand() % remove_candidates.size()];
        double loss = (double)abs(node_dscores[node]) / (double)node_weights[node];
        if (tabu_list[node]) continue;
        if (loss > current_loss) continue;
        if (loss < current_loss) {
            node_to_remove_2 = node;
            current_loss = loss;
        } else if (timestamp[node] < timestamp[node_to_remove_2]) { // tiebreak by timestamp
            node_to_remove_2 = node;
            current_loss = loss;
        }
    }
    remove_node(node_to_remove_2, step_cover_weight, step_cover);

    if (debug) print_set(uncovered_edges, "uncovered_edges");

    tabu_list = vb(n, false);

    while (!uncovered_edges.empty()) {
        int node_to_add = 0;
        double best_gain = 0;
        // choose a vertex to add, considering both recently removed nodes
        for (int consider : {node_to_remove_1, node_to_remove_2}) {
            // consider its neighbours
            for (pii p : AL[consider]) {
                int neighbour = p.first;
                if (step_cover[neighbour] || node_conf_change[neighbour] == false) continue;
                double gain = (double)node_dscores[neighbour] / (double)node_weights[neighbour];
                if (gain > best_gain) {
                    best_gain = gain;
                    node_to_add = neighbour;
                } else if (gain == best_gain && timestamp[neighbour] < timestamp[node_to_add]) {
                    node_to_add = neighbour;
                }
            }
            // consider the node itself
            if (node_conf_change[consider] && !step_cover[consider]) {
                int gain = (double)node_dscores[consider] / (double)node_weights[consider];
                if (gain > best_gain) {
                    best_gain = gain;
                    node_to_add = consider;
                } else if (gain == best_gain && timestamp[consider] < timestamp[node_to_add]) {
                    node_to_add = consider;
                }
            }
        }

        add_node(node_to_add, step_cover_weight, step_cover);
        if (debug) cout << "added node " << node_to_add << endl;
        if (debug) print_set(uncovered_edges, "uncovered_edges");
        tabu_list[node_to_add] = true;
        timestamp[node_to_add] = step;

        // update edge weights of uncovered edges
        for (auto &edge_idx : uncovered_edges) {
            edge edge = edge_list[edge_idx];
            edge_weights[edge_idx]++;
            node_dscores[edge.u]++;
            node_dscores[edge.v]++;
        }
    }

    // remove redundant nodes
    for (int node : remove_candidates) {
        if (step_cover[node] && node_dscores[node] == 0) {
            remove_node(node, step_cover_weight, step_cover);
        }
    }

    // update cover if improves
    if (step_cover_weight < best_weight) {
        best_weight = step_cover_weight;
        best_cover = step_cover;
    }
}

int main() {
    ios_base::sync_with_stdio(false);
    cin.tie(NULL);

    auto end_time = high_resolution_clock::now() + milliseconds(1995);
    // read input and initialise structures
    cin >> n >> e;
    AL = vector<vector<pii>>(n, vector<pii>());
    best_cover = vb(n, false);
    node_weights = vi(n, 0);
    node_dscores = vi(n, 0);
    node_conf_change = vb(n, false);
    edge_list = vector<edge>(e);
    node_degrees = vi(n, 0);
    edge_weights = vi(e, 1);
    tabu_list = vb(n, false);
    timestamp = vi(n, 0);
    threshold = (int)(0.5 * n);

    for (int i = 0; i < n; i++) {
        cin >> node_weights[i];
        node_conf_change[i] = true;
    }

    for (int i = 0; i < e; i++) {
        int u, v;
        cin >> u >> v;
        AL[u].push_back(make_pair(v, i));
        AL[v].push_back(make_pair(u, i));
        edge_list[i].u = u;
        edge_list[i].v = v;
        node_degrees[u]++;
        node_degrees[v]++;
    }

    init_mwvc();

    generate_dscores_init();

    generate_remove_candidates();

    step_cover_weight = best_weight;
    step_cover = best_cover;

    while (high_resolution_clock::now() < end_time) {
        step_function();
        step++;
    }

    cout << best_weight << endl;
    for (int i = 0; i < n; i++) {
        if (best_cover[i]) cout << i << " ";
    }
    cout << endl;
    return 0;
}
