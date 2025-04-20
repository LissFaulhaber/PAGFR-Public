#include <iostream>
#include <queue>
#include "Dijkstra.h"
#define INFINITO 1000000

using namespace std;

Dijkstra::Dijkstra(InstanceRead *_p){
    p = _p;
}
Dijkstra::~Dijkstra(){
    delete[] adj;
}
void Dijkstra::Graph(int n_arcs, int n_vertex){
    adj = new list<pair<int, double> >[n_vertex];
    int label, label2;
    for(int a = 0; a < n_arcs; ++a){
        for(int i = 0; i < n_vertex; ++i){
            if(p->v_arcs[a].i == p->v_node[i].id){
                label = i;
            }
            if(p->v_arcs[a].j == p->v_node[i].id){
                label2 = i;
            }
        }
        adj[label].push_back(make_pair(label2, p->v_arcs[a].dist));
    }
}
double Dijkstra::Minpath(int n_vertex, int origem, int destino){
    double cost[n_vertex];
    int visit[n_vertex];
    int orig, dest;
    priority_queue <pair<int, double>, vector<pair<int, double> >, greater<pair<int, double> > > pq;
    for (int i = 0; i < n_vertex; ++i){
        cost[i] = INFINITO;
        visit[i] = false;
        if(p->v_node[i].id == origem){
            orig = i;
        }
        if(p->v_node[i].id == destino){
            dest = i;
        }
    }
    cost[orig] = 0.0;
    pq.push(make_pair(cost[orig], orig));
    while (!pq.empty())
    {
        pair<int, double> p = pq.top();
        int u = p.second;
        pq.pop();
        if(visit[u] == false){
            visit[u] = true;
            list<pair<int, double> >::iterator it;
            for(it = adj[u].begin(); it != adj[u].end(); it++){
                int v = it->first;
                double dist = it->second;
                if(cost[v] > (cost[u] + dist)){
                    cost[v] = cost[u] + dist;
                    pq.push(make_pair(cost[v], v));
                }
            }
        }
    }
    return cost[dest];
}

void Dijkstra::Path(double* spl){
    int n_vertex = p->nvertex;
    int n_arcs = p->narcs;
    int n_OD = p->nOD;
    Graph(n_arcs, n_vertex);
    for(int k = 0; k < n_OD; ++k){
        double result = Minpath(n_vertex, p->v_OD[k].orig, p->v_OD[k].dest);
        spl[k] = result;
        //cout << "Custo menor caminho " << k << ": " << spl[k] << endl;
    }
}
