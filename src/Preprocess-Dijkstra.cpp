#include <iostream>
#include <fstream>
#include "Preprocess-Dijkstra.h"
#include "EVCSLP.h"
#include "EVCSLP-SPL.h"

using namespace std;

Preprocess::Preprocess(InstanceRead *_p){
    p = _p;
}
void Preprocess::ShortestPath(int n_vertex, int n_arcs){
    Dijkstra newinstance(p);
    newinstance.Graph(n_arcs, n_vertex);
    for (int i = 0; i < n_vertex; ++i){
        for (int j = 0; j < n_vertex; ++j){
            if (i != j){
                sp[i][j] = newinstance.Minpath(n_vertex, i, j);
            }
            if (i == j){
                sp[i][j] = 0.0;
            }
        }
    }
}
void Preprocess::ArcsUsed(){
    int n_vertex = p->nvertex;
    int n_arcs = p->narcs;
    int n_OD = p->nOD;
    sp = new double* [n_vertex];
	for(int i = 0; i<n_vertex;++i){
		sp[i] = new double [n_vertex];
	}
    ShortestPath(n_vertex, n_arcs);
    double teste = 0.0;
    set = new int*[n_OD];
    for (int k = 0; k < n_OD; ++k){
        set[k] = new int[n_arcs];
        //cout << "Caminho" << k << ": " << endl;
        for (int a = 0; a < n_arcs; ++a){
            double summation = 0.0;
            summation = sp[p->v_OD[k].orig][p->v_arcs[a].i] + p->v_arcs[a].dist + sp[p->v_arcs[a].j][p->v_OD[k].dest];
            if (fabs(summation - p->spl[k]) < 1.0e-8){
                set[k][a] = 1;
                //cout << a << " ";
            }else{
                set[k][a] = 0;
            }
        }
        //cout << endl;
    }
    ofstream file;
    file.open("graph.txt",ios::trunc);
    for(int k = 0; k < n_OD-1; ++k){
        for(int p=k+1; p< n_OD; ++p){
            for (int a = 0; a < n_arcs; ++a){
                if (set[k][a] == 1 && set[p][a] == 1){
                    file << k << " " << p << endl;
                    break;
                }
            }
        }
    }
    calcCobKruskal(n_OD);
}
int Preprocess::find_set(int v){
    if (v != parent[v])
        parent[v] = find_set(parent[v]);
    return parent[v];
}
void Preprocess::calcCobKruskal(int nOD){
    parent = new int[nOD];
    int rank_[nOD];
    for(int i=0; i<nOD; ++i){
        parent[i] = i;
        rank_[i] = 0;
    }
    int a,b, i, j;
    ifstream file("graph.txt");
    while(!file.eof()){
        file >> i;
        file >> j;
        a = find_set(i);
        b = find_set(j);
        if (a == b) continue;
        if (rank_[a] > rank_[b]) parent[b] = a;
        else parent[a] = b;
        if (rank_[a] == rank_[b]) ++rank_[b];
    }
    for (int k = 0; k < nOD; ++k){
        group.push_back(parent[k]);
    }
    group.sort();
    group.unique();
    cout << group.size() << endl;
}
void Preprocess::NewSolve(){
    list<int>::iterator k;
    InstanceRead *p1 = new InstanceRead(*p);
    int aux = p->nOD;
    int aux2 = p->narcs;
    int aux3 = p->nvertex;
    for (k = group.begin(); k != group.end(); ++k){
        cout << "interator " << *k << endl;
        int origem = 0;
        int contador[aux2];
        int arcos = 0;
        int vert = 0;
        for (int a = 0; a < aux2; ++a){
            contador[a] = 0;
        }
        for (int i = 0; i < aux; ++i){
            if(parent[i] == *k){
                p1->v_OD[origem] = p->v_OD[i];
                p1->spl[origem] = p->spl[i];
                origem++;
                for (int a = 0; a < aux2; ++a){
                    if(set[i][a] > 1E-5){
                        contador[a] += 1;
                    }
                }
            }
        }
        for (int a = 0; a < aux2; ++a){
            if (contador[a] > 1E-5){
                p1->v_arcs[arcos] = p->v_arcs[a];
                arcos++;
            }
        }
        p1->nOD = origem;
        p1->narcs = arcos;
        for (int i = 0; i < aux3; ++i){
            for(int a = 0; a < arcos; ++a){
                if(p->v_node[i].id == p1->v_arcs[a].i || p->v_node[i].id == p1->v_arcs[a].j){
                    p1->v_node[vert] = p->v_node[i];
                    vert++;
                    break;
                }
            }
        }
        p1->nvertex = vert;
        cout << "OD " << origem << " Arcs " << p1->narcs << " Vertex " << p1->nvertex << endl;
        /* EVCSLP j(p1);
        j.solveEVCSLP(); */
        EVCSLPSPL j(p1);
        j.solveSPL();
    }
}
