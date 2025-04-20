#pragma once
#include <vector>
#include <list>
#include "InstanceRead.h"
#include "Dijkstra.h"
using namespace std;

class Preprocess{
    private:
        InstanceRead *p;
        int** set;
        double** sp;
        list<int> group;
        void ShortestPath(int n_vertex, int n_arcs);
        int find_set(int v);
        void calcCobKruskal(int nOD);
    public:
        Preprocess(InstanceRead * _p);
        void ArcsUsed();
        void NewSolve();
        int* parent;
};