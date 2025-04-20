#pragma once
#include <list>
#include "InstanceRead.h"

class Dijkstra{
    private:
        std::list<std::pair<int, double> > * adj;
        InstanceRead *p;
    public:
        Dijkstra(InstanceRead * _p);
	~Dijkstra();
        void Graph(int n_arcs, int n_vertex);
        double Minpath(int n_vertex, int origem, int destino);
        void Path(double* spl); 
};
