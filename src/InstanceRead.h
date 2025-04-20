#pragma once
#define MAX_NOS 2500 // numero maximo de nos
#define MAX_ARE 2600 // numero maximo de arestas
#define MAX_OD 325000 // numero maximo de pares OD


class InstanceRead{
    public:

        // ----------------- ESTRUTURAS ----------------
        struct edge{
	    int i, j, id, prl;
	    double dist;
        };

        struct arcs{
	    int i, j, id;
	    double dist;
        };

        struct OD{
            int id, orig, dest;
        };

        struct vertex{
            int id, label;
        };

        // ----- Dados de entrada
        int nvertex;
        int nedges;
        int ncities;
        int nOD;
        int narcs;
        double C;
        //double M;

        //Mapear ID x ID real
        vertex* v_node;
        edge* v_edge;
        arcs* v_arcs;
        OD* v_OD;
        int* v_cities; 
        // ----- Rede
        int** neighbor_edge;
        int** neighbor_node;
        int* neighbor_size;

        InstanceRead();
        InstanceRead(const InstanceRead& v);
        ~InstanceRead();

        void usage(char *argv[]);
        void load_data(const char *const file_name);
        void print_masked_data();

        double* spl;
        double* lpl;
        double** sumM;
        double multiM;
        
};