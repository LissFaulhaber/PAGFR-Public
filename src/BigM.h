#include "ilcplex/ilocplex.h"
#include "InstanceRead.h"

class BigM{
    public:
        BigM(InstanceRead *_p);
        void Results(double* lpl, double** sumM, double multiM);
    private:
        InstanceRead *p;
        void addObjective(IloEnv &env, IloModel &model, IloIntVarArray &x, int n_arcs);
        void addConstraintFlow(IloEnv &env, IloModel &model, IloIntVarArray &x, int origem, int destino, int n_vertex, int n_arcs);
        void addConstraintLongestPath (IloEnv &env, IloModel &model, IloIntVarArray &x, int n_vertex, int n_arcs);
        double LongestPath(int orig, int dest);
        double MultiplierM(int n_edges);
        double SummationM(int i, int j, int n_arcs);     
};