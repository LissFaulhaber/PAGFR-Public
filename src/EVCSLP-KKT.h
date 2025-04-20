#include "ilcplex/ilocplex.h"
#include "InstanceRead.h"

typedef IloArray<IloNumVarArray> IloNumVarMatrix;
typedef IloArray<IloIntVarArray> IloIntVarMatrix;
typedef IloArray<IloIntArray> IloIntMatrix;

class EVCSLPKKT{
    public:
        EVCSLPKKT(InstanceRead *_p);
        void solveKKT();
    private:
        void addObjective(IloEnv &env, IloModel &model, IloIntVarArray &y,int n_vertex);
        void addConstraintBatteryorig(IloEnv &env, IloModel &model, IloNumVarMatrix &z, double C, int n_OD, int n_arcs, int n_vertex);
        void addConstraintBatterylevel(IloEnv &env, IloModel &model, IloNumVarMatrix &z, IloIntVarMatrix &x, IloNumVarMatrix &w, IloNumVarMatrix &omega, double C, int n_OD, int n_arcs, int n_vertex);
        void addConstraintNonnegative (IloEnv &env, IloModel &model, IloNumVarMatrix &z, IloIntVarMatrix &x, int n_OD, int n_arcs);
        void addConstraintOvercharge (IloEnv &env, IloModel &model, IloNumVarMatrix &z, IloIntVarMatrix &x, double C, int n_OD, int n_arcs);
        void addConstraintBatterylinear (IloEnv &env, IloModel &model, IloIntVarArray &y, IloNumVarMatrix &w, double C, int n_OD, int n_vertex);
        void addConstraintBatteryLimsup (IloEnv &env, IloModel &model, IloNumVarMatrix &z, IloIntVarMatrix &x, IloNumVarMatrix &w, int n_OD, int n_arcs, int n_vertex);
        void addConstraintBatteryLiminf (IloEnv &env, IloModel &model,  IloNumVarMatrix &z, IloIntVarMatrix &x, IloIntVarArray &y, IloNumVarMatrix &w, double C, int n_OD, int n_arcs, int n_vertex);
        void addConstraintBatteryOmegalinY (IloEnv &env, IloModel &model, IloNumVarMatrix &omega, IloIntVarArray &y, int n_OD, int n_vertex);
        void addConstraintBatteryOmegalinX (IloEnv &env, IloModel &model, IloNumVarMatrix &omega, IloIntVarMatrix &x, int n_OD, int n_vertex, int n_arcs);
        void addConstraintBatteryOmegalinear (IloEnv &env, IloModel &model, IloNumVarMatrix &omega, IloIntVarArray &y, IloIntVarMatrix &x, int n_OD, int n_vertex, int n_arcs);
        void addConstraintFlowbalance (IloEnv &env, IloModel &model, IloIntVarMatrix &x, int n_OD, int n_arcs, int n_vertex);
        void addConstraintPath (IloEnv &env, IloModel &model, IloIntVarMatrix &x, int n_OD, int n_arcs, int n_vertex);
        void addConstraintKKTMi (IloEnv &env, IloModel &model, IloNumVarMatrix &v, IloNumVarMatrix &lamb, IloNumVarMatrix &mi, int n_OD, int n_arcs);
        void addConstraintKKTV (IloEnv &env, IloModel &model, IloIntVarMatrix &x, IloNumVarMatrix &v, int n_OD, int n_arcs);
        void addConstraintKKTLamb (IloEnv &env, IloModel &model, IloIntVarMatrix &x, IloNumVarMatrix &lamb, int n_OD, int n_arcs);
        InstanceRead *p;
        //IloConstraintArray constraints;
        //bool chkConflicts(IloCplex& Problem);
};