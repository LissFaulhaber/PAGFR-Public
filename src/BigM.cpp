#include <iostream>
#include "BigM.h"

using namespace std;

BigM::BigM(InstanceRead *_p){
    p = _p;
}
void BigM::addObjective(IloEnv &env, IloModel &model, IloIntVarArray &x, int n_arcs){
    IloExpr summation(env);
    for (int a = 0; a < n_arcs; ++a){
        summation += x[a]*p->v_arcs[a].dist;
    }
    model.add(IloMaximize(env,summation));
}
void BigM::addConstraintFlow(IloEnv &env, IloModel &model, IloIntVarArray &x, int origem, int destino, int n_vertex, int n_arcs){
    for (int i = 0; i < n_vertex; ++i){
        IloExpr summation (env);
        for(int a=0; a<n_arcs; ++a){
            if(p->v_arcs[a].i == i){
                summation += x[a];
            }
            else if(p->v_arcs[a].j == i){
                summation -= x[a];
            }
        }
        if (i == origem){
            IloRange constraintFloworigin(env, 1, summation, 1);
            stringstream constr;
            constr << "Flow_origin(" << i << ")";
            constraintFloworigin.setName(constr.str().c_str());
            model.add (constraintFloworigin);
        }
        else if (i == destino){
            IloRange constraintFlowdest(env, -1, summation, -1);
            stringstream constr;
            constr << "Flow_dest(" << i << ")";
            constraintFlowdest.setName(constr.str().c_str());
            model.add (constraintFlowdest);
        }
        else{
            IloRange constraintFlowother(env, 0, summation, 0);
            stringstream constr;
            constr << "Flow_other(" << i << ")";
            constraintFlowother.setName(constr.str().c_str());
            model.add (constraintFlowother);
        }
    }
}
void BigM::addConstraintLongestPath (IloEnv &env, IloModel &model, IloIntVarArray &x, int n_vertex, int n_arcs){
    for(int i=0; i<n_vertex; ++i){
        for(int a=0; a<n_arcs; ++a){
            IloExpr summation (env);
            int j;
            if(p->v_arcs[a].i < p->v_arcs[a].j){
            j = a + 1;
            }
            if (p->v_arcs[a].j < p->v_arcs[a].i){
                j = a - 1;
            }
            if (p->v_arcs[a].i == i){
                summation += x[a] + x[j];
                IloRange constraintPath(env, -IloInfinity, summation, 1);
                stringstream constr;
                constr << "Path(" << i << ")(" << a <<")";
                constraintPath.setName(constr.str().c_str());
                model.add (constraintPath);
            }
        }
    }
}
double BigM::LongestPath(int orig, int dest){
    double lp = 0;
    IloEnv (env);
    try{
        IloModel model(env);

        int origem = orig;
        int destino = dest;
        int n_arcs = p->narcs;
        int n_vertex = p->nvertex;
        IloIntVarArray x(env, n_arcs, 0, 1);
        for (int a = 0; a < n_arcs; ++a){
            stringstream varx;
            varx << "x(" << a << ")";
            x[a].setName(varx.str().c_str());
            model.add(x[a]);
        }


        addConstraintFlow(env, model, x, origem, destino, n_vertex, n_arcs);
        addConstraintLongestPath(env, model, x, n_vertex, n_arcs);
        addObjective(env, model, x, n_arcs);

        IloCplex cplex(model);

        cplex.exportModel ("modelBigM.lp");
        
        ofstream logfile("Output.log");
        cplex.setOut(logfile);
        
        cplex.solve();

        lp = cplex.getObjValue();

    }
    catch (IloException& e) {
      cerr << "Concert exception caught: " << e << endl;
   }
   catch (...) {
      cerr << "Unknown exception caught" << endl;
   }
   env.end();
   return lp;
}
double BigM::MultiplierM(int n_edges){
    double summation = 0;
    double bigest = 0;
    for (int a = 0; a < n_edges; ++a){
        summation += p->v_edge[a].dist;
        if (p->v_edge[a].dist > bigest){
            bigest = p->v_edge[a].dist;
        }
    }
    summation = (summation*bigest);
    return summation;
}
double BigM::SummationM(int i, int j, int n_arcs){
    double summation = 0;
    for (int a = 0; a < n_arcs; ++a){
        summation += p->v_arcs[a].dist;
        if (p->v_arcs[a].i == i && p->v_arcs[a].j == j){
            summation -= p->v_arcs[a].dist;
        }
    }
    return summation;
}
void BigM::Results(double* lpl, double** sumM, double multiM){
    int n_OD = p->nOD;
    int n_arcs = p->narcs;
    int n_vertex = p->nvertex;
    int n_edges = p->nedges;
    for (int k = 0; k <  n_OD; ++k){
        lpl[k] = LongestPath(p->v_OD[k].orig, p->v_OD[k].dest);
    }
    /* for (int k = 0; k < n_OD; ++k){
        cout << lpl[k] << endl;
    } */
    for (int i = 0; i < n_vertex-1; ++i){
        for (int j = i+1; j < n_vertex; ++j){
            sumM[i][j] = SummationM(i, j, n_arcs);
            //cout << sumM[i][j] << endl;
        }
    }
    multiM = MultiplierM(n_edges);
    //cout << multiM << endl;
}
