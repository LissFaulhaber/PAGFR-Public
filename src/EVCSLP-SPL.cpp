#include <iostream>
#include "EVCSLP-SPL.h"

using namespace std;

EVCSLPSPL::EVCSLPSPL(InstanceRead *_p){
    p = _p;
}
void EVCSLPSPL::addObjective(IloEnv &env, IloModel &model, IloIntVarArray &y,int n_vertex){
    
    IloExpr summation(env);
    for(int i=0;i<n_vertex; ++i){
        summation +=  y[i];
    }
    model.add(IloMinimize(env,summation));
}
void EVCSLPSPL::addConstraintBatteryorig (IloEnv &env, IloModel &model, IloNumVarMatrix &z, double C, int n_OD, int n_arcs, int n_vertex){
    for(int k=0; k<n_OD; ++k){
        IloExpr summation(env);
        for(int a=0; a<n_arcs; ++a){
            if(p->v_arcs[a].i == p->v_OD[k].orig){
                summation += z[k][a];
            }
        }
        summation -= C;
        IloRange constraintBatteryorig(env, 0, summation, 0);
        stringstream constr;
        constr << "Batteryorig(" << k << ")("<< p->v_OD[k].orig << ")";
        constraintBatteryorig.setName(constr.str().c_str());
        model.add(constraintBatteryorig);
    }
}
void EVCSLPSPL::addConstraintBatterylevel(IloEnv &env, IloModel &model, IloNumVarMatrix &z, IloIntVarMatrix &x, IloNumVarMatrix &w, IloNumVarMatrix &omega, double C, int n_OD, int n_arcs, int n_vertex)
{
    for(int k=0; k<n_OD; ++k){
        for(int i=0; i<n_vertex; ++i){
            if (p->v_node[i].id != p->v_OD[k].orig && p->v_node[i].id != p->v_OD[k].dest)
            {
                IloExpr summation(env);
                IloExpr constraint(env);
                IloExpr charge(env);
                for(int a=0; a<n_arcs; ++a){
                    if(p->v_arcs[a].i == p->v_node[i].id){
                        constraint += z[k][a];
                    }
                    else if(p->v_arcs[a].j == p->v_node[i].id){
                        charge += (z[k][a] - p->v_arcs[a].dist*x[k][a]);
                    }
                }
                summation += charge;
                summation += (omega[k][i]*C);
                summation -= w[k][i];
                constraint -= summation;
                IloRange constraintBatterylevel(env, 0, constraint, 0);
                stringstream constr;
                constr << "Batterylevel(" << k << ")(" << p->v_node[i].id << ")";
                constraintBatterylevel.setName(constr.str().c_str());
                model.add(constraintBatterylevel);
            }
        }
    }
}
void EVCSLPSPL::addConstraintNonnegative (IloEnv &env, IloModel &model, IloNumVarMatrix &z, IloIntVarMatrix &x, int n_OD, int n_arcs, int n_vertex){
    for(int k=0; k<n_OD; ++k){
        for(int a=0; a<n_arcs; ++a){
            if(p->v_arcs[a].j != p->v_OD[k].orig){
                IloExpr constraint(env);
                constraint += z[k][a];
                constraint -= p->v_arcs[a].dist*x[k][a];
                IloRange constraintNonnegative(env, 0, constraint, IloInfinity);
                stringstream constr;
                constr << "Nonnegativity(" << k << ")(" << a << ")";
                constraintNonnegative.setName(constr.str().c_str());
                model.add(constraintNonnegative);
            }
        }
    }
}
void EVCSLPSPL::addConstraintOvercharge (IloEnv &env, IloModel &model, IloNumVarMatrix &z, IloIntVarMatrix &x, double C, int n_OD, int n_arcs, int n_vertex){
    for(int k=0; k<n_OD; ++k){
        for(int a=0; a<n_arcs; ++a){
            IloExpr constraint(env);
            constraint += z[k][a];
            constraint -= C * x[k][a];
            IloRange constraintOvercharge(env, - IloInfinity, constraint, 0);
            stringstream constr;
		    constr << "OverCharge(" << k << ")(" << a << ")";
            constraintOvercharge.setName(constr.str().c_str());
            model.add (constraintOvercharge);
        }
    }
}
void EVCSLPSPL::addConstraintBatterylinear (IloEnv &env, IloModel &model, IloIntVarArray &y, IloNumVarMatrix &w, double C, int n_OD, int n_vertex){
    for(int k=0; k<n_OD; ++k){
        for(int i=0; i<n_vertex; ++i){
            if (p->v_node[i].id != p->v_OD[k].orig && p->v_node[i].id != p->v_OD[k].dest){
                IloExpr constraint(env);
                constraint += w[k][i];
                constraint -= C * y[i];
                IloRange constraintBatterylinear(env, - IloInfinity, constraint, 0);
                stringstream constr;
                constr << "Battery_linear(" << k << ")(" << p->v_node[i].id << ")";
                constraintBatterylinear.setName(constr.str().c_str());
                model.add(constraintBatterylinear);
            }
        }
    }
}
void EVCSLPSPL::addConstraintBatteryLimsup(IloEnv &env, IloModel &model, IloNumVarMatrix &z, IloIntVarMatrix &x, IloNumVarMatrix &w, int n_OD, int n_arcs, int n_vertex)
{
    for(int k=0; k<n_OD; ++k){
        for(int i=0; i<n_vertex; ++i){
            if (p->v_node[i].id != p->v_OD[k].orig && p->v_node[i].id != p->v_OD[k].dest){
                IloExpr summation(env);
                IloExpr constraint (env);
                for (int a=0; a<n_arcs; ++a){
                    if(p->v_arcs[a].j == p->v_node[i].id){
                        summation += z[k][a];
                        summation -= p->v_arcs[a].dist*x[k][a];
                    }
                }
                constraint = w[k][i];
                constraint -= summation;
                IloRange constraintBatteryLimsup(env, - IloInfinity, constraint, 0);
                stringstream constr;
                constr << "Battery_Limsup(" << k << ")(" << p->v_node[i].id << ")";
                constraintBatteryLimsup.setName(constr.str().c_str());
                model.add (constraintBatteryLimsup);
            }
        }
    }
}
void EVCSLPSPL::addConstraintBatteryLiminf (IloEnv &env, IloModel &model,  IloNumVarMatrix &z, IloIntVarMatrix &x, IloIntVarArray &y, IloNumVarMatrix &w, double C, int n_OD, int n_arcs, int n_vertex){
    for(int k=0; k<n_OD; ++k){
        for(int i=0; i<n_vertex; ++i){
            if(p->v_node[i].id != p->v_OD[k].orig && p->v_node[i].id != p->v_OD[k].dest){
                IloExpr charge (env);
                IloExpr summation (env);
                IloExpr constraint (env);
                for(int a=0; a<n_arcs; ++a){
                    if(p->v_arcs[a].j == p->v_node[i].id){
                        summation += z[k][a];
                        summation -= p->v_arcs[a].dist*x[k][a];
                    }
                }
                charge += (1 - y[i]) * C;
                constraint += w[k][i];
                constraint += charge;
                constraint -= summation;
                IloRange constraintBatteryLiminf(env, 0, constraint, IloInfinity);
                stringstream constr;
                constr << "Battery_Liminf(" << k << ")(" << p->v_node[i].id << ")";
                constraintBatteryLiminf.setName(constr.str().c_str());
                model.add (constraintBatteryLiminf);
            }
        }
    }
}
void EVCSLPSPL::addConstraintBatteryOmegalinY (IloEnv &env, IloModel &model, IloNumVarMatrix &omega, IloIntVarArray &y, int n_OD, int n_vertex){
    for(int k=0; k < n_OD; ++k){
        for (int i=0; i<n_vertex; ++i){
            if(p->v_node[i].id!=p->v_OD[k].orig && p->v_node[i].id!=p->v_OD[k].dest){
                IloExpr constraint(env);
                constraint += omega[k][i];
                constraint -= y[i];
                IloRange constraintBatteryOmegalinY(env, -IloInfinity, constraint, 0);
                stringstream constr;
                constr << "Battery_OmegaY(" << k << ")(" << p->v_node[i].id << ")";
                constraintBatteryOmegalinY.setName(constr.str().c_str());
                model.add (constraintBatteryOmegalinY);
            }
        }
    }
}
void EVCSLPSPL::addConstraintBatteryOmegalinX (IloEnv &env, IloModel &model, IloNumVarMatrix &omega, IloIntVarMatrix &x, int n_OD, int n_vertex, int n_arcs){
    for(int k=0; k < n_OD; ++k){
        for (int i=0; i<n_vertex; ++i){
            if(p->v_node[i].id!=p->v_OD[k].orig && p->v_node[i].id!=p->v_OD[k].dest){
                IloExpr constraint(env);
                IloExpr summation(env);
                for(int a=0; a<n_arcs; ++a){
                    if(p->v_arcs[a].j == p->v_node[i].id){
                        summation += x[k][a];
                    }
                }
                constraint += omega[k][i];
                constraint -= summation;
                IloRange constraintBatteryOmegalinX(env, -IloInfinity, constraint, 0);
                stringstream constr;
                constr << "Battery_OmegaX(" << k << ")(" << p->v_node[i].id << ")";
                constraintBatteryOmegalinX.setName(constr.str().c_str());
                model.add (constraintBatteryOmegalinX);
            }
        }
    }
}
void EVCSLPSPL::addConstraintBatteryOmegalinear (IloEnv &env, IloModel &model, IloNumVarMatrix &omega, IloIntVarArray &y, IloIntVarMatrix &x, int n_OD, int n_vertex, int n_arcs){
    for(int k=0; k < n_OD; ++k){
        for (int i=0; i<n_vertex; ++i){
            if(p->v_node[i].id!=p->v_OD[k].orig && p->v_node[i].id!=p->v_OD[k].dest){
                IloExpr constraint(env);
                IloExpr summation(env);
                for(int a=0; a<n_arcs; ++a){
                    if(p->v_arcs[a].j == p->v_node[i].id){
                        summation += x[k][a];
                    }
                }
                constraint += omega[k][i];
                constraint -= (-1 + y[i] + summation);
                IloRange constraintBatteryOmegalinear(env, 0, constraint, IloInfinity);
                stringstream constr;
                constr << "Battery_Omega(" << k << ")(" << p->v_node[i].id << ")";
                constraintBatteryOmegalinear.setName(constr.str().c_str());
                model.add (constraintBatteryOmegalinear);
            }
        }
    }
}
void EVCSLPSPL::addConstraintFlowbalance (IloEnv &env, IloModel &model, IloIntVarMatrix &x, int n_OD, int n_arcs, int n_vertex){
    for(int k=0; k<n_OD; ++k){
        for (int i = 0; i < n_vertex; ++i){
            IloExpr summation (env);
            for(int a=0; a<n_arcs; ++a){
                if(p->v_arcs[a].i == p->v_node[i].id){
                    summation += x[k][a];
                }
                else if(p->v_arcs[a].j == p->v_node[i].id){
                    summation -= x[k][a];
                }
            }
            if (p->v_node[i].id == p->v_OD[k].orig){
                IloRange constraintFloworigin(env, 1, summation, 1);
                stringstream constr;
                constr << "Flow_origin(" << k << ")(" << p->v_node[i].id << ")";
                constraintFloworigin.setName(constr.str().c_str());
                model.add (constraintFloworigin);
            }
            else if (p->v_node[i].id == p->v_OD[k].dest){
                IloRange constraintFlowdest(env, -1, summation, -1);
                stringstream constr;
                constr << "Flow_dest(" << k << ")(" << p->v_node[i].id << ")";
                constraintFlowdest.setName(constr.str().c_str());
                model.add (constraintFlowdest);
            }
            else{
                IloRange constraintFlowother(env, 0, summation, 0);
                stringstream constr;
                constr << "Flow_other(" << k << ")(" << p->v_node[i].id << ")";
                constraintFlowother.setName(constr.str().c_str());
                model.add (constraintFlowother);
            }
        }
    }
}
void EVCSLPSPL::addConstraintPath (IloEnv &env, IloModel &model, IloIntVarMatrix &x, int n_OD, int n_arcs, int n_vertex){
    for(int k=0; k<n_OD; ++k){
        for(int i=0; i<n_vertex; ++i){
            for(int a=0; a<n_arcs; ++a){
                IloExpr summation (env);
                if (p->v_arcs[a].i == p->v_node[i].id){
                    summation += x[k][a];
                    for(int j = 0; j < n_arcs; ++j){
                        if((p->v_arcs[j].j == p->v_arcs[a].i) && (p->v_arcs[j].i == p->v_arcs[a].j)){
                            summation += x[k][j];
                        }
                    }
                    IloRange constraintPath(env, -IloInfinity, summation, 1);
                    stringstream constr;
                    constr << "Path(" << k << ")(" << p->v_node[i].id << ")(" << a <<")";
                    constraintPath.setName(constr.str().c_str());
                    model.add (constraintPath);
                }
            }
        }
    }
}
void EVCSLPSPL::addConstraintShortPath(IloEnv &env, IloModel &model, IloIntVarMatrix &x, int n_OD, int n_arcs){
    for (int k = 0; k < n_OD; ++k){
        IloExpr summation(env);
        for (int a = 0; a < n_arcs; ++a){
            summation += x[k][a] * p->v_arcs[a].dist;
        }
        summation -= p->spl[k];
        IloRange constraintShortPath(env, 0, summation, 0);
        stringstream constr;
        constr << "ShortPath(" << k << ")";
        constraintShortPath.setName(constr.str().c_str());
        model.add (constraintShortPath);
    }
}
void EVCSLPSPL::solveSPL(){
    IloEnv env;
    try{
        IloModel modelSPL(env);

        int n_vertex = p->nvertex;
        int n_OD = p->nOD;
        int n_arcs = p->narcs;
        IloIntVarArray y(env, n_vertex, 0, 1);
        IloNumVarMatrix z(env, n_OD);
        IloIntVarMatrix x(env, n_OD);
        for (int i = 0; i < n_vertex; ++i){
            stringstream vary;
		    vary << "y(" << i << ")";
            y[i].setName(vary.str().c_str());
	        modelSPL.add(y[i]);          
        }
        for(int k = 0; k < n_OD; ++k){
            z[k] = IloNumVarArray (env,n_arcs, 0, IloInfinity,ILOFLOAT);
            x[k] = IloIntVarArray(env,n_arcs, 0, 1);
            for (int a = 0; a < n_arcs; ++a){
                stringstream varz;
                varz << "z(" << k << ")(" << a << ")";
                z[k][a].setName(varz.str().c_str());
                modelSPL.add(z[k][a]);
                stringstream varx;
                varx << "x(" << k << ")(" << a << ")";
                x[k][a].setName(varx.str().c_str());
                modelSPL.add(x[k][a]);
            }
        }
        IloNumVarMatrix w(env, n_OD);
        IloNumVarMatrix omega(env, n_OD);
        for (int k = 0; k < n_OD; ++k){
            w[k] = IloNumVarArray(env, n_vertex,0,IloInfinity,ILOFLOAT);
            omega[k] = IloNumVarArray(env, n_vertex,0,IloInfinity,ILOFLOAT);
            for (int i = 0; i < n_vertex; ++i){
                stringstream varw;
		        varw << "w(" << k << ")(" << i << ")";
                w[k][i].setName(varw.str().c_str());
	            modelSPL.add(w[k][i]);
                stringstream varomega;
		        varomega << "omega(" << k << ")(" << i << ")";
                omega[k][i].setName(varomega.str().c_str());
	            modelSPL.add(omega[k][i]);
            }
        }
        double C = p->C;

        
        addConstraintBatteryorig(env, modelSPL, z, C, n_OD, n_arcs, n_vertex);
        addConstraintBatterylevel(env, modelSPL, z, x, w, omega, C, n_OD, n_arcs, n_vertex);
        addConstraintNonnegative(env, modelSPL, z, x, n_OD, n_arcs, n_vertex);
        addConstraintOvercharge(env, modelSPL, z, x, C, n_OD, n_arcs, n_vertex);
        addConstraintBatterylinear(env, modelSPL, y, w, C, n_OD, n_vertex);
        addConstraintBatteryLimsup(env, modelSPL, z, x, w, n_OD, n_arcs, n_vertex);
        addConstraintBatteryLiminf(env, modelSPL, z, x, y, w, C, n_OD, n_arcs, n_vertex);
        addConstraintBatteryOmegalinY (env, modelSPL, omega, y, n_OD, n_vertex);
        addConstraintBatteryOmegalinX (env, modelSPL, omega, x, n_OD, n_vertex, n_arcs);
        addConstraintBatteryOmegalinear (env, modelSPL, omega, y, x, n_OD, n_vertex, n_arcs);
        addConstraintShortPath(env, modelSPL, x, n_OD, n_arcs);
        addConstraintFlowbalance(env, modelSPL, x, n_OD, n_arcs, n_vertex);
        addConstraintPath(env, modelSPL, x, n_OD, n_arcs, n_vertex);
        addObjective(env, modelSPL, y, n_vertex);

        IloCplex cplex(modelSPL);

        cplex.setParam (IloCplex::Param::TimeLimit, 10800);

        cplex.exportModel ("modelSPL.lp");

        cplex.solve();

        /* for (int k = 0; k < n_OD; ++k){
            double custo = 0;
            for (int a = 0; a < n_arcs; ++a){
                int arest = cplex.getValue(x[k][a]);
                custo = custo + (arest * p->v_arcs[a].dist);
            }
            cout << "Custo " << custo << endl;
        } */        
        cout << "Obj " << cplex.getObjValue() << endl;
        for (int i = 0; i < n_vertex; ++i){
            if (cplex.getValue(y[i]) > 1E-5){
                cout << "Posto " << p->v_node[i].id << endl;
            }
        }
        /* for(int k = 0; k < n_OD; ++k){
            cout << "Arcos usados caminho " << p->v_OD[k].id << endl;
            for (int a = 0; a < n_arcs; ++a){
                if (cplex.getValue(x[k][a]) > 1E-5){
                cout << a << " ";
                }
            }
            cout << endl;
        } */
    }
   catch (IloException& e) {
      cerr << "Concert exception caught: " << e << endl;
   }
   catch (...) {
      cerr << "Unknown exception caught" << endl;
   }
   env.end();
}

