#include <iostream>
#include <list>
#include "EVCSLP.h"
using namespace std;

/* bool EVCSLP::chkConflicts(IloCplex &Problem)
{
	if ((Problem.getStatus() == IloAlgorithm::Infeasible) || (Problem.getStatus() == IloAlgorithm::InfeasibleOrUnbounded))
	{
		cout << "SolveProblem::chkConflicts()" << "Checking conflict..." << endl;
		cout << "SolveProblem::chkConflicts()" << "No solution!!! Starting conflict refinement..." << endl;

		IloNumArray preferences(Problem.getEnv());
		for (IloInt i = 0; i < constraints.getSize(); i++)
			preferences.add(1.0);

		if (Problem.refineConflict(constraints, preferences))
		{
			IloCplex::ConflictStatusArray conflict = Problem.getConflict(constraints);
			(*Problem.getEnv().getImpl()).useDetailedDisplay(IloTrue);
			cout << "SolveProblem::chkConflicts()" << " Conflict :" << endl;
			for (IloInt i = 0; i < constraints.getSize(); i++)
			{
				if (conflict[i] == IloCplex::ConflictMember)
					cout << "SolveProblem::chkConflicts()" << "Proved  : " << constraints[i]<<"\n";
				if (conflict[i] == IloCplex::ConflictPossibleMember)
					cout << "SolveProblem::chkConflicts()" << "Possible: " << constraints[i]<<"\n";
			}
		}
		else
		{
			cout << "SolveProblem::chkConflicts()" << " Conflict could not be refined." << endl;
		}

		return true;
	}
	return false;
} */

EVCSLP::EVCSLP(InstanceRead *_p){
    p = _p;
}
void EVCSLP::addObjective(IloEnv &env, IloModel &model, IloIntVarArray &y,int n_vertex){
    
    IloExpr summation(env);
    for(int i=0;i<n_vertex; ++i){
        summation +=  y[i];
    }
    model.add(IloMinimize(env,summation));
}
void EVCSLP::addConstraintBatteryorig (IloEnv &env, IloModel &model, IloNumVarMatrix &z, double C, int n_OD, int n_arcs, int n_vertex){
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
        //constraints.add(constraintBatteryorig);
    }
}
void EVCSLP::addConstraintBatterylevel(IloEnv &env, IloModel &model, IloNumVarMatrix &z, IloIntVarMatrix &x, IloNumVarMatrix &w, IloNumVarMatrix &omega, double C, int n_OD, int n_arcs, int n_vertex)
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
                //constraints.add(constraintBatterylevel);
            }
        }
    }
}
void EVCSLP::addConstraintNonnegative (IloEnv &env, IloModel &model, IloNumVarMatrix &z, IloIntVarMatrix &x, int n_OD, int n_arcs, int n_vertex){
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
                //constraints.add(constraintNonnegative);
            }
        }
    }
}
void EVCSLP::addConstraintOvercharge (IloEnv &env, IloModel &model, IloNumVarMatrix &z, IloIntVarMatrix &x, double C, int n_OD, int n_arcs, int n_vertex){
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
            //constraints.add(constraintOvercharge);
        }
    }
}
void EVCSLP::addConstraintBatterylinear (IloEnv &env, IloModel &model, IloIntVarArray &y, IloNumVarMatrix &w, double C, int n_OD, int n_vertex){
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
                //constraints.add(constraintBatterylinear);
            }
        }
    }
}
void EVCSLP::addConstraintBatteryLimsup(IloEnv &env, IloModel &model, IloNumVarMatrix &z, IloIntVarMatrix &x, IloNumVarMatrix &w, int n_OD, int n_arcs, int n_vertex)
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
                //constraints.add(constraintBatteryLimsup);
            }
        }
    }
}
void EVCSLP::addConstraintBatteryLiminf (IloEnv &env, IloModel &model,  IloNumVarMatrix &z, IloIntVarMatrix &x, IloIntVarArray &y, IloNumVarMatrix &w, double C, int n_OD, int n_arcs, int n_vertex){
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
                //constraints.add(constraintBatteryLiminf);
            }
        }
    }
}
void EVCSLP::addConstraintBatteryOmegalinY (IloEnv &env, IloModel &model, IloNumVarMatrix &omega, IloIntVarArray &y, int n_OD, int n_vertex){
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
                //constraints.add(constraintBatteryOmegalinY);
            }
        }
    }
}
void EVCSLP::addConstraintBatteryOmegalinX (IloEnv &env, IloModel &model, IloNumVarMatrix &omega, IloIntVarMatrix &x, int n_OD, int n_vertex, int n_arcs){
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
                //constraints.add(constraintBatteryOmegalinX);
            }
        }
    }
}
void EVCSLP::addConstraintBatteryOmegalinear (IloEnv &env, IloModel &model, IloNumVarMatrix &omega, IloIntVarArray &y, IloIntVarMatrix &x, int n_OD, int n_vertex, int n_arcs){
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
                //constraints.add(constraintBatteryOmegalinear);
            }
        }
    }
}
void EVCSLP::addConstraintFlowbalance (IloEnv &env, IloModel &model, IloIntVarMatrix &x, int n_OD, int n_arcs, int n_vertex){
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
                //constraints.add(constraintFloworigin);
            }
            else if (p->v_node[i].id == p->v_OD[k].dest){
                IloRange constraintFlowdest(env, -1, summation, -1);
                stringstream constr;
                constr << "Flow_dest(" << k << ")(" << p->v_node[i].id << ")";
                constraintFlowdest.setName(constr.str().c_str());
                model.add (constraintFlowdest);
                //constraints.add(constraintFlowdest);
            }
            else{
                IloRange constraintFlowother(env, 0, summation, 0);
                stringstream constr;
                constr << "Flow_other(" << k << ")(" << p->v_node[i].id << ")";
                constraintFlowother.setName(constr.str().c_str());
                model.add (constraintFlowother);
                //constraints.add(constraintFlowother);
            }
        }
    }
}
void EVCSLP::addConstraintPath (IloEnv &env, IloModel &model, IloIntVarMatrix &x, int n_OD, int n_arcs, int n_vertex){
    for(int k=0; k<n_OD; ++k){
        for(int i=0; i<n_vertex; ++i){
            for(int a=0; a<n_arcs; ++a){
                IloExpr summation (env);
                /* int j;
                if(p->v_arcs[a].i < p->v_arcs[a].j){
                    j = a + 1;
                }
                if (p->v_arcs[a].j < p->v_arcs[a].i){
                    j = a - 1;
                } */
                if (p->v_arcs[a].i == p->v_node[i].id){
                    summation += x[k][a];
                    for(int j = 0; j < n_arcs; ++j){
                        if((p->v_arcs[j].j == p->v_arcs[a].i) && (p->v_arcs[j].i == p->v_arcs[a].j)){
                            summation += x[k][j];
                        }
                    }
                    IloRange constraintPath(env, -IloInfinity, summation, 1);
                    stringstream constr;
                    constr << "Path(" << k << ")(" << i << ")(" << a <<")";
                    constraintPath.setName(constr.str().c_str());
                    model.add (constraintPath);
                    //constraints.add(constraintPath);
                }
            }
        }
    }
}
void EVCSLP::addConstraintFME (IloEnv &env, IloModel &model, IloIntVarMatrix &x, IloNumVarMatrix &pi, int n_OD, int n_arcs, int n_vertex){
    for(int k=0; k<n_OD; ++k){
        for(int a=0; a<n_arcs; ++a){
            IloExpr summation (env);
            IloExpr constraint(env);
            /* int f;
            if(p->v_arcs[a].i < p->v_arcs[a].j){
                f = a + 1;
            }
            if (p->v_arcs[a].j < p->v_arcs[a].i){
                f = a - 1;
            } */
            summation += p->v_arcs[a].dist;
            for(int j = 0; j < n_arcs; ++j){
                if(p->v_arcs[j].j == p->v_arcs[a].i && p->v_arcs[j].i == p->v_arcs[a].j){
                    summation -=(2*p->v_arcs[a].dist*x[k][j]);
                }
            }
            for(int i = 0; i < n_vertex; ++i){
                if(p->v_node[i].id == p->v_arcs[a].i){
                    constraint += pi[k][i];
                }
                if(p->v_node[i].id == p->v_arcs[a].j){
                    constraint -= pi[k][i];
                }
            }
            constraint -= summation;
            IloRange constraintFME(env, -IloInfinity, constraint, 0);
            stringstream constr;
            constr << "FMELin(" << k << ")(" << a << ")";
            constraintFME.setName(constr.str().c_str());
            model.add (constraintFME);
            //constraints.add(constraintFME);
        }
    }
}
void EVCSLP::addConstraintFMEPi (IloEnv &env, IloModel &model, IloNumVarMatrix &pi, int n_OD, int n_vertex){
    for(int k=0; k<n_OD; ++k){
        IloExpr constraint (env);
        for(int i = 0; i < n_vertex; ++i){
            if(p->v_node[i].id == p->v_OD[k].dest){
                constraint = pi[k][i];
            }
        }
        IloRange constraintFMEPi(env, 0, constraint, 0);
        stringstream constr;
        constr << "Pidest(" << k << ")(" << p->v_OD[k].dest << ")";
        constraintFMEPi.setName(constr.str().c_str());
        model.add (constraintFMEPi);
        //constraints.add(constraintFMEPi);
    }   
}

void EVCSLP::solveEVCSLP(){
    IloEnv env;
    try{
        IloModel model(env);
        //constraints = IloConstraintArray(env);

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
	        model.add(y[i]);          
        }
        for(int k = 0; k < n_OD; ++k){
            z[k] = IloNumVarArray (env,n_arcs, 0, IloInfinity,ILOFLOAT);
            x[k] = IloIntVarArray(env,n_arcs, 0, 1);
            for (int a = 0; a < n_arcs; ++a){
                stringstream varz;
                varz << "z(" << k << ")(" << a << ")";
                z[k][a].setName(varz.str().c_str());
                model.add(z[k][a]);
                stringstream varx;
                varx << "x(" << k << ")(" << a << ")";
                x[k][a].setName(varx.str().c_str());
                model.add(x[k][a]);
            }
        }
        IloNumVarMatrix pi(env, n_OD);
        IloNumVarMatrix w(env, n_OD);
        IloNumVarMatrix omega(env, n_OD);
        for (int k = 0; k < n_OD; ++k){
            pi[k] = IloNumVarArray(env, n_vertex, 0, IloInfinity, ILOFLOAT);
            w[k] = IloNumVarArray(env, n_vertex, 0, IloInfinity, ILOFLOAT);
            omega[k] = IloNumVarArray(env, n_vertex, 0, IloInfinity, ILOFLOAT);
            for (int i = 0; i < n_vertex; ++i){
                stringstream varpi;
		        varpi << "pi(" << k << ")(" << i << ")";
                pi[k][i].setName(varpi.str().c_str());
	            model.add(pi[k][i]);
                stringstream varw;
		        varw << "w(" << k << ")(" << i << ")";
                w[k][i].setName(varw.str().c_str());
	            model.add(w[k][i]);
                stringstream varomega;
		        varomega << "omega(" << k << ")(" << i << ")";
                omega[k][i].setName(varomega.str().c_str());
	            model.add(omega[k][i]);
            }
        }
        double C = p->C;

        
        addConstraintBatteryorig(env, model, z, C, n_OD, n_arcs, n_vertex);
        addConstraintBatterylevel(env, model, z, x, w, omega, C, n_OD, n_arcs, n_vertex);
        addConstraintNonnegative(env, model, z, x, n_OD, n_arcs, n_vertex);
        addConstraintOvercharge(env, model, z, x, C, n_OD, n_arcs, n_vertex);
        addConstraintBatterylinear(env, model, y, w, C, n_OD, n_vertex);
        addConstraintBatteryLimsup(env, model, z, x, w, n_OD, n_arcs, n_vertex);
        addConstraintBatteryLiminf(env, model, z, x, y, w, C, n_OD, n_arcs, n_vertex);
        addConstraintBatteryOmegalinY (env, model, omega, y, n_OD, n_vertex);
        addConstraintBatteryOmegalinX (env, model, omega, x, n_OD, n_vertex, n_arcs);
        addConstraintBatteryOmegalinear (env, model, omega, y, x, n_OD, n_vertex, n_arcs);
        addConstraintFlowbalance(env, model, x, n_OD, n_arcs, n_vertex);
        addConstraintPath(env, model, x, n_OD, n_arcs, n_vertex);
        addConstraintFME(env, model, x, pi, n_OD, n_arcs, n_vertex);
        addConstraintFMEPi(env, model, pi, n_OD, n_vertex);
        addObjective(env, model, y, n_vertex);

        IloCplex cplex(model);

        cplex.setParam (IloCplex::Param::TimeLimit, 10800);

        
        cplex.exportModel("model.lp");
        //cplex.setParam(IloCplex::Param::Threads, 1);

        cplex.solve();
        //chkConflicts(cplex);
   
        cout << "Obj " << cplex.getObjValue() << endl;
        for (int i = 0; i < n_vertex; ++i){
            if (cplex.getValue(y[i]) > 1E-5){
                cout << "Posto " << p->v_node[i].id << endl;
            }
        }
        /* list <int> arco;
        for(int k = 0; k < n_OD; ++k){
            //double count = 0.0;
            //cout << "Arcos usados caminho " << p->v_OD[k].id << endl;
            for (int a = 0; a < n_arcs; ++a){
                if (cplex.getValue(x[k][a]) > 1E-5){
                arco.push_back(p->v_arcs[a].id);
                //count += p->v_arcs[a].dist;
                }
            }
            //cout << "Total caminho " << count << endl;
        }
        arco.sort();
        arco.unique();
        cout << arco.size() << endl; 
        list<int>::iterator k;
        for (k = arco.begin(); k != arco.end(); ++k){
            cout << *k << " ";
        }
        cout << endl; */
    }
   catch (IloException& e) {
      cerr << "Concert exception caught: " << e << endl;
   }
   catch (...) {
      cerr << "Unknown exception caught" << endl;
   }
   env.end();
}
