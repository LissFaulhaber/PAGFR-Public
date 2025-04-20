#include <iostream>
#include "EVCSLP-KKT.h"
using namespace std;

/* bool EVCSLPKKT::chkConflicts(IloCplex &Problem)
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

EVCSLPKKT::EVCSLPKKT(InstanceRead *_p){
    p = _p;
}
void EVCSLPKKT::addObjective(IloEnv &env, IloModel &model, IloIntVarArray &y,int n_vertex){
    
    IloExpr summation(env);
    for(int i=0;i<n_vertex; ++i){
        summation +=  y[i];
    }
    model.add(IloMinimize(env,summation));
}
void EVCSLPKKT::addConstraintBatteryorig (IloEnv &env, IloModel &model, IloNumVarMatrix &z, double C, int n_OD, int n_arcs, int n_vertex){
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
void EVCSLPKKT::addConstraintBatterylevel(IloEnv &env, IloModel &model, IloNumVarMatrix &z, IloIntVarMatrix &x, IloNumVarMatrix &w, IloNumVarMatrix &omega, double C, int n_OD, int n_arcs, int n_vertex)
{
    for(int k=0; k<n_OD; ++k){
        for(int i=0; i<n_vertex; ++i){
            if (i != p->v_OD[k].orig && i != p->v_OD[k].dest)
            {
                IloExpr summation(env);
                IloExpr constraint(env);
                IloExpr charge(env);
                for(int a=0; a<n_arcs; ++a){
                    if(p->v_arcs[a].i == i){
                        constraint += z[k][a];
                    }
                    else if(p->v_arcs[a].j == i){
                        charge += (z[k][a] - p->v_arcs[a].dist*x[k][a]);
                    }
                }
                summation += charge;
                summation += (omega[k][i]*C);
                summation -= w[k][i];
                constraint -= summation;
                IloRange constraintBatterylevel(env, 0, constraint, 0);
                stringstream constr;
                constr << "Batterylevel(" << k << ")(" << i << ")";
                constraintBatterylevel.setName(constr.str().c_str());
                model.add(constraintBatterylevel);
                //constraints.add(constraintBatterylevel);
            }
        }
    }
}
void EVCSLPKKT::addConstraintNonnegative (IloEnv &env, IloModel &model, IloNumVarMatrix &z, IloIntVarMatrix &x, int n_OD, int n_arcs){
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
void EVCSLPKKT::addConstraintOvercharge (IloEnv &env, IloModel &model, IloNumVarMatrix &z, IloIntVarMatrix &x, double C, int n_OD, int n_arcs){
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
void EVCSLPKKT::addConstraintBatterylinear (IloEnv &env, IloModel &model, IloIntVarArray &y, IloNumVarMatrix &w, double C, int n_OD, int n_vertex){
    for(int k=0; k<n_OD; ++k){
        for(int i=0; i<n_vertex; ++i){
            if (i != p->v_OD[k].orig && i != p->v_OD[k].dest){
                IloExpr constraint(env);
                constraint += w[k][i];
                constraint -= C * y[i];
                IloRange constraintBatterylinear(env, - IloInfinity, constraint, 0);
                stringstream constr;
                constr << "Battery_linear(" << k << ")(" << i << ")";
                constraintBatterylinear.setName(constr.str().c_str());
                model.add(constraintBatterylinear);
                //constraints.add(constraintBatterylinear);
            }
        }
    }
}
void EVCSLPKKT::addConstraintBatteryLimsup(IloEnv &env, IloModel &model, IloNumVarMatrix &z, IloIntVarMatrix &x, IloNumVarMatrix &w, int n_OD, int n_arcs, int n_vertex)
{
    for(int k=0; k<n_OD; ++k){
        for(int i=0; i<n_vertex; ++i){
            if (i != p->v_OD[k].orig && i != p->v_OD[k].dest){
                IloExpr summation(env);
                IloExpr constraint (env);
                for (int a=0; a<n_arcs; ++a){
                    if(p->v_arcs[a].j == i){
                        summation += z[k][a];
                        summation -= p->v_arcs[a].dist*x[k][a];
                    }
                }
                constraint = w[k][i];
                constraint -= summation;
                IloRange constraintBatteryLimsup(env, - IloInfinity, constraint, 0);
                stringstream constr;
                constr << "Battery_Limsup(" << k << ")(" << i << ")";
                constraintBatteryLimsup.setName(constr.str().c_str());
                model.add (constraintBatteryLimsup);
                //constraints.add(constraintBatteryLimsup);
            }
        }
    }
}
void EVCSLPKKT::addConstraintBatteryLiminf (IloEnv &env, IloModel &model,  IloNumVarMatrix &z, IloIntVarMatrix &x, IloIntVarArray &y, IloNumVarMatrix &w, double C, int n_OD, int n_arcs, int n_vertex){
    for(int k=0; k<n_OD; ++k){
        for(int i=0; i<n_vertex; ++i){
            if(i != p->v_OD[k].orig && i != p->v_OD[k].dest){
                IloExpr charge (env);
                IloExpr summation (env);
                IloExpr constraint (env);
                for(int a=0; a<n_arcs; ++a){
                    if(p->v_arcs[a].j == i){
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
                constr << "Battery_Liminf(" << k << ")(" << i << ")";
                constraintBatteryLiminf.setName(constr.str().c_str());
                model.add (constraintBatteryLiminf);
                //constraints.add(constraintBatteryLiminf);
            }
        }
    }
}
void EVCSLPKKT::addConstraintBatteryOmegalinY (IloEnv &env, IloModel &model, IloNumVarMatrix &omega, IloIntVarArray &y, int n_OD, int n_vertex){
    for(int k=0; k < n_OD; ++k){
        for (int i=0; i<n_vertex; ++i){
            if(i!=p->v_OD[k].orig && i!=p->v_OD[k].dest){
                IloExpr constraint(env);
                constraint += omega[k][i];
                constraint -= y[i];
                IloRange constraintBatteryOmegalinY(env, -IloInfinity, constraint, 0);
                stringstream constr;
                constr << "Battery_OmegaY(" << k << ")(" << i << ")";
                constraintBatteryOmegalinY.setName(constr.str().c_str());
                model.add (constraintBatteryOmegalinY);
                //constraints.add(constraintBatteryOmegalinY);
            }
        }
    }
}
void EVCSLPKKT::addConstraintBatteryOmegalinX (IloEnv &env, IloModel &model, IloNumVarMatrix &omega, IloIntVarMatrix &x, int n_OD, int n_vertex, int n_arcs){
    for(int k=0; k < n_OD; ++k){
        for (int i=0; i<n_vertex; ++i){
            if(i!=p->v_OD[k].orig && i!=p->v_OD[k].dest){
                IloExpr constraint(env);
                IloExpr summation(env);
                for(int a=0; a<n_arcs; ++a){
                    if(p->v_arcs[a].j == i){
                        summation += x[k][a];
                    }
                }
                constraint += omega[k][i];
                constraint -= summation;
                IloRange constraintBatteryOmegalinX(env, -IloInfinity, constraint, 0);
                stringstream constr;
                constr << "Battery_OmegaX(" << k << ")(" << i << ")";
                constraintBatteryOmegalinX.setName(constr.str().c_str());
                model.add (constraintBatteryOmegalinX);
                //constraints.add(constraintBatteryOmegalinX);
            }
        }
    }
}
void EVCSLPKKT::addConstraintBatteryOmegalinear (IloEnv &env, IloModel &model, IloNumVarMatrix &omega, IloIntVarArray &y, IloIntVarMatrix &x, int n_OD, int n_vertex, int n_arcs){
    for(int k=0; k < n_OD; ++k){
        for (int i=0; i<n_vertex; ++i){
            if(i!=p->v_OD[k].orig && i!=p->v_OD[k].dest){
                IloExpr constraint(env);
                IloExpr summation(env);
                for(int a=0; a<n_arcs; ++a){
                    if(p->v_arcs[a].j == i){
                        summation += x[k][a];
                    }
                }
                constraint += omega[k][i];
                constraint -= (-1 + y[i] + summation);
                IloRange constraintBatteryOmegalinear(env, 0, constraint, IloInfinity);
                stringstream constr;
                constr << "Battery_Omega(" << k << ")(" << i << ")";
                constraintBatteryOmegalinear.setName(constr.str().c_str());
                model.add (constraintBatteryOmegalinear);
                //constraints.add(constraintBatteryOmegalinear);
            }
        }
    }
}
void EVCSLPKKT::addConstraintFlowbalance (IloEnv &env, IloModel &model, IloIntVarMatrix &x, int n_OD, int n_arcs, int n_vertex){
    for(int k=0; k<n_OD; ++k){
        for (int i = 0; i < n_vertex; ++i){
            IloExpr summation (env);
            for(int a=0; a<n_arcs; ++a){
                if(p->v_arcs[a].i == i){
                    summation += x[k][a];
                }
                else if(p->v_arcs[a].j == i){
                    summation -= x[k][a];
                }
            }
            if (i == p->v_OD[k].orig){
                IloRange constraintFloworigin(env, 1, summation, 1);
                stringstream constr;
                constr << "Flow_origin(" << k << ")(" << i << ")";
                constraintFloworigin.setName(constr.str().c_str());
                model.add (constraintFloworigin);
                //constraints.add(constraintFloworigin);
            }
            else if (i == p->v_OD[k].dest){
                IloRange constraintFlowdest(env, -1, summation, -1);
                stringstream constr;
                constr << "Flow_dest(" << k << ")(" << i << ")";
                constraintFlowdest.setName(constr.str().c_str());
                model.add (constraintFlowdest);
                //constraints.add(constraintFlowdest);
            }
            else{
                IloRange constraintFlowother(env, 0, summation, 0);
                stringstream constr;
                constr << "Flow_other(" << k << ")(" << i << ")";
                constraintFlowother.setName(constr.str().c_str());
                model.add (constraintFlowother);
                //constraints.add(constraintFlowother);
            }
        }
    }
}
void EVCSLPKKT::addConstraintPath (IloEnv &env, IloModel &model, IloIntVarMatrix &x, int n_OD, int n_arcs, int n_vertex){
    for(int k=0; k<n_OD; ++k){
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
                    summation += x[k][a] + x[k][j];
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
void EVCSLPKKT::addConstraintKKTMi (IloEnv &env, IloModel &model, IloNumVarMatrix &v, IloNumVarMatrix &lamb, IloNumVarMatrix &mi, int n_OD, int n_arcs){
    for(int k=0; k<n_OD; ++k){
        for(int a=0; a<n_arcs; ++a){
            IloExpr constraint (env);
            constraint = p->v_arcs[a].dist - mi[k][p->v_arcs[a].i] + mi[k][p->v_arcs[a].j] - v[k][a] + lamb[k][a];
            IloRange constraintKKT(env, 0, constraint, 0);
            stringstream constr;
            constr << "KKTLinMi(" << k << ")(" << a << ")";
            constraintKKT.setName(constr.str().c_str());
            model.add (constraintKKT);
            //constraints.add(constraintKKT);
        }
    }
}
void EVCSLPKKT::addConstraintKKTV (IloEnv &env, IloModel &model, IloIntVarMatrix &x, IloNumVarMatrix &v, int n_OD, int n_arcs){
    for(int k=0; k<n_OD; ++k){
        for(int a=0; a<n_arcs; ++a){
            IloExpr summation (env);
            summation += v[k][a];
            summation -= p->lpl[k]*(1-x[k][a]);
            IloRange constraintKKTV(env, -IloInfinity, summation, 0);
            stringstream constr;
            constr << "KKTLinV(" << k << ")(" << a << ")";
            constraintKKTV.setName(constr.str().c_str());
            model.add (constraintKKTV);
            //constraints.add(constraintKKTV);
        }
    }
}
void EVCSLPKKT::addConstraintKKTLamb (IloEnv &env, IloModel &model, IloIntVarMatrix &x, IloNumVarMatrix &lamb, int n_OD, int n_arcs){
    for(int k=0; k<n_OD; ++k){
        for(int a=0; a<n_arcs; ++a){
            IloExpr summation (env);
            summation += lamb[k][a];
            summation -= p->lpl[k]*(1-(1-x[k][a]));
            IloRange constraintKKTLamb(env, -IloInfinity, summation, 0);
            stringstream constr;
            constr << "KKTLinLamb(" << k << ")(" << a << ")";
            constraintKKTLamb.setName(constr.str().c_str());
            model.add (constraintKKTLamb);
            //constraints.add(constraintKKTLamb);
        }
    }
}
void EVCSLPKKT::solveKKT(){
   IloEnv env;
    try{
        IloModel modelkkt(env);
        //constraints = IloConstraintArray(env);
        int n_vertex = p->nvertex;
        int n_OD = p->nOD;
        int n_arcs = p->narcs;
        IloIntVarArray y(env, n_vertex, 0, 1);
        IloIntVarMatrix x(env, n_OD);
        IloNumVarMatrix z(env, n_OD);
        IloNumVarMatrix lamb(env, n_OD);
        IloNumVarMatrix v(env, n_OD);
        for (int i = 0; i < n_vertex; ++i){
            stringstream vary;
		    vary << "y(" << i << ")";
            y[i].setName(vary.str().c_str());
	        modelkkt.add(y[i]);          
        }
        for(int k = 0; k < n_OD; ++k){
            x[k] = IloIntVarArray(env,n_arcs, 0, 1);
            z[k] = IloNumVarArray (env,n_arcs, 0, IloInfinity,ILOFLOAT);
            lamb[k] = IloNumVarArray (env, n_arcs, 0, IloInfinity);
            v[k] = IloNumVarArray (env, n_arcs, 0, IloInfinity);
            for (int a = 0; a < n_arcs; ++a){
                stringstream varx;
                varx << "x(" << k << ")(" << a << ")";
                x[k][a].setName(varx.str().c_str());
                modelkkt.add(x[k][a]);
                stringstream varz;
                varz << "z(" << k << ")(" << a << ")";
                z[k][a].setName(varz.str().c_str());
                modelkkt.add(z[k][a]);
                stringstream varlamb;
                varlamb << "lamb(" << k << ")(" << a << ")";
                lamb[k][a].setName(varlamb.str().c_str());
                modelkkt.add(lamb[k][a]);
                stringstream varv;
                varv << "v(" << k << ")(" << a << ")";
                v[k][a].setName(varv.str().c_str());
                modelkkt.add(v[k][a]);
            }
        }
        IloNumVarMatrix mi(env, n_OD);
        IloNumVarMatrix w(env, n_OD);
        IloNumVarMatrix omega(env, n_OD);
        for (int k = 0; k < n_OD; ++k){
            mi[k] = IloNumVarArray(env, n_vertex,-IloInfinity,IloInfinity,ILOFLOAT);
            w[k] = IloNumVarArray(env, n_vertex,0,IloInfinity,ILOFLOAT);
            omega[k] = IloNumVarArray(env, n_vertex,0,IloInfinity,ILOFLOAT);
            for (int i = 0; i < n_vertex; ++i){
                stringstream varmi;
		        varmi << "mi(" << k << ")(" << i << ")";
                mi[k][i].setName(varmi.str().c_str());
	            modelkkt.add(mi[k][i]);
                stringstream varw;
		        varw << "w(" << k << ")(" << i << ")";
                w[k][i].setName(varw.str().c_str());
	            modelkkt.add(w[k][i]);
                stringstream varomega;
		        varomega << "omega(" << k << ")(" << i << ")";
                omega[k][i].setName(varomega.str().c_str());
	            modelkkt.add(omega[k][i]);
            }
        }
        double C = p->C;
        
        addConstraintBatteryorig(env, modelkkt, z, C, n_OD, n_arcs, n_vertex);
        addConstraintBatterylevel(env, modelkkt, z, x, w, omega, C, n_OD, n_arcs, n_vertex);
        addConstraintNonnegative(env, modelkkt, z, x, n_OD, n_arcs);
        addConstraintOvercharge(env, modelkkt, z, x, C, n_OD, n_arcs);
        addConstraintBatterylinear(env, modelkkt, y, w, C, n_OD, n_vertex);
        addConstraintBatteryLimsup(env, modelkkt, z, x, w, n_OD, n_arcs, n_vertex);
        addConstraintBatteryLiminf(env, modelkkt, z, x, y, w, C, n_OD, n_arcs, n_vertex);
        addConstraintBatteryOmegalinY (env, modelkkt, omega, y, n_OD, n_vertex);
        addConstraintBatteryOmegalinX (env, modelkkt, omega, x, n_OD, n_vertex, n_arcs);
        addConstraintBatteryOmegalinear (env, modelkkt, omega, y, x, n_OD, n_vertex, n_arcs);
        addConstraintFlowbalance(env, modelkkt, x, n_OD, n_arcs, n_vertex);
        addConstraintPath(env, modelkkt, x, n_OD, n_arcs, n_vertex);
        addConstraintKKTMi (env, modelkkt, v, lamb, mi, n_OD, n_arcs);
        addConstraintKKTV (env, modelkkt, x, v, n_OD, n_arcs);
        addConstraintKKTLamb (env, modelkkt, x, lamb, n_OD, n_arcs);
        addObjective(env, modelkkt, y, n_vertex);

        IloCplex cplex(modelkkt);

        cplex.setParam (IloCplex::Param::TimeLimit, 10800);

        cplex.exportModel ("modelKKT.lp");

        cplex.solve();
        //chkConflicts(cplex);
        cout << "Obj " << cplex.getObjValue() << endl;
        for (int i = 0; i < n_vertex; ++i){
            if (cplex.getValue(y[i]) > 0){
                cout << "Posto " << i << endl;
            }
        }
        /* for (int k = 0; k < n_OD; ++k){
            double custo = 0;
            for (int a = 0; a < n_arcs; ++a){
                int arest = cplex.getValue(x[k][a]);
                custo = custo + (arest * p->v_arcs[a].dist);
                //cout << cplex.getValue(z[k][a]) << ",";
            }
            //cout << endl;
            cout << "Custo " << custo << endl;
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
