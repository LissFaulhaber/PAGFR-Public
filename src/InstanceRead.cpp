#include <bits/stdc++.h>
#include "InstanceRead.h"

using namespace std;

//---------------------------------------------------------------------------
InstanceRead::InstanceRead(){}
InstanceRead::~InstanceRead()
{
	//cout << "Objeto destruído" << endl;
	for(int i = 0; i<nvertex;++i){
		delete[] neighbor_edge[i], neighbor_node[i];
	}
	delete[] v_node, v_edge, v_arcs, v_OD, v_cities, neighbor_size, neighbor_edge, neighbor_node, spl;
}
InstanceRead::InstanceRead(const InstanceRead& v){
	nOD = v.nOD;
	nvertex = v.nvertex;
	narcs = v.narcs;
	v_OD = new OD[nOD];
	memcpy(v_OD, v.v_OD, sizeof(OD)*nOD);
	v_arcs = new arcs[narcs];
	memcpy(v_arcs, v.v_arcs, sizeof(arcs)*narcs);
	v_node = new vertex[nvertex];
	memcpy(v_node, v.v_node, sizeof(vertex)*nvertex);
	C = v.C;
	spl = new double[nOD];
	memcpy(spl, v.spl, sizeof(nOD));
	//lpl = v.lpl;
}
void InstanceRead::usage(char *argv[])
{
	cout << "Usage:" <<endl;
	cout << "\t"<<argv[0]<<" <input_instance_name> <tempo_limite> <L> <B>"<<endl;
	cout << "\t"<<"<input_instance_name>: nome do arquivo de entrada"<<endl;
	cout << "\t"<<"<L>: limite maximo de contadores | (0,1) percentual do total de arestas"<<endl;
	cout << "\t"<<"<B>: limite maximo de faixas | (0,1) percentual do total de faixas"<<endl;
}
void InstanceRead::load_data(const char *const file_name){

	ifstream f_inst(file_name);
	if (!f_inst.is_open())
	{
		cout << "ERROR: File " << file_name << " not found!" << endl;
		exit(0);
	}


	int add_edge, add_vert, id, ei, ej, prl, city;
	double dist;
	bool nf_ei, nf_ej;

	f_inst >> nvertex;
	f_inst >> nedges;
	f_inst >> ncities;
	nOD = (ncities*(ncities-1))/2;
	narcs = nedges*2;

	v_node = new vertex [nvertex];
	v_edge = new edge [nedges];
	v_arcs = new arcs [narcs];
	v_OD = new OD [nOD];
	v_cities = new int [nvertex];

	// ----- Rede
	neighbor_edge = new int* [nvertex];
	for(int i = 0; i<nvertex;++i){
		neighbor_edge[i] = new int [nedges];
	}
	neighbor_node = new int* [nvertex];
	for(int i = 0; i<nvertex;++i){
		neighbor_node[i] = new int [nvertex];
	}
	neighbor_size = new int [nedges];

	spl = new double [nOD];
	lpl = new double [nOD];
    	sumM = new double* [nvertex];
	for(int i = 0; i<nvertex;++i){
		sumM[i] = new double [nvertex];
	}
	
	for (int c = 0; c < ncities; ++c){
		f_inst >> city;
		v_cities[c] = city;
	}

	add_vert = 0;
	add_edge = 0;
	C = 0.0;

	for (int e = 0; e < nedges; ++e)
	{
		f_inst >> id;
		f_inst >> ei;
		f_inst >> ej;
		f_inst >> prl;
		f_inst >> dist;

		if(ei == ej) continue;

		nf_ei = nf_ej = true;

		for (int i = 0; i < add_vert; ++i){
			if(ei == v_node[i].label){
				nf_ei = false;
			}
			if(ej == v_node[i].label){
				nf_ej = false;
			}
		}

		if (nf_ei == true){
			v_node[add_vert].label = ei;
			v_node[add_vert].id = add_vert;
			add_vert++;
		}

		if (nf_ej == true){
			v_node[add_vert].label = ej;
			v_node[add_vert].id = add_vert;
			add_vert++;
		}

		int idx, idy;
		for (int i = 0; i < add_vert; ++i)
		{
			if(ei == v_node[i].label){
				idx = i;
			}else if (ej == v_node[i].label){
				idy = i;
			}
		}

		if(idx<idy){
			v_edge[add_edge].i = idx;
			v_edge[add_edge].j = idy;
		}else{
			v_edge[add_edge].i = idy;
			v_edge[add_edge].j = idx;
		}

		neighbor_edge[idx][neighbor_size[idx]] = add_edge;
		neighbor_edge[idy][neighbor_size[idy]] = add_edge;

		neighbor_node[idx][neighbor_size[idx]] = idy;
		neighbor_node[idy][neighbor_size[idy]] = idx;

		neighbor_size[idx]++;
		neighbor_size[idy]++;

		v_edge[add_edge].id = id;
		v_edge[add_edge].prl = prl;
		v_edge[add_edge].dist = dist;
		if (v_edge[add_edge].dist > C){
			C = v_edge[add_edge].dist;
		}
		add_edge++;

	}

	//cout<<"Check: "<<add_vert<<" "<<nvertex<<endl;

	f_inst.close();

	for (int a = 0; a < narcs; ++a){
		int e = a/2;
		v_arcs[a].id = a;
		v_arcs[a].dist = v_edge[e].dist;
		if (a%2 == 0){
			v_arcs[a].i = v_edge[e].i;
			v_arcs[a].j = v_edge[e].j;
		}else{
			v_arcs[a].i = v_edge[e].j;
			v_arcs[a].j = v_edge[e].i;
		}
	}

	int k = 0;
	for (int i = 0; i < ncities-1; ++i){
		for (int j = i+1; j < ncities; ++j){
			for (int v = 0; v < nvertex; ++v){
				if (v_node[v].label == v_cities[i]){
					v_OD[k].orig = v;
				}
				else if (v_node[v].label == v_cities[j]){
					v_OD[k].dest = v;
				}
			}
			v_OD[k].id = k;
			k++;
		}
	}

	//M = 3*C;

//Fim da leitura
}

void InstanceRead::print_masked_data(){
	cout << nvertex << " " << nedges << " " << narcs << " " << nOD << " " << endl;
	//cout << "C" << C << endl;
	for (int i = 0; i < nvertex; ++i)
	{
		cout << v_node[i].id << " " << v_node[i].label << endl;
	}
	/* for (int e = 0; e < nedges; ++e)
	{
		cout << v_edge[e].id << " " << v_edge[e].i << " "
			 << v_edge[e].j << " " << v_edge[e].prl << " "
			 << v_edge[e].dist << endl;
	}
	for (int a = 0; a < narcs; ++a){
		cout << v_arcs[a].id << " " << v_arcs[a].i << " "
			 << v_arcs[a].j << " " << v_arcs[a].dist << endl;
	}
	for (int k = 0; k < nOD; ++k){
		cout << v_OD[k].orig << " " << v_OD[k].dest << endl;
	} */
}
/* void del_dynamic_data(){

} */
