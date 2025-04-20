#include <iostream>
#include <sys/time.h> 
#include "InstanceRead.h"
#include "BigM.h"
#include "Dijkstra.h"
#include "Preprocess-Dijkstra.h"
#include "EVCSLP-KKT.h"
#include "EVCSLP.h"
#include "EVCSLP-SPL.h"


using namespace std;
int main (int argc,char *argv[]){ 

    InstanceRead *p = new InstanceRead();
    if(argc!=3){
		p->usage(argv);
	}else{ 
		const char *datafile = argv[1];
		int model = atoi(argv[2]);

		p->load_data(datafile);
		cout << datafile << endl;
		/*cout<< "Masked data" << endl;
		p->print_masked_data();
		cout << "Dijkstra resultados" << endl; */
		if(model == 1){
			Dijkstra *d = new Dijkstra(p);
			d->Path(p->spl);
			delete d;
		}
		if(model == 0){
			BigM *g = new BigM(p);
			g->Results(p->lpl, p->sumM, p->multiM);
			delete g;
		}
		if(model == 3){
		Preprocess i(p);
		cout << "Caminhos conexos" << endl;
		i.ArcsUsed();
		i.NewSolve();
		}		
		//cout<< "Teste" << endl;
		clock_t hI = clock();
		struct timeval start, end; 

    	gettimeofday(&start, NULL);
		//cout<< "Original data" << endl;
		//print_original_data();
		if(model == 0){
			EVCSLPKKT f(p);
			cout << "Modelo KKT" << endl;
			f.solveKKT();
		}if(model == 1){
			EVCSLPSPL s(p);
			cout << "Modelo SPL" << endl;
			s.solveSPL();
		}if(model == 2){
			EVCSLP e(p);
			cout << "Modelo Bellmans" << endl;
			e.solveEVCSLP();
		}
		gettimeofday(&end, NULL);
		double time_taken = (end.tv_sec - start.tv_sec) * 1e6;
		time_taken = (time_taken + (end.tv_usec -  start.tv_usec)) * 1e-6;
		cout << "Execution Time: " << time_taken << endl;
		clock_t hF = clock();
		cout<<"Tempo de Execucao: "<<((double)hF - hI) / CLOCKS_PER_SEC<<endl;
		delete p;
		//delete i;
	} 

	return 0;

}
