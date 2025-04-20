####################################################
###################### MAKE ########################
####################################################

EXECUTAVEL = EVCSLP
PATHEXEC = ./bin
PATHSRC= ./src
PATHTEMP = ./.temp

all:
	mkdir -p $(PATHEXEC)
	mkdir -p $(PATHTEMP)
	make $(EXECUTAVEL)

#Juntando todos os objetos e gerando o programa
$(EXECUTAVEL): $(PATHTEMP)/InstanceRead.o $(PATHTEMP)/BigM.o $(PATHTEMP)/Dijkstra.o $(PATHTEMP)/Preprocess-Dijkstra.o $(PATHTEMP)/EVCSLP-KKT.o $(PATHTEMP)/EVCSLP.o $(PATHTEMP)/EVCSLP-SPL.o $(PATHTEMP)/main.o
	$(CPP) $(CCFLAGS) $(PATHTEMP)/InstanceRead.o $(PATHTEMP)/BigM.o $(PATHTEMP)/Dijkstra.o $(PATHTEMP)/Preprocess-Dijkstra.o $(PATHTEMP)/EVCSLP-KKT.o $(PATHTEMP)/EVCSLP.o $(PATHTEMP)/EVCSLP-SPL.o $(PATHTEMP)/main.o $(CCLNFLAGS) -o $(PATHEXEC)/$(EXECUTAVEL)

$(PATHTEMP)/InstanceRead.o: $(PATHSRC)/InstanceRead.cpp
	$(CPP) $(CCFLAGS) -c $(PATHSRC)/InstanceRead.cpp -o $(PATHTEMP)/InstanceRead.o
$(PATHTEMP)/BigM.o: $(PATHSRC)/BigM.cpp
	$(CPP) $(CCFLAGS) -c $(PATHSRC)/BigM.cpp -o $(PATHTEMP)/BigM.o
$(PATHTEMP)/Dijkstra.o: $(PATHSRC)/Dijkstra.cpp
	$(CPP) $(CCFLAGS) -c $(PATHSRC)/Dijkstra.cpp -o $(PATHTEMP)/Dijkstra.o
$(PATHTEMP)/Preprocess-Dijkstra.o: $(PATHSRC)/Preprocess-Dijkstra.cpp
	$(CPP) $(CCFLAGS) -c $(PATHSRC)/Preprocess-Dijkstra.cpp -o $(PATHTEMP)/Preprocess-Dijkstra.o
$(PATHTEMP)/EVCSLP-KKT.o: $(PATHSRC)/EVCSLP-KKT.cpp
	$(CPP) $(CCFLAGS) -c $(PATHSRC)/EVCSLP-KKT.cpp -o $(PATHTEMP)/EVCSLP-KKT.o
$(PATHTEMP)/EVCSLP.o: $(PATHSRC)/EVCSLP.cpp
	$(CPP) $(CCFLAGS) -c $(PATHSRC)/EVCSLP.cpp -o $(PATHTEMP)/EVCSLP.o
$(PATHTEMP)/EVCSLP-SPL.o: $(PATHSRC)/EVCSLP-SPL.cpp
	$(CPP) $(CCFLAGS) -c $(PATHSRC)/EVCSLP-SPL.cpp -o $(PATHTEMP)/EVCSLP-SPL.o
$(PATHTEMP)/main.o: $(PATHSRC)/main.cpp
	$(CPP) $(CCFLAGS) -c $(PATHSRC)/main.cpp -o $(PATHTEMP)/main.o


####################################################
###################### CLEAN #######################
####################################################
clean:
	rm -rf $(PATHEXEC)
	rm -rf $(PATHTEMP)

####################################################
####################### CPLEX ######################
####################################################

##### CPLEX CONFIGURATION's
# System architecture
SYSTEM = x86-64_linux
# Static library format for Cplex
LIBFORMAT = static_pic
# Cplex directory
CPLEXDIR = /opt/ibm/ILOG/CPLEX_Studio201/cplex
# Concert directory
CONCERTDIR = /opt/ibm/ILOG/CPLEX_Studio201/concert

##### CPLEX DIRECTIVE's
# Cplex static libraries directory
CPLEXLIBDIR = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
# Concert static libraries directory
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
# Cplex header's directory
CPLEXINCDIR = $(CPLEXDIR)/include
# Concert header's directory
CONCERTINCDIR = $(CONCERTDIR)/include

####################################################
##################### COMPILER #####################
####################################################

##### COMPILER CONFIGURATION's
# Compiler
CPP = g++
# Compilation parameters
CCOPT = -std=c++17 -O3 -g -fPIC -fexceptions -DIL_STD 
## Include libraries identifiers
CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert  -m64 -lm -pthread -ldl
# Header's include path
CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR) 
