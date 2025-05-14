This repository contains the code for all algorithms discussed in the paper "Solving the strip packing problem with a decomposition framework and a generic solver: implementation, tuning, and reinforcement-learning-based hybridization" by Fatih Burak Akcay and Maxence Delorme. 

Our algorithms are coded in C++ and, unless stated otherwise, use the commercial solver Gurobi for the ILP models and IBM CPLEX CP SOLVER for the CP models. 
The code is divided over 2 main folders (one for the direct models and one for the decomposition related models), each containing the subfolders for each method discussed.

'Direct Models' main folder:
- CC_CP					| Model CC-CP
- CC_ILP				| Model CC-ILP
- D_ILP					| Model D-ILP
- DC_CP 				| Model DC-CP
- DC_ILP				| Model DC-ILP

'Decomposition folder' has three subfolders: 'Secondary Problem', 'Rejecting Solutions' and 'Cutting an MP solution'. 
For 'Secondary Problem' subfolder	|Includes a full decomposition model where MP is formulated as a CP model and secondary problem is formulated using the model SP-ILP and SP-CP in 					|the slave_functions.cpp
For 'Rejecting Solutions' subfolder:
- CPLEX_ALTERNATE			| Model MP-ILP+DB+ALT (solver: IBM CPLEX) 
- CPLEX_ILP_MIN_REJ			| Model MP-ILP+REJ (solver: IBM CPLEX)
- CPLEX_ILP_MIN_REJ_UB			| Model MP-ILP+REJ with UB (solver: IBM CPLEX)
- CPLEX_REJECT				| Model MP-ILP+DB+REJ (solver: IBM CPLEX)
- CP_ALT				| Model MP-CP+ALT 
- ILP_BIN_ALT				| Model MP-ILP+ALT 
- ILP_INT_ALT				| Model MP-ILP+ALT+INT
For Cutting an MP solution subfolder:
- MP_CP_NG				| Model MP-CP+DEL-NGC
- MP_ILP_NG				| Model MP-ILP+DEL-NGC
- CP_DEL_LCBC				| Model MP-CP+DEL-LCBC
- ILP_DEL_LCBC				| Model MP-ILP+DEL-LCBC
- ILP_OTF				| Model MP-ILP+OTF
- ILP_FAKE_OTF				| Model MP-ILP+DUMCB
- ILP_OBJ_DIR_A				| Model MP-ILP+A+DEL
- ILP_OBJ_DIR_W				| Model MP-ILP+W+DEL
- ILP_OBJ_DIR_H				| Model MP-ILP+H+DEL
- RL					| REINFORCEMENT LEARNING MODEL with FT-ILP + FT-CP

Each subfolder contains the same substructure:
- helper_functions.cpp		| Contains a number of secondary functions (this file is usually the same for each subfolder)
- helper_functions.h		| The header file corresponding to helper_functions.cpp (this file is usually the same for each subfolder)	| 
- main.cpp			| The front-end code for using the method  
- main.h			| The header file corresponding to main.cpp 
- LB.cpp			| Contains the lower bound models
- LB.h				| The header file corresponding to LB.cpp
- model.cpp			| Contains the model of the tested method
- model.h			| The header file corresponding to model.cpp 
- makefile			| Used for compiling under linux (it needs to be written by the user)
Note: for decomposition related models, slave_functions.cpp and slave_functions.h are included, containing secondary problem related functions.

********
Once compiled, the following command can be used to run the algorithm:
	./PROGRAM "./PATH_INSTANCE" "NAME_INSTANCE" "./PATH_AND_NAME_OUTPUT_GENERAL" 
where
- PROGRAM is the name of the compiled software 
- ./PATH_INSTANCE is the relative path of the folder where the instance to solve is located
- NAME_INSTANCE is the name of the instance to solve
- ./PATH_AND_NAME_OUTPUT_GENERAL is the name of the file (together with its relative path) where performance metrics (such as the optimality status, the CPU time required, or the number of variables) are stored after solving an instance
********

Moreover, "_INPUT.rar" contains a txt-file for each of our test instances. 
Each txt-file is structured as follows:
- the first line contains the number of items (n)
- the second line contains the width of the strip (W)
- the remaining (n) lines all contain, for each item:
    	- the item dimensions (width and height) and its demand (where demand is set to 1) (w_j h_j 1) 
