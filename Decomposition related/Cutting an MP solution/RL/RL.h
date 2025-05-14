#ifndef RL_H
#define RL_H

#include "model.h"
void RL(Instance& inst, Solution& sol, int& idx, int& idx2, int& a);
void createStory(int& idx, int& idx2, int& a, double& time_Val);
void train(Instance& inst, Solution& sol, Agent& agent);
#endif 

