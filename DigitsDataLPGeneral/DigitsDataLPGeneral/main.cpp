//
//  main.cpp
// MNIST Digits Data
// This is the all w formulation. Not a good choice. This version is for comparison purposes.
//
//  Created by Stephan Patterson.
//  Copyright © 2018 Stephan Patterson. All rights reserved.

#include "SuppPt.hpp"
#include "/Library/gurobi801/mac64/include/gurobi_c++.h"
#include <iostream>
#include <cmath>
#include <sstream>
#include <ctime>
#include <fstream>

bool validateDist(const std::vector<SuppPt>::iterator, const std::vector<SuppPt>::iterator);

int main(int argc, const char * argv[]) {
   
   int digit = 1;
   const int Pnum = 4;
   std::cout << "Exact using Digit: " << digit << " Pnum: " << Pnum << std::endl;
   
   auto start = std::chrono::steady_clock::now();
   auto t = start;
   
   std::ifstream indata;
   std::ostringstream fname;
   fname << "/Users/spatterson/Documents/Barycenter/digit" << digit << "p.txt";
   indata.open(fname.str());
   double lambda[Pnum];
   int Psnum[Pnum];
   for (int i = 0; i < Pnum;++i)
   {
      lambda[i]= 1.0/Pnum;
      Psnum[i] = 0;
   }
   
   const int xnum = 16;
   const int ynum = 16;
   const int gridsize = xnum*ynum;
   double tempx;
   double tempy;
   double tempmass;
   double xstep = 1.0/xnum;
   double ystep = 1.0/ynum;
   int index = 0;
   std::vector<SuppPt> Psupp(gridsize*Pnum);
   
   int swwd = 0;// Which Digit to start with. Default of 0 starts from first available.
   for (int i = 0; i < swwd; ++i)
   {
      for (int j = 0; j < xnum; ++j)
      {
         for (int k = 0; k < ynum; ++k)
         {
            indata >> tempmass;
         }
      }
   }
   for (int i = 0; i < Pnum; ++i)
   {
      tempx = 0;
      tempy = 0;
      for (int j = 0; j < xnum; ++j)
      {
         tempx = 0;
         for (int k = 0; k < ynum; ++k)
         {
            indata >> tempmass;
            if (tempmass > 0)
            {
               Psupp[index] = SuppPt(tempx, tempy, tempmass);
               ++index;
               ++Psnum[i];
            }
            tempx += xstep;
         }
         tempy += ystep;
      }
   }
   
   indata.close();
   for (int i = 0; i < Pnum; ++i)
      std::cout << "Size of S_" << i << ": " <<Psnum[i] << std::endl;
   
   // remove extra support points from Psupp
   const int totalsupp = index;
   Psupp.resize(totalsupp);
   // ensure P_i sum to exactly 1
   std::vector<SuppPt>::iterator Psupptest = Psupp.begin();
   for (int i = 0; i < Pnum; ++i)
   {
      validateDist(Psupptest, Psupptest+Psnum[i]);
      Psupptest += Psnum[i];
   }
   
   //Calculate number of possible supp pts S0
   long int S0 = 1;
   for (int i = 0; i < Pnum; ++i)
   {
      S0 *= Psnum[i];
   }
   std::cout << "Size of S0: " << S0 << std::endl;
   
   int startindices[Pnum];
   int endindices[Pnum];
   int indices[Pnum];
   
   startindices[0] = 0;
   for (int i = 0; i < Pnum; ++i)
   {
      if (i > 0)
      {
         startindices[i] = endindices[i-1]+1;
      }
      endindices[i] = startindices[i]+Psnum[i]-1;
      indices[i] = startindices[i];
   }
   
   
   std::vector<SuppPt> Pbar0(S0);
   std::vector<SuppPt>::iterator it = Pbar0.begin();
   
   GRBEnv* env = new GRBEnv();
   GRBModel* model = new GRBModel(*env);
   model->set(GRB_IntParam_Method, 0);
   
   std::vector< GRBVar > w(S0);
   GRBLinExpr  exp[totalsupp];
   
   index = 0;
   for (unsigned long int j = 0; j < S0; ++j)
   {
      double sum1 = 0;
      double sum2 = 0;
      double cost = 0;
      for (int i = 0; i < Pnum; ++i)
      {
         index = indices[i];
         sum1 += Psupp[index].loc1;
         sum2 += Psupp[index].loc2;
      }
      sum1 /= Pnum;
      sum2 /= Pnum;
      
      Pbar0[j] = SuppPt(sum1, sum2, 0.0);
      for (int i = 0; i < Pnum; ++i)
      {
         index = indices[i];
         cost +=lambda[i]*((Pbar0[j].loc1-Psupp[index].loc1)*(Pbar0[j].loc1-Psupp[index].loc1)+(Pbar0[j].loc2-Psupp[index].loc2)*(Pbar0[j].loc2-Psupp[index].loc2));
      }
      w[j] = model->addVar(0.0, 1.0, cost, GRB_CONTINUOUS);//This updates the objective
      for (int i = 0; i < Pnum; ++i)
      {
         index = indices[i];
         exp[index] += w[j];
      }
      
      int k = Pnum-1;
      if (indices[k] < endindices[k])
      {
         ++indices[k];
      }
      else
      {
         int temp = k-1;
         while (temp >= 0)
         {
            if (indices[temp] == endindices[temp])
            {
               temp -= 1;
            }
            else
            {
               ++indices[temp];
               break;
            }
         }
         for (int l = k; l > temp; --l)
         {
            indices[l] = startindices[l];
         }
      }
   }
   for (int j = 0; j < totalsupp; ++j)
   {
      model->addConstr(exp[j] == Psupp[j].mass);
   }
   
   model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
   model->set(GRB_IntParam_Presolve, 0);
   
   t = std::chrono::steady_clock::now();
   std::cout << "Setup time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-start).count() << "ms" << std::endl;
   
   model->optimize();
   
   auto told = t;
   t = std::chrono::steady_clock::now();
   std::cout << "Model Solve time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count() << "ms" << std::endl;
   
   std::vector<SuppPt> Pbar01(S0);
   int S01=0;
   if (model->get(GRB_IntAttr_Status) == 2)
   {
      it = Pbar0.begin();
      for (int j = 0; j < S0; ++j)
      {
         tempmass = w[j].get(GRB_DoubleAttr_X);
         if (tempmass != 0)
         {
            Pbar01[S01] = SuppPt(Pbar0[j].loc1, Pbar0[j].loc2, tempmass);
            ++S01;
         }
      }
      
      Pbar01.resize(S01);
      
      std::ostringstream filename2;
      filename2 << "/Users/spatterson/Documents/testouttrue" << digit << "_" << Pnum << ".txt";
      std::ofstream outresults;
      outresults.open(filename2.str());
      for (it = Pbar01.begin(); it != Pbar01.end(); ++it)
      {
         outresults << *it << '\n';
      }
      told = t;
      t = std::chrono::steady_clock::now();
   }
   else
   {
      std::cout << "Optimal Solution not reached." << std::endl;
   }
   
   model->terminate();
   model->reset();
   delete model;
   delete env;
   
   told = t;
   t = std::chrono::steady_clock::now();
   std::cout << "Total run time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count() << "ms"<< std::endl;
   
   return 0;
}

//Checks that the mass sums to 1 for the given points
//Add check for duplicate support points?
bool validateDist(const std::vector<SuppPt>::iterator it1, const std::vector<SuppPt>::iterator it2)
{
   std::vector<SuppPt>::iterator it = it1;
   double total = 0;
   while (it != it2)
   {
      total += it->mass;
      ++it;
   }
   
   if (total < 0.95 || total > 1.05)
   {
      std::cout << "Warning: Far from 1. Consider alternative." << total << std::endl;
   }
   if (total != 1)
   {
      it1->mass = it1->mass + (1.0-total);
   }
   return true;
}
