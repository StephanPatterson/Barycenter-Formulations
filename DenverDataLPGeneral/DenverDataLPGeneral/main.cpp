//
//  main.cpp
//  DenverGen
//
//  Created by Stephan Patterson on 10/7/18.
//  Copyright © 2018 Stephan Patterson. All rights reserved.
//

#include "/Library/gurobi801/mac64/include/gurobi_c++.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include "SuppPt.hpp"

bool validateDist(const std::vector<SuppPt>::iterator, const std::vector<SuppPt>::iterator);

int main(int argc, const char * argv[]) {
   
   auto start = std::chrono::steady_clock::now();
   auto t = start;
   
   int startyear = 12;
   std::string startmonth = "Feb";//Jan 2012 not in data set. Start: Feb 2012.
   int Pnum =11;//If NumMoMeas = 1, Number of months to process. Will go consecutively from start month, skipping any empty months.
   //Otherwise, Number of measures that will be created. Will combine NumMoMeas of months into each measure
   int NumMoMeas = 1;//Number of Months in a single Measure. Use 1 for original approach of each month is its own measure
   
   std::vector<SuppPt> Psupp;
   int totalsupp = 0;
   int startindices[Pnum];
   int indices[Pnum];
   int endindices[Pnum];
   int Psnum[Pnum];
   double lambda[Pnum];
   int dotoy = 0; //If 0, runs denver data. if 1, runs toy example
   if (dotoy ==1)
   {
      if (Pnum >3)
      {
         std::cout << "Adjust hard coded sizes." << std::endl;
         return 0;
      }
      Psnum[0] = 2;
      Psnum[1] = 2;
      Psnum[2] = 3;
      
      double loc1 = 0.0;
      double loc2 = 0.0;
      for (int i = 0; i < Pnum; ++i)
      {
         loc2 = 0.0;
         startindices[i] = totalsupp;
         indices[i] = startindices[i];
         if (i > 0)
         {
            endindices[i-1] = startindices[i]-1;
         }
         for (int j = 0; j < Psnum[i]; ++j)
         {
            Psupp.push_back(SuppPt(Pnum*loc1, Pnum*loc2 , 1.0/Psnum[i], totalsupp));
            ++totalsupp;
            ++loc2;
         }
         ++loc1;
      }
      endindices[Pnum-1] = totalsupp;
      --endindices[Pnum-1];
      
      // Hard code a swap of P_2 (index 1) so that we begin non-optimal
      std::swap(Psupp[startindices[1]], Psupp[startindices[1]+1]);
      --Psupp[startindices[1]].combos;
      ++Psupp[startindices[1]+1].combos;
      double sum = 0;
      for (int i = 0; i < Pnum-1; ++i)
      {
         //      lambda[i] = (double)Psnum[i]/totalsupp;
         lambda[i] = 1.0/Pnum;
         sum += lambda[i];
      }
      lambda[Pnum-1] = 1-sum;
   }
   else
   {
      std::string temp;
      
      std::ifstream indata;
      std::ostringstream fname;
      int year;
      std::string month;
      fname << "/Users/spatterson/Documents/Column Generation/murdercleaned.csv";
      indata.open(fname.str());
      getline(indata, temp, '\n');// Advance past headers
      indata >> year;
      if (startyear < year)
      {
         std::cout << "Invalid Year Selection." << std::endl;
         return 1;
      }
      while (year < startyear)
      {
         getline(indata, temp, '\n');// Advance to next entry
         indata >> year;
      }
      indata >> month;
      while (month != startmonth)
      {
         getline(indata, temp, '\n');
         indata >> year;
         indata >> month;
         if (year != startyear)
         {
            std::cout << "Selected Month/Year not in data." << std::endl;
            return 1;
         }
         if (indata.eof())
         {
            indata.close();
            std::cout << "Invalid Input. End of file reached." << std::endl;
            return 1;
         }
      }
      double loc1;
      double loc2;
      indata >> loc1;
      indata >> loc2;
      std::string currentmonth = startmonth;
      
      for (int i = 0; i < Pnum; ++i)
      {
         Psnum[i] = 0;
         startindices[i] = totalsupp;
         indices[i] = totalsupp;
         for (int j = 0; j < NumMoMeas; ++j)
         {
            while (month == currentmonth )
            {
               Psupp.push_back(SuppPt(loc1, loc2 , 1.0));
               ++Psnum[i];
               ++totalsupp;
               
               indata >> year;
               indata >> month;
               indata >> loc1;
               indata >> loc2;
               if (indata.eof())
               {
                  indata.close();
                  std::cout << "Invalid Input. End of file reached." << std::endl;
                  return 1;
               }
            }
            currentmonth = month;
         }
         double totalmass = 0;
         for (int j = 0; j < Psnum[i]-1; ++j)
         {
            Psupp[startindices[i]+j].mass /= Psnum[i];
            totalmass += Psupp[startindices[i]+j].mass;
         }
         Psupp[totalsupp-1].mass = 1-totalmass;
         currentmonth = month;
         endindices[i] = totalsupp-1;
      }
      Psupp.resize(totalsupp);
      indata.close();
      
      double sum = 0;
      for (int i = 0; i < Pnum-1; ++i)
      {
         lambda[i] = (double)Psnum[i]/totalsupp;
         sum += lambda[i];
      }
      lambda[Pnum-1] = 1-sum;
   }
   
   int index = 0;
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
   
   std::vector<SuppPt> Pbar0(S0);
   std::vector<SuppPt>::iterator it = Pbar0.begin();
   
   GRBEnv* env = new GRBEnv();
   GRBModel* model = new GRBModel(*env);
   model->set(GRB_IntParam_Method, 0);
   model->set(GRB_IntParam_Presolve, 0);
   model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
   
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
         sum1 += lambda[i]*Psupp[index].loc1;
         sum2 += lambda[i]*Psupp[index].loc2;
      }
      
      Pbar0[j] = SuppPt(sum1, sum2, 0.0);
      for (int i = 0; i < Pnum; ++i)
      {
         index = indices[i];
         cost +=lambda[i]*((Pbar0[j].loc1-Psupp[index].loc1)*(Pbar0[j].loc1-Psupp[index].loc1)+(Pbar0[j].loc2-Psupp[index].loc2)*(Pbar0[j].loc2-Psupp[index].loc2));
         //         std::cout << index << " " ;
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
   t = std::chrono::steady_clock::now();
   std::cout << "Setup Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-start).count()  <<"ms" << std::endl;
   /*
    t = clock();
    std::cout << "Setup time: " << (double)(t-t0)/CLOCKS_PER_SEC << std::endl;
    */
   model->optimize();
   auto t0 = t;
   t = std::chrono::steady_clock::now();
   std::cout << "Solve Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-t0).count()  <<"ms" << std::endl;
   
   std::vector<SuppPt> Pbar01(S0);
   int S01=0;
   double tempmass;
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
      filename2 << "/Users/spatterson/Documents/Column Generation/GenObj" << startyear << "_" << Pnum << ".txt";
      std::ofstream outresults;
      outresults.open(filename2.str());
      
      outresults << std::setprecision(10) << model->get(GRB_DoubleAttr_ObjVal) << '\n';
      for (it = Pbar01.begin(); it != Pbar01.end(); ++it)
      {
         outresults << *it << '\n';
      }
   }
   else
   {
      std::cout << "Optimal Solution not reached." << std::endl;
   }
   
   model->terminate();
   model->reset();
   delete model;
   delete env;
   t = std::chrono::steady_clock::now();
   std::cout << "Total Running time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-start).count()  <<"ms" << std::endl;
   /*
    t = clock()-t0;
    std::cout << "Total run time: " << (double)t/CLOCKS_PER_SEC << std::endl;
    */
   return 0;
}

bool validateDist(const std::vector<SuppPt>::iterator it1, const std::vector<SuppPt>::iterator it2)
{
   std::vector<SuppPt>::iterator it = it1;
   double total = 0;
   while (it != it2)
   {
      total += it->mass;
      ++it;
   }
   //   std::cout << total << std::endl;
   
   if (total < 0.95 || total > 1.05)
   {
      std::cout << "Warning: Far from 1. Consider alternative." << total << std::endl;
   }
   //   std::cout << total << std::endl;
   if (total != 1)
   {
      it1->mass = it1->mass + (1.0-total);
   }
   return true;
}

