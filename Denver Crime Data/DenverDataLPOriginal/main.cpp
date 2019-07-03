//
//  main.cpp
//  Denver Data with LP Original
//
//  Created by Stephan Patterson on 1/7/19.
//  Copyright © 2019 Stephan Patterson. All rights reserved.
//

#include "/Library/gurobi801/mac64/include/gurobi_c++.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include "SuppPt.hpp"

bool validateDist(const std::vector<SuppPt>::iterator, const std::vector<SuppPt>::iterator);

int main(int argc, const char * argv[]) {
   
   auto start = std::chrono::steady_clock::now();
   auto t = start;
   
   int startyear = 14;
   std::string startmonth = "Jan";
   int Pnum =12;//If NumMoMeas = 1, Number of months to process. Will go consecutively from start month, skipping any empty months.
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
            //            std::cout << Psupp[totalsupp] << std::endl;
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
      fname << "/DenverCrime.csv";
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
      while (month != startmonth)//Advance to start month
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
   
   //Create possible support set S.
   int index2 = 0;
   std::vector<SuppPt> Pbar0(S0);
   index = 0;
   for (int j = 0; j < S0; ++j)
   {
      double sum1 = 0;
      double sum2 = 0;
      for (int i = 0; i < Pnum; ++i)
      {
         index = indices[i];
         sum1 += lambda[i]*Psupp[index].loc1;
         sum2 += lambda[i]*Psupp[index].loc2;
      }
      //Removing the check for uniqueness for the Denver data, as it is expected to be in general position with no duplicate x_j
      /*      bool unique = 0;
       for (int l = 0; l < index2; ++l)
       {
       if (Pbar0[l].loc1 == sum1)
       {
       if (Pbar0[l].loc2 == sum2)
       {
       unique = 1;
       break;
       }
       }
       }
       
       if (unique == 0)
       {*/
      Pbar0[index2] = SuppPt(sum1, sum2, 0.0);
      ++index2;
      //      }
      
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
   
   S0 = index2;
   std::cout<< "Size of no-duplicate possible support " << S0 << std::endl;
   
   GRBEnv* env = new GRBEnv();
   GRBModel* model = new GRBModel(*env);
   model->set(GRB_IntParam_Method, 2);
   model->set(GRB_IntParam_Presolve, 0);
   
   GRBVar* z;
   //If S0 is too large, this will need to be adjusted to add variables in batches.
   z = model->addVars(S0);
   
   // Optional code for naming the variables. Nice if printing the LP.
   //   for (int j = 0; j < S0; ++j)
   //   {
   //      std::ostringstream vname;
   //      vname << "w_" << j+1;
   //      w[j].set(GRB_StringAttr_VarName, vname.str());
   //         w[j].set(GRB_DoubleAttr_UB, 1.0);
   //   }
   
   std::vector<SuppPt>::iterator it = Pbar0.begin();
   int Psmax = *std::max_element(Psnum, Psnum+Pnum);
   std::vector< std::vector <GRBVar> > y(S0,std::vector<GRBVar> (Psmax));
   std::vector< std::vector<double> > cost(S0, std::vector<double>(Psmax));
   
   for (int i = 0; i < Pnum; ++i)
   {
      GRBLinExpr exp2[Psnum[i]];
      // Calculate distances
      index2 = startindices[i];
      for (int j = 0; j < S0; ++j)
      {
         for (int k = 0; k < Psnum[i]; ++k)
         {
            cost[j][k] = (Pbar0[j].loc1-Psupp[k+index2].loc1)*(Pbar0[j].loc1-Psupp[k+index2].loc1)+(Pbar0[j].loc2-Psupp[k+index2].loc2)*(Pbar0[j].loc2-Psupp[k+index2].loc2);
         }
      }
      for (int j = 0; j < S0; ++j)
      {
         GRBLinExpr exp;
         for (int k = 0; k < Psnum[i]; ++k)
         {
            //            std::ostringstream vname;
            //            vname << "y" << i+1 << "_" << j+1  << "_" << k+1;
            y[j][k] = model->addVar(0.0, 1.0, lambda[i]*cost[j][k], GRB_CONTINUOUS);
            //            y[j][k].set(GRB_StringAttr_VarName, vname.str());
            exp += y[j][k];
            exp2[k] += y[j][k];
         }
         model->addConstr(exp == z[j]);
      }
      for (int k = 0; k < Psnum[i]; ++k)
      {
         model->addConstr(exp2[k] == Psupp[index2+k].mass);
      }
   }
   model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
   
   //      model->update();
   //   std::ostringstream filename;
   //   filename << "/Users/spatterson/Documents/test.lp";
   //   model->write(filename.str());
   t = std::chrono::steady_clock::now();
   std::cout << "Setup Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-start).count()  <<"ms" << std::endl;
   
   model->optimize();
   
   auto t0 = t;
   t = std::chrono::steady_clock::now();
   std::cout << "Solve Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-t0).count()  <<"ms" << std::endl;
   
   //This reads out the solution, pulling only those support points from Pbar0 (possible support) and assigning them to Pbar01 with appropriate mass.
   //Depending on what is desired with the barycenter, other forms of output and storage may be preferred.
   std::vector<SuppPt> Pbar01(S0);
   int S01=0;
   double tempmass;
   if (model->get(GRB_IntAttr_Status) == 2)
   {
      it = Pbar0.begin();
      for (int j = 0; j < S0; ++j)
      {
         tempmass = z[j].get(GRB_DoubleAttr_X);
         if (tempmass != 0)
         {
            Pbar01[S01] = SuppPt(Pbar0[j].loc1, Pbar0[j].loc2, tempmass);
            ++S01;
         }
         //   std::cout << w[j].get(GRB_DoubleAttr_X) << std::endl;
      }
      
      Pbar01.resize(S01);
   }
   else
   {
      std::cout << "Optimal Solution not reached." << std::endl;
   }
   
   delete [] z;
   model->terminate();
   model->reset();
   delete model;
   delete env;
   
   t = std::chrono::steady_clock::now();
   std::cout << "Total Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-start).count()  <<"ms" << std::endl;
   
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
