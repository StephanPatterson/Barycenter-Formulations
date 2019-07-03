//
//  main.cpp
//  DigitsDataLPOriginal: This will calculate the exact barycenter for a set of digits data using Grid Knowledge
//  Now called Original Formulation with Knowledge - using set G
//
//  This should be the fastest grid version that uses LP Original due to the short setup time.
//  Created by Stephan Patterson on 11/22/17.
//  Copyright © 2017 Stephan Patterson. All rights reserved.


#include "SuppPt.hpp"
#include "/Library/gurobi801/mac64/include/gurobi_c++.h"
#include <iostream>
#include <cmath>
#include <sstream>
#include <chrono>
#include <fstream>

bool validateDist(const std::vector<SuppPt>::iterator, const std::vector<SuppPt>::iterator);

int main(int argc, const char * argv[]) {
   
   auto start = std::chrono::steady_clock::now();
   auto t = start;
   
   //std::cout.precision(17);
   
   //which digit to be analyzed
   int digit = 0;
   
   //opens file containing data for the chosen digit
   std::ifstream indata;
   std::ostringstream fname;
   fname << "/Users/spatterson/Documents/Barycenter/digit" << digit << "p.txt";
   indata.open(fname.str());
   
   // How many digits to include.
   const int Pnum = 3;
   std::cout << "Exact using Digit: " << digit << " Pnum: " << Pnum << std::endl;
   double lambda[Pnum];
   int Psnum[Pnum];
   for (int i = 0; i < Pnum;++i)
   {
      lambda[i]= 1.0/Pnum;
      Psnum[i] = 0;
   }
   
   //If xnum/ynum not both powers of 2, there will be numerical error with using steps as coded.
   const int xnum = 16;
   const int ynum = 16;
   const int gridsize = xnum*ynum;
   double tempx;
   double tempy;
   double tempmass;
   double xstep = 1.0/xnum;
   double ystep = 1.0/ynum;
   int index = 0;
   double minx = 1.0;
   double miny = 1.0;
   double maxx = 0.0;
   double maxy = 0.0;
   std::vector<SuppPt> Psupp(gridsize*Pnum);
   //These coordinates begin at the bottom left corner
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
               if (tempx < minx)
               {
                  minx = tempx;
               }
               if (tempy < miny)
               {
                  miny = tempy;
               }
               if (tempx > maxx)
               {
                  maxx = tempx;
               }
               if (tempy > maxy)
               {
                  maxy = tempy;
               }
            }
            tempx += xstep;
         }
         tempy += ystep;
      }
   }
   
   indata.close();
   
   auto dataendtime = std::chrono::steady_clock::now();
   std::cout << "Total Data Read-in Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(dataendtime-start).count() << "ms" << std::endl;
   
   for (int i = 0; i < Pnum; ++i)
      std::cout << "Size of S_" << i <<": " <<Psnum[i] << std::endl;
   std::cout << "Grid generated between: " << minx << " " << maxx << " " << miny << " " << maxy <<std::endl;
   int gridxdlength = (maxx-minx)/xstep;
   int gridydlength = (maxy-miny)/ystep;
   int gridxlength = Pnum*(gridxdlength);
   int gridylength = Pnum*(gridydlength);
   ++gridxdlength;
   ++gridydlength;
   ++gridxlength;
   ++gridylength;
   
   const int totalsupp = index;
   Psupp.resize(totalsupp);
   
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
   
   //This saves the starting index for each measure. Should not be modified after setting.
   int startindices[Pnum];
   
   startindices[0] = 0;
   for (int i = 1; i < Pnum; ++i)
   {
      startindices[i] = startindices[i-1]+Psnum[i-1];
   }
   int index2 = 0;
   
   std::vector<SuppPt> Pbar0(S0);
   double xgstep = xstep/Pnum;
   double ygstep = ystep/Pnum;
   
   //This is the old version of generating the grid steps. New version better handles Pnum that is a non-multiple of 2
   //   double xgstep[Pnum];
   //   double ygstep[Pnum];
   //   double xtotal = 0;
   //   double ytotal = 0;
   //   for (int i = 0; i < Pnum; ++i)
   //   {
   //      xgstep[i] = xstep-xtotal;
   //      xgstep[i] /= (Pnum-i);
   //      xtotal += xgstep[i];
   //      ygstep[i] = ystep-ytotal;
   //      ygstep[i] /= (Pnum-i);
   //      ytotal += ygstep[i];
   //   }
   //   std::cout << xstep << " " << xtotal << std::endl;
   ////   xgstep[Pnum-1] = xstep-xtotal;
   ////   ygstep[Pnum-1] = xstep-ytotal;
   //   for(int i = 0; i < Pnum; ++i)
   //   {
   //      std::cout << xgstep[i] << " " << ygstep[i] << std::endl;
   //   }
   //   xtotal  = 0;
   //   ytotal = 0;
   //   for (int i = 0; i < Pnum; ++i)
   //   {
   //      xtotal += xgstep[i];
   //      ytotal += ygstep[i];
   //   }
   //   std::cout << xtotal << " " <<ytotal << std::endl;
   tempy = miny;
   
   for (int j = 0; j < gridylength; ++j)
   {
      tempx = minx;
      int ygindex = j % Pnum;
      if (tempy >= miny && tempy <= maxy)
      {
         for (int k = 0; k < gridxlength; ++k)
         {
            int xgindex = k % Pnum;
            if (tempx >= minx && tempx <= maxx)
            {
               Pbar0[index2] = SuppPt(tempx, tempy, 0.0);
               ++index2;
            }
            tempx += xgstep;
            if (xgindex == Pnum-1)
            {
               tempx = minx + (k+1)/Pnum*xstep;
            }
         }
      }
      tempy += ygstep;
      if (ygindex == Pnum-1)
      {
         tempy = miny + (j+1)/Pnum*ystep;
      }
   }
   
   S0 = index2;
   std::cout<< "Size of no-duplicate possible support " << S0 << std::endl;
   Pbar0.resize(S0);
   
   t = std::chrono::steady_clock::now();
   std::cout << "Total setup time, prior to model instantiation: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-dataendtime).count() << "ms" << std::endl;
   auto told = t;
   
   GRBEnv* env = new GRBEnv();
   GRBModel* model = new GRBModel(*env);
   model->set(GRB_IntParam_Method, 0);
   model->set(GRB_IntParam_Presolve, 0);
   
   GRBVar* z;
   //If S0 is too big this needs to be modified to add in batches.
   z = model->addVars(S0);
   
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
      //         model->update();
      //         std::cout << y[0][0].get(GRB_StringAttr_VarName) << std::endl;
   }
   //      model.setObjective(obj);
   model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
   //   model->update();
   //   std::ostringstream filename;
   //   filename << "/Users/spatterson/Documents/test.lp";
   //   model->write(filename.str());

   t = std::chrono::steady_clock::now();
   std::cout << "Total model setup time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count() << "ms" << std::endl;
   
   model->optimize();
   
   t = std::chrono::steady_clock::now();
   std::cout << "Total model solution time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count() << "ms" << std::endl;
   
   //Read out the solution. This is then written to file; display code is separate.
   std::vector<SuppPt> Pbar01(S0);
   int S01=0;
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
   std::cout << "Total time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-start).count() << "ms" << std::endl;

   return 0;
}

//Checks that the mass sums to 1 for the given points
//Possible Future Modification: Add check for duplicate support points?
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

