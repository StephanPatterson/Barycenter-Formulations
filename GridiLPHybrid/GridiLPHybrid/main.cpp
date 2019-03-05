//
//  main.cpp
//  Hybrid using "i"
//  Dots are sample, hard-coded generation.
//
//  Created by Stephan Patterson on 3/12/18.
//  Copyright © 2018 Stephan Patterson. All rights reserved.
//

#include "SuppPt.hpp"
#include "/Library/gurobi801/mac64/include/gurobi_c++.h"
#include <iostream>
#include <fstream>
#include <cmath>
#include <sstream>
#include <ctime>
#include <list>
#include <vector>
#include <random>

bool validateDist(const std::vector<SuppPt>::iterator, const std::vector<SuppPt>::iterator);

template< int n>
struct costcoord
{
   double ch;
   unsigned long int ind[n];
   //   int xjindex;
};

int main(int argc, const char * argv[]) {
   
   clock_t t = clock();
   clock_t t0 = t;
   
   const int Pnum = 4;
   int digit = 1;
   int finer = 1; // Grid step for i-dot
   
   int Psnum[Pnum];
   double lambda[Pnum];
   for (int i = 0; i < Pnum;++i)
   {
      lambda[i]= 1.0/Pnum;
      Psnum[i]=0;
   }
   
   std::ifstream indata;
   std::ostringstream fname;
   fname << "/Users/spatterson/Documents/Barycenter/digit" << digit << "p.txt";
   indata.open(fname.str());
   const int xnum = 16;
   const int ynum = 16;
   const int gridsize = xnum*ynum;
   double tempx;
   double tempy;
   double tempmass;
   double xstep = 1.0/xnum;
   double ystep = 1.0/ynum;
   unsigned long int index = 0;
   double minx = 1.0;
   double miny = 1.0;
   double maxx = 0.0;
   double maxy = 0.0;
   double maximass = 0;
   unsigned long int startindices[Pnum];
   unsigned long int indices[Pnum];
   unsigned long int endindices[Pnum];
   std::vector<SuppPt> Psupp(gridsize*Pnum+8*Pnum);

   
   //These coordinates read in assuming bottom left corner of each pixel
   for (int i = 0; i < Pnum; ++i)
   {
      startindices[i] = index;
      indices[i] = index;
      tempx = 0;
      tempy = 0;
      maximass = 0;
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
               if (tempmass > maximass)
               {
                  maximass = tempmass;
               }
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
      
      double xstep2 = xstep/finer;
      double ystep2 = ystep/finer;
      //Finer Grid code
      if (i == 0)
      {
         tempy = 29*ystep;
         tempx = 6*xstep;
         Psupp[index] = SuppPt(tempx, tempy-ystep2, maximass);
         ++Psnum[i];
         ++index;

         Psupp[index] = SuppPt(tempx-xstep2, tempy, maximass/3);
         ++Psnum[i];
         ++index;
         Psupp[index] = SuppPt(tempx, tempy, maximass);
         ++Psnum[i];
         ++index;
      }
      else if (i == 1)
      {
         tempx = 6*xstep;
         tempy = 36*ystep;
         Psupp[index] = SuppPt(tempx, tempy, maximass);
         ++Psnum[i];
         ++index;
         tempy += ystep2;
         Psupp[index] = SuppPt(tempx, tempy, maximass);
         ++Psnum[i];
         ++index;
         tempy += ystep2;
         Psupp[index] = SuppPt(tempx, tempy, maximass);
         ++Psnum[i];
         tempy += ystep2;
         ++index;
         Psupp[index] = SuppPt(tempx-xstep2, tempy, maximass);
         ++Psnum[i];
         tempy += ystep2;
         ++index;
      }
      else if ( i == 2)
      {
         tempx = 9*xstep;
         tempy = 35*ystep;
         
         Psupp[index] = SuppPt(tempx, tempy, maximass);
         ++Psnum[i];
         ++index;
         Psupp[index] = SuppPt(tempx, tempy+ystep2, maximass*1.1);
         ++Psnum[i];
         ++index;
         Psupp[index] = SuppPt(tempx, tempy-ystep2, maximass/3);
         ++Psnum[i];
         ++index;
         Psupp[index] = SuppPt(tempx+xstep2, tempy, maximass);
         ++Psnum[i];
         ++index;
         Psupp[index] = SuppPt(tempx-xstep2, tempy, maximass);
         ++Psnum[i];
         ++index;
      }
      //      else if (i == 3)
      //      {
      //         tempx = 7*xstep;
      //         tempy = 27*ystep;
      //         Psupp[index] = SuppPt(tempx, tempy, maximass);
      //         ++Psnum[i];
      //         ++index;
      //         Psupp[index] = SuppPt(tempx, tempy+ystep, maximass*2);
      //         ++Psnum[i];
      //         ++index;
      //         Psupp[index] = SuppPt(tempx, tempy-ystep, maximass/3);
      //         ++Psnum[i];
      //         ++index;
      //         Psupp[index] = SuppPt(tempx+xstep, tempy, maximass);
      //         ++Psnum[i];
      //         ++index;
      //         Psupp[index] = SuppPt(tempx-xstep, tempy, maximass);
      //         ++Psnum[i];
      //         ++index;
      //      }
      else
      {
         tempx = 11*xstep;
         tempy = 24*ystep;
         Psupp[index] = SuppPt(tempx, tempy-ystep2, maximass);
         ++Psnum[i];
         ++index;
         tempx-=xstep2;
         Psupp[index] = SuppPt(tempx, tempy-ystep2, maximass);
         ++Psnum[i];
         ++index;
         tempx-=xstep2;
         Psupp[index] = SuppPt(tempx, tempy, maximass);
         ++Psnum[i];
         ++index;
         tempx-=xstep2;
         Psupp[index] = SuppPt(tempx, tempy, maximass/3);
         ++Psnum[i];
         ++index;

      }
      
      //Scale to probability distribution
      tempmass = 0;
      for (int j = startindices[i]; j < index; ++j)
      {
         tempmass += Psupp[j].mass;
      }
      for (int j = startindices[i]; j < index; ++j)
      {
         Psupp[j].mass /= tempmass;
      }
      
      endindices[i] = index-1;
   }
   
   indata.close();
   /*
    const int totalsupp = index;
    int Pbarmaxsize = Pnum-1;
    Pbarmaxsize += totalsupp;
    
    
    std::vector<SuppPt>::iterator Psupptest = Psupp.begin();
    for (int i = 0; i < Pnum; ++i)
    {
    validateDist(Psupptest, Psupptest+Psnum[i]);
    Psupptest += Psnum[i];
    }
    
    
    index = 0;
    for (int i = 0; i < Pnum; ++i)
    {
    std::cout << "P_" << i << std::endl;
    for (int j = 0; j < Psnum[i]; ++j)
    {
    std::cout << Psupp[index] << std::endl;
    ++index;
    }
    }
    */
   index = 0;
   for (int i = 0; i < Pnum; ++i)
   {
      std::ostringstream filename2;
      filename2 << "/Users/spatterson/Documents/testi_" << i << "_" << finer <<".txt";
      std::ofstream outresults;
      outresults.open(filename2.str());
      
      for (int j = 0; j < Psnum[i]; ++j)
      {
         outresults << Psupp[index] << '\n';
         ++index;
      }
   }
   
   // remove extra support points from Psupp
   const int totalsupp = index;
   //   Psupp.resize(totalsupp);
   int Pbarmaxsize = Pnum-1;
   Pbarmaxsize += totalsupp;
   
   // ensure P_i sum to exactly 1
   std::vector<SuppPt>::iterator Psupptest = Psupp.begin();
   for (int i = 0; i < Pnum; ++i)
   {
      validateDist(Psupptest, Psupptest+Psnum[i]);
      Psupptest += Psnum[i];
   }
   
   //Calculate number of possible supp pts S0
   unsigned long int S0 = 1;
   for (int i = 0; i < Pnum; ++i)
   {
      S0 *= Psnum[i];
   }
   std::cout << S0 << std::endl;
   
   //Run all combinations of S_i to get possible support points of barycenter
   unsigned long int index2 = 0;
   
   std::vector<SuppPt> Pbar0(S0);
   unsigned long int maxcombos = 1;
   
   index = 0;
   
   std::vector< std::vector<bool> > Pgraph(S0, std::vector<bool>(totalsupp));
   //   std::vector< std::list<costcoord<Pnum> > > Xjcombocost(S0);
   std::vector< std::vector<costcoord<Pnum> > > Xjcombocost(S0, std::vector<costcoord<Pnum> >());
   
   for (int j = 0; j < S0; ++j)
   {
      unsigned long int dupPbarindex = 0;
      double sum1 = 0;
      double sum2 = 0;
      for (int i = 0; i < Pnum; ++i)
      {
         index = indices[i];
         sum1 += Psupp[index].loc1;
         sum2 += Psupp[index].loc2;
      }
      sum1 /= Pnum;
      sum2 /= Pnum;
      
      bool unique = 0;
      for (int l = 0; l < index2; ++l)
      {
         if (Pbar0[l].loc1 == sum1)
         {
            if (Pbar0[l].loc2 == sum2)
            {
               unique = 1;
               ++Pbar0[l].combos;
               if (Pbar0[l].combos > maxcombos)
               {
                  maxcombos = Pbar0[l].combos;
               }
               dupPbarindex = l;
               break;
            }
         }
      }
      
      
      if (unique == 0)
      {
         Pbar0[index2] = SuppPt(sum1, sum2, 0.0, 1);
         dupPbarindex = index2;
         ++index2;
      }
      
      costcoord<Pnum> coord;
      coord.ch = 0;
      for (int i = 0; i < Pnum; ++i)
      {
         index = indices[i];
         Pgraph[dupPbarindex][index] = true;
         double diff1= (Pbar0[dupPbarindex].loc1-Psupp[index].loc1);
         double diff2= (Pbar0[dupPbarindex].loc2-Psupp[index].loc2);
         coord.ch +=lambda[i]*(diff1*diff1+diff2*diff2);
         coord.ind[i] = index;
      }
      
      //      coord.xjindex = dupPbarindex;
      Xjcombocost[dupPbarindex].emplace_back(coord);
      
      //Logic for iterating through all combinations
      int k = Pnum;
      --k;
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
   Pbar0.resize(S0);
   //   Xjcombocost.resize(S0);
   std::cout<< "Size of no-duplicate possible support " << S0 << std::endl;
   //
   //Set up environment
   /*   try { GRBEnv env = GRBEnv();}
    catch(GRBException e){std::cout << "Error Code = " << e.getErrorCode() << std::endl;
    std::cout << e.getMessage() << std::endl;
    }*/
   GRBEnv* env = new GRBEnv();
   GRBModel* model = new GRBModel(*env);
   model->set(GRB_IntParam_Method, 0);
   model->set(GRB_IntParam_Presolve, 0);
   
   GRBVar y;
   //exp2(ik) == mass(ik) will have both w and y
   std::vector<GRBLinExpr> exp2(totalsupp);
   //exp(j,i) == z(j) will only have y
   std::vector< std::vector<GRBLinExpr > > exp(S0, std::vector<GRBLinExpr>(Pnum));
   std::vector< GRBVar > w(S0*maxcombos);
   unsigned long int wnum = 0;
   
   index = 0;
   index2 = 0;
   unsigned long int Syik[S0];
   unsigned long int Syindex = 0;
   std::vector<unsigned long int> Swik(S0*maxcombos);
   for (unsigned long int j = 0; j < S0; ++j)
   {
      unsigned long int yz = 1;
      for (unsigned long int xik = 0; xik < totalsupp; ++xik)
      {
         if (Pgraph[j][xik] == true)
         {
            ++yz;
         }
      }
      if (yz < Pbar0[j].combos)
      {
         index = 0;
         Syik[Syindex] = j;
         ++Syindex;
         // here is where we introduce either y,z or w. This is the y, z code
         for (int i = 0; i < Pnum; ++i)
         {
            while (index <= endindices[i])
            {
               if (Pgraph[j][index] == true)
               {
                  double cost = (Pbar0[j].loc1-Psupp[index].loc1)*(Pbar0[j].loc1-Psupp[index].loc1)+(Pbar0[j].loc2-Psupp[index].loc2)*(Pbar0[j].loc2-Psupp[index].loc2);
                  y = model->addVar(0.0, 1.0, lambda[i]*cost, GRB_CONTINUOUS);
                  exp[j][i] += y;
                  exp2[index] += y;
               }
               ++index;
            }
         }
      }
      else
      {
         // here is the w code
         //         unsigned long int listlength = Pbar0[j].combos;
         //         std::cout << listlength << " " << Xjcombocost[j].size() << std::endl;
         //         std::list<costcoord<Pnum> >::iterator Xjit = Xjcombocost[j].begin();
         std::vector<costcoord<Pnum> >::iterator Xjit = Xjcombocost[j].begin();
         
         //         while (Xjit != Xjcombocost[j].end())
         //         {
         for (int t = 0; t < Pbar0[j].combos; ++t)
         {
            w[wnum] = model->addVar(0.0, 1.0, Xjit->ch, GRB_CONTINUOUS);
            for (int i = 0; i < Pnum; ++i)
            {
               exp2[Xjit->ind[i]] += w[wnum];
            }
            ++Xjit;
            Swik[wnum] = j;
            ++wnum;
         }
         //         unsigned long int listlength = Xjcombocost[j].size();
         //         double costs[listlength];
         //         std::list<costcoord<Pnum> >::iterator Xjit = Xjcombocost[j].begin();
         //         unsigned long int d=0;
         //         while (Xjit != Xjcombocost[j].end())
         //         {
         //            costs[d] = Xjit->ch;
         //            ++Xjit;
         //            ++d;
         //         }
         //         Xjit = Xjcombocost[j].begin();
         //         GRBVar* w;
         //         w = model->addVars(NULL, NULL, costs, NULL, NULL, listlength);
         //         d = 0;
         //         while (Xjit != Xjcombocost[j].end())
         //         {
         //            for (int i = 0; i < Pnum; ++i)
         //            {
         //               exp2[Xjit->ind[i]] += w[d];
         //            }
         //            ++Xjit;
         //            ++d;
         //         }
      }
   }
   
   //   Pgraph.clear(); //This is a minor slowdown that would improve memory efficiency
   //Add z variables.
   
   std::cout << wnum << std::endl;
   std::cout << Syindex << std::endl;
   
   GRBVar* z;
   z = model->addVars(Syindex);
   int zindex;
   for (int i = 0; i < Pnum; ++i)
   {
      zindex = 0;
      for (int jindex = 0; jindex < Syindex; ++jindex)
      {
         model->addConstr(exp[Syik[jindex]][i] == z[zindex]);
         ++zindex;
      }
   }
   for (int i = 0; i < totalsupp; ++i)
   {
      model->addConstr(exp2[i] == Psupp[i].mass);
   }
   
   model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
   
   //      model->update();
   //      std::ostringstream filename;
   //      filename << "/Users/spatterson/Documents/test.lp";
   //      model->write(filename.str());
   t = clock();
   std::cout << "Setup time: " << (double)(t-t0)/CLOCKS_PER_SEC << std::endl;
   
   model->optimize();
   
   t = clock()-t;
   std::cout << "Model Solve time: " << (double)t/CLOCKS_PER_SEC << std::endl;
   
   std::vector<SuppPt>::iterator it = Pbar0.begin();
   std::vector<SuppPt> Pbar01(Pbarmaxsize);
   int S01=0;
   if (model->get(GRB_IntAttr_Status) == 2)
   {
      it = Pbar0.begin();
      for (int j = 0; j < Syindex; ++j)
      {
         tempmass = z[j].get(GRB_DoubleAttr_X);
         if (tempmass != 0)
         {
            Pbar01[S01] = SuppPt(Pbar0[Syik[j]].loc1, Pbar0[Syik[j]].loc2, tempmass);
            ++S01;
         }
      }
      for (int j = 0; j < wnum; ++j)
      {
         tempmass = w[j].get(GRB_DoubleAttr_X);
         if (tempmass != 0)
         {
            //            index = ;
            Pbar01[S01] = SuppPt(Pbar0[Swik[j]].loc1, Pbar0[Swik[j]].loc2, tempmass);
            ++S01;
         }
      }
      
      Pbar01.resize(S01);
      
      //      t = clock()-t;
      //      std::cout << "Get W time: " << (double)t/CLOCKS_PER_SEC << std::endl;
      
      //      S0 = Pbar0.size();*/
      /*
       std::ostringstream filename2;
       filename2 << "/Users/spatterson/Documents/testouttrue" << digit << "_" << Pnum << ".txt";
       std::ofstream outresults;
       outresults.open(filename2.str());
       for (it = Pbar01.begin(); it != Pbar01.end(); ++it)
       {
       outresults << *it << '\n';
       }*/
      //      t = clock()-t;
      //      std::cout << "Write time: " << (double)t/CLOCKS_PER_SEC << std::endl;
      
      
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
   
   t = clock()-t0;
   std::cout << "Total run time: " << (double)t/CLOCKS_PER_SEC << std::endl;
   
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
