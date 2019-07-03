//
//  main.cpp
//  ColGenBothStart
//  This introduces both initial vertices: Greedy-start and 2-app-start
//  No noticable improvement to begin with two.
//  Not fully updated for efficiency.
//
//  Created by Stephan Patterson on 12/7/18.
//  Copyright © 2018 Stephan Patterson. All rights reserved.
//
// Start: Greedy and 2-approx = 2 variables in initial RMP

#include "/Library/gurobi801/mac64/include/gurobi_c++.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <chrono>
#include <iomanip>
#include "SuppPt.hpp"

bool validateDist(const std::vector<SuppPt>::iterator, const std::vector<SuppPt>::iterator);
void mvmult(const int , const int long , const int long , const std::vector< double> &, double * const, const int * const , const int );
void mvmult(const int , const int long , const int long , const std::vector< double> &, std::vector<int> &, double * const, const int * const , const int );
std::vector< double> mvmult(std::vector< std::vector<int> >::iterator,const int, const std::vector< double> &);
std::vector< double> vmmult(const std::vector< std::vector<int> > & , const int , const long int , const double* const );
int vmmult(const std::vector< std::vector<int> > &, const int , const long int , const double* const , std::vector<double> &);
std::vector<double> dualtotals(const double *const, const int, const long int, const int *);
int adjusttotals(std::vector<double> &, const double * const, const int, const long int, const int * const);
int makegreedy(std::vector<SuppPt> &, std::vector<int> & , int *, const int &, std::vector<double> &, const long int &);

int main(int argc, const char * argv[]) {
   
   auto start = std::chrono::steady_clock::now();
   auto t = start;
   
   
   int Pnum =12;//If NumMoMeas = 1, Number of months to process. Will go consecutively from start month, skipping any empty months.
      //Otherwise, Number of measures that will be created. Will combine NumMoMeas of months into each measure
   int startyear = 13;
   std::string startmonth = "Jan";

   std::cout.precision(10);
   int Psnum[Pnum];
   std::vector<SuppPt> Psupp;
   int totalsupp = 0;
   int startindices[Pnum];
   int indices[Pnum];
   int endindices[Pnum];
   int Aprn=0;
   double lambda[Pnum];

   int warmstart = 1; //If 1, the previous solution is saved and reintroduced before ->optimize()
   const double pricetol = -1e-6; // -1e-7 will usually allow duplication before termination. -1e-6 cuts off just a couple iterations with no apparent loss of accuracy in master objective
   int NumMoMeas = 1;//Number of Months in a single Measure. Use 1 for original approach of each month is its own measure
   int comptoTrue = 1;//Whether to do computations to measure error per iter vs. final answer from LP General
   double Trueobj = 0; //Part of comptoTrue. Unneeded otherwise
   std::ofstream outerror; //Part of comptoTrue. Unneeded otherwise
   double currenttol = 1.0;
      
   std::string temp;
   
   std::ifstream indata;
   std::ostringstream fname;
   int year;
   std::string month;
   fname << "DenverCrime.csv";
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
            Psupp.push_back( SuppPt(loc1, loc2 , 1.0, totalsupp));// Using combos to store original index
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
      std::cout << "Month measure: " << i+1 << " Size: " << Psnum[i] << std::endl;
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
   //   Psupp.resize(totalsupp);
   indata.close();
   
   if (comptoTrue == 1)
   {
      std::ostringstream fname2;
      fname2 << "/Users/spatterson/Documents/Column Generation/GenObj" << startyear << "_" << Pnum << ".txt";
      indata.open(fname2.str());
      indata >> Trueobj;
      indata.close();
      std::cout <<Trueobj << std::endl;
      std::ostringstream fname3;
      fname3 << "/Users/spatterson/Documents/Column Generation/ColGenBothYr" << startyear << "_" << Pnum << ".txt";
      outerror.open(fname3.str());
      if (Trueobj == 0)
      {
         std::cout << "Error with reading in true answer." << std::endl;
         return 1;
      }
   }
   
   for (int i = 0; i < 2; ++i)
   {
      Aprn += Psnum[i];
   }
   double sum = 0;
   for (int i = 0; i < Pnum-1; ++i)
   {
      lambda[i] = (double)Psnum[i]/totalsupp;
      sum += lambda[i];
   }
   lambda[Pnum-1] = 1-sum;
   
   
   int Amrn = totalsupp - Aprn;
   
   //std::cout << totalsupp << " " <<Aprn << " " << Amrn << std::endl;
   
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
   long int blockreplength = S0;
   blockreplength /= Psnum[0];
   blockreplength /= Psnum[1];
   
   //Pbar0 contains the possible support points x_j
   std::vector<SuppPt> Pbar0(S0);
   //c constains the cost associated with each grouping in the general formulation
   std::vector<double> c(S0);
   //Amaster contains only those columns which are being used for the master problem
   //this contains the columns from the general formulation, not the DW formulation
   
   //This section generates all c, x_j and does not require repeating
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
      }
      c[j] = cost;
      
      
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
   
   //   std::vector<SuppPt>::iterator Psuppit = Psupp.begin();
   std::vector<SuppPt> Psupp2(totalsupp);
   //   std::copy(Psuppit, Psuppit+totalsupp, Psupp2.begin() ); //Now doing this in loop
   
   GRBEnv* env = new GRBEnv();
   GRBModel* master = new GRBModel(*env);
   //   master->set(GRB_IntParam_Method, 0);
   master->set(GRB_IntParam_Presolve, 0);
   master->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
   master->set(GRB_IntParam_OutputFlag, 0); //Turning off display on screen. Edit to see iterations, objective value, etc.
   
   std::vector< GRBVar > lamk;
   //   std::vector< double > lamkcurrentopt;
   //   std::vector< std::vector<double> > winit;
   //   std::vector< double> cmas;
   
   //Start 2-approx
   GRBModel* model = new GRBModel(*env);
   //   model->set(GRB_IntParam_Method, 0);
   model->set(GRB_IntParam_Presolve, 0);
   
   GRBVar* z;
   z = model->addVars(totalsupp);
   
   int Psmax = *std::max_element(Psnum, Psnum+Pnum);
   std::vector< std::vector< std::vector <GRBVar> >> y(Pnum, std::vector<std::vector<GRBVar> > (totalsupp,std::vector<GRBVar> (Psmax)));
   std::vector< std::vector<double> > costtemp(totalsupp, std::vector<double>(Psmax));
   int index2;
   for (int i = 0; i < Pnum; ++i)
   {
      GRBLinExpr exp2[Psnum[i]];
      // Calculate distances
      index2 = startindices[i];
      for (int j = 0; j < totalsupp; ++j)
      {
         for (int k = 0; k < Psnum[i]; ++k)
         {
            costtemp[j][k] = (Psupp[j].loc1-Psupp[k+index2].loc1)*(Psupp[j].loc1-Psupp[k+index2].loc1)+(Psupp[j].loc2-Psupp[k+index2].loc2)*(Psupp[j].loc2-Psupp[k+index2].loc2);
         }
      }
      for (int j = 0; j < totalsupp; ++j)
      {
         GRBLinExpr exp;
         for (int k = 0; k < Psnum[i]; ++k)
         {
            //            std::ostringstream vname;
            //            vname << "y" << i+1 << "_" << j+1  << "_" << k+1;
            y[i][j][k] = model->addVar(0.0, 1.0, lambda[i]*costtemp[j][k], GRB_CONTINUOUS);
            //            y[j][k].set(GRB_StringAttr_VarName, vname.str());
            exp += y[i][j][k];
            exp2[k] += y[i][j][k];
         }
         model->addConstr(exp == z[j]);
      }
      for (int k = 0; k < Psnum[i]; ++k)
      {
         model->addConstr(exp2[k] == Psupp[index2+k].mass);
      }
   }
   
   model->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
   model->optimize();
   
   std::vector<double> winit(S0,0);
   for (int j = 0; j < totalsupp; ++j)
   {
      if (z[j].get(GRB_DoubleAttr_X) > 0)
      {
         std::vector<SuppPt> xik;
         int ikindex = 0;
         std::vector<int> psnum2(Pnum,0);
         for (int i = 0; i < Pnum; ++i)
         {
            for (int k = 0; k < Psnum[i]; ++k)
            {
               if (y[i][j][k].get(GRB_DoubleAttr_X) > 1e-06)
               {
                  xik.push_back(Psupp[k+startindices[i]]);
                  xik[ikindex].setMass(y[i][j][k].get(GRB_DoubleAttr_X));
                  ++ikindex;
                  ++psnum2[i];
               }
            }
         }
         makegreedy( xik, psnum2, Psnum, Pnum, winit, S0);
      }
   }
   
   //Declare constraints
   GRBLinExpr  exp[Amrn+1];
   //Used for coefficient return from matrix vector product
   double expcoeff[Amrn+1];;
   
   //Calculate coefficient for each numfeas in DW
   double cost = 0;
   int k = 0;
   for (std::vector<double>::iterator wit = winit.begin(); wit != winit.end(); ++wit)
   {
      cost += c[k]*(*wit);
      ++k;
   }
   lamk.push_back(master->addVar(0.0, GRB_INFINITY, cost, GRB_CONTINUOUS));
   mvmult(Amrn, blockreplength, S0, winit, expcoeff, Psnum, Pnum);
   
   for (int j = 0; j < Amrn; ++j)
   {
      exp[j] += expcoeff[j]*lamk[0];
   }
   exp[Amrn] += lamk[0];
   // End of setup for adding 2-approx solution to 1st master problem
   // Now run greedy
   //Amaster contains only those columns which are being used for the master problem
   // currently this contains the columns from the general formulation, not the DW formulation
   std::vector< std::vector<int> > Amaster;
   
   //This section generates c, x_j
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
      }
      c[j] = cost;
      
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
   
   int Amcn = 0;
   std::vector<SuppPt>::iterator Psuppit = Psupp.begin();
   
   //   auto greedyendtime = std::chrono::steady_clock::now();
   //   std::cout << "Total Greedy Construction Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(greedyendtime-dataendtime).count() << "ms" << std::endl;
   
   std::vector< double > expcoeff2;
   std::vector< std::vector<double> > winit2;
   int jindex = 0;
   std::vector< double> cmas;
   for (int n = 0; n < 1; ++n)
   {
      std::vector<double> wtemp;
      std::copy(Psuppit, Psuppit+totalsupp, Psupp2.begin() );
      for (int i = 0; i < Pnum; ++i)
      {
         indices[i] = startindices[i];
      }
      
      while (Psupp2[Psnum[0]-1].mass >1e-14)
      {
         std::vector< int > Amascol(Amrn,0);
         
         double minmass = 1;
         for (int i = 0; i < Pnum; ++i)
         {
            index = indices[i];
            if (Psupp2[index].mass < minmass)
            {
               minmass = Psupp2[index].mass;
            }
         }
         
         wtemp.push_back(minmass);
         
         long int loocur = S0;
         int unki = 0;
         for (int i = 0; i < Pnum; ++i)
         {
            unsigned long int index2 = Psupp2[indices[i]].combos-startindices[i];
            loocur /= Psnum[i];
            unki += loocur*index2;
         }
         cmas.push_back(c[unki]);
         
         for (int i = 2; i < Pnum; ++i)
         {
            Amascol[Psupp2[indices[i]].combos-Aprn] = 1;
         }
         Amaster.push_back(Amascol);
         ++Amcn;
         
         for (int i = 0; i < Pnum; ++i)
         {
            index = indices[i];
            Psupp2[index].mass -= minmass;
            if (Psupp2[index].mass <= 0)
            {
               ++indices[i];
            }
         }
         ++jindex;
         
      }
      winit2.push_back(wtemp);
   }
   
   //Calculate coefficient for each numfeas in DW
   cost = 0;
   std::vector<std::vector<int> >::iterator Amasit = Amaster.begin();
   k = 0;
   cost = 0;
   for (std::vector<double>::iterator wit = winit2[0].begin(); wit != winit2[0].end(); ++wit)
   {
      cost += cmas[k]*(*wit);
      ++k;
   }
   lamk.push_back(master->addVar(0.0, GRB_INFINITY, cost, GRB_CONTINUOUS));
   
   expcoeff2 = mvmult(Amasit, Amcn, winit2[0]);
   Amasit += winit2[0].size();
   
   for (int j = 0; j < Amrn; ++j)
   {
      exp[j] += expcoeff2[j]*lamk[1];
   }
   exp[Amrn] += lamk[1];
   //End of Greedy

   for (int j = 0; j < Amrn; ++j)
   {
      master->addConstr(exp[j] == Psupp[Aprn+j].mass);
   }
   master->addConstr(exp[Amrn] == 1);

   /*
    std::ostringstream filename;
    filename << "/Users/spatterson/Documents/Column Generation/masterfull0.lp";
    master->write(filename.str());
    */
   
   //   std::cout << "Setup time: " << (double)(t-t0)/CLOCKS_PER_SEC << std::endl;
   t = std::chrono::steady_clock::now();
   std::cout << "Total setup time, prior to first solve: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-start).count() << "ms" << std::endl;
   
   
   //   std::cout << "Master Solve: " << std::endl;
   master->optimize();
   double masterobj = master->get(GRB_DoubleAttr_ObjVal);
   
   if (comptoTrue == 1)
   {
      double diff =std::abs(masterobj -Trueobj);
      outerror <<  std::setprecision(10) << diff << ' ';
      if (currenttol > diff)
      {
         //Record times of accuracy
         t = std::chrono::steady_clock::now();
         std::cout << "Time since start: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-start).count() << "ms for accuracy " << diff << " at Iteration number: "<< 0 <<std::endl;
         
      }
      while (currenttol > diff)
      {
         //Lower current accuracy estimate
         currenttol /= 10;
      }
   }
   
   GRBConstr* constrs =master->getConstrs();
   int mastersolvenum = 1;
   auto told = t;
   t = std::chrono::steady_clock::now();
   /*
    std::cout << "Master " << mastersolvenum << " solve time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count() << "ms with objective " << masterobj << std::endl;*/
   
   //Should add a check here that we got an optimal solution
   
   GRBEnv* env2 = new GRBEnv();
   GRBModel* price = new GRBModel(*env2);
   //   price->set(GRB_IntParam_Method, 0);
   price->set(GRB_IntParam_Presolve, 0);
   price->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
   price->set(GRB_IntParam_OutputFlag,0);
   
   GRBLinExpr  expprice[Aprn];
   
   int numuniquecol = Psnum[0]*Psnum[1];
   std::vector<int> mincindices(numuniquecol,0);
   /*std::cout << constrs[0].get(GRB_DoubleAttr_Pi) << std::endl;
    std::cout << constrs[1].get(GRB_DoubleAttr_Pi) << std::endl;
    
    std::cout << constrs[2].get(GRB_DoubleAttr_Pi) << std::endl;*/
   /*
    //This is for writing the dual solution to a file.
    std::ofstream outdata;
    std::ostringstream yfname;
    
    yfname << "/Users/spatterson/Documents/Column Generation/ystore.txt";
    outdata.open(yfname.str());
    
    */
   double* yhat = master->get(GRB_DoubleAttr_Pi, constrs, Amrn+1);
   std::vector<double> yhatprev(Amrn);
   std::vector<double> yhatAmas= dualtotals( yhat,  Amrn, blockreplength, Psnum);
   std::vector<double> cnew;
   
   GRBLinExpr objective = -1*yhat[Amrn];
   price->setObjective(objective);
   
   std::vector<GRBVar> x(numuniquecol);
   
   //Search and use no-duplicate columns in APrice
   int cnewindex = -1;
   int yindex = 0;
   
   for (int kindex = 0; kindex < S0; ++kindex)
   {
      if (kindex % blockreplength == 0)
      {
         yindex = 0;
         cnew.push_back(c[kindex]-yhatAmas[yindex]);
         ++cnewindex;
         mincindices[cnewindex] = kindex;
      }
      else
      {
         if (c[kindex]-yhatAmas[yindex] < cnew[cnewindex])
         {
            cnew[cnewindex] =c[kindex]-yhatAmas[yindex];
            mincindices[cnewindex] = kindex;
         }
      }
   }
   
   
   int jp1 = startindices[0];
   int jp2 = startindices[1];
   
   for (int i = 0; i < numuniquecol; ++i)
   {
      x[i] = price->addVar(0.0, GRB_INFINITY, cnew[i], GRB_CONTINUOUS);
      //Add to constraints without referencing Aprice
      expprice[jp1] += x[i];
      expprice[jp2] += x[i];
      ++jp2;
      if (jp2 == startindices[2])
      {
         jp2 = startindices[1];
         ++jp1;
      }
   }
   
   for (int i = 0; i < Aprn; ++i)
   {
      price->addConstr(expprice[i] == Psupp[i].mass);
   }/*
     std::ostringstream pricefilename;
     pricefilename << "/Users/spatterson/Documents/Column Generation/pricefull0.lp";
     price->write(pricefilename.str());*/
   
   //   std::cout << "Price Solve: " << std::endl;
   price->optimize();
   
   int pricesolvenum = 1;
   told = t;
   t = std::chrono::steady_clock::now();
   
   double pobjvalue = price->get(GRB_DoubleAttr_ObjVal);
   
   //   std::cout << "Price " << pricesolvenum << " solve time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count()  <<"ms with objective" << pobjvalue << std::endl;
   
   //   std::cout<<pobjvalue <<std::endl;
   std::vector<double> xsol(numuniquecol,0);
   
   for (int i = 0; i < numuniquecol; ++i)
   {
      xsol[i] = x[i].get(GRB_DoubleAttr_X);
   }
   
   int iter = 0;
   const long int maxiter = 10000;
   
   std::vector<int> vbases;
   for (int i = 0; i < lamk.size(); ++i)
   {
      vbases.push_back(lamk[i].get(GRB_IntAttr_VBasis));
   }
   vbases.push_back(-1);
   int cbases[Amrn+1];
   for (int i = 0; i <= Amrn; ++i)
   {
      cbases[i] = constrs[i].get(GRB_IntAttr_CBasis);
   }
   
   while (pobjvalue < pricetol)
   {
      if (iter > maxiter)
      {
         std::cout << "Max iterations reached." << std::endl;
         break;
      }
      double cost = 0;
      for (int k = 0; k < numuniquecol; ++k)
      {
         cost += c[mincindices[k]]*xsol[k];
      }
      /*
       told = t;
       t = std::chrono::steady_clock::now();
       std::cout << "cx calculation time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count()  <<"ms" << std::endl;
       */
      
      double expcoeff2[Amrn+1];
      mvmult( Amrn, blockreplength, numuniquecol, xsol, mincindices, expcoeff2, Psnum, Pnum);
      expcoeff2[Amrn] = 1.0;
      
      /*
       told = t;
       t = std::chrono::steady_clock::now();
       std::cout << "Constraint coefficient calculation time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count()  <<"ms" << std::endl;*/
      
      /*
       lamkcurrentopt.clear();
       for (int i = 0; i < lamk.size(); ++i)
       {
       lamkcurrentopt.push_back(lamk[i].get(GRB_DoubleAttr_X));
       }*/
      
      GRBColumn newmastercol;
      newmastercol.addTerms(expcoeff2, master->getConstrs(), Amrn+1);
      lamk.push_back(master->addVar(0.0, GRB_INFINITY, cost, GRB_CONTINUOUS, newmastercol));//adding this variable is when the auto warm start is lost
      
      /*
       told = t;
       t = std::chrono::steady_clock::now();
       std::cout << "Variable add to model time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count()  <<"ms" << std::endl;
       */
      /*
       std::ostringstream masfilename2;
       masfilename2 << "/Users/spatterson/Documents/Column Generation/masterfull" << iter+1 <<".lp";
       
       master->write(masfilename2.str());
       */
      //      std::cout << "Master Solve: " << std::endl;
      
      //This code should warm start
      if (warmstart == 1 )
      {
         master->update();
         for (int i = 0; i < lamk.size(); ++i)
         {
            lamk[i].set(GRB_IntAttr_VBasis, vbases[i]);
         }
         for (int i = 0; i <=Amrn; ++i)
         {
            constrs[i].set(GRB_IntAttr_CBasis, cbases[i]);
         }
         /*
          master->update();
          
          for (int i = 0; i < lamk.size()-1; ++i)
          {
          lamk[i].set(GRB_DoubleAttr_PStart, lamkcurrentopt[i]);
          }
          lamk[lamk.size()-1].set(GRB_DoubleAttr_PStart,0);
          constrs =master->getConstrs();
          for (int i = 0; i <= Amrn; ++i )
          {
          constrs[i].set(GRB_DoubleAttr_DStart,yhat[i]);
          }*/
      }
      /*
       told = t;
       t = std::chrono::steady_clock::now();
       std::cout << "Reintroduce warm solution time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count()  <<"ms" << std::endl;
       */
      master->optimize();
      masterobj = master->get(GRB_DoubleAttr_ObjVal);
      /*
       told = t;
       t = std::chrono::steady_clock::now();
       std::cout << "Master solve time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count()  <<"ms" << std::endl;
       */
      if ((comptoTrue == 1) && (iter % 5 == 0))
      {
         double diff =std::abs(masterobj -Trueobj);
         outerror <<  std::setprecision(10) << diff << ' ';
         if (currenttol > diff)
         {
            //Record times of accuracy
            t = std::chrono::steady_clock::now();
            std::cout << "Time since start: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-start).count() << "ms for accuracy " << diff << " at Iteration number: "<< iter << std::endl;
            std::cout << "Previous Price Objective was: " <<pobjvalue << std::endl;
            
         }
         while (currenttol > diff)
         {
            //Lower current accuracy estimate
            currenttol /= 10;
         }
      }
      
      if (warmstart == 1)
      {
         for (int i = 0; i < lamk.size(); ++i)
         {
            vbases[i] = (lamk[i].get(GRB_IntAttr_VBasis));
         }
         vbases.push_back(-1);
         for (int i = 0; i <= Amrn; ++i)
         {
            cbases[i] = constrs[i].get(GRB_IntAttr_CBasis);
         }
      }
      
      ++mastersolvenum;
      /*
      told = t;
      t = std::chrono::steady_clock::now();
      std::cout << "Master " << mastersolvenum << " solve time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count() << "ms with objective " << masterobj << std::endl;*/
      for (int i = 0; i < Amrn; ++i)
      {
         yhatprev[i] = yhat[i];
      }
      
      yhat = master->get(GRB_DoubleAttr_Pi, constrs, Amrn+1);
      double ychange[Amrn];
      for (int i = 0; i < Amrn; ++i)
      {
         ychange[i] =yhat[i];
         ychange[i] -= yhatprev[i];
      }
      /*
       This bit of code writes the dual solutions to a file.
       for (int i = 0; i < Amrn+1; ++i)
       {
       outdata << yhat[i] << " ";
       }
       outdata << std::endl;
       outdata << std::endl;
       */
      
      
      //      told = t;
      //      t = std::chrono::steady_clock::now();
      
      adjusttotals(yhatAmas, ychange, Amrn, blockreplength,Psnum);
      //      told = t;
      //      t = std::chrono::steady_clock::now();
      //      std::cout << "Calculate product yTA new " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count()  <<"ms" << std::endl;
 
      /*
       told = t;
       t = std::chrono::steady_clock::now();
       std::cout << "Calculate product yTA " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count()  <<"ms" << std::endl;
       */
      /*
       if (yhatPrev == yhatAmas)
       {
       std::cout << "Equivalent dual solution will now produce same pricing solution. Algorithm stalls." << std::endl;
       break;
       }*/
      if (iter > 1)
      {
         int equal = 1;
         for (int i = 0; i < Amrn; ++i)
         {
            if (ychange[i] != 0)
            {
               equal = 0;
               break;
            }
         }
         if (equal == 1)
         {
            std::cout << "Equivalent dual solution will now produce same pricing solution. Algorithm stalls." << std::endl;
            break;
         }
      }
      cnew.clear();
      
      GRBLinExpr objective = -1*yhat[Amrn];
      /*
       told = t;
       t = std::chrono::steady_clock::now();
       std::cout << "Check duplicate & start objective " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count()  <<"ms" << std::endl;
       */
      
      
      //Search and use no-duplicate columns in APrice
      int cnewindex = -1;
      int yindex = 0;
      for (int kindex = 0; kindex < S0; ++kindex)
      {
         if (kindex % blockreplength == 0)
         {
            yindex = 0;
            cnew.push_back(c[kindex]-yhatAmas[yindex]);
            ++cnewindex;
            mincindices[cnewindex] = kindex;
         }
         else
         {
            if (c[kindex]-yhatAmas[yindex] < cnew[cnewindex])
            {
               cnew[cnewindex] =c[kindex]-yhatAmas[yindex];
               mincindices[cnewindex] = kindex;
            }
         }
         ++yindex;
      }
      for (int j = 0; j < numuniquecol; ++j)
      {
         objective += cnew[j]*x[j];
      }
      
      
      price->setObjective(objective);
      /*
       told = t;
       t = std::chrono::steady_clock::now();
       std::cout << "Set objective " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count()  <<"ms" << std::endl;
       */
      /*
       std::ostringstream prifilename2;
       prifilename2 << "/Users/spatterson/Documents/Column Generation/pricefull" << iter+1 <<".lp";
       
       price->write(prifilename2.str());
       */
      ++pricesolvenum;
      //      std::cout << "Price Solve: " << std::endl;
      price->optimize();
      pobjvalue = price->get(GRB_DoubleAttr_ObjVal);
      /*
       told = t;
       t = std::chrono::steady_clock::now();
       std::cout << "Price " << pricesolvenum << " solve time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count()  <<"ms with objective " << pobjvalue << std::endl;*/
      
      for (int i = 0; i < numuniquecol; ++i)
      {
         xsol[i] = x[i].get(GRB_DoubleAttr_X);
         /*
        if (xsol[i] != 0)
        {
        std::cout << i << " " << xsol[i] << " " << mincindices[i] << std::endl;
        }*/
      }
      
      ++iter;
   }
   told = t;
   t = std::chrono::steady_clock::now();
   std::cout << "Total Running time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-start).count()  <<"ms" << std::endl;
   std::cout << "Last price obj: " << pobjvalue << std::endl;
   std::cout << "Optimal value: " << masterobj << " found after " << iter << " solves." << std::endl;
   
   std::ostringstream masfilename2;
   masfilename2 << "/Users/spatterson/Documents/Column Generation/masterfull" << Pnum << "_"<< iter+1 <<".lp";
   
   master->write(masfilename2.str());
   master->terminate();
   delete master;
   
   std::ostringstream prifilename2;
   prifilename2 << "/Users/spatterson/Documents/Column Generation/pricefull" << Pnum << "_"<< iter+1 <<".lp";
   
   price->write(prifilename2.str());
   price->terminate();
   delete price;
   delete env;
   
   outerror.close();
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

std::vector< double> mvmult(const std::vector< std::vector<int> > &A, const int Ars, const int long blocklength, const int long Acs, const std::vector< double> &x)
{
   std::vector<double > Ax(Ars);
   for (int i = 0; i < Ars; ++i)
   {
      Ax.push_back(0);
      for (int j = 0; j < Acs; ++j)
      {
         if (A[j][i] != 0)
         {
            Ax[i] += (double)A[j][i]*x[j];
         }
      }
   }
   return Ax;
}

std::vector< double> mvmult(std::vector<std::vector<int> >::iterator A, const int colnum, const std::vector< double> &x)
{
   //This gives the row number
   unsigned long int retlength = A->size();
   std::vector<double > Ax(retlength,0);
   std::vector<std::vector<int> >::iterator A2temp = A;
   
   for (int i = 0; i < colnum; ++i)
   {
      int j = 0;
      for (std::vector<int>::iterator Atemp = A2temp->begin(); Atemp != A2temp->end(); ++Atemp)
      {
         if (*Atemp != 0)
         {
            Ax[j] += x[i];
         }
         ++j;
      }
      ++A2temp;
   }
   return Ax;
}

std::vector< double> vmmult(const std::vector< std::vector<int> > &A, const int Ars, const long int Acs, const double* const yT)
{
   std::vector<double > yA(Acs,0);
   for (int j = 0; j < Ars; ++j)
   {
      if (yT[j] != 0)
      {
         for (int i = 0; i < Acs; ++i)
         {
            if (A[i][j] != 0)
            {
               yA[i] += yT[j];
            }
         }
      }
   }
   return yA;
}

int vmmult(const std::vector< std::vector<int> > &A, const int Ars, const long int Acs, const double* const yT, std::vector<double> &yA)
{
   for (int i = 0; i < Acs; ++i)
   {
      yA[i] = 0;
   }
   for (int j = 0; j < Ars; ++j)
   {
      if (yT[j] != 0)
      {
         for (int i = 0; i < Acs; ++i)
         {
            if (A[i][j] != 0)
            {
               yA[i] += yT[j];
            }
         }
      }
   }
   return 0;
}

void mvmult(const std::vector< std::vector<int> > &A, const int Ars, const int long blocklength, const int long Acs, const std::vector< double> &x, double * Ax)
{
   int colindex = 0;
   for (int i = 0; i < Ars; ++i)
   {
      Ax[i] = 0;
   }
   for (int j = 0; j < Acs; ++j)
   {
      if (j % blocklength == 0)
      {
         colindex = 0;
      }
      if (x[j] != 0)
      {
         for (int i = 0; i < Ars; ++i)
         {
            if (A[colindex][i] != 0)
            {
               Ax[i] += x[j];
            }
         }
      }
      ++colindex;
   }
}

void mvmult(const int Ars, const int long blocklength, const int long Acs, const std::vector< double> &x, double * const Ax, const int * const Psnum, const int Pnum)
{
   unsigned long int index = 0;
   unsigned long int startindex = 0;
   unsigned long int ucn = 0;
   unsigned long int blockilength = blocklength;
   int blockinum = 0;
   for (int i = 0; i < Ars; ++i)
   {
      Ax[i] = 0;
   }
   for (int j = 0; j < Acs; ++j)
   {
      if (x[j] != 0)
      {
         ucn = j-blocklength*blockinum;
         blockilength = blocklength;
         startindex = 0;
         for (int i = 2; i < Pnum; ++i)
         {
            blockilength /= Psnum[i];
            index = ucn / blockilength;
            Ax[startindex + index] += x[j];
            ucn -= index*blockilength;
            startindex += Psnum[i];
         }
      }
      if ((j+1)%blocklength == 0)
      {
         ++blockinum;
      }
   }
}

void mvmult(const int Ars, const int long blocklength, const int long Acs, const std::vector< double> &x, std::vector<int> &mcind, double * const Ax, const int * const Psnum, const int Pnum)
{
   unsigned long int index = 0;
   unsigned long int startindex = 0;
   unsigned long int ucn = 0;
   unsigned long int blockilength = blocklength;
   for (int i = 0; i < Ars; ++i)
   {
      Ax[i] = 0;
   }
   for (int j = 0; j < Acs; ++j)
   {
      if (x[j] != 0)
      {
         ucn = mcind[j]%blocklength;
         blockilength = blocklength;
         startindex = 0;
         
         for (int i = 2; i < Pnum; ++i)
         {
            blockilength /= Psnum[i];
            index = ucn / blockilength;
            Ax[startindex + index] += x[j];
            ucn -= index*blockilength;
            startindex += Psnum[i];
         }
      }
   }
}

std::vector<double> dualtotals(const double * const yhat, const int Amrn, const long int blocklength, const int * Psnum)
{
   std::vector<double> ytotals(blocklength,0);
   int i = 2;
   int totalPs = 0;
   long int startindex = 0;
   long int blockinum = 1;
   long int blockilength = blocklength;
   long int blocki1length = blockilength/Psnum[i];
   long int j = 0;
   long int index = 0;
   while (j < Amrn)
   {
      if (yhat[j] != 0)
      {
         index = 0;
         for (int k = 0; k < blockinum; ++k)
         {
            for (int l = 0; l < blocki1length; ++l)
            {
               ytotals[startindex+index] += yhat[j];
               ++index;
            }
            index = blockilength*(k+1);
         }
      }
      ++j;
      startindex += blocki1length;
      if (j-totalPs == Psnum[i])
      {
         blockinum *= Psnum[i];
         blockilength /= Psnum[i];
         startindex = 0;
         totalPs += Psnum[i];
         ++i;
         blocki1length /= Psnum[i];
      }
   }
   return ytotals;
}

int adjusttotals(std::vector<double> &yhatAmas, const double * const ychange, const int Amrn, const long int blocklength, const int * const Psnum)
{
   int i = 2;
   int totalPs = 0;
   long int startindex = 0;
   long int blockinum = 1;
   long int blockilength = blocklength;
   long int blocki1length = blockilength/Psnum[i];
   long int j = 0;
   long int index = 0;
   while (j < Amrn)
   {
      if (ychange[j] != 0)
      {
         index = 0;
         for (int k = 0; k < blockinum; ++k)
         {
            for (int l = 0; l < blocki1length; ++l)
            {
               yhatAmas[startindex+index] += ychange[j];
               ++index;
            }
            index = blockilength*(k+1);
         }
      }
      ++j;
      startindex += blocki1length;
      if (j-totalPs == Psnum[i])
      {
         blockinum *= Psnum[i];
         blockilength /= Psnum[i];
         startindex = 0;
         totalPs += Psnum[i];
         ++i;
         blocki1length /= Psnum[i];
      }
   }
   return 0;
}

int makegreedy(std::vector<SuppPt> &Psupp2, std::vector<int> &Psnum2, int * Psnum, const int &Pnum, std::vector<double> &winit,const long int &S0)
{
   std::vector<int> indices(Pnum,0);
   std::vector<int> startindices(Pnum,0);
   
   for (int i = 1; i < Pnum; ++i)
   {
      indices[i] = indices[i-1];
      indices[i] += Psnum2[i-1];
      startindices[i] = indices[i];
   }
   
   while (Psupp2[Psnum2[0]-1].mass >1e-12)
   {
      double minmass = 1;
      for (int i = 0; i < Pnum; ++i)
      {
         if (Psupp2[indices[i]].mass < minmass)
         {
            minmass = Psupp2[indices[i]].mass;
         }
      }
      
      long int loocur = S0;
      int unki = 0;
      int index = 0;
      for (int i = 0; i < Pnum; ++i)
      {
         unsigned long int index2 = Psupp2[indices[i]].combos-index;
         loocur /= Psnum[i];
         unki += loocur*index2;
         index += Psnum[i];
         Psupp2[indices[i]].mass -= minmass;
         if (Psupp2[indices[i]].mass <= 1e-12)
         {
            ++indices[i];
         }
      }
      
      winit[unki] = minmass;
   }
   return 0;
}
