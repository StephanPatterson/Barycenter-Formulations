//
//  main.cpp
//  Greedy Start
//  No Pricing Problem Sorting
//
//  Created by Stephan Patterson on 12/4/18.
//  Copyright © 2018 Stephan Patterson. All rights reserved.
//

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
//std::vector< double> mvmult(const std::vector< std::vector<int> > &, const int, const long int, const long int, const std::vector< double> &);
//void mvmult(const std::vector< std::vector<int> > &, const int, const long int, const long int, const std::vector< double> &, double *);
void mvmult(const int , const int long , const int long , const std::vector< double> &, double * , const int * const , const int );
void mvmult(const int , const int long , const int long , const std::vector< double> &, std::vector<int> &, double * const, const int * const , const int );
std::vector< double> mvmult(std::vector< std::vector<int> >::iterator,const int, const std::vector< double> &);
std::vector< double> vmmult(const std::vector< std::vector<int> > & , const int , const long int , const double* const );
int vmmult(const std::vector< std::vector<int> > &, const int , const long int , const double* const , std::vector<double> &);
std::vector<double> dualtotals(const double *const, const int, const long int, const int *, const int &);
int adjusttotals(std::vector<double> &, const double * const, const int, const long int, const int * const, const int &);
int MeasToTXTFile(std::vector<SuppPt> &,const int &, const int &, const int &, const int &, const int &);
int MeasToCommaFile(std::vector<SuppPt> &,const int &, const int &, const int &, const int &, const int &);

int main(int argc, const char * argv[]) {
   
   auto start = std::chrono::steady_clock::now();
   auto t = start;
   std::vector<double> mastert;
   std::vector<double> pricet;
   
   std::cout.precision(10);
   const int Pnum =12;//If NumMoMeas = 1, Number of months to process. Will go consecutively from start month, skipping any empty months.
   //Otherwise, Number of measures that will be created. Will combine NumMoMeas of months into each measure
   int startyear = 15;
   std::string startmonth = "Jan";
   int NumMoMeas = 1;//Number of Months in a single Measure. Use 1 for original approach of each month is its own measure
   
   int Psnum[Pnum];
   int numinprice = 2; // Numbers not = 1,2 have not been fully implemented
   std::vector<SuppPt> Psupp;
   int totalsupp = 0;
   int startindices[Pnum];
   int indices[Pnum];
   int endindices[Pnum];
   int Aprn=0;
   double lambda[Pnum];
   //   std::vector<std::string> monthnames = {"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};
   
   int warmstart = 1; //If 1, the previous solution is saved and reintroduced before ->optimize()
   //   int dotoy = 0; //If 0, runs denver data. if 1, runs toy example
   const double pricetol = -1e-6; // -1e-7 will usually allow duplication before termination. -1e-6 cuts off just a couple iterations with no apparent loss of accuracy in master objective
   
   int comptoTrue = 0;// If 1, will read in the objective value produced by the general formulation for the given year, and calculate error at each iteration
   double Trueobj = 0;
   std::ofstream outerror;
   double currenttol = 1.0;
   
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

   for (int i = 0; i < Pnum-2; ++i)
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
      //      std::cout << "Month measure: " << i+1 << " Size: " << Psnum[i] << std::endl;
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
   //   NumMoMeas = 1;
   for (int i = Pnum-2; i < Pnum; ++i)
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
      //      std::cout << "Month measure: " << i+1 << " Size: " << Psnum[i] << std::endl;
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
   
   indata.close();
   
   if (comptoTrue == 1)
   {
      std::ostringstream fname2;
      fname2 << "/GenObj" << startyear << "_" << Pnum << ".txt";
      indata.open(fname2.str());
      indata >> Trueobj;
      indata.close();
      std::cout <<Trueobj << std::endl;
      std::ostringstream fname3;
      fname3 << "/ColGenGRandomYr" << startyear << "_" << Pnum << ".txt";
      outerror.open(fname3.str());
      if (Trueobj == 0)
      {
         std::cout << "Error with reading in true answer." << std::endl;
         return 1;
      }
   }
   
   for (int i = 0; i < numinprice; ++i)
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

   /*
    auto dataendtime = std::chrono::steady_clock::now();
    std::cout << "Total Data Read-in Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(dataendtime-start).count() << "ms" << std::endl;
    */
   int Amrn = totalsupp - Aprn;
   
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
   // currently this contains the columns from the general formulation, not the DW formulation
   std::vector< std::vector<int> > Amaster;
   
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
         //         std::cout << index << " " ;
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
   std::vector<SuppPt> Psupp2(totalsupp);
   //   std::copy(Psuppit, Psuppit+totalsupp, Psupp2.begin() ); //Now doing this in loop
   /*
    auto greedyendtime = std::chrono::steady_clock::now();
    std::cout << "Total Greedy Construction Time: " << std::chrono::duration_cast<std::chrono::milliseconds>(greedyendtime-dataendtime).count() << "ms" << std::endl;
    */
   
   GRBEnv* env = new GRBEnv();
   GRBModel* master = new GRBModel(*env);
   master->set(GRB_IntParam_Method, 0);
   master->set(GRB_IntParam_Presolve, 0);
   master->set(GRB_IntAttr_ModelSense, GRB_MINIMIZE);
   //   master->set(GRB_IntParam_UpdateMode,0); //why did this break adding constraints?
   master->set(GRB_IntParam_OutputFlag, 0); //Turning off display on screen. Disable to see iterations, objective value, etc.
   
   std::vector< GRBVar > lamk;
   //   std::vector< double > lamkcurrentopt;
   std::vector< std::vector<double> > winit;
   int jindex = 0;
   std::vector< double> cmas;
//   std::ofstream outdest;
//   outdest.precision(10);
//   std::ostringstream fnamedest;
//   fnamedest << "/DenDataGreedyTo" << startyear << "_" << Pnum << "_" << 1 << ".txt";
//   outdest.open(fnamedest.str());
//   for (int i = 0; i < Pnum; ++i)
//   {
//      outdest <<"LOCX" <<i<<",LOCY"<<i;
//      if (i < Pnum-1)
//      { outdest<<",";
//      }
//   }
//   outdest <<'\n';
   
      std::vector<double> wtemp;
//      std::vector<SuppPt> Wsupp;//having both of these doesn't make sense, long term. Need full support for visualizations only. Otherwise vector of doubles sufficient.
      std::copy(Psuppit, Psuppit+totalsupp, Psupp2.begin() );
      for (int i = 0; i < Pnum; ++i)
      {
         indices[i] = startindices[i];
      }
      
      while (Psupp2[totalsupp-1].mass >1e-14)
      {
         //         std::cout << "mass: " << Psupp2[Psnum[0]-1].mass << std::endl;
         std::vector< int > Amascol(Amrn,0);
         
         double minmass = 1;
         for (int i = 2; i < Pnum; ++i)
         {
            index = indices[i];
//            outdest <<Psupp2[index].loc1<<","<<Psupp2[index].loc2;
//            if (i < Pnum-1)
//            {
//               outdest<<",";
//            }
            if (Psupp2[index].mass < minmass)
            {
               minmass = Psupp2[index].mass;
            }
         }
//         outdest << '\n';
         
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
//         Wsupp.push_back(Pbar0[unki]);
//         Wsupp[Amcn].mass = minmass;
         
         for (int i = 2; i < Pnum; ++i)
         {
            Amascol[Psupp2[indices[i]].combos-Aprn] = 1;
         }
         Amaster.push_back(Amascol);
         ++Amcn;
         
         for (int i = 2; i < Pnum; ++i)
         {
            index = indices[i];
            Psupp2[index].mass -= minmass;
            if (Psupp2[index].mass < 1e-14)
            {
               ++indices[i];
            }
         }
         ++jindex;
         
      }

      winit.push_back(wtemp);
//      MeasToCommaFile(Wsupp,Amcn,startyear,Pnum,-1,0);
//   outdest.close();
   
   //Declare constraints
   GRBLinExpr  exp[Amrn+1];
   //Used for coefficient return from matrix vector product
   std::vector< double > expcoeff;
   
   //Calculate coefficient for each numfeas in DW
   double cost = 0;
   std::vector<std::vector<int> >::iterator Amasit = Amaster.begin();
   int k = 0;
   for (std::vector<double>::iterator wit = winit[0].begin(); wit != winit[0].end(); ++wit)
   {
      cost += cmas[k]*(*wit);
      ++k;
   }
   lamk.push_back(master->addVar(0.0, GRB_INFINITY, cost, GRB_CONTINUOUS));
   
   expcoeff = mvmult(Amasit, Amcn, winit[0]);
   Amasit += winit[0].size();
   
   for (int j = 0; j < Amrn; ++j)
   {
      exp[j] += expcoeff[j]*lamk[0];
   }
   exp[Amrn] += lamk[0];
   for (int j = 0; j < Amrn; ++j)
   {
      master->addConstr(exp[j] == Psupp[Aprn+j].mass);
   }
   master->addConstr(exp[Amrn] == 1);
//   std::ostringstream filename;
//   filename << "/masterfull0.lp";
//   master->write(filename.str());
   
   /*
    t = std::chrono::steady_clock::now();
    std::cout << "Total setup time, prior to first solve: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-start).count() << "ms" << std::endl;
    */
   
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
   /*
    auto told = t;
    t = std::chrono::steady_clock::now();
    
    mastert.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count());
    std::cout << "Master " << mastersolvenum << " solve time: " << mastert[0] << "ms with objective " << masterobj << std::endl;
    */
   
   //Should add a check here that we got an optimal solution
   
   GRBEnv* env2 = new GRBEnv();
   GRBModel* price = new GRBModel(*env2);
   price->set(GRB_IntParam_Method, 0);
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
    
    yfname << "/ystore.txt";
    outdata.open(yfname.str());
    
    */
   double* yhat = master->get(GRB_DoubleAttr_Pi, constrs, Amrn+1);
   std::vector<double> yhatprev(Amrn);
   /*
    for (int i = 0; i < Amrn; ++i)
    {
    yhatprev[i] = yhat[i];
    }
    */
   
   /*
    for (int i = 0; i < Amrn+1; ++i)
    {
    outdata << yhat[i] << " ";
    }
    outdata << std::endl;
    
    outdata << std::endl;
    */
   std::vector<double> yhatAmas= dualtotals( yhat,  Amrn, blockreplength, Psnum, Pnum);
   
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
      ++yindex;
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
     pricefilename << "/pricefull0.lp";
     price->write(pricefilename.str());*/
   
   //   std::cout << "Price Solve: " << std::endl;
   price->optimize();
   
   int pricesolvenum = 1;
   double pobjvalue = price->get(GRB_DoubleAttr_ObjVal);
   
   /*
    told = t;
    t = std::chrono::steady_clock::now();
    
    pricet.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count());
    std::cout << "Price " << pricesolvenum << " solve time: " << pricet[0] << "ms with objective " << pobjvalue << std::endl;
    */
   
   //   std::cout<<pobjvalue <<std::endl;
   std::vector<double> xsol(numuniquecol,0);
   
   for (int i = 0; i < numuniquecol; ++i)
   {
      xsol[i] = x[i].get(GRB_DoubleAttr_X);
   }
   /*
    told = t;
    t = std::chrono::steady_clock::now();
    std::cout << "Extract Price Solution time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count()  <<"ms" << std::endl;
    */
   
   int iter = 0;
   const long int maxiter = 10000;
   
   std::vector<int> vbases;
   vbases.push_back(0);
   vbases.push_back(-1);
   int cbases[Amrn+1];
   for (int i = 0; i <= Amrn; ++i)
   {
      cbases[i] = constrs[i].get(GRB_IntAttr_CBasis);
   }
   
   /*
    told = t;
    t = std::chrono::steady_clock::now();
    std::cout << "Save First Sol for Warm Start time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count()  <<"ms" << std::endl;*/
   
   while (pobjvalue < pricetol)
   {
      if (iter >= maxiter)
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
       std::cout << "cx calculation time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count()  <<"ms" << std::endl;*/
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
       std::cout << "Variable add to model time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count()  <<"ms" << std::endl;*/
      
      /*
       std::ostringstream masfilename2;
       masfilename2 << "/masterfull" << iter+1 <<".lp";
       
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
       std::cout << "Reintroduce warm solution time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count()  <<"ms" << std::endl;*/
      
      master->optimize();
      masterobj = master->get(GRB_DoubleAttr_ObjVal);
      
      /*
       told = t;
       t = std::chrono::steady_clock::now();
       
       mastert.push_back(std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count());
       std::cout << "Master solve time: " << mastert[iter]  <<"ms on iteration " << iter << " with objective: " << masterobj << std::endl;*/
      
      if (comptoTrue == 1 && (iter %5 == 0))
      {
         double diff =std::abs(masterobj -Trueobj);
         outerror <<  std::setprecision(10) << diff << ' ';
         if (currenttol > std::abs(diff))
         {
            //Record times of accuracy
            auto t2 = std::chrono::steady_clock::now();
            std::cout << "Time since start: " << std::chrono::duration_cast<std::chrono::milliseconds>(t2-start).count() << "ms for accuracy " << diff << " at Iteration number: "<< iter << std::endl;
            std::cout << "Previous Price Objective was: " <<pobjvalue << std::endl;
            
         }
         while (currenttol > std::abs(diff))
         {
            //Lower current accuracy estimate
            currenttol /= 10;
         }
      }
      
      if (warmstart == 1)
      {
         //         constrs =master->getConstrs();
         //         vbases.clear();
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
       std::cout << "Warm Start Save Time:" << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count() << "ms."<< std::endl;*/
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
      
      adjusttotals(yhatAmas, ychange, Amrn, blockreplength, Psnum, Pnum);
      
      /*      std::cout << "Start of yhatA comparison. New, then previous." << std::endl;
       for (int i = 0; i <blockreplength; ++i)
       {
       std::cout << yhatAmas[i] << " " << yhatAmasprev[i] << std::endl;
       }*/
      /*
       told = t;
       t = std::chrono::steady_clock::now();
       std::cout << "Calculate product yTA new " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count()  <<"ms" << std::endl;*/
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
       std::cout << "Check duplicate & start objective " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count()  <<"ms" << std::endl;*/
      //Search and use no-duplicate columns in APrice
      //Start of original price objective construction
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
            double temp =c[kindex]-yhatAmas[yindex];
            if ( temp < cnew[cnewindex])
            {
               cnew[cnewindex] =temp;
               mincindices[cnewindex] = kindex;
            }
         }
         ++yindex;
      }
      for (int j = 0; j < numuniquecol; ++j)
      {
         objective += cnew[j]*x[j];
      }
      //End of original price objective construction
      /*
       std::cout << "Start of c compares. New, then previous. " << std::endl;
       for (int j = 0; j < numuniquecol; ++j)
       {
       std::cout << cnew[j] << " " << cprev[j] << std::endl;
       }
       */
      price->setObjective(objective);
      /*
       told = t;
       t = std::chrono::steady_clock::now();
       std::cout << "Set new price objective " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count()  <<"ms" << std::endl;
       */
      
      
      /*
       std::ostringstream prifilename2;
       prifilename2 << "/pricefull" << iter+1 <<".lp";
       
       price->write(prifilename2.str());
       */
      ++pricesolvenum;
      //      std::cout << "Price Solve: " << std::endl;
      price->optimize();
      pobjvalue = price->get(GRB_DoubleAttr_ObjVal);
      
      /*
       told = t;
       t = std::chrono::steady_clock::now();
       std::cout << "Price " << pricesolvenum << " solve time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count()  <<"ms with objective " << pobjvalue << std::endl;
       */
      
      for (int i = 0; i < numuniquecol; ++i)
      {
         xsol[i] = x[i].get(GRB_DoubleAttr_X);
      }
      
      /*
       told = t;
       t = std::chrono::steady_clock::now();
       std::cout << "Extract Price Solution time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-told).count()  <<"ms" << std::endl;
       */
      ++iter;
   }
   
   
   //   told = t;
   t = std::chrono::steady_clock::now();
   std::cout << "Total Running time: " << std::chrono::duration_cast<std::chrono::milliseconds>(t-start).count()  <<"ms" << std::endl;
   std::cout << "Last price obj: " << pobjvalue << std::endl;
   std::cout << "Optimal value: " << masterobj << " found after " << iter << " solves." << std::endl;
   
//   std::ostringstream masfilename2;
//   masfilename2 << "/masterfull" << Pnum << "_"<< iter+1 <<".lp";
//   master->write(masfilename2.str());
   master->terminate();
   delete master;
   
//   std::ostringstream prifilename2;
//   prifilename2 << "/pricefull" << Pnum << "_"<< iter+1 <<".lp";
//   price->write(prifilename2.str());
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
         //         std::cout << *Atemp << std::endl;
         if (*Atemp != 0)
         {
            Ax[j] += x[i];
            //            std::cout << j << " " << *Atemp << " " << x[i] << " ";
         }
         ++j;
      }
      ++A2temp;
   }/*
     std::cout << std::endl;
     for (int j = 0; j < retlength; ++j)
     {
     std::cout << Ax[j] << std::endl;
     }*/
   return Ax;
}
/*
 std::vector< double> vmmult(const std::vector< std::vector<int> > &A, const int Ars, const long int Acs, const double* yT)
 {
 std::vector<double > yA(Acs);
 for (int i = 0; i < Acs; ++i)
 {
 yA.push_back(0);
 for (int j = 0; j < Ars; ++j)
 {
 if (A[i][j] != 0)
 {
 yA[i] += yT[j];
 }
 }
 //      std::cout << yA[i] << std::endl;
 }
 return yA;
 }*/

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


/*
 void mvmult(const std::vector< std::vector<int> > A, const int Ars, const int long Acs, const std::vector< double> x, double * Ax)
 {
 for (int i = 0; i < Ars; ++i)
 {
 Ax[i] = 0;
 for (int j = 0; j < Acs; ++j)
 {
 if (A[j][i] != 0)
 {
 Ax[i] += x[j];
 }
 }
 }
 }
 */
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

/*
 void mvmult(const int Ars, const int long blocklength, const int long Acs, const std::vector< double> &x, double * Ax, const int * const Psnum, const int Pnum)
 {
 int indices[Pnum];
 int endindices[Pnum];
 int startindices[Pnum];
 indices[2] = 0;
 startindices[2] = 0;
 for (int i = 3; i < Pnum; ++i)
 {
 indices[i] = indices[i-1]+Psnum[i-1];
 startindices[i] = indices[i];
 endindices[i-1] = indices[i];
 --endindices[i-1];
 }
 endindices[Pnum-1] = indices[Pnum-2];
 endindices[Pnum-1] += Psnum[Pnum-2];
 --endindices[Pnum-1];
 for (int i = 0; i < Ars; ++i)
 {
 Ax[i] = 0;
 }
 for (int j = 0; j < Acs; ++j)
 {
 if (x[j] != 0)
 {
 //         std::cout << j << " jindex" << std::endl;
 for (int i = 2; i < Pnum; ++i)
 {
 //            std::cout << indices[i] << std::endl;
 Ax[indices[i]] += x[j];
 }
 }
 if ((j+1) % blocklength != 0)
 {
 int k = Pnum-1;
 if (indices[k] < endindices[k])
 {
 ++indices[k];
 }
 else
 {
 int temp = k-1;
 while (temp >= 2)
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
 else
 {
 for (int i = 2; i < Pnum; ++i)
 {
 indices[i] = startindices[i];
 }
 }
 }
 }*/

void mvmult(const int Ars, const int long blocklength, const int long Acs, const std::vector< double> &x, double * Ax, const int * const Psnum, const int Pnum)
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
         //         std::cout << j << std::endl;
         ucn = j-blocklength*blockinum;
         blockilength = blocklength;
         startindex = 0;
         for (int i = 2; i < Pnum; ++i)
         {
            blockilength /= Psnum[i];
            //            std::cout << ucn << " " << blockilength << std::endl;
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
         //         std::cout << j << std::endl;
         ucn = mcind[j]%blocklength;
         blockilength = blocklength;
         startindex = 0;
         
         for (int i = 2; i < Pnum; ++i)
         {
            blockilength /= Psnum[i];
            //            std::cout << ucn << " " << blockilength << std::endl;
            index = ucn / blockilength;
            Ax[startindex + index] += x[j];
            ucn -= index*blockilength;
            startindex += Psnum[i];
         }
      }
   }
}

std::vector<double> dualtotals(const double * const yhat, const int Amrn, const long int blocklength, const int * Psnum, const int &Pnum)
{
   std::vector<double> ytotals(blocklength,0);
   int i = 2; //Should really pass this in. Related to numinprice
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
         //         std::cout << yhat[j] << std::endl;
         index = 0;
         for (int k = 0; k < blockinum; ++k)
         {
            for (int l = 0; l < blocki1length; ++l)
            {
               //               std::cout << "index: " << startindex+index << " yhat: " << yhat[j] << std::endl;
               ytotals[startindex+index] += yhat[j];
               ++index;
            }
            index = blockilength*(k+1);
         }
      }
      ++j;
      //      std::cout << j << std::endl;
      startindex += blocki1length;
      if (j-totalPs == Psnum[i])
      {
         blockinum *= Psnum[i];
         blockilength /= Psnum[i];
         startindex = 0;
         totalPs += Psnum[i];
         ++i;
         if ( i < Pnum)
         {
            blocki1length /= Psnum[i];
         }
      }
   }
   return ytotals;
}

int adjusttotals(std::vector<double> &yhatAmas, const double * const ychange, const int Amrn, const long int blocklength, const int * const Psnum, const int & Pnum)
{
   int i = 2; //Should really pass this in. Related to numinprice
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
         //         std::cout << yhat[j] << std::endl;
         index = 0;
         for (int k = 0; k < blockinum; ++k)
         {
            for (int l = 0; l < blocki1length; ++l)
            {
               //               std::cout << "index: " << startindex+index << " yhat: " << yhat[j] << std::endl;
               yhatAmas[startindex+index] += ychange[j];
               ++index;
            }
            index = blockilength*(k+1);
         }
      }
      ++j;
      //      std::cout << j << std::endl;
      startindex += blocki1length;
      if (j-totalPs == Psnum[i])
      {
         blockinum *= Psnum[i];
         blockilength /= Psnum[i];
         startindex = 0;
         totalPs += Psnum[i];
         ++i;
         if (i < Pnum)
         {
            blocki1length /= Psnum[i];
         }
      }
   }
   return 0;
}
int MeasToTXTFile(std::vector<SuppPt> &Psupp,const int &Psnum, const int &startyear, const int &Pnum, const int &filenum, const int &startindex)
{
   std::ofstream outmeasure;
   std::ostringstream fname2;
   fname2 << "/DenDataOrigMeas" << startyear << "_" << Pnum << "_" << filenum << ".txt";
   outmeasure.open(fname2.str());
   
   for (int k = 0; k < Psnum; ++k)
   {
      outmeasure << Psupp[startindex+k] << '\n';
   }
   outmeasure.close();
   
   return 0;
}

int MeasToCommaFile(std::vector<SuppPt> &Psupp,const int &Psnum, const int &startyear, const int &Pnum, const int &filenum, const int &startindex)
{
   std::ofstream outmeasure;
   outmeasure.precision(10);
   std::ostringstream fname2;
   fname2 << "/DenDataGreedyMeas" << startyear << "_" << Pnum << "_" << filenum << ".txt";
   outmeasure.open(fname2.str());
   outmeasure <<"LOC1,LOC2,MASS"<<'\n';
   
   for (int k = 0; k < Psnum; ++k)
   {
      if (Psupp[startindex+k].mass>1e-12)
      {
         outmeasure <<Psupp[startindex+k].loc1 <<','<<Psupp[startindex+k].loc2<<','<<Psupp[startindex+k].mass<<'\n';
      }
   }
   outmeasure.close();
   
   return 0;
}
