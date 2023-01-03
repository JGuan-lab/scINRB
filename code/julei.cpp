#include <Rcpp.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <set>
#include <string>
#include <algorithm>

using namespace std;

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

double g_nodes_num;
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). 
NumericVector outereva(CharacterVector cluster,CharacterVector realcluster,int length) {
  int i=0;
  double a=0;
  double b=0;
  double c=0;
  double d=0;
  NumericVector result(3);
  while (i<length) {
    int j=i+1;
    while (j<length) {
      if(cluster[i]==cluster[j]){
        if(realcluster[i]==realcluster[j])
          a=a+1;
        else b=b+1;
      }else{
        if(realcluster[i]==realcluster[j])
          c=c+1;
        else d=d+1;
      }
      j=j+1;
    }
    i=i+1;
  }
  double jc=a/(a+b+c);
  double fmi=sqrt((a/(a+b))*(a/(a+c)));
  double ri=2*(a+d)/(length*(length-1));
  result[0]=jc;
  result[1]=fmi;
  result[2]=ri;
  
  return result;
}

// [[Rcpp::export]]
List innereva(NumericMatrix distE,List cluster){
  int clustertype=cluster.size();
  NumericVector avg(clustertype);
  NumericVector diam(clustertype);
  int i=0;
  for(i=0;i<clustertype;i++){
    NumericVector test=cluster[i];
    double indexlength=test.length();
    int j=0;
    int k=0;
    double sum=0;
    double max=0;
    for(j=0;j<indexlength;j++){
      for(k=j+1;k<indexlength;k++){
        sum=sum+distE(test[j]-1,test[k]-1);
        if(max<distE(test[j]-1,test[k]-1))
          max=distE(test[j]-1,test[k]-1);
      }
    }
    avg[i]=2*sum/(indexlength*(indexlength-1));
    diam[i]=max;
  }
  
  NumericMatrix dmin(clustertype,clustertype);
  for(i=0;i<clustertype;i++){
    int j=0;
    for(j=i;j<clustertype;j++){
      if(i==j)
        dmin(i,j)=0;
      else{
        int k=0;
        int h=0;
        double min=100;
        NumericVector c1=cluster[i];
        NumericVector c2=cluster[j];
        for(k=0;k<c1.length();k++){
          for(h=0;h<c2.length();h++){
            if(min>distE(c1[k]-1,c2[h]-1))
              min=distE(c1[k]-1,c2[h]-1);
          }
        }
        dmin(i,j)=min;
        dmin(j,i)=min;
      }
    }
  }
  return List::create(Named("avg")=avg,Named("diam")=diam,Named("dmin")=dmin);
}


double juleilog2(double x)
{
  return log(x) / log(2);
}
double juleih(double x)
{
  if (x > 0)
    return -1 * x * juleilog2(x);
  else
    return 0;
}
//Hi(Xi)
double juleiHi(vector<int> & Xi)
{
  double p1 = Xi.size() / g_nodes_num;
  double p0 = 1 - p1;
  return juleih(p0) + juleih(p1);
}
//H(X)

double juleiH(List X)
{
  double res = 0;
  
  for (int i = 0; i < X.size(); i++)
  {
    std::vector<int> xi=X[i];
    res += juleiHi(xi);
  }
  return res;
}

vector<int> intersection(vector<int>& a, vector<int>& b)
{
  sort(a.begin(), a.end());
  sort(b.begin(), b.end());
  
  vector<int> res(max(a.size(), b.size()));
  auto iter = std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), res.begin());
  res.resize(iter - res.begin());
  
  return res;
}

vector<int> difference(vector<int>& a, vector<int>& b)
{
  sort(a.begin(), a.end());
  sort(b.begin(), b.end());
  
  vector<int> res(max(a.size(), b.size()));
  auto iter = std::set_difference(a.begin(), a.end(), b.begin(), b.end(), res.begin());
  res.resize(iter - res.begin());
  
  return res;
}
double H_Xi_joint_Yj(vector<int> & Xi, vector<int> & Yj)
{
  double P11 = intersection(Xi, Yj).size() / g_nodes_num;
  double P10 = difference(Xi, Yj).size() / g_nodes_num;
  double P01 = difference(Yj, Xi).size() / g_nodes_num;
  double P00 = 1 - P11 - P10 - P01;
  
  if (juleih(P11) + juleih(P00) >= juleih(P01) + juleih(P10))
    return juleih(P11) + juleih(P10) + juleih(P01) + juleih(P00);
  else
    return juleiHi(Xi) + juleiHi(Yj);
}
double H_Xi_given_Yj(vector<int> & Xi, vector<int> & Yj)
{
  return H_Xi_joint_Yj(Xi, Yj) - juleiHi(Yj);
}

double H_Xi_given_Y(vector<int> & Xi, List Y)
{
  std::vector<int> y0=Y[0];
  double res = H_Xi_given_Yj(Xi, y0);
  for (int i = 1; i < Y.size(); i++)
  {
    std::vector<int> yi=Y[i];
    res = min(res, H_Xi_given_Yj(Xi, yi));
  }
  return res;
}

double H_Xi_given_Y_norm(vector<int> & Xi, List Y)
{
  return H_Xi_given_Y(Xi, Y) / juleiHi(Xi);
}
double H_X_given_Y(List X, List Y)
{
  double res = 0;
  for (int i = 0; i < X.size(); i++)
  {
    std::vector<int> xi=X[i];
    res += H_Xi_given_Y(xi, Y);
  }
  
  return res;
}

double H_X_given_Y_norm(List X, List Y)
{
  double res = 0;
  for (int i = 0; i < X.size(); i++)
  {
    std::vector<int> xi=X[i];
    res += H_Xi_given_Y_norm(xi, Y);
  }
  
  return res / X.size();
}

int getNodesNum(List X, List Y)
{
  set<int> s;
  for (int i = 0; i < X.size(); i++)
  {
    std::vector<int> xi=X[i];
    for (int j = 0; j < xi.size(); j++)
      s.insert(xi[j]);
  }
  for (int i = 0; i < Y.size(); i++)
  {
    std::vector<int> yi=Y[i];
    for (int j = 0; j < yi.size(); j++)
      s.insert(yi[j]);
  }
  return s.size();
}

// [[Rcpp::export]]

double NMI(List X, List Y)
{
  if (X.size() == 0 || Y.size() == 0)
    return 0;
  g_nodes_num = getNodesNum(X, Y);
  
  return 1 - 0.5 * (H_X_given_Y_norm(X, Y) + H_X_given_Y_norm(Y, X));
}


//I(X:Y)
double juleiI(List X, List Y)
{
  double x=0.5 * (juleiH(X) + juleiH(Y) - H_X_given_Y(X, Y) - H_X_given_Y(Y, X));
  return x;
}

// [[Rcpp::export]]
double NMI_max(List X, List Y)
{
  if (X.size() == 0 || Y.size() == 0)
    return 0;
  g_nodes_num = getNodesNum(X, Y);
  return juleiI(X,Y) / max(juleiH(X), juleiH(Y));
}
// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
*/
