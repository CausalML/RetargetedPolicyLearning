#include <Rcpp.h>
using namespace Rcpp;
#include "gurobi_c++.h"

// [[Rcpp::export]]
NumericVector optlinearpol(DataFrame data, int n, int d, int timelimit, int verbose) {
  NumericVector wts1 = data["wts1"];
  NumericVector wts2 = data["wts2"];
  NumericVector wts3 = data["wts3"];
  
  GRBEnv env = GRBEnv(true);
  env.start();
  
  GRBModel model = GRBModel(env);
  model.getEnv().set(GRB_DoubleParam_TimeLimit, timelimit);
  model.getEnv().set(GRB_IntParam_Threads, 1);
  model.getEnv().set(GRB_IntParam_LogToConsole, verbose);
  
  std::vector<GRBVar> b1(d);
  std::vector<GRBVar> b2(d);
  std::vector<GRBVar> b3(d);
  for(int j = 0; j < d; j++) {
    b1[j] = model.addVar(-1.,1.,0.,GRB_CONTINUOUS,"b1_"+std::to_string(j));
    b2[j] = model.addVar(-1.,1.,0.,GRB_CONTINUOUS,"b2_"+std::to_string(j));
    b3[j] = model.addVar(-1.,1.,0.,GRB_CONTINUOUS,"b3_"+std::to_string(j));
  }
  GRBVar a1   = model.addVar(-GRB_INFINITY,GRB_INFINITY,0.,GRB_CONTINUOUS,"a1");
  GRBVar a2   = model.addVar(-GRB_INFINITY,GRB_INFINITY,0.,GRB_CONTINUOUS,"a2");
  GRBVar a3   = model.addVar(-GRB_INFINITY,GRB_INFINITY,0.,GRB_CONTINUOUS,"a3");
  std::vector<GRBVar> z1(n);
  std::vector<GRBVar> z2(n);
  std::vector<GRBVar> z3(n);
  for(int i = 0; i < n; i++) {
    z1[i] = model.addVar(0,1,-wts1[i],GRB_BINARY,"z1_"+std::to_string(i));
    z2[i] = model.addVar(0,1,-wts2[i],GRB_BINARY,"z2_"+std::to_string(i));
    z3[i] = model.addVar(0,1,-wts3[i],GRB_BINARY,"z3_"+std::to_string(i));
  }

  std::vector<NumericVector> x(d);
  for(int j = 0; j < d; j++) {
    x[j] = data["x."+std::to_string(j+1)];
  }
  
  GRBLinExpr bx1,bx2,bx3;
  for(int i = 0; i < n; i++) {
    bx1 = GRBLinExpr(a1,1.0);
    bx2 = GRBLinExpr(a2,1.0);
    bx3 = GRBLinExpr(a3,1.0);
    for(int j = 0; j < d; j++) {
      bx1+=GRBLinExpr(b1[j],x[j][i]);
      bx2+=GRBLinExpr(b2[j],x[j][i]);
      bx3+=GRBLinExpr(b3[j],x[j][i]);
    }
    model.addGenConstrIndicator(z1[i],1, bx1-bx2,GRB_GREATER_EQUAL,0.001);
    model.addGenConstrIndicator(z1[i],1, bx1-bx3,GRB_GREATER_EQUAL,0.001);
    model.addGenConstrIndicator(z2[i],1, bx2-bx1,GRB_GREATER_EQUAL,0.001);
    model.addGenConstrIndicator(z2[i],1, bx2-bx3,GRB_GREATER_EQUAL,0.001);
    model.addGenConstrIndicator(z3[i],1, bx3-bx1,GRB_GREATER_EQUAL,0.001);
    model.addGenConstrIndicator(z3[i],1, bx3-bx2,GRB_GREATER_EQUAL,0.001);
    model.addConstr(z1[i]+z2[i]+z3[i],GRB_EQUAL,1);
  }
  
  model.optimize();
  
  NumericVector out(3*(d+1)+1);
  
  if(model.get(GRB_IntAttr_SolCount) > 0) {
    out[0] = a1.get(GRB_DoubleAttr_X);
    for(int j = 0; j < d; j++) {
      out[1+j] = b1[j].get(GRB_DoubleAttr_X);
    }
    out[1+d] = a2.get(GRB_DoubleAttr_X);
    for(int j = 0; j < d; j++) {
      out[2+d+j] = b2[j].get(GRB_DoubleAttr_X);
    }
    out[2+2*d] = a3.get(GRB_DoubleAttr_X);
    for(int j = 0; j < d; j++) {
      out[3+2*d+j] = b3[j].get(GRB_DoubleAttr_X);
    }
    out[3*(d+1)] = 1; // success code
  } else {
    for(int j = 0; j < 3*(d+1); j++) {
      out[j] = 0;
    }
    out[3*(d+1)] = -1; // fail code
  }
  
  return out;
}