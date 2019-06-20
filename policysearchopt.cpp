#include <Rcpp.h>
# define PIC              3.1415926535897932384626433832795028841971693993751
# define SQRTTWO           1.4142135623730950488016887242096980785696718753769
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector policysearch(DataFrame data, int ntheta) {
  NumericVector x1 = data["x.1"];
  NumericVector x2 = data["x.2"];
  NumericVector wts0 = data["wts0"];
  NumericVector wts1 = data["wts1"];
  NumericVector out(2);
  unsigned int n = wts0.length();
  float theta;
  float suml,sumr = 0;
  float theta_inc = PIC * 2. / ntheta;
  NumericVector vals;
  std::vector<std::tuple<float,float,float>> zipped(n);
  float totwts0 = sum(wts0);
  float totwts1 = sum(wts1);
  out[0] = 0.0;
  out[1] = totwts1>totwts0 ? std::numeric_limits<double>::max() : -std::numeric_limits<double>::max();
  float best = totwts1>totwts0 ? totwts1 : totwts0;
  for (int it=0; it<=ntheta; it++){
    theta = it==ntheta ? 5.0*PIC/4.0 : theta_inc*it;
    vals  = cos(theta)*x1 + sin(theta)*x2;
    for(size_t i=0; i<n; ++i){std::get<0>(zipped[i])=vals[i];std::get<1>(zipped[i])=wts0[i];std::get<2>(zipped[i])=wts1[i];}
    std::sort(zipped.begin(), zipped.end());
    suml = 0.0;
    sumr = totwts0;
    for(size_t ib=0; ib<n-1; ib++){
      sumr -= std::get<1>(zipped[ib]);
      suml += std::get<2>(zipped[ib]);
      if(suml+sumr > best){
        best = suml+sumr;
        out[0] = theta;
        out[1] = (std::get<0>(zipped[ib]) + std::get<0>(zipped[ib+1])) / 2.0;
      }
    }
  }
  return out;
}
