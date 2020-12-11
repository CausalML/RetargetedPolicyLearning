// find best policy I[x>theta]

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector policysearch(DataFrame data) {
  NumericVector x = data["x"];
  NumericVector wts0 = data["wts0"];
  NumericVector wts1 = data["wts1"];
  NumericVector out(1);
  unsigned int n = wts0.length();
  float suml,sumr = 0;
  std::vector<std::tuple<float,float,float>> zipped(n);
  float totwts0 = sum(wts0);
  float totwts1 = sum(wts1);
  out[0] = 0.0;
  out[1] = totwts1>totwts0 ? std::numeric_limits<double>::max() : -std::numeric_limits<double>::max();
  float best = totwts1>totwts0 ? totwts1 : totwts0;
  for(size_t i=0; i<n; ++i){std::get<0>(zipped[i])=x[i];std::get<1>(zipped[i])=wts0[i];std::get<2>(zipped[i])=wts1[i];}
  std::sort(zipped.begin(), zipped.end());
  suml = 0.0;
  sumr = totwts0;
  for(size_t ib=0; ib<n-1; ib++){
    sumr -= std::get<1>(zipped[ib]);
    suml += std::get<2>(zipped[ib]);
    if(suml+sumr > best){
      best = suml+sumr;
      out[0] = (std::get<0>(zipped[ib]) + std::get<0>(zipped[ib+1])) / 2.0;
    }
  }
  return out;
}
