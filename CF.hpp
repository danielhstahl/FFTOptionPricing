Complex BMWithDrift(
  Complex& u, //u
  double drift, //drift (eg, r-sigma*sigma*.5)
  double vol,//volatility (eg, sigma)
  double t //time to maturity
){
  //Complex ui(0, u);
  u=u*Complex(0, 1);
  return exp((u*drift+u*u*vol*vol*.5)*t);
}
