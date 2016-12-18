FSTS::FSTS(int numSteps_, double xmin_, double xmax_){
  numSteps=numSteps_;
  xmin=xmin_;
  xmax=xmax_;
  dx=(xmax-xmin)/(double)(numSteps-1);
  vmax=M_PI/dx;
  du=2*vmax/numSteps;

}
std::vector<priceAndUnderlying> FSTS::OptionPrice(
  double t, //time to maturity
  double discount,//e^{-rt} in BS
  auto& payoff,//function payoff from payoffs.h
  auto& CF //function characteristic function from CF.h
){
  std::vector<Complex> payoffV(numSteps);
  double vmin=-vmax+du;
  for(int i=0; i<numSteps;i++){
    payoffV[i]=payoff(xmin+i*dx)*exp(Complex(0, -vmin*i*dx));
  }
  fft(payoffV);
  for(int i=0; i<numSteps; i++){
    Complex u=Complex(vmin+i*du, 0);
    payoffV[i]=payoffV[i]*CF(u, t);//*exp(Complex(0, -vmin*i*dx));
  }
  ifft(payoffV);
  std::vector<priceAndUnderlying> price(numSteps);
  for(int i=0;i<numSteps;++i){
    price[i].price=discount*(payoffV[i]*exp(Complex(0, vmin*i*dx))).getReal()/numSteps; //*cos(vmin*i*dx)
    price[i].underlying=xmin+i*dx;
  }
  return price;
}
