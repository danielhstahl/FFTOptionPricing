CarrMadan::CarrMadan(int numSteps_){
  setValues(numSteps_, .25, 1.5);
}
CarrMadan::CarrMadan(int numSteps_, double ada_, double alpha_){
  setValues(numSteps_, ada_, alpha_);
}
void CarrMadan::setValues(int numSteps_, double ada_, double alpha_){
  numSteps=numSteps_;
  ada=ada_;
  alpha=alpha_;
  lambda=(2*M_PI/numSteps)/ada;
  b=.5*numSteps*lambda;
}
std::vector<priceAndStrike> CarrMadan::OptionPrice(
  double t, //time to maturity
  double discount,//e^{-rt} in BS
  auto& augCF, //augmented CF
  auto& CF //function characteristic function from CF.h
){
  std::vector<Complex> cmpl(numSteps);
  int everyOther=1;
  cmpl[0]=discount*augCF(0, t, alpha, CF);
  for(int i=1; i<numSteps; ++i){
    cmpl[i]=discount*augCF(i*ada, t, alpha, CF)*exp(Complex(0, b*i*ada))*(3+everyOther);
    everyOther=everyOther*(-1);
  }
  fft(cmpl);
  std::vector<priceAndStrike> ck(numSteps);
  for(int i=0; i<numSteps; ++i){
    ck[i].strike=-b+lambda*i;
    ck[i].price=cmpl[i].getReal()*exp(-alpha*ck[i].strike)*ada/(M_PI*3.0);
  }
  return ck;
}
