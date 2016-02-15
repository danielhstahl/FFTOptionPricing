double Call(double s0, double y, double k){
  double st=s0*exp(y);
  if(st>k){
    return st-k;
  }
  else{
    return 0;
  }
}
double Put(double s0, double y, double k){
  double st=s0*exp(y);
  if(st<k){
    return k-st;
  }
  else{
    return 0;
  }
}
Complex Call(double v, double t, double alpha, auto& cf){ //used for Carr-Madan approach
  Complex u=v-Complex(0, alpha+1);
  return cf(u, t)/(alpha*alpha+alpha-v*v+Complex(0, (2*alpha+1)*v));
}
