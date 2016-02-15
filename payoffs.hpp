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
