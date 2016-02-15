
// Better optimized but less intuitive
void ifft(std::vector<Complex>& x){
	// DFT
	//std::reverse(x.begin(), x.end()); //reverse is O(n)
	unsigned int N = x.size(), k = N, n;
	double thetaT = M_PI / N;
	Complex phiT = Complex(cos(thetaT), sin(thetaT));
  Complex T;
	while (k > 1)	{
		n = k;
		k >>= 1;
		phiT = phiT * phiT;
		T = Complex(1.0, 0);
		for (unsigned int l = 0; l < k; l++)
		{
			for (unsigned int a = l; a < N; a += n)
			{
				unsigned int b = a + k;
				Complex t = x[a] - x[b];
				x[a] = x[a]+x[b];
				x[b] = t * T;
			}
			T = T*phiT;
		}
	}
	// Decimate
	unsigned int m = (unsigned int)log2(N);
	for (unsigned int a = 0; a < N; a++)
	{
		unsigned int b = a;
		// Reverse bits
		b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
		b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
		b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
		b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
		b = ((b >> 16) | (b << 16)) >> (32 - m);
		if (b > a)
		{
			Complex t = x[a];
			x[a] = x[b];
			x[b] = t;
		}
	}
}

// fft (in-place)
void fft(std::vector<Complex>& x){
	//std::reverse(x.begin(), x.end()); //reverse is O(n)
  unsigned int N = x.size(), k = N, n;
  double thetaT = M_PI / N;
  Complex phiT = Complex(cos(thetaT), sin(-thetaT)); //conjugate going in
  Complex T;
  while (k > 1)	{
    n = k;
    k >>= 1;
    phiT = phiT * phiT;
    T = Complex(1.0, 0);
    for (unsigned int l = 0; l < k; l++)
    {
      for (unsigned int a = l; a < N; a += n)
      {
        unsigned int b = a + k;
        Complex t = x[a] - x[b];
        x[a] = x[a]+x[b];
        x[b] = t * T;
      }
      T = T*phiT;
    }
  }
  // Decimate
  unsigned int m = (unsigned int)log2(N);
  for (unsigned int a = 0; a < N; a++)
  {
    unsigned int b = a;
    // Reverse bits
    b = (((b & 0xaaaaaaaa) >> 1) | ((b & 0x55555555) << 1));
    b = (((b & 0xcccccccc) >> 2) | ((b & 0x33333333) << 2));
    b = (((b & 0xf0f0f0f0) >> 4) | ((b & 0x0f0f0f0f) << 4));
    b = (((b & 0xff00ff00) >> 8) | ((b & 0x00ff00ff) << 8));
    b = ((b >> 16) | (b << 16)) >> (32 - m);
    if (b > a)
    {
      Complex t = x[a];
      x[a] = x[b];
      x[b] = t;
    }
  }
}

/*void fft(std::vector<Complex>& x){
  fft(1, x);
}
void ifft(std::vector<Complex>& x){
  fft(-1, x);
}
void fft(int dir,std::vector<Complex>& x){
   long n,i,i1,j,k,i2,l,l1,l2;
   double c1,c2;//,tx,ty,t1,t2,u1,u2,z;
    n=x.size();
   // Do the bit reversal
   i2 = n >> 1;
   j = 0;
   for (i=0;i<n-1;i++) {
      if (i < j) {
        std::swap(x[i],x[j]);
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }
	 unsigned int m = (unsigned int)log2(n);


   l2 = 1;
   Complex c(-1, 0);
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      Complex u(1, 0);
      //u1 = 1.0;
      //u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<n;i+=l2) {
            i1 = i + l1;

            Complex t=x[i1]*u;
            x[i1]=x[i]-t;
            x[i]=x[i]+t;
         }
         u=u*c;
      }
      c2 = sqrt((1.0 - c.getReal()) / 2.0);
      if (dir == 1){
         c2 = -c2;
       }
      c1 = sqrt((1.0 + c.getReal()) / 2.0);
      c.setValues(c1, c2);
      //c=Complex(c1, c2);
   }

}
*/
