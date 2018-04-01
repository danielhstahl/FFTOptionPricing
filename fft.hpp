/**
Note that neither fft or ifft include the division by N.  The division by
N needs to be done outside this function.  
*/
auto genericfft(CArray&& x, Complex&& phiT, int N){
  
	unsigned int k=N;
  unsigned int n;
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
			T *= phiT;
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
  return std::move(x);
}


auto ifft(CArray&& x){
  unsigned int N = x.size();
  double thetaT = M_PI / N;
  return std::move(genericfft(std::move(x), Complex(cos(thetaT), sin(thetaT), N)));
}

auto fft(CArray&& x){
  unsigned int N = x.size();
  double thetaT = M_PI / N;
  return std::move(genericfft(std::move(x), Complex(cos(thetaT), sin(-thetaT), N)));
}