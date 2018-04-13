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
  x=genericfft(std::move(x), Complex(cos(thetaT), sin(thetaT)), N);
  return std::move(x);
}

auto fft(CArray&& x){
  unsigned int N = x.size();
  double thetaT = M_PI / N;
  x=genericfft(std::move(x), Complex(cos(thetaT), sin(-thetaT)), N);
  return std::move(x);
}
constexpr auto cmp=std::complex<double>(1.0, 0.0);
constexpr auto cmpi=std::complex<double>(0.0, 1.0);
template<typename FN, typename Method, typename EditOutput>
auto genericIntegration(double outputMin, double inputMin, double inputMax, FN&& fn, int N, Method&& f, EditOutput&& fnOutput){
	const double integrateRange=inputMax-inputMin;
	const double inputDx=integrateRange/(N-1); //this is ada in carr madan
	const double outputMax=outputMin+(N-1)*2.0*M_PI/integrateRange;
	const double du=(outputMax-outputMin)/N;
	return futilities::for_each_parallel(
		f(futilities::for_each_parallel(0, N, [&](const auto& index){
			const double currInput=inputMin+inputDx*index;		
			const auto fnAtX=fn(currInput)*inputDx*exp(outputMin*inputDx*index*cmpi)/3.0;
			//simpsons rule
			if(index==0||index==N-1){
				return fnAtX;
			} 
			else{
				return index%2==0?2.0*fnAtX:4.0*fnAtX;
			}
		})),
		[&](const auto& v, const auto& index){
			const double u=outputMin+index*du;
			return fnOutput(u, v*exp(u*inputMin*cmpi));
		}
	);
}
template<typename FN, typename Method>
auto genericIntegration(double outputMin, double inputMin, double inputMax, FN&& fn, int N, Method&& f){
	return genericIntegration(outputMin, inputMin, inputMax, fn, N, f, [](const auto& u, const auto& v){
		return v;
	});
}

template<typename FN>
auto integrateFFT(double xMin, double uMin, double uMax, FN&& fn, int N){
	return genericIntegration(xMin, uMin, uMax, fn, N, fft);
}
template<typename FN,  typename EditOutput>
auto integrateFFT(double xMin, double uMin, double uMax, FN&& fn, int N, EditOutput&& fnOutput){
	return genericIntegration(xMin, uMin, uMax, fn, N, fft, fnOutput);
}
template<typename FN, typename Method>
auto integrateIFFT(double uMin, double xMin, double xMax, FN&& fn, int N){
	return genericIntegration(uMin, xMin, xMax, fn, N, ifft);
}
template<typename FN, typename EditOutput>
auto integrateIFFT(double uMin, double xMin, double xMax, FN&& fn, int N, EditOutput&& fnOutput){
	return genericIntegration(uMin, xMin, xMax, fn, N, ifft, fnOutput);
}
