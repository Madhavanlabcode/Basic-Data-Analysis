package util.fourier;

public class FFT2D_Subset {

	double[][][] data;
	
	int N, M;
	int s; //we assume that the contour region is a square.
	
	double[][][] fHat;
	public double[][][] fHat2;
	
	FFT1D fft;
	
	//This takes the fft of a subset of the data from i = rx to rx+sx; same for y.
	//sx and sy must be powers of 2.

	
	//real
	public FFT2D_Subset(double[][] data, int size)
	{
		N = data.length; M = data[0].length;
		this.data = new double[N][N][2];
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
			{
				this.data[i][j][0] = data[i][j];
				this.data[i][j][1] = 0;
			}
		s = size;
		fHat = new double[s][s][2];
		fHat2 = new double[s][s][2];
		
		fft = new FFT1D(s);
	}
	
	public void doFFT(int rx, int ry)
	{
		//First, FFT the rows:
		for (int i = 0; i < s; i++)
		{
			fft.doFFT(data[i + rx], fHat[i], ry);
		}
		//Then the columns
		for (int i = 0; i < s; i++)
		{
			fft.doFFT(fHat, fHat, i);
		}

		//place the "centered" values in fHat2
		setFHat2();

	}

	void setFHat2()
	{
		for (int i = 0; i < s; i++)
			for (int j = 0; j < s; j++)
			{
				fHat2[i][j][0] = fHat[(i + s/2) % s][(j + s/2) % s][0];
				fHat2[i][j][1] = fHat[(i + s/2) % s][(j + s/2) % s][1];
			}
		
	}
	
}
	
