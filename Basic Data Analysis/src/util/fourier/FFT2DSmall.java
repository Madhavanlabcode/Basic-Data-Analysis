package util.fourier;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_2D;
import util.FieldOps;

//This class performs a 2D FFT in the simplest possible way. We use the 1D FFT code on the 
//rows of the matrix. Then, we transpose the matrix and use 1D FFT on the rows again. 
//Then we transpose back.

//On 8 July 2013, I updated this to be able to handle non-powers of 2 data using the methods from Jtransforms.
//The code therefore works in two different ways, depending on if the data size is a power of 2 or not.
//We assume the data is square 2**m X 2**m. ---> now we assume it is merely square
public class FFT2DSmall {
	
	int N;
	int M;
	public double[][][] f, fHat;
	public double[][] falt;

	FFT1D oneDFFT;
	
	boolean isPower2;
	DoubleFFT_2D jtFFT = null;
	
	//does not copy f.
	public FFT2DSmall(double[][] f)
	{	
		N = f.length;
		M = f[0].length;
		this.falt = f;
		fHat = new double[N][M][2];
		isPowerOfTwo();
		if (isPower2)
			oneDFFT = new FFT1D(falt);
		else
			jtFFT = new DoubleFFT_2D(N, M);
	}
	public FFT2DSmall(double[][] f, double[][][] fHat)
	{	
		N = f.length;
		M = f[0].length;
		this.falt = f;
		this.fHat = fHat;
		isPowerOfTwo();
		if (isPower2)
			oneDFFT = new FFT1D(falt);
		else
			jtFFT = new DoubleFFT_2D(N, M);
	}
	
	private void isPowerOfTwo()
	{
		isPower2 = false;
		boolean isN = false, isM = false;
		for (int i = 0; i < 20 && !isN; i++)
			isN = (Math.pow(2, i) == N);
		for (int j = 0; j < 20 && ! isM; j++)
			isM = Math.pow(2, j) == M;
		isPower2 = isN && isM;
		isPower2 = false;
//		System.out.println("It is " + isPower2 + " that the data size is a power of 2.");
	}
	public void setF(double[][] f)
	{
		int oldN = N;
		int oldM = M;
		N = f.length;
		M = f[0].length;
		this.falt = f;
		if (oldN != N || oldM != M)
			fHat = new double[N][M][2];
		isPowerOfTwo();
		if (isPower2)
			oneDFFT = new FFT1D(falt);
		else
			jtFFT = new DoubleFFT_2D(N, M);
	}

	//This copies f by reference
	public FFT2DSmall(double[][][] f)
	{
		N = f.length;
		M = f[0].length;
		this.f = f;
		fHat = new double[N][M][2];
		isPowerOfTwo();
		if (isPower2)
			oneDFFT = new FFT1D(f[0]);
		else
			jtFFT = new DoubleFFT_2D(N, M);
	}
	//This copies f by reference
	public FFT2DSmall(double[][][] fHat, boolean centered)
	{
		N = fHat.length;
		M = fHat[0].length;
		if(!centered) this.fHat = fHat;
		else {this.fHat = new double [N][M][2]; FieldOps.shift(fHat, this.fHat);}
		f = new double[N][M][2];
		isPowerOfTwo();
		if (isPower2)
			oneDFFT = new FFT1D(falt);
		else
			jtFFT = new DoubleFFT_2D(N, M);
	}
	//assuming the size remains the same.
	public void reset(double[][][] fHat, boolean centered)
	{
		N = fHat.length;
		M = f[0].length;
		if(!centered) this.fHat = fHat;
		else {this.fHat = new double [N][M][2]; FieldOps.shift(fHat, this.fHat);}
//		f = new double[N][N][2];
//		oneDFFT = new FFT1D(this.fHat[0]);
	}
	public void reset(double[][] f)
	{
		this.falt = f;
		this.f = null;
	}
	
	public void doFFT()
	{
		
		if (isPower2){
			//First, FFT the rows:
			if (f != null){
				for (int i = 0; i < N; i++)
				{
					oneDFFT.doFFT(f[i], fHat[i]);
				}
				//Then the columns
				for (int i = 0; i < N; i++)
				{
					oneDFFT.doFFT(fHat, fHat, i);
				}
			}
			else{
				for (int i = 0; i < N; i++)
				{
					oneDFFT.doFFTRe(falt[i], fHat[i]);
				}
				//Then the columns
				for (int i = 0; i < N; i++)
				{
					oneDFFT.doFFT(fHat, fHat, i);
				}
			}
		}
		else
		{
			//invoke methods of JTransform;
			double[] temp = new double [N*M*2];
			if (f != null)
				for (int i = 0; i < N; i++)
					for (int j = 0; j < M; j++)
					{
						temp[2*(i*M + j)] = f[i][j][0];
						temp[2*(i*M + j) + 1] = f[i][j][1];
					}
			else
				for (int i = 0; i < N; i++)
					for (int j = 0; j < M; j++)
					{
						temp[2*(i*M + j)] = falt[i][j];
						temp[2*(i*M + j) + 1] = 0;
					}
				
			jtFFT.complexForward(temp);
			for (int i = 0; i < N; i++)
				for (int j = 0; j < M; j++){
					fHat[i][j][0] = temp[2*(i*M + j)];
					fHat[i][j][1] = temp[2*(i*M + j) + 1];
				}
		}
	}
	
	//we start from fHat, and fft back.
	//The non-square FFT only works if !isPower2;
	public void doIFFT()
	{
		if (f == null)
			f = new double[N][M][2];
		if (isPower2){
			//Take the star:
			FieldOps.star(fHat, fHat);
			for (int i = 0; i < N; i++)
			{
				oneDFFT.doFFT(fHat[i], f[i]);
			}
			
			//Then the columns
			for (int i = 0; i < N; i++)
			{
				oneDFFT.doFFT(f, f, i);
			}
			FieldOps.star(f, f);
			FieldOps.star(fHat, fHat);
			for (int i = 0; i < N; i++)
				for(int j = 0; j < N; j++)
				{
					f[i][j][0] /= N*N;
					f[i][j][1] /= N*N;
				}
		}
		else
		{
			//invoke methods of JTransform;
			double[] temp = new double [N*M*2];
			for (int i = 0; i < N; i++)
				for (int j = 0; j < M; j++)
				{
					temp[2*(i*M + j)] = fHat[i][j][0];
					temp[2*(i*M + j) + 1] = fHat[i][j][1];
				}
			jtFFT.complexInverse(temp, true);
			for (int i = 0; i < N; i++)
				for (int j = 0; j < M; j++){
					f[i][j][0] = temp[2*(i*M + j)];
					f[i][j][1] = temp[2*(i*M + j) + 1];
				}
		}
//		//then, transpose the temp1;
//		rotate(temp1, temp2);
//		//FFT on the transpose.
//		for (int i = 0; i < N; i++)
//		{
//			oneDFFT.doFFT(temp1[i], fHat[i]);
//		}
//		//Transpose back for the answer
//		rotateBack(fHat, temp2);
		
		//place the "centered" values in fHat2
	}
	
	public double getFHat2Re(int i, int j)
	{
		return fHat[(i + N/2) % N][(j + M/2) % M][0];
	}
	public double getFHat2Im(int i, int j)
	{
		return fHat[(i + N/2) % N][(j + M/2) % M][1];
	}
	public double getFHat2Mag(int i, int j)
	{
		return Math.sqrt(fHat[(i + N/2) % N][(j + M/2) % M][1]*fHat[(i + N/2) % N][(j + M/2) % M][1] + fHat[(i + N/2) % N][(j + M/2) % M][0]*fHat[(i + N/2) % N][(j + M/2) % M][0]);
	}
	
	//should not be used to save memory.
	public double[][][] fHat2()
	{
		double[][][] fHat2 = new double [N][M][2];
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
			{
				fHat2[i][j] = fHat[(i + N/2) % N][(j + M/2) % M];
			}
		return fHat2;
	}
	public void fHat2(double[][][] target)
	{
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
			{
				target[i][j] = fHat[(i + N/2) % N][(j + M/2) % M];
			}
	}
	void transpose(double[][][] source, double[][][] target)
	{
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M ; j++)
			{
				target[i][j][0] = source[j][i][0]; target[i][j][1] = source[j][i][1];
			}
	}
	
	//this is a transpose with a vertical flip
	void rotate(double[][][] sourceTarget, double[][][] temp)
	{
		transpose(sourceTarget, temp);
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
			{
				sourceTarget[i][j][0] = temp[i][M-1 - j][0];
				sourceTarget[i][j][1] = temp[i][M-1 - j][1];
			}
	}
	void rotateBack(double[][][] sourceTarget, double[][][] temp)
	{
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
			{
				temp[i][j][0] = sourceTarget[i][M-1 - j][0];
				temp[i][j][1] = sourceTarget[i][M-1 - j][1];
			}
		transpose(temp, sourceTarget);
	}
}
