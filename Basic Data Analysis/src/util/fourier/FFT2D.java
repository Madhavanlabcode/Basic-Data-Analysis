package util.fourier;

import util.FieldOps;

//This class performs a 2D FFT in the simplest possible way. We use the 1D FFT code on the 
//rows of the matrix. Then, we transpose the matrix and use 1D FFT on the rows again. 
//Then we transpose back.
//We assume the data is square 2**m X 2**m.
public class FFT2D {

	int N;
	public double[][][] f, temp1, temp2, fHat, fHat2;
	
	FFT1D oneDFFT;
	
	//real fft
	public FFT2D(double[][] f)
	{
		N = f.length;
		this.f = new double[N][N][2];
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				this.f[i][j][0] = f[i][j]; this.f[i][j][1] = 0;
			}
		temp1 = new double[N][N][2];
		fHat = new double[N][N][2];
		fHat2 = new double[N][N][2];
		oneDFFT = new FFT1D(this.f[0]);
	}

	//This copies f by reference
	public FFT2D(double[][][] f)
	{
		N = f.length;
		this.f = f;
		temp1 = new double[N][N][2];
		fHat = new double[N][N][2];
		fHat2 = new double[N][N][2];
		oneDFFT = new FFT1D(this.f[0]);
	}
	
	public void doFFT()
	{
		//First, FFT the rows:
		for (int i = 0; i < N; i++)
		{
			oneDFFT.doFFT(f[i], temp1[i]);
		}
		//Then the columns
		for (int i = 0; i < N; i++)
		{
			oneDFFT.doFFT(temp1, fHat, i);
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
		setFHat2();
	}
	
	//we start from fHat, and fft back.
	public void doIFFT()
	{
		//Take the star:
		FieldOps.star(fHat, temp1);
		for (int i = 0; i < N; i++)
		{
			oneDFFT.doFFT(temp1[i], temp1[i]);
		}
		
		//Then the columns
		for (int i = 0; i < N; i++)
		{
			oneDFFT.doFFT(temp1, f, i);
		}
		FieldOps.star(f, f);
		for (int i = 0; i < N; i++)
			for(int j = 0; j < N; j++)
			{
				f[i][j][0] /= N*N;
				f[i][j][1] /= N*N;
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
		setFHat2();
	}
	
	void setFHat2()
	{
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				fHat2[i][j][0] = fHat[(i + N/2) % N][(j + N/2) % N][0];
				fHat2[i][j][1] = fHat[(i + N/2) % N][(j + N/2) % N][1];
			}
		
	}
	void transpose(double[][][] source, double[][][] target)
	{
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N ; j++)
			{
				target[i][j][0] = source[j][i][0]; target[i][j][1] = source[j][i][1];
			}
	}
	
	//this is a transpose with a vertical flip
	void rotate(double[][][] sourceTarget, double[][][] temp)
	{
		transpose(sourceTarget, temp);
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				sourceTarget[i][j][0] = temp[i][N-1 - j][0];
				sourceTarget[i][j][1] = temp[i][N-1 - j][1];
			}
	}
	void rotateBack(double[][][] sourceTarget, double[][][] temp)
	{
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
			{
				temp[i][j][0] = sourceTarget[i][N-1 - j][0];
				temp[i][j][1] = sourceTarget[i][N-1 - j][1];
			}
		transpose(temp, sourceTarget);
	}
}
