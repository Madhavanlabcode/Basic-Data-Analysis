package util.fourier;

import util.ArrayOps;

public class FFT1D {

	int N;
	//This is the complex data. the first index is i, the second is 0 for real, 1 for imag.
	double[][] f;
	
	//the FT of f.
	public double[][] fHat;
	
	//the values of k which are possible, in the usual increasing order. This is with array spacing = unit length.
	double[] ki;
	double[][] phasei;
	
	//In order to do this without allocating new memory we need a block of complex vectors, this size
	//of which is log_2(N);
	int logN;
	
	double[][][] storage;
	
	public FFT1D(double[][] f)
	{
		this.f = f;
		N = f.length;
		
		ki = new double[N];
		phasei = new double[N][2];
		fHat = new double[N][2];
		
		for (int i = 0; i < N; i++){
			ki[i] = -2*Math.PI*i/N;
			phasei[i][0] = Math.cos(ki[i]);
			phasei[i][1] = Math.sin(ki[i]);
		}
		
		logN = 0;
		int n = N;
		while(n > 1)
		{
			n/= 2; logN++;
		}
		storage = new double[logN][N][2];
	}
	
	public FFT1D(int m)
	{
		N = m;
		
		f = new double[N][2];
		ki = new double[N];
		phasei = new double[N][2];
		fHat = new double[N][2];
		
		for (int i = 0; i < N; i++){
			ki[i] = -2*Math.PI*i/N;
			phasei[i][0] = Math.cos(ki[i]);
			phasei[i][1] = Math.sin(ki[i]);
		}
		
		logN = 0;
		int n = N;
		while(n > 1)
		{
			n/= 2; logN++;
		}
		storage = new double[logN][N][2];
	}
	
	//The algorithm is not ruined if f and fHat are the same array.
	void doFFT(double[][] f, double[][] fHat)
	{
		//populate the storage array.
		int i, j, k, m, p;
		int n = N/2, npieces = 1;
		for (i = 0; i < n; i++)
		{
			m = 2*i;
			storage[0][i][0] = f[m][0];
			storage[0][i][1] = f[m][1];
			storage[0][i + n][0] = f[m + 1][0];
			storage[0][i + n][1] = f[m + 1][1];
		}
		for(j = 1; j < logN; j++)
		{
			n /= 2; npieces *= 2;
			for (k = 0; k < npieces; k++)
				for (i = 0; i < n; i++)
				{
					m = 2*k*n;
					storage[j][i + m][0] = storage[j-1][2*i + m][0];
					storage[j][i + m][1] = storage[j-1][2*i + m][1];
					storage[j][i + n + m][0] = storage[j-1][2*i + 1 + m][0];
					storage[j][i + n + m][1] = storage[j-1][2*i + 1 + m][1];

				}
//			System.out.println(n + ", " + npieces + ", " + j + ", " + logN);
		}
		
		//Now we have in the last row of storage, N/2 even/odd pairs. We must first
		//re-arrange each pair, then go back up rearranging each row.
		//pair rearrangement:
		double[] temp1 = new double[2];
		double[] temp2 = new double[2];
		for (i = 0; i < npieces; i++)
		{
			m = 2*i;
			temp1[0] = storage[logN-1][m][0] + storage[logN-1][m + 1][0];
			temp1[1] = storage[logN-1][m][1] + storage[logN-1][m + 1][1];
			temp2[0] = storage[logN-1][m][0] - storage[logN-1][m + 1][0];
			temp2[1] = storage[logN-1][m][1] - storage[logN-1][m + 1][1];
			storage[logN-1][m][0] = temp1[0];
			storage[logN-1][m][1] = temp1[1];
			storage[logN-1][m + 1][0] = temp2[0];
			storage[logN-1][m + 1][1] = temp2[1];
		}
		//Now, recursively working back up the tree, we reach the top:
		//We are first combining N/2 pairs into N/4 quads: n = 2, npieces = 128
		for (j = logN - 1; j >= 1; j--)
		{
			n *= 2; npieces /= 2;
			for (k = 0; k < npieces; k++)
				for (i = 0; i < n; i++)
				{
					p = i*N/(2*n); m = 2*k*n;
					//phase is -2 pi i/n, the index of which is i*N/n
					temp1[0] = phasei[p][0]*storage[j][i + m + n][0] - phasei[p][1]*storage[j][i + m + n][1];
					temp1[1] = phasei[p][1]*storage[j][i + m + n][0] + phasei[p][0]*storage[j][i + m + n][1];
					storage[j-1][i + m][0] = temp1[0] + storage[j][i + m][0];
					storage[j-1][i + m][1] = temp1[1] + storage[j][i + m][1];
					storage[j-1][i + m + n][0] = storage[j][i + m][0] - temp1[0];
					storage[j-1][i + m + n][1] = storage[j][i + m][1] - temp1[1];
				}
//			System.out.println(n + ", " + npieces + ", " + j + ", " + logN);
		}
//		for (i = 0; i < N; i++){
//			for (j = 0; j < logN; j++)
//			{
//				System.out.print(Math.sqrt(storage[j][i][0]*storage[j][i][0] + storage[j][i][1]*storage[j][i][1]) + "\t");
//			}
//		System.out.println();
//		}
		
		//Lastly we do the final step transferring the stuff in storage[0] to fHat.
		for (i = 0; i < N; i++)
		{
			fHat[i][0] = storage[0][i][0];
			fHat[i][1] = storage[0][i][1];
//			System.out.println(fHat[i][0]*fHat[i][0] + fHat[i][1]*fHat[i][1]);
		}
	}
	void doFFTRe(double[] f, double[][] fHat)
	{
		//populate the storage array.
		int i, j, k, m, p;
		int n = N/2, npieces = 1;
		for (i = 0; i < n; i++)
		{
			m = 2*i;
			storage[0][i][0] = f[m];
			storage[0][i][1] = 0;
			storage[0][i + n][0] = f[m + 1];
			storage[0][i + n][1] = 0;
		}
		for(j = 1; j < logN; j++)
		{
			n /= 2; npieces *= 2;
			for (k = 0; k < npieces; k++)
				for (i = 0; i < n; i++)
				{
					m = 2*k*n;
					storage[j][i + m][0] = storage[j-1][2*i + m][0];
					storage[j][i + m][1] = storage[j-1][2*i + m][1];
					storage[j][i + n + m][0] = storage[j-1][2*i + 1 + m][0];
					storage[j][i + n + m][1] = storage[j-1][2*i + 1 + m][1];

				}
//			System.out.println(n + ", " + npieces + ", " + j + ", " + logN);
		}
		
		//Now we have in the last row of storage, N/2 even/odd pairs. We must first
		//re-arrange each pair, then go back up rearranging each row.
		//pair rearrangement:
		double[] temp1 = new double[2];
		double[] temp2 = new double[2];
		for (i = 0; i < npieces; i++)
		{
			m = 2*i;
			temp1[0] = storage[logN-1][m][0] + storage[logN-1][m + 1][0];
			temp1[1] = storage[logN-1][m][1] + storage[logN-1][m + 1][1];
			temp2[0] = storage[logN-1][m][0] - storage[logN-1][m + 1][0];
			temp2[1] = storage[logN-1][m][1] - storage[logN-1][m + 1][1];
			storage[logN-1][m][0] = temp1[0];
			storage[logN-1][m][1] = temp1[1];
			storage[logN-1][m + 1][0] = temp2[0];
			storage[logN-1][m + 1][1] = temp2[1];
		}
		//Now, recursively working back up the tree, we reach the top:
		//We are first combining N/2 pairs into N/4 quads: n = 2, npieces = 128
		for (j = logN - 1; j >= 1; j--)
		{
			n *= 2; npieces /= 2;
			for (k = 0; k < npieces; k++)
				for (i = 0; i < n; i++)
				{
					p = i*N/(2*n); m = 2*k*n;
					//phase is -2 pi i/n, the index of which is i*N/n
					temp1[0] = phasei[p][0]*storage[j][i + m + n][0] - phasei[p][1]*storage[j][i + m + n][1];
					temp1[1] = phasei[p][1]*storage[j][i + m + n][0] + phasei[p][0]*storage[j][i + m + n][1];
					storage[j-1][i + m][0] = temp1[0] + storage[j][i + m][0];
					storage[j-1][i + m][1] = temp1[1] + storage[j][i + m][1];
					storage[j-1][i + m + n][0] = storage[j][i + m][0] - temp1[0];
					storage[j-1][i + m + n][1] = storage[j][i + m][1] - temp1[1];
				}
//			System.out.println(n + ", " + npieces + ", " + j + ", " + logN);
		}
//		for (i = 0; i < N; i++){
//			for (j = 0; j < logN; j++)
//			{
//				System.out.print(Math.sqrt(storage[j][i][0]*storage[j][i][0] + storage[j][i][1]*storage[j][i][1]) + "\t");
//			}
//		System.out.println();
//		}
		
		//Lastly we do the final step transferring the stuff in storage[0] to fHat.
		for (i = 0; i < N; i++)
		{
			fHat[i][0] = storage[0][i][0];
			fHat[i][1] = storage[0][i][1];
//			System.out.println(fHat[i][0]*fHat[i][0] + fHat[i][1]*fHat[i][1]);
		}
	}
	
	//This is exactly the same as doFFT, except we do a column
	void doFFT(double[][][] f, double[][][] fHat, int a)
	{
		//populate the storage array.
		int i, j, k, m, p;
		int n = N/2, npieces = 1;
		for (i = 0; i < n; i++)
		{
			m = 2*i;
			storage[0][i][0] = f[m][a][0];
			storage[0][i][1] = f[m][a][1];
			storage[0][i + n][0] = f[m + 1][a][0];
			storage[0][i + n][1] = f[m + 1][a][1];
		}
		for(j = 1; j < logN; j++)
		{
			n /= 2; npieces *= 2;
			for (k = 0; k < npieces; k++)
				for (i = 0; i < n; i++)
				{
					m = 2*k*n;
					storage[j][i + m][0] = storage[j-1][2*i + m][0];
					storage[j][i + m][1] = storage[j-1][2*i + m][1];
					storage[j][i + n + m][0] = storage[j-1][2*i + 1 + m][0];
					storage[j][i + n + m][1] = storage[j-1][2*i + 1 + m][1];

				}
//			System.out.println(n + ", " + npieces + ", " + j + ", " + logN);
		}
		
		//Now we have in the last row of storage, N/2 even/odd pairs. We must first
		//re-arrange each pair, then go back up rearranging each row.
		//pair rearrangement:
		double[] temp1 = new double[2];
		double[] temp2 = new double[2];
		for (i = 0; i < npieces; i++)
		{
			m = 2*i;
			temp1[0] = storage[logN-1][m][0] + storage[logN-1][m + 1][0];
			temp1[1] = storage[logN-1][m][1] + storage[logN-1][m + 1][1];
			temp2[0] = storage[logN-1][m][0] - storage[logN-1][m + 1][0];
			temp2[1] = storage[logN-1][m][1] - storage[logN-1][m + 1][1];
			storage[logN-1][m][0] = temp1[0];
			storage[logN-1][m][1] = temp1[1];
			storage[logN-1][m + 1][0] = temp2[0];
			storage[logN-1][m + 1][1] = temp2[1];
		}
		//Now, recursively working back up the tree, we reach the top:
		//We are first combining N/2 pairs into N/4 quads: n = 2, npieces = 128
		for (j = logN - 1; j >= 1; j--)
		{
			n *= 2; npieces /= 2;
			for (k = 0; k < npieces; k++)
				for (i = 0; i < n; i++)
				{
					p = i*N/(2*n); m = 2*k*n;
					//phase is -2 pi i/n, the index of which is i*N/n
					temp1[0] = phasei[p][0]*storage[j][i + m + n][0] - phasei[p][1]*storage[j][i + m + n][1];
					temp1[1] = phasei[p][1]*storage[j][i + m + n][0] + phasei[p][0]*storage[j][i + m + n][1];
					storage[j-1][i + m][0] = temp1[0] + storage[j][i + m][0];
					storage[j-1][i + m][1] = temp1[1] + storage[j][i + m][1];
					storage[j-1][i + m + n][0] = storage[j][i + m][0] - temp1[0];
					storage[j-1][i + m + n][1] = storage[j][i + m][1] - temp1[1];
				}
//			System.out.println(n + ", " + npieces + ", " + j + ", " + logN);
		}
//		for (i = 0; i < N; i++){
//			for (j = 0; j < logN; j++)
//			{
//				System.out.print(Math.sqrt(storage[j][i][0]*storage[j][i][0] + storage[j][i][1]*storage[j][i][1]) + "\t");
//			}
//		System.out.println();
//		}
		
		//Lastly we do the final step transferring the stuff in storage[0] to fHat.
		for (i = 0; i < N; i++)
		{
			fHat[i][a][0] = storage[0][i][0];
			fHat[i][a][1] = storage[0][i][1];
//			System.out.println(fHat[i][0]*fHat[i][0] + fHat[i][1]*fHat[i][1]);
		}
	}
	void doFFTRe(double[][] f, double[][][] fHat, int a)
	{
		//populate the storage array.
		int i, j, k, m, p;
		int n = N/2, npieces = 1;
		for (i = 0; i < n; i++)
		{
			m = 2*i;
			storage[0][i][0] = f[m][a];
			storage[0][i][1] = 0;
			storage[0][i + n][0] = f[m + 1][a];
			storage[0][i + n][1] = 0;
		}
		for(j = 1; j < logN; j++)
		{
			n /= 2; npieces *= 2;
			for (k = 0; k < npieces; k++)
				for (i = 0; i < n; i++)
				{
					m = 2*k*n;
					storage[j][i + m][0] = storage[j-1][2*i + m][0];
					storage[j][i + m][1] = storage[j-1][2*i + m][1];
					storage[j][i + n + m][0] = storage[j-1][2*i + 1 + m][0];
					storage[j][i + n + m][1] = storage[j-1][2*i + 1 + m][1];

				}
//			System.out.println(n + ", " + npieces + ", " + j + ", " + logN);
		}
		
		//Now we have in the last row of storage, N/2 even/odd pairs. We must first
		//re-arrange each pair, then go back up rearranging each row.
		//pair rearrangement:
		double[] temp1 = new double[2];
		double[] temp2 = new double[2];
		for (i = 0; i < npieces; i++)
		{
			m = 2*i;
			temp1[0] = storage[logN-1][m][0] + storage[logN-1][m + 1][0];
			temp1[1] = storage[logN-1][m][1] + storage[logN-1][m + 1][1];
			temp2[0] = storage[logN-1][m][0] - storage[logN-1][m + 1][0];
			temp2[1] = storage[logN-1][m][1] - storage[logN-1][m + 1][1];
			storage[logN-1][m][0] = temp1[0];
			storage[logN-1][m][1] = temp1[1];
			storage[logN-1][m + 1][0] = temp2[0];
			storage[logN-1][m + 1][1] = temp2[1];
		}
		//Now, recursively working back up the tree, we reach the top:
		//We are first combining N/2 pairs into N/4 quads: n = 2, npieces = 128
		for (j = logN - 1; j >= 1; j--)
		{
			n *= 2; npieces /= 2;
			for (k = 0; k < npieces; k++)
				for (i = 0; i < n; i++)
				{
					p = i*N/(2*n); m = 2*k*n;
					//phase is -2 pi i/n, the index of which is i*N/n
					temp1[0] = phasei[p][0]*storage[j][i + m + n][0] - phasei[p][1]*storage[j][i + m + n][1];
					temp1[1] = phasei[p][1]*storage[j][i + m + n][0] + phasei[p][0]*storage[j][i + m + n][1];
					storage[j-1][i + m][0] = temp1[0] + storage[j][i + m][0];
					storage[j-1][i + m][1] = temp1[1] + storage[j][i + m][1];
					storage[j-1][i + m + n][0] = storage[j][i + m][0] - temp1[0];
					storage[j-1][i + m + n][1] = storage[j][i + m][1] - temp1[1];
				}
//			System.out.println(n + ", " + npieces + ", " + j + ", " + logN);
		}
//		for (i = 0; i < N; i++){
//			for (j = 0; j < logN; j++)
//			{
//				System.out.print(Math.sqrt(storage[j][i][0]*storage[j][i][0] + storage[j][i][1]*storage[j][i][1]) + "\t");
//			}
//		System.out.println();
//		}
		
		//Lastly we do the final step transferring the stuff in storage[0] to fHat.
		for (i = 0; i < N; i++)
		{
			fHat[i][a][0] = storage[0][i][0];
			fHat[i][a][1] = storage[0][i][1];
//			System.out.println(fHat[i][0]*fHat[i][0] + fHat[i][1]*fHat[i][1]);
		}
	}
	
	//In this method we fft f starting from a given index x.
	//Other than that, no difference
	void doFFT(double[][] f, double[][] fHat, int x)
	{
		//populate the storage array.
		int i, j, k, m, p;
		int n = N/2, npieces = 1;
		for (i = 0; i < n; i++)
		{
			m = 2*i;
			storage[0][i][0] = f[x + m][0];
			storage[0][i][1] = f[x + m][1];
			storage[0][i + n][0] = f[x + m + 1][0];
			storage[0][i + n][1] = f[x + m + 1][1];
		}
		for(j = 1; j < logN; j++)
		{
			n /= 2; npieces *= 2;
			for (k = 0; k < npieces; k++)
				for (i = 0; i < n; i++)
				{
					m = 2*k*n;
					storage[j][i + m][0] = storage[j-1][2*i + m][0];
					storage[j][i + m][1] = storage[j-1][2*i + m][1];
					storage[j][i + n + m][0] = storage[j-1][2*i + 1 + m][0];
					storage[j][i + n + m][1] = storage[j-1][2*i + 1 + m][1];

				}
//			System.out.println(n + ", " + npieces + ", " + j + ", " + logN);
		}
		
		//Now we have in the last row of storage, N/2 even/odd pairs. We must first
		//re-arrange each pair, then go back up rearranging each row.
		//pair rearrangement:
		double[] temp1 = new double[2];
		double[] temp2 = new double[2];
		for (i = 0; i < npieces; i++)
		{
			m = 2*i;
			temp1[0] = storage[logN-1][m][0] + storage[logN-1][m + 1][0];
			temp1[1] = storage[logN-1][m][1] + storage[logN-1][m + 1][1];
			temp2[0] = storage[logN-1][m][0] - storage[logN-1][m + 1][0];
			temp2[1] = storage[logN-1][m][1] - storage[logN-1][m + 1][1];
			storage[logN-1][m][0] = temp1[0];
			storage[logN-1][m][1] = temp1[1];
			storage[logN-1][m + 1][0] = temp2[0];
			storage[logN-1][m + 1][1] = temp2[1];
		}
		//Now, recursively working back up the tree, we reach the top:
		//We are first combining N/2 pairs into N/4 quads: n = 2, npieces = 128
		for (j = logN - 1; j >= 1; j--)
		{
			n *= 2; npieces /= 2;
			for (k = 0; k < npieces; k++)
				for (i = 0; i < n; i++)
				{
					p = i*N/(2*n); m = 2*k*n;
					//phase is -2 pi i/n, the index of which is i*N/n
					temp1[0] = phasei[p][0]*storage[j][i + m + n][0] - phasei[p][1]*storage[j][i + m + n][1];
					temp1[1] = phasei[p][1]*storage[j][i + m + n][0] + phasei[p][0]*storage[j][i + m + n][1];
					storage[j-1][i + m][0] = temp1[0] + storage[j][i + m][0];
					storage[j-1][i + m][1] = temp1[1] + storage[j][i + m][1];
					storage[j-1][i + m + n][0] = storage[j][i + m][0] - temp1[0];
					storage[j-1][i + m + n][1] = storage[j][i + m][1] - temp1[1];
				}
//			System.out.println(n + ", " + npieces + ", " + j + ", " + logN);
		}
//		for (i = 0; i < N; i++){
//			for (j = 0; j < logN; j++)
//			{
//				System.out.print(Math.sqrt(storage[j][i][0]*storage[j][i][0] + storage[j][i][1]*storage[j][i][1]) + "\t");
//			}
//		System.out.println();
//		}
		
		//Lastly we do the final step transferring the stuff in storage[0] to fHat.
		for (i = 0; i < N; i++)
		{
			fHat[i][0] = storage[0][i][0];
			fHat[i][1] = storage[0][i][1];
//			System.out.println(fHat[i][0]*fHat[i][0] + fHat[i][1]*fHat[i][1]);
		}
	}

	public void doFFT()
	{
		doFFT(f, fHat);
	}
	
	public void doIFFT()
	{
		//Take the star:
		ArrayOps.star(fHat, f);
		doFFT(f, f);
		ArrayOps.star(f, f);
		
		//Then the columns
		for (int i = 0; i < N; i++)
		{
			f[i][0] /= N;
			f[i][1] /= N;
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
	
	
	
	public static void main(String[] args)
	{
		double[][] x = new double[512][2];
		for (int i = 0; i < 512; i++){
			x[i][0] = i == 201 ? 1 : 0; x[i][1] = 0;
		}
		x[52][0] = 1;
		FFT1D fft = new FFT1D(x);
		fft.doFFT();
	}
	
}
