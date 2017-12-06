package util.fourier;

//This class is a "slow fourier transform."
public class SFT {

	int N;
	double[][] sin;
	double[][] cos;
	double[] kn;
	
	public SFT(int N)
	{
		this.N = N;
		sin = new double[N][N];
		cos = new double[N][N];
		kn  = new double[N];
		
		for (int i = 0; i < N; i++)
		{
			kn[i] = (2*Math.PI/N)*(i - N/2); 
			for (int j = 0; j < N; j++)
			{
				sin[i][j] = Math.sin(kn[i]*j);
				cos[i][j] = Math.cos(kn[i]*j);
			}
		}
		
	}
	
	//we assume everything is size N*N;
	public void sft(double[][] data, double[][][] sft)
	{
//		sft[0] is real. sft[1] is imaginary.
		//A_k = sum(x) sum(y) e^i*k*y data(x, y) e^i*k*x
		int i, j, k, m;
		double[][] temp = new double[2][N];
		for (k = 0; k < N; k++){
		System.out.println("" + k + " ");
			for (m = 0; m < N; m++)
			{
				System.out.print("" + m + " ");
				sft[0][k][m] = 0; sft[1][k][m] = 0;
				for (i = 0; i < N; i++){
					temp[0][i] = 0; temp[1][i] = 0;
					for (j = 0; j < N; j++)
					{
						temp[0][i] += cos[m][j]*data[i][j];
						temp[1][i] += sin[m][j]*data[i][j];
					}
					sft[0][k][m] += temp[0][i]*cos[k][i] - temp[1][i]*sin[k][i];
					sft[1][k][m] += temp[0][i]*sin[k][i] + temp[1][i]*cos[k][i];
				}
			}
		}
	}
}
