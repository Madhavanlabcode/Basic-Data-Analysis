package util.scalar.functions;

public class TwoDArray {

	public static double[][] getGaussMask(double L, int N)
	{
		double[][] gauss = new double [N][N];
		for (int k = 0; k < N; k++)
			for (int l = 0; l < N; l++)
				gauss[k][l] = Math.exp((-((double)(k*k) + (l*l))/(2*L*L)));
		return gauss;
	}
	public static void getGaussMask(double L, int N, double[][] gauss)
	{
		for (int k = 0; k < N; k++)
			for (int l = 0; l < N; l++)
				gauss[k][l] = Math.exp((-((double)(k*k) + (l*l))/(2*L*L)));
	}

	public static void getCrossMask(double L, double xi, int N, double[][] gauss)
	{
		for (int k = 0; k < N; k++)
			for (int l = 0; l < N; l++)
				gauss[k][l] = Math.exp((-((double)(k*k) + (l*l) - xi*(k-l)*(k-l))/(2*L*L)));
	}
}
