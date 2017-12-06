package main;

public class CentroidField {

	//assume everything is square
	static double[][] centroidField(double[][] data, int border, int pieceSize, int csize)
	{
		int N = data.length, b = border, s = pieceSize;
		int i, j, k, m;
		double[][] field = new double [N-2*b][N-2*b];
		boolean checkedge = b <= s;

		double[] centroid;
		int[] round = new int[2];
		boolean edge;
		for (i = b; i < N-b; i++)
			for (j = b; j < N-b; j++)
			{
				centroid = centroid(data, i-s/2, s, j-s/2, s, checkedge);
				round[0] = (int)centroid[0]; round[1] = (int)centroid[1];
				edge = centroid[0] < csize/2 || centroid[0] >= N-1-csize/2 || centroid[1] < csize/2 || centroid[1] >= N-1-csize/2;
				if (!edge){
					for (k = round[0] - csize/2; k < round[0] + csize/2 - 1; k++)
						for (m = round[1] - csize/2; m < round[1] + csize/2 - 1; m++)
							field[k][m]++;
//					field[round[0]][round[1]] += 9;
//					field[round[0]+1][round[1]] += 4;
//					field[round[0]-1][round[1]] += 4;
//					field[round[0]][round[1]+1] += 4;
//					field[round[0]][round[1]-1] += 4;
//					field[round[0]+1][round[1]+1] += 1;
//					field[round[0]-1][round[1]+1] += 1;
//					field[round[0]+1][round[1]-1] += 1;
//					field[round[0]-1][round[1]-1] += 1;
				}
			}
		return field;
	}

	static double[][] centroidField2(double[][] data, int border, int pieceSize, int csize)
	{
		int N = data.length, b = border, s = pieceSize;
		int i, j, k, m;
		double[][] field = new double [N-2*b][N-2*b];
		boolean checkedge = b <= s;

		double[] centroid;
		double[] sigma = new double[2], censq;
		double weight;
		int[] round = new int[2];
		boolean edge;
		for (i = b; i < N-b; i++)
			for (j = b; j < N-b; j++)
			{
				centroid = centroid(data, i-s/2, s, j-s/2, s, checkedge);
				censq = cenSq(data, i-s/2, s, j-s/2, s, checkedge);
				sigma[0] = censq[0] - centroid[0]*centroid[0];
				sigma[1] = censq[1] - centroid[1]*centroid[1];
				weight = Math.min(sigma[0]/sigma[1], sigma[1]/sigma[0]);
				round[0] = (int)centroid[0]; round[1] = (int)centroid[1];
				edge = centroid[0] < csize/2 || centroid[0] >= N-1-csize/2 || centroid[1] < csize/2 || centroid[1] >= N-1-csize/2;
				if (!edge){
					for (k = round[0] - csize/2; k < round[0] + csize/2 - 1; k++)
						for (m = round[1] - csize/2; m < round[1] + csize/2 - 1; m++)
							field[k][m] += weight/Math.max(sigma[0], sigma[1]);
//					field[round[0]][round[1]] += 9;
//					field[round[0]+1][round[1]] += 4;
//					field[round[0]-1][round[1]] += 4;
//					field[round[0]][round[1]+1] += 4;
//					field[round[0]][round[1]-1] += 4;
//					field[round[0]+1][round[1]+1] += 1;
//					field[round[0]-1][round[1]+1] += 1;
//					field[round[0]+1][round[1]-1] += 1;
//					field[round[0]-1][round[1]-1] += 1;
				}
			}
		return field;
	}
	
	static void smooth(double[][] data, int smoothness)
	{
		//we assume smoothness is 2*n + 1, so that (smoothness - 1)/2 = n;
		int n = (smoothness - 1)/2;
		int N = data.length, M = data[0].length;
		double temp;
		int area = smoothness*smoothness;
		int i, j, k, p;
		//first the easy part
		for (i = n; i < N - n; i++)
			for (j = n; j < M - n; j++){
				temp = 0;
				for (k = -n; k < n+1; k++)
					for (p = -n; p < n+1; p++)
					{
						temp += data[i+k][j+p];
					}
				data[i][j] = temp/area;
			}
		//ignore the hard part.
	}
	
	public static double[] centroid(double[][] field, int x, int dx, int y, int dy, boolean checkedge)
	{
		double sum = 0;
		
		double[] centroid = new double[2];
		if (checkedge)
			for (int i = x; i < x + dx; i++)
				for (int j = y; j < y + dy; j++)
				{
					if (i >= 0 && i < field.length && j >= 0 && j < field.length){
						centroid[0] += i*field[i][j];
						centroid[1] += j*field[i][j];
						sum += field[i][j];
					}
				}
		else
			for (int i = x; i < x + dx; i++){
				for (int j = y; j < y + dy; j++)
				{
					centroid[0] += i*field[i][j];
					centroid[1] += j*field[i][j];
					sum += field[i][j];
//					System.out.print(Math.log(field[i][j]) + "\t");
				}
//				System.out.println();
			}
		
		centroid[0] /= sum;
		centroid[1] /= sum;
		
		return centroid;
	}
	static double[] cenSq(double[][] field, int x, int dx, int y, int dy, boolean checkedge)
	{
		double sum = 0;
		
		double[] centroid = new double[2];
		if (checkedge)
			for (int i = x; i < x + dx; i++)
				for (int j = y; j < y + dy; j++)
				{
					if (i >= 0 && i < field.length && j >= 0 && j < field.length){
						centroid[0] += i*i*field[i][j];
						centroid[1] += j*j*field[i][j];
						sum += field[i][j];
					}
				}
		else
			for (int i = x; i < x + dx; i++)
				for (int j = y; j < y + dy; j++)
				{
					centroid[0] += i*i*field[i][j];
					centroid[1] += j*j*field[i][j];
					sum += field[i][j];
			}
		
		centroid[0] /= sum;
		centroid[1] /= sum;
		
		return centroid;
	}
}
