package util;

public class MatrixOps {

	public static double[][] inv2x2(double[][] mat)
	{
		double det = mat[0][0]*mat[1][1] - mat[0][1]*mat[1][0];
		return new double[][] {{mat[1][1]/det, -mat[0][1]/det}, {-mat[1][0]/det, mat[0][0]/det}}; 
	}
}
