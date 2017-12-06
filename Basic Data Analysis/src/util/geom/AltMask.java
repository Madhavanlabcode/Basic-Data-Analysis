package util.geom;

import util.FieldOps;
import util.matrix.Matrix;

public class AltMask {
	
	
	public static void putGaussSidedSquareMaskRotated(int nx, int ny, int s, double[] c, double theta, double thickness, double[][] mask, double[][] temp, double[][] rotMatrix)
	{
		Matrix.putRotationMatrix(theta, rotMatrix);
//		putGaussSidedSquareMask(nx, nx, s, c, thickness, temp);
		FieldOps.applyLinearTransformation(temp, rotMatrix, c, mask);
	}
	public static void putSolidGaussSidedSquareMaskRotated(int nx, int ny, int s, double[] c, double theta, double thickness, double[][] mask, double[][] temp, double[][] rotMatrix)
	{
		Matrix.putRotationMatrix(theta, rotMatrix);
//		putSolidGaussSidedSquareMask(nx, nx, s, c, thickness, temp);
		FieldOps.applyLinearTransformation(temp, rotMatrix, c, mask);
	}
	
	
	
	
	public static void putRectMaskRotated(int nx, int ny, int side, double[] c, double theta, double[][] mask, double[][] rotMatrix)
	{
		Matrix.putRotationMatrix(theta, rotMatrix);
		double[] r = new double [2];
		double[] rPrime = new double [2];
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				r[0] = i; r[1] = j;
				TwoD.vectorChangeOrigin(r, rotMatrix, rPrime, c);
				if (Distance.distanceToSquareBody(side, c[0], c[1], rPrime[0], rPrime[1]) == 0)
					mask[i][j] = 1;
				else
					mask[i][j] = 0;
			}
	}
	
}
