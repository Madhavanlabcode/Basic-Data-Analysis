package util.geom;

import util.FieldOps;
import util.matrix.Matrix;

public class Mask {

	public static double[][] rectMaskCent(int nx, int ny, int lx, int ly)
	{
		double[][] mask = new double [nx][ny];
		int marginx = (nx-lx)/2;
		int marginy = (ny-ly)/2;
		for (int i = 0; i < lx; i++)
			for (int j = 0; j < ly; j++)
				mask[i+marginx][j+marginy] = 1;
		return mask;
	}
	
	public static double[][] rectMask(int nx, int ny, int lx, int ly, int[] c)
	{
		double[][] mask = new double [nx][ny];
		for (int i = 0; i < lx; i++)
			for (int j = 0; j < ly; j++)
				if (i + c[0] >= 0 && i + c[0] < nx && j + c[1] >= 0 && j + c[1] < ny)
					mask[i+c[0]][j+c[1]] = 1;
		return mask;
	}
	
	//center c, side length s, side thickness t;
	public static double[][] gaussSidedSquareMask(int nx, int ny, int s, int[] c, double thickness)
	{
		double[][] mask = new double[nx][ny];
		
		double d = 0;
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				d = Distance.distanceToSquareEdge(s, c[0], c[1], i, j);
				mask[i][j] = Math.exp(-(d*d/(thickness*thickness)));
			}
		mask[0][0] = 1;
		return mask;
	}
	public static void putGaussSidedSquareMask(int nx, int ny, int s, double[] c, double thickness, double[][] mask)
	{
		double d = 0;
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				d = Distance.distanceToSquareEdge(s, c[0], c[1], i, j);
				mask[i][j] = Math.exp(-(d*d/(thickness*thickness)));
			}
		mask[0][0] = 1;
	}
	public static double[][] solidGaussSidedSquareMask(int nx, int ny, int s, int[] c, double thickness)
	{
		double[][] mask = new double[nx][ny];
		
		double d = 0;
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				d = Distance.distanceToSquareBody(s, c[0], c[1], i, j);
				mask[i][j] = Math.exp(-(d*d/(thickness*thickness)));
			}
		return mask;
	}
	public static void putSolidGaussSidedSquareMask(int nx, int ny, int s, double[] c, double thickness, double[][] mask)
	{
		double d = 0;
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				d = Distance.distanceToSquareBody(s, c[0], c[1], i, j);
				mask[i][j] = Math.exp(-(d*d/(thickness*thickness)));
			}
	}
	
	//returns simple gaussian e^-r^2/x^2 without other factors
	public static double[][] getGaussianMask(double x)
	{
		int L = (int)Math.round((4*x+2));
		double d;
		double[][] g = new double [L][L];
		for (int i = 0; i < L; i++)
			for (int j = 0; j < L; j++)
			{
				d = Distance.distance(i, j);
				g[i][j] = Math.exp(-(d*d)/(x*x));
			}
		return g;
	}
	
	
	public static void putRectMask(int nx, int ny, int lx, int ly, double[] c, double[][] mask)
	{
		FieldOps.zero(mask);
		int x, y;
		for (int i = 0; i < lx; i++)
			for (int j = 0; j < ly; j++)
			{
				x = i + (int)Math.round(c[0]) - lx/2; y = j + (int)Math.round(c[1]) - ly/2;
				if (x >= 0 && x < nx && y >= 0 && y < ny)
					mask[x][y] = 1;
			}
	}
	public static void putGaussSidedSquareMaskRotated(int nx, int ny, int s, double[] c, double theta, double thickness, double[][] mask, double[][] temp, double[][] rotMatrix)
	{
		Matrix.putRotationMatrix(theta, rotMatrix);
		putGaussSidedSquareMask(nx, nx, s, c, thickness, temp);
		FieldOps.applyLinearTransformation(temp, rotMatrix, c, mask, 0);
	}
	public static void putSolidGaussSidedSquareMaskRotated(int nx, int ny, int s, double[] c, double theta, double thickness, double[][] mask, double[][] temp, double[][] rotMatrix)
	{
		Matrix.putRotationMatrix(theta, rotMatrix);
		putSolidGaussSidedSquareMask(nx, nx, s, c, thickness, temp);
		FieldOps.applyLinearTransformation(temp, rotMatrix, c, mask, 0);
	}
	public static void putRectMaskRotated(int nx, int ny, int lx, int ly, double[] c, double theta, double[][] temp, double[][] mask, double[][] rotMatrix)
	{
		Matrix.putRotationMatrix(theta, rotMatrix);
		putRectMask(nx, ny, lx, ly, c, temp);
		FieldOps.applyLinearTransformation(temp, rotMatrix, c, mask, 0);
	}
	
}
