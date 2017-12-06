package util.matrix;

import util.FieldOps;
import util.geom.AtomicCoordinatesSet;
import util.geom.Distance;

public class Matrix {
	
	public static double[] temp = new double [2];
	private static double[][] refX = new double[][] {{-1, 0}, {0, 1}};
	
	
	public static void putRotationMatrix(double theta, double[][] matrix)
	{
		matrix[0][0] = Math.cos(theta);
		matrix[0][1] = Math.sin(theta);
		matrix[1][0] = -Math.sin(theta);
		matrix[1][1] = Math.cos(theta);
	}
	public static double[][] getRotationMatrix(double theta)
	{
		double[][] matrix = new double [2][2];
		putRotationMatrix(theta, matrix);
		return matrix;
	}
	public static double[] getProductWith(double[][] mat, double[] vec)
	{
		double[] ans = new double [2];
		putProductWith(mat, vec, ans);
		return ans;
	}
	public static double[] getProductWithOrigin(double[][] mat, double[] vec, double[] origin)
	{
		double[] ans = new double [2];
		putProductWithOrigin(mat, vec, ans, origin);
		return ans;
	}
	public static void putProductWith(double[][] mat, double[] vec, double[] ans) {
		ans[0] = mat[0][0]*vec[0] + mat[1][0]*vec[1];
		ans[1] = mat[0][1]*vec[0] + mat[1][1]*vec[1];
	}
	public static void putProductWith(double[][] mat1, double[][] mat2, double[][] ans) {
		ans[0][0] = mat1[0][0]*mat2[0][0] + mat1[1][0]*mat2[0][1];
		ans[0][1] = mat1[0][1]*mat2[0][0] + mat1[1][1]*mat2[0][1];
		ans[1][0] = mat1[0][0]*mat2[1][0] + mat1[1][0]*mat2[1][1];
		ans[1][1] = mat1[0][1]*mat2[1][0] + mat1[1][1]*mat2[1][1];
	}
	public static double[][] getProductWith(double[][] mat1, double[][] mat2) {
		double[][] ans = new double [2][2];
		putProductWith(mat1, mat2, ans);
		return ans;
	}
	public static void putProductWithOrigin(double[][] mat, double[] vec, double[] ans, double[] origin) {
		temp[0] = vec[0] - origin[0];
		temp[1] = vec[1] - origin[1];
		ans[0] = (mat[0][0]*temp[0] + mat[1][0]*temp[1]) + origin[0];
		ans[1] = (mat[0][1]*temp[0] + mat[1][1]*temp[1]) + origin[1];
	}

	//gets the matrix for reflection about a vector
	public static double[][] getReflectionMatrix(double[] vector)
	{
		double[] unit = Distance.unitVector(vector);
		double theta = FieldOps.atan(unit[0], unit[1]);
		double[][] rot = new double [2][2], rotInv = new double [2][2];
		putRotationMatrix(-theta, rotInv);
		putRotationMatrix(theta, rot);
		double[][] ans = getProductWith(getProductWith(rot, refX), rotInv);	
		return ans;
	}
	
	public static double[] rotateVector(double[] vector, double angle)
	{
		return getProductWith(getRotationMatrix(angle), vector);
	}
	/**
	 * This method rotates a vector so that the unit vector in the direction of latt.getA() will point in the plus-X direction.
	 * 
	 */
	public static double[][] getRotationMatrix(AtomicCoordinatesSet latt)
	{
		double[] unitA = Distance.unitVector(latt.getA());
		double[] unitB = rotateVector(unitA, Math.PI/2);
		unitA[0] = -unitA[0];
		unitA[1] = -unitA[1];
		unitB[0] = -unitB[0];
		unitB[1] = -unitB[1];
		return new double[][] {unitA, unitB};
	}
}
