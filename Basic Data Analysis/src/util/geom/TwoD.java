package util.geom;

public class TwoD {

	//the rotation matrix IS the new basis.
	private static double[][] temp = new double [2][2];
	
	public static void vectorChange(double[] vector, double[][] matrix, double[] result)
	{
		result[0] = vector[0]*matrix[0][0] + vector[1]*matrix[1][0];
		result[1] = vector[0]*matrix[0][1] + vector[1]*matrix[1][1];
	}
	public static void vectorChangeOrigin(double[] vector, double[][] matrix, double[] result, double[] origin)
	{
		temp[0][0] = vector[0]-origin[0]; temp[0][1] = vector[1]-origin[1];
		vectorChange(temp[0], matrix, temp[1]);
		result[0] = temp[1][0]+origin[0];
		result[1] = temp[1][1]+origin[1];
	}
}
