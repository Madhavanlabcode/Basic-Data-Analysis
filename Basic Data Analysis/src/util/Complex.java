package util;

//This performs complex arithmetic with double[] arrays.
public class Complex {

	public static double[] product(double[] a, double[] b)
	{
		return new double[] {a[0]*b[0] - a[1]*b[1], a[0]*b[1] + a[1]*b[0]};
	}
	public static double[] getQuotient(double[] a, double[] b)
	{
		double denom = b[0]*b[0] + b[1]*b[1];
		return new double[] {(a[0]*b[0] + a[1]*b[1])/denom, (-a[0]*b[1]+a[1]*b[0])/denom};
	}
	public static double[] getQuotient(double[] a, double[] b, double lambda)
	{
		double denom = b[0]*b[0] + b[1]*b[1] + lambda;
		return new double[] {(a[0]*b[0] + a[1]*b[1])/denom, (-a[0]*b[1]+a[1]*b[0])/denom};
	}
	public static double mag(double[] a)
	{
		return Math.sqrt(a[0]*a[0] + a[1]*a[1]);
	}
	public static double phase(double[] a)
	{
		return FieldOps.atan(a[0], a[1]);
	}
	public static double phaseCenteredZero(double[] a)
	{
		return FieldOps.atanCenteredZero(a[0], a[1]);
	}
	public static double mag(double x, double y)
	{
		return Math.sqrt(x*x + y*y);
	}
	public static void product(double[] a, double[] b, double[] target)
	{
		target[0] = a[0]*b[0] - a[1]*b[1];
		target[1] = a[0]*b[1] + a[1]*b[0];
	}
	public static void product(double are, double aim, double bre, double bim, double[] target)
	{
		target[0] = are*bre - aim*bim;
		target[1] = are*bim + aim*bre;
	}
	public static double[] expi(double x)
	{
		return new double[] {Math.cos(x), Math.sin(x)};
	}
	public static double[] expmi(double x)
	{
		return new double[] {Math.cos(x), -Math.sin(x)};
	}
	public static void expi(double x, double[] target)
	{
		target[0] = Math.cos(x);
		target[1] = Math.sin(x);
	}
	public static void expmi(double x, double[] target)
	{
		target[0] = Math.cos(x);
		target[1] = -Math.sin(x);
	}
	public static void sum(double[] a, double[] b, double[] target)
	{
		target[0] = a[0]+b[0];
		target[1] = a[1]+b[1];
	}
	public static void times(double[] a, double real, double[] target)
	{
		target[0] = a[0]*real;
		target[1] = a[1]*real;
	}
}
