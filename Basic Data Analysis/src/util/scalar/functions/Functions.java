package util.scalar.functions;

import flanagan.interpolation.CubicSpline;
import util.scalar.Function;

public class Functions {

	public static class Gaussian implements Function
	{
		private double x0, dx;
		private double n, twoDx2;
		
		public Gaussian(double X0, double Dx)
		{
			x0 = X0;
			dx = Dx;
			
			n = 1/Math.sqrt(2*Math.PI*dx*dx);
			twoDx2 = 2*dx*dx;
		}
		
		public double of(double x)
		{
			double x1 = x - x0;
			return n*Math.exp(-(x1*x1)/twoDx2);
		}
	}
	/**
	 * Symmetric in the dimensions: x -> r
	 * @author madhavanlab2011
	 *
	 */
	public static class Gaussian2D implements Function
	{
		private double x0, dx;
		private double n, twoDx2;
		
		public Gaussian2D(double X0, double Dx)
		{
			x0 = X0;
			dx = Dx;
			
			n = 1/(2*Math.PI*dx*dx);
			twoDx2 = 2*dx*dx;
		}
		
		public double of(double x)
		{
			double x1 = x - x0;
			return n*Math.exp(-(x1*x1)/twoDx2);
		}
	}
	public static class Line implements Function
	{
		private double m, b;
		
		public Line(double m, double b)
		{
			this.m = m;
			this.b = b;
		}
		
		public double of(double x)
		{
			return m*x + b;
		}
	}
	
	/**
	 * Wraps the CubicSpline of Flanagan providing default constant values on either side of the interpolation range.
	 * @author madhavanlab2011
	 *
	 */
	public static class CubicSplineWrapper implements Function
	{
		double xmin, xmax;
		double yBelow, yAbove;
		CubicSpline f;
		public CubicSplineWrapper(CubicSpline f, double yBelow, double yAbove)
		{
			this.f = f;
			this.yBelow = yBelow;
			this.yAbove = yAbove;
			xmin = f.getXmin();
			xmax = f.getXmax();
		}
		public CubicSplineWrapper(double[] x, double[] y)
		{
			f = new CubicSpline(x, y);
			yBelow = y[0];
			yAbove = y[y.length-1];
			xmin = f.getXmin();
			xmax = f.getXmax();
		}
		public CubicSplineWrapper(double[] x, double[] y, double yBelow, double yAbove)
		{
			f = new CubicSpline(x, y);
			this.yBelow = yBelow;
			this.yAbove = yAbove;
			xmin = f.getXmin();
			xmax = f.getXmax();
		}
		@Override
		public double of(double x) {
			// TODO Auto-generated method stub
			if (x < xmin) return yBelow;
			if (x > xmax) return yAbove;
			return f.interpolate(x);
		}
	}
	
	public static class LagrangeInteroplationWrapper implements Function
	{
		double[] x, y;
		public LagrangeInteroplationWrapper(double[] x, double[] y)
		{
			this.x = x;
			this.y = y;
		}
		@Override
		public double of(double x) {
			return aitken(x, this.x, y);
		}
	}
	
	/**
	 * This method performs Lagrange interpolation to get the value of the function defined by the two arrays xi and fi, at the point x.
	 * @author Tao Peng (2006)
	 * @param x
	 * @param xi
	 * @param fi
	 * @return
	 */
	  public static double aitken(double x, double xi[],
			    double fi[]) {
			    int n = xi.length-1;
			    double ft[] = (double[]) fi.clone();
			    for (int i=0; i<n; ++i) {
			      for (int j=0; j<n-i; ++j) {
			        ft[j] = (x-xi[j])/(xi[i+j+1]-xi[j])*ft[j+1]
			               +(x-xi[i+j+1])/(xi[j]-xi[i+j+1])*ft[j];
			      }
			    }
			    return ft[0];
			  }

}
