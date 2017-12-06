package util;

import util.regression.ACM_CustomFunctions;
import flanagan.analysis.Regression;

public class ArrayUtil {

	/**
	 * This returns the array of best-estimate fitting parameters for the function type named. 
	 * @param x	The x-axis data
	 * @param y	The y-axis data
	 * @param function	A string naming the type of function to be used (see CustomFunctions for a list of acceptable strings).
	 * @return
	 */
	public static double[] fitToFunction(double[] x, double[] y, String function)
	{
		ACM_CustomFunctions r; double[] start, step; Regression reg;
		r = ACM_CustomFunctions.getNew(function);
		double cent = (ArrayOps.max(x)+ArrayOps.min(x))/2;
		double xp = ArrayOps.max(x), xm = ArrayOps.min(x), ym = ArrayOps.min(y), yp = ArrayOps.max(y);
		double dx = xp-xm;
		double dy = yp-ym;

		double centroid = x[(int)Math.round(ArrayOps.centIndex(y))];
		double spread = ArrayOps.sigmaIndex(y)*(dx/x.length);
		double area = ArrayOps.sum(y)*(dx/x.length);

		if (function.equalsIgnoreCase("Lorentzian"))
		{
			//we have to cleverly select the initial values:
			//the total peak height is 2A0/(width*pi)
			//therefore A0 = height*width*pi/2
			start = new double [4];
			start[1] = ArrayOps.min(y);
			start[0] = cent;
			start[2] = dx/2;
			start[3] = dy*start[2]*Math.PI/2;
			
			step = new double [4];
			step[1] = dy/10;
			step[0] = dx/10;
			step[2] = dx/10;
			step[3] = start[3]/5;
			
//			double[] step = new double [start.length];
//			for (int i = 0; i < start.length; i++)
//				step[i] = start[i]/20;
			reg = new Regression(x, y);
			reg.simplex(r, start, step);
			return reg.getBestEstimates();
		}
		else if (function.equalsIgnoreCase("TwoLorentzian"))
		{
			start = new double [7];
			start[0] = cent - dx*0.1;
			start[1] = dx/3;
			start[2] = dy*start[1]*Math.PI/4;
			start[3] = cent+dx*0.1;
			start[4] = dx/3;
			start[5] = start[2];
			start[6] = ArrayOps.min(y);
			step = new double [7];
			step[0] = dx/10;
			step[1] = dx/15;
			step[2] = start[2]/10;
			step[3] = dx/10;
			step[4] = dx/15;
			step[5] = step[2];
			step[6] = dy/10;
			reg = new Regression(x, y);
			reg.simplex(r, start, step);
			return reg.getBestEstimates();
		}
		else if (function.equalsIgnoreCase("Fermi"))
		{
			start = new double [4];
			start[0] = cent;
			start[1] = dx/10;
			start[2] = dy;
			if (y[0] < y[x.length-1]) start[2] *= -1;
			start[3] = ArrayOps.min(y);
			
			step = new double [4];
			step[0] = dx/10;
			step[1] = dx/50;
			step[2] = dy/10;
			step[3] = dy/20;
			//will attempt step procedure now;
			reg = new Regression(x, y);
			reg.simplex(r, start, step);
			return reg.getBestEstimates();
		}
		else if (function.equalsIgnoreCase("FermiLine"))
		{
			start = new double [5];
			start[0] = cent;
			start[1] = dx/10;
			start[2] = dy;
			if (y[0] < y[x.length-1]) start[2] *= -1;
			start[3] = ArrayOps.min(y);
			start[4] = 0;
			
			step = new double [5];
			step[0] = dx/10;
			step[1] = dx/50;
			step[2] = dy/10;
			step[3] = dy/20;
			step[4] = (dy/dx)/10;
			//will attempt step procedure now;
			reg = new Regression(x, y);
			reg.simplex(r, start, step);
			return reg.getBestEstimates();
		}
		else if (function.equalsIgnoreCase("1PeakLorentzian_Line"))
		{
			start = new double [5];
			start[0] = cent;
			start[1] = dx/4;
			start[2] = dy*start[1]*Math.PI/4;
			start[3] = ym;
			start[4] = 0;
			step = new double [5];
			step[0] = dx/10;
			step[1] = dx/20;
			step[2] = start[2]/10;
			step[3] = dy/20;
			step[4] = (dy/dx)/10;
			reg = new Regression(x, y);
			reg.simplex(r, start, step);
			return reg.getBestEstimates();
		}
		else if (function.equalsIgnoreCase("TwoGauss"))
		{
			start = new double [6];
			start[0] = centroid - spread*0.5;
			start[1] = spread/1.5;
			start[2] = area/2;
			start[3] = centroid + spread*0.5;
			start[4] = spread/1.5;
			start[5] = area/2;
			step = new double [6];
			step[0] = spread/10;
			step[1] = spread/15;
			step[2] = start[2]/10;
			step[3] = spread/10;
			step[4] = spread/15;
			step[5] = start[5]/10;
			reg = new Regression(x, y);
			reg.simplexPlot(r, start, step);
			return reg.getBestEstimates();
		}
		return null;
	}
	/**
	 * This returns the calculated y-values for the function defined by the given parameter list and function name.
	 * @param x	X-axis data
	 * @param params	Function parameter values
	 * @param function	Function name
	 * @return
	 */
	public static double[] getExpectedValues(double[] x, double[] params, String function)
	{
		ACM_CustomFunctions r = ACM_CustomFunctions.getNew(function);
		double[] y = new double [x.length];
		
		
		
		for (int i = 0; i < x.length; i++)
			y[i] = r.function(params, new double [] {x[i]});
		return y;
	}
}
