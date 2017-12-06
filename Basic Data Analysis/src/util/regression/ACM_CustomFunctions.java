package util.regression;

import util.Printer;

/**
 * This is supposed to contain wrappers for fitting custom functions in the Apache Commons Math library.
 * @author madhavanlab2011
 *
 */
public abstract class ACM_CustomFunctions {
	
	public String[] plist;
	//in all of these the peak center should be [0] so that a satter plot table can be written standardly.
	
	public abstract double function(double[] param, double[] x);

	public static class TwoLorentzian  extends ACM_CustomFunctions{
		/*parameters p:
			0 = center, 1
			1 = width, 1
			2 = height, 1
			3 = center, 2
			4 = width, 2
			5 = height, 2
			6 = additive constant
		 */
		public TwoLorentzian()
		{
			super();
			plist = new String[] {"Center_1", "Width_1", "Height_1", "Center_2", "Width_2", "Height_2", "Additive Const"};
		}
		public double function(double[] param, double[] x) {
			// TODO Auto-generated method stub
			return param[6] + (param[2]/Math.PI)*((param[1]/2)/((x[0]-param[0])*(x[0]-param[0]) + (param[1]*param[1]/4))) +
					(param[5]/Math.PI)*((param[4]/2)/((x[0]-param[3])*(x[0]-param[3]) + (param[4]*param[4]/4)));
		}
		
	}
	public static class LorentzianPlusConst extends ACM_CustomFunctions{
		/*parameters p:
			0 = center
			1 = additive constant
			2 = width
			3 = height
		 */
		public LorentzianPlusConst()
		{
			super();
			plist = new String[] {"Center", "Additive Const", "Width", "Height"};
		}
		public double function(double[] param, double[] x) {
			// TODO Auto-generated method stub
			return param[1] + (param[3]/Math.PI)*((param[2]/2)/((x[0]-param[0])*(x[0]-param[0]) + (param[2]*param[2]/4)));
		}
		
	}
	public static class LorentzianPlusLine extends ACM_CustomFunctions{
		/*parameters p:
			0 = center
			1 = width
			2 = height
		 	3 = offset
		 	4 = slope
		 */
		public LorentzianPlusLine()
		{
			super();
			plist = new String[] {"Center", "Width", "Height", "Offset", "Slope"};
		}
		public double function(double[] param, double[] x) {
			// TODO Auto-generated method stub
			return param[3]+param[4]*x[0] + (param[2]/Math.PI)*((param[1]/2)/((x[0]-param[0])*(x[0]-param[0]) + (param[1]*param[1]/4)));
		}
		
	}
	public static class LorentzianPlusLine_NPeaks extends ACM_CustomFunctions{
		/*parameters p:
			0 = center
			1 = width
			2 = height
			...
		 	n*3 = offset
		 	n*3 + 1 = slope
		 */
		
		final int n;
		
		public LorentzianPlusLine_NPeaks(int n)
		{
			super();
			this.n = n;
			plist = new String[3*n+2];
			for (int i = 0; i < n; i++)
			{
				plist[3*i] = "Center " + i;
				plist[3*i+1] = "Width " + i;
				plist[3*i+2] = "Height " + i;
			}
			plist[plist.length-1] = "Slope";
			plist[plist.length-1] = "Offset";
		}
		public double function(double[] param, double[] x) {
			// TODO Auto-generated method stub
			double ans = param[3*n]+param[3*n+1]*x[0];
			for (int i = 0; i < n; i++)
				ans += (param[3*i+2]/Math.PI)*((param[3*i+1]/2)/((x[0]-param[3*i])*(x[0]-param[3*i]) + (param[3*i+1]*param[3*i+1]/4)));
			return ans;
		}
		
	}
	
	public static class Fermi extends ACM_CustomFunctions{
		
		public Fermi()
		{
			super();
			plist = new String[] {"Center", "Temperature", "Height", "Additive Const"};
		}
		public double function(double[] param, double[] x) {
			// TODO Auto-generated method stub
			double p = (x[0] - param[0])/param[1];
			if (x[0] < -0.514){
				System.out.print(Printer.arrayLnHorizontal(param));
			}
			if (p > 100)
			{
				p = 100;
			}
			if (p < -100) {
				p = -100; 
			}
			return param[3] + param[2]*(1/(1+Math.exp(p)));
		}
		
	}
	public static class FermiLine extends ACM_CustomFunctions{
		
		public FermiLine()
		{
			super();
			plist = new String[] {"Center", "Temperature", "Height", "Additive Const", "Slope"};
		}
		public double function(double[] param, double[] x) {
			// TODO Auto-generated method stub
			double p = (x[0] - param[0])/param[1];
//			if (x[0] < -0.514){
//				System.out.print(Printer.arrayLnHorizontal(param));
//			}
//			if (p > 100)
//			{
//				p = 100;
//			}
//			if (p < -100) {
//				p = -100; 
//			}
			return param[3] + param[2]*(1/(1+Math.exp(p))) + param[4]*x[0];
		}
		
	}
	
	public static class TwoGauss extends ACM_CustomFunctions{

		public TwoGauss()
		{
			super();
			plist = new String[] {"Mean1", "Sigma1", "Weight1","Mean2", "Sigma2","Weight2"};
		}
		
		@Override
		public double function(double[] param, double[] x) {
			return param[2]*normalDist(x[0], param[0], param[1]) + param[5]*normalDist(x[0], param[3], param[4]);
		}
	}
	
	public static class Gauss_2D extends ACM_CustomFunctions{

		public Gauss_2D()
		{
			super();
			plist = new String[] {"Mean_x", "Mean_y", "Sigma", "Weight"};
		}
		@Override
		public double function(double[] param, double[] x) {
			return (param[3]/(param[2]*Math.sqrt(2*Math.PI))) * Math.exp(-(Math.pow((x[0]-param[0])/param[2], 2) + Math.pow((x[1]-param[1])/param[2], 2))/2);
		}
		
	}
	public static class Gauss_2D_Off extends ACM_CustomFunctions{

		public Gauss_2D_Off()
		{
			super();
			plist = new String[] {"Mean_x", "Mean_y", "Sigma", "Weight", "Offset"};
		}
		@Override
		public double function(double[] param, double[] x) {
			return (param[3]/(param[2]*Math.sqrt(2*Math.PI))) * Math.exp(-(Math.pow((x[0]-param[0])/param[2], 2) + Math.pow((x[1]-param[1])/param[2], 2))/2) + param[4];
		}
		
	}
	
	/**
	 * A hemisphere elongated in the z-direction, rising from a constant background.
	 * Codename "Dome".
	 * @author madhavanlab2011
	 *
	 */
	public static class Dome_2D extends ACM_CustomFunctions{
		public Dome_2D()
		{
			super();
			plist = new String[] {"Cent_x", "Cent_y", "Radius", "Height", "Offset"};
		}
		@Override
		public double function(double[] param, double[] x) {
			double rsq = (x[0]-param[0])*(x[0]-param[0]) + (x[1]-param[1])*(x[1]-param[1]);
			double hsq = (param[2]*param[2] - rsq)/(param[2]*param[2]);
			if (hsq <= 0) return param[4];
			else return param[4] + param[3]*Math.sqrt(hsq);
			
		}
	}
	
	/**
	 * A dome that looks like cos((pi/2)*(r/radius))^2 for r < radius, and 0 for r > radius.
	 * Codename CosSqDimple 
	 * @author madhavanlab2011
	 *
	 */
	public static class CosSq_2D_Dimple extends ACM_CustomFunctions
	{
		public CosSq_2D_Dimple()
		{
			super();
			plist = new String[] {"Cent_x", "Cent_y", "Radius", "Height", "Offset"};
		}

		@Override
		public double function(double[] param, double[] x) {
			// TODO Auto-generated method stub
			double r = Math.sqrt((x[0]-param[0])*(x[0]-param[0]) + (x[1]-param[1])*(x[1]-param[1]));
			if (r > param[2]) return param[4];
			else return param[4] + param[3] * Math.pow(Math.cos((Math.PI/2)*(r/param[2])), 2);
		}
		
	}
	public static class ShiftedExponential extends ACM_CustomFunctions
	{
		public ShiftedExponential()
		{
			super();
			plist = new String[] {"DecayConst", "FrontFactor", "AdditiveConst"};
		}

		@Override
		public double function(double[] param, double[] x) {
			return param[1]*Math.exp(param[0]*x[0]) + param[2];
		}
	}

	public static ACM_CustomFunctions getNew(String name)
	{
		if (name.equalsIgnoreCase("Lorentzian"))
			return new LorentzianPlusConst();
		if (name.equalsIgnoreCase("TwoLorentzian"))
			return new TwoLorentzian();
		if (name.equalsIgnoreCase("Fermi"))
			return new Fermi();
		if (name.equalsIgnoreCase("FermiLine"))
			return new FermiLine();
		for (int i = 1; i < 5; i++)
			if(name.equalsIgnoreCase("" + i + "PeakLorentzian_Line"))
				return new LorentzianPlusLine_NPeaks(i);
		if (name.equalsIgnoreCase("TwoGauss"))
			return new TwoGauss();
		if (name.equalsIgnoreCase("Gauss2D"))
			return new Gauss_2D();
		if (name.equalsIgnoreCase("Gauss2DConst"))
			return new Gauss_2D_Off();
		if (name.equalsIgnoreCase("Dome"))
			return new Dome_2D();
		if (name.equalsIgnoreCase("CosSqDimple"))
			return new CosSq_2D_Dimple();
		if (name.equalsIgnoreCase("ShiftedExp"))
			return new ShiftedExponential();
		
		return null;
	}
	
	public static int getNParameters(String name)
	{
		if (name.equalsIgnoreCase("Lorentzian"))
			return 4;
		if (name.equalsIgnoreCase("TwoLorentzian"))
			return 5;
		if (name.equalsIgnoreCase("Fermi"))
			return 4;
		if (name.equalsIgnoreCase("FermiLine"))
			return 5;
		for (int i = 1; i < 5; i++)
			if(name.equalsIgnoreCase("" + i + "PeakLorentzian_Line"))
				return i*3 + 2;
		if (name.equalsIgnoreCase("TwoGauss"))
			return 6;
		if (name.equalsIgnoreCase("Gauss2D"))
			return 4;
		if (name.equalsIgnoreCase("Gauss2DConst"))
			return 5;
		if (name.equalsIgnoreCase("Dome"))
			return 5;
		if (name.equalsIgnoreCase("CosSqDimple"))
			return 5;
		if (name.equalsIgnoreCase("ShiftedExp"))
			return 3;
		return -1;
		
	}
	public static double[] getExpectedValues(double[] x, double[] params, String function)
	{
		ACM_CustomFunctions r = ACM_CustomFunctions.getNew(function);
		double[] y = new double [x.length];
		
		
		
		for (int i = 0; i < x.length; i++)
			y[i] = r.function(params, new double [] {x[i]});
		return y;
	}
	public static double parabola(double[] param, double x)
	{
		return param[1]*(x-param[0])*(x-param[0]) + param[2];
	}

	public static double normalDist(double x, double mean, double sigma)
	{
		return Math.exp(-Math.pow((x-mean)/sigma,2)/2)/(sigma*Math.sqrt(2*Math.PI));
	}
	


}
