package util.regression;

public abstract class FFTFittingCustomFunctions {
	
	public abstract double function(double[] param, double[] x);
	public String[] plist;

	public static class TwoDGaussian  extends FFTFittingCustomFunctions {
		/*parameters p:
			0 = center K
			1 = width in line
			2 = transverse width
			3 = height
			4 = additive constant
		 */
		public TwoDGaussian()
		{
			super();
			plist = new String[] {"Center K", "Width_In_Line", "Width Transverse", "Height", "Background"};
		}
		public double function(double[] param, double[] x) {
			// TODO Auto-generated method stub
			return param[4] + param[3]*Math.exp(-((Math.pow(x[0]-param[0], 2)/(param[1]*param[1]) + (x[1]*x[1])/(param[2]*param[2]))));
		}
	}
	
	public static class TwoDGaussianFreeXY  extends FFTFittingCustomFunctions {
		/*parameters p:
			0 = center x
			1 = center y
			2 = sigma_x
			3 = sigma_y
			4 = amplitude
			5 = additive constant
		 */
		public TwoDGaussianFreeXY()
		{
			super();
			plist = new String[] {"Center X", "Center Y", "Sigma X", "Sigma Y", "Height", "Background"};
		}
		public double function(double[] param, double[] x) {
			// TODO Auto-generated method stub
			return param[5] + param[4]*Math.exp(-((Math.pow(x[0]-param[0], 2)/(param[2]*param[2]) + (Math.pow(x[1]-param[1], 2))/(param[3]*param[3]))));
		}
	}

}
