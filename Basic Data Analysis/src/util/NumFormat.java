package util;

public class NumFormat {

	public static String scientific(double x, int nSigFigs)
	{
		double logd = Math.log10(Math.abs(x));
		int log = logd >= 0 ? (int)logd : (int)(logd - 1);
		double x1 = x/(Math.pow(10, log));
		if (x != 0 && !(Math.abs(log) < 3))
			return String.format("%." + (nSigFigs - 1) + "f", x1) + "E" + log;
		else
			return String.format("%." + (nSigFigs - 1) + "f", x);
	}
	
	public static String voltage(double v)
	{
		if (v < 1)
			return "" + FieldOps.round(v*1000) + " mV";
		else return String.format("%.2f", v) + " V";
	}
}
