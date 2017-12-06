package util.scalar;

import util.Constants;

public interface Function {
	
	public double of(double x);
	
	static class BlackBodyLambda implements Function
	{
		static final double hc = Constants.h * Constants.c;
		static final double eightpihc = (8 * Math.PI * hc);
		
		double kt;
		
		public BlackBodyLambda(double temperature)
		{
			this.kt = temperature * Constants.kB;
		}
		
		public double of(double x)
		{
			return (eightpihc/Math.pow(x, 5)) * 1/(Math.exp(hc/(x*kt)));
		}
	}
	static class Constant implements Function
	{
		double a;
		public Constant(double a)
		{
			this.a = a;
		}
		public double of(double x)
		{
			return a;
		}
	}
	static class Step implements Function
	{
		double h;
		public Step(double h)
		{
			this.h = h;
		}
		
		public double of(double x)
		{
			if (x < 0) return 0;
			else if (x > 0) return h;
			else return h/2;
		}
		
	}
}
