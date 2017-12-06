package util.scalar.functions;

import util.scalar.Function;

public class AdHocFunctions {

	public static class LogisticMap implements Function
	{
		double lambda;
			
		public LogisticMap(double lambda) {
			this.lambda = lambda;
		}


		public double of(double x) {
			return 4*lambda*x*(1-x);
		}
	}
	
	public static class gLogistic implements Function
	{
		int n;
		LogisticMap m;
		double temp;
		public gLogistic(int n, LogisticMap m)
		{
			this.n = n;
			this.m = m;
		}
		
		public double of(double x)
		{
			temp = x;
			for (int i = 0; i < n; i++)
				temp = m.of(temp);
			return temp;
		}
	}
	public static class Lambdak implements Function
	{
		double g;

		public Lambdak(double g) {
			super();
			this.g = g;
		}

		public double of(double x) {
			return Math.sqrt(1 + g*g + 2*g*Math.cos(x));
		}
		
		
	}
	public static class LOfN implements Function
	{
		double g;
		int nsteps = 100000;
		double dk = Math.PI/nsteps;
		public Lambdak lambda;
		public LOfN(double g) {
			this.g = g;
			this.lambda = new Lambdak(g);
		}

		public double of(double x) {
			double sum = 0;
			double k;
			for (int i = 0; i < nsteps; i++)
			{
				k = i*dk + dk/2;
				sum += Math.cos(k*x)*dk/lambda.of(k);
			}
			return sum/Math.PI;
		}
		
		
	}
	public static class GOfN implements Function
	{
		double g;
		int nsteps = 10000;
		double dk = Math.PI/nsteps;
		public Lambdak lambda;
		public LOfN ln;
		public GOfN(double g) {
			this.g = g;
			this.lambda = new Lambdak(g);
			ln = new LOfN(g);
		}

		public double of(double x) {
			return ln.of(x) + g*ln.of(x+1);
		}
		
		
	}
}
