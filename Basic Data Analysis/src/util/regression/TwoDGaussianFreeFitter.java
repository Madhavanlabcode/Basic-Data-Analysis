package util.regression;

import java.util.ArrayList;

import impurity.PointImp;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;

import util.FieldOps;
import util.regression.ACM_NonLinearFitter.PointComp;

public class TwoDGaussianFreeFitter {
	static final int MAX_ITERATIONS = 65535;
	
	public double[][] fftmag;
	public double[][] values;	//The two-dimensional function which stores the value of the fitting function

	int nx, ny;
	ChiSquared chisq;
	public int userStage = 0;
	
	public double[] initGuesses = new double[6];
	public double[] fitToFunction(boolean useWholeMap)
	{
		FFTFittingCustomFunctions.TwoDGaussianFreeXY func = new FFTFittingCustomFunctions.TwoDGaussianFreeXY();
		chisq = new ChiSquared(func, initGuesses, fftmag, useWholeMap);
		NelderMeadSimplex simp = getSimplex(initGuesses);
		PointComp comp = new PointComp();
		boolean done = false; //if the simplex holds its value for a greater number of steps than twice the dimension of the parameter space, 
		//it will be regarded as done.
		int nstepsStuck = 0;
		
		double error = 0, lastError = Double.MAX_VALUE;
		
		double[] bestParam;
		
		bestParam = getMinParams(simp.getPoints(), chisq);
		error = chisq.value(bestParam);
		
//		ArrayList<String> logLines = new ArrayList<String>();
		
		int nparams = 5;
		
//		ColorScale c = ColorScales.getNew(data, 1);
//		File base = FileOps.selectSave(null);

		System.out.print(chisq.pointsToVisit.size() + " ");
		int i;
		int n = 0;
		for (i = 0; i < MAX_ITERATIONS && !done; i++)
		{
			if (lastError == error) nstepsStuck++;
			else nstepsStuck = 0;
			lastError = error;
			simp.iterate(chisq, comp);
			
			bestParam = getMinParams(simp.getPoints(), chisq);
			error = chisq.value(bestParam);
			
//			SRAW.writeImage(base.toString() + MovieMaker.fromInt(i), ImageEditing.getBufferedImage(getExpectedData(data.length, data[0].length, bestParam, fname), c));//			logLines.add("" + i + "\t" + error + "\t" + lastError);
			done = nstepsStuck > nparams*2;
			n = i;
		}
		
		for (i = 0; i < bestParam.length; i++)
			System.out.print(" " + String.format("%.2f", bestParam[i]) + "\t");
		System.out.println("The total number of iterations was " + n + "  " + error);

//		String[] log = new String [logLines.size()];
//		for (i = 0; i < log.length; i++)
//			log[i] = logLines.get(i);
//		ColumnIO.writeLines(log, FileOps.selectSave(null).toString());
		
		values = new double [nx][ny];
		for (i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				values[i][j] = func.function(bestParam, new double[] {i, j});
			}
		
		return bestParam;
	}

	/**
	 * We assume the lattice is already symmetrized, so that 90 degrees or 60 degrees (or 120) is the angle between the vectors.
	 */
	public TwoDGaussianFreeFitter(double[][] fftmag)
	{
		this.fftmag = fftmag;
		nx = fftmag.length;
		ny = fftmag[0].length;
		//In either case, we form the two lattices, on and off directions, by taking first a unit vector in the direction of a, and secondly that unit vector rotated by
		//the necessary angle over two.
	}
	public NelderMeadSimplex getSimplex(double[] initGuesses)
	{
		NelderMeadSimplex simp;
		
		double[] step;
		
		step = new double[initGuesses.length];
		step[0] = initGuesses[2]/3;
		step[1] = initGuesses[3]/3;
		for (int i = 2; i < initGuesses.length; i++)
		{
			step[i] = initGuesses[i]/3;
		}
		simp = new NelderMeadSimplex(step);
		simp.build(initGuesses);
		return simp;
	}

	public static class ChiSquared implements MultivariateFunction
	{
		//the function to be fit
		FFTFittingCustomFunctions f;
		double[][] data; //The two-dimensional function to be fitted
		double[][] weights;
		
		boolean periodic = false;
		boolean useWholeMap = false;
		ArrayList<int[]> pointsToVisit;
		
		static double widthConst_Weight = 10000000;
		static double widthConst_Bool = 2;
		
		
		
		int nx, ny;
		
		/**
		 * This is for 2d data. the coordinates will be entered as [i,j] and the values are dataset[i][j].
		 * 
		 * With the boolean we specify the direction of the cut. This will determine the mask used on the data.
		 * Note that k and width guesses are taken in pixel units
		 * @param dataset
		 * @param fname
		 */
		public ChiSquared(FFTFittingCustomFunctions f, double[] initGuesses, double[][] data, boolean useWholeMap){
			this.f = f;
			this.data = data;
			nx = data.length;
			ny = data[0].length;
			this.weights = new double [nx][ny];
			this.useWholeMap = useWholeMap;
			pointsToVisit = new ArrayList<int[]>();
			//Now, populate the weights array. The trick is to use periodic boundary conditions so that the edge of the map is no obstacle.
			int ireal, jreal;
			double[] kCent = new double[2];
			kCent[0] = initGuesses[0];
			kCent[1] = initGuesses[1];
			int ki = FieldOps.round(kCent[0]);
			int kj = FieldOps.round(kCent[1]);
			double xW = initGuesses[2]*widthConst_Weight;
			double yW = initGuesses[3]*widthConst_Weight;
			double xwb = initGuesses[2]*widthConst_Bool;
			double ywb = initGuesses[3]*widthConst_Bool;
			int nptsToVisit = 0;
			double weightedDistSq;
			double boolDistSq;
			if (!useWholeMap){
				if (periodic)
					for (int i = -nx/2; i < nx/2; i++)
						for (int j = -ny/2; j < ny/2; j++)
						{
							ireal = (ki + i + nx) % nx;
							jreal = (kj + j + ny) % ny;
							weightedDistSq = ( (i*i)/(xW*xW) + (j*j)/(yW*yW));
							boolDistSq = ( (i*i)/(xwb*xwb) + (j*j)/(ywb*ywb));
							weights[ireal][jreal] = Math.exp(-weightedDistSq);
							if (boolDistSq <= 1)
								pointsToVisit.add(new int[] {ireal, jreal});
						}
				else
					for (int i = -nx/2; i < nx/2; i++)
						for (int j = -ny/2; j < ny/2; j++)
						{
							ireal = (ki + i);
							jreal = (kj + j);
							weightedDistSq = ( (i*i)/(xW*xW) + (j*j)/(yW*yW));
							boolDistSq = ( (i*i)/(xwb*xwb) + (j*j)/(ywb*ywb));
							if (FieldOps.withinBounds(ireal, jreal, data)){
								weights[ireal][jreal] = Math.exp(-weightedDistSq);
								if (boolDistSq <= 1)
									pointsToVisit.add(new int[] {ireal, jreal});
							}
						}
			}
			
//			Layer.writeBIN(Layer.getFreeLayer(weights));
		}
		
		@Override
		/**
		 * returns the sum of squared residuals between the data and the stored function
		 */
		public double value(double[] param) {
			double sum = 0;
			double[] r;
			double[] rSym;
//			System.out.println(Printer.vectorP(param));
//			for (int i = 0; i < nx; i++)
//				for (int j = 0; j < ny; j++)
//				{
//					//First, set the value of r.
//					r = new double[] {i, j};
//					sum += weights[i][j]*Math.pow(f.function(param, r) - data[i][j], 2);
//				}
			int i, j;
			if (!useWholeMap)
				for (int k = 0; k < pointsToVisit.size(); k++)
					{
						//First, set the value of r.
						i = pointsToVisit.get(k)[0];
						j = pointsToVisit.get(k)[1];
						r = new double[] {i, j};
						sum += weights[i][j]*Math.pow(f.function(param, r) - data[i][j], 2);
					}
			else
				for (i = 0; i < nx; i++)
					for (j = 0; j < ny; j++)
					{
						//First, set the value of r.
						r = new double[] {i, j};
						sum += Math.pow(f.function(param, r) - data[i][j], 2);
					}
				
			//			System.out.print(sum + "\t" + Printer.arrayLnHorizontal(param));
			return sum;
		}
		
	}
	public static double[] getMinParams(PointValuePair[] points, ChiSquared chisq)
	{
		int mini = 0;
		for (int i = 0; i < points.length; i++)
		{
//			if (comp.compare(points[mini], points[i]) > 0)
			if (chisq.value(points[i].getPoint()) < chisq.value(points[mini].getPoint()))
				mini = i;
		}
		return points[mini].getPoint();
	}
	
	/**
	 * This is intended to fit impurities to gaussians using the impurity's location as a guess and the indicated width,
	 * in the usual data field.
	 * 
	 * return field [0] is a list of imps [imps.length][nparams]
	 * return field [1] is the sum of all values [nx][ny]
	 * return field [2] is the sum of all guessing weights [nx][ny];
	 * @param imps
	 * @param data
	 * @return
	 */
	public static double[][][] fitListOfImpurities(PointImp[] imps, double[][] data, double widthGuess, boolean expectBrightPeak)
	{
		double[][] ans = new double[imps.length][];
		TwoDGaussianFreeFitter f = new TwoDGaussianFreeFitter(data);
		double heightGuess = FieldOps.sigma(data) * (expectBrightPeak ? 1 : -1);
		double background = expectBrightPeak ? FieldOps.min(data) : FieldOps.max(data);
		System.out.println();
		
		int nx = data.length, ny = data[0].length;
		double[][] values = new double [nx][ny];
		double[][] weights = new double [nx][ny];
//		double[][] weightedValues = new double [nx][ny];
		for (int i = 0; i < imps.length; i++)
		{			
			System.out.print(i + " ");

			f.initGuesses = new double [] {imps[i].pixelPos[0], imps[i].pixelPos[1], widthGuess, widthGuess, heightGuess, background};
			ans[i] = f.fitToFunction(false);
			
			for (int j = 0; j < nx; j++)
				for (int k = 0; k < ny; k++){
					values[j][k] += f.values[j][k];
//					weightedValues[j][k] += f.values[j][k]*f.chisq.weights[j][k];
				}
			for (int p = 0; p < f.chisq.pointsToVisit.size(); p++)
			{
				int j = f.chisq.pointsToVisit.get(p)[0];
				int k = f.chisq.pointsToVisit.get(p)[1];
				weights[j][k] += f.chisq.weights[j][k];
			}
		}
		return new double[][][] {ans, values, weights};
	}
	

}
