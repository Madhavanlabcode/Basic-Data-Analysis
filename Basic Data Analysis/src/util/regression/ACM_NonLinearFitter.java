package util.regression;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.analysis.UnivariateFunction;
import org.apache.commons.math3.analysis.solvers.BisectionSolver;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;

import util.ArrayOps;
import util.Complex;
import util.FieldOps;

public class ACM_NonLinearFitter {
	
	static final int MAX_ITERATIONS = 65535;
	
	/**
	 * Fits an (x, y) data set to the function defined as fname.
	 * All of the steps needed to iterate the simplex until the minimum is achieved are carried out.
	 * @param x
	 * @param y
	 * @param fname
	 * @return the set of function parameters corresponding to the best value achieved.
	 */
	public static double[] fitToFunction(double[] x, double[] y, String fname)
	{
		NelderMeadSimplex simp = getSimplex(x, y, fname);
		ChiSquared chisq = new ChiSquared(x, y, fname);
		PointComp comp = new PointComp();
		boolean done = false; //if the simplex holds its value for a greater number of steps than twice the dimension of the parameter space, 
		//it will be regarded as done.
		int nstepsStuck = 0;
		
		double error = 0, lastError = Double.MAX_VALUE;
		
		double[] bestParam;
		
		bestParam = getMinParams(simp.getPoints(), chisq);
		error = chisq.value(bestParam);
		
//		ArrayList<String> logLines = new ArrayList<String>();
		
		int nparams = ACM_CustomFunctions.getNParameters(fname);
		
		int i;
		for (i = 0; i < MAX_ITERATIONS && !done; i++)
		{
			if (lastError == error) nstepsStuck++;
			else nstepsStuck = 0;
			lastError = error;
			simp.iterate(chisq, comp);
			
//			bestParam = getMinParams(simp.getPoints(), comp);
			bestParam = getMinParams(simp.getPoints(), chisq);
			error = chisq.value(bestParam);
			
//			logLines.add("" + i + "\t" + error + "\t" + lastError);
			done = nstepsStuck > nparams*2;
		}
		
		System.out.println(" " + i);

//		String[] log = new String [logLines.size()];
//		for (i = 0; i < log.length; i++)
//			log[i] = logLines.get(i);
//		ColumnIO.writeLines(log, FileOps.selectSave(null).toString());
		return bestParam;
	}
	public static FittingResult fitToFunctionFull(double[] x, double[] y, String fname)
	{
		NelderMeadSimplex simp = getSimplex(x, y, fname);
		ChiSquared chisq = new ChiSquared(x, y, fname);
		PointComp comp = new PointComp();
		boolean done = false; //if the simplex holds its value for a greater number of steps than twice the dimension of the parameter space, 
		//it will be regarded as done.
		int nstepsStuck = 0;
		
		double error = 0, lastError = Double.MAX_VALUE;
		
		double[] bestParam;
		
		bestParam = getMinParams(simp.getPoints(), chisq);
		error = chisq.value(bestParam);
		
		
		int nparams = ACM_CustomFunctions.getNParameters(fname);
		
		int i;
		for (i = 0; i < MAX_ITERATIONS && !done; i++)
		{
			if (lastError == error) nstepsStuck++;
			else nstepsStuck = 0;
			lastError = error;
			simp.iterate(chisq, comp);
			
			bestParam = getMinParams(simp.getPoints(), chisq);
			error = chisq.value(bestParam);
			
			done = nstepsStuck > nparams*2;
		}
		double[] ycalc = getExpectedY(x, bestParam, fname);
		return new FittingResult(bestParam, x, y, ycalc, error);
	}
	public static FittingResult fitToFunctionFull(double[] x, double[] y, ACM_CustomFunctions f, double[] start, double[] step)
	{
		NelderMeadSimplex simp = new NelderMeadSimplex(step);
		simp.build(start);

		ChiSquared chisq = new ChiSquared(x, y, f);
		PointComp comp = new PointComp();
		boolean done = false; //if the simplex holds its value for a greater number of steps than twice the dimension of the parameter space, 
		//it will be regarded as done.
		int nstepsStuck = 0;
		double error = 0, lastError = Double.MAX_VALUE;
		
		double[] bestParam;
		
		bestParam = getMinParams(simp.getPoints(), chisq);
		error = chisq.value(bestParam);
		
		
		int nparams = f.plist.length;
		ArrayList<double[]> bestParamsHist = new ArrayList<double[]>();
		int i;
		for (i = 0; i < MAX_ITERATIONS && !done; i++)
		{
			if (lastError == error) nstepsStuck++;
			else nstepsStuck = 0;
			lastError = error;
			simp.iterate(chisq, comp);
			
			bestParam = getMinParams(simp.getPoints(), chisq);
			error = chisq.value(bestParam);
			
			double[] bestParamTemp = new double [nparams+1];
			for (int k = 0; k < nparams; k++)
				bestParamTemp[k] = bestParam[k];
			bestParamTemp[nparams] = error;
			bestParamsHist.add(bestParamTemp.clone());
			
			done = nstepsStuck > nparams*2;
		}
		double[] ycalc = new double[y.length];
		for (i = 0; i < x.length; i++)
			ycalc[i] = f.function(bestParam, new double[] {x[i]});
		FittingResult ans = new FittingResult(bestParam, x, y, ycalc, error);
		ans.paramHistory = bestParamsHist;
		return ans;
	}
	public static FittingResult fitToFunctionFull(double[][] x, double[] y, ACM_CustomFunctions f, double[] start, double[] step)
	{
		NelderMeadSimplex simp = new NelderMeadSimplex(step);
		simp.build(start);

		ChiSquared chisq = new ChiSquared(x, y, f);
		PointComp comp = new PointComp();
		boolean done = false; //if the simplex holds its value for a greater number of steps than twice the dimension of the parameter space, 
		//it will be regarded as done.
		int nstepsStuck = 0;
		double error = 0, lastError = Double.MAX_VALUE;
		
		double[] bestParam;
		
		bestParam = getMinParams(simp.getPoints(), chisq);
		error = chisq.value(bestParam);
		
		
		int nparams = f.plist.length;
		ArrayList<double[]> bestParamsHist = new ArrayList<double[]>();
		int i;
		for (i = 0; i < MAX_ITERATIONS && !done; i++)
		{
			if (lastError == error) nstepsStuck++;
			else nstepsStuck = 0;
			lastError = error;
			simp.iterate(chisq, comp);
			
			bestParam = getMinParams(simp.getPoints(), chisq);
			error = chisq.value(bestParam);
			
			double[] bestParamTemp = new double [nparams+1];
			for (int k = 0; k < nparams; k++)
				bestParamTemp[k] = bestParam[k];
			bestParamTemp[nparams] = error;
			bestParamsHist.add(bestParamTemp.clone());
			
			done = nstepsStuck > nparams*4;
		}
		System.out.println("There were " + i + " iterations.");
		double[] ycalc = new double[y.length];
		for (i = 0; i < x.length; i++)
			ycalc[i] = f.function(bestParam, x[i]);
		FittingResult ans = new FittingResult(bestParam, x[0], y, ycalc, error);
		ans.paramHistory = bestParamsHist;
		return ans;
	}
	
	/**
	 * The FittingResult uses a parabola y = a(x-h)^2 + k. param[0] is h: the apex of the parabola.
	 * @param x
	 * @param y
	 * @return
	 */
	public static FittingResult fitToParabola(double[] x, double[] y)	{
		int nParams = 3;
		
		double[] abce = ArrayOps.fitToParabola(x, y);
		double[] hak = new double [3];
		
		hak[1] = abce[0];
		if (abce[0] == 0) hak[0] = 0;
		else hak[0] = -abce[1]/(2*abce[0]);
		hak[2] = abce[2]-abce[0]*hak[0]*hak[0];
		
		double[] ycalc = new double [x.length];
		for (int i = 0; i < x.length; i++)
			ycalc[i] = abce[0]*x[i]*x[i] + abce[1]*x[i] + abce[2];
		
		return new FittingResult(hak, x, y, ycalc, abce[3]);
	}
	
	public static double[] fitToFunction(double[][] data, String fname)
	{
		NelderMeadSimplex simp = getSimplex(data, fname);
		ChiSquared chisq = new ChiSquared(data, fname);
		PointComp comp = new PointComp();
		boolean done = false; //if the simplex holds its value for a greater number of steps than twice the dimension of the parameter space, 
		//it will be regarded as done.
		int nstepsStuck = 0;
		
		double error = 0, lastError = Double.MAX_VALUE;
		
		double[] bestParam;
		
		bestParam = getMinParams(simp.getPoints(), chisq);
		error = chisq.value(bestParam);
		
//		ArrayList<String> logLines = new ArrayList<String>();
		
		int nparams = ACM_CustomFunctions.getNParameters(fname);
		
//		ColorScale c = ColorScales.getNew(data, 1);
//		File base = FileOps.selectSave(null);

		
		int i;
		for (i = 0; i < MAX_ITERATIONS && !done; i++)
		{
			if (lastError == error) nstepsStuck++;
			else nstepsStuck = 0;
			lastError = error;
			simp.iterate(chisq, comp);
			
//			bestParam = getMinParams(simp.getPoints(), comp);
			bestParam = getMinParams(simp.getPoints(), chisq);
			error = chisq.value(bestParam);
			
//			SRAW.writeImage(base.toString() + MovieMaker.fromInt(i), ImageEditing.getBufferedImage(getExpectedData(data.length, data[0].length, bestParam, fname), c));//			logLines.add("" + i + "\t" + error + "\t" + lastError);
			done = nstepsStuck > nparams*2;
		}
		
		System.out.println("The total number of iterations was " + i);

//		String[] log = new String [logLines.size()];
//		for (i = 0; i < log.length; i++)
//			log[i] = logLines.get(i);
//		ColumnIO.writeLines(log, FileOps.selectSave(null).toString());
		return bestParam;
	}
	
	
	/**
	 * returns the simplex for an 1d (x, y) array for the given function name
	 * @param x
	 * @param y
	 * @param fname
	 * @return
	 */
	public static NelderMeadSimplex getSimplex(double[] x, double [] y, String fname)
	{
		NelderMeadSimplex simp = null;
		
		double[] start = null, step = null;
		double cent = (ArrayOps.max(x)+ArrayOps.min(x))/2;
		double xp = ArrayOps.max(x), xm = ArrayOps.min(x), ym = ArrayOps.min(y), yp = ArrayOps.max(y);
		double dx = xp-xm;
		double dy = yp-ym;

		double centroid = x[(int)Math.round(ArrayOps.centIndex(y))];
		double spread = ArrayOps.sigmaIndex(y)*(dx/x.length);
		double area = ArrayOps.sum(y)*(dx/x.length);

		if (fname.equalsIgnoreCase("Lorentzian"))
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
		}
		else if (fname.equalsIgnoreCase("TwoLorentzian"))
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
		}
		else if (fname.equalsIgnoreCase("Fermi"))
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
		}
		else if (fname.equalsIgnoreCase("FermiLine"))
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
		}
		else if (fname.equalsIgnoreCase("1PeakLorentzian_Line"))
		{
			start = new double [5];
			start[0] = cent;
			start[1] = dx/2;
			start[2] = dy*start[1]*Math.PI/4;
			start[3] = ym;
			start[4] = 0;
			step = new double [5];
			step[0] = dx/5;
			step[1] = dx/10;
			step[2] = start[2]/10;
			step[3] = dy/20;
			step[4] = (dy/dx)/10;
		}
		else if (fname.equalsIgnoreCase("TwoGauss"))
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
		}
		else if (fname.equalsIgnoreCase("ShiftedExp"))
		{
			start = new double [3];
			start[0] = Math.abs(yp) > Math.abs(ym) ? (yp/ym)/dx : (ym/yp)/dx;
			start[1] = dy;
			start[2] = ym;
			step = new double [3];
			step[0] = start[0]/10;
			step[1] = start[1]/10;
			step[2] = start[2]/10;
		}
		simp = new NelderMeadSimplex(step);
		simp.build(start);
		return simp;
	}
	
	
	/**
	 * returns the simplex for a 2d dataset. Dataset is assumed to be square. Right now the only supported type is the 2d Gaussian
	 * @param data
	 * @param fname
	 * @return
	 */
	public static NelderMeadSimplex getSimplex(double[][] data, String fname)
	{
		NelderMeadSimplex simp;
		
		double[] start, step;
		double[] cent = new double[] {data.length/2, data[0].length/2};
		double xm = 0, xp = data.length, ym = 0, yp = data[0].length;
		double dx = xp-xm;
		double dy = yp-ym;
		
		double zm = FieldOps.min(data);
		double zx = FieldOps.max(data);
		double dz = zx-zm;

//		double[] centroid = CentroidField.centroid(data, 0, data.length, 0, data[0].length, false);
		double[] centroid = new double[] {FieldOps.centXSq(data), FieldOps.centYSq(data)};
		double[] spread = FieldOps.sigmaR(data);
		double area = ArrayOps.sum(data) - data.length*data[0].length*zm;

		if (fname.equalsIgnoreCase("Gauss2D"))
		{
			start = new double [4];
			start[0] = centroid[0];
			start[1] = centroid[1];
			start[2] = Complex.mag(spread)/Math.sqrt(2);
			start[3] = area;
			step = new double [4];
			step[0] = Complex.mag(spread)/10;
			step[1] = Complex.mag(spread)/10;
			step[2] = start[2]/5;
			step[3] = area/10;
			simp = new NelderMeadSimplex(step);
			simp.build(start);
			return simp;
		}
		if (fname.equalsIgnoreCase("Gauss2DConst"))
		{
			start = new double [5];
			start[0] = centroid[0];
			start[1] = centroid[1];
			start[2] = Complex.mag(spread)/Math.sqrt(2);
			start[3] = (zx-zm)*(start[2]*Math.sqrt(2*Math.PI));
			start[4] = zm;
			step = new double [5];
			step[0] = Complex.mag(spread)/10;
			step[1] = Complex.mag(spread)/10;
			step[2] = start[2]/5;
			step[3] = start[3]/10;
			step[4] = dz/10;
			simp = new NelderMeadSimplex(step);
			simp.build(start);
			return simp;
		}
		if (fname.equalsIgnoreCase("Dome"))
		{
			start = new double [5];
			start[0] = centroid[0];
			start[1] = centroid[1];
			start[2] = Complex.mag(spread);
			start[3] = (zx-zm);
			start[4] = zm;
			step = new double [5];
			step[0] = Complex.mag(spread)/10;
			step[1] = Complex.mag(spread)/10;
			step[2] = start[2]/5;
			step[3] = start[3]/10;
			step[4] = dz/10;
			simp = new NelderMeadSimplex(step);
			simp.build(start);
			return simp;
		}
		if (fname.equalsIgnoreCase("CosSqDimple"))
		{
			start = new double [5];
			start[0] = centroid[0];
			start[1] = centroid[1];
			start[2] = Complex.mag(spread)*2;
			start[3] = (zx-zm);
			start[4] = zm;
			step = new double [5];
			step[0] = Complex.mag(spread)/10;
			step[1] = Complex.mag(spread)/10;
			step[2] = start[2]/5;
			step[3] = start[3]/10;
			step[4] = dz/10;
			simp = new NelderMeadSimplex(step);
			simp.build(start);
			return simp;
		}
		return null;
	}

	
	/**
	 * This gives the y values of the fitting.
	 * @param x
	 * @param params
	 * @param fname
	 * @return
	 */
	public static double[] getExpectedY(double[] x, double[] params, String fname)
	{
		ACM_CustomFunctions f = ACM_CustomFunctions.getNew(fname);
		double[] y = new double [x.length];
		
		
		
		for (int i = 0; i < x.length; i++)
			y[i] = f.function(params, new double [] {x[i]});
		return y;
	}
	/**
	 * The expected values of z(x, y) for 2D data.
	 * @param x
	 * @param params
	 * @param fname
	 * @return
	 */
	public static double[][] getExpectedData(int nx, int ny, double[] params, String fname)
	{
		ACM_CustomFunctions f = ACM_CustomFunctions.getNew(fname);
		double[][] data = new double [nx][ny];
		
		
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				data[i][j] = f.function(params, new double [] {i, j});
		return data;
	}
	public static class ChiSquared implements MultivariateFunction
	{
		//the function to be fit
		ACM_CustomFunctions f;
		
		double[][] coordinates; //the coordinates array is [npts][d] where d is the dimension. npts goes as size^d.
		double[] values;	//the array of data values [npts]
		
		int npts, d;
		
		/**
		 * This constructor is for 1d data consisting of ordered pairs (x, y) fit to the function.
		 * @param x
		 * @param y
		 */
		public ChiSquared(double[] x, double[] y, String fname)
		{
			f = ACM_CustomFunctions.getNew(fname);
			npts = x.length;
			coordinates = new double [npts][1];
			values = y.clone();
			
			for (int i = 0; i < npts; i++)
				coordinates[i] = new double[] {x[i]};
		}
		public ChiSquared(double[] x, double[] y, ACM_CustomFunctions f)
		{
			this.f = f;
			npts = x.length;
			coordinates = new double [npts][1];
			values = y.clone();
			
			for (int i = 0; i < npts; i++)
				coordinates[i] = new double[] {x[i]};
		}
		
		/**
		 * This is for 2d data. the coordinates will be entered as [i,j] and the values are dataset[i][j].
		 * @param dataset
		 * @param fname
		 */
		public ChiSquared(double[][] dataset, String fname){
			f = ACM_CustomFunctions.getNew(fname);
			npts = dataset.length*dataset[0].length;
			
			values = new double [npts];
			coordinates = new double [npts][2];
			for (int i = 0; i < dataset.length; i++)
				for (int j = 0; j < dataset.length; j++)
				{
					values[i*dataset.length + j] = dataset[i][j];
					coordinates[i*dataset.length + j] = new double [] {i, j};
				}
		}
		/**
		 * This is for data where x has more than one parameter (e.g. 2, N and B for LL's)
		 * @param x
		 * @param y
		 * @param f
		 */
		public ChiSquared(double[][] x, double[] y, ACM_CustomFunctions f)
		{
			npts = x.length;
			this.f = f;
			
			coordinates = new double [npts][x[0].length];
			for (int i = 0; i < npts; i++)
				coordinates[i] = x[i].clone();
			values = y.clone();
		}
		
		@Override
		/**
		 * returns the sum of squared residuals between the data and the stored function
		 */
		public double value(double[] param) {
			double sum = 0;
			for (int i = 0; i < npts; i++)
				sum += Math.pow(f.function(param, coordinates[i]) - values[i], 2);
			return sum;
		}
		
	}
	
	public static class PointComp implements Comparator<PointValuePair>
	{
		@Override
		public int compare(PointValuePair arg0, PointValuePair arg1) {
			if (arg0.getValue() > arg1.getValue()) return 1;
			else if (arg0.getValue() < arg1.getValue()) return -1;
			return 0;
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
	
	
	public static class FittingResult implements Comparable
	{
		public double[] fitParams;
		public double[] x, y;//input data for storage
		public double[] yCalc;
		public double chsqPerDOF;
		
		public ArrayList<double[]> paramHistory;
		
		public int nParams;
		
		public FittingResult(double[] fitParams, double[] x, double[] y,
				double[] yCalc, double chisq) {
			this.fitParams = fitParams;
			this.x = x.clone();
			this.y = y.clone();
			this.yCalc = yCalc;
			this.chsqPerDOF = chisq/(x.length-1);
			nParams = fitParams.length;
		}
		@Override
		public int compareTo(Object arg0) {
			// TODO Auto-generated method stub
			if (chsqPerDOF < ((FittingResult)arg0).chsqPerDOF)
				return (int)(Math.min(-1, ((FittingResult)arg0).chsqPerDOF - chsqPerDOF));
			else if (chsqPerDOF > ((FittingResult)arg0).chsqPerDOF)
				return (int)(Math.max(1, chsqPerDOF - ((FittingResult)arg0).chsqPerDOF));
			return 0;
		}
		
		public static double[] getParamArray(int i, FittingResult[] r)
		{
			double[] ans = new double [r.length];
			for (int j = 0; j < r.length; j++)
				ans[j] = r[j].fitParams[i];
			return ans;
		}
		public static double[][] getPTable(FittingResult[] r)
		{
			int maxR = 0;
			for (int i = 0; i < r.length; i++)
				maxR = Math.max(maxR, r[i].nParams);
			double[][] ans = new double [maxR][r.length];
			for (int j = 0; j < r.length; j++)
				for (int i = 0; i < r[j].nParams; i++)
					ans[i][j] = r[j].fitParams[i];
			return ans;
		}
		
		
		public int getByteLength()
		{
			return 4 + 4 + 8*(fitParams.length + x.length + y.length + yCalc.length + 1);
		}
		
		/**
		 * Writes to a dataOutputStream Hopefully this is the buffered stream provided
		 * by FittingResult.writeBIN([])
		 * @param outd
		 */
		public void writeBin(DataOutputStream outd)
		{
			try {
				outd.writeInt(getByteLength());
				outd.writeInt(nParams);
				for (int i = 0; i < nParams; i++)
					outd.writeDouble(fitParams[i]);
				for (int i = 0; i < x.length; i++)
					outd.writeDouble(x[i]);
				for (int i = 0; i < y.length; i++)
					outd.writeDouble(y[i]);
				for (int i = 0; i < y.length; i++)
					outd.writeDouble(yCalc[i]);
				outd.writeDouble(chsqPerDOF);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		public static FittingResult readBIN(DataInputStream ind)
		{
			try {
				int bytes = ind.readInt();
				int nParams = ind.readInt();
				int nDoubles = (bytes - 8 - 8 - nParams*8)/3;
				double[] x = new double [nDoubles];
				double[] y = new double [nDoubles];
				double[] yCalc = new double[nDoubles];
				double[] fitParams = new double [nParams];
				for (int i = 0; i < nParams; i++)
					fitParams[i] = ind.readDouble();
				for (int i = 0; i < x.length; i++)
					x[i] = ind.readDouble();
				for (int i = 0; i < y.length; i++)
					y[i] = ind.readDouble();
				for (int i = 0; i < y.length; i++)
					yCalc[i] = ind.readDouble();
				double chisq = ind.readDouble();
				return new FittingResult(fitParams, x, y, yCalc, chisq);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			return null;
		}
		public static void writeBIN(FittingResult[] results, String filepath)
		{
			int length = 4;
			for (int i = 0; i < results.length; i++)
				length += results[i].getByteLength();
			
			FileOutputStream outf = null;
			BufferedOutputStream outbuff = null;
			DataOutputStream outd = null;
			
			try {
				outf = new FileOutputStream(filepath);
				outbuff = new BufferedOutputStream(outf, length);
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			outd = new DataOutputStream(outbuff);
			try {
				outd.writeInt(results.length);
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			for (int i = 0; i < results.length; i++)
				results[i].writeBin(outd);
			try {
//				outbuff.flush();
//				outf.flush();
				outd.close();
				outbuff.close();
				outf.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		public static FittingResult[] readBIN(String filepath)
		{
			File file = new File(filepath);
			FileInputStream inf = null;
			BufferedInputStream inbuff = null;
			DataInputStream ind = null;
			try {
				inf = new FileInputStream(file);
				inbuff = new BufferedInputStream(inf, (int)file.length());
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				return null;
			}
			ind = new DataInputStream(inbuff);
			int n = 0;
			try {
				n = ind.readInt();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			FittingResult[] ans = new FittingResult[n];
			for (int i = 0; i < n; i++)
				ans[i] = readBIN(ind);
		try{
			ind.close();
			inbuff.close();
			inf.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return ans;

			
		}
		public double getYAt(double x, String fname) {
			return ACM_CustomFunctions.getNew(fname).function(fitParams, new double [] {x});
		}
	}		
	/**
	 * Finds the actual maximum of a 1PeakLorentzianLine(if any) by setting the derivative equal to zero. The maximum is assumed to lie within +/- 10 times the width.
	 * @param f
	 * @return
	 */
	public static double getPeakMaximum(FittingResult f)
	{
		
		double H = f.fitParams[1]*f.fitParams[2]/Math.PI;
		double g2 = f.fitParams[1]*f.fitParams[1];
		double m = f.fitParams[4];

		if (m == 0) return f.fitParams[0];
		
		BisectionSolver bs = new BisectionSolver();
		Quartic dLdchi = new Quartic(1, 0, g2/2, -H/m, g2*g2/16);
		double upperValue = 0;
		double lowerValue = 0;
		try{
//			lowerValue = bs.solve(MAX_ITERATIONS, dLdchi, -100*f.fitParams[1], 0, -99*f.fitParams[1], null);
			lowerValue = bs.solve(MAX_ITERATIONS, dLdchi, -100*f.fitParams[1], 0, -99*f.fitParams[1]);
		}
		catch (Exception tmme)
		{
			tmme.printStackTrace();
		}
		try{
			upperValue = bs.solve(MAX_ITERATIONS, dLdchi, 0, 100*f.fitParams[1], 99*f.fitParams[1]);
		}
		catch (Exception tmme)
		{
			tmme.printStackTrace();
		}
		
		if (upperValue == Double.NaN && lowerValue == Double.NaN){
			System.out.println("Solver returned NaN");
			return f.fitParams[0];
		}
		else if (upperValue != Double.NaN && lowerValue == Double.NaN)
			return upperValue + f.fitParams[0];
		else if (upperValue == Double.NaN && lowerValue != Double.NaN)
			return lowerValue + f.fitParams[0];
		else
			return (f.getYAt(upperValue, "1PeakLorentzian_Line") > f.getYAt(lowerValue, "1PeakLorentzian_Line")) ? upperValue + f.fitParams[0] : lowerValue + f.fitParams[0]; 
		
	}
	
	private static class Quartic implements UnivariateFunction
	{
		double a, b, c, d, e;
		public Quartic(double a, double b, double c, double d, double e) {
			super();
			this.a = a;
			this.b = b;
			this.c = c;
			this.d = d;
			this.e = e;
		}
		@Override
		public double value(double x) {
			return a*x*x*x*x + b*x*x*x + c*x*x + d*x + e;
		}
		
	}


}
