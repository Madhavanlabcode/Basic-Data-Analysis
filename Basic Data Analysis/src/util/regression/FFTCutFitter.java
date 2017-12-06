package util.regression;

import java.util.ArrayList;


import javax.swing.JFileChooser;

import org.apache.commons.math3.analysis.MultivariateFunction;
import org.apache.commons.math3.optim.PointValuePair;
import org.apache.commons.math3.optim.nonlinear.scalar.noderiv.NelderMeadSimplex;

import util.FieldOps;
import util.fileops.ColumnIO;
import util.fileops.FileOps;
import util.fileops.Topomap;
import util.geom.AtomicCoordinatesSet;
import util.geom.Distance;
import util.matrix.Matrix;
import util.regression.ACM_NonLinearFitter.PointComp;

public class FFTCutFitter {

	static final int MAX_ITERATIONS = 65535;
	
	int LATTICE_TYPE;
	static final int TRIANGULAR = 0;
	static final int SQUARE = 1;
	double halfAngle;
	
	double[][][] fftmag;
	int nx, ny;
	AtomicCoordinatesSet latt;

	public AtomicCoordinatesSet onLatt, offLatt;

	boolean[][][] angleMask;
	
	public UserGeneratedFFTFit user = null;
	
	public double[] fitToFunction(int index, double kguess, double widthguess, boolean onTrue)
	{
		ChiSquared chisq = new ChiSquared(index, this, new FFTFittingCustomFunctions.TwoDGaussian(), onTrue, kguess, widthguess);
		NelderMeadSimplex simp = getSimplex(fftmag[index], kguess, widthguess, chisq.mask);
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

		
		int i;
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
		}
		
		System.out.println("The total number of iterations was " + i);

//		String[] log = new String [logLines.size()];
//		for (i = 0; i < log.length; i++)
//			log[i] = logLines.get(i);
//		ColumnIO.writeLines(log, FileOps.selectSave(null).toString());
		return bestParam;
	}

	/**
	 * We assume the lattice is already symmetrized, so that 90 degrees or 60 degrees (or 120) is the angle between the vectors.
	 */
	public FFTCutFitter(double[][][] fftmag, AtomicCoordinatesSet symLatt)
	{
		this.fftmag = fftmag;
		nx = fftmag[0].length;
		ny = fftmag[0][0].length;
		this.latt = symLatt;
		double angle = latt.getAngleBetweenVectors();
		if (angle == Math.PI/2) { LATTICE_TYPE = SQUARE; halfAngle = Math.PI/4;}
		else	{LATTICE_TYPE = TRIANGULAR; halfAngle = Math.PI/6;}
		
		//In either case, we form the two lattices, on and off directions, by taking first a unit vector in the direction of a, and secondly that unit vector rotated by
		//the necessary angle over two.
		double[] ua = Distance.unitVector(latt.getA());
		onLatt = new AtomicCoordinatesSet(ua, Matrix.rotateVector(ua, Math.PI/2), new double[] {nx/2, ny/2});
		ua = Matrix.rotateVector(ua, halfAngle);
		offLatt = new AtomicCoordinatesSet(ua, Matrix.rotateVector(ua, Math.PI/2), new double[] {nx/2, ny/2});
		
		//Now, we deploy the boolean field to cut off only the strip
		angleMask = new boolean [2][nx][ny];
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				angleMask[0][i][j] = onLatt.getAngleBetween(new double[] {i-nx/2, j-ny/2}, 0) <= halfAngle; 
				angleMask[1][i][j] = offLatt.getAngleBetween(new double[] {i-nx/2, j-ny/2}, 0) <= halfAngle; 
			}
	}
	
	public NelderMeadSimplex getSimplex(double[][] data, double kguess, double widthguess, boolean[][] mask)
	{
		NelderMeadSimplex simp;
		
		double[] start, step;
		
		double zm = FieldOps.min(data, mask);
		double zx = FieldOps.max(data, mask);
		double dz = zx-zm;

		start = new double [5];
		start[0] = kguess;
		start[1] = widthguess;
		start[2] = widthguess;
		start[3] = (zx-zm);
		start[4] = zm;
		step = new double [5];
		step[0] = widthguess/2;
		step[1] = widthguess/3;
		step[2] = widthguess/3;
		step[3] = start[3]/4;
		step[4] = dz/4;
		simp = new NelderMeadSimplex(step);
		simp.build(start);
		return simp;
	}

	
	
	public static class ChiSquared implements MultivariateFunction
	{
		//the function to be fit
		FFTFittingCustomFunctions f;
		FFTCutFitter parent;
		double[][] fftmag; //The two-dimensional function to be fitted
		double[][] values;	//The two-dimensional function which stores the value of the fitting function
		boolean[][] mask;
		
		static double widthConst = 3;
		
		AtomicCoordinatesSet latt;
		
		int nx, ny;
		
		boolean onTrue;
		
		/**
		 * This is for 2d data. the coordinates will be entered as [i,j] and the values are dataset[i][j].
		 * 
		 * With the boolean we specify the direction of the cut. This will determine the mask used on the data.
		 * Note that k and width guesses are taken in pixel units
		 * @param dataset
		 * @param fname
		 */
		public ChiSquared(int index, FFTCutFitter parent, FFTFittingCustomFunctions f, boolean onTrue, double kguess, double widthGuess){
			this.f = f;
			this.fftmag = parent.fftmag[index];
			nx = fftmag.length;
			ny = fftmag[0].length;
			this.mask = new boolean [nx][ny];
			this.values = new double [nx][ny];
			this.onTrue = onTrue;
			
			//Now, populate the mask for later use
			boolean[][] parentMask = onTrue ? parent.angleMask[0] : parent.angleMask[1];
			latt = onTrue ? parent.onLatt : parent.offLatt;
			
			double[] r, rSym;
			double[] kCent = new double[2];
			kCent[0] = latt.getA()[0]*kguess;
			kCent[1] = latt.getA()[1]*kguess;
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
				{
					r = new double[] {i, j};
					rSym = latt.getAtomicCoords(r);
					mask[i][j] = parentMask[i][j] && Distance.distance(rSym, kCent) < widthConst*widthGuess;
				}
			
			

		}
		
		@Override
		/**
		 * returns the sum of squared residuals between the data and the stored function
		 */
		public double value(double[] param) {
			double sum = 0;
			double[] r;
			double[] rSym;
			
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
					if (mask[i][j])
					{
						//First, set the value of r.
						r = new double[] {i-nx/2, j-ny/2};
						rSym = parent.latt.getAtomicCoords(r);
						sum += Math.pow(f.function(param, rSym) - values[i][j], 2);
					}
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
	 * This is to store guesses supplied by the user 
	 * @author madhavanlab2011
	 *
	 */
	public static class EQPair
	{
		double e;
		double q;
		boolean onTrue; //on or off direction.
		public EQPair(double e, double q, boolean onTrue) {
			super();
			this.e = e;
			this.q = q;
			this.onTrue = onTrue;
		}
		public double getE() {
			return e;
		}
		public double getQ() {
			return q;
		}
		public boolean isOnTrue() {
			return onTrue;
		}
		
		public String toString()
		{
			return "" + q + "\t" + e + "\t" + (onTrue ? "on" : "off");
		}
	}
	
	public static class UserGeneratedFFTFit{
		ArrayList<ArrayList<EQPair>> modes;
		
		int currentMode = -1;
		
		public UserGeneratedFFTFit()
		{
			modes = new ArrayList<ArrayList<EQPair>> ();
			addNewMode();
		}
		
		public void addNewMode()
		{
			modes.add(new ArrayList<EQPair>());
			currentMode++;
		}
		public void removeMode()
		{
			modes.remove(modes.size()-1);
			currentMode--;
		}
		
		public void addEQPair(EQPair eq)
		{
			modes.get(currentMode).add(eq);
		}
		public void removeLastEQPair()
		{
			modes.get(currentMode).remove(modes.get(currentMode).size()-1);
		}

		public void printModes() {
			for (int i = 0; i < modes.size(); i++){System.out.println("Mode " + i);
				for (int j = 0; j < modes.get(i).size(); j++){
					System.out.println(modes.get(i).get(j));
				}
				System.out.println();
			}
		}
	}
	public void processKeyStroke(char c)
	{
		switch (c){
		case 'u':
			if (user == null) user = new UserGeneratedFFTFit();
			else user.addNewMode();
			break;
		case 'U':
			if (user != null) user.removeMode();
			break;
		case 'N':
			if (user != null) user.printModes();
			break;
		case 'F': //do a fit.
			if (user != null);
			
		}
	}
	public static void main(String[] args)
	{
		Topomap.setStdDir();
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		Topomap t = Topomap.open(fc);
		AtomicCoordinatesSet latt = new AtomicCoordinatesSet(ColumnIO.readString(FileOps.selectOpen(fc).toString()));
		FFTCutFitter it = new FFTCutFitter(t.data, latt);
		
//		ImageEditing.copyToClipboard(ImageEditing.weldHorizontal(ImageEditing.getBufferedImage(it.angleMask[0]), ImageEditing.getBufferedImage(it.angleMask[1])));
	}


}
