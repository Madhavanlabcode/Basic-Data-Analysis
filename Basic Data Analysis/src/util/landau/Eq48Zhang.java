package util.landau;

import drawing.GraphDrawerCart;
import util.Printer;
import util.SpectraUtil;
import util.SpectraUtil.LandauLevelPeak;
import util.regression.ACM_CustomFunctions;
import util.regression.ACM_NonLinearFitter;

/**
 * This class attempts to implement Eq. 48 of PRB 82 045122
 * @author madhavanlab2011
 *
 */
public class Eq48Zhang extends ACM_CustomFunctions implements LandauLevelModel {

	//The model parameters:
	//Values are taken from Table IV of the paper for Bi2Te3 
	//except for c0 which is made 0.315 to compare with Wenwen's 1% data (was 0.180 in the paper)
	double c0 = -0.315, c2 = -200, gsz = -50.34, A = 2.87;
	//The value of A-tilde is given as the value of B-tildea
	//c0 is C_0 tilde, same with C2, A is Atilde.
	static final double MUB = 5.7883818066e-5; //This is the Bohr magneton in eV/T.
	static final double HBAROVERE = 2.067833758e5*(1/(Math.PI)); //This is in T*A^2.
	static final double EOH = 1/HBAROVERE;
	
	static final double c0Def = -0.315, c2Def = 49.68, gszDef = -50.34, ADef = 2.87;

	public Eq48Zhang()
	{
		super();
		plist = new String[] {"c0", "c2", "gsz", "Atwiddle"};

	}
	
	public double energy(double N, double B)
	{
		double s = Math.signum(N);
		if (N == 0) {
			return c0 + EOH*B*c2 - (MUB/2)*gsz*B;
		}
		
		return c0 + 2*N*B*EOH*c2 + s*Math.sqrt(Math.pow(-c2*EOH*B + MUB*gsz*B/2, 2) + 2*Math.abs(N)*EOH*B*A*A);
	}
	
	public void setDefaults()
	{
		c0 = c0Def;
		c2 = c2Def;
		gsz = gszDef;
		A = ADef;
	}
	

	
	@Override
	/**
	 * This is intended to leave the parameters c0, c2, gsa, and A in the best-fit values.
	 */
	public void fit(LandauLevelPeak[] peaks) {
		double[][] x = new double [peaks.length][2];
		for (int i = 0; i < peaks.length; i++)
			x[i] = new double[] {peaks[i].n, peaks[i].B};
		double[] y = new double [peaks.length];
		for (int i = 0; i < y.length; i++)
			y[i] = peaks[i].max;
		
		double[] start = new double[] {c0Def, c2Def, gszDef, ADef};
		double[] step = new double[] {c0Def, c2Def, gszDef, ADef};
		
		ACM_NonLinearFitter.FittingResult ans = ACM_NonLinearFitter.fitToFunctionFull(x, y, this, start, step);
		
		double[][] orderedPairs = new double[peaks.length][];
		double[] orderedX = new double [ peaks.length];
		double[] orderedY = new double [peaks.length];
		
		for (int i = 0; i < peaks.length; i++)
		{
			orderedPairs[i] = peaks[i].getOrderedPair();
			orderedX[i] = orderedPairs[i][0];
			orderedY[i] = orderedPairs[i][1];
			System.out.println(peaks[i].n + "\t" + orderedPairs[i][0] + "\t" + orderedPairs[i][1] + "\t" + ans.yCalc[i]);
		}
		
		System.out.println(Printer.arrayLnHorizontal(ans.fitParams));
//		GraphDrawerCart.plotGraphs(orderedX, new double[][] {orderedY, ans.yCalc});
		
		
	}

	@Override
	public void fit(LandauLevelPeak[][] peaks) {
		int npx = 0;
		int[] lengths = new int [peaks.length];
		for (int i = 0; i < peaks.length; i++){
			lengths[i] = peaks[i].length;
			for (int j = 0; j < peaks[i].length; j++)
				npx++;
		}
		
		double[][] x = new double [npx][2];
		double[] y = new double [npx];
		int n = 0;
		for (int i = 0; i < peaks.length; i++)
			for (int j = 0; j < peaks[i].length; j++){
				x[n] = new double[] {peaks[i][j].n, peaks[i][j].B};
				y[n++] = peaks[i][j].max;
		}
		double[] start = new double[] {c0Def, c2Def, gszDef, ADef};
		double factor = 0.01;
		double[] step = new double[] {factor*c0Def, factor*c2Def, factor*gszDef, factor*ADef};
		
		ACM_NonLinearFitter.FittingResult ans = ACM_NonLinearFitter.fitToFunctionFull(x, y, this, start, step);

		double[][] peakX = new double [peaks.length][];
		double[][] peakY = new double [peaks.length][];
		double[][] peakYCalc = new double [peaks.length][];
		n = 0;
		for (int i = 0; i < peaks.length; i++){
			peakX[i] = new double [peaks[i].length];
			peakY[i] = new double [peaks[i].length];
			peakYCalc[i] = new double [peaks[i].length];
				
			for (int j = 0; j < peaks[i].length; j++){
				peakX[i][j] = peaks[i][j].getOrderedPair()[0];
				peakY[i][j] = peaks[i][j].getOrderedPair()[1];
				peakYCalc[i][j] = ans.yCalc[n++];
			}
		}
		
		GraphDrawerCart.plotGraphs(peakX, peakY);
		Printer.printColumnSeries(peakX, peakY);
		GraphDrawerCart.plotGraphs(peakX, peakYCalc);
		Printer.printColumnSeries(peakX, peakYCalc);
		
		
		System.out.println(Printer.arrayVertical(ans.fitParams));
	}
	/**
	 * Here x[0] is N and x[1] is B.
	 */
	public double function(double[] param, double[] x) {
		double[] oldPara = new double[] {c0, c2, gsz, A};
		setParam(param);
		
		double ans = energy(x[0], x[1]);
		setParam(oldPara);
		return ans;
	}
	
	public void setParam(double[] param)
	{
		c0 = param[0];
		c2 = param[1];
		gsz = param[2];
		A = param[3];
	}

	public static void main(String[] args)
	{
		LandauLevelPeak[][] filePeaks = SpectraUtil.getFromACertainFolder(SpectraUtil.LandauLevelFittingFolder.get1PercentFolder(), 6, null);
		
		double B;
		Eq48Zhang zh = new Eq48Zhang();
		int nmax = 8;
		int nB = 20;
		LandauLevelPeak[][] peaksTr = new LandauLevelPeak[nB][nmax];
		LandauLevelPeak[][] peaks = new LandauLevelPeak[nB][nmax];
		LandauLevelPeak[] peakList = new LandauLevelPeak[nB*nmax];
		int n = 0;
		for (int j = 0; j < nB; j++){
			B = (j)/2.0 + 1;
			for (int i = 0; i < nmax; i++)
			{
				int nll = i;
				peakList[n++] = new LandauLevelPeak(zh.energy(nll, B), nll, B);
				peaksTr[j][i] = new LandauLevelPeak(zh.energy(nll, B), nll, B);
			}
		}
		peaks = SpectraUtil.splitByIndex(peakList);
		double[][][] stuff = SpectraUtil.LandauLevelPeak.getOrderedPairs(peaks);
		double[][][] stuff2 = SpectraUtil.LandauLevelPeak.getOrderedPairs(peaksTr);

		Printer.printColumnSeries(stuff[0], stuff[1]);
		
		GraphDrawerCart.plotGraphs(stuff[0], stuff[1]);
		GraphDrawerCart.plotGraphs(stuff2[0], stuff2[1]);
	
		new Eq48Zhang_SliderPanel_2(filePeaks, zh);
	}
	
	public static void spawnSliderPanelFromNB_File()
	{
		
	}

}
