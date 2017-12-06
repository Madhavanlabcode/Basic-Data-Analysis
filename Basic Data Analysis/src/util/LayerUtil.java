package util;

import image.ImageEditing;
import impurity.PointImp;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import org.apache.commons.math3.analysis.interpolation.BicubicSplineInterpolatingFunction;
import org.apache.commons.math3.analysis.interpolation.BicubicSplineInterpolator;

import drawing.GraphDrawerCart;
import drawing.LayerViewer;
import drawing.TopomapViewer;
import main.SRAW;
import schrodinger.MovieMaker;
import util.UnitCellUtil.UnitCell2D;
import util.color.ColorScale;
import util.color.ColorScales;
import util.fileops.ColumnIO;
import util.fileops.FileOps;
import util.fileops.Layer;
import util.fileops.RHKFileOps;
import util.fileops.Topomap;
import util.fourier.AtomicCoordinatesGenerator;
import util.fourier.FFT2DSmall;
import util.fourier.FFTOps;
import util.geom.AtomicCoordinatesSet;
import util.geom.Distance;
import util.matrix.Matrix;
import util.regression.ACM_NonLinearFitter;
import util.regression.TwoDGaussianFreeFitter;

public class LayerUtil {
	static JFileChooser fc;
	
	public static void main(String[] args)
	{
		Topomap.setStdDir();
		fc = new JFileChooser(Topomap.stddir);
		testDeconvolution();
//		Topomap map = Topomap.open(fc);
//		Layer t = map.getLayer(125);
//		Layer t = Layer.openFree(fc);
//		AtomicCoordinatesSet latt = new AtomicCoordinatesSet(FileOps.openText(fc));
//		boolean[][] mask = LayerUtil.BooleanFields.getFilledCutEllipseAroundBraggPeak(t.nx, t.ny, latt, 120, 60, 280, 0);
//		LayerViewer.show(Layer.getFreeLayer(ArrayOps.toDouble(mask)), 1000, true);
//		double[][] data = ColumnIO.readNColumns(FileOps.selectOpen(fc), 3, 2);
//		double[] ans = FieldOps.fitToPlane(data[0], data[1], data[2]);
//		Printer.printlnHorizontal(ans);
	}	
	
	public static void testDeconvolution(){
		Layer h0 = Layer.openFree(fc);
		Layer g0 = Layer.openFree(fc);
		double[][] h = FFTOps.obtainFFTmagCent(h0.data);
		LayerViewer.show(Layer.newLayer(h0, h), 1000, true);
		double[][] g = FFTOps.obtainFFTmagCent(g0.data);
		LayerViewer.show(Layer.newLayer(g0, g), 1000, true);
//		double[][] f = FieldOps.getSimpleFourierDeconvolution(h, g);
		double[][] f = FieldOps.getWienerDeconvolution(h, g, 10000, 10000, 2000000000000.0);
		LayerViewer.show(Layer.newLayer(g0, f), 1000, true);
	}
	public static void contourIntegrateADislocation(){
		Layer t = Layer.openFree(fc);
		Layer u = Layer.openFree(fc);
		AtomicCoordinatesSet latt = new AtomicCoordinatesSet(FileOps.openText(fc));
		double[][] rotMat = Matrix.getRotationMatrix(latt);
		//		double[] rotatedA = Matrix.getProductWith(Matrix.getRotationMatrix(latt), Distance.unitVector(latt.getA()));
//		System.out.println(Printer.arrayLnHorizontal(rotatedA));
		int nlayers = 40;
		double[][][] data = new double [nlayers][][];
		int nx = 240, ny = 240;
		double[][][] vector = new double[2][nx][ny];
//		vector = FieldOps.gradient(TopomapUtil.MapCuttingMethods.getRadialOrigin(nx, ny));
		vector = new double[][][] {t.data, u.data};
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				double[] temp = Matrix.getProductWith(rotMat, new double[] {vector[0][i][j], vector[1][i][j]});
				vector[0][i][j] = temp[0];
				vector[1][i][j] = temp[1];
			}
//		Matrix.rotateVector(vector, angle)
//		Layer t = Layer.getFreeLayer(vector[0]);
//		Layer u = Layer.getFreeLayer(vector[1]);
//		new TopomapViewer_complex2(Topomap.newTopomap(new Layer[] {t}), Topomap.newTopomap(new Layer[] {u}), fc.getCurrentDirectory().toString() + "\\", 512);
		int[] center = Printer.getTwoInts();
		for (int i = 0; i < nlayers; i++)
		{
			ArrayList<int[]> contour = FieldOps.getSquareContour(center, i+1);
			boolean[][] pixels = TopomapUtil.FourierFilterMethods.getSuppressionField(contour, nx, ny);
			data[i] = ArrayOps.toDouble(pixels);
			double c = FieldOps.getLineIntegral(vector, contour);
			System.out.println(i + "\t" + c);
		}
		Topomap m = new Topomap(data, ArrayOps.generateArrayInclBoth(1, nlayers, nlayers), t.x, t.y, null);

	}
	public static void copyImageOfBrillouinZone(Layer t, AtomicCoordinatesSet latt)
	{
		BufferedImage ks = new BufferedImage(t.nx, t.ny, BufferedImage.TYPE_INT_ARGB);
		ImageEditing.color_ARGB(ks, new java.awt.Color(0, 0, 0, 0));
		ImageEditing.drawBrillouinZone(ks, latt, Color.RED);
		ImageEditing.copyToClipboard(ks);
	}
	public static Layer makeSymmetrizedFFT(Layer t, AtomicCoordinatesSet latt, boolean square)
	{
		double[][] newData = new double [t.nx][t.ny];
		boolean log = JOptionPane.showConfirmDialog(null, "Use log scale?") == JOptionPane.YES_OPTION;
//		if (!square) ColumnIO.writeString(LayerUtil.getSymmetrizedTriangleLattice(latt, t.getLayer(0)).toString(), FileOps.selectSave(null).toString());
			if (square)
				newData = LayerUtil.symmetrizeFFT_1(latt, t, log);
			else
				newData = LayerUtil.symmetrizeFFTTriang(latt, t, log);
//			FieldOps.log(newData[i]);
		
		double dx = 2*Math.PI/t.xLength;
		double dy = 2*Math.PI/t.yLength;
		double[] x = new double [t.nx];
		double[] y = new double [t.ny];
		for (int i = 0; i < t.nx; i++)
			x[i] = (i-t.nx/2)*dx;
		for (int j = 0; j < t.ny; j++)
			y[j] = (j-t.ny/2)*dy;
		
		return new Layer(newData, x, y, t.v, t.current);
	}

	public static void rotate180(Layer t)
	{
		ArrayOps.flip(t.x);
		ArrayOps.flip(t.y);
		ArrayOps.flipX(t.data);
		ArrayOps.flipY(t.data);
	}
	public static void printTwoScatterPlotsRt2(PointImp[] imps, AtomicCoordinatesSet latt, double offsetXScatter, double offsetYScatter, double offsetXRt2, double offsetYRt2)
	{
		PointImp[][] split = ImpurityUtil.splitImps(imps, ImpurityUtil.splitImpsWRTRoot2(imps, latt, offsetXRt2, offsetYRt2));
		double[][][] impLattPos = new double[][][] {ImpurityUtil.getLatticePositionsModuloN(split[0], 1, offsetXScatter, offsetYScatter, latt), ImpurityUtil.getLatticePositionsModuloN(split[1], 1, offsetXScatter, offsetYScatter, latt)};
//		ImageEditing.copyToClipboard(ImageEditing.getImpurityModuloMap(impLattPos, imps, 512, 256));
		String[] lines = new String[Math.max(split[0].length, split[1].length)];
		for (int i = 0; i < lines.length; i++)
			lines[i] = "" + (i < split[0].length ? impLattPos[0][i][0] + "\t" + impLattPos[0][i][1] : "\t") + "\t" + (i < split[1].length ? impLattPos[1][i][0] + "\t" + impLattPos[1][i][1] : "\t");
		
		for (int i = 0; i < lines.length; i++)
			System.out.println(lines[i]);
	}
	public static void fitEachImpurityToAGaussian(Layer t, PointImp[] imps)
	{
		//Fit each impurity to a gaussian.
		PointImp[] fitted = new PointImp[imps.length];
		double radius = Double.parseDouble(JOptionPane.showInputDialog("Radius of fitting region? (pixels)"));
		double[][][] temp = TwoDGaussianFreeFitter.fitListOfImpurities(imps, t.data, radius, false);
		double[][] table = temp[0];
		double[][] value = temp[1];
		for (int i = 0; i < fitted.length; i++){
			fitted[i] = new PointImp(new double [] {table[i][0], table[i][1]});
		}
		PointImp.writeToFile(fitted, FileOps.selectSave(fc));
		LayerViewer.show(Layer.newLayer(t, value), 512, false);
		LayerViewer.show(Layer.newLayer(t, temp[2]), 512, false);
		double[][] visibleTopo = new double [t.nx][t.ny];
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
				visibleTopo[i][j] = temp[2][i][j]*t.data[i][j];
		LayerViewer.show(Layer.newLayer(t, visibleTopo), 512, false);
	}
	public static void splitLayerByBinsOfLayer(Layer t, String dir, String name, int nbins)
	{
		Layer dist = Layer.open(fc);
		int[][] bins = FieldOps.getPercentileBinsForField(dist.data, nbins);
		Layer[] split = splitLayerBinned(t, bins);
		for (int i = 0; i < split.length; i++)
			Layer.writeBIN(split[i], dir + name + "split_" + (i+1) + "_of_" + split.length + ".bin");
	}

	public static void makeAtomFittingOutput(Layer t, AtomicCoordinatesSet latt, int npix, int nbins, double fraction, String fname, String outdir, boolean absColorScale)
	{
		UnitCell2D[] cells = UnitCellUtil.getUnitCells(fraction, t, latt, npix, 8);
		double[][] fitParams_0 = UnitCellUtil.fitAtomsToPeakWithPicture(latt, t, fname, outdir, absColorScale, cells);
//		double[][] fitParams = FileOps.openTable(fc);
		boolean[] wanted = chooseWantedAtoms(fitParams_0, npix);
		double[][] fitParams = purgeUnwantedAtoms(wanted, fitParams_0);
		//		double[][] fitParams = FileOps.openTable(fc);
		ColumnIO.writeTable(fitParams, outdir + "params edited.txt");
		ColumnIO.writeTable(fitParams_0, outdir + "params full.txt");
		ColumnIO.writeBin(wanted, outdir + "wanted list.dat");
		
		
		fitLatticeSitesHistogramToGauss(t, latt, nbins, fitParams[3], null, outdir + "height or area.txt");
		fitLatticeSitesHistogramToGauss(t, latt, nbins, fitParams[2], null, outdir + "width.txt");
//		double[] height = new double [fitParams[0].length];
//		for (int i = 0; i < height.length; i++)
//			height[i] = fitParams[3][i]/fitParams[2][i];
//		fitLatticeSitesHistogramToGauss(t, latt, 256, height, null);
		fitLatticeSitesHistogramToGauss(t, latt, nbins, fitParams[4], null, outdir + "offset.txt");
		
		double[] dh_bicubic = new double [cells.length];
		double[] fit_dh = new double [cells.length];
		double[][] fit;
		double[][] params = FieldOps.transpose(fitParams_0);
		for (int i = 0; i < cells.length; i++)
		{
			dh_bicubic[i] = FieldOps.max(cells[i].data) - FieldOps.min(cells[i].data);
			fit = ACM_NonLinearFitter.getExpectedData(cells[i].nx, cells[i].ny, params[i], fname);
			fit_dh[i] = FieldOps.max(fit) - FieldOps.min(fit);
		}
		fitLatticeSitesHistogramToGauss(t, latt, nbins, dh_bicubic, null, outdir + "bicubic dh.txt");
		fitLatticeSitesHistogramToGauss(t, latt, nbins, fit_dh, null, outdir + "fit dh.txt");
		
		
		
		
	}
	/**
	 * This takes a list of fitting parameters (from different unit cells) and purges ones which meet specific requirements.
	 * It is assumed that params[0] and [1] are the (x, y) coordinates of the center and that [2] is the (isotropic) spread.
	 * @param fitParams
	 * @param size
	 * @return
	 */
	public static boolean[] chooseWantedAtoms(double[][] fitParams, int size) {
		boolean[] wanted = new boolean [fitParams[0].length];
		if (size == -1){ 
			for (int i = 0; i < wanted.length; i++)
				wanted[i] = true;
			return wanted;
		}
		
		int nwanted = 0;
		for (int i = 0; i < wanted.length; i++)
		{
			if (fitParams[0][i] > 0 && fitParams[1][i] > 0 && fitParams[0][i] < size && fitParams[1][i] < size) //if the "mean" is within the window
				if (Math.abs(fitParams[2][i]) < 5*size) //if the gaussian spread is not freaking huge
				{
					wanted[i] = true;
					nwanted++;
				}
				else
					System.out.println("Rejected atom " + i);
		}
		return wanted;
	}
	public static double[][] purgeUnwantedAtoms(boolean[] wanted, double[][] fitParams) {
		int nwanted = 0;
		for (int i =0; i < wanted.length; i++)
			if (wanted[i]) nwanted++;
		double[][] paramsLeft = new double [fitParams.length][nwanted];
		int n = 0;
		for (int i = 0; i < wanted.length; i++)
			if (wanted[i])
			{
				for (int j = 0; j < fitParams.length; j++)
					paramsLeft[j][n] = fitParams[j][i];
				n++;
			}
		System.out.println("" + (wanted.length - nwanted) + " atoms were purged.");
		return paramsLeft;
	}

	public static void makeLatticeSitesHistograms(Layer t, AtomicCoordinatesSet latt, int nbins)
	{
//		LatticeSite[] all = getLatticeSites(t, latt);
//		double[] allval = new double [all.length];
//		for (int i = 0; i < all.length; i++)
//		{
//			allval[i] = all[i].value;
////			if (Math.abs(allval[i]) < 1e-13) System.out.println(i + "\t" + all[i]);
//		}
		
//		double[] allval = UnitCellUtil.getLatticeSiteMeans(latt, t, 16, 0.5);
//		double[] allval = UnitCellUtil.getLatticeSiteMeansLiteral(latt, t, 1);
		double[] allval = UnitCellUtil.getLatticeSiteMaxima(latt, t, 0.3);
//		double[] allval = UnitCellUtil.fitAtomsToPeak(latt, t, 32, "Gauss2DConst", 0)[3];
		
		
		
		double[] binmins = ArrayOps.generateArrayNotInclUpper(ArrayOps.min(allval), ArrayOps.max(allval), nbins);
		double mean = FieldOps.mean(t.data);
//		mean = -3.308116467660996e-19;
//		System.out.println(mean);
		int[] hist = ArrayOps.getHistogramExcludeValue(allval, binmins, mean);
		for (int i = 0; i < nbins; i++)
		{
			System.out.println(binmins[i] + "\t" + hist[i]);
		}
	
//		GraphDrawerCart.plotGraph(binmins, hist);
		
//		double[] para = ArrayUtil.fitToFunction(binmins, ArrayOps.toDouble(hist), "TwoGauss");
//		double[] fitCurve = ArrayUtil.getExpectedValues(binmins, para, "TwoGauss");
		double[] para = ACM_NonLinearFitter.fitToFunction(binmins, ArrayOps.toDouble(hist), "TwoGauss");
		double[] fitCurve = ACM_NonLinearFitter.getExpectedY(binmins, para, "TwoGauss");
		
		double[] pg1 = new double[] {para[0], para[1], para[2], 0, 1, 0};
		double[] pg2 = new double[] {0, 1, 0, para[3], para[4], para[5]};
		double[] g1 = ArrayUtil.getExpectedValues(binmins, pg1, "TwoGauss");
		double[] g2 = ArrayUtil.getExpectedValues(binmins, pg2, "TwoGauss");
		double[] freq = new double [g1.length];
		for (int i = 0; i < hist.length; i++)
			freq[i] = hist[i]/(double)(allval.length);
		System.out.println(Printer.arrayVertical(para));
		System.out.println("Total atoms:");
		double g1s = ArrayOps.sum(g1);
		double g2s = ArrayOps.sum(g2);
		System.out.println(g1s);
		System.out.println(g2s);
		System.out.println("Doping percentage:\r\n" + (g1s/(g1s+g2s)));
		double[][] matrix = new double[][] {binmins, ArrayOps.toDouble(hist), fitCurve, g1, g2, freq};
		FileOps.writeTableASCII(fc, matrix);
		
		GraphDrawerCart.plotGraph(matrix, 0);
		
		//		double[] perc = ArrayOps.getPercentiles(allval);
//		for (int i = 0 ; i < 100; i++)
//			System.out.println(i + "\t" + perc[i]);
	}
	/**
	 * This requires no intervention from the user. The pictures written are a graph of the histogram, the layer, and its FFT.
	 * The histogram text files are also written.
	 * @param t
	 * @param latt
	 * @param nbins
	 */
	public static void makeLatticeSitesHistogramsPictures(Layer t, AtomicCoordinatesSet latt, int nbins, String fileSuffix)
	{
//		LatticeSite[] all = getLatticeSites(t, latt);
//		double[] allval = new double [all.length];
//		for (int i = 0; i < all.length; i++)
//		{
//			allval[i] = all[i].value;
////			if (Math.abs(allval[i]) < 1e-13) System.out.println(i + "\t" + all[i]);
//		}
//		double[] allval = UnitCellUtil.getValuesAtInt(latt, t);
//		double[] allval = UnitCellUtil.getLatticeSiteMeans(latt, t, 16, 0.5);
//		double[] allval = UnitCellUtil.getLatticeSiteMeansLiteral(latt, t, 1);
		double[] allval = UnitCellUtil.getLatticeSiteMaxima(latt, t, 0.25);
//		double[] allval = UnitCellUtil.fitAtomsToPeak(latt, t, 32, "Gauss2DConst", 0)[3];
		ArrayList<String> lines = new ArrayList<String>();
		String dir = fc.getCurrentDirectory().toString() + "\\histoutput\\";
		if (!new File(dir).exists()) new File(dir).mkdir();
		
		double[] binmins = ArrayOps.generateArrayNotInclUpper(ArrayOps.min(allval), ArrayOps.max(allval), nbins);
		double mean = FieldOps.mean(t.data);
//		mean = -3.308116467660996e-19;
//		System.out.println(mean);
		int[] hist = ArrayOps.getHistogramExcludeValue(allval, binmins, mean);
		for (int i = 0; i < nbins; i++)
		{
			System.out.println(binmins[i] + "\t" + hist[i]);
		}
	
//		GraphDrawerCart.plotGraph(binmins, hist);
		
//		double[] para = ArrayUtil.fitToFunction(binmins, ArrayOps.toDouble(hist), "TwoGauss");
//		double[] fitCurve = ArrayUtil.getExpectedValues(binmins, para, "TwoGauss");
		double[] para = ACM_NonLinearFitter.fitToFunction(binmins, ArrayOps.toDouble(hist), "TwoGauss");
		double[] fitCurve = ACM_NonLinearFitter.getExpectedY(binmins, para, "TwoGauss");
		
		double[] pg1 = new double[] {para[0], para[1], para[2], 0, 1, 0};
		double[] pg2 = new double[] {0, 1, 0, para[3], para[4], para[5]};
		double[] g1 = ArrayUtil.getExpectedValues(binmins, pg1, "TwoGauss");
		double[] g2 = ArrayUtil.getExpectedValues(binmins, pg2, "TwoGauss");

		String parameters = Printer.arrayVertical(para);
		lines.add(parameters);
		lines.add("Total atoms:");
		double g1s = ArrayOps.sum(g1);
		double g2s = ArrayOps.sum(g2);
		lines.add("" + g1s);
		lines.add("" + g2s);
		lines.add("Doping Percentage:\r\n" + (g1s/(g1s+g2s)));
		for (int i = 0; i < lines.size(); i++)
			System.out.println(lines.get(i));
		double[][] matrix = new double[][] {binmins, fitCurve, g1, g2, ArrayOps.toDouble(hist)};
		ColumnIO.writeTable(matrix, dir + "hist" + fileSuffix + ".txt");
		
		ColumnIO.writeLines(lines, dir + "summary_" + fileSuffix + ".txt");
		
		BufferedImage graph = GraphDrawerCart.getPlot(matrix, 0);
		BufferedImage layer = ImageEditing.getBufferedImage(t.data, ColorScales.NSCALES-1);
		double[][] fftd = FFTOps.obtainFFTmagCent(t.data);
		FieldOps.log(fftd);
		BufferedImage fft = ImageEditing.getBufferedImage(fftd, ColorScales.NSCALES-1);
		
		SRAW.writeImage(dir + "Hist" + fileSuffix, graph);		
		SRAW.writeImage(dir + "layer" + fileSuffix, layer);		
		SRAW.writeImage(dir + "fft" + fileSuffix, fft);		
		
		//Finally, write pictures of the lattice sites:
		File latticeDir = new File(dir + "lattice sites " + fileSuffix + "\\");
		if (!latticeDir.exists()) latticeDir.mkdir();
		UnitCell2D[] cells = UnitCellUtil.getUnitCells(latt, t, 0.25, 2);
		ColorScale scale = ColorScales.getNew(t.data, ColorScales.NSCALES-1);
		for (int i = 0; i < cells.length; i++)
		{
			SRAW.writeImage(latticeDir.toString() + "\\" + MovieMaker.fromInt(i) + "_" + cells[i].getA() + "_" + cells[i].getB(), cells[i].getLiteralImage(scale, t, latt, 1));
		}
				//		double[] perc = ArrayOps.getPercentiles(allval);
//		for (int i = 0 ; i < 100; i++)
//			System.out.println(i + "\t" + perc[i]);
	}
	
	public static Layer cropLayer(Layer t, int xi, int xf, int yi, int yf)
	{
		int nx = (xf-xi), ny = (yf-yi);
		double[][] data = new double [nx][ny];
		double[] x = new double [nx];
		double[] y = new double [ny];
		for (int i = 0; i < nx; i++)
			x[i] = t.x[i+xi];
		for (int j = 0; j < ny; j++)
			y[j] = t.y[j+yi];
		
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				data[i][j] = t.data[i+xi][j+yi];
		
		return new Layer(data,x,y,t.v,t.current);
	}
	public static void fitLatticeSitesHistogramToGauss(Layer t, AtomicCoordinatesSet latt, int nbins, double[] allval, double[] bounds, String outpath)
	{
//		LatticeSite[] all = getLatticeSites(t, latt);
//		double[] allval = new double [all.length];
//		for (int i = 0; i < all.length; i++)
//		{
//			allval[i] = all[i].value;
////			if (Math.abs(allval[i]) < 1e-13) System.out.println(i + "\t" + all[i]);
//		}
		
//		double[] allval = UnitCellUtil.getLatticeSiteMeans(latt, t, 16, 0.5);
//		double[] allval = UnitCellUtil.getLatticeSiteMeansLiteral(latt, t, 1);
//		double[] allval = UnitCellUtil.getLatticeSiteMaxima(latt, t);
		
		
		double[] binmins;
		if (bounds == null)
			binmins = ArrayOps.generateArrayNotInclUpper(ArrayOps.min(allval), ArrayOps.max(allval), nbins);
		else
			binmins = ArrayOps.generateArrayNotInclUpper(bounds[0], bounds[1], nbins);
		double mean = FieldOps.mean(t.data);
//		mean = -3.308116467660996e-19;
//		System.out.println(mean);
		int[] hist = ArrayOps.getHistogramExcludeValue(allval, binmins, mean);
		for (int i = 0; i < nbins; i++)
		{
			System.out.println(binmins[i] + "\t" + hist[i]);
		}
	
//		GraphDrawerCart.plotGraph(binmins, hist);
		
//		double[] para = ArrayUtil.fitToFunction(binmins, ArrayOps.toDouble(hist), "TwoGauss");
//		double[] fitCurve = ArrayUtil.getExpectedValues(binmins, para, "TwoGauss");
		double[] para = ACM_NonLinearFitter.fitToFunction(binmins, ArrayOps.toDouble(hist), "TwoGauss");
		double[] fitCurve = ACM_NonLinearFitter.getExpectedY(binmins, para, "TwoGauss");
		
		double[] pg1 = new double[] {para[0], para[1], para[2], 0, 1, 0};
		double[] pg2 = new double[] {0, 1, 0, para[3], para[4], para[5]};
		double[] g1 = ArrayUtil.getExpectedValues(binmins, pg1, "TwoGauss");
		double[] g2 = ArrayUtil.getExpectedValues(binmins, pg2, "TwoGauss");
		System.out.println(Printer.arrayVertical(para));
		System.out.println("Total atoms:");
		double g1s = ArrayOps.sum(g1);
		double g2s = ArrayOps.sum(g2);
		System.out.println(g1s);
		System.out.println(g2s);
		System.out.println("Doping percentage:\r\n" + (g1s/(g1s+g2s)));
		double[][] matrix = new double[][] {binmins, ArrayOps.toDouble(hist), fitCurve, g1, g2};
		if (outpath == null)
			FileOps.writeTableASCII(fc, matrix);
		else
			ColumnIO.writeTable(matrix, outpath);
		
		GraphDrawerCart.plotGraph(matrix, 0);
		
		//		double[] perc = ArrayOps.getPercentiles(allval);
//		for (int i = 0 ; i < 100; i++)
//			System.out.println(i + "\t" + perc[i]);
	}
	

	/**
	 * This is intended solely for square lattices.
	 * @param latt
	 * @param t
	 * @return
	 */
	public static double[][] symmetrizeFFT_1(AtomicCoordinatesSet latt, Layer t, boolean log)
	{
		double[][] braggVec = latt.getReciprocal();
		double meanTemp = 0;
		if (!log) {
			meanTemp = FieldOps.mean(t.data);
			FieldOps.plusEquals(t.data, -meanTemp);
		}
		double[][] fftmag = FFTOps.obtainFFTmagCent(t.data);
		double mean = FieldOps.mean(fftmag);
		double[] origin = new double[] {t.nx/2, t.ny/2};
		
		double[][][] reflection = new double [2][t.nx][t.ny];
		FieldOps.applyLinearTransformation(fftmag, Matrix.getReflectionMatrix(braggVec[0]), origin, reflection[0], mean);
		FieldOps.applyLinearTransformation(fftmag, Matrix.getReflectionMatrix(braggVec[1]), origin, reflection[1], mean);
		double[][][] ans = new double [3][t.nx][t.ny];
		FieldOps.add(fftmag, reflection[0], ans[0]);
//		ans[1] = FieldOps.rotatePlus90(ans[0]);
		ans[1] = FieldOps.rotatePlus90_aboutPixel(ans[0], t.nx/2, t.ny/2);
		FieldOps.add(ans[0], ans[1], ans[2]);

//		FieldOps.log(fftmag);
//		FieldOps.log(reflection[0]);
//		FieldOps.log(reflection[1]);
//		FieldOps.log(ans[0]);
//		FieldOps.log(ans[2]);
//		String s = FileOps.selectSave(fc).toString();
//		SRAW.writeImage(s + "orig", fftmag);
//		SRAW.writeImage(s + "refl_1", reflection[0]);
//		SRAW.writeImage(s + "refl_2", reflection[1]);
//		SRAW.writeImage(s + "avg_1", ans[0]);
//		SRAW.writeImage(s + "final_answer", ans[2]);
//				FileOps.writeImage(fc, fftmag);
//		FileOps.writeImage(fc, reflection[0]);
//		FileOps.writeImage(fc, reflection[1]);
		if (log) FieldOps.log(ans[2]);
		else FieldOps.plusEquals(t.data, meanTemp);

		return ans[2];
	}
	/**
	 * This method is for triangular lattices where 3 Bragg peaks are equally visible.
	 * It is not to be assumed that they form a perfect hexagon. The hexagon is formed in the method itself.
	 * 
	 * The symmetrization is six-fold in the sense that the six hextants repeat themselves. The hextants themselves
	 * have unique left and right halves.
	 * @param latt
	 * @param t
	 * @return
	 */
	public static double[][] symmetrizeFFTTriang(AtomicCoordinatesSet latt, Layer t, boolean log)
	{
		double[][] braggVec = latt.getReciprocal();
		
		double meanreal = FieldOps.mean(t.data);
		if (!log)
			FieldOps.minusEquals(t.data, meanreal);
		
		double angle = Math.toDegrees(Math.acos((braggVec[0][0]*braggVec[1][0] + braggVec[0][1]*braggVec[1][1])/(Complex.mag(braggVec[0])*Complex.mag(braggVec[1]))));
		double angle_b1 = Math.toDegrees(Complex.phase(braggVec[0]));
		double angle_b2 = angle_b1 + angle;
		
		double avg = (angle_b1 + angle_b2 - (angle > 90 ? 120 : 60))/2;
		System.out.println(angle + "\t" + angle_b1 + "\t" + avg);
		
		avg = Math.toRadians(avg);
		
		double[][] fftmag = FFTOps.obtainFFTmagCent(t.data);
		double mean = FieldOps.mean(fftmag);
		double[] origin = new double[] {t.nx/2, t.ny/2};
		
		double[][] unitBragg = new double [3][2];
		
		unitBragg[0] = new double [] {Math.cos(avg), Math.sin(avg)};
		unitBragg[1] = new double [] {Math.cos(avg+Math.PI/3), Math.sin(avg+Math.PI/3)};
		unitBragg[2] = new double [] {Math.cos(avg+2*Math.PI/3), Math.sin(avg+2*Math.PI/3)};

		
		double[][][] reflection = new double [3][t.nx][t.ny];
		FieldOps.applyLinearTransformation(fftmag, Matrix.getReflectionMatrix(unitBragg[0]), origin, reflection[0], mean);
		FieldOps.applyLinearTransformation(fftmag, Matrix.getReflectionMatrix(unitBragg[1]), origin, reflection[1], mean);
		FieldOps.applyLinearTransformation(fftmag, Matrix.getReflectionMatrix(unitBragg[2]), origin, reflection[2], mean);
		double[][][] ans = new double [3][t.nx][t.ny];
		
		FieldOps.add(reflection[0], reflection[1], ans[0]);
		FieldOps.add(ans[0], reflection[2], ans[0]);
		
		double min = FieldOps.min(fftmag);
		
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
				if (ans[0][i][j] < 3*min)
					ans[0][i][j] = 3*min;

		if (log) FieldOps.log(ans[0]);
		else FieldOps.plusEquals(t.data, meanreal);
		return ans[0];
	}
	public static AtomicCoordinatesSet getSymmetrizedTriangleLattice(AtomicCoordinatesSet latt, Layer t)
	{
		double[][] braggVec = latt.getReciprocal();
		
		double angle = Math.toDegrees(Math.acos((braggVec[0][0]*braggVec[1][0] + braggVec[0][1]*braggVec[1][1])/(Complex.mag(braggVec[0])*Complex.mag(braggVec[1]))));
		double angle_b1 = Math.toDegrees(Complex.phase(braggVec[0]));
		double angle_b2 = angle_b1 + angle;
		
		double avg = (angle_b1 + angle_b2 - (angle > 90 ? 120 : 60))/2;
		System.out.println(angle + "\t" + angle_b1 + "\t" + avg);
		
		avg = Math.toRadians(avg);
		
		double[][] fftmag = FFTOps.obtainFFTmagCent(t.data);
		double mean = FieldOps.mean(fftmag);
		double[] origin = new double[] {t.nx/2, t.ny/2};
		
		double[][] unitBragg = new double [3][2];
		
		unitBragg[0] = new double [] {Math.cos(avg), Math.sin(avg)};
		unitBragg[1] = new double [] {Math.cos(avg+Math.PI/3), Math.sin(avg+Math.PI/3)};
		unitBragg[2] = new double [] {Math.cos(avg+2*Math.PI/3), Math.sin(avg+2*Math.PI/3)};
		
		return AtomicCoordinatesSet.generateCentered(unitBragg, t.nx);
	}
	public static Layer[] getLatticeDerivatives(Layer t, AtomicCoordinatesSet latt)
	{
		double[] aHat = Distance.unitVector(latt.getA());
		double[] bHat = Distance.unitVector(latt.getB());
		
		Layer dfda = Layer.newLayer(t, FieldOps.directionalDerivative2D(t.data, aHat));

		Layer dfdb = Layer.newLayer(t, FieldOps.directionalDerivative2D(t.data, bHat));
		return new Layer[] {dfda, dfdb};
	}
	public static Layer getDerivativeNematicity(Layer t, AtomicCoordinatesSet latt)
	{
		double[] aHat = Distance.unitVector(latt.getA());
		double[] bHat = Distance.unitVector(latt.getB());
		
		Layer dfda = Layer.newLayer(t, FieldOps.getSecondDirectionalDerivative2D(t.data, aHat));

		Layer dfdb = Layer.newLayer(t, FieldOps.getSecondDirectionalDerivative2D(t.data, bHat));
		Layer ans = Layer.newLayer(t, FieldOps.minus(dfda.data, dfdb.data));
		return ans;
	}
	public static Layer fourierFilter(Layer t, ArrayList<int[]> includedPts)
	{
		double[][] newData = new double [t.nx][t.ny];
		FFT2DSmall fft;
		fft = FFTOps.obtainFFT(t.data);
		
		int[] x = new int [includedPts.size()], y = new int [includedPts.size()];
		for (int i = 0; i < includedPts.size(); i++)
		{
			x[i] = includedPts.get(i)[0];
			y[i] = includedPts.get(i)[1];
		}
		double minMag = FieldOps.magMin(fft.fHat);
//		double minMag = Math.exp(-20);
		
		fft = FFTOps.obtainFFT(FieldOps.copy(t.data));
		FFTOps.supressElseModesTo(fft, x, y, minMag);
		fft.doIFFT();
		for (int j = 0; j < t.nx; j++)
			for (int k = 0; k < t.ny; k++)
				newData[j][k] = fft.f[j][k][0];
			
		return Layer.newLayer(t, newData);
	}
	public static void fourierFilter(Layer t, boolean[][] suppress)
	{
		double[][] newData = new double [t.nx][t.ny];
		FFT2DSmall fft;
		fft = FFTOps.obtainFFT(t.data);
		
		double minMag = FieldOps.magMin(fft.fHat);
//		Systme.out.println(FieldOps.magMin(fft.fHat));
//		double minMag;
		if (minMag == 0) minMag = Math.exp(-20);
		System.out.println(Math.log(minMag));
		
		fft = FFTOps.obtainFFT(FieldOps.copy(t.data));
		FFTOps.supressModesTo(fft, suppress, minMag);
		fft.doIFFT();
		for (int j = 0; j < t.nx; j++)
			for (int k = 0; k < t.ny; k++)
				newData[j][k] = fft.f[j][k][0];

		t.data = newData;
	}
	
	/**
	 * Returns a layer with half the number of pixels in each direction as the input layer. (Useful for 
	 * reducing memory/time requirements of U-field.)
	 * @return the new layer
	 */
	public static Layer coarsen2(Layer t)
	{
		double[][] data = FieldOps.reduce(2, t.data);
		double[] x = new double [t.nx/2], y = new double [t.ny/2]; 
		for (int i = 0; i < x.length; i++)
			x[i] = t.x[2*i];
		for (int i = 0; i < y.length; i++)
			y[i] = t.y[2*i];
		
		return new Layer (data, x, y, t.v, t.current);
		
	}
	
	public static Layer coarsenN(Layer t, int N)
	{
		double[][] data = FieldOps.reduce(N, t.data);
		double[] x = new double [t.nx/N], y = new double [t.ny/N]; 
		for (int i = 0; i < x.length; i++)
			x[i] = t.x[N*i];
		for (int i = 0; i < y.length; i++)
			y[i] = t.y[N*i];
		
		return new Layer (data, x, y, t.v, t.current);
		
	}
	
	/**
	 * This is supposed to obtain a list of "Lattice Sites" including the value of the topography at those points.
	 * @author madhavanlab2011
	 *
	 */
	public static LatticeSite[] getLatticeSites(Layer t, AtomicCoordinatesSet latt)
	{
		double[][][] latticeSites = AtomicCoordinatesGenerator.getLatticeSites(latt, t.nx, null);
		int la = latticeSites.length;
		int lb = latticeSites[0].length;
		boolean[][] inArea = new boolean [la][lb];
		int n = 0;
		int iorigin = 0, jorigin = 0;
		double[] tempat;
		for (int i = 0; i < la; i++)
			for (int j = 0; j < lb; j++){
				inArea[i][j] = (latticeSites[i][j][0] <= t.nx && latticeSites[i][j][0] > 0 && latticeSites[i][j][1] <= t.ny && latticeSites[i][j][1] > 0);
				tempat = latt.getAtomicCoords(latticeSites[i][j]);
				if (FieldOps.round(tempat[0]) == 0 && FieldOps.round(tempat[1]) == 0)
				{
					iorigin = i;
					jorigin = j;
					System.out.println("Origin at (" + iorigin + ", " + jorigin + ")");
				}
				if (inArea[i][j]) n++;
			}
		System.out.println(n);
		LatticeSite[] ls = new LatticeSite[n];
		double[] offset = new double [] {0, 0};
		int m = 0;
		for (int i = 0; i < la; i++)
			for (int j = 0; j < lb; j++)
			{	
				if (inArea[i][j]){
					ls[m++] = new LatticeSite(offset, new int[] {i-iorigin, j-jorigin}, t, latticeSites[i][j]);
				}
			}
		
		return ls;

	}
	/**
	 * This method is garaunteed to produce the two-dimensional array of adjacent pts. 
	 * @param t
	 * @param latt
	 * @return
	 */
	public static Topomap convertToTopomap(Layer l)
	{
		return Topomap.newTopomap(new Layer[] {l});
	}
	
	/**
	 * This computes the "rt2 strength" in each rt2 unit cell to be found in the map. The output is a "tiled" 
	 * picture in which the tiles are rt2 unit cells. The "rt2 strength" is defined as the layer value on the rt2 site
	 * minus the average value of the 4 adjacent 1x1 atoms, divided by the summed magnitudes.
	 * @param t
	 * @param oneXone
	 * @return
	 */
	public static Layer makeRt2Map(Layer t, AtomicCoordinatesSet oneXone)
	{
		double[][] rt2Field = new double[t.nx][t.ny];
		
		AtomicCoordinatesSet rt2 = oneXone.getRt2Lattice();
		TwoDLattice rt2latt = new TwoDLattice(t, rt2);
		
		double[] rt2coord;
		int rt2x, rt2y;
		double[] rt2PointPix;
		double[] atcoord;
		double[] adjacent = new double [4];
		double adjacentAvg;
		int atx, aty;
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
			{
				rt2coord = rt2.getAtomicCoords(i, j);
				rt2x = FieldOps.round(rt2coord[0]);
				rt2y = FieldOps.round(rt2coord[1]);
				rt2PointPix = rt2.getPixelCoords(rt2x, rt2y);
				atcoord = oneXone.getAtomicCoords(rt2PointPix);
				atx = FieldOps.round(atcoord[0]);
				aty = FieldOps.round(atcoord[1]);
				adjacent[0] = t.evaluateAt(oneXone.getPixelCoords(atx-1, aty));
				adjacent[1] = t.evaluateAt(oneXone.getPixelCoords(atx+1, aty));
				adjacent[2] = t.evaluateAt(oneXone.getPixelCoords(atx, aty-1));
				adjacent[3] = t.evaluateAt(oneXone.getPixelCoords(atx, aty+1));
				adjacentAvg = ArrayOps.sum(adjacent)/4;
				rt2Field[i][j] = (t.evaluateAt(rt2PointPix) - adjacentAvg);//(Math.abs(t.evaluateAt(rt2PointPix)) + Math.abs(adjacentAvg));
			}
		return Layer.newLayer(t, rt2Field);
	}
	
	public static void markAllAtomsWithinPercentileRange(Layer t, Layer forPic, int minp, int maxp, AtomicCoordinatesSet latt)
	{
		LatticeSite[] all = getLatticeSites(t, latt);
		double[] allval = new double [all.length];
		for (int i = 0; i < all.length; i++)
			allval[i] = all[i].value;
		
		double[] perc = ArrayOps.getPercentiles(allval);
		
		BufferedImage im = ImageEditing.getBufferedImage(forPic.data, 0);
		Graphics g = im.getGraphics();
		int x, y;
		for (int i = 0; i < all.length; i++)
		{
			if (all[i].value >= perc[minp] && all[i].value < perc[maxp])
			{
				x = FieldOps.round(all[i].rpix[0]);
				y = FieldOps.round(all[i].rpix[1]);
				g.setColor(java.awt.Color.BLUE);
				LayerViewer.drawPlus(g, x, y, 2);
			}
		}
		SRAW.writeImage(FileOps.selectSave(null).toString(), im);
	}
	public static class LatticeSite
	{
		public double[] offset; //the offset in lattice units. Default (0, 0).
		public double value; //the value of the field at this point.
		public int[] rlatt; //the position of the point in lattice units.
		public double[] rpix; //the pixel position corresponding to the point.
		
		public LatticeSite(double[] offset, int[] rlatt, Layer l, double[] rpix) {
			super();
			this.offset = offset;
			this.rlatt = rlatt;
			this.rpix = rpix;
			if (l == null)
				value = 0;
			else
				value = l.evaluateAt(rpix);
		}
		
		public String toString()
		{
			return "Lattice: " + Printer.vectorP(rlatt) + "\tPixels: " + Printer.vectorP(rpix) + "\t Value = " + value;
		}
	}
	
	/**
	 * Implementing a 2d array of LatticeSites with information of whether they are or are not in a set of specified bounds, as well as a shift of lattice origin.
	 * @author madhavanlab2011
	 *
	 */
	public static class TwoDLattice
	{
		public LatticeSite[][] sites;
		public int na, nb;
		public boolean[][] inArea;
		public int iorigin = 0, jorigin = 0;
		
		public TwoDLattice(Layer t, AtomicCoordinatesSet latt)
		{
			double[][][] latticeSites = AtomicCoordinatesGenerator.getLatticeSites(latt, t.nx, null);
			na = latticeSites.length;
			nb = latticeSites[0].length;
			inArea = new boolean [na][nb];
			sites = new LatticeSite[na][nb];
			double[] tempat;
			for (int i = 0; i < na; i++)
				for (int j = 0; j < nb; j++){
					inArea[i][j] = (latticeSites[i][j][0] <= t.nx && latticeSites[i][j][0] > 0 && latticeSites[i][j][1] <= t.ny && latticeSites[i][j][1] > 0);
					tempat = latt.getAtomicCoords(latticeSites[i][j]);
					if (FieldOps.round(tempat[0]) == 0 && FieldOps.round(tempat[1]) == 0)
					{
						iorigin = i;
						jorigin = j;
						System.out.println("Origin at (" + iorigin + ", " + jorigin + ")");
					}
				}
			double[] offset = new double [] {0, 0};
			for (int i = 0; i < na; i++)
				for (int j = 0; j < nb; j++){
					sites[i][j] = new LatticeSite(offset, new int[] {i-iorigin, j-jorigin}, t, latticeSites[i][j]);
				}
		}
		
		public LatticeSite get(int a, int b)
		{
			return sites[a+iorigin][b+jorigin];
		}
	}
	
	public static void writeGradPicture(Layer t)
	{
		double[][][] grad = FieldOps.gradientNM2(t.data);
		BufferedImage out = ImageEditing.getBufferedImage(grad, null);
		String s = FileOps.selectSave(fc).toString();
		SRAW.writeImage(s, out);
		double[][] mag = FieldOps.magnitude(grad);
		SRAW.writeImage(s + "mag", mag);
	}
	
	/**
	 * Writes a boolean map according to whether or not each pixel belongs to a "step edge".
	 * The simple criterion is whether the magnitude of the gradient at that pixel exceeds a cutoff.
	 * @param t
	 */
	public static void writeSharpnessMap(Layer t, double cutoff)
	{
		double[][][] grad = FieldOps.gradientNM2(t.data);
		double[][] mag = FieldOps.magnitude(grad);
		boolean[][] isStep = FieldOps.isGreaterThan(mag, cutoff);
		SRAW.writeImage(FileOps.selectSave(fc).toString(), isStep);
	}
	
	/**
	 * This writes a boolean map according to which a "sharp" region is defined as a step edge provided
	 * its area exceeds a certain number of pixels.
	 * @param t
	 * @param cutoff
	 */
	public static ArrayList<ArrayList<int[]>> getStepEdges(Layer t, double cutoff, int sizeCutoff)
	{
		double[][][] grad = FieldOps.gradientNM2(t.data);
		double[][] mag = FieldOps.magnitude(grad);
		boolean[][] isStep = FieldOps.isGreaterThan(mag, cutoff);
		
		
		
		//Go through each step edge and decide if it is large enough:
		ArrayList<ArrayList<int[]>> edges = new ArrayList<ArrayList<int[]>>();
		
		boolean[][] visited = new boolean[t.nx][t.ny];
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
				visited[i][j] = !isStep[i][j];
		
		ArrayList<int[]> queue = new ArrayList<int[]>();
		ArrayList<int[]> thisStep = new ArrayList<int[]>();
		int count = 0;
		int[] currentPoint;
		int x, y;
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
				if (!visited[i][j])
				{
					System.out.println(count);
					count++;
					thisStep = new ArrayList<int[]>();
					queue.add(new int[] {i,j});
					thisStep.add(new int[] {i,j});
					while(queue.size() > 0)
					{
						currentPoint = queue.remove(0);
						x = currentPoint[0]; y = currentPoint[1];
						visited[x][y] = true;
						if (x > 0 && !visited[x-1][y]){
							queue.add(new int[] {x-1,y});
							thisStep.add(new int[] {x-1,y});
							visited[x-1][y] = true;
						}
						if (y > 0 && !visited[x][y-1]){
							queue.add(new int[] {x,y-1});
							thisStep.add(new int[] {x,y-1});
							visited[x][y-1] = true;
						}
						if (x < t.nx-1 && !visited[x+1][y]){
							queue.add(new int[] {x+1,y});
							thisStep.add(new int[] {x+1,y});
							visited[x+1][y] = true;
						}
						if (y < t.ny-1 && !visited[x][y+1]){
							queue.add(new int[] {x,y+1});
							thisStep.add(new int[] {x,y+1});
							visited[x][y+1] = true;
						}
					}
					
					//Now, we only add step edges above the cutoff:
					if (thisStep.size() >= sizeCutoff)
						edges.add(thisStep);
				}

		//Now we go through and switch isStep to the actual isStep:
		return edges;
	}
	
	public static Layer fitEachPlateauToPlaneAndZeroTheEdges(Layer t, boolean[][] isStepEdge)
	{
		double[][] data = t.data.clone();
		int[][] blobIDs = FieldOps.splitByStepEdges(isStepEdge);
		FieldOps.subtractPlaneFits(data, blobIDs);
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++)
				if (blobIDs[i][j] == -1)
					data[i][j] = 0;
		return Layer.newLayer(t, data);
	}
	public static Layer fitEachPlateauToParabolaAndZeroTheEdges(Layer t, boolean[][] isStepEdge)
	{
		double[][] data = t.data.clone();
		int[][] blobIDs = FieldOps.splitByStepEdges(isStepEdge);
		FieldOps.subtractPlaneFitsNoMean(data, blobIDs);
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++)
				if (blobIDs[i][j] == -1)
					data[i][j] = 0;
		
		FieldOps.changeZeroToAverage(data);
		return Layer.newLayer(t, data);
	}
	public static boolean[][] getStepEdgeMap(Layer t, double cutoff, int sizeCutoff)
	{
		ArrayList<ArrayList<int[]>> edges = getStepEdges(t, cutoff, sizeCutoff);
		boolean[][] isStep = new boolean [t.nx][t.ny];
		for (int i = 0; i < edges.size(); i++)
			for (int j = 0; j < edges.get(i).size(); j++)
			{
				isStep[edges.get(i).get(j)[0]][edges.get(i).get(j)[1]] = true;
			}
		return isStep;
	}
	public static void writeStepEdgeMap(Layer t, double cutoff, int sizeCutoff)
	{
		boolean[][] isStep = getStepEdgeMap(t, cutoff, sizeCutoff);
		File f = FileOps.selectSave(fc);
		SRAW.writeImage(f.toString(), isStep);
//		ColumnIO.writeBin(isStep, f.toString() + ".dat");
	}
	
	/**
	 * This makes a map of the step edge heights. Each pixel has a "step edge height"
	 * which is defined as zero (if it is not a step edge) or a value equal to the height
	 * at some point forward along the gradient vector, minus the height at some point
	 * backward along the gradient vector.
	 * 
	 * @param t
	 * @param cutoff
	 * @param sizeCutoff
	 * @param vectorLength
	 */
	public static double[][] getStepEdgeHeightMap(Layer t, double cutoff, int sizeCutoff, double vectorLength)
	{
		double[][][] grad = FieldOps.gradient(t.data);
		double[] locGrad;
		boolean[][] isStep = getStepEdgeMap(t, cutoff, sizeCutoff);
		double[][] stepHeight = new double [t.nx][t.ny];
		double[] unitVector = new double [2];
		double[] forePoint = new double [2], aftPoint = new double [2];
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
				if (isStep[i][j])
				{
					locGrad = new double[] {grad[0][i][j], grad[1][i][j]};
					unitVector = Distance.unitVector(locGrad);
					forePoint[0] = i + vectorLength*unitVector[0]/2;
					forePoint[1] = j + vectorLength*unitVector[1]/2;
					aftPoint[0] = i - vectorLength*unitVector[0]/2;
					aftPoint[1] = j - vectorLength*unitVector[1]/2;
					stepHeight[i][j] = t.evaluateAt(forePoint) - t.evaluateAt(aftPoint);
					stepHeight[i][j] = Math.max(stepHeight[i][j], 0);
				}
		
		return stepHeight;
	}
	
	public static void writeStepEdgeHeightStuff(Layer t, double gradcut, int sizecut)
	{
		int lineLength = 15;
		RHKFileOps.doFitting(t, 10);
		Layer stepHeight = Layer.newLayer(t, getStepEdgeHeightMap(t, gradcut, sizecut, lineLength));
		double[] heights = FieldOps.getArray(stepHeight.data);
		double[] bins = ArrayOps.generateArrayInclBoth(FieldOps.min(stepHeight.data) + 1e-11, FieldOps.max(stepHeight.data) + 1e-11, 256);
		int[] hist = ArrayOps.getHistogram(heights, bins);
		
		Layer.writeBIN(stepHeight, fc);
		ColumnIO.writeLines(Printer.getHistLines(bins, hist), fc.getSelectedFile().toString() + "_hist.txt");
		String note = "Gradient cutoff  = " + gradcut + "\r\n"+
				"Size cutoff = " + sizecut + " pixels\r\n"+
				"Gradient Line length " + lineLength;
		ColumnIO.writeString(note, fc.getSelectedFile().toString() + "_data note.txt");
		GraphDrawerCart.plotGraph(bins, hist);
		for (int i = 0; i < hist.length; i++)
			System.out.println(bins[i] + "\t" + hist[i]);
		new LayerViewer(stepHeight, fc.getCurrentDirectory().toString(), 512);
	}
	
	public static Layer getStepHeightMapDefined(Layer t, ArrayList<StepEdgePixel> defined)
	{
		double[][] data = new double [t.nx][t.ny];
		for (int i = 0; i < defined.size(); i++)
		{
			data[defined.get(i).pixel[0]][defined.get(i).pixel[1]] = defined.get(i).getHeight(t);
		}
		return Layer.newLayer(t, data);
	}
	public static void writeStepEdgeHeightStuffDefined(Layer t, double gradcut, int sizecut, double minLength)
	{
		RHKFileOps.doFitting(t, 10);
		boolean[][] isStep = LayerUtil.getStepEdgeMap(t, gradcut, sizecut);
		ArrayList<LayerUtil.StepEdgePixel> defined = LayerUtil.getStepEdgePixelsWithUpperLower(isStep, t, minLength);
		Layer stepHeight = getStepHeightMapDefined(t, defined);
		double[] heights = FieldOps.getArray(stepHeight.data);
		double[] bins = ArrayOps.generateArrayInclBoth(1e-11, FieldOps.max(stepHeight.data) + 1e-11, 256);
		int[] hist = ArrayOps.getHistogram(heights, bins);
		
		Layer.writeBIN(stepHeight, fc);
		ColumnIO.writeLines(Printer.getHistLines(bins, hist), fc.getSelectedFile().toString() + "_hist.txt");
		String note = "Gradient cutoff  = " + gradcut + "\r\n"+
				"Size cutoff = " + sizecut + " pixels\r\n"+
				"Gradient Line length " + minLength;
		SRAW.writeImage(fc.getSelectedFile().toString() + " Step edge pixels", LayerUtil.StepEdgePixel.getImage(defined, t));
		ColumnIO.writeString(note, fc.getSelectedFile().toString() + "_data note.txt");
		GraphDrawerCart.plotGraph(bins, hist);
		for (int i = 0; i < hist.length; i++)
			System.out.println(bins[i] + "\t" + hist[i]);
		new LayerViewer(stepHeight, fc.getCurrentDirectory().toString(), 512);
	}
	
	/**
	 * Part of a series on k-space drift correction. The FFT is NOT considered to be centered.
	 * @param t
	 * @param bragg
	 * @param mask
	 * @param expU
	 */
	public static void putExpUFourier(double[][][] shiftedFFTZ, double[][] mask, double[][][] expU, double[][][] expUfft)
	{
		FieldOps.multiply(shiftedFFTZ, mask, expUfft);
		FFTOps.putIFFT(expUfft, expU, false);
	}
	
	/**
	 * This method must do everything including output the results.
	 * @param t
	 * @param mask - the array [nlayers][nx][ny]
	 */
	public static void doDriftCorrectFTSeries(Layer t, double[][][] mask, int[][] braggi, boolean includeLinear, String outDir, boolean writePictures)
	{
		File f = new File(outDir);
		if (!f.exists()) f.mkdirs();
		String imDir = outDir + "Image output\\";
		File fim = new File(imDir);
		if (!fim.exists() && writePictures) fim.mkdirs();
		
		int nlayers = mask.length;
		double[][][] fftz = new double [t.nx][t.ny][2];
		double[][][] expU = new double [t.nx][t.ny][2];
		double[][][] expUfft = new double [t.nx][t.ny][2];
		double[][][][] shiftedFFTZ = new double [2][t.nx][t.ny][2];
		
		double[][][][] expUPhase = new double [2][nlayers][t.nx][t.ny];
		BufferedImage expUIM = new BufferedImage(t.nx, t.ny, BufferedImage.TYPE_INT_RGB);
		BufferedImage expUFFTZIM = new BufferedImage(t.nx, t.ny, BufferedImage.TYPE_INT_RGB);
		BufferedImage phaseImage;
		FFTOps.putFFT(t.data, fftz, false);
		for (int i = 0; i < 2; i++)
			FieldOps.shift(fftz, shiftedFFTZ[i], -braggi[i][0], -braggi[i][1]);
		double temp;
		
		double naturalMin = FieldOps.magMin(fftz);
		//First loop: get expU for each bragg and mask, and write the output before application of the ufield
		Topomap phase;
		double[] indices = ArrayOps.generateArrayInclBoth(1, mask.length, mask.length);
		for (int n = 0; n < 2; n++){
			for (int i = 0; i < mask.length; i++)
			{
				putExpUFourier(shiftedFFTZ[n], mask[i], expU, expUfft);
				FieldOps.phase(expU, expUPhase[n][i]);
				
				//write output:
				if (writePictures){
					expUIM = ImageEditing.getBufferedImage(expU, null);
					SRAW.writeImage(imDir + "expu_b" + n + "_" + MovieMaker.fromInt(i), expUIM);
					//write the fft image
					expUFFTZIM = FFTOps.getImageCent(expUfft, true, false, naturalMin);
					SRAW.writeImage(imDir + "expufft_b" + n + "_" + MovieMaker.fromInt(i), expUFFTZIM);
					phaseImage = ImageEditing.getBufferedImage(expUPhase[n][i]);
					SRAW.writeImage(imDir + "phase_b" + n + "_" + MovieMaker.fromInt(i), phaseImage);
				}
			}
			phase = new Topomap(expUPhase[n], indices, t.x, t.y, null);
			Topomap.writeBIN(phase, outDir + "phase_b" + n + ".bin");
		}
		
		//The above information is enough to apply any combination of u-field to the topography.
		//Now we will attempt to apply the u-field at each mask.
		int[][][] phaseN = new int [2][t.nx][t.ny];
		double[][][] phaseCont = new double [2][t.nx][t.ny];
		double[][] braggTrue = new double [2][2];
		for (int i = 0; i < 2; i++)
		{
			braggTrue[i][0] = braggi[i][0]*2*Math.PI/t.nx;
			braggTrue[i][1] = braggi[i][1]*2*Math.PI/t.nx;
		}
		double[][][] u = new double [t.nx][t.ny][2];
		BufferedImage uIM;
		double[][][] after = new double[mask.length][t.nx][t.ny];
		BicubicSplineInterpolatingFunction interp = null;
		BicubicSplineInterpolator erp = new BicubicSplineInterpolator();
		interp = erp.interpolate(ArrayOps.generateArrayInclBoth(0, t.nx-1, t.nx), ArrayOps.generateArrayInclBoth(0, t.ny-1, t.ny), t.data);
		double mean = FieldOps.mean(t.data);
		for (int k = 0; k < mask.length; k++)
		{
			for (int i = 0; i < 2; i++)
			{
				FieldOps.putPhaseSteps(expUPhase[i][k], t.nx/2, t.ny/2, phaseN[i]);
				FieldOps.putAddedPhaseSteps(expUPhase[i][k], phaseN[i], 2*Math.PI, phaseCont[i]);
				if (writePictures){
					phaseImage = ImageEditing.getBufferedImage(phaseCont[i]);
					SRAW.writeImage(imDir + "cont_phase_b" + i + "_" + MovieMaker.fromInt(k), phaseImage);
				}
			}
			FieldOps.putU(phaseCont[0], phaseCont[1], braggTrue, -1, u);
			FieldOps.subtractAvg(u);
			if (writePictures){
				uIM = ImageEditing.getBufferedImage(u, null);
				SRAW.writeImage(imDir + "u_" + MovieMaker.fromInt(k), uIM);
			}
			FieldOps.applyUFieldBiCubic(interp, u, after[k], mean);
		}
		
		Topomap done = new Topomap(after, indices, t.x, t.y, null);
		Topomap.writeBIN(done, outDir + "after.bin");
		new TopomapViewer(done, outDir, 1024);
//		System.exit(0);
	}
	public static ArrayList<StepEdgePixel> getStepEdgePixels(boolean[][] isStep, Layer t)
	{
		double[][][] grad = FieldOps.gradient(t.data);
		ArrayList<StepEdgePixel> p = new ArrayList<StepEdgePixel>();
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
				if (isStep[i][j])
					p.add(new StepEdgePixel(new int[] {i, j}, grad));
		return p;
	}
	public static Layer[] splitLayerBinned(Layer t, int[][] bins)
	{
		int min = FieldOps.min(bins);
		int max = FieldOps.max(bins);
		int n = (max-min)+1;
		double[][][] data = new double [n][t.nx][t.ny];
		data = FieldOps.splitByBins(t.data, bins);
		for (int j = 0; j < n; j++)
		{
			FieldOps.changeZeroToAverage(data[j]);
		}

		Layer[] ans = new Layer[n];
		for (int i = 0; i < n; i++)
			ans[i] = Layer.newLayer(t, data[i]);
		return ans;
	}

	public static ArrayList<StepEdgePixel> getStepEdgePixelsWithUpperLower(boolean[][] isStep, Layer t, double minLength){
		double[][][] grad = FieldOps.gradient(t.data);
		ArrayList<StepEdgePixel> p = new ArrayList<StepEdgePixel>();
		StepEdgePixel temp;
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
				if (isStep[i][j])
				{
					temp = new StepEdgePixel (new int[] {i, j}, grad);
					temp.setUpperLower(isStep, minLength);
					if (temp.pixDown != null)
						p.add(temp);
				}
			
		return p;
	}
	public static class StepEdgePixel{
		public int[] pixel;
		public double[] grad;
		
		//The nearest pixel above and below the step edge, in the direction of the gradient.
		public int[] pixUp = new int [2], pixDown = new int [2];
		public StepEdgePixel(int[] pixel, double[][][] grad)
		{
			this.pixel = pixel;
			this.grad = new double[] {grad[0][pixel[0]][pixel[1]], grad[1][pixel[0]][pixel[1]]};
		}
		
		public double getHeight(Layer t) {
			return (t.data[pixUp[0]][pixUp[1]] - t.data[pixDown[0]][pixDown[1]]);
		}

		public void setUpperLower(boolean[][] isStep, double minLength)
		{
			double[] unitVector = new double [2];
			unitVector = Distance.unitVector(grad);
			boolean done = false;
			double length = minLength;
			while(!done)
			{
				pixUp[0] = FieldOps.round(pixel[0] + length*unitVector[0]/2);
				pixUp[1] = FieldOps.round(pixel[1] + length*unitVector[1]/2);
				pixDown[0] = FieldOps.round(pixel[0] - length*unitVector[0]/2);
				pixDown[1] = FieldOps.round(pixel[1] - length*unitVector[1]/2);

				if (pixUp[0] < 0 || pixUp[0] >= isStep.length ||
						pixUp[1] < 0 || pixUp[1] >= isStep[0].length ||
						pixDown[0] < 0 || pixDown[0] >= isStep.length ||
						pixDown[1] < 0 || pixDown[1] >= isStep[0].length)
				{
					pixUp = null;
					pixDown = null;
					return;
				}
				done = !isStep[pixUp[0]][pixUp[1]] && !isStep[pixDown[0]][pixDown[1]]; 
				length++;
			}
		}
		public void setUpperLowerFixed(double length)
		{
			double[] unitVector = new double [2];
			unitVector = Distance.unitVector(grad);
			pixUp[0] = FieldOps.round(pixel[0] + length*unitVector[0]/2);
			pixUp[1] = FieldOps.round(pixel[1] + length*unitVector[1]/2);
			pixDown[0] = FieldOps.round(pixel[0] - length*unitVector[0]/2);
			pixDown[1] = FieldOps.round(pixel[1] - length*unitVector[1]/2);
		}
		
		public static BufferedImage getImage(ArrayList<StepEdgePixel> defined, Layer t)
		{
			BufferedImage ans = new BufferedImage(t.nx, t.ny, BufferedImage.TYPE_INT_RGB);
			ImageEditing.blacken(ans);
			for (int i = 0; i < defined.size(); i++)
			{
				ans.setRGB(defined.get(i).pixel[0], defined.get(i).pixel[1], Color.white.getRGB());
				ans.setRGB(defined.get(i).pixUp[0], defined.get(i).pixUp[1], Color.red.getRGB());
				ans.setRGB(defined.get(i).pixDown[0], defined.get(i).pixDown[1], Color.blue.getRGB());
			}
			return ans;
		}
	}
	
	public static class BooleanFields{
		/**
		 * This produces a cut, solid ellipse whose center is at the Bragg peak.
		 * 
		 * @param nx
		 * @param ny
		 * @param latt
		 * @param a - the semi-axis directed along the bragg vector
		 * @param b - the transverse semi-axis
		 * @param cutoff - if a pixel is further along the bragg vector from the origin than this, false.
		 * 
		 * @return
		 */
		public static boolean[][] getFilledCutEllipseAroundBraggPeak(int nx, int ny, AtomicCoordinatesSet latt, double a, double b, double cutoff, int braggIndex){
			int[][] bragg = AtomicCoordinatesSet.generateBragg(latt, nx);
			AtomicCoordinatesSet recip = latt.getReciprocalLattice();
			
			double[] aHat, bHat;
			if (braggIndex == 0){
				aHat = Distance.unitVector(recip.getA());
				bHat = Distance.unitVector(recip.getB());
			}
			else{
				aHat = Distance.unitVector(recip.getB());
				bHat = Distance.unitVector(recip.getA());
			}
			
			int[] cent = bragg[braggIndex];
			cent[0] += nx/2;
			cent[1] += ny/2;
			
			boolean[][] ans = new boolean [nx][ny];
			int xc, yc, xa, ya;
			double acent, abragg, bbragg;
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++){
					xc = i-nx/2;
					yc = j-ny/2;
					xa = i-cent[0];
					ya = j-cent[1];
					acent = Distance.dot(aHat[0], aHat[1], xc, yc);
					abragg = Distance.dot(aHat[0], aHat[1], xa, ya);
					bbragg = Distance.dot(bHat[0], bHat[1], xa, ya);
					ans[i][j] = (abragg*abragg)/(a*a) + (bbragg*bbragg)/(b*b) <= 1 && acent < cutoff;
				}
			return ans;
		}
	}
	
}
