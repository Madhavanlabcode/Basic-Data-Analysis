package util;

import image.GifSequenceWriter;
import image.ImageEditing;
import impurity.PointImp;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import org.apache.commons.math3.stat.regression.RegressionResults;
import org.apache.commons.math3.stat.regression.SimpleRegression;

import drawing.GraphDrawerCart;
import drawing.LayerViewer;
import drawing.SpectraDrawer;
import drawing.TopomapViewer;
import main.SRAW;
import misc.FermiVelocityCalculator;
import schrodinger.MovieMaker;
import util.calc.StripCut;
import util.color.ColorScales;
import util.fileops.ColumnIO;
import util.fileops.FileOps;
import util.fileops.Layer;
import util.fileops.LineOfSpectra;
import util.fileops.PointSpectra;
import util.fileops.RHKFileOps;
import util.fileops.Topomap;
import util.fourier.AtomicCoordinatesGenerator;
import util.fourier.FFT2DSmall;
import util.fourier.FFT3D_Wrapper;
import util.fourier.FFTOps;
import util.fourier.ImpurityListEditor;
import util.geom.AtomicCoordinatesSet;
import util.geom.Distance;
import util.matrix.Matrix;
import util.regression.ACM_CustomFunctions;
import util.regression.ACM_NonLinearFitter;
import util.regression.ACM_NonLinearFitter.FittingResult;
import util.robot.NanonisTopomapTaker.NanonisNumber;
import util.scalar.functions.Functions;
import flanagan.analysis.Regression;

public class TopomapUtil {
	static JFileChooser fc;
	
	public static void main(String[] args)
	{
		
		Topomap.setStdDir();
		fc = new JFileChooser(Topomap.stddir);
		
//		doAdvancedGaussianBinning(Topomap.open(fc));
		doAdvancedGaussianBinning(Topomap.open(fc), 10);
		System.exit(0);
//		mergeSomeMaps();
//		String basename = "_layer_";
//		Topomap tm = Topomap.open(fc);
//		String dirm = FileOps.selectDir(fc);
//		String[] mapFiles = new String [tm.nlayers];
//		for (int i = 0; i < tm.nlayers; i++)
//			mapFiles[i] = dirm + basename + i + ".bin";
//		AtomicCoordinatesSet lattm = new AtomicCoordinatesSet(FileOps.openText(fc));
		
//		writeCutsOfAllMapsInAFolder(mapFiles, tm.v, lattm, 32, 1, FileOps.selectSave(fc).toString());
//		System.exit(0);
		
		File f = FileOps.selectOpen(fc);
//		Layer t = Layer.readBIN(f.toString());
		Topomap t = Topomap.readBIN(f.toString());
//		String s = FileOps.selectSave(fc).toString();
//		doPseudofieldStuff(s);
//		Topomap t = Topomap.openFromLayersAuto(fc);
		String name = fc.getSelectedFile().getName().substring(0, fc.getSelectedFile().getName().length()-4);
		String dir = fc.getCurrentDirectory().toString() + "\\";
		AtomicCoordinatesSet latt = new AtomicCoordinatesSet(FileOps.openText(fc));
//		MapCuttingMethods.writeHighSymmetryCuts(t, latt, 32, 1, dir, name);
//		MapCuttingMethods.writeRadialCutRt2(t, latt, 40, 1, dir, name, 0);
//		MapCuttingMethods.writeWedgeCutsRt2(t, latt, 40, 1, dir, name, 0);
//		MapCuttingMethods.writeRadialCutsOrigin(t, 2, dir, name);
		MapCuttingMethods.writeCutsToBothBraggPeaks(t, latt, 32, 1, dir, name, false, false);
//		doCompleteAdvancedBinning(t, s);
//		doCompleteAdvancedBinning_SinglePath(t, s);
//		getAverageOfAdvancedBinning(t, s);
//		PointSpectra.writeBIN(ps, fc); 
//		System.exit(0);
		//I assume i have opened uaa:
//		subtractPolynomialFitFrom2ndDerivativeFolder(dir, name.substring(0,name.length()-3), 4);
//		System.out.println(s.substring(0, s.length() - 7));
//		writeTroughTopomaps(t, latt, s);
		
		//Now let us create the thing which will represent the bins: the absolute value
		//of the remainder when phi is interegerd. At the integer values this should be zero.
//		double[][][] phiTroughs = new double [16][t.nx][t.ny];
//		for (int k = 0; k < 16; k++)
//		for (int i = 0; i < t.nx; i++)
//			for (int j = 0; j < t.ny; j++)
//			{
//				phiTroughs[k][i][j] = Math.abs(phi[i][j] - FieldOps.round(phi[i][j]+k/16.0) + k/16.0);
//			}
//		LayerViewer.show(Layer.getFreeLayer(phiTroughs), 1000, true);
//		new TopomapViewer(new Topomap(phiTroughs, ArrayOps.generateArrayNotInclUpper(0, 15, 15), t.x, t.y, null), Topomap.stddir, 1000);
//		GraphDrawerCart.plotGraph(xAxis, averageAlongA);
//		GraphDrawerCart.plotGraph(xAxis, averageSmoothed);
//		System.out.println("" + maxas.length + " maxima found.");
		
//		boolean[][] troughsA = FieldOps.isDirectionalExtremum(t.data, Distance.unitVector(latt.getB()), 0, 0);
//		ImageEditing.copyToClipboard(ImageEditing.getBufferedImage(troughsA));
//		double[][] dispA = ArrayOps.toDouble(troughsA);
//		Topomap tA = Topomap.newTopomap(t, dispA);
//		LayerViewer.show(Layer.getFreeLayer(dispA), 1000, true);
//		new TopomapViewer(tA, dir, 1000);
		//		Topomap.writeBIN(t, fc);
//		Layer shift = Layer.open(fc);
//		double[][] temp = new double [t.nx][t.ny];
//		for (int i = 0; i < t.nx; i++)
//			for (int j = 0; j < t.ny; j++)
//				temp[i][j] = (i-250)*((t.v[1]-t.v[0])/250);
//		FieldOps.negate(shift.data);
//		Topomap shifted = TopomapUtil.shiftEachSpectrum(t, shift.data);
//		Topomap shifted = TopomapUtil.getFourierSpectrumLineFit(t);
//		Topomap.writeBIN(shifted, fc);
//		t.resizeSimple(1/Math.sqrt(2));
//		Topomap.writeBIN(t, f.toString());
//		System.exit(0);
//		double[][] data = new double [t.nx][t.ny];
//		for (int i = 0; i < t.nx; i++)
//			for (int j = 0; j < t.ny; j++)
//				data[i][j] = ArrayOps.min(t.getSpectrum(i, j));
//		Layer.writeBIN(new Layer(data, t.x, t.y, 1, 1), fc);
//		Topomap minSource = Topomap.open(fc);
//		writeShiftedAccordingToMinimum(t, minSource, fc);
//		Topomap t = Topomap.openFromLayers(fc);
//		flipCertainLayersY(t, new int[] {0, 2, 4, 6});
//		for (int i = 0; i < t.nlayers; i++)
//			FieldOps.changeZeroToAverage(t.data[i]);
//		double[][] min = new double [t.nx][t.ny];
//		for (int i = 0; i < t.nx; i++)
//			for (int j = 0; j < t.ny; j++)
//				min[i][j] = t.v[ArrayOps.maxIndex(t.getSpectrum(i, j))];
		
//		Layer.writeBIN(new Layer(min, t.x, t.y, 1, 1), fc);
//		Topomap.writeBIN(t, fc);
//		System.exit(0);
//		
//		double[] tmv = new double[20];
//		for (int i = 0; i < 10; i++)
//			tmv[i] = 0.1 - i*0.01;
//		for (int i = 10; i < tmv.length; i++)
//			tmv[i] = -0.01 - (i-10)*0.01;
//		ArrayOps.flip(tmv);
//		writeNanonisNumberAveragesCommas(t, Double.parseDouble(JOptionPane.showInputDialog("Multiplicitave factor?")), true, tmv);
//		for (int i = 0; i < t.nlayers; i++)
//			System.out.println(FieldOps.mean(t.data[i]));
//		for (int i = 0; i < t.nlayers; i++)
//			System.out.println("" + t.v[i] + "\t" + FieldOps.mean(t.data[i]));// + "\t"+ NanonisNumber.parseDouble(t.data[i]).toString());
//		for (int i = 0; i <= t.nlayers; i++)
//			System.out.println("" +(t.v[i+1]- t.v[i]));
		
//		Layer[] fits = fitAllSpectraToParabola_Double(t, -0.221, 0.008, 1);
//		LayerViewer.show(fits[0], 512, true);
//		LayerViewer.show(fits[1], 512, true);
//		doLandauLevelFitting(t, -0.2097, 0.006);
		//		Layer[] fits = adaptivelyFitAllSpectraToAParabola(t, -0.2206, 0.005);
//		LayerViewer.show(fits[0], 512, true);
//		LayerViewer.show(fits[1], 512, true);
//		LayerViewer.show(fits[2], 512, true);
		
//		Layer[] last = new Layer[fits.length-2];
//		for (int i = 0; i < fits.length-2; i++)
//			last[i] = fits[i+2];
//		new TopomapViewer(Topomap.newTopomap(last), fc.getCurrentDirectory().toString() + "\\",512);
//		Topomap z = Topomap.open(fc);
//		Topomap zadj = TopomapUtil.getZAdjusteddIdV(t, z, 1.91e10);
//		Topomap.writeBIN(zadj, fc);
		
//		Topomap.writeBIN(getZMap(t), fc);
//		Topomap.writeBIN(getWindowMapAboutV(t, 0.02), fc);
//		Topomap.writeBIN(getCurrentMapAboutV(t, 0), fc);
//		Layer[] minima = getSpectrumMinimum(t);
//		LayerViewer.show(minima[0], 1024);
//		LayerViewer.show(minima[1], 1024);
		
//		removeLockinZeroesAtExtremes(t, -0.262, true, 0.02);
//		invert(t);
//		Topomap.writeBIN(t, fc);
		
//		wrapIterativelySuppress(t);
		
//		saveCropped(t);
		
//		Topomap.writeBIN(Topomap.newTopomap(new Layer[] {Layer.open(fc)}), fc);
		
//		int nlayers = Integer.parseInt(JOptionPane.showInputDialog("How many layers?", 1));
//		Layer [] things = new Layer [nlayers];
//		for (int i =0; i < nlayers; i++)
//			things[i] = Layer.open(fc);
//		Topomap.writeBIN(Topomap.newTopomap(things), fc);
		
//		rotate180(t);
//		Topomap.writeBIN(t, fc);
		
//		doFourierFilterUserSpecified(t);
//		t.y = t.x;
//		Topomap.writeBIN(t, fc);
////		Layer height = Layer.open(fc);
//		

//		double[][] dist = MapCuttingMethods.getRadialRt2Point(t, latt, 1);
//		writeCutsToBothBraggPeaks(t, latt, 12, 2, dir, name);
//		MapCuttingMethods.writeHighSymmetryCuts(t, latt, 32, 1, dir, name);
//		MapCuttingMethods.writeRadialCutRt2(t, latt, 40, 1, dir, name, 0);
//		MapCuttingMethods.writeWedgeCutsRt2(t, latt, 40, 1, dir, name, 0);
//		MapCuttingMethods.writeRadialCutsOrigin(t, 2, dir, name);
//		MapCuttingMethods.writeCutsToBothBraggPeaks(t, latt, 32, 2, dir, name, false);
//		int npts = 80, maxLength = 80;
//		int braggIndex = 0;
//		double[][] distance = MapCuttingMethods.getStripInwardFromBraggPeak(t, latt, braggIndex, 32);//MapCuttingMethods.getInwardRightAngleRt2Point(t, latt, 0);
//		LayerViewer.show(Layer.getFreeLayer(distance), 1024, true);
//		boolean[][][] map = MapCuttingMethods.getMasks(distance, npts, maxLength);//MapCuttingMethods.getRadialRt2Masks(t, latt, 0, maxLength, npts);
//		MapCuttingMethods.maskMap(t, map, ArrayOps.generateArrayNotInclLower(0, maxLength*(t.xLength/t.nx), npts), dir, name, "bragg" + (braggIndex+1));
//		double[][][] mapd = new double [map.length][][];
//		for (int i = 0; i < map.length; i++) mapd[i] = ArrayOps.toDouble(map[i]);
		
//		new TopomapViewer(new Topomap(mapd, ArrayOps.generateArrayNotInclLower(0, maxLength, npts), t.x, t.y, null), null, 512);
		
//		splitTopomapByBinsOfLayer(t, dir, name, 4, fc);
//		splitTopomapSpectraByBinsOfLayer(t, dir, name, 4, t.nx, fc);
		
//		Topomap.writeBIN(split[1], dir + name + "split1.bin");
//		Topomap.writeBIN(shiftLeftOnePixel(Topomap.open()));
//		doRadialStripCut(fc.getCurrentDirectory().toString(), "\\", 2, 1000, t, "1PeakLorentzian_Line");
//		PointSpectra ps = getRadialAverageStripCut(2, t);
//		PointSpectra.writeBIN(ps, fc);
//		Layer.writeBIN(ps.toLayer(), fc);
		
		//		doRadialCutPeakFitting(fc.getCurrentDirectory().toString() + "\\new fits\\", "", s, "1PeakLorentzian_Line", 9, 2);		
//		doRadialCutPeakTwoPass(fc.getCurrentDirectory().toString(), "\\iterative fits\\2\\", s, "Lorentzian");		

//		boolean[][][] filter = doRadialCutPeakTwoPass(s, "Lorentzian");
//		Layer t = Layer.open(fc);
//		double[][] dist = ImpurityUtil.getYukawa_Power_FromAll(PointImp.readFromGaussSquareFile(FileOps.selectOpen(fc)), t.nx, t.ny, 0.25, 1, 12);
//		double[][] dist = ImpurityUtil.getDistanceTo(PointImp.readFromGaussSquareFile(FileOps.selectOpen(fc)), t.nx, t.ny);

//		double[][] dist = ImpurityUtil.countImpuritiesGauss(PointImp.readFromGaussSquareFile(FileOps.selectOpen(fc)), t.nx, t.ny, 24);
//		double[][] dist = ImpurityUtil.countImpuritiesDensity(PointImp.readFromGaussSquareFile(FileOps.selectOpen(fc)), t.nx, t.ny, 24);
//		Layer.writeBIN(Layer.newLayer(t, dist), fc);
		
		//Make a topomap of density maps with different lengths
//		String input = JOptionPane.showInputDialog("Enter radii min and max, spearated by commas.");
//		double rmin = Double.parseDouble(input.split(",")[0].trim());
//		double rmax = Double.parseDouble(input.split(",")[1].trim());
//		double[] radii = ArrayOps.generateArrayInclBoth(rmin, rmax, Integer.parseInt(JOptionPane.showInputDialog("How many layers?")));
//		if (radii.length > 1)
//			Topomap.writeBIN(new Topomap(ImpurityUtil.countImpuritiesGauss(PointImp.readFromGaussSquareFile(FileOps.selectOpen(fc)), t.nx, t.ny, radii), radii, t.x, t.y, null));
//		else
//			Layer.writeBIN(new Topomap(ImpurityUtil.countImpuritiesGauss(PointImp.readFromGaussSquareFile(FileOps.selectOpen(fc)), t.nx, t.ny, radii), radii, t.x, t.y, null).getLayer(0));
		
		//Apply u-field, generated from DriftCorrectionMethods' AUTO function, to a topomap.
//		Topomap tu = applyUField(DriftCorrectionAnalysis.getUFromAutoDriftCorrFolder(null, fc), t);
//		Topomap.writeBIN(tu, fc);
		
		//Apply u-field to a list of impurities
//		PointImp[] imps = PointImp.readFromGaussSquareFile(FileOps.selectOpen(fc));
//		int[] energyIndices = PointImp.readEnergyIndices(fc.getSelectedFile());
//		PointImp[] shift = ImpurityUtil.shiftWithUField(imps, DriftCorrectionAnalysis.getUFromAutoDriftCorrFolder(null, fc));
//		PointImp.writeToFile(shift, FileOps.selectSave(fc));
//		String mapname = fc.getSelectedFile().getName();
//		mapname = mapname.substring(0, mapname.length()-4);
//	
//		ImageEditing.writeForAVIMovie(ImageEditing.getImpurityModuloMapMovie(ImpurityUtil.getLatticePositionsModuloN(imps, 1, 0, 0, latt), imps, t, energyIndices), mapname);
	
//		PointSpectra ps = ImpurityUtil.getSpectraAt(imps, t, true, 5);
//		PointSpectra ps = ImpurityUtil.getLocalDeviationSpectraNorm(imps, t, 1, 2);
//		PointSpectra ps = ImpurityUtil.getLocalPercentiles(imps, t, 1.5);
//		PointSpectra ps = PointSpectra.open(fc);
//		double[] maxima = new double [ps.nspec];
//		int[] maxi = new int [ps.nspec];
//		for (int i = 0; i < ps.nspec; i++)
//		{
//			maxi[i] = ArrayOps.maxIndex(ps.data[i]);
//			maxima[i] = t.v[maxi[i]];
//			System.out.println(maxima[i]);
//		}
//		double de = t.v[1]-t.v[0];
//		double[] bins = ArrayOps.generateArrayNotInclUpper(t.v[0] - de/2, t.v[t.nlayers-1] - de/2, t.nlayers);
//		int[] hist = ArrayOps.getHistogram(maxima, ArrayOps.generateArrayNotInclUpper(ArrayOps.min(maxima), ArrayOps.max(maxima), bins.length));
////		boolean[] greaterThanZero = new boolean [imps.length];
////		for (int i = 0; i < imps.length; i++)
////			greaterThanZero[i] = maxima[i] > 0;
////		PointImp[][] split = ImpurityUtil.splitImps(imps, greaterThanZero);
////		PointImp.writeToFile(split[0], new File(fc.getSelectedFile().toString() + "N_A.txt"));
////		PointImp.writeToFile(split[1], new File(fc.getSelectedFile().toString() + "N_B.txt"));
//	
//		for (int i = 0; i < bins.length; i++)
//			System.out.println("" + bins[i] + "\t" + hist[i]);
//		GraphDrawerCart.plotGraph(bins, hist);
//		PointSpectra.writeBIN(ps, fc);
//		
//		//Write the energy-dependent list of impurities
//		String[] lines = new String [imps.length+3];
//		lines[0] = "Impurity List";
//		lines[1] = "";
//		lines[2] = "X\tY\tMaxIndex";
//		for (int i = 0; i < imps.length; i++)
//			lines[i+3] = "" + imps[i].pixelPos[0] + "\t" + imps[i].pixelPos[1] + "\t" + maxi[i];
//		
//		FileOps.writeLines(fc, lines);
//		PointImp[] fitted = new PointImp[imps.length];
//		double[][][] temp = TwoDGaussianFreeFitter.fitListOfImpurities(imps, t.data[24], 2);
//		double[][] table = temp[0];
//		double[][] value = temp[1];
//		for (int i = 0; i < fitted.length; i++){
//			fitted[i] = new PointImp(new double [] {table[i][0], table[i][1]});
//		}
//		PointImp.writeToFile(fitted, FileOps.selectSave(fc));
//		LayerViewer.show(Layer.newLayer(t, value), 512);
//		LayerViewer.show(Layer.newLayer(t, temp[2]), 512);
		//		double[][] ints = new double [t.nx][t.ny];
//		int[][] bins = FieldOps.getPercentileBinsForField(field, nbins)
		//		Topomap t = Topomap.open(fc);
		
//		Topomap[] mask = impurityMaskPercentileHighFreq(t, 0.02, 0.98, 30);
//		Topomap.writeBIN(mask[0], fc);
//		Topomap.writeBIN(mask[1], fc);
		
		//		for (int i = 0; i < t.nlayers; i++)
//			FieldOps.changeZeroToAverage(t.data[i]);
//		Topomap.writeBIN(t);
//		Topomap tp = applyUField(FileOps.open2Comp(true), t);
//		Topomap.writeBIN(tp);
//		Topomap.writeBIN(getJDOS(t), fc);
		//		Layer l = Layer.openFree(fc);
//		Topomap.writeBIN(addTopographyTo(t, l));
//		Topomap subt = subtractBackground(t, 5)[0];
//		Topomap.writeBIN(subt, fc);
//		Topomap.writeBIN(truncateBias(t, 16, t.nlayers-12));
//		fitEachSpectrum(t, fc.getCurrentDirectory().toString() + "\\Curve Fit\\", "ShiftedExp");
//		fitEachSpectrumExponential(t, fc.getCurrentDirectory().toString() + "\\Curve Fit\\");
//		Topomap.writeBIN(getLocalNematicity(t, latt), fc);
		PointSpectra[] onoff = doHighSymmetryStripCut(latt, dir + "\\", name, Double.parseDouble(JOptionPane.showInputDialog("Width of line?", "" + 2)), t);
		
//		PointSpectra stripCut = doStripCutAcross_Rt2(latt, dir + "\\", name + "_rt2", Double.parseDouble(JOptionPane.showInputDialog("Width of line?", "" + 2)), t);
		
//		PointSpectra ivsk = PointSpectra.open(fc);
//		doIvsKPeakFittingStandard(fc.getCurrentDirectory().toString() + "\\", fc.getSelectedFile().getName().substring(0, fc.getSelectedFile().getName().length()-4), ivsk, 20, 56, true, 0);
//				FileOps.writeString(fc, latt.getRt2Lattice().toString());
//		Topomap.writeBIN(makeSymmetrizedFFTs(t, latt));
//		Topomap s = makeSymmetrizedFFTs(t, latt, false);
//		Topomap.writeBIN(s, fc);
//		new TopomapViewer(s, fc.getCurrentDirectory().toString(), 512);
//		Topomap.writeBINsImages(t, false);
//		Topomap tm = expandBi(t, 2);
//		Topomap.writeBIN(tm);
//		ColumnIO.writeTable(getSpectraAtLatticeAndOff(t, latt), FileOps.selectSave(fc).toString());
//		PointSpectra p = getSpectraAtLattice(t, latt, new double[] {0, 0});
//		PointSpectra.writeBIN(p, fc);
		
//		LineOfSpectra[] ls = splitTopomapIntoLinesOfSpectra(t.nx > t.ny, t, height);
//		String filepath = FileOps.selectSave(fc).toString();
//		for (int i = 0; i < ls.length; i++)
//		{
//			LineOfSpectra.writeBIN(ls[i], filepath + "_" + i + ".bin");
//		}
	}
	
	public static void doAdvancedGaussianBinning(Topomap t){
		Topomap phia = Topomap.open(fc);
		Topomap phib = Topomap.open(fc);
		String output = FileOps.selectSave(fc).toString();
		int nx = phia.nx, ny = phia.ny;
		double[][] distanceFirst = new double [nx][ny];//FieldOps.getQuadratureSum(phia.data[0], phib.data[0]);
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				distanceFirst[i][j] = Distance.distance(phia.data[0][i][j] - FieldOps.round(phia.data[0][i][j]),
				phib.data[0][i][j] - FieldOps.round(phib.data[0][i][j]));
//		LayerViewer.show(Layer.getFreeLayer(distanceFirst), 1000, true);
		double L = FieldOps.getCorrectGaussianLength(distanceFirst, 30, 0.2);
//		double[][] gauss= FieldOps.getGaussianWeights(distanceFirst, L);
//		LayerViewer.show(Layer.getFreeLayer(gauss), 1000, true);
		int nbins = 65;
		double[] divOffset = ArrayOps.generateArrayInclBoth(0, 2, nbins);
		double[] anisOffsetA = ArrayOps.generateArrayInclBoth(-0.5, 1.5, nbins);
		double[] anisOffsetB = ArrayOps.generateArrayInclBoth(0, 2, nbins);
		for (int k = 0; k < t.nlayers; k++)
		{
//			Topomap weights = getGaussianWeightsForLayer(phia, phib, k, divOffset, divOffset, L);
			Topomap weights = getGaussianWeightsForLayer(phia, phib, k, anisOffsetA, anisOffsetB, L);
			Topomap smoothed = doGeneralAdvancedBinning(k, t, weights);
			Topomap.writeBIN(smoothed, output + "_layer_" + k + ".bin");
		}
	}
	public static void doAdvancedGaussianBinning(Topomap t, int k){
		Topomap phia = Topomap.open(fc);
		Topomap phib = Topomap.open(fc);
		String output = FileOps.selectSave(fc).toString();
		boolean saveWeights = JOptionPane.showConfirmDialog(null, "Save the weights?") == JOptionPane.YES_OPTION;
		boolean saveAvg = JOptionPane.showConfirmDialog(null, "Save the average spectrum?") == JOptionPane.YES_OPTION;
		int nx = phia.nx, ny = phia.ny;
		double[][] distanceFirst = new double [nx][ny];//FieldOps.getQuadratureSum(phia.data[0], phib.data[0]);
		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
				distanceFirst[i][j] = Distance.distance(phia.data[0][i][j] - FieldOps.round(phia.data[0][i][j]),
				phib.data[0][i][j] - FieldOps.round(phib.data[0][i][j]));
//		LayerViewer.show(Layer.getFreeLayer(distanceFirst), 1000, true);
		double L = FieldOps.getCorrectGaussianLength(distanceFirst, 30, 0.2);
//		double[][] gauss= FieldOps.getGaussianWeights(distanceFirst, L);
//		LayerViewer.show(Layer.getFreeLayer(gauss), 1000, true);
		int nbins = 65;
		double[] divOffset = ArrayOps.generateArrayInclBoth(0, 2, nbins);
		double[] anisOffsetA = ArrayOps.generateArrayInclBoth(-0.5, 1.5, nbins);
		double[] anisOffsetB = ArrayOps.generateArrayInclBoth(0, 2, nbins);
//		Topomap weights = getGaussianWeightsForLayer(phia, phib, k, divOffset, divOffset, L);
		Topomap weights = getGaussianWeightsForLayer(phia, phib, k, anisOffsetA, anisOffsetB, L);
		Topomap smoothed = doGeneralAdvancedBinning(k, t, weights);
		if (saveWeights) Topomap.writeBIN(weights, output + "_weigh_" + k + ".bin");
		Topomap.writeBIN(smoothed, output + "_layer_" + k + ".bin");
		if (saveAvg){
			String[] lines = new String [smoothed.nlayers];
			double[] smoothSpec = smoothed.getAverageSpectrum();
			for (int i = 0; i < lines.length; i++)
				lines[i] = "" + smoothed.v[i] + "\t" + smoothSpec[i];
			ColumnIO.writeLines(lines, output + "_averageSpec_" + k + ".txt");
		}
	}
	
	/**
	 * We are going to assume that the mapFiles and mapVoltages arrays are already sorted in ascending order of voltage.
	 * @param mapFiles
	 * @param mapVoltages
	 * @param latt
	 * @param dir
	 * @param width
	 * @param thickness
	 */
	public static void writeCutsOfAllMapsInAFolder(String[] mapFiles, double[] mapVoltages, AtomicCoordinatesSet latt, int width, double thickness, String output){
//		ArrayList<Layer> ls = new ArrayList<Layer>();
//		File[] f = new File(dir).listFiles();
		Topomap t = null; boolean firstTime = true;
		boolean[][][][] braggMasks = null;
		int nk = 0, ne = mapVoltages.length, nlayers = 0;
		double[] k = null;
		double[][][][] ans = null;
//		ArrayList<Double> energies = new ArrayList<Double>();
		for (int i = 0; i < mapFiles.length; i++)
		{
			System.out.println();
			System.out.println(mapFiles[i]);
			t = Topomap.readBIN(mapFiles[i].toString());
			Topomap fft; double[][][] fftmag = new double [t.nlayers][][];
			for (int r = 0; r < t.nlayers; r++) fftmag[r] = FFTOps.obtainFFTmagCent(t.data[r]);
			fft = Topomap.newTopomap(t, fftmag);
			if (firstTime){
//				double[][][] bragDist = {TopomapUtil.MapCuttingMethods.getStripOutwardTowardsBraggPeak(t, latt, 0, width),
//						TopomapUtil.MapCuttingMethods.getStripOutwardTowardsBraggPeak(t, latt, 1, width)};
				braggMasks = new boolean [2][][][];
//				braggMasks[0] = TopomapUtil.MapCuttingMethods.getMasksThickness(bragDist[0], thickness);
//				braggMasks[1] = TopomapUtil.MapCuttingMethods.getMasksThickness(bragDist[1], thickness);
				braggMasks = TopomapUtil.MapCuttingMethods.getMapsToBothBraggPeaks(t, latt, width, (int) thickness);
				nk = braggMasks[0].length;
				nlayers = t.nlayers;
				k = t.getKValuesReal(nk);
				ans = new double [2][nlayers][nk][ne];
				firstTime = false;
			}
			//Now let us operate upon the topomap which will give us a single horizontal slab of the final product.
			for (int p = 0; p < nk; p++){System.out.print(" " + p);
				double[][] spectrum = new double [2][];
				for (int m = 0; m < 2; m++){
					spectrum[m] = TopomapUtil.MapCuttingMethods.getSpectrumAveragedOverBooleans(fft, braggMasks[m][p]);
					for (int q = 0; q < t.nlayers; q++){
						ans[m][q][p][i] = spectrum[m][q];
					}
				}
			}
		}
		
		Topomap[] maps = new Topomap[2];
		for (int i = 0; i < maps.length; i++){
			maps[i] = new Topomap(ans[i], t.v, k, mapVoltages, null);
			Topomap.writeBIN(maps[i], output + "_bragg_" + i + ".bin");
		}
	}
	public static void printCentroidOfSnTeEllipse(Topomap t, AtomicCoordinatesSet latt){
		double[] aHat = Distance.unitVector(latt.getA()), bHat = Distance.unitVector(latt.getB());
		//		String s = FileOps.selectSave(fc).toString();
//		int degree = Printer.getAnInt("Enter the degree of polynomial.");
		boolean[][] mask = LayerUtil.BooleanFields.getFilledCutEllipseAroundBraggPeak(t.nx, t.ny, latt, 120, 60, 280, 0);
		double[][] centroids = MapCuttingMethods.getFFTCentroidEachLayer(t, mask);
		for (int i = 0; i < t.nlayers; i++){
			double a = Distance.dot(aHat[0], aHat[1], centroids[i][0]-t.nx/2, centroids[i][1]-t.ny/2);
			double b = Distance.dot(bHat[0], bHat[1], centroids[i][0]-t.nx/2, centroids[i][1]-t.ny/2);
			System.out.println("" + i + "\t" + a + "\t" + b);
		}

	}
	public static void doCompleteAdvancedBinning(Topomap t, String s){
		Topomap phia = Topomap.open(fc);
		Topomap phib = Topomap.open(fc);
//		String dira = s + "phia\\";
//		String dirb = s + "phib\\";
		String dirdiv = s + "phithroughdiv\\";
		String diranis = s + "phiacrossanis\\";
//		if (!new File(dira).exists()) new File(dira).mkdir();
//		if (!new File(dirb).exists()) new File(dirb).mkdir();
		if (!new File(dirdiv).exists()) new File(dirdiv).mkdir();
		if (!new File(diranis).exists()) new File(diranis).mkdir();
		
		int nbins = 65;
		double gaussL = 2;
//		doRidiculouslyAdvancedBinning(t, phia, nbins, 6, 1, 5, dira);
//		doRidiculouslyAdvancedBinning(t, phib, nbins, 6, 1, 5, dirb);
		
		double[] divOffset = ArrayOps.generateArrayInclBoth(0, 2, nbins);
		doEvenMoreRidiculouslyAdvancedBinning(t, phia, phib, gaussL, divOffset, divOffset, 5, dirdiv);
		double[] anisOffsetA = ArrayOps.generateArrayInclBoth(-0.5, 1.5, nbins);
		double[] anisOffsetB = ArrayOps.generateArrayInclBoth(0, 2, nbins);
		doEvenMoreRidiculouslyAdvancedBinning(t, phia, phib, gaussL, anisOffsetA, anisOffsetB, 5, diranis);
	}
	public static void doCompleteAdvancedBinning_SinglePath(Topomap t, String s){
		Topomap phia = Topomap.open(fc);
		Topomap phib = Topomap.open(fc);
		String dirpath = s + "thepath\\";
		if (!new File(dirpath).exists()) new File(dirpath).mkdir();
		
		double gaussL = 6;
		
		double[] pathElement1A = ArrayOps.generateArrayNotInclUpper(0, 1, 32);
		double[] pathElement1B = ArrayOps.generateArrayNotInclUpper(0, 0, 32);
		double[] pathElement2A = ArrayOps.generateArrayNotInclUpper(1, 1, 32);
		double[] pathElement2B = ArrayOps.generateArrayNotInclUpper(0, 1, 32);
		double[] pathElement3A = ArrayOps.generateArrayInclBoth(1, 2, 41);
		double[] pathElement3B = ArrayOps.generateArrayInclBoth(1, 2, 41);
		double[] pathA = ArrayOps.concatenate(new double[][] {pathElement1A, pathElement2A, pathElement3A});
		double[] pathB = ArrayOps.concatenate(new double[][] {pathElement1B, pathElement2B, pathElement3B});
		//		double[] divOffset = ArrayOps.generateArrayInclBoth(0, 2, nbins);
		doEvenMoreRidiculouslyAdvancedBinning(t, phia, phib, gaussL, pathA, pathB, 5, dirpath);
	}
	public static void getAverageOfAdvancedBinning(Topomap t, String s){
		Topomap phia = Topomap.open(fc);
		Topomap phib = Topomap.open(fc);
		int nbins = 65;
//		doRidiculouslyAdvancedBinning(t, phia, nbins, 6, 1, 5, dira);
//		doRidiculouslyAdvancedBinning(t, phib, nbins, 6, 1, 5, dirb);
		
//		double[] divOffset = ArrayOps.generateArrayInclBoth(0, 2, nbins);
//		doEvenMoreRidiculouslyAdvancedBinning(t, phia, phib, 6, divOffset, divOffset, 5, dirdiv);
		double[] anisOffsetA = ArrayOps.generateArrayInclBoth(-0.5, 1.5, nbins);
		double[] anisOffsetB = ArrayOps.generateArrayInclBoth(0, 2, nbins);
		PointSpectra ps =  getAverageOfEvenMoreRidiculouslyAdvancedBinning(t, phia, phib, 2, anisOffsetA, anisOffsetB, 5);
		PointSpectra.writeBIN(ps, s);
	}
	public static void doRidiculouslyAdvancedBinning(Topomap t, Topomap phi, int nbins, double smL, double offsetMax, int recipPerc, String s){
		double[] offset = ArrayOps.generateArrayInclBoth(0, offsetMax, nbins);
		for (int k = 0; k < t.nlayers; k++){System.out.println(k);
			double[][][] phiTroughs = new double [nbins][t.nx][t.ny];
			for (int q = 0; q < nbins; q++){
				for (int i = 0; i < t.nx; i++)
					for (int j = 0; j < t.ny; j++)
						phiTroughs[q][i][j] = Math.abs(phi.data[k][i][j] - FieldOps.round(phi.data[k][i][j]+offset[q]) + offset[q]);
				int[][] bins = FieldOps.getPercentileBinsForField(phiTroughs[q], recipPerc);
				double[][] smoothWeight = FieldOps.getSmoothedWeightingFunctions(bins, smL, 0);
				phiTroughs[q] = FieldOps.spatialFilter(t.data[k], smoothWeight);
			}
			Topomap.writeBIN(new Topomap(phiTroughs, offset, t.x, t.y, null), s + "_layer_" + k + ".bin");
		}
	}
	public static void doEvenMoreRidiculouslyAdvancedBinning(Topomap t, Topomap phia, Topomap phib, double smL, double[] offsetA, double[] offsetB, int recipPerc, String s){
		int nbins = offsetA.length;
		double[][][] phiTroughs = new double [nbins][t.nx][t.ny];
//		for (int k = 0; k < t.nlayers; k++){System.out.println(k);
		for (int k = 6; k < 7; k++){System.out.println(k);
			for (int q = 0; q < nbins; q++){
				for (int i = 0; i < t.nx; i++)
					for (int j = 0; j < t.ny; j++)
						phiTroughs[q][i][j] = Math.abs(phia.data[k][i][j] - FieldOps.round(phia.data[k][i][j]+offsetA[q]) + offsetA[q])
						+ Math.abs(phib.data[k][i][j] - FieldOps.round(phib.data[k][i][j]+offsetB[q]) + offsetB[q]);
				int[][] bins = FieldOps.getPercentileBinsForField(phiTroughs[q], recipPerc);
				double[][] smoothWeight = FieldOps.getSmoothedWeightingFunctions(bins, smL, 0);
				phiTroughs[q] = FieldOps.spatialFilter(t.data[k], smoothWeight);
			}
			Topomap.writeBIN(new Topomap(phiTroughs, ArrayOps.generateArrayInclBoth(0, nbins-1, nbins), t.x, t.y, null), s + "_layer_" + k + ".bin");
			System.gc();
		}
	}
	/**
	 * This splits t into nlayers topomaps each of which is one particular layer binned by the weights given in the Weights topomap.
	 * @param t
	 * @param weights
	 * @param s
	 */
	public static void doGeneralAdvancedBinning(Topomap t, Topomap weights, String s){
		int nbins = weights.nlayers;
		for (int k = 0; k < t.nlayers; k++){//System.out.println(k);
			double[][][] binned = new double [nbins][][];
			for (int q = 0; q < nbins; q++){
				binned[q] = FieldOps.spatialFilter(t.data[k], weights.data[q]);
			}
			Topomap.writeBIN(new Topomap(binned, ArrayOps.generateArrayInclBoth(0, nbins-1, nbins), t.x, t.y, null), s + "_layer_" + k + ".bin");
			System.gc();
		}
	}
	public static Topomap doGeneralAdvancedBinning(int k, Topomap t, Topomap weights){
		int nbins = weights.nlayers;
		double[][][] binned = new double [nbins][][];
		for (int q = 0; q < nbins; q++){
			binned[q] = FieldOps.spatialFilter(t.data[k], weights.data[q]);
		}
		return Topomap.newTopomap(weights, binned);
	}
	public static Topomap getGaussianWeightsForLayer(Topomap phia, Topomap phib, int k, double[] offsetA, double[] offsetB, double L){
	
		int nbins = offsetA.length;
		int nx = phia.nx, ny = phib.ny;
		//First find the gaussian length;
//		double L = -1;
		double[][] distance = new double [nx][ny];
		double[][][] weights = new double [nbins][nx][ny];
		for (int q = 0; q < nbins; q++){
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
					distance[i][j] = Distance.distance(phia.data[k][i][j] - FieldOps.round(phia.data[k][i][j]+offsetA[q]) + offsetA[q],
					phib.data[k][i][j] - FieldOps.round(phib.data[k][i][j]+offsetB[q]) + offsetB[q]);
			weights[q] = FieldOps.getGaussianWeights(distance, L);
		}
		return new Topomap(weights, ArrayOps.generateArray(0, 1, nbins), phia.x, phia.y, null);
	}
	public static PointSpectra getAverageOfEvenMoreRidiculouslyAdvancedBinning(Topomap t, Topomap phia, Topomap phib, double smL, double[] offsetA, double[] offsetB, int recipPerc){
		int nbins = offsetA.length;
		//Here nspec will be the number of layers of the topograph,
		//and each spectrum will represent a spectrum over Q.
		double[][] specdata = new double [t.nlayers][nbins];
		for (int k = 0; k < t.nlayers; k++){System.out.println(k);
			double[][][] phiTroughs = new double [nbins][t.nx][t.ny];
			for (int q = 0; q < nbins; q++){
				for (int i = 0; i < t.nx; i++)
					for (int j = 0; j < t.ny; j++)
						phiTroughs[q][i][j] = Math.abs(phia.data[k][i][j] - FieldOps.round(phia.data[k][i][j]+offsetA[q]) + offsetA[q])
						+ Math.abs(phib.data[k][i][j] - FieldOps.round(phib.data[k][i][j]+offsetB[q]) + offsetB[q]);
				int[][] bins = FieldOps.getPercentileBinsForField(phiTroughs[q], recipPerc);
				double[][] smoothWeight = FieldOps.getSmoothedWeightingFunctions(bins, smL, 0);
				specdata[k][q] = FieldOps.mean(t.data[k], smoothWeight);
			}
		}
		return new PointSpectra(specdata, offsetA, t.v, new double [t.nlayers]);
	}
	public static void writeTroughTopomaps(Topomap t, AtomicCoordinatesSet latt, String s)
	{
		
		for (int p = 0; p < 2; p++){
			double[] vHat = Distance.unitVector(p == 0 ? latt.getA() : latt.getB());
			double[][][] phi = new double [t.nlayers][t.nx][t.ny];
			double[][][] intsd = new double [t.nlayers][t.nx][t.ny];
			for (int q = 0; q < t.nlayers; q++){System.out.println("" + p + "\t" + q);
				double[][] faa = FieldOps.getSecondDirectionalDerivative2D(t.data[q], vHat);
				double thickness = 2;
				double smooth = 3;
				double[] averageAlongA = FieldOps.getSpectrumWithBools(faa, TopomapUtil.MapCuttingMethods.getMasksThickness(TopomapUtil.MapCuttingMethods.getUnrestrictedDotWithDirection(t.nx, t.ny, latt, p), thickness));
				int npts = averageAlongA.length; int margin = 16;
				double[] averageSmoothed = ArrayOps.gaussSmooth(averageAlongA, smooth);
				double[] derivativeSmoothed = ArrayOps.getDerivative(averageSmoothed);
				double[] xAxis = ArrayOps.generateArray(-averageAlongA.length*thickness/2, thickness, averageAlongA.length);
				ArrayList<double[]> maxima = new ArrayList<double[]>();
				for (int i = margin; i < npts - margin; i++)
					if (derivativeSmoothed[i] > 0 && derivativeSmoothed[i+1] < 0 && averageSmoothed[i] > 0 && averageSmoothed[i+1] > 0){
						//Now we fit to a parabola the thing, and place the maximum there.
		//				double[] tempX = ArrayOps.copyWithin(xAxis, i-margin/2, i+margin/2);
		//				double[] tempY = ArrayOps.copyWithin(averageSmoothed, i-margin/2, i+margin/2);
						SpectraUtil.SpectraFittingRange sfr = new SpectraUtil.SpectraFittingRange(i-margin/2, margin, xAxis, averageSmoothed);
						FittingResult r = sfr.fitToParabola();
						double[] added = new double[] {r.fitParams[0], r.fitParams[2]};
//						System.out.print(Printer.arrayLnHorizontal(added));
						maxima.add(added);
					}
				
				//Now I have the list of maxima as a function of a. I may now proceed to label each pixel
				//both by the index of the nearest maximum below it, and by its phase on the way to the next one,
				//the maxima being regarded as 1 apart.
				double[] maxas = new double [maxima.size()];
				double[] intervals = new double [maxima.size()-1];
				for (int i = 0; i < maxas.length; i++){
					maxas[i] = maxima.get(i)[0];
					if (i < maxas.length-1) intervals[i] = maxima.get(i+1)[0] - maxima.get(i)[0];
				}
				double avgInterval = ArrayOps.mean(intervals);
				double[] expandedmaxas = new double [maxas.length+4];//This is so that we can always find a "maximum" greater and less than a. 
				expandedmaxas[0] = maxas[0] - 2*avgInterval;
				expandedmaxas[1] = maxas[0] - avgInterval;
				for (int i = 0; i < maxas.length; i++)
					expandedmaxas[i+2] = maxas[i];
				expandedmaxas[maxas.length+2] = maxas[maxas.length-1]+avgInterval;
				expandedmaxas[maxas.length+3] = maxas[maxas.length-1]+2*avgInterval;
				int[][] ints = new int [t.nx][t.ny];
//				double[][] phi = new double [t.nx][t.ny];
				for (int i = 0; i < t.nx; i++)
					for (int j = 0; j < t.ny; j++)
					{
						double a = Distance.dot(vHat[0], vHat[1], i-t.nx/2, j-t.ny/2);
		//				ints[i][j] = ArrayOps.indexOf(maxas, a, true);
						ints[i][j] = ArrayOps.indexOf(expandedmaxas, a, true);
						intsd[q][i][j] = ints[i][j];
						//				if (ints [i][j] == 0) System.out.println("" + i + "\t" + j);
						//Now determine the distance to the abscissae defined by ints[i][j] and ints[i][j]-1:
						phi[q][i][j] = (ints[i][j]-1) + (a - expandedmaxas[ints[i][j]-1])/(expandedmaxas[ints[i][j]] - expandedmaxas[ints[i][j]-1]);
					}
				
//				LayerViewer.show(Layer.getFreeLayer(ArrayOps.toDouble(ints)), 1000, true);
//				LayerViewer.show(Layer.getFreeLayer(phi), 1000, true);
			}
			Topomap tp = Topomap.newTopomap(t, phi);
			Topomap ti = Topomap.newTopomap(t, intsd);
			Topomap.writeBIN(ti, s + "ints_" + p + ".bin");
			Topomap.writeBIN(tp, s + "phi_" + p + ".bin");
		}
	}
	public static void writeSomeDistanceMaps_troughsAB(Topomap t, AtomicCoordinatesSet latt)
	{
		String s = FileOps.selectSave(fc).toString();
		boolean[][][] troughsA = FieldOps.isDirectionalExtremum(t.data, Distance.unitVector(latt.getA()), 0);
		boolean[][][] troughsB = FieldOps.isDirectionalExtremum(t.data, Distance.unitVector(latt.getB()), 0);
		double[][][] distanceA = new double [t.nlayers][t.nx][t.ny];
		double[][][] distanceB = new double [t.nlayers][t.nx][t.ny];
		double[][][] distanceApB = new double [t.nlayers][t.nx][t.ny];
		double[][][] distanceAmB = new double [t.nlayers][t.nx][t.ny];
		for (int i = 0; i < t.nlayers; i++){
			System.out.println(i);
			Distance.putDistanceFromTrueSparse(troughsA[i], distanceA[i]);
			Distance.putDistanceFromTrueSparse(troughsB[i], distanceB[i]);
			distanceApB[i] = FieldOps.add(distanceA[i], distanceB[i]);
			distanceAmB[i] = FieldOps.minus(distanceA[i], distanceB[i]);
		}
		Topomap tda = Topomap.newTopomap(t, distanceA);
		Topomap tdb = Topomap.newTopomap(t, distanceB);
		Topomap tdapb = Topomap.newTopomap(t, distanceApB);
		Topomap tdamb = Topomap.newTopomap(t, distanceApB);
		Topomap.writeBIN(tda, s + "dist_a.bin");
		Topomap.writeBIN(tdb, s + "dist_b.bin");
		Topomap.writeBIN(tdapb, s + "dist_apb.bin");
		Topomap.writeBIN(tdamb, s + "dist_amb.bin");

	}
	public static void mergeSomeMaps()
	{
		int n = Integer.parseInt(JOptionPane.showInputDialog("How many maps to merge?"));
		Topomap[] them = new Topomap[n];
		for (int i = 0; i < n; i++)
			them[i] = Topomap.open(fc);
		Topomap merged = mergeMaps(them);
		Topomap.writeBIN(merged, fc);
		
		System.exit(0);
	}
	public static void doPseudofieldStuff(String s)
	{
		Topomap comp = Topomap.open(fc);
		Topomap anis = Topomap.open(fc);
		Topomap e12 = Topomap.open(fc);
		double aC = 3.8, aU = 1.2, a3 = 0;
		double pixelsPerAngstrom = 1024.0/1300.0;
		double constant = FermiVelocityCalculator.vectorPotential_teslaAngstroms(1);
		int nl = comp.nlayers, nx = comp.nx, ny = comp.ny;
		
		double[][][] A10 = new double [nl][nx][ny];
		double[][][] A21 = new double [nl][nx][ny];
		double[][][] B0 = new double [nl][nx][ny];;
		double[][][] B1 = new double [nl][nx][ny];;
		
		for (int i = 0; i < comp.nlayers; i++){
			double[][][] temp = FieldOps.getPseudofieldStuff_alternate(comp.data[i], anis.data[i], e12.data[i], aC, aU, a3);
			A10[i] = temp[0];
			A21[i] = temp[3];
			B0[i] = temp[4];
			B1[i] = temp[5];
			FieldOps.timesEquals(B0[i], pixelsPerAngstrom*constant);
			FieldOps.timesEquals(B1[i], pixelsPerAngstrom*constant);
		}
		Topomap A10m = Topomap.newTopomap(comp, A10);
		Topomap A21m = Topomap.newTopomap(comp, A21);
		Topomap B0m = Topomap.newTopomap(comp, B0);
		Topomap B1m = Topomap.newTopomap(comp, B1);
		
		Topomap.writeBIN(A10m, s + "A1_x.bin");
		Topomap.writeBIN(A21m, s + "A2_y.bin");
		Topomap.writeBIN(B0m, s + "B0.bin");
		Topomap.writeBIN(B1m, s + "B1.bin");
		
		System.exit(0);
	}
	public static void doLandauLevelFitting(Topomap t, double vCent, double hWidth)
	{
		Layer[] fits = adaptivelyFitAllSpectraToAParabola(t, vCent, hWidth);
		LayerViewer.show(fits[0], 512, true);
		LayerViewer.show(fits[1], 512, true);
		LayerViewer.show(fits[2], 512, true);
	}
	public static Topomap getFourierSpectrumLineFit(Topomap t)
	{
		double[][][] data = FieldOps.copy(t.data);
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
			{
				double[] initSpec = t.getSpectrum(i, j);
				double m = (initSpec[t.nlayers-1] - initSpec[0])/(t.nlayers-1);
				for (int k = 0; k < t.nlayers; k++)
					data[k][i][j] = m*k;
			}
		return Topomap.newTopomap(t, data);
		
	}
	public static Layer[] adaptivelyFitAllSpectraToAParabola(Topomap t, double initCent, double hWidth)
	{
		int sx = t.nx/2, sy = t.ny/2;
		double[][] maxPos = new double [t.nx][t.ny];
		double[][] maxCurv = new double [t.nx][t.ny];
		double[][] maxRaw = new double [t.nx][t.ny];
		double[][] maxOffset = new double [t.nx][t.ny];
		boolean[][] visited = new boolean [t.nx][t.ny];
		SpectraUtil.SpectraFittingRange sfr = new SpectraUtil.SpectraFittingRange(initCent-hWidth, initCent+hWidth, t.v, t.getSpectrum(sx, sy));
		FittingResult r = sfr.fitToParabola();
		maxPos[sx][sy] = r.fitParams[0];
		maxCurv[sx][sy] = r.fitParams[1];
		maxRaw[sx][sy] = r.x[ArrayOps.maxIndex(r.y)];
		maxOffset[sx][sy] = r.fitParams[2];
		boolean tryingToGetMinimum = maxCurv[sx][sy] > 0;

		visited[sx][sy] = true;
		ArrayList<int[]> queue = new ArrayList<int[]>();
		queue.add(new int[] {sx, sy, sx+1, sy});
		queue.add(new int[] {sx, sy, sx-1, sy});
		queue.add(new int[] {sx, sy, sx, sy+1});
		queue.add(new int[] {sx, sy, sx, sy-1});
		
		int[] points;
		int counter = 0;
		while(queue.size() > 0)
		{
//			if (counter % 10000 == 0) System.out.println(counter + "\t" + queue.size());
			counter++;
			points = queue.remove(0);
			visited[points[2]][points[3]] = true;
			double guess = maxPos[points[0]][points[1]];
			sfr = new SpectraUtil.SpectraFittingRange(guess-hWidth, guess+hWidth, t.v, t.getSpectrum(points[2], points[3]));
			r = sfr.fitToParabola();
			maxPos[points[2]][points[3]] = r.fitParams[0];
			maxCurv[points[2]][points[3]] = r.fitParams[1];
			maxRaw[points[2]][points[3]] = r.x[tryingToGetMinimum ? ArrayOps.maxIndex(r.y) : ArrayOps.maxIndex(r.y)];
			maxOffset[points[2]][points[3]] = r.fitParams[2];
//			System.out.print("" + points[2] + "\t" + points[3] + "\t" + Printer.arrayLnHorizontal(sfr.fitY));
			
			if (points[2] > 0 && !visited[points[2]-1][points[3]]){
				queue.add(new int[] {points[2], points[3], points[2]-1, points[3]}); visited[points[2]-1][points[3]] = true;}
			if (points[3] > 0 && !visited[points[2]][points[3]-1]){
				queue.add(new int[] {points[2], points[3], points[2], points[3]-1}); visited[points[2]][points[3]-1] = true;}
			if (points[2] < t.nx-1 && !visited[points[2]+1][points[3]]){
				queue.add(new int[] {points[2], points[3], points[2]+1, points[3]}); visited[points[2]+1][points[3]] = true;}
			if (points[3] < t.ny-1 && !visited[points[2]][points[3]+1]){
				queue.add(new int[] {points[2], points[3], points[2], points[3]+1}); visited[points[2]][points[3]+1] = true;}
		}
		Layer[] ans = new Layer[4];
		ans[0] = Layer.getFreeLayer(maxPos);
		ans[1] = Layer.getFreeLayer(maxCurv);
		ans[2] = Layer.getFreeLayer(maxRaw);
		ans[3] = Layer.getFreeLayer(maxOffset);
		return ans;
	}
	
	
	public static Layer[] fitAllSpectraToParabola(Topomap t, double center, double hWidth, boolean showGraph)
	{
		SpectraUtil.SpectraFittingRange range = new SpectraUtil.SpectraFittingRange(center-hWidth, center+hWidth, t.v, t.getSpectrum(0, 0));
//		int kOff = ArrayOps.indexOf(t.v, range.y[0], true);
		double[][] peak = new double [t.nx][t.ny], curvature = new double [t.nx][t.ny];
		double[][][] calcData = new double [range.fitY.length][t.nx][t.ny];
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
			{
				range = new SpectraUtil.SpectraFittingRange(center-hWidth, center+hWidth, t.v, t.getSpectrum(i, j));
				FittingResult parab = range.fitToParabola();
				peak[i][j] = parab.fitParams[0];
				curvature[i][j] = parab.fitParams[1];
				for (int k = 0; k < parab.yCalc.length; k++)
					calcData[k][i][j] = parab.yCalc[k];//-parab.y[k];
					
			}
		Layer[] ans = new Layer[range.fitY.length+2];
		ans[0] = Layer.getFreeLayer(peak);
		ans[1] = Layer.getFreeLayer(curvature);
		for (int k = 0; k < calcData.length; k++){
			ans[k+2] = Layer.getFreeLayer(calcData[k]);
			ans[k+2].v = range.fitX[k];
		}
		
		if (showGraph){
		double[][] gx = new double [2*t.nx*t.ny][];
		double[][] gy = new double [2*t.nx*t.ny][];
		double dh = 0.005;
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
			{
				int n = 2*(i*t.ny+j);
				gx[n] = t.v;
				gy[n] = ArrayOps.add(t.getSpectrum(i, j), n*dh);
				gx[n+1] = range.fitX;
				gy[n+1] = new double [calcData.length];
				for (int k = 0; k < calcData.length; k++)
					gy[n+1][k] = calcData[k][i][j] + n*dh;
				}
		GraphDrawerCart g = GraphDrawerCart.getNew(gx, gy);

		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
			{
				int n = 2*(i*t.ny+j);
				g.setColor(Color.RED, n);
				g.setColor(Color.BLUE, n+1);				
			}
		g.showWindow();
		}
		return ans;
	}
	/**
	 * This uses an initial guess as the starting point then smooths it and uses that as the next guess.
	 * @param t
	 * @param center
	 * @param hWidth
	 * @return
	 */
	public static Layer[] fitAllSpectraToParabola_Double(Topomap t, double center, double hWidth, double reductionRatio)
	{
		SpectraUtil.SpectraFittingRange range = new SpectraUtil.SpectraFittingRange(center-hWidth, center+hWidth, t.v, t.getSpectrum(0, 0));
//		int kOff = ArrayOps.indexOf(t.v, range.y[0], true);
		double[][] peak = new double [t.nx][t.ny], curvature = new double [t.nx][t.ny];
		double[][][] calcData = new double [range.fitY.length][t.nx][t.ny];
		double[][] guesses = fitAllSpectraToParabola(t, center, hWidth, false)[0].data;
		FieldOps.cutOffExtremes(guesses, 0.005, 0.995);
		guesses = FieldOps.gaussSmooth(guesses, 4);
		double[][] gx = new double [2*t.nx*t.ny][];
		double[][] gy = new double [2*t.nx*t.ny][];
		double dh = 0.005;
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
			{
				range = new SpectraUtil.SpectraFittingRange(guesses[i][j]-hWidth/reductionRatio, guesses[i][j]+hWidth/reductionRatio, t.v, t.getSpectrum(i, j));
				FittingResult parab = range.fitToParabola();
				peak[i][j] = parab.fitParams[0];
				curvature[i][j] = parab.fitParams[1];
				int n = 2*(i*t.ny+j);
				gx[n] = t.v;
				gy[n] = ArrayOps.add(t.getSpectrum(i, j), n*dh);
				gx[n+1] = range.fitX;
				gy[n+1] = new double [parab.yCalc.length];
				for (int k = 0; k < parab.yCalc.length; k++){
					gy[n+1][k] = parab.yCalc[k] + n*dh;
//					calcData[k][i][j] = parab.yCalc[k] - parab.y[k];
				}
				
			}
		Layer[] ans = new Layer[2];
		ans[0] = Layer.getFreeLayer(peak);
		ans[1] = Layer.getFreeLayer(curvature);
//		for (int k = 0; k < calcData.length; k++){
//			ans[k+2] = Layer.getFreeLayer(calcData[k]);
//			ans[k+2].v = range.fitX[k];
//		}

		GraphDrawerCart g = GraphDrawerCart.getNew(gx, gy);

		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
			{
				int n = 2*(i*t.ny+j);
				g.setColor(Color.RED, n);
				g.setColor(Color.BLUE, n+1);				
			}
		g.showWindow();
		return ans;
	}
	
	public static void saveCropped(Topomap t)
	{
		int xi; int xf; int yi; int yf;
		String input = JOptionPane.showInputDialog("Enter the minimum and maximum of x, comma separated.", "" + 0 + "," + t.nx);
		String[] token = input.split(",");
		xi = Integer.parseInt(token[0]); xf = Integer.parseInt(token[1]);
		input = JOptionPane.showInputDialog("Enter the minimum and maximum of y, comma separated.", "" + 0 + "," + t.ny);
		token  = input.split(",");
		yi = Integer.parseInt(token[0]); yf = Integer.parseInt(token[1]);
		
		Topomap.writeBIN(Topomap.getCropped(t, xi, xf, yi, yf), FileOps.selectSave(fc).toString());
	}
	
	public static void writeNanonisNumberAveragesCommas(Topomap t, double factor, boolean flip, int interval)
	{
		double[] mean = new double [t.nlayers];
		for (int i = 0; i < t.nlayers; i++)
			mean[i] = FieldOps.mean(t.data[i]);
		
		String it = "";
		int j;
		for (int i = 0; i < t.nlayers; i+= interval){
			j = flip ? t.nlayers-i-1 : i;
			String tok = NanonisNumber.parseDouble(mean[j]*factor).toString();
			it += tok + ",";
//			System.out.println(t.v[i] +"\t"+ tok);
			System.out.println(tok);
		}
		
		System.out.println(it);
		it = "";
		for (int i = 0; i < t.nlayers; i+= interval){
			j = flip ? t.nlayers-i-1 : i;
			String tok = NanonisNumber.parseDouble(mean[j]*factor).toString();
			it += tok + ",";
			System.out.println(t.v[j] +"\t"+ tok);
//			System.out.println(tok);
		}
		
		System.out.println(it);
	}
	public static void flipCertainLayersY(Topomap t, int[] layers)
	{
		for (int i = 0; i < layers.length; i++)
			ArrayOps.flipY(t.data[layers[i]]);
	}
	public static void writeNanonisNumberAveragesCommas(Topomap t, double factor, boolean flip, double[] v)
	{
		double[] mean = new double [t.nlayers];
		for (int i = 0; i < t.nlayers; i++)
			mean[i] = FieldOps.mean(t.data[i]);
		
		Functions.CubicSplineWrapper f = new Functions.CubicSplineWrapper (t.v, mean);
		
		String it = "";
		for (int i = 0; i < v.length; i++){
			String tok = NanonisNumber.parseDouble(Math.abs(f.of(v[i]) - f.of(0))*factor).toString();
			it += tok + ",";
//			System.out.println(t.v[i] +"\t"+ tok);
			System.out.println(tok);
		}
		
		System.out.println(it);
		it = "";
		for (int i = 0; i < v.length; i++){
			String tok = NanonisNumber.parseDouble(f.of(v[i]) - f.of(0)*factor).toString();
			it += tok + ",";
			System.out.println(v[i] +"\t"+ tok);
//			System.out.println(tok);
		}
		for (int i = 0; i < v.length; i++){
			String tok = NanonisNumber.parseDouble(v[i]).toString();
			it += tok + ",";
			System.out.println(tok);
//			System.out.println(tok);
		}
		
		System.out.println(it);
	}
	
	public static void rotate180(Topomap t)
	{
		ArrayOps.flip(t.x);
		ArrayOps.flip(t.y);
		for (int i = 0; i < t.nlayers; i++)
		{
			ArrayOps.flipX(t.data[i]);
			ArrayOps.flipY(t.data[i]);
		}
	}
	public static void flipEnergy(Topomap t)
	{
		ArrayOps.flip(t.v);
		double[][][] temp = new double [t.nlayers][][];
		for (int i = 0; i < t.nlayers; i++)
		{
			temp[i] = t.data[t.nlayers-1-i];
		}
		t.data = temp;
	}
	public static void flipATopomap(JFileChooser fc)
	{
		File open = FileOps.selectOpen(fc);
		Topomap t = Topomap.readBIN(open.toString());
		flipEnergy(t);
		Topomap.writeBIN(t, open.toString());
	}
	
	public static void wrapIterativelySuppress(Topomap t)
	{
		double li = Double.parseDouble(JOptionPane.showInputDialog("Enter the real-space smoothing length in pixels."));
		double lk = Double.parseDouble(JOptionPane.showInputDialog("Enter the energy smoothing length in pixels."));
		boolean saveEachIter = JOptionPane.showConfirmDialog(null, "Save each iteration?") == JOptionPane.YES_OPTION;
		boolean saveSheet = JOptionPane.showConfirmDialog(null, "Save each sheet?") == JOptionPane.YES_OPTION;
		boolean saveDiff = JOptionPane.showConfirmDialog(null, "Save each collection of differences?") == JOptionPane.YES_OPTION;
		double minperc = Double.parseDouble(JOptionPane.showInputDialog("Minimum distrubition cutoff (0 to 1)?"));
		double maxperc = Double.parseDouble(JOptionPane.showInputDialog("Maximum distrubition cutoff (0 to 1)?"));
		String basepath = FileOps.selectSave(fc).toString();
		iterativelySuppressBadPixels(t, minperc, maxperc, li, lk, saveDiff, saveSheet, saveEachIter, basepath, Integer.parseInt(JOptionPane.showInputDialog("How many iterations?")));
	}
	
	/**
	 * For use with the 3D gaussian smoothing methods.
	 * @param t
	 * @param lower
	 * @param upper
	 * @param li
	 * @param lk
	 */
	public static void iterativelySuppressBadPixels(Topomap t, double lower, double upper, double li, double lk, boolean saveEachDiff, boolean saveEachSheet, boolean saveEachIteration, String basepath, int niter)
	{
		double[][][] smooth = null;
		File sheet = null, differences = null, iteration = null;
		for (int i = 0; i < niter; i++)
		{ 	System.out.println(i);
			if (saveEachSheet) sheet = new File(basepath + "sheet_" + i + ".bin");
			if (saveEachDiff) differences = new File(basepath + "diff_" + i + ".bin");
			if (saveEachIteration) iteration = new File(basepath + "iteration_" + i + ".bin");
			smooth = FieldOps.getGaussianSmoothing3D(t.data, lk, li, li);
			FieldOps.minusEquals(t.data, smooth);
			double[][][] diff = FieldOps.cutOffExtremes3D(t.data, lower, upper, false, saveEachDiff);
			FieldOps.plusEquals(t.data, smooth);
			
			if (saveEachSheet)
				Topomap.writeBIN(Topomap.newTopomap(t, smooth), sheet.toString());
			if (saveEachDiff)
				Topomap.writeBIN(Topomap.newTopomap(t, diff), differences.toString());
			if (saveEachIteration)
				Topomap.writeBIN(t, iteration.toString());
		}
		Topomap.writeBIN(t, basepath + "Final iteration.bin");
		Topomap.writeBIN(Topomap.newTopomap(t, smooth), basepath + "Final sheet.bin");
	}
	public static int[] getHistogram(Topomap t, int nbins)
	{
		double[] allData = FieldOps.getArray(t.data);
		int[] ans = ArrayOps.getHistogram(allData, ArrayOps.generateArrayNotInclUpper(FieldOps.min(t.data), ArrayOps.max(t.data), nbins));
		return ans;
	}
	
	public static Topomap getLocalNematicity(Topomap t, AtomicCoordinatesSet latt)
	{
		double[][][] data = new double [t.nlayers][t.nx][t.ny];
		for (int i = 0; i < t.nlayers; i++)
			data[i] = LayerUtil.getDerivativeNematicity(t.getLayer(i), latt).data;
		return Topomap.newTopomap(t, data);
	}
	
	public static Layer getSimplestGapMap(Topomap t, double cutoff)
	{
		double[][] data = new double [t.nx][t.ny];
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
				for (int k = 0; k < t.nlayers; k++)
					if (t.data[k][i][j] < cutoff) data[i][j]++;
		
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
				data[i][j] *= t.v[1]-t.v[0];
		
		Layer ans = Layer.newLayer(t, data);
		ans.v = cutoff;
		return ans;
	}
	public static Layer[] getSpectrumMinimum(Topomap t)
	{
		double[][] min = new double [t.nx][t.ny];
		double[][] minV = new double [t.nx][t.ny];
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
			{
				min[i][j] = ArrayOps.min(t.getSpectrum(i, j));
				minV[i][j] = t.v[ArrayOps.minIndex(t.getSpectrum(i,j))];
			}
		Layer minl = Layer.newLayer(t, min);

		Layer minv = Layer.newLayer(t, minV);
		return new Layer[] {minl, minv};
	}
	public static LineOfSpectra[] splitTopomapIntoLinesOfSpectra(boolean horizontal, Topomap t, Layer z)
	{
		//if !horizontal, vertical. 
		LineOfSpectra[] ls;
		
		double[][] data;
		double[] x; double[] y; double [] zz;
		if (horizontal) //then each line cut is horizontal, and they are spread out vertically.
		{
			ls = new LineOfSpectra[t.ny];
			for (int i = 0; i < t.ny; i++)
			{
				y = ArrayOps.clone(t.y[i], t.nx);
				x = t.x;
				data = new double[t.nx][t.nlayers];
				zz = new double[t.nx];
				for (int k = 0; k < t.nx; k++)
				{
					zz[k] = z.data[k][i]*(0.346/0.189);
					data[k] = t.getSpectrum(k, i);
				}
				ls[i] = new LineOfSpectra(data, t.v, x, y, zz);
			}
			return ls;
		}
		else
		{
			ls = new LineOfSpectra[t.nx];
			for (int i = 0; i < t.nx; i++)
			{
				y = t.y;
				x = ArrayOps.clone(t.x[i], t.ny);
				data = new double[t.ny][t.nlayers];
				zz = new double[t.ny];
				for (int k = 0; k < t.ny; k++)
				{
					zz[k] = z.data[i][k]*(0.346/0.189);
					data[k] = t.getSpectrum(i, k);
				}
				ls[i] = new LineOfSpectra(data, t.v, x, y, zz);
			}
			return ls;
		}
	}

	public static void splitTopomapByBinsOfLayer(Topomap t, String dir, String name, int nbins, JFileChooser fc)
	{
		Layer dist = Layer.open(fc);
		String toponame = fc.getSelectedFile().getName().substring(0, fc.getSelectedFile().getName().length()-4);
		int[][] bins = FieldOps.getPercentileBinsForField(dist.data, nbins);
		int min = FieldOps.min(bins);
		int max = FieldOps.max(bins);
		int n = (max-min)+1;
		Topomap split;
		String outdir = dir + name + "split_" + toponame + "_" + n + "bins\\";
		if (!new File (outdir).exists()) new File (outdir).mkdir();
		double[][] spec = new double [n][t.nlayers];
		boolean writeMaps = JOptionPane.showConfirmDialog(null, "Write the separate topomaps?", "Hard disk usage", JOptionPane.YES_NO_OPTION) == JOptionPane.YES_OPTION;
		for (int i = 0; i < n; i++){
			split = splitTopomapBinned(t, bins, i);
			if (writeMaps) 
				Topomap.writeBIN(split, outdir + name + (i+1) + "_of_" + n + ".bin");
			spec[i] = split.getAverageSpectrum();
		}
		String header = "bin";
		for (int i = 0; i < nbins; i++)
			header += "\t" + (i+1);
		header += "\r\n";
		for (int i = 0; i < t.nlayers; i++){
			header += t.v[i];
			for (int j = 0; j < n; j++)
				header += "\t" + spec[j][i];
			header += "\r\n";
		}
		ColumnIO.writeString(header, outdir + name + "table.txt");
		double[][] binsd = ArrayOps.toDouble(bins);
		Layer.writeBIN(Layer.newLayer(dist, binsd), outdir + name + "bins.bin");
		GraphDrawerCart.plotGraphs(t.v, spec);

	}
	public static void splitTopomapByBinsOfLayerFast(Topomap t, String dir, String name, int nbins, JFileChooser fc)
	{
		Layer dist = Layer.open(fc);
		String toponame = fc.getSelectedFile().getName().substring(0, fc.getSelectedFile().getName().length()-4);
		int[][] bins = FieldOps.getPercentileBinsForField(dist.data, nbins);
		int min = FieldOps.min(bins);
		FieldOps.minusEquals(bins, min); //the bin array will now start at zero, the same as the array of sepctra
		
		int max = FieldOps.max(bins);
		int[] npixPerBin = new int [nbins];
		String outdir = dir + name + "split_" + toponame + "_" + nbins + "bins\\";
		if (!new File (outdir).exists()) new File (outdir).mkdir();
		double[][] spec = new double [nbins][t.nlayers];
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
				for (int k = 0; k < t.nlayers; k++)
				{
					spec[bins[i][j]][k] += t.data[k][i][j];
					if (k == 0) npixPerBin[bins[i][j]]++;
				}
		
		for (int i = 0; i < nbins; i++)
		{
			for (int k = 0; k < t.nlayers; k++)
				spec[i][k] /= npixPerBin[i];
			
			System.out.println("" + i + "\t" + npixPerBin[i]);
		}
		String[] lines = new String[t.nlayers+1];
		lines[0] = "bin";
		for (int i = 0; i < nbins; i++)
			lines[0] += "\t" + (i+1);
		for (int i = 0; i < t.nlayers; i++){
			lines[i+1] = "" + t.v[i];
			for (int j = 0; j < nbins; j++)
				lines[i+1] += "\t" + spec[j][i];
		}
		ColumnIO.writeLines(lines, outdir + name + "table.txt");
		double[][] binsd = ArrayOps.toDouble(bins);
		Layer.writeBIN(Layer.newLayer(dist, binsd), outdir + name + "bins.bin");
		GraphDrawerCart.plotGraphs(t.v, spec);

	}
	public static Topomap splitTopomapByBinsOfLayerFast_autoHist(Topomap t, String dir, String name, int nbins, JFileChooser fc)
	{
		Layer dist = Layer.open(fc);
		String toponame = fc.getSelectedFile().getName().substring(0, fc.getSelectedFile().getName().length()-4);
		int[][] bins = FieldOps.getPercentileBinsForField(dist.data, nbins);
		int min = FieldOps.min(bins);
		FieldOps.minusEquals(bins, min); //the bin array will now start at zero, the same as the array of sepctra
		
		int max = FieldOps.max(bins);
		int[] npixPerBin = new int [nbins];
		String outdir = dir + "binned\\" + name + "split_" + toponame + "_" + nbins + "bins\\";
		if (!new File (outdir).exists()) new File (outdir).mkdirs();
		double[][] spec = new double [nbins][t.nlayers];
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
				for (int k = 0; k < t.nlayers; k++)
				{
					spec[bins[i][j]][k] += t.data[k][i][j];
					if (k == 0) npixPerBin[bins[i][j]]++;
				}
				
		
		for (int i = 0; i < nbins; i++)
		{
			for (int k = 0; k < t.nlayers; k++)
				spec[i][k] /= npixPerBin[i];
			
			System.out.println("" + i + "\t" + npixPerBin[i]);
		}
		String[] lines = new String[t.nlayers+1];
		lines[0] = "bin";
		for (int i = 0; i < nbins; i++)
			lines[0] += "\t" + (i+1);
		for (int i = 0; i < t.nlayers; i++){
			lines[i+1] = "" + t.v[i];
			for (int j = 0; j < nbins; j++)
				lines[i+1] += "\t" + spec[j][i];
		}
		ColumnIO.writeLines(lines, outdir + name + "table.txt");
		double[][] binsd = ArrayOps.toDouble(bins);
		Layer.writeBIN(Layer.newLayer(dist, binsd), outdir + name + "bins.bin");
		GraphDrawerCart.plotGraphs(t.v, spec);
		Topomap ans = t.nlayers > 1 ? getSpectralDistributionBasicBinned(1, t.nlayers, t, bins) : null;
		if (t.nlayers > 1 && JOptionPane.showConfirmDialog(null, "Save the spectra distribution map?") == JOptionPane.YES_OPTION)
			Topomap.writeBIN(ans, outdir+name+"spec.bin");
		return ans;
		
	}
	
	/**
	 * Returns a topomap in which each layer is the selected bin in white, other layers being black.
	 * @param bins
	 * @param t
	 * @return
	 */
	public static Topomap getBinTopomap(int[][] bins, Topomap t)
	{
		int max = FieldOps.max(bins);
		return null;
	}
	
	public static Layer splitTopomapByBinsOfLayer(Topomap t, String dir, String name, boolean[][] use, int nbins, JFileChooser fc)
	{
		Layer dist = Layer.open(fc);
		String toponame = fc.getSelectedFile().getName().substring(0, fc.getSelectedFile().getName().length()-4);
		int[][] bins = FieldOps.getPercentileBinsForField(dist.data, nbins, use);
		int min = FieldOps.min(bins);
		int max = FieldOps.max(bins);
		int n = (max-min);
		Topomap split;
		String outdir = dir + name + "_bool_split_" + toponame + "_" + n + "bins\\";
		if (!new File (outdir).exists()) new File (outdir).mkdir();
		double[][] spec = new double [n][t.nlayers];
		boolean writeMaps = JOptionPane.showConfirmDialog(null, "Write the separate topomaps?", "Hard disk usage", JOptionPane.YES_NO_OPTION) == JOptionPane.YES_OPTION;
		for (int i = 0; i < n; i++){
			split = splitTopomapBinned(t, bins, i+1);
			if (writeMaps) 
				Topomap.writeBIN(split, outdir + name + (i+1) + "_of_" + n + ".bin");
			spec[i] = split.getAverageSpectrum();
		}
		String header = "bin";
		for (int i = 0; i < nbins; i++)
			header += "\t" + (i+1);
		header += "\r\n";
		for (int i = 0; i < t.nlayers; i++){
			header += t.v[i];
			for (int j = 0; j < n; j++)
				header += "\t" + spec[j][i];
			header += "\r\n";
		}
		ColumnIO.writeString(header, outdir + name + "table.txt");
		double[][] binsd = ArrayOps.toDouble(bins);
		Layer ans = Layer.newLayer(dist, binsd);
		Layer.writeBIN(ans, outdir + name + "bins.bin");
		GraphDrawerCart.plotGraphs(t.v, spec);
		return ans;
	}
	public static void splitTopomapSpectraByBinsOfLayer(Topomap t, String dir, String name, int nbins, int mapsize, JFileChooser fc)
	{
		if (!new File(dir).exists())
			new File(dir).mkdirs();
		
		Layer dist = Layer.open(fc);
		String layerName = fc.getSelectedFile().getName().substring(0, fc.getSelectedFile().getName().length()-4);
		String fdir = dir + layerName + "\\";
		if (!new File(fdir).exists())
			new File(fdir).mkdirs();
		int[][] bins = FieldOps.getPercentileBinsForField(dist.data, nbins);
		
		Layer[] layers = new Layer[nbins];

		double[][][] isABin = new double [nbins][t.nx][t.ny];
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
				for (int k = 0; k < nbins; k++)
					isABin[k][i][j] = bins[i][j] == k ? 1 : 0;
		Topomap.writeBIN(new Topomap(isABin, ArrayOps.generateArrayNotInclUpper(0, nbins, nbins), t.x, t.y, null), fdir + name + "pointInEachBin.bin");

		for (int i = 0; i < nbins; i++)
			layers[i] = TopomapUtil.getSpectralDistributionFancy(mapsize, mapsize, t, bins, i);
				
		Topomap ans = Topomap.newTopomap(layers);
		Topomap.writeBIN(ans, fdir + name + "spectraSplit.bin");
		double[][] sum = new double [t.nx][t.ny];
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
				for (int k = 0; k < nbins; k++)
					sum[i][j] += layers[k].data[i][j];
		
		Layer.writeBIN(Layer.newLayer(layers[0], sum), fdir + name + "totalSpectralWeight.bin");
	}
	
	public static double[][] getSpectraAtLatticeAndOff(Topomap t, AtomicCoordinatesSet latt)
	{
		double[][][] latticeSites = AtomicCoordinatesGenerator.getLatticeSites(latt, t.nx, null);
		double[][][] halfLatticeSites = AtomicCoordinatesGenerator.getLatticeSites(latt, t.nx, new double[] {0.5, 0.5});
		int la = latticeSites.length;
		int lb = latticeSites[0].length;
		boolean[][] inArea = new boolean [la][lb];
		int n = 0;
		for (int i = 0; i < la; i++)
			for (int j = 0; j < lb; j++){
				inArea[i][j] = (latticeSites[i][j][0] <= t.nx && latticeSites[i][j][0] > 0 && latticeSites[i][j][1] <= t.ny && latticeSites[i][j][1] > 0);
				if (inArea[i][j]) n++;
			}
		System.out.println(n);
		double[][] spectra = new double [3][t.nlayers];
		for (int k = 0; k < t.nlayers; k++){
			spectra[0][k] = t.v[k];
		
			for (int i = 0; i < la; i++)
				for (int j = 0; j < lb; j++)
				{	
					if (inArea[i][j]){
						spectra[1][k] += t.getValueAt(latticeSites[i][j], k);
						spectra[2][k] += t.getValueAt(halfLatticeSites[i][j], k);
					}
				}
			spectra[1][k] /= n;
			spectra[2][k] /= n;
		}
		return spectra;
	}
	public static PointSpectra getSpectraAtLattice(Topomap t, AtomicCoordinatesSet latt, double[] offset)
	{
		double[][][] latticeSites = AtomicCoordinatesGenerator.getLatticeSites(latt, t.nx, offset);
		int la = latticeSites.length;
		int lb = latticeSites[0].length;
		boolean[][] inArea = new boolean [la][lb];
		int n = 0;
		for (int i = 0; i < la; i++)
			for (int j = 0; j < lb; j++){
				inArea[i][j] = (latticeSites[i][j][0] <= t.nx && latticeSites[i][j][0] > 0 && latticeSites[i][j][1] <= t.ny && latticeSites[i][j][1] > 0);
				if (inArea[i][j]) n++;
			}
		System.out.println(n);
		double[][] spectra = new double [n][t.nlayers];
		double[] x = new double [n];
		double[] y = new double [n];
		int index = 0;
		for (int i = 0; i < la; i++)
			for (int j = 0; j < lb; j++){
				for (int k = 0; k < t.nlayers; k++)
					if (inArea[i][j]){
						spectra[index][k] = t.getValueAt(latticeSites[i][j], k);
						x[index] = latticeSites[i][j][0];
						y[index] = latticeSites[i][j][1];
					}
				if (inArea[i][j]) index++;
			}
		
		return new PointSpectra(spectra, t.v, x, y);
	}
	public static PointSpectra[] getSpectraAtLatticeIfSignalBetween(Topomap t, AtomicCoordinatesSet latt, double[] offset, double[] sBounds, int biasIndex)
	{
		double[][][] latticeSites = AtomicCoordinatesGenerator.getLatticeSites(latt, t.nx, offset);
		int la = latticeSites.length;
		int lb = latticeSites[0].length;
		boolean[][] inArea = new boolean [la][lb];
		int[] nSpecEach = new int [sBounds.length];
		int n = 0;
		boolean grouptag = false;
		for (int i = 0; i < la; i++)
			for (int j = 0; j < lb; j++){
				inArea[i][j] = (latticeSites[i][j][0] <= t.nx && latticeSites[i][j][0] > 0 && latticeSites[i][j][1] <= t.ny && latticeSites[i][j][1] > 0);
				if (inArea[i][j]){
					n++;
				grouptag = false;
				for (int k = 0; k < sBounds.length && !grouptag; k++)
					if (t.getValueAt(latticeSites[i][j], biasIndex) < sBounds[k]){
						nSpecEach[k]++; grouptag = true;
					}
				}
				
			}
		System.out.println(n);
		double[][][] spectra = new double [sBounds.length][][];
		double[][] x = new double [sBounds.length][];
		double[][] y = new double [sBounds.length][];
		for (int i = 0; i < sBounds.length; i++)
		{
			spectra[i] = new double [nSpecEach[i]][t.nlayers];
			x[i] = new double [nSpecEach[i]];
			y[i] = new double [nSpecEach[i]];
		}
		
		int[] index = new int [sBounds.length];
		int currentGroup = 0;
		for (int i = 0; i < la; i++)
			for (int j = 0; j < lb; j++){
				grouptag = false;
				for (int k = 0; k < sBounds.length && !grouptag; k++)
					if (t.getValueAt(latticeSites[i][j], biasIndex) < sBounds[k]){
						currentGroup = k; grouptag = true;
					}
				
				for (int k = 0; k < t.nlayers; k++)
					if (inArea[i][j]){
						spectra[currentGroup][index[currentGroup]][k] = t.getValueAt(latticeSites[i][j], k);
						x[currentGroup][index[currentGroup]] = latticeSites[i][j][0];
						y[currentGroup][index[currentGroup]] = latticeSites[i][j][1];
					}
				if (inArea[i][j]) index[currentGroup]++;
			}
		
		PointSpectra[] p = new PointSpectra[sBounds.length];
		for (int i = 0; i < p.length; i++)
			p[i] = new PointSpectra(spectra[i], t.v, x[i], y[i]);
		
		return p;
	}
	
	public static boolean[][] getSignalBelow(Topomap t, double cutoff, int index)
	{
		boolean[][] ans = new boolean[t.nx][t.ny];
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
				ans[i][j] = t.data[index][i][j] < cutoff;
		return ans;
	}
	
	public static void makeCorrelationMovie()
	{
		JFileChooser fc = new JFileChooser("C:\\data\\");
		Topomap t = Topomap.open(fc);
		double[][] imp = FileOps.openBin(fc);
		double[][][] c = TopomapUtil.correlationMapEachEnergy(t, imp);
		Topomap t2 = Topomap.newTopomap(t, c);
		Topomap.writeBIN(t2, fc);
		MovieMaker.makeMovie(t2, true);
	}
	public static void makeIntensityMovie()
	{
		JFileChooser fc = new JFileChooser("C:\\data\\");
		Topomap t = Topomap.open(fc);
		double[][] imp = FieldOps.normalizeWaveFunction(FileOps.openBin(fc));
		double[][][] c = TopomapUtil.unnormalizedCorrelationMapEachEnergy(t, imp);
		Topomap t2 = Topomap.newTopomap(t, c);
		Topomap.writeBIN(t2, fc);
		MovieMaker.makeMovie(t2, true);
	}
	public static void makeTopomapMovie(boolean singleScale)
	{
		JFileChooser fc = new JFileChooser("C:\\data\\");
		Topomap t = Topomap.open(fc);
		MovieMaker.makeMovie(t, singleScale);
	}
	
	//This now also plots a scatter plot of E vs K, assuming that maximum intensity occurs at the momentum K.
	//It also fits the area around the maximum of the spectrum to a certain function which is specified within.
	public static void doRadialStripCut(String dir, String ending, double width, int npts, Topomap t, String functionName)//, boolean normalizeEachLine)
	{
		String outdir = dir + ending + "_rcut\\";// + (normalizeEachLine ? "norm_layer" : "") + "\\";
		String outdirn = dir + ending + "_rcut_norm_layer\\";
		File outd = new File (outdir);
		if (!outd.exists()) outd.mkdir();
		outd = new File (outdirn);
		if (!outd.exists()) outd.mkdir();
		
		int layerpixheight = 6;
		
		if (t == null)
			t = RHKFileOps.getTopoTopo(dir, ending);
		
		double[][] data = radialFFTStripCutEachEnergy(t, width, npts);
		double[][] fittingCurveData = new double [data.length][data[0].length];
		
		//obtain E vs K.
		double[] k = t.getKValuesReal(data[0].length);
		int bmin = 2, bmax = data[0].length-1;
		double [] kmax = new double [t.nlayers];
		for (int i = 0; i < t.nlayers; i++)
		{
			kmax[i] = k[ArrayOps.maxIndex(data[i], bmin, bmax)];
		}
		int peakmin, peakmax;
		int peaksize = 7;
		double[] peakRegion, peakRegionK;
		double[] bkgroundRegion;
		double backmean;
		int bkmin = bmax-11, bkmax = bmax - 1;
		double pkmean;
		double[][] fitPtable = null;
		double[][] fittingPoints = new double [2*peaksize+1][t.nlayers];
		double[][] fittingPointsK = new double [2*peaksize+1][t.nlayers];
		double[] fitParams;
		for (int i = 0; i < t.nlayers; i++)
		{
			peakmin = FieldOps.round(ArrayOps.maxIndex(data[i], bmin, bmax))-peaksize;
			peakmax = peakmin+2*peaksize+1;
			peakmax = Math.min(peakmax, (t.nx/2)-1);
			peakRegion = new double [peakmax-peakmin];
			peakRegionK = new double [peakmax-peakmin];
			bkgroundRegion = new double [bkmax - bkmin];
			for (int b = 0; b < bkgroundRegion.length; b++)
				bkgroundRegion[b] = data[i][b+bkmin];
			
			backmean = ArrayOps.mean(bkgroundRegion);
			
			for (int j = 0; j < peakRegion.length; j++)
			{
				peakRegion[j] = data[i][j+peakmin];
				peakRegionK[j] = k[j+peakmin];
			}
			pkmean = ArrayOps.mean(peakRegion);
//			reg = new Regression(peakRegionK, peakRegion);
//			reg.linear();
			fitParams = ACM_NonLinearFitter.fitToFunction(peakRegionK, peakRegion, functionName);
			fittingCurveData[i] = ArrayUtil.getExpectedValues(k, fitParams, functionName);
			
			if (fitPtable == null) fitPtable = new double [1+fitParams.length][t.nlayers];
			
			fitPtable[0][i] = t.v[i];
			for (int j = 0; j < fitParams.length; j++)
				fitPtable[j+1][i] = fitParams[j];
			
			
			
			for (int j = 0; j < peakRegion.length; j++)
			{
				fittingPoints[j][i] = peakRegion[j];
				fittingPointsK[j][i] = peakRegionK[j];
			}
			kmax[i] = fitParams[0];
		}
		ColumnIO.writeTwoColumns(kmax, t.v, "EvsK_scatter.txt", outdir);
		ColumnIO.writeTable(fitPtable, outdir + "Lorenztian Ps mean width height.txt");
		//Note re above: the Lorenztian is (Height/pi)*(Width/2)/((x-mean)^2 + (width/2)^2)
		ColumnIO.writeTable(fittingPointsK, outdir + "Fitting region K values.txt");
		ColumnIO.writeTable(fittingPoints, outdir + "Fitting region intensity values.txt");
//		ColumnIO.writeTable(FieldOps.transpose(fittingPoints), outdir + "Fitting region intensity values Transpose.txt");
//		ColumnIO.writeTable(FieldOps.transpose(fittingPointsK), outdir + "Fitting region K values Transpose.txt");
		
		PointSpectra.writeBIN(new PointSpectra(data, k, new double[t.nlayers], t.v), outdir + "Intensity spectra I vs k.bin");
		PointSpectra.writeBIN(new PointSpectra(fittingCurveData, k, new double[t.nlayers], t.v), outdir + "Curve fit I vs k.bin");
		//
		
		
		layerpixheight = (bmax-bmin)/(t.nlayers+1);
		System.out.println(layerpixheight);
		double[][] imageout = FieldOps.truncate(data, 0, data.length, bmin, bmax);
		double[][] tableout = FieldOps.augment(data, t.v, k);
		SRAW.writeInflatedImage(outdir + "Evsk.bmp", imageout, layerpixheight, 1);
		ColumnIO.writeTable(tableout, outdir + "RadialCuts.txt");
		ColumnIO.writeTable(FieldOps.augment(fittingCurveData, t.v, k), outdir + "RadialCut_FitCurves.txt");
		ColumnIO.writeTable(imageout, outdir + "RadialCuts_T.txt");
		ColumnIO.writeString("[" + k[bmin] + ", " + k[bmax-1] + "]", outdir + "Kbounds.txt");

		double[][] normalized = FieldOps.normalizeEachRow(imageout, false);
		SRAW.writeInflatedImage(outdirn + "Evsk.bmp", normalized, layerpixheight, 1);
		ColumnIO.writeTable(normalized, outdirn + "RadialCuts_Picture.txt");
//		imageout = FieldOps.truncate(normalized, 0, data.length, 2, data[0].length/2);
		normalized = FieldOps.normalizeEachRow(data, false);
		tableout = FieldOps.augment(normalized, t.v, t.getKValuesReal(data[0].length));
		ColumnIO.writeTable(tableout, outdirn + "RadialCuts.txt");
		
		double[][][] polar = polarCoordFFTEachEnergy(t, width*5, npts/5);
		for (int i = 0; i < polar.length; i++)
		{
			imageout = FieldOps.truncate(polar[i], 0, polar[i].length, bmin, bmax);
			SRAW.writeImage(outdir + "polar_" + MovieMaker.fromInt(i), imageout);
		}
	}
	public static PointSpectra[] doHighSymmetryStripCut(AtomicCoordinatesSet unitBragg, String dir, String ending, double width, Topomap t)//, boolean normalizeEachLine)
	{

		String outdir = dir + ending + "_cut_" + "w"+(int)width+"\\";// + (normalizeEachLine ? "norm_layer" : "") + "\\";
		String outdirn = dir + ending + "_cut_norm_layer_" + "w"+(int)width+"\\";
		File outd = new File (outdir);
		if (!outd.exists()) outd.mkdir();
		outd = new File (outdirn);
		if (!outd.exists()) outd.mkdir();
		
		int layerpixheight = 6;
		
		if (t == null)
			t = RHKFileOps.getTopoTopo(dir, ending);
		
		double braggAngleDeg = Math.toDegrees(unitBragg.getAngleBetweenVectors());
		
		double[][] onData = stripCutEachEnergy(t, width, unitBragg.getA());
		double[] otherVector = Matrix.getProductWith(Matrix.getRotationMatrix(Math.toRadians(braggAngleDeg/2)), unitBragg.getA());
		double[][] offData = stripCutEachEnergy(t, width, otherVector);
	
		//obtain E vs K. For K, we remember that the dataset is the symmetrized FFT so that it's real space is actually k-space.
		double[] k = new double [t.nx/2];
		for (int i = 0; i < t.nx/2; i++)
			k[i] = Math.abs(t.x[i+t.nx/2]);
		int bmin = 2, bmax = onData[0].length-1;
//		ColumnIO.writeTable(FieldOps.transpose(fittingPoints), outdir + "Fitting region intensity values Transpose.txt");
//		ColumnIO.writeTable(FieldOps.transpose(fittingPointsK), outdir + "Fitting region K values Transpose.txt");
		
		PointSpectra on = new PointSpectra(onData, k, new double[t.nlayers], t.v);
		PointSpectra off = new PointSpectra(offData, k, new double[t.nlayers], t.v);
		PointSpectra.writeBIN(on, outdir + "Intensity spectra I vs k on direction.bin");
		PointSpectra.writeBIN(off, outdir + "Intensity spectra I vs k off direction.bin");
		Layer onl = on.toLayer();
		Layer offl = off.toLayer();
		Layer.writeBIN(onl, outdir + "Intensity spectra I vs k on directionlayer.bin");
		PointSpectra.writeBIN(off, outdir + "Intensity spectra I vs k off directionlayer.bin");
		//
		layerpixheight = (bmax-bmin)/(t.nlayers+1);
		System.out.println(layerpixheight);
		double[][] imageouton = FieldOps.truncate(onData, 0, onData.length, bmin, bmax);
		double[][] tableouton = FieldOps.augment(onData, t.v, k);
		double[][] imageoutoff = FieldOps.truncate(offData, 0, onData.length, bmin, bmax);
		double[][] tableoutoff = FieldOps.augment(offData, t.v, k);
		SRAW.writeInflatedImage(outdir + "Evsk_on.bmp", imageouton, layerpixheight, 1);
		ColumnIO.writeTable(tableouton, outdir + "RadialCuts_on.txt");
		SRAW.writeInflatedImage(outdir + "Evsk_off.bmp", imageoutoff, layerpixheight, 1);
		ColumnIO.writeTable(tableoutoff, outdir + "RadialCuts_off.txt");
		ColumnIO.writeString("[" + k[bmin] + ", " + k[bmax-1] + "]", outdir + "Kbounds.txt");

		double[][] normalized = FieldOps.normalizeEachRow(imageouton, false);
		SRAW.writeInflatedImage(outdirn + "Evsk_on.bmp", normalized, layerpixheight, 1);
		ColumnIO.writeTable(normalized, outdirn + "RadialCuts_Picture_on.txt");
		normalized = FieldOps.normalizeEachRow(imageoutoff, false);
		SRAW.writeInflatedImage(outdirn + "Evsk_off.bmp", normalized, layerpixheight, 1);
		ColumnIO.writeTable(normalized, outdirn + "RadialCuts_Picture_off.txt");
//		imageout = FieldOps.truncate(normalized, 0, data.length, 2, data[0].length/2);
		writeGIFWithDirections(unitBragg.getA(), otherVector, t, outdir + "Movie.gif");
		new SpectraDrawer(off, null);
		return new PointSpectra[] {on, off};
	
	}
	public static PointSpectra getRadialAverageStripCut(double width, Topomap t)//, boolean normalizeEachLine)
	{
		double[][] data = radialFFTStripCutEachEnergy(t, width, 1000);
		double[] k = t.getKValuesReal(data[0].length);
		return new PointSpectra(data, k, new double[t.nlayers],  t.v);
	}
	public static PointSpectra doStripCutAcross_Rt2(AtomicCoordinatesSet latt, String dir, String ending, double width, Topomap t)//, boolean normalizeEachLine)
	{

		String outdir = dir + ending + "_cut_" + "w"+(int)width+"\\";// + (normalizeEachLine ? "norm_layer" : "") + "\\";
		File outd = new File (outdir);
		if (!outd.exists()) outd.mkdir();
		
		int layerpixheight = 6;
		
		if (t == null)
			t = RHKFileOps.getTopoTopo(dir, ending);
		
		double braggAngleDeg = Math.toDegrees(latt.getAngleBetweenVectors());
		
		//Now we must ascertain the positions of the Rt2 peaks:
		AtomicCoordinatesSet rt2 = latt.getRt2Lattice();
		System.out.println(rt2);
		int[][] braggPoints = AtomicCoordinatesSet.generateBragg(rt2, t.nx);
		double[][] data = stripCutEachEnergy(t, width, braggPoints[0], braggPoints[1]);
//		double[] otherVector = Matrix.getProductWith(Matrix.getRotationMatrix(Math.toRadians(braggAngleDeg/2)), unitBragg.getA());
//		double[][] offData = stripCutEachEnergy(t, width, otherVector);
//	
//		//obtain E vs K. For K, we remember that the dataset is the symmetrized FFT so that it's real space is actually k-space.
		double[] k = new double [t.nx/2];
		for (int i = 0; i < t.nx/2; i++)
			k[i] = Math.abs(t.x[i+t.nx/2]);
		int bmin = 2, bmax = data[0].length-1;
		PointSpectra on = new PointSpectra(data, k, new double[t.nlayers], t.v);
		PointSpectra.writeBIN(on, outdir + "Intensity spectra I vs k on direction.bin");
		//
		layerpixheight = (bmax-bmin)/(t.nlayers+1);
		System.out.println(layerpixheight);
		double[][] imageouton = FieldOps.truncate(data, 0, data.length, bmin, bmax);
		double[][] tableouton = FieldOps.augment(data, t.v, k);
		SRAW.writeInflatedImage(outdir + "Evsk_on.bmp", imageouton, layerpixheight, 1);
		ColumnIO.writeTable(tableouton, outdir + "RadialCuts_on.txt");
		ColumnIO.writeString("[" + k[bmin] + ", " + k[bmax-1] + "]", outdir + "Kbounds.txt");

		double[][] normalized = FieldOps.normalizeEachRow(imageouton, false);
//		imageout = FieldOps.truncate(normalized, 0, data.length, 2, data[0].length/2);
		writeGIFWithLineSeg(braggPoints[0], braggPoints[1], t, outdir + "Movie.gif");
		new SpectraDrawer(on, null);
		return null;
	
	}
	public static void writeGIFWithDirections(double[] onDir, double[] offDir, Topomap t, String path)
	{
		//We also write a gif showing the directions as a function of energy:
		BufferedImage[] gifstack = new BufferedImage[t.nlayers];

		java.awt.Graphics g;
		double[] ondR = new double [2], offdR = new double [2];
		ondR[0] = onDir[0]*t.nx/(2*Complex.mag(onDir));
		ondR[1] = onDir[1]*t.nx/(2*Complex.mag(onDir));
		offdR[0] = offDir[0]*t.nx/(2*Complex.mag(offDir));
		offdR[1] = offDir[1]*t.nx/(2*Complex.mag(offDir));
		for (int i = 0; i < t.nlayers; i++)
		{
			gifstack[i] = ImageEditing.getBufferedImage(t.data[i]);
			g = gifstack[i].getGraphics();
			g.setColor(java.awt.Color.BLUE);
			g.drawLine(t.nx/2, t.ny/2, (int)(t.nx/2+ondR[0]), (int)(t.ny/2+ondR[1]));
			g.drawString("on", (int)(t.nx/2+ondR[0])-10, (int)(t.ny/2+ondR[1]));
			g.drawString(NumFormat.voltage(t.v[i]), 5, t.ny-5);
			g.setColor(java.awt.Color.MAGENTA);
			g.drawLine(t.nx/2, t.ny/2, (int)(t.nx/2+offdR[0]), (int)(t.ny/2+offdR[1]));
			g.drawString("off", (int)(t.nx/2+offdR[0])-10, (int)(t.ny/2+offdR[1]));
			
		}
		GifSequenceWriter.writeGifSequence(gifstack, 500, path);
	}
	public static void writeGIFWithLineSeg(int[] begin, int[] end, Topomap t, String path)
	{
		//We also write a gif showing the directions as a function of energy:
		BufferedImage[] gifstack = new BufferedImage[t.nlayers];

		java.awt.Graphics g;
		for (int i = 0; i < t.nlayers; i++)
		{
			gifstack[i] = ImageEditing.getBufferedImage(t.data[i]);
			g = gifstack[i].getGraphics();
			g.setColor(java.awt.Color.BLUE);
			g.drawLine((int)t.nx/2+begin[0], (int)t.ny/2+begin[1], (int)t.nx/2+end[0], (int)t.ny/2+end[1]);
			g.drawString("start", (int)(t.nx/2+begin[0])-10, (int)(t.ny/2+begin[1]));
			g.drawString(NumFormat.voltage(t.v[i]), 5, t.ny-5);
			g.setColor(java.awt.Color.MAGENTA);
		}
		GifSequenceWriter.writeGifSequence(gifstack, 500, path);
	}
	public static void writeGIFWithVoltage(Topomap t, String path, int sizeFactor, int cscale)
	{
		//We also write a gif showing the directions as a function of energy:
		BufferedImage[] gifstack = new BufferedImage[t.nlayers];
		BufferedImage temp;
		java.awt.Graphics g;
		for (int i = 0; i < t.nlayers; i++)
		{
			temp = ImageEditing.getBufferedImage(t.data[i], cscale);
			
			gifstack[i] = ImageEditing.getEnlarged(temp, sizeFactor);	
			g = gifstack[i].getGraphics();
			g.setColor(ColorScales.getUnusedColor(cscale));
			g.drawString(NumFormat.voltage(t.v[i]), 5, t.ny*sizeFactor-5);
		}
		GifSequenceWriter.writeGifSequence(gifstack, 500, path);
	}
	/**
	 * This tries to do a standard fitting to a single dispersing QPI mode. We start by eyeballing k-value of the peak and fitting its neithborhood to Lorentzian + line.
	 * Then we go to adjacent energies and use the value obtained previously as the guess for the next one (or close). 
	 * @param dir 	Output directory
	 * @param ending	
	 * @param t
	 * @param peaksize	Half-size of the region of k-space in which to conduct the fitting
	 * @param firstPeakGuess	Location of the peak in one layer as established by the user
	 * @param upTrue	Whether in the course of going through the energies we should go up or down.
	 * $
	 */
	public static ACM_NonLinearFitter.FittingResult[] doIvsKPeakFittingStandard(String dir, String ending, PointSpectra t, int peaksize, int firstPeakGuess, boolean upTrue, int firstLayer)//, boolean normalizeEachLine)
	{
		String outdir = dir + ending + "LorentzLines_" + peaksize +"_" + "\\";// + (normalizeEachLine ? "norm_layer" : "") + "\\";
		String suboutdir = outdir + "Curves\\";
		File outd = new File (outdir);
		if (!outd.exists()) outd.mkdirs();
		File subd = new File (suboutdir);
		if (!subd.exists()) subd.mkdirs();

		int layerpixheight = 6;
		
		double[][] data = t.data;
		double[][] fittingCurveData = new double [data.length][data[0].length];
		
		//obtain E vs K.
		double[] k = t.v;
		int bmin = firstPeakGuess-peaksize, bmax = firstPeakGuess+peaksize;
		double [] kmax = new double [t.y.length];
		for (int i = 0; i < t.y.length; i++)
		{
			kmax[i] = k[ArrayOps.maxIndex(data[i], bmin, bmax)];
		}
		double[] peakRegion, peakRegionK;
		double[] peakRegionFit;
		double[][] fitPtable = null;
		
		
		ACM_NonLinearFitter.FittingResult[] fit = new FittingResult[t.y.length];
		for (int i = 0; i < t.y.length; i++)
		{
			peakRegion = new double [bmax-bmin];
			peakRegionK = new double [bmax-bmin];
			peakRegionFit = new double [bmax-bmin];
			
			for (int j = 0; j < peakRegion.length; j++)
			{
				peakRegion[j] = data[i][j+bmin];
				peakRegionK[j] = k[j+bmin];
			}
			
			fit[i] = ACM_NonLinearFitter.fitToFunctionFull(peakRegionK, peakRegion, "1PeakLorentzian_Line");
			
			if (fitPtable == null) fitPtable = new double [1+fit[i].fitParams.length][t.y.length];
			
			fitPtable[0][i] = t.y[i];
			for (int j = 0; j < fit[i].fitParams.length; j++)
				fitPtable[j+1][i] = fit[i].fitParams[j];
			
			ColumnIO.writeTable(new double[][] {peakRegionK, peakRegion, peakRegionFit}, suboutdir + "Curve " + MovieMaker.fromInt(i) + " " + t.y[i] + ".txt");

			bmin = (int)(kmax[i]/(t.v[1]-t.v[0]) - peaksize);
			bmax = (int)(kmax[i]/(t.v[1]-t.v[0]) + peaksize);
			kmax[i] = fit[i].fitParams[0];
			
			for (int j = 0; j < k.length; j++)
				fittingCurveData[i][j] = fit[i].getYAt(k[j], "1PeakLorentzian_Line");
		}
//		fittingCurveData = FieldOps.transpose(fittingCurveData);
		ColumnIO.writeTwoColumns(kmax, t.y, "EvsK_scatter.txt", outdir);
		ColumnIO.writeTable(fitPtable, outdir + "Fitting Ps mean width height.txt");
		//Note re above: the Lorenztian is (Height/pi)*(Width/2)/((x-mean)^2 + (width/2)^2)
//		ColumnIO.writeTable(FieldOps.transpose(fittingPoints), outdir + "Fitting region intensity values Transpose.txt");
//		ColumnIO.writeTable(FieldOps.transpose(fittingPointsK), outdir + "Fitting region K values Transpose.txt");
		PointSpectra.writeBIN(new PointSpectra(fittingCurveData, k, new double[t.y.length], t.y), outdir + "Curve fit I vs k.bin");
		//
		
		
		layerpixheight = (bmax-bmin)/(t.y.length+1);
		System.out.println(layerpixheight);
		SRAW.writeInflatedImage(outdir + "Evsk_fitting_curves.bmp", fittingCurveData, layerpixheight, 1);
		
		BufferedImage spectra = ImageEditing.getBufferedImage(t.data);
		int sx, sy;
		//in the image, sx is the energy factor and sy is the k factor
		if (t.y.length < t.v.length)
		{
			sy = 2;
			sx = sy*(t.v.length/t.y.length);
		}
		else
		{
			sx = 2;
			sy = sx*(t.y.length/t.v.length);
		}
		BufferedImage enlargedSpectra = new BufferedImage(spectra.getWidth()*sx, spectra.getHeight()*sy, BufferedImage.TYPE_INT_RGB);
		ImageEditing.enlargeBasicStretch(spectra, enlargedSpectra, sx, sy);
		SRAW.writeImage(outdir + "Original spectra small", spectra);
		
		int x, y;
		for (int i = 0; i < t.y.length; i++)
		{
			x = sx*i + sx/2;
			y = ArrayOps.indexOf(t.v, fit[i].fitParams[0], true)*sy + sy/2;
			if (t.y[i] >= -0.22)
			ImageEditing.drawPlus(enlargedSpectra, x, y, 4, java.awt.Color.BLUE);
		}
		SRAW.writeImage(outdir + "Spectra", enlargedSpectra);
		
		ColumnIO.writeTable(FieldOps.augment(fittingCurveData, t.v, k), outdir + "RadialCut_FitCurves.txt");
		return fit;
	}
	
	//In this method we first go through a fitting as in doRadialCutPeakFitting, but we then use the results of this fitting to set the center and width of the fitting region
	//for the second iteration of peak fitting, which is expected to be more accurate. The second peak fitting function is user-defined with the first being "Lorentzian" for now.
	public static void doRadialCutPeakTwoPass(String dir, String ending, PointSpectra t, String fname2)//, boolean normalizeEachLine)
	{
		String fname1 = "Lorentzian";
		String outdir1 = dir + ending + "peak2p1" + fname1 + "\\";// + (normalizeEachLine ? "norm_layer" : "") + "\\";
		String outdir2 = dir + ending + "peak2p2" + fname2 + "\\";// + (normalizeEachLine ? "norm_layer" : "") + "\\";
		String outdcurve = outdir2 + "curves\\";
		File outd = new File (outdir1);
		if (!outd.exists()) outd.mkdirs();
		outd = new File (outdir2);
		if (!outd.exists()) outd.mkdirs();
		outd = new File (outdcurve);
		if (!outd.exists()) outd.mkdirs();
		
		int layerpixheight = 6;
		
		double[][] data = t.data;
		double[][] fittingCurveData = new double [data.length][data[0].length];
		
		//obtain E vs K.
		double[] k = t.v;
		int bmin = 2, bmax = data[0].length-1;
		double [] kmax = new double [t.y.length];
		for (int i = 0; i < t.y.length; i++)
		{
			kmax[i] = k[ArrayOps.maxIndex(data[i], bmin, bmax)];
		}
		int peakmin, peakmax;
		int peaksize = 7;
		double[] peakRegion, peakRegionK;
		double[] bkgroundRegion;
		double backmean;
		int bkmin = bmax-11, bkmax = bmax - 1;
		double pkmean;
		double[][] fitPtable = null;
		double[][] fittingPoints = new double [2*peaksize+1][t.y.length];
		double[][] fittingPointsK = new double [2*peaksize+1][t.y.length];
		double[][] fitCurveSub;
		double[] fitParams = null;
		for (int i = 0; i < t.y.length; i++)
		{
			peakmin = FieldOps.round(ArrayOps.maxIndex(data[i], bmin, bmax))-peaksize;
			peakmax = peakmin+2*peaksize+1;
			peakRegion = new double [peakmax-peakmin];
			peakRegionK = new double [peakmax-peakmin];
			bkgroundRegion = new double [bkmax - bkmin];
			for (int b = 0; b < bkgroundRegion.length; b++)
				bkgroundRegion[b] = data[i][b+bkmin];
			
			backmean = ArrayOps.mean(bkgroundRegion);
			
			for (int j = 0; j < peakRegion.length; j++)
			{
				peakRegion[j] = data[i][j+peakmin];
				peakRegionK[j] = k[j+peakmin];
			}
			pkmean = ArrayOps.mean(peakRegion);
//			reg = new Regression(peakRegionK, peakRegion);
//			reg.linear();
			fitParams = ArrayUtil.fitToFunction(peakRegionK, peakRegion, fname1);
			fittingCurveData[i] = ArrayUtil.getExpectedValues(k, fitParams, fname1);
			if (fitPtable == null) fitPtable = new double [1+fitParams.length][t.y.length];
			
			fitPtable[0][i] = t.y[i];
			for (int j = 0; j < fitParams.length; j++)
				fitPtable[j+1][i] = fitParams[j];
				
			for (int j = 0; j < peakRegion.length; j++)
			{
				fittingPoints[j][i] = peakRegion[j];
				fittingPointsK[j][i] = peakRegionK[j];
			}
			kmax[i] = fitParams[0];
		}
		ColumnIO.writeTwoColumns(kmax, t.v, "EvsK_scatter.txt", outdir1);
		ColumnIO.writeTable(fitPtable, outdir1 + "Fitting Ps mean width height.txt");
		//Note re above: the Lorenztian is (Height/pi)*(Width/2)/((x-mean)^2 + (width/2)^2)
		ColumnIO.writeTable(fittingPointsK, outdir1 + "Fitting region K values.txt");
		ColumnIO.writeTable(fittingPoints, outdir1 + "Fitting region intensity values.txt");
//		ColumnIO.writeTable(FieldOps.transpose(fittingPoints), outdir + "Fitting region intensity values Transpose.txt");
//		ColumnIO.writeTable(FieldOps.transpose(fittingPointsK), outdir + "Fitting region K values Transpose.txt");
		PointSpectra.writeBIN(new PointSpectra(fittingCurveData, k, new double[t.y.length], t.y), outdir1 + "Curve fit I vs k.bin");
		//
		layerpixheight = (bmax-bmin)/(t.y.length+1);
		System.out.println(layerpixheight);
		double[][] imageout = FieldOps.truncate(fittingCurveData, 0, data.length, bmin, bmax);
		SRAW.writeInflatedImage(outdir1 + "Evsk_fitting_curves.bmp", imageout, layerpixheight, 1);
		ColumnIO.writeTable(FieldOps.augment(fittingCurveData, t.y, k), outdir1 + "RadialCut_FitCurves.txt");
		//SECOND ITERATION:
		double dk = k[1] - k[0];
		double[][] ptablecopy = FieldOps.copy(fitPtable);
		fitPtable = new double [ACM_CustomFunctions.getNew(fname2).plist.length+1][t.x.length];
		for (int i = 0; i < t.y.length; i++)
		{
			peaksize = (int)(((ptablecopy[3][i]/dk)+1)*2.0/3);
			System.out.println(peaksize);
			peakmin = FieldOps.round(ptablecopy[1][i]/dk)-peaksize;
			peakmax = peakmin+2*peaksize+1;
			peakRegion = new double [peakmax-peakmin];
			peakRegionK = new double [peakmax-peakmin];
			bkgroundRegion = new double [bkmax - bkmin];
			for (int b = 0; b < bkgroundRegion.length; b++)
				bkgroundRegion[b] = data[i][b+bkmin];
			
			backmean = ArrayOps.mean(bkgroundRegion);
			
			for (int j = 0; j < peakRegion.length; j++)
			{
				peakRegion[j] = data[i][j+peakmin];
				peakRegionK[j] = k[j+peakmin];
			}
			pkmean = ArrayOps.mean(peakRegion);
//			reg = new Regression(peakRegionK, peakRegion);
//			reg.linear();
			fitParams = ArrayUtil.fitToFunction(peakRegionK, peakRegion, fname2);
			fittingCurveData[i] = ArrayUtil.getExpectedValues(k, fitParams, fname2);
			//write the curve subsets:
			fitCurveSub = new double [3][peakRegionK.length];
			for (int j = 0; j < peakRegionK.length; j++)
			{
				fitCurveSub[0][j] = peakRegionK[j];
				fitCurveSub[1][j] = peakRegion[j];
				fitCurveSub[2] = ArrayUtil.getExpectedValues(peakRegionK, fitParams, fname2);
			}
			ColumnIO.writeTable(fitCurveSub, outdcurve + "Fitting Curve " + i + " " + t.y[i] + ".txt");
			//
			fitPtable[0][i] = t.y[i];
			for (int j = 0; j < fitParams.length; j++)
				fitPtable[j+1][i] = fitParams[j];
				
//			for (int j = 0; j < peakRegion.length; j++)
//			{
//				fittingPoints[j][i] = peakRegion[j];
//				fittingPointsK[j][i] = peakRegionK[j];
//			}
			kmax[i] = fitParams[0];
		}
		ColumnIO.writeTwoColumns(kmax, t.y, "EvsK_scatter.txt", outdir2);
		ColumnIO.writeTable(fitPtable, outdir2 + "Fitting Ps mean width height.txt");
		//Note re above: the Lorenztian is (Height/pi)*(Width/2)/((x-mean)^2 + (width/2)^2)
//		ColumnIO.writeTable(fittingPointsK, outdir1 + "Fitting region K values.txt");
//		ColumnIO.writeTable(fittingPoints, outdir1 + "Fitting region intensity values.txt");
//		ColumnIO.writeTable(FieldOps.transpose(fittingPoints), outdir + "Fitting region intensity values Transpose.txt");
//		ColumnIO.writeTable(FieldOps.transpose(fittingPointsK), outdir + "Fitting region K values Transpose.txt");
		PointSpectra.writeBIN(new PointSpectra(fittingCurveData, k, new double[t.y.length], t.y), outdir2 + "Curve fit I vs k.bin");
		//
		layerpixheight = (bmax-bmin)/(t.y.length+1);
		System.out.println(layerpixheight);
		imageout = FieldOps.truncate(fittingCurveData, 0, data.length, bmin, bmax);
		SRAW.writeInflatedImage(outdir2 + "Evsk_fitting_curves.bmp", imageout, layerpixheight, 1);
		ColumnIO.writeTable(FieldOps.augment(fittingCurveData, t.y, k), outdir2 + "RadialCut_FitCurves.txt");
	}
	
	/**
	 * This method does the same thing as RadialCutPeakTwoPass but does not write any files and instead returns a boolean array indicating the fourier filter
	 * which would include the fitting region of the 2nd fitting pass.
	 * @param dir
	 * @param ending
	 * @param t
	 * @param fname2
	 * @return
	 */
	public static boolean[][][] doRadialCutPeakTwoPassFFilter(PointSpectra t, String fname2)//, boolean normalizeEachLine)
	{
		String fname1 = "Lorentzian";
			
		int layerpixheight = 6;
		
		double[][] data = t.data;
		double[][] fittingCurveData = new double [data.length][data[0].length];
		
		//obtain E vs K.
		double[] k = t.v;
		int bmin = 2, bmax = data[0].length-1;
		double [] kmax = new double [t.y.length];
		for (int i = 0; i < t.y.length; i++)
		{
			kmax[i] = k[ArrayOps.maxIndex(data[i], bmin, bmax)];
		}
		int peakmin, peakmax;
		int peaksize = 7;
		double[] peakRegion, peakRegionK;
		double[] bkgroundRegion;
		double backmean;
		int bkmin = bmax-11, bkmax = bmax - 1;
		double pkmean;
		double[][] fitPtable = null;
		double[][] fittingPoints = new double [2*peaksize+1][t.y.length];
		double[][] fittingPointsK = new double [2*peaksize+1][t.y.length];
		double[][] fitCurveSub;
		double[] fitParams = null;
		for (int i = 0; i < t.y.length; i++)
		{
			peakmin = FieldOps.round(ArrayOps.maxIndex(data[i], bmin, bmax))-peaksize;
			peakmax = peakmin+2*peaksize+1;
			peakRegion = new double [peakmax-peakmin];
			peakRegionK = new double [peakmax-peakmin];
			bkgroundRegion = new double [bkmax - bkmin];
			for (int b = 0; b < bkgroundRegion.length; b++)
				bkgroundRegion[b] = data[i][b+bkmin];
			
			backmean = ArrayOps.mean(bkgroundRegion);
			
			for (int j = 0; j < peakRegion.length; j++)
			{
				peakRegion[j] = data[i][j+peakmin];
				peakRegionK[j] = k[j+peakmin];
			}
			pkmean = ArrayOps.mean(peakRegion);
//			reg = new Regression(peakRegionK, peakRegion);
//			reg.linear();
			fitParams = ArrayUtil.fitToFunction(peakRegionK, peakRegion, fname1);
			fittingCurveData[i] = ArrayUtil.getExpectedValues(k, fitParams, fname1);
			if (fitPtable == null) fitPtable = new double [1+fitParams.length][t.y.length];
			
			fitPtable[0][i] = t.y[i];
			for (int j = 0; j < fitParams.length; j++)
				fitPtable[j+1][i] = fitParams[j];
				
			for (int j = 0; j < peakRegion.length; j++)
			{
				fittingPoints[j][i] = peakRegion[j];
				fittingPointsK[j][i] = peakRegionK[j];
			}
			kmax[i] = fitParams[0];
		}
		//
		layerpixheight = (bmax-bmin)/(t.y.length+1);
		System.out.println(layerpixheight);
		//SECOND ITERATION:
		double dk = k[1] - k[0];
		double[][] ptablecopy = FieldOps.copy(fitPtable);
		fitPtable = new double [ACM_CustomFunctions.getNew(fname2).plist.length+1][t.x.length];
		boolean[][][] filter = new boolean [t.nspec][t.nlayers][t.nlayers];
		for (int i = 0; i < t.y.length; i++)
		{
			peaksize = (int)(((ptablecopy[3][i]/dk)+1)*2.0/3);
			System.out.println(peaksize);
			peakmin = FieldOps.round(ptablecopy[1][i]/dk)-peaksize;
			peakmax = peakmin+2*peaksize+1;
			peakRegion = new double [peakmax-peakmin];
			peakRegionK = new double [peakmax-peakmin];
			bkgroundRegion = new double [bkmax - bkmin];
			for (int b = 0; b < bkgroundRegion.length; b++)
				bkgroundRegion[b] = data[i][b+bkmin];
			
			backmean = ArrayOps.mean(bkgroundRegion);
			
			for (int j = 0; j < peakRegion.length; j++)
			{
				peakRegion[j] = data[i][j+peakmin];
				peakRegionK[j] = k[j+peakmin];
			}
			pkmean = ArrayOps.mean(peakRegion);
//			reg = new Regression(peakRegionK, peakRegion);
//			reg.linear();
			fitParams = ArrayUtil.fitToFunction(peakRegionK, peakRegion, fname2);
			fittingCurveData[i] = ArrayUtil.getExpectedValues(k, fitParams, fname2);
			//write the curve subsets:
			fitCurveSub = new double [3][peakRegionK.length];
			for (int j = 0; j < peakRegionK.length; j++)
			{
				fitCurveSub[0][j] = peakRegionK[j];
				fitCurveSub[1][j] = peakRegion[j];
				fitCurveSub[2] = ArrayUtil.getExpectedValues(peakRegionK, fitParams, fname2);
			}
			//
			fitPtable[0][i] = t.y[i];
			for (int j = 0; j < fitParams.length; j++)
				fitPtable[j+1][i] = fitParams[j];
				
//			for (int j = 0; j < peakRegion.length; j++)
//			{
//				fittingPoints[j][i] = peakRegion[j];
//				fittingPointsK[j][i] = peakRegionK[j];
//			}
			kmax[i] = fitParams[0];
			if (fname2.equalsIgnoreCase("Lorentzian"))
			{
				peaksize = (int)(((fitPtable[3][i]/dk)+1)*2.0/3);
				peakmin = FieldOps.round(fitPtable[1][i]/dk)-peaksize;
				peakmax = peakmin+2*peaksize+1;
				filter[i] = FourierFilterMethods.getSuppressionFieldRing(peakmin, peakmax-1, 2*t.nlayers, 2*t.nlayers, false);
			
			}
			else
				filter[i] = FourierFilterMethods.getSuppressionFieldRing(peakmin, peakmax-1, 2*t.nlayers, 2*t.nlayers, false);
		}
		return filter;
	}
	/**
	 * rF is [nlayers][2];
	 *
	 */
	public static Topomap getMaskingLayers(PointImp[] imps, Topomap t, int origIndex, double[][] rF, double gaussradius){
		double[][] dr = new double [t.nlayers][2];
		for (int i = 0; i < t.nlayers; i++){
			dr[i][0] = rF[i][0] - rF[origIndex][0];
			dr[i][1] = rF[i][1] - rF[origIndex][1];
		}
		
		double[][][] mask = new double [t.nlayers][t.nx][t.ny];
		for (int i  = 0; i < mask.length; i++){
			PointImp[] temp = new PointImp[imps.length];
			for (int j = 0; j < imps.length; j++)
				temp[j] = new PointImp(new double[] {imps[j].pixelPos[0] + dr[i][0], imps[j].pixelPos[1] + dr[i][1]});
			mask[i] = ImpurityUtil.getMaskingLayer(t.getLayer(i), gaussradius, temp);
		}
		
		return Topomap.newTopomap(t, mask);
			
	}
	public static double[][] radialFFTStripCutEachEnergy(Topomap t, double stripwidth, int nphipts)
	{
		StripCut[] s = new StripCut[t.nlayers];
		double[][] results = new double[t.nlayers][];
		FieldUtil.RadiallyAveragedStripCut[] rc = new FieldUtil.RadiallyAveragedStripCut[t.nlayers];
		
		System.out.print("Taking FFTs (out of " + t.nlayers + "): ");
		for (int i = 0; i < t.nlayers; i++)
		{
			s[i] = new StripCut(FFTOps.obtainFFTmagCent(t.data[i]), stripwidth);
			rc[i] = new FieldUtil.RadiallyAveragedStripCut(s[i], nphipts);
//			results[i] = new double [s[i].result.length];
			System.out.print("" + (i+1) + " ");
		}
		System.out.println();
		System.out.print("Taking cuts (out of " + t.nlayers + "): ");
		for (int i = 0; i < t.nlayers; i++)
		{
			rc[i].calculate();
			System.out.print("" + (i+1) + " ");
			results[i] = rc[i].average;
		}
		System.out.println();
		return results;
	}
	public static double[][] stripCutEachEnergy(Topomap t, double stripwidth, double[] unitVector)
	{
		StripCut[] s = new StripCut[t.nlayers];
		double[][] results = new double[t.nlayers][];
		for (int i = 0; i < t.nlayers; i++)
		{
			s[i] = new StripCut(t.data[i], stripwidth);
			s[i].setPhi(Complex.phase(unitVector));
			s[i].makeCut();
			results[i] = s[i].result;
		}
		return results;
	}
	public static double[][] stripCutEachEnergy(Topomap t, double stripwidth, int[] origin, int[] destination)
	{
		StripCut[] s = new StripCut[t.nlayers];
		double[][] results = new double[t.nlayers][];
		for (int i = 0; i < t.nlayers; i++)
		{
			s[i] = new StripCut(t.data[i], stripwidth, origin, destination);
			s[i].makeCut();
			results[i] = s[i].result;
		}
		return results;
	}
	public static double[][][] polarCoordFFTEachEnergy(Topomap t, double stripwidth, int nphipts)
	{
		StripCut[] s = new StripCut[t.nlayers];
		double[][][] results = new double[t.nlayers][][];
		FieldUtil.RadiallyAveragedStripCut[] rc = new FieldUtil.RadiallyAveragedStripCut[t.nlayers];
		
		System.out.print("Taking FFTs (out of " + t.nlayers + "): ");
		for (int i = 0; i < t.nlayers; i++)
		{
			s[i] = new StripCut(FFTOps.obtainFFTmagCent(t.data[i]), stripwidth);
			rc[i] = new FieldUtil.RadiallyAveragedStripCut(s[i], nphipts);
//			results[i] = new double [s[i].result.length];
			System.out.print("" + (i+1) + " ");
		}
		System.out.println();
		System.out.print("Taking cuts (out of " + t.nlayers + "): ");
		for (int i = 0; i < t.nlayers; i++)
		{
			rc[i].calculate();
			System.out.print("" + (i+1) + " ");
			results[i] = rc[i].allresults;
		}
		System.out.println();
		return results;
	}

	
	public static void impurityMaskPercentile(Topomap t, String homedir, double lowerperc, double upperperc, String append)
	{
		String outdir = homedir.substring(0, homedir.length() - 1) + append + "\\";
		if (!new File(outdir).exists())
			new File(outdir).mkdir();
		double[][] data;
		for (int i = 0; i < t.nlayers; i++)
		{
			data = FieldOps.copy(t.data[i]);
			FieldOps.cutOffExtremes(data, lowerperc, upperperc);
			RHKFileOps.write3Files(outdir, data, t.names[i]);
		}
	}
	/**
	 * Returns two topomaps. The first is the impurity-masked topomap; the second is the cutoff-pixel topomap
	 * @param t
	 * @param homedir
	 * @param lowerperc
	 * @param upperperc
	 * @param append
	 * @return
	 */
	public static Topomap[] impurityMaskPercentile(Topomap t, double lowerperc, double upperperc)
	{
		double[][][][] data = new double [2][t.nlayers][t.nx][t.ny];
		for (int i = 0; i < t.nlayers; i++)
		{
			data[0][i] = FieldOps.copy(t.data[i]);
			data[1][i] = FieldOps.cutOffExtremesPixels(data[0][i], lowerperc, upperperc);
		}
		return new Topomap[] {Topomap.newTopomap(t, data[0]), Topomap.newTopomap(t, data[1])};
	}
	/**
	 * Same as above, but fourier-filters out the long-wavelength modes before eliminating the extreme pixels.
	 * This should tend to eliminate noisy pixels rather than bright patches.
	 * @param t
	 * @param lowerperc
	 * @param upperperc
	 * @param ftRadius
	 * @return
	 */
	public static Topomap[] impurityMaskPercentileHighFreq(Topomap t, double lowerperc, double upperperc, double ftRadius)
	{
		double[][][][] data = new double [3][t.nlayers][t.nx][t.ny];
		for (int i = 0; i < t.nlayers; i++)
		{
			data[0][i] = FFTOps.getFourierFilteredIFFT(t.data[i], FourierFilterMethods.getSuppressionFieldCircle(ftRadius, t.nx, t.ny));
			data[2][i] = FFTOps.getFourierFilteredIFFT(t.data[i], FourierFilterMethods.getSuppressionFieldNotCircle(ftRadius, t.nx, t.ny));
			data[1][i] = FieldOps.cutOffExtremesPixels(data[0][i], lowerperc, upperperc);
			FieldOps.add(data[2][i], data[0][i], data[0][i]);
		}
		return new Topomap[] {Topomap.newTopomap(t, data[0]), Topomap.newTopomap(t, data[1])};
	}
	
	public static double[][][] correlationMapEachEnergy(Topomap t, double[][] imp)
	{
		double[][][] corrmap = new double[t.nlayers][t.nx][t.ny];
		for (int i = 0; i < t.nlayers; i++)
		{
			corrmap[i] = FieldOps.getCorrelationMapCent(imp, t.data[i]);
			System.out.print(i + " ");
		}
		System.out.println();
		return corrmap;
	}
	public static double[][][] unnormalizedCorrelationMapEachEnergy(Topomap t, double[][] imp)
	{
		double[][][] corrmap = new double[t.nlayers][t.nx][t.ny];
		for (int i = 0; i < t.nlayers; i++)
		{
			corrmap[i] = FieldOps.getCorrelationMapCent(imp, t.data[i]);
			System.out.print(i + " ");
		}
		System.out.println();
		return corrmap;
	}
	
	//returns a new topomap where the data block is replaced by the layer files composing the series.
	public static Topomap loadSeries(Topomap t)
	{
		JFileChooser fc = new JFileChooser("C:\\data\\");
		File f = null;
		JOptionPane.showMessageDialog(null, "Select first layer.");
		if (fc.showOpenDialog(null) == JFileChooser.APPROVE_OPTION) 
			f = fc.getSelectedFile();
		double[][] field = ColumnIO.readSquareTable(f.toString());
		//we will read them all in now:
		String s = f.toString();
		String dir = s.substring(0, s.lastIndexOf("\\")+1);
		
		int firstIndex = Integer.parseInt(s.substring(dir.length(), s.indexOf(".dat")));
		int index = firstIndex;
		while (new File(dir + MovieMaker.fromInt(index)+".dat").exists())
			index++;
		double[][][] dataset = new double[index-firstIndex][field.length][field[0].length];
		
		dataset[0] = field;
		System.out.print("Read (out of " + dataset.length + ")");
		for (int i = 1; i < dataset.length; i++){ System.out.print("" + i + " ");
			dataset[i] = ColumnIO.readSquareTable(dir + MovieMaker.fromInt(firstIndex+i)+".dat");
		}
		System.out.println();
		
		return Topomap.newTopomap(t, dataset);
	}
	
	//this returns a set of curves of intensity vs energy evaluated at each point in the points array, which is [n][2].
	public static double[][] getIntensityAroundPoints(Topomap t, double[][] points, double radius)
	{
		double[][] intensity = new double [points.length][t.nlayers];
		for (int i = 0; i < points.length; i++)
			for (int j = 0; j < t.nlayers; j++)
			{
				intensity[i][j] = FieldOps.getLocalAvgCircle(t.data[j], points[i][0], points[i][1], radius);
			}
		return intensity;
	}
	public static double[][] getCorrelationAtPoints(Topomap t, double[][] points, double[][] imp)
	{
		double[][] intensity = new double [points.length][t.nlayers];
		for (int i = 0; i < points.length; i++)
			for (int j = 0; j < t.nlayers; j++)
			{
				intensity[i][j] = FieldOps.correlation(imp, t.data[j], (int)(points[i][0]+0.5), (int)(points[i][1]+0.5));
			}
		return intensity;
	}
	
	public static void writeIntensityAtImpurities(double radius)
	{
		Topomap t = Topomap.open();
		ArrayList<ImpurityListEditor.Impurity> list = ImpurityListEditor.getImpurities(1);
		double[][] pts = new double [list.size()][2];
		for (int i = 0; i < list.size(); i++)
		{
			pts[i] = list.get(i).getPosition();
		}
		double[][] intensityCurves = getIntensityAroundPoints(t, pts, radius);
		FileOps.writeTableASCII(intensityCurves);
	}
	public static void writeCorrelationAtImpurities()
	{
		Topomap t = Topomap.open();
		double[][] imp = FileOps.openBin(Topomap.stddir);
		ArrayList<ImpurityListEditor.Impurity> list = ImpurityListEditor.getImpurities(1);
		double[][] pts = new double [list.size()][2];
		for (int i = 0; i < list.size(); i++)
		{
			pts[i] = list.get(i).getPosition();
		}
		double[][] intensityCurves = getCorrelationAtPoints(t, pts, imp);
		FileOps.writeTableASCII(intensityCurves);
	}

	//this assumes that the ufield is already the same size as a layer
	public static Topomap applyUField(double[][][] u, Topomap t)
	{
		double[][][] after = new double[t.nlayers][t.nx][t.ny];
		
		double[][][] uActual = u;
//		double[][] temp = new double [t.nx][t.ny];
		
		if (u.length != t.nx){
			uActual = new double [t.nx][t.ny][2];
			double ratio = u.length/t.nx;
			if (ratio > 1)
				FieldOps.reduce((int)ratio, u, uActual);
		}
		else System.out.println("same size");
		
		for (int i = 0; i < t.nlayers; i++)
		{
//			FieldOps.applyUFieldSmooth(t.data[i], uActual, 2, 2, after[i]);
			FieldOps.applyUFieldBiCubic(t.data[i], uActual, after[i]);
			FieldOps.changeZeroToAverage(after[i]);
		}
		Topomap ta = Topomap.newTopomap(t, after);
		return ta;
	}
	public static Topomap applyUField(boolean fancy, double[][][] u, Topomap t)
	{
		double[][][] after = new double[t.nlayers][t.nx][t.ny];
		
		double[][][] uActual = u;
//		double[][] temp = new double [t.nx][t.ny];
		
		if (u.length != t.nx){
			uActual = new double [t.nx][t.ny][2];
			double ratio = u.length/t.nx;
			if (ratio > 1)
				FieldOps.reduce((int)ratio, u, uActual);
		}
		else System.out.println("same size");
		
		for (int i = 0; i < t.nlayers; i++)
		{
//			FieldOps.applyUFieldSmooth(t.data[i], uActual, 2, 2, after[i]);
			if (!fancy){
			FieldOps.applyUFieldBiCubic(t.data[i], uActual, after[i]);
			FieldOps.changeZeroToAverage(after[i]);
			}
			else{
				after[i] = DriftCorrectionMethods.applyUFieldSpecial(u, t.data[i]);
				boolean[][] outsidePixels = FieldOps.getOutsidePixels(u);
				FieldOps.zero(after[i], outsidePixels);
				FieldOps.changeZeroToAverage(after[i]);

			}
		}
		
		Topomap ta = Topomap.newTopomap(t, after);
		return ta;
	}
	
	/**
	 * This assumes u has same number of layers as the topomap.
	 * @param u
	 * @param t
	 * @return
	 */
	public static Topomap applyUField(double[][][][] u, Topomap t)
	{
		double[][][] after = new double[t.nlayers][t.nx][t.ny];
		
		double[][][][] uActual = u;
//		double[][] temp = new double [t.nx][t.ny];
		
//		if (u[0].length != t.nx){
//			uActual = new double [t.nlayers][t.nx][t.ny][2];
//			double ratio = u[0].length/t.nx;
//			if (ratio > 1)
//				FieldOps.reduce((int)ratio, u, uActual);
//		}
//		else System.out.println("same size");
		
		for (int i = 0; i < t.nlayers; i++)
		{
//			FieldOps.applyUFieldSmooth(t.data[i], uActual, 2, 2, after[i]);
			FieldOps.applyUFieldBiCubic(t.data[i], uActual[i], after[i]);
			FieldOps.changeZeroToAverage(after[i]);
		}
		Topomap ta = Topomap.newTopomap(t, after);
		return ta;
	}
	public static Topomap shiftLeftOnePixel(Topomap t)
	{
		double[][][] data = new double [t.nlayers][t.nx][t.ny];
		for (int i = 0; i < t.nlayers; i++)
			for (int j = 0; j < t.nx-1; j++)
				for (int k = 0; k < t.ny; k++)
				{
					data[i][j][k] = t.data[i][j+1][k];
				}
		int j = t.nx-1;
		for (int i = 0; i < t.nlayers; i++)
				for (int k = 0; k < t.ny; k++)
				{
					data[i][j][k] = t.data[i][j][k];
				}
		return new Topomap(data, t.v, t.x, t.y, null);
	}
	
//	public static void addOrigin(Topomap t)
//	{
//		double xc = Double.parseDouble(JOptionPane.showInputDialog("Enter the x-coordinate of the center."));
//		double yc = Double.parseDouble(JOptionPane.showInputDialog("Enter the y-coordinate of the center."));
//		t.makeOrigin(new double[] {xc, yc});
//	}
	public static Topomap expandBi(Topomap t, int factor)
	{
		int nnx = t.nx*factor, nny = t.ny*factor;
		double[] x = new double [nnx];
		double[] y = new double [nny];
		for (int i = 0; i < nnx; i++)
			x[i] = t.x[0] + (t.xLength/nnx)*i;
		for (int i = 0; i < nny; i++)
			y[i] = t.y[0] + (t.yLength/nny)*i;
		
		double[][][] data = new double[t.nlayers][][];
		for (int i = 0; i < t.nlayers; i++)
			data[i] = FieldOps.expandBi(t.data[i], factor);
		return new Topomap(data, t.v, x, y, null);
	}
	public static Topomap expandDouble(Topomap t)
	{
		int nnx = t.nx*2, nny = t.ny*2;
		double[] x = new double [nnx];
		double[] y = new double [nny];
		for (int i = 0; i < nnx; i++)
			x[i] = t.x[0] + (t.xLength/nnx)*i;
		for (int i = 0; i < nny; i++)
			y[i] = t.y[0] + (t.yLength/nny)*i;
		
		double[][][] data = new double[t.nlayers][][];
		for (int i = 0; i < t.nlayers; i++)
			data[i] = FieldOps.expand(t.data[i]);
		return new Topomap(data, t.v, x, y, null);
	}
	/**
	 * If lattice is square
	 * @param t
	 * @param latt
	 * @return
	 */
	public static Topomap makeSymmetrizedFFTs(Topomap t, AtomicCoordinatesSet latt, boolean square)
	{
		double[][][] newData = new double [t.nlayers][t.nx][t.ny];
		boolean log = JOptionPane.showConfirmDialog(null, "Use log scale?") == JOptionPane.YES_OPTION;
//		if (!square) ColumnIO.writeString(LayerUtil.getSymmetrizedTriangleLattice(latt, t.getLayer(0)).toString(), FileOps.selectSave(null).toString());
		for (int i = 0; i < t.nlayers; i++){
			if (square)
				newData[i] = LayerUtil.symmetrizeFFT_1(latt, t.getLayer(i), log);
			else
				newData[i] = LayerUtil.symmetrizeFFTTriang(latt, t.getLayer(i), log);
//			FieldOps.log(newData[i]);
		}
		
		double dx = 2*Math.PI/t.xLength;
		double dy = 2*Math.PI/t.yLength;
		double[] x = new double [t.nx];
		double[] y = new double [t.ny];
		for (int i = 0; i < t.nx; i++)
			x[i] = (i-t.nx/2)*dx;
		for (int j = 0; j < t.ny; j++)
			y[j] = (j-t.ny/2)*dy;
		
		return new Topomap(newData, t.v, x, y, null);
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
	//This sorts the spectra by signal at the indicated bias, and obtains the cutoffs between the groups.
	public static double[] obtainSignalRegionGroups(int biasIndex, int ngroups, Topomap t)
	{
		double[] signals = new double [t.nx*t.ny];
		for (int i = 0; i < signals.length; i++)
			signals[i] = t.data[biasIndex][i%t.nx][i/t.nx];
		ArrayOps.quicksort(signals);
		
		int[] bounds = new int [ngroups];
		for (int i = 0; i < bounds.length; i++)
			bounds[i] = (((i+1)*signals.length)/ngroups)-1;
		
		double[] values = new double [ngroups];
		for (int i = 0; i < bounds.length; i++)
			values[i] = signals[bounds[i]];
		return values;
	}
	public static Topomap mergeMaps(Topomap[] t)
	{
		double[] vRaw, vFinal;
		int nlayers = 0;
		for (int i = 0; i < t.length; i++)
			nlayers += t[i].nlayers;
		
		vRaw = new double [nlayers];
		vFinal = new double [nlayers];
		int nx = t[0].nx, ny = t[0].ny;
		int n = 0;
		int[] is = new int [nlayers];
		int[] js = new int [nlayers];
		for (int i = 0; i < t.length; i++)
		for (int j = 0; j < t[i].nlayers; j++){
			vRaw[n] = t[i].v[j];
			vFinal[n] = vRaw[n];
			is[n] = i;
			js[n] = j;
			n++;
		}
		
		ArrayOps.quicksort(vFinal);
		double[][][] data =new double[nlayers][][];
		for (int i = 0; i < nlayers; i++)
			for (int j = 0; j < nlayers; j++)
				if (vFinal[i] == vRaw[j])
					data[i] = t[is[j]].data[js[j]];
		return new Topomap(data, vFinal, t[0].x, t[0].y, null);
	}
	public static void invert(Topomap t)
	{
		double[] newV = new double[t.nlayers];
		double[][][] newData = new double [t.nlayers][][];
		for (int i = 0; i < t.nlayers; i++)
		{
			newV[i] = t.v[t.nlayers-1-i];
			newData[i] = t.data[t.nlayers-1-i];
		}
		t.v = newV;
		t.data = newData;
	}
	public static Topomap fourierFilter(Topomap t, ArrayList<int[]> includedPts, boolean include)
	{
		double[][][] newData = new double [t.nlayers][t.nx][t.ny];
		FFT2DSmall fft;
		fft = FFTOps.obtainFFT(t.data[0]);
		
		int[] x = new int [includedPts.size()], y = new int [includedPts.size()];
		for (int i = 0; i < includedPts.size(); i++)
		{
			x[i] = includedPts.get(i)[0];
			y[i] = includedPts.get(i)[1];
		}
		
		double minMag = FieldOps.magMin(fft.fHat);
//		Systme.out.println(FieldOps.magMin(fft.fHat));
//		double minMag;
		if (minMag == 0) minMag = Math.exp(-20);
		System.out.println(Math.log(minMag));
		
		for (int i = 0; i < t.nlayers; i++)
		{
			fft = FFTOps.obtainFFT(FieldOps.copy(t.data[i]));
			if (include)	
				FFTOps.supressElseModesTo(fft, x, y, minMag);
			else
				FFTOps.supressModesTo(fft, x, y, minMag);
				
			fft.doIFFT();
			for (int j = 0; j < t.nx; j++)
				for (int k = 0; k < t.ny; k++)
					newData[i][j][k] = fft.f[j][k][0];
		}
		return Topomap.newTopomap(t, newData);
	}
	public static Topomap fourierFilter(Topomap t, boolean[][] suppress)
	{
		double[][][] newData = new double [t.nlayers][t.nx][t.ny];
		FFT2DSmall fft;
		fft = FFTOps.obtainFFT(t.data[0]);
		
		double minMag = FieldOps.magMin(fft.fHat);
//		Systme.out.println(FieldOps.magMin(fft.fHat));
//		double minMag;
		if (minMag == 0) minMag = Math.exp(-20);
		System.out.println(Math.log(minMag));
		
		for (int i = 0; i < t.nlayers; i++)
		{
			fft = FFTOps.obtainFFT(FieldOps.copy(t.data[i]));
			FFTOps.supressModesTo(fft, suppress, minMag);
			fft.doIFFT();
			for (int j = 0; j < t.nx; j++)
				for (int k = 0; k < t.ny; k++)
					newData[i][j][k] = fft.f[j][k][0];
		}
		return Topomap.newTopomap(t, newData);
	}
	/**
	 * This one uses a 3D suppression array so that the suppression region can "disperse."
	 * @param t
	 * @param suppress
	 * @return
	 */
	public static Topomap fourierFilter(Topomap t, boolean[][][] suppress)
	{
		double[][][] newData = new double [t.nlayers][t.nx][t.ny];
		FFT2DSmall fft;
		fft = FFTOps.obtainFFT(t.data[0]);
		
		double minMag = FieldOps.magMin(fft.fHat);
//		Systme.out.println(FieldOps.magMin(fft.fHat));
//		double minMag;
		if (minMag == 0) minMag = Math.exp(-20);
		System.out.println(Math.log(minMag));
		
		for (int i = 0; i < t.nlayers; i++)
		{
			fft = FFTOps.obtainFFT(FieldOps.copy(t.data[i]));
			FFTOps.supressModesTo(fft, suppress[i], minMag);
			fft.doIFFT();
			for (int j = 0; j < t.nx; j++)
				for (int k = 0; k < t.ny; k++)
					newData[i][j][k] = fft.f[j][k][0];
		}
		return Topomap.newTopomap(t, newData);
	}
	
	public static Topomap spatialFilter(Topomap t, boolean[][] preserve)
	{
		double mean = 0;
		double[][][] ans = new double[t.nlayers][t.nx][t.ny];
		for (int i = 0; i < t.nlayers; i++)
		{
			mean = FieldOps.mean(t.data[i], preserve);
			for (int j = 0; j < t.nx; j++)
				for (int k = 0; k < t.ny; k++)
					if (!preserve[j][k])
						ans[i][j][k] = mean;
					else
						ans[i][j][k] = t.data[i][j][k];
		}
		return Topomap.newTopomap(t, ans);
	}
	public static Topomap spatialFilter(Topomap t, double[][] weight)
	{
		double mean = 0;
		double[][][] ans = new double[t.nlayers][t.nx][t.ny];
		for (int i = 0; i < t.nlayers; i++)
		{
			mean = FieldOps.mean(t.data[i], weight);
			for (int j = 0; j < t.nx; j++)
				for (int k = 0; k < t.ny; k++)
				{
					ans[i][j][k] = weight[j][k]*t.data[i][j][k] + (1-weight[j][k])*mean;
				}
		}
		return Topomap.newTopomap(t, ans);
	}
	
	//fits each spectrum in the map to the function specified by fname, and write various necessary files into the directory dir.
	public static void fitEachSpectrum(Topomap t, String dir, String fname)
	{
		ACM_CustomFunctions f = ACM_CustomFunctions.getNew(fname);
		
		File outd = new File (dir);
		if (!outd.exists()) outd.mkdir();

		double[][][] fitParams = new double [f.plist.length][t.nx][t.ny];
		double[][][] result = new double [t.nlayers][t.nx][t.ny];
		double[] spectrum, fitPs, estimate;
		for (int i = 0; i < t.nx; i++){//System.out.println();
			for (int j = 0; j < t.ny; j++)
			{
				if (j == 0) System.out.println(Printer.vectorP(i, j) + "\t");
				fitPs = ACM_NonLinearFitter.fitToFunction(t.v, t.getSpectrum(i, j), fname);
				estimate = ACM_CustomFunctions.getExpectedValues(t.v, fitPs, fname);
				for (int k = 0; k < t.nlayers; k++)
					result[k][i][j] = estimate[k];
				for (int k = 0; k < fitPs.length; k++)
					fitParams[k][i][j] = fitPs[k];
			}
		}
		
		//output:
		Topomap.writeBIN(Topomap.newTopomap(t, result), dir + "Fitting Curves" + ".bin");
		
		for (int i = 0; i < t.nlayers; i++)
			result[i] = FieldOps.minus(t.data[i], result[i]);
		
		Topomap.writeBIN(Topomap.newTopomap(t, result), dir + "Difference values" + ".bin");
		
		
		for (int i = 0; i < fitParams.length; i++)
		{
			Layer.writeBIN(Layer.newLayer(t, fitParams[i]), dir + "Fit " + f.plist[i] + ".bin");
		}
		
	}
	//This is for a linear regression to a simple exponential
	public static void fitEachSpectrumExponential(Topomap t, String dir)
	{
		//take the absolute value
		Topomap copy = Topomap.newTopomap(t, t.data.clone());
		for (int i = 0; i < t.nlayers; i++)
		{
			FieldOps.abs(copy.data[i]);
			FieldOps.log(copy.data[i]);
		}
		
		double[][] tempData = new double [t.nlayers][2];
		for (int i = 0; i < t.nlayers; i++)
			tempData[i][0] = t.v[i];
		
		File outd = new File (dir);
		if (!outd.exists()) outd.mkdir();

		double[][][] fitParams = new double [2][t.nx][t.ny];
		double[][][] result = new double [t.nlayers][t.nx][t.ny];
		double[] spectrum, fitPs, estimate;
		SimpleRegression reg = new SimpleRegression(true);
		RegressionResults results;
		double[] regRes;
		for (int i = 0; i < t.nx; i++){//System.out.println();
			for (int j = 0; j < t.ny; j++)
			{
				reg.clear();
				for (int k = 0; k < t.nlayers; k++)
					tempData[k][1] = copy.data[k][i][j];
				reg.addData(tempData);
				
				results = reg.regress();
				fitPs = results.getParameterEstimates();
//				Printer.printlnHorizontal(regRes);
				for (int k = 0; k < t.nlayers; k++)
					result[k][i][j] = Math.exp(fitPs[0] + fitPs[1]*t.v[k]);
				fitParams[0][i][j] = fitPs[0];
				fitParams[1][i][j] = fitPs[1];
				//				for (int k = 0; k < t.nlayers; k++)
//					result[k][i][j] = estimate[k];
//				for (int k = 0; k < fitPs.length; k++)
//					fitParams[k][i][j] = fitPs[k];
			}
		}
		
		//output:
		if (!dir.equals("")){
			Topomap.writeBIN(Topomap.newTopomap(t, result), dir + "Fitting Curves" + ".bin");
			
			for (int i = 0; i < t.nlayers; i++)
				result[i] = FieldOps.minus(FieldOps.getAbs(t.data[i]), result[i]);
			
			Topomap.writeBIN(Topomap.newTopomap(t, result), dir + "Difference values" + ".bin");
			
			
			Layer.writeBIN(Layer.newLayer(t, fitParams[0]), dir + "Fit Front Factor.bin");
			Layer.writeBIN(Layer.newLayer(t, fitParams[1]), dir + "Fit decay length Factor.bin");
		}
		double hbar = 6.58211928e-16;
		double c = 299792458;
		double me = 0.510998910e6/(c*c);
		double[][] workFunction = new double [t.nx][t.ny];
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
			{
				workFunction[i][j] = hbar*hbar*fitParams[1][i][j]*fitParams[1][i][j]/(8*me);
			}
		Layer wf = Layer.newLayer(t, workFunction);
		if (!dir.equals(""))
			Layer.writeBIN(wf, dir + "fit Work function.bin");
		else
			Layer.writeBIN(wf);
		
		new LayerViewer(wf, dir, 512);
	}
	
	public static Topomap[] subtractBackground(Topomap topo, int npoly)
	{
		Regression reg;
		double[] g = new double[topo.nlayers];
		double[] a;
		double[][][] poly = new double [topo.nlayers][topo.nx][topo.ny];
		double[][][] subt = new double [topo.nlayers][topo.nx][topo.ny];
		double[][] powers = new double[topo.nlayers][npoly+1];
		for (int k = 0; k < topo.nlayers; k++)
			for (int i = 0; i < npoly+1; i++)
				powers[k][i] = Math.pow(topo.v[k], i);
			
		
		for (int j = 0; j < topo.ny; j++)
			for (int i = 0; i < topo.nx; i++)
			{
				for (int k = 0; k < topo.nlayers; k++)
					g[k] = topo.data[k][i][j];
				reg = new Regression(topo.v,g);
				reg.polynomial(npoly);
				a = reg.getBestEstimates();
//				System.out.print(i + ",\t" + j + ",\t");
//				for (int k = 0; k  < a.length; k++)
//					System.out.print(a[k] + ",\t");
//				System.out.println();
				for (int k = 0; k < topo.nlayers; k++)
				{
					for (int m = 0; m < npoly+1; m++)
						poly[k][i][j] += a[m]*powers[k][m];
					subt[k][i][j] = topo.data[k][i][j]-poly[k][i][j];
				}
			}
		return new Topomap[] {new Topomap(subt, topo.v, topo.x, topo.y, null), new Topomap(poly, topo.v, topo.x, topo.y, null)};
	}
//	public static double[][] getFFTRadialAveragedSpec(Topomap t)
//	{
//		double[][][] fft = t.getFFT(false);
//		StripCut s = new StripCut(fft[0], 1);
//		FieldUtil.RadiallyAveragedStripCut rc = new FieldUtil.RadiallyAveragedStripCut(s, 1000);
//		
//		double[][] ans = new double [t.nlayers][s.npts];
//		for (int i = 0; i < t.nlayers; i++)
//		{
//			s.changeMap(fft[i]);
//			rc.calculate();
//			ans[i] = FieldOps.copy(rc.average);
//		}
//		return ans;
//		
//	}
	public static Topomap getZAdjusteddIdV(Topomap t, Layer z, double k)
	{
		double[][][] ans = new double[t.nlayers][t.nx][t.ny];
		for (int i = 0; i < t.nlayers; i++)
		{
			ans[i] = FieldOps.getZAdjusteddIdV(t.data[i], z.data, k, true);
		}
		return Topomap.newTopomap(t, ans);
	}
	/**
	 * This one is for actual "topomaps" where the z has as many layers as the di/dV. We won't try to get an absolute reference for Z,
	 * but simply do each layer separately.
	 * @param t
	 * @param z
	 * @param k
	 * @return
	 */
	public static Topomap getZAdjusteddIdV(Topomap t, Topomap z, double k)
	{
		double[][][] ans = new double[t.nlayers][t.nx][t.ny];
		for (int i = 0; i < t.nlayers; i++)
		{
			ans[i] = FieldOps.getZAdjusteddIdV(t.data[i], z.data[i], k, true);
		}
		return Topomap.newTopomap(t, ans);
	}

	
	public static Topomap addTopographyTo(Topomap t, Layer topo)
	{
		Layer[] combo = new Layer[t.nlayers+1];
		combo[0] = topo;
		for (int i = 1; i < t.nlayers; i++)
		{
			combo[i] = t.getLayer(i-1);
		}
		return Topomap.newTopomap(combo);
	}
	public static Topomap truncateBias(Topomap t, int imin, int imax)
	{
		Layer[] combo = new Layer[imax-imin];
		for (int i = 0; i < combo.length; i++)
		{
			combo[i] = t.getLayer(i+imin);
		}
		return Topomap.newTopomap(combo);
	}
	
	/**
	 * This arranges for all maps to have the same minimum. For our purposes, we will shift all layers UP.
	 * @param t
	 * @param minSource
	 */
	public static void writeShiftedAccordingToMinimum(Topomap t, Topomap minSource, JFileChooser fc)
	{
		int[][] min = new int [t.nx][t.ny];
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
				min[i][j] = ArrayOps.minIndex(minSource.getSpectrum(i, j));
		
		int minmin = ArrayOps.max(min);
		FieldOps.minusEquals(min, minmin);
		
		Topomap ans = shiftEachSpectrum(t, min);
		Topomap.writeBIN(ans, fc);

	}
	public static Topomap shiftEachSpectrum(Topomap t, int[][] shifts)
	{
		double[][][] data = new double [t.nlayers][t.nx][t.ny];
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
			{
				for (int k = 0; k < t.nlayers; k++)
				{
					if (k + shifts[i][j] >= 0 && k + shifts[i][j] < t.nlayers)
						data[k][i][j] = t.data[k+shifts[i][j]][i][j];
					else
						data[k][i][j] = 0;
				}
			}
		
		for (int k = 0; k < t.nlayers; k++)
			FieldOps.changeZeroToAverage(data[k]);
		
		return Topomap.newTopomap(t, data);
	}
	public static Topomap shiftEachSpectrum(Topomap t, double[][] shifts_metric)
	{
		double[][] shifts_pixel = FieldOps.copy(shifts_metric);
		FieldOps.timesEquals(shifts_pixel, 1/(t.v[1]-t.v[0]));//being t.ny;
//		Layer.writeBIN(Layer.getFreeLayer(shifts_pixel), fc);
		double[][][] data = new double [t.nlayers][t.nx][t.ny];
		double[] q = new double [t.nlayers];
		double[][] expqr = new double [t.nlayers][2];
		for (int i = 0; i < t.nlayers; i++){
			q[i] = 2*Math.PI*i/t.nlayers;
			if (q[i] >= Math.PI) q[i] -= 2*Math.PI;
		}
		
		double[] tempZ = new double [2];
		double[] tempFFTZ = new double[2];
//		double[][][] tempNonsense = new double [2*t.nlayers][t.nx][t.ny];
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
			{
				double[] fftz = FFTOps.get1DFFTComplex(t.getSpectrum(i, j));
				for (int k = 0; k < t.nlayers; k++)
				{
					expqr[k][0] = Math.cos(q[k]*shifts_pixel[i][j]);
					expqr[k][1] = Math.sin(-q[k]*shifts_pixel[i][j]);
					tempFFTZ[0] = fftz[2*k];
					tempFFTZ[1] = fftz[2*k+1];
					Complex.product(expqr[k], tempFFTZ, tempZ);
					fftz[2*k] = tempZ[0];
					fftz[2*k+1] = tempZ[1];
//					tempNonsense[2*k][i][j] = fftz[2*k];
//					tempNonsense[2*k+1][i][j] = fftz[2*k+1];
				}
				double[] newspec = FFTOps.getIFFTReal(fftz);
				for (int k = 0; k < t.nlayers; k++)
					data[k][i][j] = newspec[k];
			}
				
//		for (int k = 0; k < t.nlayers; k++)
//			FieldOps.changeZeroToAverage(data[k]);
//		Topomap.writeBIN(new Topomap(tempNonsense, ArrayOps.generateArray(0, 1, 2*t.nlayers), t.x, t.y, null), fc);
		return Topomap.newTopomap(t, data);
	}
	public static Topomap shiftEachSpectrum_fancy(Topomap t, double[][] shifts_metric)
	{
		double[][] shifts_pixel = FieldOps.copy(shifts_metric);
		FieldOps.timesEquals(shifts_pixel, 1/(t.v[1]-t.v[0]));//being t.ny;
//		Layer.writeBIN(Layer.getFreeLayer(shifts_pixel), fc);
		double[][][] data = new double [t.nlayers][t.nx][t.ny];
		double[] q = new double [t.nlayers];
		double[][] expqr = new double [t.nlayers][2];
		for (int i = 0; i < t.nlayers; i++){
			q[i] = 2*Math.PI*i/t.nlayers;
			if (q[i] >= Math.PI) q[i] -= 2*Math.PI;
		}
		double[] initSpec = new double [t.nlayers];
		double m, b;
		
		double[] tempZ = new double [2];
		double[] tempFFTZ = new double[2];
//		double[][][] tempNonsense = new double [2*t.nlayers][t.nx][t.ny];
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
			{
				initSpec = t.getSpectrum(i, j);
				m = (initSpec[t.nlayers-1] - initSpec[0])/(t.nlayers-1);
				b = initSpec[0];
				for (int k = 0; k < t.nlayers; k++)
					initSpec[k] -= m*k;// + b;
				double[] fftz = FFTOps.get1DFFTComplex(initSpec);
				for (int k = 0; k < t.nlayers; k++)
				{
					expqr[k][0] = Math.cos(q[k]*shifts_pixel[i][j]);
					expqr[k][1] = Math.sin(-q[k]*shifts_pixel[i][j]);
					tempFFTZ[0] = fftz[2*k];
					tempFFTZ[1] = fftz[2*k+1];
					Complex.product(expqr[k], tempFFTZ, tempZ);
					fftz[2*k] = tempZ[0];
					fftz[2*k+1] = tempZ[1];
////					tempNonsense[2*k][i][j] = fftz[2*k];
////					tempNonsense[2*k+1][i][j] = fftz[2*k+1];
				}
				double[] newspec = FFTOps.getIFFTReal(fftz);
				for (int k = 0; k < t.nlayers; k++){
					data[k][i][j] = newspec[k]+ m*(k);//+shifts_pixel[i][j]);// + b;
					
				}
			}
				
//		for (int k = 0; k < t.nlayers; k++)
//			FieldOps.changeZeroToAverage(data[k]);
//		Topomap.writeBIN(new Topomap(tempNonsense, ArrayOps.generateArray(0, 1, 2*t.nlayers), t.x, t.y, null), fc);
		return Topomap.newTopomap(t, data);
	}
	
	/**
	 * This creates a k-space "topomap" which is equal to the spectral function ( = -2*Im G) calculated using the dispersion relation and re and im parts
	 * of the self energy supplied as parameters. The re and im parts are assumed functions of the energy only.
	 * 
	 * The dispersion relation will be supplied as internal-method code by the user.
	 * @param v
	 * @param kx
	 * @param ky
	 * @param dispersion
	 * @param re
	 * @param im
	 * @return
	 */
	public static Topomap getSpectralFunction(double[] v, double[] kx, double[] ky, double[] re, double[] im){
		double[][][] data = new double [v.length][kx.length][ky.length];
		
		double bareEnergy;
		
		for (int k = 0; k < v.length; k++)
			for (int i = 0; i < kx.length; i++)
				for (int j = 0; j < ky.length; j++)
				{
					//Calculate the energy associated with this k-point (dispersion relation): 
					bareEnergy = v[0] + 2*0.5*(kx[i]*kx[i] + ky[j]*ky[j]);
					data[k][i][j] = FieldOps.getSpectralFunction(v[k], re[k], im[k], bareEnergy);
				}
		return new Topomap(data, v, kx, ky, null);
	}
	public static void subtractParabolicFit(Topomap[] u)
	{
		for (int i = 0; i < u.length; i++)
			for (int j = 0; j < u[0].nlayers; j++)
				FieldOps.subtractParabolicFit(u[i].data[j]);
	}

	/**
	 * This returns a topomap each layer of which is a layer of the JDOS according to a certain FieldOps method.
	 * @param dos
	 * @return
	 */
	public static Topomap getJDOS(Topomap dos)
	{
		double[][][] data = new double [dos.nlayers][][];
		for (int i = 0; i < dos.nlayers; i++)
		{
			System.out.println(i);
			data[i] = FieldOps.getJDOS_sameSize(dos.data[i]);
		}
		return Topomap.newTopomap(dos, data);
	}
	
	/**
	 * This takes a map and creates the Z map. The Z map is g(r, +E)/g(r, -E) and is defined in supplementary info of doi:10.1038/nature09169
	 * .
	 * @param map
	 * @return
	 */
	public static Topomap getZMap(Topomap t)
	{
		int zeroIndex = ArrayOps.indexOf(t.v, 0, true);
		int nlayers;
		for (nlayers = 0; (zeroIndex+nlayers < t.nlayers && zeroIndex-nlayers >= 0); nlayers++);
		System.out.println(nlayers);
		double[][][] ansData = new double [nlayers][t.nx][t.ny];
		for (int k = 0; k < nlayers; k++)
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
					ansData[k][i][j] = t.data[zeroIndex+k][i][j]/t.data[zeroIndex-k][i][j];
		
		double[] v = new double [nlayers];
		for (int i = 0; i < nlayers; i++)
			v[i] = t.v[i+zeroIndex];
		
		return new Topomap(ansData, v, t.x, t.y, null);
	}
	/**
	 * This gives a signal equal to (integral dE of T(r,E) from E-window to E+window) / 2*window. The ends of the map have truncated windows.
	 * @param t
	 * @param energyWindow
	 * @return
	 */
	public static Topomap getWindowMapAboutV(Topomap t, double energyWindow)
	{
		double dv = t.v[1]-t.v[0];
		int nSpaces = (int)(energyWindow/dv);
		int[] lowerLimit = new int [t.nlayers];
		int[] upperLimit = new int [t.nlayers];
		for (int i = 0; i < nSpaces; i++){
			lowerLimit[i] = 0;
			upperLimit[i] = 2*i;
		}
		for (int i = nSpaces; i < t.nlayers - nSpaces; i++)
		{
			lowerLimit[i] = i-nSpaces;
			upperLimit[i] = i+nSpaces;
		}
		for (int i = t.nlayers-nSpaces; i < t.nlayers; i++)
		{
			int delta = (t.nlayers-1)-i;
			lowerLimit[i] = (t.nlayers-1) - 2*delta;
			upperLimit[i] = t.nlayers-1;
		}
		double[][][] ansData = new double [t.nlayers][t.nx][t.ny];
		for (int k = 0; k < t.nlayers; k++)
		{
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
				{
					for (int p = lowerLimit[k]; p <= upperLimit[k]; p++)
						ansData[k][i][j] += t.data[p][i][j];
					ansData[k][i][j] /= (upperLimit[k]-lowerLimit[k]+1);
				}
		}
		return Topomap.newTopomap(t, ansData);
	}
	public static Topomap getCurrentMapAboutV(Topomap t, double V)
	{
		int zeroIndex = ArrayOps.indexOf(t.v, V, true);
		double[][][] ansData = new double [t.nlayers][t.nx][t.ny];
		for (int k = 0; k < zeroIndex; k++)
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++){
					for (int p = k; p <= zeroIndex; p++)
						ansData[k][i][j] += t.data[p][i][j];
					ansData[k][i][j] /= Math.abs(k-zeroIndex)+1;
				}
		for (int k = zeroIndex; k < t.nlayers; k++)
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++){
					for (int p = zeroIndex; p <= k; p++)
						ansData[k][i][j] += t.data[p][i][j];
					ansData[k][i][j] /= Math.abs(k-zeroIndex)+1;
				}
		
		return Topomap.newTopomap(t, ansData);
	}
	
	public static Layer getSpectralDistributionBasic(int sx, int ny, Topomap t)
	{
		double xmin = ArrayOps.min(t.v);
		double xmax = ArrayOps.max(t.v);
		double ymin = ArrayOps.min(t.data);
		double ymax = ArrayOps.max(t.data);
		
		double[][] histogram = new double[t.nlayers][ny];
		double[] y = ArrayOps.generateArrayInclBoth(ymax, ymin, ny);
		double[] coord;
		Layer l = new Layer(histogram, t.v, y, 1, 1);
		for (int k = 0; k < t.nlayers; k++)
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
				{
					coord = l.getPixelCoords(0, t.data[k][i][j]);
					histogram[k][FieldOps.round(coord[1])]++;
				}
		
		double epsilon = Double.parseDouble(JOptionPane.showInputDialog("Add epsilon?", "" + 0));
		//Now resize it by sx:
		double[] x = ArrayOps.generateArrayInclBoth(xmin, xmax, t.nlayers*sx);
		double[][] finalHistogram = new double [x.length][ny];
		if (sx == 1)
			for (int i = 0; i < x.length; i++)
				for (int j = 0; j < ny; j++)
					finalHistogram[i][j] = l.evaluateAtMetric(x[i], y[j]) + epsilon;
		else
			for (int i = 0; i < x.length; i++)
				for (int j = 0; j < ny; j++)
					finalHistogram[i][j] = l.data[i][j] + epsilon;
		return new Layer(finalHistogram, x, y, 1, 1);
				
	}
	public static Layer getSpectralDistributionBasicLimited(int sx, int ny, double lower, double upper, Topomap t)
	{
		double xmin = ArrayOps.min(t.v);
		double xmax = ArrayOps.max(t.v);
		double ymin = ArrayOps.min(t.data);
		double ymax = ArrayOps.max(t.data);
		
		double[][] histogram = new double[t.nlayers][ny];
		double[] y = ArrayOps.generateArrayInclBoth(upper, lower, ny);
		double[] coord;
		Layer l = new Layer(histogram, t.v, y, 1, 1);
		int p;
		for (int k = 0; k < t.nlayers; k++)
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
				{
					coord = l.getPixelCoords(0, t.data[k][i][j]);
					p = FieldOps.round(coord[1]);
					if (p >= 0 && p < ny) histogram[k][p]++;
				}
		
		double epsilon = Double.parseDouble(JOptionPane.showInputDialog("Add epsilon?", "" + 0));
		//Now resize it by sx:
		double[] x = ArrayOps.generateArrayInclBoth(xmin, xmax, t.nlayers*sx);
		double[][] finalHistogram = new double [x.length][ny];
		if (sx == 1)
			for (int i = 0; i < x.length; i++)
				for (int j = 0; j < ny; j++)
					finalHistogram[i][j] = l.evaluateAtMetric(x[i], y[j]) + epsilon;
		else
			for (int i = 0; i < x.length; i++)
				for (int j = 0; j < ny; j++)
					finalHistogram[i][j] = l.data[i][j] + epsilon;
		return new Layer(finalHistogram, x, y, 1, 1);
				
	}
	
	/**
	 * Here the bins array starts from zero and the topomap is an array of layers, one for each bin.
	 * @param sx
	 * @param ny
	 * @param t
	 * @param bins
	 * @return
	 */
	public static Topomap getSpectralDistributionBasicBinned(int sx, int ny, Topomap t, int[][] bins)
	{
		double xmin = ArrayOps.min(t.v);
		double xmax = ArrayOps.max(t.v);
		double ymin = ArrayOps.min(t.data);
		double ymax = ArrayOps.max(t.data);
		
		int nbins = FieldOps.max(bins)+1;
		double[] v = ArrayOps.generateArray(0, 1, nbins);
		double[][][] histogram = new double[nbins][t.nlayers][ny];
		double[] y = ArrayOps.generateArrayInclBoth(ymax, ymin, ny);
		double[] coord;
		Topomap l = new Topomap(histogram, v, t.v, y, null);
		for (int k = 0; k < t.nlayers; k++)
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
				{
					coord = l.getPixelCoords(0, t.data[k][i][j]);
					histogram[bins[i][j]][k][FieldOps.round(coord[1])]++;
				}
		
//		double epsilon = Double.parseDouble(JOptionPane.showInputDialog("Add epsilon?", "" + 0));
		//Now resize it by sx:
//		double[] x = ArrayOps.generateArrayInclBoth(xmin, xmax, t.nlayers*sx);
//		double[][][] finalHistogram = new double [nbins][x.length][ny];
//		if (sx == 1)
//			for (int i = 0; i < x.length; i++)
//				for (int j = 0; j < ny; j++)
//					finalHistogram[bins[i][j]][i][j] = l.getValueAt(x[i], y[j], bins[i][j]);
//		else
//			for (int i = 0; i < x.length; i++)
//				for (int j = 0; j < ny; j++)
//					for (int k = 0; k < nbins; k++)
//						finalHistogram[k][i][j] = l.data[k][i][j];
//		return new Topomap(finalHistogram, v, x, y, null);
		return l;
	}
	public static Layer getSpectralDistributionFancy(int nx, int ny, Topomap t)
	{
		double xmin = ArrayOps.min(t.v);
		double xmax = ArrayOps.max(t.v);
		double ymin = ArrayOps.min(t.data);
		double ymax = ArrayOps.max(t.data);
		
		double[] y = ArrayOps.generateArrayInclBoth(ymax, ymin, ny);
		double[] x = ArrayOps.generateArrayInclBoth(xmin, xmax, nx);
		
		double[] spectrum, spectrumInt;
		spectrumInt = new double [t.nlayers];
		double[] tvInt = new double [t.nlayers];
		double[][] finalHistogram = new double [x.length][ny];
		Layer ans = new Layer(finalHistogram, x, y, 1, 1);
		for (int i = 0; i < t.nlayers; i++)
			tvInt[i] = ans.getPixelCoords(t.v[i], 0)[0];
		
		double distance, gaussCutoff = 1;
		
		for (int i = 0; i < t.nx; i++){System.out.print("\r\n" + i);
			for (int j = 0; j < t.ny; j++)
			{System.out.print(j + " ");
				
				spectrum = t.getSpectrum(i, j);
				for (int p = 0; p < spectrum.length; p++)
					spectrumInt[p] = ans.getPixelCoords(0, spectrum[p])[1];
				for (int m = 0; m < nx; m++)
					for (int n = 0; n < ny; n++)
					{
						distance = Distance.minimumDistanceToFunctionLinear(tvInt, spectrumInt, m, n);
						if (distance < 8*gaussCutoff)
							finalHistogram[m][n] += Math.exp(-(distance*distance)/(2*gaussCutoff*gaussCutoff));
					}
			}
		}
		return ans;
				
	}
	public static Layer getSpectralDistributionFancy(int nx, int ny, Topomap t, int[][] bins, int selectedBin)
	{
		double xmin = ArrayOps.min(t.v);
		double xmax = ArrayOps.max(t.v);
		double ymin = ArrayOps.min(t.data);
		double ymax = ArrayOps.max(t.data);
		
		double[] y = ArrayOps.generateArrayInclBoth(ymax, ymin, ny);
		double[] x = ArrayOps.generateArrayInclBoth(xmin, xmax, nx);
		
		double[] spectrum, spectrumInt;
		spectrumInt = new double [t.nlayers];
		double[] tvInt = new double [t.nlayers];
		double[][] finalHistogram = new double [x.length][ny];
		Layer ans = new Layer(finalHistogram, x, y, 1, 1);
		for (int i = 0; i < t.nlayers; i++)
			tvInt[i] = ans.getPixelCoords(t.v[i], 0)[0];
		
		double distance, gaussCutoff = 1;
		
		for (int i = 0; i < t.nx; i++){System.out.print(" " + i);
			for (int j = 0; j < t.ny; j++)
			{
				if (bins[i][j] == selectedBin){
					spectrum = t.getSpectrum(i, j);
					for (int p = 0; p < spectrum.length; p++)
						spectrumInt[p] = ans.getPixelCoords(0, spectrum[p])[1];
					for (int m = 0; m < nx; m++)
						for (int n = 0; n < ny; n++)
						{
							distance = Distance.minimumDistanceToFunctionLinear(tvInt, spectrumInt, m, n);
							if (distance < 8*gaussCutoff)
								finalHistogram[m][n] += Math.exp(-(distance*distance)/(2*gaussCutoff*gaussCutoff));
						}
				}
			}
		}
		System.out.println();
		return ans;
				
	}
	public static class FourierFilterMethods
	{
		public static boolean[][] getSuppressionField(ArrayList<int[]> points, int nx, int ny)
		{
			boolean[][] ans = new boolean [nx][ny];
			for (int i = 0; i < points.size(); i++)
			{
				ans[(points.get(i)[0]+nx)%nx][(points.get(i)[1]+ny)%ny] = true;
			}
			return ans;
		}
		
		public static boolean[][] getSuppressionFieldCircle(double radius, int nx, int ny)
		{
			int x = 0, y = 0;
			boolean[][] ans = new boolean [nx][ny];
			for (x = (int)(-radius-1); x <= radius; x++)
				for (y = (int)(-radius-1); y <= radius; y++)
					if (x*x + y*y <= radius*radius)
						ans[(x+nx)%nx][(y+ny)%ny] = true;

			return ans;
		}
		public static boolean[][] getSuppressionFieldRing(int rmin, int rmax, int nx, int ny, boolean supressTheRing)
		{
			int x = 0, y = 0;
			boolean[][] ans = new boolean [nx][ny];
			boolean condition;
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
				{
					x = i-rmax;
					y = j-rmax;
					condition = ((x*x + y*y >= rmin*rmin) && (x*x + y*y <= rmax*rmax));
					ans[(x+nx)%nx][(y+ny)%ny] = supressTheRing ? condition : !condition;
				}

			return ans;
		}
		/**
		 * Suppresses everything on the lines x = 0, y = 0, if suppress cross = true. Otherwise suppresses everything else.
		 * This is an alternative to "line slope subtract" in both directions.
		 * @param nx
		 * @param ny
		 * @param suppressCross
		 * @return
		 */
		public static boolean[][] getSupressionFieldCross(int nx, int ny, boolean supressCross)
		{
			boolean[][] ans = new boolean [nx][ny];
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
					ans[i][j] = (i == 0 || j == 0) ? supressCross : !supressCross;
			return ans;
		}
		public static boolean[][] getSupressionFieldLine(int nx, int ny, boolean vertical, boolean suppress)
		{
			boolean[][] ans = new boolean [nx][ny];
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
					ans[i][j] = (vertical ? (i == 0) : (j== 0)) ? suppress : !suppress;
			return ans;
		}
		public static boolean[][] getSupressionFieldStripThroughOrigin(int nx, int ny, double thickness, double[] unitVector, boolean suppress)
		{
			
			boolean[][] ans = new boolean [nx][ny];
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++){
					double dist = Distance.distanceToLineSegment(nx/2, nx/2+unitVector[0], ny/2, ny/2+unitVector[1], i, j, false)[0];
					ans[i][j] = (dist > thickness/2) ? !suppress : suppress;
				}
			boolean[][] realAns = new boolean[nx][ny];
			FieldOps.shift(ans, realAns);
			return realAns;
		}
		public static boolean[][] getSupressionFieldCrossNotO(int nx, int ny, boolean supressCross)
		{
			boolean[][] ans = new boolean [nx][ny];
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
					ans[i][j] = ((i == 0 || j == 0) && (i != 0 || j != 0)) ? supressCross : !supressCross;
			return ans;
		}
		public static boolean[][] getSuppressionFieldNotCircle(double radius, int nx, int ny)
		{
			int x = 0, y = 0;
			boolean[][] ans = new boolean [nx][ny];
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
					ans[i][j] = true;
			
			for (x = -(int)(radius-1); x <= radius; x++)
				for (y = -(int)(radius-1); y <= radius; y++)
					if (x*x + y*y <= radius*radius)
						ans[(x+nx)%nx][(y+ny)%ny] = false;

			return ans;
		}
		
		/**
		 * This draws 1 circle of radius r0 centered at the origin, and n circles of radius r1 centered at the point pts[i] for i < n.
		 * 
		 * @param r0
		 * @param r1
		 * @param pts
		 * @param nx
		 * @param ny
		 * @return
		 */
		public static boolean[][] getSuppressionFieldCircles(int r0, int r1, int[][] pts, int nx, int ny)
		{
			boolean[][] ans = new boolean [nx][ny];
			int x, y; int k;
			boolean inacircle = false;
			for (x = -nx; x < 2*nx; x++)
				for (y = -nx; y < ny; y++)
				{
					inacircle = false;
					inacircle = x*x + y*y < r0*r0;
					for (k = 0; !inacircle && k < pts.length; k++)
						inacircle = inacircle || ((x-pts[k][0])*(x-pts[k][0]) + (y-pts[k][1])*(y-pts[k][1]) < r1*r1);
					ans[(x+nx)%nx][(y+ny)%ny] = true;
				
				}
			return ans;
		}
		public static boolean[][] getSuppressionFieldCircles(int[] r, int[][] pts, int nx, int ny)
		{
			boolean[][] ans = new boolean [nx][ny];
			int x, y; int k;
			boolean inacircle = false;
			for (x = -nx; x < 2*nx; x++){// System.out.println();
				for (y = -ny; y < 2*ny; y++)
//			for (x = 0; x < nx; x++){ System.out.println();
//				for (y = 0; y < ny; y++)
				{
					inacircle = false;
					for (k = 0; !inacircle && k < pts.length; k++)
						inacircle = inacircle || ((x-pts[k][0])*(x-pts[k][0]) + (y-pts[k][1])*(y-pts[k][1]) < r[k]*r[k]);
					if (inacircle){
						ans[(x+nx)%nx][(y+ny)%ny] = true;
//						System.out.println("WTF");
					}

//					System.out.print(inacircle ? "t" : "f");
					
				}
			}
			return ans;
		}
		public static void putSuppressionFieldCircles(int r0, int r1, int[][] pts, int nx, int ny, boolean[][] ans, boolean defalt)
		{
			int x, y; int k;
			boolean inacircle = false;
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
					ans[i][j] = defalt;
			for (x = -nx; x < 2*nx; x++)
				for (y = -nx; y < ny; y++)
				{
//					ans[(x+nx)%nx][(y+ny)%ny] = false;
					inacircle = false;
					inacircle = x*x + y*y < r0*r0;
					for (k = 0; !inacircle && k < pts.length; k++)
						inacircle = inacircle || ((x-pts[k][0])*(x-pts[k][0]) + (y-pts[k][1])*(y-pts[k][1]) < r1*r1);
					if (inacircle) ans[(x+nx)%nx][(y+ny)%ny] = !defalt;
				
				}
		}
		/**
		 * The previous method did not default to all points within r1 of the pts. This one does that, and also not default
		 * to all points OUTSIDE r2 of the points.
		 * @param r0
		 * @param r1
		 * @param r2
		 * @param pts
		 * @param nx
		 * @param ny
		 * @param ans
		 * @param defalt
		 */
		public static void putSuppressionFieldCircles(int r0, int r1, int r2, int[][] pts, int nx, int ny, boolean[][] ans, boolean defalt)
		{
			int x, y; int k;
			boolean inacircle = false;
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
					ans[i][j] = defalt;
			for (x = -nx; x < 2*nx; x++)
				for (y = -nx; y < ny; y++)
				{
//					ans[(x+nx)%nx][(y+ny)%ny] = false;
					inacircle = false;
					inacircle = x*x + y*y < r0*r0;
					for (k = 0; !inacircle && k < pts.length; k++){
						inacircle = inacircle || ((x-pts[k][0])*(x-pts[k][0]) + (y-pts[k][1])*(y-pts[k][1]) < r1*r1);
						inacircle = inacircle || ((x-pts[k][0])*(x-pts[k][0]) + (y-pts[k][1])*(y-pts[k][1]) > r2*r2);
					
					}
					if (inacircle) ans[(x+nx)%nx][(y+ny)%ny] = !defalt;
				
				}
		}
		
		/**
		 * This switches the boolean from default to !default if any of the conditions is met. The conditions is that
		 * the value of scalar[k] should be greater than or less than cutoff[k] according to the greater[k] boolean.
		 * @param scalar
		 * @param cutoff
		 * @param greater
		 * @param ans
		 * @param defalt
		 */
		public static void putSuppressionFieldGeneral(double[][][] scalar, double[] cutoff, boolean[] greater, boolean[][] ans)
		{
			boolean inacircle = false;
			int nx = scalar[0].length;
			int ny = scalar[0][0].length;
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++){
					ans[i][j] = false;
					for (int k = 0; k < scalar.length && !ans[i][j]; k++)
					{
						ans[i][j] = ans[i][j] || greater[k] ? (scalar[k][i][j] > cutoff[k]) : (scalar[k][i][j] < cutoff[k]);
//						switchit = switchit || (scalar[k][i][j] > cutoff[k]);
//						switchit = true;
					}
//					if (switchit){
//						ans[i][j] = !ans[i][j];
//					}
				}
		}
		public static void putSuppressionFieldGeneral(double[][][] scalar, double[] cutoff, boolean[] greater, boolean[][][] ans)
		{
			boolean inacircle = false;
			int nx = scalar[0].length;
			int ny = scalar[0][0].length;
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
					for (int k = 0; k < scalar.length; k++){
						ans[k][i][j] = false;
						if (greater[k] ? (scalar[k][i][j] > cutoff[k]) : (scalar[k][i][j] < cutoff[k]))
							ans[k][i][j] = true;
					}
		}
		/**
		 * This method sets all pixels within R of the 0, 0 corner of the cirlce to 1, else zero.
		 * @param r
		 * @param mask
		 */
		public static void putCircleMaskOrigin(double r, double[][] mask)
		{
			for (int i = 0; i < mask.length; i++)
				for (int j = 0; j < mask[0].length; j++)
				{
					if (Distance.distancePeriodic(i, j, mask.length, mask[0].length) <= r)
						mask[i][j] = 1;
					else mask[i][j] = 0;
				}
		}
		/**
		 * This one puts a radial Fermi function mask. m = 1/(exp((r-r0)/T) + 1)
		 * @param r
		 * @param T
		 * @param mask
		 */
		public static void putFermiMaskOrigin(double r0, double T, double[][] mask)
		{
			double r;
			for (int i = 0; i < mask.length; i++)
				for (int j = 0; j < mask[0].length; j++)
				{
					r = Distance.distancePeriodic(i, j, mask.length, mask[0].length);
					mask[i][j] = 1/(Math.exp((r-r0)/T) + 1);
				}
		}
		public static void putGaussMaskOrigin(double sigma, double[][] mask)
		{
			double r;
			for (int i = 0; i < mask.length; i++)
				for (int j = 0; j < mask[0].length; j++)
				{
					r = Distance.distancePeriodic(i, j, mask.length, mask[0].length);
					mask[i][j] = Math.exp(-(r*r)/(2*sigma*sigma));
				}
		}
		
		/**
		 * This gets a series of masks. the area increases linearly from rmin to rmax
		 * @param rmin
		 * @param rmax
		 * @param n
		 * @return
		 */
		public static double[][][] getCircleMasks(double rmin, double rmax, int n, int N)
		{
			double sqmin = rmin*rmin;
			double sqmax = rmax*rmax;
			double[] num = ArrayOps.generateArrayInclBoth(sqmin, sqmax, n);
			double[][][] mask = new double [n][N][N];
			for (int i = 0; i < n; i++)
				putCircleMaskOrigin(Math.sqrt(num[i]), mask[i]);
			return mask;
		}
		public static double[][][] getFermiMasks(double rmin, double rmax, int n, int N, double T)
		{
			double sqmin = rmin*rmin;
			double sqmax = rmax*rmax;
			double[] num = ArrayOps.generateArrayInclBoth(sqmin, sqmax, n);
			double[][][] mask = new double [n][N][N];
			for (int i = 0; i < n; i++)
				putFermiMaskOrigin(Math.sqrt(num[i]), T, mask[i]);
			return mask;
		}
	}

	/**
	 * This outputs the average spectra of each region in "bins". A "region" consists of all points with the same integer value of bin.
	 * 
	 * @param t
	 * @param bins
	 * @return
	 */
	public static double[][] getAverageSpectraBinned(Topomap t, int[][] bins)
	{
		int min = FieldOps.min(bins);
		int max = FieldOps.max(bins);
		int n = (max-min)+1;
		
		double[][] ans = new double [n][t.nlayers];
		for (int i =0; i < n; i++)
			for (int j = 0; j < t.nlayers; j++)
				ans[i][j] = FieldOps.mean(t.data[j], bins, min+i);
		return ans;
	}
	
	public static Topomap[] splitTopomapBinned(Topomap t, int[][] bins)
	{
		int min = FieldOps.min(bins);
		int max = FieldOps.max(bins);
		int n = (max-min)+1;
		double[][][] temp = new double [n][t.nx][t.ny];
		double[][][][] data = new double [n][t.nlayers][t.nx][t.ny];
		
		for (int i = 0; i < t.nlayers; i++)
		{
			temp = FieldOps.splitByBins(t.data[i], bins);
			for (int j = 0; j < n; j++)
			{
				FieldOps.changeZeroToAverage(temp[j]);
				data[j][i] = FieldOps.copy(temp[j]);
			}
		}

		Topomap[] ans = new Topomap[n];
		for (int i = 0; i < n; i++)
			ans[i] = Topomap.newTopomap(t, data[i]);
		return ans;
	}
	public static Topomap splitTopomapBinned(Topomap t, int[][] bins, int binIndex)
	{
		int min = FieldOps.min(bins);
		int max = FieldOps.max(bins);
		double[][] temp = new double [t.nx][t.ny];
		double[][][] data = new double [t.nlayers][t.nx][t.ny];
		
		for (int i = 0; i < t.nlayers; i++)
		{
			temp = FieldOps.splitByBins(t.data[i], bins)[binIndex];
			FieldOps.changeZeroToAverage(temp);
			data[i] = FieldOps.copy(temp);
		}

		Topomap ans = Topomap.newTopomap(t, data);
		return ans;
	}
	public static Topomap splitTopomapBinned(Topomap t, int[][] bins, int binIndex, boolean[][] use)
	{
		int min = FieldOps.min(bins);
		int max = FieldOps.max(bins);
		double[][] temp = new double [t.nx][t.ny];
		double[][][] data = new double [t.nlayers][t.nx][t.ny];
		
		for (int i = 0; i < t.nlayers; i++)
		{
			temp = FieldOps.splitByBins(t.data[i], bins)[binIndex];
			FieldOps.changeZeroToAverage(temp);
			data[i] = FieldOps.copy(temp);
		}

		Topomap ans = Topomap.newTopomap(t, data);
		return ans;
	}
	
	public static void subtractPolynomialFitEachSpectrum(Topomap t, int n)
	{
		double[] spec, subfit;
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.nx; j++)
			{
				spec = t.getSpectrum(i, j);
				subfit = ArrayOps.subtractPolynomialFit(t.v, spec, n)[0];
				for (int k = 0; k < t.nlayers; k++)
					t.data[k][i][j] = subfit[k];
			}
	}
	
	public static double[] getAverageFFTOfSpectra(Topomap t)
	{
		double[] ans = new double [t.nlayers];
		double[] temp = new double [t.nlayers];
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
			{
				temp = FFTOps.get1DFFTMag(t.getSpectrum(i, j));
				for (int k = 0; k < t.nlayers; k++)
					ans[k] += temp[k]/(t.nx*t.ny);
			}
		return ans;
	}
	public static void replaceWithFFTOfSpectra(Topomap t)
	{
		boolean log = JOptionPane.showConfirmDialog(null, "Take the log?") == JOptionPane.YES_OPTION;

		double[] ans = new double [t.nlayers];
		double[] temp = new double [t.nlayers];
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
			{
				temp = FFTOps.get1DFFTMag(t.getSpectrum(i, j));
				if (log) FieldOps.log(temp);
				for (int k = 0; k < t.nlayers; k++)
					t.data[k][i][j] = temp[k];
			}
	}
	/**
	 * Replaces lockin-generated zero values with the maximum value in each spectrum, to reduce how ugly the image looks.
	 * @param t
	 * @param vcutoff
	 * @param lessThan
	 * @param amplitudecutoff
	 */
	public static void removeLockinZeroesAtExtremes(Topomap t, double vcutoff, boolean lessThan, double amplitudecutoff)
	{
		
		double[][] max = new double [t.nx][t.ny];
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
				max[i][j] = ArrayOps.max(t.getSpectrum(i, j));
		
		if (lessThan)
			for (int k = 0; t.v[k] <  vcutoff; k++){System.out.print(" " + k);
				for (int i = 0; i < t.nx; i++)
					for (int j = 0; j < t.ny; j++)
					{
						if (t.data[k][i][j] < amplitudecutoff)
							t.data[k][i][j] = max[i][j];
					}
			}
		else
		for (int k = t.nlayers-1; t.v[k] > vcutoff; k--){System.out.print(" " + k);
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
				{
					if (t.data[k][i][j] < amplitudecutoff)
						t.data[k][i][j] = max[i][j];
				}
		}

	}
	public static void writeAverageSpectraAroundImps(Topomap t, double gaussradius, JFileChooser fc)
	{
		PointImp[] imps = PointImp.readFromGaussSquareFile(FileOps.selectOpen(fc));
		int nspec = imps.length;
		double[][] averageSpec = new double [nspec][t.nlayers];
		double[][] correlations = new double [nspec][t.nlayers];
		double[][] map = new double [t.nx][t.ny];
		double[][] pixelWeight = ImpurityUtil.getPixelWeight(t.nx, t.ny, gaussradius);
		for (int i = 0; i < imps.length; i++)
		{
			System.out.println(i);
			map = ImpurityUtil.makeImpurityGaussMaskNormal(imps, t.nx, t.ny, gaussradius, pixelWeight, i);
			for (int j = 0; j < t.nlayers; j++)
			{
				averageSpec[i][j] = FieldOps.innerProduct(map, t.data[j]);
				correlations[i][j] = FieldOps.correlation(map, t.data[j]);
			}
		}
		
		double[] x = new double [nspec];
		double[] y = new double [nspec];
		for (int i = 0; i < nspec; i++)
		{
			x[i] = imps[i].pixelPos[0];
			y[i] = imps[i].pixelPos[1];
		}
		PointSpectra average = new PointSpectra(averageSpec, t.v, x, y);
		PointSpectra corr = new PointSpectra(correlations, t.v, x, y);
		
		File f = FileOps.selectSave(fc);
		
		PointSpectra.writeBIN(average, f.toString() + "average.bin");
		PointSpectra.writeBIN(corr, f.toString() + "corr.bin");
		new SpectraDrawer(average, null);
		new SpectraDrawer(corr, null);
		
		ColumnIO.writeLines(average.toTextTable(), f.toString() + "average.txt");
		ColumnIO.writeLines(corr.toTextTable(), f.toString() + "corr.txt");
	}
	
	public static Layer sumOneDimension(Topomap t, int dimension)
	{
		if (dimension == 0)
		{
			double[][] data = new double [t.nx][t.ny];
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
					for (int k = 0; k < t.nlayers; k++)
					{
						data[i][j] += t.data[k][i][j];
					}
			return new Layer(data, t.x, t.y, 1, 1);
		}
		if (dimension == 1)
		{
			double[][] data = new double [t.nlayers][t.ny];
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
					for (int k = 0; k < t.nlayers; k++)
					{
						data[k][j] += t.data[k][i][j];
					}
			return new Layer(data, t.v, t.y, 1, 1);
		}
		else //if (dimension == 2)
		{
			double[][] data = new double [t.nx][t.nlayers];
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
					for (int k = 0; k < t.nlayers; k++)
					{
						data[i][k] += t.data[k][i][j];
					}
			return new Layer(data, t.x, t.v, 1, 1);
		}
	}
	
	/**
	 * The box array consitutes the area of the summed-over-Y fft which we will suppress in the 
	 * fft of the topography before IFFT-ing it
	 * @param t
	 * @param boxz
	 * 
	 * @return
	 */
	public static Topomap fourierFilterTheMap(Topomap t, ArrayList<FilterBox> box)
	{
		FFT3D_Wrapper fft = new FFT3D_Wrapper(FieldOps.copy(t.data));
		fft.doFFT();
		double minMag = FieldOps.min(t.data);
		double mag = 0;
		int kstar, istar;
		int itemp;
		for (int p = 0; p < box.size(); p++)
		{
			for (int k = box.get(p).jmin; k < box.get(p).jmin + box.get(p).height; k++) //k is the layer index
				for (int i = box.get(p).imin; i < box.get(p).imin + box.get(p).width; i++){
					itemp = (i+t.nx/2)%t.nx;
					kstar = t.nlayers-k;
					istar = t.nx-itemp;//This assumes the map has an even number of pixels horizontally, which should always be the case
					if (istar == t.nx) istar = 0;
					for (int j = 0; j < t.ny; j++)
					{
						mag = Complex.mag(fft.output[k][itemp][j]);
						fft.output[k][itemp][j][0] *= minMag/(mag);
						fft.output[k][itemp][j][1] *= minMag/(mag);
//						kstar = (k-(t.nlayers/2)+t.nlayers)%t.nlayers;
//						istar = (i+t.nx/2)%t.nx;
						if (istar != itemp || kstar != k){
							fft.output[kstar][istar][j][0] *= minMag/(mag);
							fft.output[kstar][istar][j][1] *= minMag/(mag);
						}
						System.out.println("" + i + "\t" + istar + "\t" + ((i+t.nx/2)%t.nx) + "\t" + k + "\t" + kstar);
					}
				}
		}
		
		fft.doIFFT();
		
		return Topomap.newTopomap(t, fft.data);
	}
	
	public static class FilterBox
	{
		public FilterBox(int imin, int jmin, int width, int height) {
			super();
			this.imin = imin;
			this.jmin = jmin;
			this.width = width;
			this.height = height;
		}
		int imin, jmin;
		int width, height;
	}
	
	public static void doFourierFilterUserSpecified(Topomap t)
	{
		ArrayList<FilterBox> box = new ArrayList<FilterBox>();
		box.add(new FilterBox(148,13,34,3));
		Topomap filt = fourierFilterTheMap(t, box);
		Topomap.writeBIN(filt, fc);
	}

	public static Topomap getFromConvolutionSeriesTest(Layer t, AtomicCoordinatesSet latt)
	{
		int nl = 64;
		double[][][] mask = new double [nl][t.nx][t.ny];
		double[] length = new double [nl];
		for (int i = 0; i < nl; i++)
		{
//			mask[i][i][i] = 1;
			length[i] = (i+1.0)/4;
			//mask[i] = FieldOps.getGaussSmoothConvolve(t.data, length[i]);
			mask[i] = FieldOps.getDWaveConvolve_sharper(t.data, length[i], latt);
		}
		
		double[][][] ans = new double [nl][][];
		for (int i = 0; i < nl; i++)
			ans[i] = FieldOps.convolveFourier(t.data, mask[i]);
		
		new TopomapViewer(new Topomap(mask, ArrayOps.generateArrayNotInclUpper(0, nl, nl), t.x, t.y, null), "", 512);
		
		return new Topomap(ans, ArrayOps.generateArrayNotInclUpper(0, nl, nl), t.x, t.y, null);
		
	}

	/**
	 * 
	 * @param t
	 * @param direction
	 * @return
	 */
	public static PointSpectra getTimeOrderedPointSpectra(Topomap t, int direction) {
		double[] x = new double [t.nx*t.ny], y = new double [t.nx*t.ny];
		double[][] data = new double [t.nx*t.ny][t.nlayers];
		if (direction == 0)
		{
			int n = 0;
			for (int j = 0; j < t.ny; j++)
				for (int i = 0; i < t.nx; i++){
					for (int k = 0; k < t.nlayers; k++) data[n][k] = t.data[k][i][j];
					x[n] = t.x[i];
					y[n] = t.y[j];
					n++;
				}
		}
		return new PointSpectra(data, t.v, x, y);
	}
	
	public static class MapCuttingMethods
	{
		//These methods all have the same basic form. They average all pixels in the map, which is usually
		//a symmetrized FFT, using a boolean filter. The boolean filter is different for each energy. The number of true values must be at least one for this to work.
		/**
		 * This computes the map using a radial average about the Rt2 points. The Rt2 poiunts
		 * are computed assuming latt is the 1x1 lattice in real space: we take latt.getRt2() then 
		 * get the reciprocal and round to the nearest int. Then we generate a distance map
		 * and look through it with various sizes to get the boolean array.
		 * 
		 * @param t - a map of the same dimensiont 
		 * @param latt - the 1x1 lattice in real space
		 * @param rt2pt - 0 for one, 1 for the other 
		 * @return
		 */
		public static double[][] getRadialRt2Point(int nx, int ny, AtomicCoordinatesSet latt, int rt2pt)
		{
			double[][] distance = new double [nx][ny];
//			AtomicCoordinatesSet rt2 = latt.getRt2Lattice();
			int cx, cy;
			int[] rt2 = getRt2Point(latt, rt2pt, nx);
			cx = rt2[0]; cy = rt2[1];
			System.out.println("" + cx + "\t" + cy);
			int x, y;
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
				{
					x = i-nx/2;
					y = j-ny/2;
					distance[i][j] = Distance.distance(x-cx, y-cy);
				}
			return distance;
		}
		public static double[][] getRadialRt2Point(Topomap t, AtomicCoordinatesSet latt, int rt2pt){
			return getRadialRt2Point(t.nx, t.ny, latt, rt2pt);
		}
		public static double[][] getRadialRt2Point(Layer t, AtomicCoordinatesSet latt, int rt2pt){
			return getRadialRt2Point(t.nx, t.ny, latt, rt2pt);
		}
		public static double[][] getRadialRt2Point(double[][] t, AtomicCoordinatesSet latt, int rt2pt){
			return getRadialRt2Point(t.length, t[0].length, latt, rt2pt);
		}
		
		public static double[][] getUnrestrictedDotWithDirection(int nx, int ny, AtomicCoordinatesSet latt, int direction)
		{
			double[][] distance = new double [nx][ny];
//			AtomicCoordinatesSet rt2 = latt.getRt2Lattice();
			double[] aHat = Distance.unitVector(direction == 0 ? latt.getA() : latt.getB());
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
				{
					distance[i][j] = FieldOps.dot(aHat, i-nx/2, j-ny/2);
				}
			return distance;
		}
		public static double[][] getRadialOrigin(Topomap t)
		{
			return getRadialOrigin(t.nx, t.ny);
		}
		public static double[][] getRadialOrigin(int nx, int ny)
		{
			double[][] distance = new double [nx][ny];
			int x, y;
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
				{
					x = i-nx/2;
					y = j-ny/2;
					distance[i][j] = Distance.distance(x, y);
				}
			return distance;
		}
		/**
		 * This computes a "distance" map where the equi-distance surfaces are wedges with right angled leading edges.
		 * @param t
		 * @param latt
		 * @param rt2pt
		 * @return
		 */
		public static double[][] getInwardRightAngleRt2Point(Topomap t, AtomicCoordinatesSet latt, int rt2pt)
		{
			double[][] distance = new double [t.nx][t.ny];
			int cx, cy;
			int[] rt2 = getRt2Point(latt, rt2pt, t);
			cx = rt2[0]; cy = rt2[1];
			double[] unitVectorToOrigin = Distance.unitVector(new double[] {-cx, -cy});
			double[] uVRight = Matrix.getProductWith(Matrix.getRotationMatrix(-Math.PI/4), unitVectorToOrigin);
			double[] uVLeft = Matrix.getProductWith(Matrix.getRotationMatrix(Math.PI/4), unitVectorToOrigin);
			System.out.println("" + cx + "\t" + cy);
			System.out.println(Printer.vectorP(unitVectorToOrigin));
			System.out.println(Printer.vectorP(uVLeft));
			System.out.println(Printer.vectorP(uVRight));
			double x, y;
			double xr, yr;
			double in;
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
				{
					x = i-t.nx/2;
					y = j-t.ny/2;
					xr = x-cx;
					yr = y-cy;
					distance[i][j] = Math.max(Distance.dot(xr, yr, uVRight[0], uVRight[1]), Distance.dot(xr, yr, uVLeft[0], uVLeft[1]));
					in = Distance.dot(xr, yr, unitVectorToOrigin[0], unitVectorToOrigin[1]);
					if (in < 0)
						distance[i][j] = 2*in;
				}
			return distance;
		}
		public static double[][] getStripInwardFromBraggPeak(Topomap t, AtomicCoordinatesSet latt, int bragg, int width)
		{
			double[][] distance = new double [t.nx][t.ny];
			int cx, cy;
			int[][] braggi = AtomicCoordinatesSet.generateBragg(latt, t.nx);
			cx = braggi[bragg][0]; cy = braggi[bragg][1];
			double[] unitVectorToOrigin = Distance.unitVector(new double[] {-cx, -cy});
			
			double[] uVRight = Matrix.getProductWith(Matrix.getRotationMatrix(-Math.PI/2), unitVectorToOrigin);
			System.out.println("" + cx + "\t" + cy);
			System.out.println(Printer.vectorP(unitVectorToOrigin));
			System.out.println(Printer.vectorP(uVRight));
			double x, y;
			double xr, yr;
			double in;
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
				{
					x = i-t.nx/2;
					y = j-t.ny/2;
					xr = x-cx;
					yr = y-cy;
					distance[i][j] = Distance.dot(xr, yr, unitVectorToOrigin[0], unitVectorToOrigin[1]);
					in = Math.abs(Distance.dot(xr, yr, uVRight[0], uVRight[1]));
					if (in > width/2)
						distance[i][j] = -1;
				}
			return distance;
			
		}
		public static double[][] getStripOutwardTowardsBraggPeak(Topomap t, AtomicCoordinatesSet latt, int bragg, int width)
		{
			double[][] distance = new double [t.nx][t.ny];
			int cx, cy;
			int[][] braggi = AtomicCoordinatesSet.generateBragg(latt, t.nx);
			cx = 0; cy = 0;
			double[] unitVectorToBragg = Distance.unitVector(new double[] {braggi[bragg][0], braggi[bragg][1]});
			
			double[] uVRight = Matrix.getProductWith(Matrix.getRotationMatrix(-Math.PI/2), unitVectorToBragg);
			System.out.println("" + cx + "\t" + cy);
			System.out.println(Printer.vectorP(unitVectorToBragg));
			System.out.println(Printer.vectorP(uVRight));
			double x, y;
			double xr, yr;
			double in;
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
				{
					x = i-t.nx/2;
					y = j-t.ny/2;
					xr = x-cx;
					yr = y-cy;
					distance[i][j] = Distance.dot(xr, yr, unitVectorToBragg[0], unitVectorToBragg[1]);
					in = Math.abs(Distance.dot(xr, yr, uVRight[0], uVRight[1]));
					if (in > width/2)
						distance[i][j] = -1;
				}
			return distance;
		}
		public static double[][] getStripOutwardTowardsRt2(Topomap t, AtomicCoordinatesSet latt, int bragg, int width)
		{
			double[][] distance = new double [t.nx][t.ny];
			int cx, cy;
			int[] rt2 = getRt2Point(latt, 0, t);
			cx = 0; cy = 0;
			double[] unitVectorToBragg = Distance.unitVector(new double[] {rt2[0], rt2[1]});
			
			double[] uVRight = Matrix.getProductWith(Matrix.getRotationMatrix(-Math.PI/2), unitVectorToBragg);
			System.out.println("" + cx + "\t" + cy);
			System.out.println(Printer.vectorP(unitVectorToBragg));
			System.out.println(Printer.vectorP(uVRight));
			double x, y;
			double xr, yr;
			double in;
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
				{
					x = i-t.nx/2;
					y = j-t.ny/2;
					xr = x-cx;
					yr = y-cy;
					distance[i][j] = Distance.dot(xr, yr, unitVectorToBragg[0], unitVectorToBragg[1]);
					in = Math.abs(Distance.dot(xr, yr, uVRight[0], uVRight[1]));
					if (in > width/2)
						distance[i][j] = -1;
				}
			return distance;
		}
		public static int[] getRt2Point(AtomicCoordinatesSet latt, int rt2pt, Topomap t)
		{
			return getRt2Point(latt, rt2pt, t.nx);
		}
		public static int[] getRt2Point(AtomicCoordinatesSet latt, int rt2pt, int nx)
		{
			int[][] centers = AtomicCoordinatesSet.generateBragg(latt, nx);
			int cx, cy;
			if (rt2pt == 0) {
				cx = (centers[0][0] + centers[1][0])/2;
				cy = (centers[0][1] + centers[1][1])/2;;
			}
			else
			{
				cx = (centers[0][0] - centers[1][0])/2;
				cy = (centers[0][1] - centers[1][1])/2;;
			}
			return new int [] {cx, cy};
		}
//		public static double[][] getRightAngledRt2Points(Topomap t, AtomicCoordinatesSet latt, int rt2pt)
//		{
//			
//		}
		public static boolean[][][] getRadialRt2Masks(Topomap t, AtomicCoordinatesSet latt, int rt2pt, double maxLength, int npts)
		{
			return getMasks(getRadialRt2Point(t, latt, rt2pt), npts, maxLength);
		}
		
		public static void maskMap(Topomap t, boolean[][][] filter, double[] xAxis, String dir, String ending, String cutname, boolean askSave)
		{
			String outdir = dir + ending + "_" + cutname + "_cut\\";// + (normalizeEachLine ? "norm_layer" : "") + "\\";
			File outd = new File (outdir);
			if (!outd.exists()) outd.mkdir();
			
						
			int npts = filter.length;
//			double[][] onData = stripCutEachEnergy(t, width, unitBragg.getA());
			double[][] result = new double [npts][t.nlayers];
			for (int i = 0; i < npts; i++)
				result[i] = getSpectrumAveragedOverBooleans(t, filter[i]);
			
			double[][] resultTranspose = FieldOps.transpose(result);
			//obtain E vs K. For K, we remember that the dataset is the symmetrized FFT so that it's real space is actually k-space.
			
			PointSpectra on = new PointSpectra(resultTranspose, xAxis, new double[t.nlayers], t.v);
			PointSpectra.writeBIN(on, outdir + "Intensity spectra I vs k spec.bin");
			Layer onl = on.toLayer();
			onl.flipY();
			Layer.writeBIN(onl, outdir + "Intensity spectra I vs k layer.bin");
			if (askSave)
				if (JOptionPane.showConfirmDialog(null, "Write the masks used?") == JOptionPane.YES_OPTION)
					Topomap.writeBIN(new Topomap(ArrayOps.toDouble(filter), xAxis, t.x, t.y, null), outdir + "masks.bin");
			//
		}
		public static Layer getCut(Topomap t, boolean[][][] filter, double[] xAxis, String dir, String ending, String cutname, boolean askSave)
		{
			int npts = filter.length;
//			double[][] onData = stripCutEachEnergy(t, width, unitBragg.getA());
			double[][] result = new double [npts][t.nlayers];
			for (int i = 0; i < npts; i++)
				result[i] = getSpectrumAveragedOverBooleans(t, filter[i]);
			
			double[][] resultTranspose = FieldOps.transpose(result);
			//obtain E vs K. For K, we remember that the dataset is the symmetrized FFT so that it's real space is actually k-space.
			
			PointSpectra on = new PointSpectra(resultTranspose, xAxis, new double[t.nlayers], t.v);
			Layer onl = on.toLayer();
//			onl.flipY();
//			Layer.writeBIN(onl, outdir + "Intensity spectra I vs k layer.bin");
//			if (askSave)
//				if (JOptionPane.showConfirmDialog(null, "Write the masks used?") == JOptionPane.YES_OPTION)
//					Topomap.writeBIN(new Topomap(ArrayOps.toDouble(filter), xAxis, t.x, t.y, null), outdir + "masks.bin");
			return onl;
			//
		}
		public static Layer getRadialCut(Topomap t)
		{
			boolean[][][] filter = getMasksThickness(getRadialOrigin(t), 1);
			int npts = filter.length;
			double dx = Math.abs(t.xLength/(t.nx-1));
			double[] xAxis = ArrayOps.generateArray(0, dx, npts);
//			double[][] onData = stripCutEachEnergy(t, width, unitBragg.getA());
			double[][] result = new double [npts][t.nlayers];
			for (int i = 0; i < npts; i++)
				result[i] = getSpectrumAveragedOverBooleans(t, filter[i]);
			
			double[][] resultTranspose = FieldOps.transpose(result);
			//obtain E vs K. For K, we remember that the dataset is the symmetrized FFT so that it's real space is actually k-space.
			
			PointSpectra on = new PointSpectra(resultTranspose, xAxis, new double[t.nlayers], t.v);
			Layer onl = on.toLayer();
			onl.flipY();
//			Layer.writeBIN(onl, outdir + "Intensity spectra I vs k layer.bin");
//			if (askSave)
//				if (JOptionPane.showConfirmDialog(null, "Write the masks used?") == JOptionPane.YES_OPTION)
//					Topomap.writeBIN(new Topomap(ArrayOps.toDouble(filter), xAxis, t.x, t.y, null), outdir + "masks.bin");
			return onl;
			//
		}
		public static boolean[][][] getMasks(double[][] distance, int npts, double maxLength){
			int nx = distance.length, ny = distance[0].length;
			boolean[][][] b = new boolean [npts][nx][ny];
			double dx = maxLength/npts;
			double[] length = ArrayOps.generateArrayInclBoth(0, maxLength, npts+1);
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
					for (int k = 0; k < npts; k++)
					{
						b[k][i][j] = distance[i][j] >= length[k] && distance[i][j] < length[k+1];
					}
			return b;
		}
		public static boolean[][][] getMasksThickness(double[][] distance, double thickness){
			int nx = distance.length, ny = distance[0].length;
			double minLength = ArrayOps.min(distance);
			double maxLength = ArrayOps.max(distance);
			int npts = (int)((maxLength-minLength)/thickness);
			boolean[][][] b = new boolean [npts][nx][ny];
			double[] length = ArrayOps.generateArrayInclBoth(minLength, maxLength, npts+1);
			for (int i = 0; i < nx; i++){//System.out.println(i);
				for (int j = 0; j < ny; j++)
					for (int k = 0; k < npts; k++)
					{
						b[k][i][j] = distance[i][j] >= length[k] && distance[i][j] < length[k+1];
					}
			}
			return b;
		}
		public static double[] getSpectrumAveragedOverBooleans(Topomap t, boolean[][] map)
		{
			double[] ans = new double [t.nlayers];
			for (int i = 0; i < t.nlayers; i++)
				ans[i] = FieldOps.mean(t.data[i], map);
			return ans;
		}
		/**
		 * 
		 * @param t - an unsymmetrized FFT
		 * @param latt - real space
		 * @param width - the width of the strip
		 * @param depth - the thickness of the strip. 2 means each pixel of the result will represent 2 pixels thick of the map.
		 */
		public static void writeCutsToBothBraggPeaks(Topomap t, AtomicCoordinatesSet latt, int width, int depth, String dir, String name, boolean isFourierTransform, boolean askSave)
		{
			int[][] braggi = AtomicCoordinatesSet.generateBragg(latt, t.nx);
			double maxLength = Distance.distance(braggi[0][0], braggi[0][1])*1.05;
			maxLength = (((int)(maxLength/depth))+1)*depth;
			int npts = (int)maxLength/depth;
			double[][][] fftmag = null;
			if (!isFourierTransform){
				fftmag = new double [t.nlayers][t.nx][t.ny];
				for (int i = 0; i < t.nlayers; i++)
					fftmag[i] = FFTOps.obtainFFTmagCent(t.data[i]);
			}
			Topomap use = isFourierTransform ? t : Topomap.newTopomap(t, fftmag);
			for (int braggIndex = 0; braggIndex < 2; braggIndex++){
				double[][] distance = MapCuttingMethods.getStripOutwardTowardsBraggPeak(t, latt, braggIndex, width);//MapCuttingMethods.getInwardRightAngleRt2Point(t, latt, 0);
		//		LayerViewer.show(Layer.getFreeLayer(distance), 1024, true);
				boolean[][][] map = MapCuttingMethods.getMasks(distance, npts, maxLength);//MapCuttingMethods.getRadialRt2Masks(t, latt, 0, maxLength, npts);
				MapCuttingMethods.maskMap(use, map, ArrayOps.generateArrayNotInclLower(0, maxLength*(t.xLength/t.nx), npts), dir, name, "bragg" + (braggIndex+1) + "_w" + width + "_d"+depth, askSave);
			}
		}
		public static boolean[][][][] getMapsToBothBraggPeaks(Topomap t, AtomicCoordinatesSet latt, int width, int depth)
		{
			int[][] braggi = AtomicCoordinatesSet.generateBragg(latt, t.nx);
			double maxLength = Distance.distance(braggi[0][0], braggi[0][1])*1.05;
			maxLength = (((int)(maxLength/depth))+1)*depth;
			int npts = (int)maxLength/depth;
			boolean[][][][] ans = new boolean [2][][][];
			for (int braggIndex = 0; braggIndex < 2; braggIndex++){
				double[][] distance = MapCuttingMethods.getStripOutwardTowardsBraggPeak(t, latt, braggIndex, width);//MapCuttingMethods.getInwardRightAngleRt2Point(t, latt, 0);
				ans[braggIndex] = MapCuttingMethods.getMasks(distance, npts, maxLength);//MapCuttingMethods.getRadialRt2Masks(t, latt, 0, maxLength, npts);
			}
			return ans;
		}
		/**
		 * This takes a presumably symmetrized FFT and writes cuts along the Bragg and Rt2 directions.
		 * @param t
		 * @param latt
		 * @param width
		 * @param depth
		 * @param dir
		 * @param name
		 */
		public static void writeHighSymmetryCuts(Topomap t, AtomicCoordinatesSet latt, int width, int depth, String dir, String name, boolean askSave)
		{
			int[][] braggi = AtomicCoordinatesSet.generateBragg(latt, t.nx);
			double maxLength = Distance.distance(braggi[0][0], braggi[0][1]);
			maxLength = (((int)(maxLength/depth))+1)*depth;
			int npts = (int)maxLength/depth;
			double[][] distance = MapCuttingMethods.getStripOutwardTowardsBraggPeak(t, latt, 0, width);//MapCuttingMethods.getInwardRightAngleRt2Point(t, latt, 0);
			double[][] distanceRt2 = MapCuttingMethods.getStripOutwardTowardsRt2(t, latt, 0, width);//MapCuttingMethods.getInwardRightAngleRt2Point(t, latt, 0);
//			LayerViewer.show(Layer.getFreeLayer(distance), 1024, true);
			boolean[][][] map = MapCuttingMethods.getMasks(distance, npts, maxLength);//MapCuttingMethods.getRadialRt2Masks(t, latt, 0, maxLength, npts);
			MapCuttingMethods.maskMap(t, map, ArrayOps.generateArrayNotInclLower(0, maxLength*(t.xLength/t.nx), npts), dir, name, "bragg" + "_w" + width + "_d"+depth, askSave);
			boolean[][][] maprt2 = MapCuttingMethods.getMasks(distanceRt2, npts, maxLength);//MapCuttingMethods.getRadialRt2Masks(t, latt, 0, maxLength, npts);
			MapCuttingMethods.maskMap(t, maprt2, ArrayOps.generateArrayNotInclLower(0, maxLength*(t.xLength/t.nx), npts), dir, name, "rt2" + "_w" + width + "_d"+depth, askSave);
		}
		public static void writeRadialCutsOrigin(Topomap t, int depth, String dir, String name, boolean askSave)
		{
			double maxLength = Distance.distance(t.nx/2, t.ny/2);
			maxLength = (((int)(maxLength/depth))+1)*depth;
			int npts = (int)maxLength/depth;
			double[][] distance = MapCuttingMethods.getRadialOrigin(t);//MapCuttingMethods.getInwardRightAngleRt2Point(t, latt, 0);
//			LayerViewer.show(Layer.getFreeLayer(distance), 1024, true);
			boolean[][][] map = MapCuttingMethods.getMasks(distance, npts, maxLength);//MapCuttingMethods.getRadialRt2Masks(t, latt, 0, maxLength, npts);
			MapCuttingMethods.maskMap(t, map, ArrayOps.generateArrayNotInclLower(0, maxLength*(t.xLength/t.nx), npts), dir, name, "radial" + "_d"+depth, askSave);
		}
		public static void writeRadialCutRt2(Topomap t, AtomicCoordinatesSet latt, double maxLength, int depth, String dir, String name, int rt2pt, boolean askSave)
		{
			int[][] braggi = AtomicCoordinatesSet.generateBragg(latt, t.nx);
//			double maxLength = Distance.distance(braggi[0][0], braggi[0][1]);
			maxLength = (((int)(maxLength/depth))+1)*depth;
			int npts = (int)maxLength/depth;
			double[][] distance = MapCuttingMethods.getRadialRt2Point(t, latt, rt2pt);
//			LayerViewer.show(Layer.getFreeLayer(distance), 1024, true);
			boolean[][][] map = MapCuttingMethods.getMasks(distance, npts, maxLength);//MapCuttingMethods.getRadialRt2Masks(t, latt, 0, maxLength, npts);
			MapCuttingMethods.maskMap(t, map, ArrayOps.generateArrayNotInclLower(0, maxLength*(t.xLength/t.nx), npts), dir, name, "radialRt2_" + rt2pt + "_d"+depth, askSave);
		}
		public static void writeWedgeCutsRt2(Topomap t, AtomicCoordinatesSet latt, double maxLength, int depth, String dir, String name, int rt2pt, boolean askSave)
		{
			int[][] braggi = AtomicCoordinatesSet.generateBragg(latt, t.nx);
//			double maxLength = Distance.distance(braggi[0][0], braggi[0][1]);
			maxLength = (((int)(maxLength/depth))+1)*depth;
			int npts = (int)maxLength/depth;
			double[][] distance = MapCuttingMethods.getInwardRightAngleRt2Point(t, latt, rt2pt);
//			LayerViewer.show(Layer.getFreeLayer(distance), 1024, true);
			boolean[][][] map = MapCuttingMethods.getMasks(distance, npts, maxLength);//MapCuttingMethods.getRadialRt2Masks(t, latt, 0, maxLength, npts);
			MapCuttingMethods.maskMap(t, map, ArrayOps.generateArrayNotInclLower(0, maxLength*(t.xLength/t.nx), npts), dir, name, "wedgeRt2_" + rt2pt + "_d"+depth, askSave);
		}
		
		public static double[][] getFFTCentroidEachLayer(Topomap t, boolean[][] mask){
			double[][] ans = new double [t.nlayers][2];
			for (int i = 0; i < t.nlayers; i++)
			{
				double[][] fftmag = FFTOps.obtainFFTmagCent(t.data[i]);
				ans[i] = FieldOps.centroid(fftmag, mask);
			}
			return ans;
		}

	}
	public static void subtractPolynomialFitFrom2ndDerivativeFolder(String dir, String base, int degree){
		String[] dudr = {"uaa.bin","uab.bin","uba.bin","ubb.bin"};
		String[] sums = {"tr_e.bin","e12.bin","w12.bin","uanis.bin"};
		String outdir = dir + "sub_" + degree + "_poly\\";
		FileOps.mkdir(outdir);
		Topomap uaai = Topomap.readBIN(dir + base + dudr[0]);
		Topomap uabi = Topomap.readBIN(dir + base + dudr[1]);
		Topomap ubai = Topomap.readBIN(dir + base + dudr[2]);
		Topomap ubbi = Topomap.readBIN(dir + base + dudr[3]);
		int nl = uaai.nlayers, nx = uaai.nx, ny = uaai.ny; 
		double[][][] tre = new double [nl][][];
		for (int i  = 0; i < uaai.nlayers; i++){
			FieldOps.subtractNSheetFit(uaai.data[i], degree);
			FieldOps.subtractNSheetFit(uabi.data[i], degree);
			FieldOps.subtractNSheetFit(ubai.data[i], degree);
			FieldOps.subtractNSheetFit(ubbi.data[i], degree);
		}
		Topomap.writeBIN(uaai, outdir + base + dudr[0]);
		Topomap.writeBIN(uabi, outdir + base + dudr[1]);
		Topomap.writeBIN(ubai, outdir + base + dudr[2]);
		Topomap.writeBIN(ubbi, outdir + base + dudr[3]);
		double[][][] e12 = new double [nl][nx][ny];
		double[][][] w12 = new double [nl][nx][ny];
		double[][][] tr_e = new double [nl][nx][ny];
		double[][][] anis = new double [nl][nx][ny];
		for (int k = 0; k < nl; k++)
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
				{
					e12[k][i][j] = (uabi.data[k][i][j] + ubai.data[k][i][j])/2;
					w12[k][i][j] = (uabi.data[k][i][j] - uabi.data[k][i][j])/2;
					tr_e[k][i][j] = (uaai.data[k][i][j] + ubbi.data[k][i][j])/2;
					anis[k][i][j] = (uaai.data[k][i][j] - ubbi.data[k][i][j])/2;
				}
		Topomap.writeBIN(Topomap.newTopomap(uaai, tr_e), outdir + base + sums[0]);
		Topomap.writeBIN(Topomap.newTopomap(uaai, e12), outdir + base + sums[1]);
		Topomap.writeBIN(Topomap.newTopomap(uaai, w12), outdir + base + sums[2]);
		Topomap.writeBIN(Topomap.newTopomap(uaai, anis), outdir + base + sums[3]);
	}

}

