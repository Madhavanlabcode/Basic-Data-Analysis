package util;

import image.ImageEditing;

import java.awt.Graphics;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;

import org.apache.commons.math3.analysis.interpolation.BicubicSplineInterpolatingFunction;
import org.apache.commons.math3.analysis.interpolation.BicubicSplineInterpolator;

import main.CentroidField;
import main.DataManip;
import main.SRAW;
import drawing.LayerViewer;
import drawing.TopomapViewer;
import drawing.TopomapViewer_complex2;
import schrodinger.MovieMaker;
import util.matrix.Matrix;
import util.color.ColorScale;
import util.color.ColorScales;
import util.fileops.ColumnIO;
import util.fileops.FileOps;
import util.fileops.Layer;
import util.fileops.RHKFileOps;
import util.fileops.Topomap;
import util.fourier.FFT2DSmall;
import util.fourier.FFTOps;
import util.geom.AtomicCoordinatesSet;
import util.geom.Distance;
import util.geom.Mask;

public class DriftCorrectionMethods {
	/**
	 * Notes on drift correction:
	 * 
	 * In the automatic drift correction a topography, the fancy method should shift the u-field so that
	 * the bragg peaks are real. But, this should only be done for the topography. 
	 * 
	 * 
	 */
	//Directory: Important drift correction methods are in DriftCorrectionConsole (real space)
	//LayerUtil (bare k-space)
	//DataManip (bare real space)
	
	public static final int FANCY_TIME_192x192 = 19; //the duration of a fancy application of the u-field to a single 192x192 pixel layer. This is used to estimate other times.
	//This number has now been reduced from 102 to 18 seconds by my speeding up of things.
	public static long maxPrecomputedDoubles = 1024*1024*128;
	public static int[][] braggReserve = new int [2][2];
	
	public int[][] bragg = new int [2][2];
	public double[][] braggd = new double [2][2];//this is w.r.t. origin of the FFT
	public boolean willRegularize = false, regularizationDecided = false, forceEven = true, evenDecided = false;
	public double angle = 0;
	public Layer t, copy;
	BraggPeakTaker pt;
	int N;
	public JFileChooser fc;
	
	String sourceDir;
	
	int doingItAutomatic = 0; // 0 - not doing it automatic; 1 - doing it automatic; 2 - Making a map of the U-field
	// 3 - Doing strain stuff.
	
	FourierFilterMaskSource masker;
	/**
	 * This method attempts to obtain the Bragg peaks from the user using the DriftCorrectionConsole.
	 * @param t 
	 * @return
	 */
	
	public DriftCorrectionMethods()
	{
		fc = new JFileChooser(Topomap.stddir);
	}
	public void getBraggPeaksFromUser()
	{
		if (doingItAutomatic != 1){
			this.t = Layer.open(fc);
			sourceDir = fc.getCurrentDirectory().toString();
			copy = Layer.newLayer(t, FieldOps.copy(t.data));
			RHKFileOps.doFitting(copy, 10);
			N = t.nx;
		}
		int o = JOptionPane.showConfirmDialog(null, "Get the Bragg peaks from a file?");
		if (o == JOptionPane.YES_OPTION){
			bragg = FieldOps.transpose(ColumnIO.readIntTableASCII(FileOps.selectOpen(fc).toString()));
			angle = Math.toRadians(Double.parseDouble(JOptionPane.showInputDialog("Enter the expected angle.")));
			setBraggiAndProceed(bragg, angle);
		}
		else{
			if (doingItAutomatic == 1){
				this.t = Layer.open(fc);
				sourceDir = fc.getCurrentDirectory().toString();
				copy = Layer.newLayer(t, FieldOps.copy(t.data));
				RHKFileOps.doFitting(copy, 10);
				N = t.nx;
			}
			pt = new BraggPeakTaker(copy, this);
		}
	}

	public void setBraggiAndProceed(int[][] braggi, double angle) {
		this.bragg = braggi;
		if (regularizationDecided && !willRegularize)
			setBraggToEven();
		Printer.printlnHorizontal(bragg[0]);
		Printer.printlnHorizontal(bragg[1]);
		this.angle = angle;
		if (pt != null) {
			pt.dispose();
			pt = null;
		}
		if (doingItAutomatic == 0)
			getMaskerFromUser();
		else if (doingItAutomatic == 1)
		{
			doAutomaticDriftCorrection();
		}
		else if (doingItAutomatic == 2)
		{
			//Here we call method to generate everything:
			String input = JOptionPane.showInputDialog("Enter the length scale limits in pixels.");
			double lmin = Integer.parseInt(input.split(",")[0]);
			double lmax = Integer.parseInt(input.split(",")[1]);
			int nl = Integer.parseInt(JOptionPane.showInputDialog("How many length scales?"));
			
			DriftCorrectionMethods.getUFieldTopomap_Real(t, braggi, ArrayOps.generateArrayInclBoth(lmin, lmax, nl));
		}
		else if (doingItAutomatic == 3)
		{
			String f = FileOps.selectSave(fc).toString();
			ColumnIO.writeTable(FieldOps.transpose(bragg), f + "Bragg.txt", "");
			double[][] braggTrue = getBraggTrue(bragg, N);
			AtomicCoordinatesSet kspace = new AtomicCoordinatesSet(braggTrue[0], braggTrue[1], new double[] {N/2, N/2});
			AtomicCoordinatesSet real = kspace.getReciprocalLattice();
			ColumnIO.writeString(real.toString(), f + "lattice.txt");
			int choice = Integer.parseInt(JOptionPane.showInputDialog("What would you like to do?\r\n"+
					"1 - Do a single-length-scale strain calculation on a single layer\r\n"+
					"2 - Do a multi-length-scale strain calculation on a single layer\r\n"+
					"3 - Do a single-length scale strain calcaulation on each layer of a topomap\r\n"));
				
			if (choice == 1){
				String input = JOptionPane.showInputDialog("Enter the length scale in pixels.");
				double L = Integer.parseInt(input);
				double[][][] results = measureTheStrain(t, braggi, L);
				Layer uaal = Layer.newLayer(t, results[0]);
				Layer uabl = Layer.newLayer(t, results[1]);
				Layer ubal = Layer.newLayer(t, results[2]);
				Layer ubbl = Layer.newLayer(t, results[3]);
				Layer e12l = Layer.newLayer(t, results[4]);
				Layer w12l = Layer.newLayer(t, results[5]);
				Layer tre = Layer.newLayer(t, results[6]);
				Layer anis = Layer.newLayer(t, results[7]);
				
				Layer.writeBIN(uaal, f + "uaa.bin");
				Layer.writeBIN(uabl, f + "uab.bin");
				Layer.writeBIN(ubal, f + "uba.bin");
				Layer.writeBIN(ubbl, f + "ubb.bin");
				Layer.writeBIN(e12l, f + "e12.bin");
				Layer.writeBIN(w12l, f + "w12.bin");
				Layer.writeBIN(tre, f + "tr_e.bin");
				Layer.writeBIN(anis, f + "uanis.bin");
			}
			else if (choice == 2)
			{
				double[] Llimts = Printer.getTwoDoubles("Length Scale limits?");
				int nlayers = Printer.getAnInt("How many length scales?");
				double[] L = ArrayOps.generateArrayInclBoth(Llimts[0], Llimts[1], nlayers);
				double[][][][] results = new double [8][nlayers][t.nx][t.ny];
				for (int i = 0; i < nlayers; i++)
				{
					double[][][] temp = measureTheStrain(t, braggi, L[i]);
					for (int j = 0; j < 8; j++)
						results[j][i] = temp[j];
				}
				Topomap uaal = new Topomap(results[0], L, t.x, t.y, null);
				Topomap uabl = new Topomap(results[1], L, t.x, t.y, null);
				Topomap ubal = new Topomap(results[2], L, t.x, t.y, null);
				Topomap ubbl = new Topomap(results[3], L, t.x, t.y, null);
				Topomap e12l = new Topomap(results[4], L, t.x, t.y, null);
				Topomap w12l = new Topomap(results[5], L, t.x, t.y, null);
				Topomap tre = new Topomap(results[6], L, t.x, t.y, null);
				Topomap uanis = new Topomap(results[7], L, t.x, t.y, null);
				
				Topomap.writeBIN(uaal, f + "uaa.bin");
				Topomap.writeBIN(uabl, f + "uab.bin");
				Topomap.writeBIN(ubal, f + "uba.bin");
				Topomap.writeBIN(ubbl, f + "ubb.bin");
				Topomap.writeBIN(e12l, f + "e12.bin");
				Topomap.writeBIN(w12l, f + "w12.bin");
				Topomap.writeBIN(tre, f + "tr_e.bin");
				Topomap.writeBIN(uanis, f + "uanis.bin");
			}
			if (choice == 3){
				String input = JOptionPane.showInputDialog("Enter the length scale in pixels.");
				double L = Integer.parseInt(input);
				Topomap map = Topomap.open(fc);
				double[][][][] results = new double [8][map.nlayers][t.nx][t.ny]; //God willing this memory allocation
				//request will not crash the program.
				
				for (int i = 0; i < map.nlayers; i++){
					double[][][] temp = measureTheStrain(map.getLayer(i), braggi, L);
					for (int j = 0; j < 8; j++)
						results[j][i] = temp[j];
				}
				Topomap uaal = Topomap.newTopomap(map, results[0]);
				Topomap uabl = Topomap.newTopomap(map, results[1]);
				Topomap ubal = Topomap.newTopomap(map, results[2]);
				Topomap ubbl = Topomap.newTopomap(map, results[3]);
				Topomap e12l = Topomap.newTopomap(map, results[4]);
				Topomap w12l = Topomap.newTopomap(map, results[5]);
				Topomap tre = Topomap.newTopomap(map, results[6]);
				Topomap anis = Topomap.newTopomap(map, results[7]);
				
				Topomap.writeBIN(uaal, f + "uaa.bin");
				Topomap.writeBIN(uabl, f + "uab.bin");
				Topomap.writeBIN(ubal, f + "uba.bin");
				Topomap.writeBIN(ubbl, f + "ubb.bin");
				Topomap.writeBIN(e12l, f + "e12.bin");
				Topomap.writeBIN(w12l, f + "w12.bin");
				Topomap.writeBIN(tre, f + "tr_e.bin");
				Topomap.writeBIN(anis, f + "uanis.bin");
			}
			System.exit(0);
			
		}
			
	}
	private void setBraggToEven()
	{
		for (int i = 0; i < 2; i++){
			bragg[i][0] = FieldOps.roundEven(braggd[i][0]);
			bragg[i][1] = FieldOps.roundEven(braggd[i][1]);
		}

	}
	public void doAutomaticDriftCorrection()
	{
		int choice = Integer.parseInt(JOptionPane.showInputDialog(null, "What do you want to do?\r\n"
				+ "0 - Drift correct a single layer\r\n"
				+ "1 - Drift correct a dI/dV map topography and apply it to the maps (deprecated)\r\n"
				+ "2 - Drift correct a topomap\r\n"
				+ "3 - Drift correct a topomap (save memory by writing layers separately)\r\n"
				+ "4 - Drift correct a topomap (generate a separate output folder for each layer)\r\n"
				+ "         (Option 4 uses real-space calculation and non-fancy application\r\n"
				+ "         since no other options make sense)\r\n"
				));
		if (!regularizationDecided)
			willRegularize = JOptionPane.showConfirmDialog(null, "Regularize the Bragg peaks?") == JOptionPane.YES_OPTION;

		if (!evenDecided)
			forceEven = JOptionPane.showConfirmDialog(null, "Force the final Bragg peaks to be even?") == JOptionPane.YES_OPTION;
		switch (choice)
		{
		case 0:
			automaticDriftCorrectSingleLayer(true);
			break;
		case 1:
			automaticDriftCorrectDIDVMap(true, Topomap.open(fc), fc.getCurrentDirectory().toString() + "\\");
			break;
		case 2:
			automaticDriftCorrectTopomap(true);
			break;
		case 3:
			automaticDriftCorrectLayerSeries(true, Topomap.open(fc));
			break;
		case 4:
			automaticDriftCorrectLayerSeriesSplitFolder(true, Topomap.open(fc));
			break;
		}
	}
	public void automaticDriftCorrectSingleLayer(boolean askForMemory)
	{
		if (t == null){
			JOptionPane.showMessageDialog(null, "Open the layer.");
			t = Layer.open(fc);
		}
		N = t.nx;
		String dir = fc.getCurrentDirectory().toString() + "\\";
		String name = fc.getSelectedFile().getName().substring(0, fc.getSelectedFile().getName().length()-4);
		
		int startMB = UFieldCalculation.getMinimumMemUsage(t.nx)/(1024*1024);
		int extraMB = UFieldCalculation.getExtraMemUsagePerLayer(t.nx)/(1024*1024);
		int memoryMB = startMB < 200 ? startMB + 16*extraMB : startMB + extraMB;
		
		if (askForMemory)
			memoryMB = Integer.parseInt(JOptionPane.showInputDialog("How much memory do you want to use? (in MB)\r\n"
					+ "The minimum is " + startMB + " MB.\r\n"
					+ "Each length scale incurs an extra " + extraMB + " MB.\r\n"
					+ "In the box we allocated enough for " + (startMB < 200 ? 16 : 1) + "different length scales.", memoryMB));
		
		
		int kSpace = JOptionPane.showConfirmDialog(null, "To do it in k-space click Yes.\r\n To do it in real space, click No.");

		int secEstimate = (int)(Math.pow(N/192.0, 4) * 1 * FANCY_TIME_192x192);
		boolean fancy = JOptionPane.showConfirmDialog(null, "Attempt fancy application of the u-field?\r\nFancy way estimated to take " + (secEstimate/60) + " minutes per layer.") == JOptionPane.YES_OPTION;
		
		DriftCorrectionCalculator complete;
		int code = willRegularize ? 1 : 0;
		if (kSpace == JOptionPane.YES_OPTION)
			complete = UFieldCalculation.getFinished(UFieldCalculation.getDefaultMask(bragg, N, memoryMB*1024*1024), t, bragg, angle, fancy, forceEven, code);
		else
			complete = UFieldCalculationReal.getFinished(UFieldCalculationReal.getNLayers(N, memoryMB*1024*1024), t, bragg, angle, fancy, forceEven, code);
			
		complete.writeFileOutput(dir, name, willRegularize, null);
		LayerViewer.show(Layer.newLayer(t, complete.after), 1024, true);
		
	}
	public void automaticDriftCorrectLayerSeries(boolean askForMemory, Topomap map)
	{
		N = map.nx;
		String dir = fc.getCurrentDirectory().toString() + "\\";
		String name = fc.getSelectedFile().getName().substring(0, fc.getSelectedFile().getName().length()-4);
		
		int startMB = UFieldCalculation.getMinimumMemUsage(N)/(1024*1024);
		int extraMB = UFieldCalculation.getExtraMemUsagePerLayer(N)/(1024*1024);
		int memoryMB = startMB < 200 ? startMB + 16*extraMB : startMB + extraMB;
		
		if (askForMemory)
			memoryMB = Integer.parseInt(JOptionPane.showInputDialog("How much memory do you want to use? (in MB)\r\n"
					+ "The minimum is " + startMB + " MB.\r\n"
					+ "Each length scale incurs an extra " + extraMB + " MB.\r\n"
					+ "In the box we allocated enough for " + (startMB < 200 ? 16 : 1) + "different length scales.", memoryMB));
		
		
		int kSpace = JOptionPane.showConfirmDialog(null, "To do it in k-space click Yes.\r\n To do it in real space, click No.");

		int secEstimate = (int)(Math.pow(N/192.0, 4) * 1 * FANCY_TIME_192x192);
		boolean fancy = JOptionPane.showConfirmDialog(null, "Attempt fancy application of the u-field?\r\nFancy way estimated to take " + (secEstimate/60) + " minutes per layer.") == JOptionPane.YES_OPTION;
		
		//First let us write the Topomap layers separately.
		int nlayers = map.nlayers;
		File[] layers = new File[nlayers];
		String kOrR = (kSpace == JOptionPane.YES_OPTION ? "k" : "r");
		String outdir = dir + name + "_Drift corr auto " + kOrR + (willRegularize ? " reg" : "") + "\\";
		if (!new File(outdir).exists()) new File(outdir).mkdir();
		File[] phaseCont0 = new File[nlayers];
		File[] phaseCont1 = new File[nlayers];
		File[] afterf = new File[nlayers];
		for (int i = 0; i < nlayers; i++)
		{
			layers[i] = new File(outdir + "layer_" + i + ".bin");
			Layer.writeBIN(map.getLayer(i), layers[i].toString());
		
			afterf[i] = new File(dir + name + "_after" + (willRegularize ? "Reg" : "") + kOrR + " " +  i + ".bin");
			phaseCont0[i] = new File(outdir + "phaseCont_0_" + i + ".bin");
			phaseCont1[i] = new File(outdir + "phaseCont_1_" + i + ".bin");
		}
		
		map = null;
		
		DriftCorrectionCalculator complete;
		int code = willRegularize ? 1 : 0;
		for (int i = 0; i < nlayers; i++){
			t = Layer.readBIN(layers[i].toString());
			if (kSpace == JOptionPane.YES_OPTION){
				if (kSpace == JOptionPane.YES_OPTION)
					complete = UFieldCalculation.getFinished(UFieldCalculation.getDefaultMask(bragg, N, memoryMB*1024*1024), t, bragg, angle, fancy, forceEven, code);
				else
					complete = UFieldCalculationReal.getFinished(UFieldCalculationReal.getNLayers(N, memoryMB*1024*1024), t, bragg, angle, fancy, forceEven, code);
				complete.writeFileOutput(dir, name, willRegularize, null, i, kSpace == JOptionPane.NO_OPTION);
				
//				LayerViewer.show(Layer.newLayer(t, complete.after), 1024);
			}
		}
		
		//After it's all done, read the files and write topomaps.
		Layer[] temp = new Layer[nlayers];
		for (int i = 0; i < nlayers; i++)
			temp[i] = Layer.readBIN(phaseCont0[i].toString());
		Topomap pc0 = Topomap.newTopomap(temp);
		Topomap.writeBIN(pc0, outdir + "phaseCont_0.bin");
		for (int i = 0; i < nlayers; i++)
			temp[i] = Layer.readBIN(phaseCont1[i].toString());
		Topomap pc1 = Topomap.newTopomap(temp);
		Topomap.writeBIN(pc1, outdir + "phaseCont_1.bin");
		for (int i = 0; i < nlayers; i++)
			temp[i] = Layer.readBIN(afterf[i].toString());
		Topomap after = Topomap.newTopomap(temp);
		Topomap.writeBIN(after, dir + name + "_after" + (willRegularize ? "Reg" : "") + kOrR + " " + ".bin");
		for (int i = 0; i < nlayers; i++)
		{
			layers[i].delete();
			phaseCont0[i].delete();
			phaseCont1[i].delete();
			afterf[i].delete();
		}
		
	}
	public void automaticDriftCorrectLayerSeriesSplitFolder(boolean askForMemory, Topomap map)
	{
		N = map.nx;
		String dir = fc.getCurrentDirectory().toString() + "\\";
		String name = fc.getSelectedFile().getName().substring(0, fc.getSelectedFile().getName().length()-4);
		
		int startMB = UFieldCalculation.getMinimumMemUsage(N)/(1024*1024);
		int extraMB = UFieldCalculation.getExtraMemUsagePerLayer(N)/(1024*1024);
		int memoryMB = startMB < 200 ? startMB + 16*extraMB : startMB + extraMB;
		
		if (askForMemory)
			memoryMB = Integer.parseInt(JOptionPane.showInputDialog("How much memory do you want to use? (in MB)\r\n"
					+ "The minimum is " + startMB + " MB.\r\n"
					+ "Each length scale incurs an extra " + extraMB + " MB.\r\n"
					+ "In the box we allocated enough for " + (startMB < 200 ? 16 : 1) + "different length scales.", memoryMB));
		
		
//		int kSpace = JOptionPane.showConfirmDialog(null, "To do it in k-space click Yes.\r\n To do it in real space, click No.");

		//First let us write the Topomap layers separately.
		int nlayers = map.nlayers;
		File[] layers = new File[nlayers];
		String kOrR = "r";//(kSpace == JOptionPane.YES_OPTION ? "k" : "r");
		String outdir[] = new String[map.nlayers];
		for (int i = 0; i < map.nlayers; i++){
			outdir[i] = dir + name + "_tempLayer_" + i + "_Drift corr auto " + kOrR + (willRegularize ? " reg" : "") + "\\";
			if (!new File(outdir[i]).exists()) new File(outdir[i]).mkdir();

		}
		File[] phaseCont0 = new File[nlayers];
		File[] phaseCont1 = new File[nlayers];
		File[] afterf = new File[nlayers];
		for (int i = 0; i < nlayers; i++)
		{
			layers[i] = new File(dir + name + "_layer_" + i + ".bin");
			Layer.writeBIN(map.getLayer(i), layers[i].toString());
		
			afterf[i] = new File(dir + name + "_after" + (willRegularize ? "Reg" : "") + kOrR + " " +  i + ".bin");
			phaseCont0[i] = new File(outdir[i] + "phaseCont_0.bin");
			phaseCont1[i] = new File(outdir[i] + "phaseCont_1.bin");
		}
		
		map = null;
		
		DriftCorrectionCalculator complete;
		int code = willRegularize ? 1 : 0;
		for (int i = 0; i < nlayers; i++){
			t = Layer.readBIN(layers[i].toString());
//			if (kSpace == JOptionPane.YES_OPTION){
//				UFieldCalculation complete;
//				if (regularize)
//					complete = UFieldCalculation.getFinished_Both(UFieldCalculation.getDefaultMask(bragg, N, memoryMB*1024*1024), t, bragg, angle, false);
//				else
//					complete = UFieldCalculation.getFinished_NonlinearOnly(UFieldCalculation.getDefaultMask(bragg, N, memoryMB*1024*1024), t, bragg, false);
//				
//				complete.writeFileOutput(dir, name, regularize, null, i);
////				LayerViewer.show(Layer.newLayer(t, complete.after), 1024);
//			}
//			else if (kSpace == JOptionPane.NO_OPTION){
			complete = UFieldCalculationReal.getFinished(UFieldCalculationReal.getNLayers(N, memoryMB*1024*1024), t, bragg, angle, false, forceEven, code);
				
			complete.writeFileOutputToDir(outdir[i], name, willRegularize, null);
//				LayerViewer.show(Layer.newLayer(t, complete.after), 1024);
//			}
		}
	}
	
	/**
	 * This drift corrects one stack of layers (e.g. the topography) and saves the output, so that it can later be applied to the dI/dV.
	 * This method is more complicated than just a loop over auto drift correct, because we need to write the file output slightly differently.
	 *  
	 * @param askForMemory
	 */
	public void automaticDriftCorrectTopomap(boolean askForMemory)
	{
		JOptionPane.showMessageDialog(null, "Open the map.");
		Topomap map = Topomap.open(fc);
		N = map.nx;
		String dir = fc.getCurrentDirectory().toString() + "\\";
		String name = fc.getSelectedFile().getName().substring(0, fc.getSelectedFile().getName().length()-4);
		
		int startMB = UFieldCalculation.getMinimumMemUsage(N)/(1024*1024);
		int extraMB = UFieldCalculation.getExtraMemUsagePerLayer(N)/(1024*1024);
		int memoryMB = startMB < 200 ? startMB + 16*extraMB : startMB + extraMB;
		
		if (askForMemory)
			memoryMB = Integer.parseInt(JOptionPane.showInputDialog("How much memory do you want to use? (in MB)\r\n"
					+ "The minimum is " + startMB + " MB.\r\n"
					+ "Each length scale incurs an extra " + extraMB + " MB.\r\n"
					+ "In the box we allocated enough for " + (startMB < 200 ? 16 : 1) + "different length scales.", memoryMB));
		
		
		int kSpace = JOptionPane.showConfirmDialog(null, "To do it in k-space click Yes.\r\n To do it in real space, click No.");
		
		int secEstimate = (int)(Math.pow(N/192.0, 4) * 1 * FANCY_TIME_192x192);
		boolean fancy = JOptionPane.showConfirmDialog(null, "Attempt fancy application of the u-field?\r\nFancy way estimated to take " + (secEstimate/60) + " minutes per layer.") == JOptionPane.YES_OPTION;

		File outdir = new File(dir + name + "_Drift corr auto " + (kSpace == JOptionPane.YES_OPTION ? "k" : "r") + (willRegularize ? " reg" : "") + "\\");
		if (!outdir.exists()) outdir.mkdir();
		String outputDir = outdir.toString() + "\\";

		double[][][] phaseContX = new double [map.nlayers][][];
		double[][][] phaseContY = new double [map.nlayers][][];
		double[][][] after = new double [map.nlayers][][];
		
		DriftCorrectionCalculator complete = null;
		int code = willRegularize ? 1 : 0;
			
			//Here goes!
			for (int i = 0; i < map.nlayers; i++){
				if (kSpace == JOptionPane.YES_OPTION)
					complete = UFieldCalculation.getFinished(UFieldCalculation.getDefaultMask(bragg, N, memoryMB*1024*1024), map.getLayer(i), bragg, angle, fancy, forceEven, code);
				else
					complete = UFieldCalculationReal.getFinished(UFieldCalculationReal.getNLayers(N, memoryMB*1024*1024), map.getLayer(i), bragg, angle, fancy, forceEven, code);
					
				ColumnIO.writeLines(complete.getEnergyLines(), outputDir + "Selection info_ " + i + "_.txt");
				after[i] = complete.after;
				phaseContX[i] = complete.phaseCont[complete.selectedLayer][0];
				phaseContY[i] = complete.phaseCont[complete.selectedLayer][1];
				SRAW.writeImage(outputDir + "outsidePixels_" + i + "_.bmp", complete.outsidePixels);
			}
			ColumnIO.writeString(complete.latt.toString(), outputDir + "lattice.txt");
			
			ColumnIO.writeTable(FieldOps.transpose(complete.bragg), outputDir + "Bragg.txt", "");
			

			Topomap.writeBIN(Topomap.newTopomap(map, phaseContX), outputDir + "phaseCont_0.bin");
			Topomap.writeBIN(Topomap.newTopomap(map, phaseContY), outputDir + "phaseCont_1.bin");
			Topomap.writeBIN(Topomap.newTopomap(map, after), dir + name + "_after" + (willRegularize ? "Reg" : "") + (kSpace == JOptionPane.YES_OPTION ? "k" : "r") + ".bin");
			
			if (willRegularize){
				ColumnIO.writeString(complete.lattReg.toString(), outputDir + "latticeReg.txt");
				Layer.writeBIN(Layer.newLayer(map.getLayer(0), FieldOps.getIndex(complete.uReg, 0)), outputDir + "uRegX.bin");
				Layer.writeBIN(Layer.newLayer(map.getLayer(0), FieldOps.getIndex(complete.uReg, 1)), outputDir + "uRegY.bin");
				ColumnIO.writeTable(FieldOps.transpose(complete.braggReg), outputDir + "BraggReg.txt", "");
			}
			
//			complete.writeFileOutput(dir, name, regularize, null);
//			LayerViewer.show(Layer.newLayer(t, complete.after), 1024);
		
	}
	public void automaticDriftCorrectDIDVMap(boolean askForMemory, Topomap map, String alternateDir)
	{
		String dir;
		if (t == null){
			JOptionPane.showMessageDialog(null, "Open the layer.");
			t = Layer.open(fc);
			dir = fc.getCurrentDirectory().toString() + "\\";
		}
		else
			dir = alternateDir;
		N = t.nx;
		String name = fc.getSelectedFile().getName().substring(0, fc.getSelectedFile().getName().length()-4);
		
		
		
		int startMB = UFieldCalculation.getMinimumMemUsage(t.nx)/(1024*1024);
		int extraMB = UFieldCalculation.getExtraMemUsagePerLayer(t.nx)/(1024*1024);
		int memoryMB = startMB < 200 ? startMB + 16*extraMB : startMB + extraMB;
		
		if (askForMemory)
			memoryMB = Integer.parseInt(JOptionPane.showInputDialog("How much memory do you want to use? (in MB)\r\n"
					+ "The minimum is " + startMB + " MB.\r\n"
					+ "Each length scale incurs an extra " + extraMB + " MB.\r\n"
					+ "In the box we allocated enough for " + (startMB < 200 ? 16 : 1) + "different length scales.", memoryMB));
		

		UFieldCalculation complete;
		copy = Layer.newLayer(t, FieldOps.copy(t.data));
		FieldOps.subtractNSheetFit(copy.data, 3);
		int secEstimate = (int)(Math.pow(N/192.0, 4) * 1 * FANCY_TIME_192x192);
		boolean fancy = JOptionPane.showConfirmDialog(null, "Attempt fancy application of the u-field?\r\nFancy way estimated to take " + (secEstimate/60) + " minutes per layer.") == JOptionPane.YES_OPTION;
		int code = willRegularize ? 1 : 0;
		complete = UFieldCalculation.getFinished(UFieldCalculation.getDefaultMask(bragg, N, memoryMB*1024*1024), t, bragg, angle, fancy, forceEven, code);
		
		Topomap after = TopomapUtil.applyUField(complete.u, map);
		complete.writeFileOutput(dir, name, willRegularize, t);
		Topomap.writeBIN(after, dir + "map_after.bin");
//		LayerViewer.show(Layer.newLayer(t, complete.after), 1024);
	}
	
	public void getMaskerFromUser()
	{
		int o = Integer.parseInt(JOptionPane.showInputDialog("What type of Fourier filter do you want to use?\r\n" +
				"0) No Fourier-filter; do it in real space\r\n" + 
				"1) Sharp-edged circle\r\n" +
				"2) Circle with a Fermi-Dirac tapering edge\r\n" +
				"3) Gaussian distribution\r\n" +
				"4) Custom user-defined shape\r\n" +
				"5) Pointy shape (points along x- and y-axes)\r\n"
				));
		double rmin = 0, rmax = 0, T = 0;
		
		String[] radiusMessage = new String [2];
		if (o == 1 || o == 2 || o == 5)
			radiusMessage[0] = "Enter the minimum radius to be used in the calculation, in piexls.";
		if (o == 3)
			radiusMessage[0] = "Enter the minimum gaussian width in piexls.";
		radiusMessage[1] = "Now enter the maximum.";
		int n = 0;
		if (o == 1 || o == 2 || o == 3 || o == 5)
		{
			rmin = Double.parseDouble(JOptionPane.showInputDialog(radiusMessage[0]));
			rmax = Double.parseDouble(JOptionPane.showInputDialog(radiusMessage[1]));
			n = Integer.parseInt(JOptionPane.showInputDialog("How many length scales shall we calculate?"));
		}
		if (o == 2 || o == 5)
			T = Double.parseDouble(JOptionPane.showInputDialog("Enter the width of the edge in pixels."));
		
		switch(o)
		{
		case 0:
		{
			int nl = Integer.parseInt(JOptionPane.showInputDialog("How many length scales shall we calculate?"));
			double[] lengths = new double [nl];
			if (nl <= 1){
			String line = JOptionPane.showInputDialog("Enter them, separated by commas");
			String[] tok = line.split(",");
			for (int i = 0; i < nl; i++)
				lengths[i] = Double.parseDouble(tok[i]);
			}
			else
			{
				int oo = JOptionPane.showConfirmDialog(null, "Length scales evenly spaced?");
				if (oo == JOptionPane.YES_OPTION)
				{
					String line = JOptionPane.showInputDialog("Enter the minimum and maximum, inclusive.");
					String[] tok = line.split(",");
					lengths = ArrayOps.generateArrayInclBoth(Double.parseDouble(tok[0]), Double.parseDouble(tok[1]), lengths.length);
				}
				else{
					String line = JOptionPane.showInputDialog("Enter them, separated by commas");
					String[] tok = line.split(",");
					for (int i = 0; i < nl; i++)
						lengths[i] = Double.parseDouble(tok[i]);
				}
			}
			File base = FileOps.selectSave(fc);
			
			UFieldCalculationReal real = new UFieldCalculationReal(t, null, bragg, angle, lengths, true);
			real.doExpUCalc();
			real.pickBestLayer();
			real.makeBestUField();
			real.applyUField(true);
			
			
			double[][][] phaseContX = new double [lengths.length][][];
			double[][][] phaseContY = new double [lengths.length][][];
			for (int i = 0; i < lengths.length; i++)
			{
				phaseContX[i] = real.phaseCont[i][0];
				phaseContY[i] = real.phaseCont[i][1];
			}
			
			Topomap pcx = new Topomap(phaseContX, lengths, t.x, t.y, null);
			Topomap pcy = new Topomap(phaseContY, lengths, t.x, t.y, null);
			Topomap.writeBIN(pcx, base.toString() + "phaseCont0.bin");
			Topomap.writeBIN(pcy, base.toString() + "phaseCont1.bin");
			Layer.writeBIN(Layer.newLayer(t, real.after), base.toString()+"after.bin");
			break;
		}
		case 1:
			masker = new CircleMask(rmin, rmax, n, N);
			askForLineMaskerAndExecute();
			break;
		case 2:
			masker = new FermiMask(rmin, rmax, n, N, T);
			askForLineMaskerAndExecute();
			break;
		case 3:
			masker = new GaussMask(rmin, rmax, n, N);
			askForLineMaskerAndExecute();
			break;
		case 4:
			new PeakEncircler(copy, this, angle, bragg);
			break;
		case 5:
			masker = new FourPointedShapeMask(rmin, rmax, n, N, T);
			askForLineMaskerAndExecute();
			break;
		}
	}
	
	/**
	 * We assume that the name of the file is Bragg.txt.
	 * @param dir
	 */
	public void writeBraggFile(String dir)
	{
		ColumnIO.writeTable(FieldOps.transpose(bragg), dir + "Bragg.txt", "");
	}
	public void readBraggFile(String dir)
	{
		if (dir == null)
		{
			bragg = ColumnIO.readIntTableASCII(FileOps.selectOpen(fc).toString());
		}
		else
			bragg = ColumnIO.readIntTableASCII(dir + "Bragg.txt");
	}
	
	public void writeLatticeFile(String dir)
	{
		if (N == 0) N = Integer.parseInt(JOptionPane.showInputDialog("The topography size was unknown. Enter the size in pixels of the topography."));
		String path = dir == null ? FileOps.selectSave(fc).toString() : dir + "lattice.txt";
		double[][] braggTrue = getBraggTrue(bragg, N);
		AtomicCoordinatesSet kspace = new AtomicCoordinatesSet(braggTrue[0], braggTrue[1], new double[] {N/2, N/2});
		AtomicCoordinatesSet real = kspace.getReciprocalLattice();
		ColumnIO.writeString(real.toString(), path);
	}
	public AtomicCoordinatesSet getLattice()
	{
		if (N == 0) N = Integer.parseInt(JOptionPane.showInputDialog("The topography size was unknown. Enter the size in pixels of the topography."));
		double[][] braggTrue = getBraggTrue(bragg, N);
		AtomicCoordinatesSet kspace = new AtomicCoordinatesSet(braggTrue[0], braggTrue[1], new double[] {N/2, N/2});
		AtomicCoordinatesSet real = kspace.getReciprocalLattice();
		return real;
	}
	
	
	public static double[][] getBraggTrue(int[][] bragg, int N)
	{
		double[][] ans = new double [2][2];
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				ans[i][j] = bragg[i][j]*2*Math.PI/N;
		return ans;
	}
	
	public void askForLineMaskerAndExecute()
	{
		int o = JOptionPane.showConfirmDialog(null, "Do you want to include lines along the x- and y-directions of k-space from the Bragg peaks?");
		if (o == JOptionPane.YES_OPTION)
		{
			int length = Integer.parseInt(JOptionPane.showInputDialog("How long should they be? (The layer is " + t.nx + " pixels."));
			executeWholeCalculation(true, new MaskWithExtendedStrips(masker, bragg, length, t));
		}
		else
			executeWholeCalculation(true, masker);
	}
	
	public void applyAutoUFieldToTopomap(Topomap t)
	{
		double[][][] u = //DriftCorrectionAnalysis.getUFromAutoDriftCorrFolder(null, fc);
				getULayerFromFolder(null, fc);
		double nxOver192 = ((double)t.nx)/(192.0);
		int secEstimate = (int)(Math.pow(nxOver192, 4) * t.nlayers * FANCY_TIME_192x192);
		
		boolean doItFancy = JOptionPane.showConfirmDialog(null, "Try to do it the fancy way?\r\nFancy way estimated to take " + (secEstimate/60) + " minutes.") == JOptionPane.YES_OPTION;
		if (!doItFancy)
		{
			Topomap tu = TopomapUtil.applyUField(u, t);
			Topomap.writeBIN(tu, fc);
		}
		else
		{
			//Now let the user decide whether to split the calculation into layers
			int choice = JOptionPane.showConfirmDialog(null, "This may take a long time. Do you wish to save the layers one by one\r\n"
					+ "so that if the calculation is interrupted, you don't have to restart totally?");
			if (choice == JOptionPane.YES_OPTION)
			{
				JOptionPane.showMessageDialog(null, "Enter the name of the final map file WITHOUT SUFFIX.");
				File base = FileOps.selectSave(fc);
				File[] temp = new File[t.nlayers];
				boolean filesExist = false;
				double[][][] result = new double [t.nlayers][t.nx][t.ny];
				//This part does the calcualtion while skipping, if the user desires it, those files which already exist.
				String reuseThese = "The following files already exist:\r\n";
				boolean[] fileExists = new boolean [t.nlayers];
				int nExistingFiles = 0;
				for (int i = 0; i < t.nlayers; i++)
				{
					temp[i] = new File(base + "_tempLayer_" + i + ".bin");
					if (temp[i].exists()){
						fileExists[i] = true;
						filesExist = true;
						reuseThese += temp[i].getName() + "\r\n";
						nExistingFiles++;
					}
				}
				reuseThese += "\r\nShall we skip them?";
				if (nExistingFiles > 30)
				{
					reuseThese = "The following files already exist:\r\nat least 30 files\r\n\r\nShall we skip them?";
				}
				int skip = filesExist ? JOptionPane.showConfirmDialog(null, reuseThese) : JOptionPane.NO_OPTION;
				if (skip == JOptionPane.NO_OPTION)
				{	//This option applies the u-field manually to all the files
					for (int i = 0; i < t.nlayers; i++)
					{
						System.out.println("Doing " + (i+1) + " out of " + t.nlayers); 
						double[][] tempdata = applyUFieldSpecial(u, t.data[i]);
						boolean[][] outsidePixels = FieldOps.getOutsidePixels(u);
						FieldOps.zero(tempdata, outsidePixels);
						FieldOps.changeZeroToAverage(tempdata);
						Layer.writeBIN(Layer.getFreeLayer(tempdata), temp[i].toString());
					}
				}
				else if (skip == JOptionPane.YES_OPTION)
				{
					for (int i = 0; i < t.nlayers; i++)
						if (!fileExists[i])
						{
							System.out.println("Doing " + (i+1) + " out of " + t.nlayers); 
							double[][] tempdata = applyUFieldSpecial(u, t.data[i]);
							boolean[][] outsidePixels = FieldOps.getOutsidePixels(u);
							FieldOps.zero(tempdata, outsidePixels);
							FieldOps.changeZeroToAverage(tempdata);
							Layer.writeBIN(Layer.getFreeLayer(tempdata), temp[i].toString());
						}
				}
				
				//The actual topomap will be assembled by reading from the saved files
				for (int i = 0; i < t.nlayers; i++){
					result[i] = Layer.openFree(temp[i]).data;
				}
				
				Topomap tu = Topomap.newTopomap(t, result);
				Topomap.writeBIN(tu, base.toString());
				
				for (int i = 0; i < t.nlayers; i++)
					temp[i].delete();
			
			}
			if (choice == JOptionPane.NO_OPTION)
			{	
				JOptionPane.showMessageDialog(null, "Enter the name of the final map file.");
				File f = FileOps.selectSave(fc);
				double[][][] result = new double [t.nlayers][t.nx][t.ny];
				for (int i = 0; i < t.nlayers; i++)
				{
					System.out.println("Doing " + (i+1) + " out of " + t.nlayers); 
					double[][] temp = applyUFieldSpecial(u, t.data[i]);
					boolean[][] outsidePixels = FieldOps.getOutsidePixels(u);
					FieldOps.zero(temp, outsidePixels);
					FieldOps.changeZeroToAverage(temp);
					FieldOps.copy(temp, result[i]);
				}
				Topomap tu = Topomap.newTopomap(t, result);
				Topomap.writeBIN(tu, f.toString());
			}
			
		}

	}
	public void applyAutoUFieldTopomapToTopomap(Topomap t)
	{
//		double[][][] u = DriftCorrectionAnalysis.getUFromAutoDriftCorrFolder(null, fc);
		double[][][][] u = getUTopomapFromFolder(null, fc);
		
		double nxOver192 = ((double)t.nx)/(192.0);
		int secEstimate = (int)(Math.pow(nxOver192, 4) * t.nlayers * FANCY_TIME_192x192);
		
		boolean doItFancy = JOptionPane.showConfirmDialog(null, "Try to do it the fancy way?\r\nFancy way estimated to take " + (secEstimate/60) + " minutes.") == JOptionPane.YES_OPTION;
		if (!doItFancy)
		{
			File f = FileOps.selectSave(fc);
			Topomap tu = TopomapUtil.applyUField(u, t);
			Topomap.writeBIN(tu, f.toString());
		}
		else
		{	//Now let the user decide whether to split the calculation into layers
			int choice = JOptionPane.showConfirmDialog(null, "This may take a long time. Do you wish to save the layers one by one\r\n"
					+ "so that if the calculation is interrupted, you don't have to restart totally?");
			if (choice == JOptionPane.YES_OPTION)
			{
				JOptionPane.showMessageDialog(null, "Enter the name of the final map file WITHOUT SUFFIX.");
				File base = FileOps.selectSave(fc);
				File[] temp = new File[t.nlayers];
				boolean filesExist = false;
				double[][][] result = new double [t.nlayers][t.nx][t.ny];
				//This part does the calcualtion while skipping, if the user desires it, those files which already exist.
				String reuseThese = "The following files already exist:\r\n";
				boolean[] fileExists = new boolean [t.nlayers];
				int nExistingFiles = 0;
				for (int i = 0; i < t.nlayers; i++)
				{
					temp[i] = new File(base + "_tempLayer_" + i + ".bin");
					if (temp[i].exists()){
						fileExists[i] = true;
						filesExist = true;
						reuseThese += temp[i].getName() + "\r\n";
						nExistingFiles++;
					}
				}
				reuseThese += "\r\nShall we skip them?";
				if (nExistingFiles > 30)
				{
					reuseThese = "The following files already exist:\r\nat least 30 files\r\n\r\nShall we skip them?";
				}
				int skip = filesExist ? JOptionPane.showConfirmDialog(null, reuseThese) : JOptionPane.NO_OPTION;
				if (skip == JOptionPane.NO_OPTION)
				{	//This option applies the u-field manually to all the files
					for (int i = 0; i < t.nlayers; i++)
					{
						System.out.println("Doing " + (i+1) + " out of " + t.nlayers); 
						double[][] tempdata = applyUFieldSpecial(u[i], t.data[i]);
						boolean[][] outsidePixels = FieldOps.getOutsidePixels(u[i]);
						FieldOps.zero(tempdata, outsidePixels);
						FieldOps.changeZeroToAverage(tempdata);
						Layer.writeBIN(Layer.getFreeLayer(tempdata), temp[i].toString());
					}
				}
				else if (skip == JOptionPane.YES_OPTION)
				{
					for (int i = 0; i < t.nlayers; i++)
						if (!fileExists[i])
						{
							System.out.println("Doing " + (i+1) + " out of " + t.nlayers); 
							double[][] tempdata = applyUFieldSpecial(u[i], t.data[i]);
							boolean[][] outsidePixels = FieldOps.getOutsidePixels(u[i]);
							FieldOps.zero(tempdata, outsidePixels);
							FieldOps.changeZeroToAverage(tempdata);
							Layer.writeBIN(Layer.getFreeLayer(tempdata), temp[i].toString());
						}
				}
				
				//The actual topomap will be assembled by reading from the saved files
				for (int i = 0; i < t.nlayers; i++){
					result[i] = Layer.openFree(temp[i]).data;
				}
				
				Topomap tu = Topomap.newTopomap(t, result);
				Topomap.writeBIN(tu, base.toString());
				
				for (int i = 0; i < t.nlayers; i++)
					temp[i].delete();
			
			}
			if (choice == JOptionPane.NO_OPTION)
			{	
				JOptionPane.showMessageDialog(null, "Enter the name of the final map file.");
				File f = FileOps.selectSave(fc);
				double[][][] result = new double [t.nlayers][t.nx][t.ny];
				for (int i = 0; i < t.nlayers; i++)
				{
					System.out.println("Doing " + (i+1) + " out of " + t.nlayers); 
					double[][] temp = applyUFieldSpecial(u[i], t.data[i]);
					boolean[][] outsidePixels = FieldOps.getOutsidePixels(u[i]);
					FieldOps.zero(temp, outsidePixels);
					FieldOps.changeZeroToAverage(temp);
					FieldOps.copy(temp, result[i]);
				}
				Topomap tu = Topomap.newTopomap(t, result);
				Topomap.writeBIN(tu, f.toString());
			}
		}

	}
	
	/**
	 * This assumes that the directory in question contains nothing but topomaps to be drift
	 * corrected by the same u[][][][], non-fancily.
	 * 
	 */
	public void applyAutoUFieldTopomapToDirectory(String dir)
	{
//		double[][][] u = DriftCorrectionAnalysis.getUFromAutoDriftCorrFolder(null, fc);
		JOptionPane.showMessageDialog(null, "Open the drift correction folder.");
		double[][][][] u = getUTopomapFromFolder(null, fc);
		JOptionPane.showMessageDialog(null, "Now select the target directory.");
		String outdir = FileOps.selectDir(fc);
//		double nxOver192 = ((double)t.nx)/(192.0);
//		int secEstimate = (int)(Math.pow(nxOver192, 4) * t.nlayers * FANCY_TIME_192x192);
//		
//		boolean doItFancy = JOptionPane.showConfirmDialog(null, "Try to do it the fancy way?\r\nFancy way estimated to take " + (secEstimate/60) + " minutes.") == JOptionPane.YES_OPTION;
//		if (!doItFancy)
//		{
//			File f = FileOps.selectSave(fc);
//			Topomap tu = TopomapUtil.applyUField(u, t);
//			Topomap.writeBIN(tu, f.toString());
//		}
		File[] f = new File(dir).listFiles();
		for (int i = 0; i < f.length; i++)
			if (f[i].toString().endsWith(".bin")){
			{
				Topomap t = Topomap.readBIN(f[i].toString());
				Topomap tu = TopomapUtil.applyUField(u, t);
				Topomap.writeBIN(tu, outdir + f[i].getName());
			}
			}

	}
	public static double[][][][] getUTopomapFromFolder(String dir, JFileChooser fc)
	{
		if (dir == null) dir = FileOps.selectDir(fc);
		Topomap phaseCont0 = Topomap.readBIN(dir + "phaseCont_0.bin");
		Topomap phaseCont1 = Topomap.readBIN(dir + "phaseCont_1.bin");
		AtomicCoordinatesSet latt = new AtomicCoordinatesSet(ColumnIO.readString(dir + "lattice.txt"));
		
		double[][] braggTrue = latt.getReciprocal();
		
		double[][][][] u = new double [phaseCont0.nlayers][phaseCont0.nx][phaseCont1.ny][2];
		for (int i = 0; i < phaseCont0.nlayers; i++)
			FieldOps.putU(phaseCont0.data[i], phaseCont1.data[i], braggTrue, 1, u[i]);
		
		double[][] extraTranslation = new double [phaseCont0.nlayers][2];
		File[] translation = new File[phaseCont0.nlayers];
		for (int i = 0; i < translation.length; i++)
		{	translation[i] = new File(dir + "translation_" + i + ".txt");
			if (translation[i].exists()){
				String tr = ColumnIO.readString(translation[i].toString());
				String[] trs = tr.split("\t");
				extraTranslation[i][0] = Double.parseDouble(trs[0]);
				extraTranslation[i][1] = Double.parseDouble(trs[1]);
			}
		}
		for (int k = 0; k < phaseCont0.nlayers; k++)
			for (int i = 0; i < phaseCont0.nx; i++)
				for (int j = 0; j < phaseCont0.ny; j++)
				{
					u[k][i][j][0] += extraTranslation[k][0];
					u[k][i][j][1] += extraTranslation[k][1];
				}

		
		java.io.File uRegX = new java.io.File(dir + "uRegX.bin");
		java.io.File uRegY = new java.io.File(dir + "uRegY.bin");
		
		if (uRegX.exists() && uRegY.exists()) //If it is a "regularized" u-field
		{
			braggReserve = FieldOps.transpose(ColumnIO.readIntTableASCII(dir + "BraggReg.txt"));
			double[][] urx = Layer.readBIN(uRegX.toString()).data;
			double[][] ury = Layer.readBIN(uRegY.toString()).data;
			for (int k = 0; k < phaseCont0.nlayers; k++)
				for (int i = 0; i < phaseCont0.nx; i++)
					for (int j = 0; j < phaseCont0.ny; j++)
					{
						u[k][i][j][0] += urx[i][j];
						u[k][i][j][1] += ury[i][j];
					}
		}
		else
			braggReserve = FieldOps.transpose(ColumnIO.readIntTableASCII(dir + "Bragg.txt"));

		return u;

	}
	public static double[][][] getULayerFromFolder(String dir, JFileChooser fc)
	{
		if (dir == null) dir = FileOps.selectDir(fc);
		Layer phaseCont0 = Layer.readBIN(dir + "phaseCont_0.bin");
		Layer phaseCont1 = Layer.readBIN(dir + "phaseCont_1.bin");
		AtomicCoordinatesSet latt = new AtomicCoordinatesSet(ColumnIO.readString(dir + "lattice.txt"));
		
		double[][][] u = new double [phaseCont0.nx][phaseCont1.ny][2];

		double[] extraTranslation = new double [2];
		File translation = new File(dir + "translation.txt");
		if (translation.exists()){
			String tr = ColumnIO.readString(translation.toString());
			String[] trs = tr.split("\t");
			extraTranslation[0] = Double.parseDouble(trs[0]);
			extraTranslation[1] = Double.parseDouble(trs[1]);
		}
		System.out.println(Printer.arrayLnHorizontal(extraTranslation));
		
		double[][] braggTrue = latt.getReciprocal();
		
		FieldOps.putU(phaseCont0.data, phaseCont1.data, braggTrue, 1, u);
		for (int i = 0; i < phaseCont0.nx; i++)
			for (int j = 0; j < phaseCont0.ny; j++)
			{
				u[i][j][0] += extraTranslation[0];
				u[i][j][1] += extraTranslation[1];
			}
		
		java.io.File uRegX = new java.io.File(dir + "uRegX.bin");
		java.io.File uRegY = new java.io.File(dir + "uRegY.bin");

		if (uRegX.exists() && uRegY.exists()) //If it is a "regularized" u-field
		{
			braggReserve = FieldOps.transpose(ColumnIO.readIntTableASCII(dir + "BraggReg.txt"));
			double[][] urx = Layer.readBIN(uRegX.toString()).data;
			double[][] ury = Layer.readBIN(uRegY.toString()).data;
			for (int i = 0; i < phaseCont0.nx; i++)
				for (int j = 0; j < phaseCont0.ny; j++)
				{
					u[i][j][0] += urx[i][j];
					u[i][j][1] += ury[i][j];
				}
		}
		else
			braggReserve = FieldOps.transpose(ColumnIO.readIntTableASCII(dir + "Bragg.txt"));
		
		return u;

	}
	/**
	 * This does the entire calculation at once, as in LayerUtil. The method is copied directly from
	 * LayerUtil except the business of getting the mask.
	 */
	public void executeWholeCalculation(boolean writePictures, FourierFilterMaskSource mask) {
		
		boolean regularize = JOptionPane.showConfirmDialog(null, "Regularize the Bragg peaks?") == JOptionPane.YES_OPTION;
		// TODO Auto-generated method stub
		String outDir = sourceDir + "\\Drift Corr\\"; 
		File f = new File(outDir);
		if (!f.exists()) f.mkdirs();
		String imDir = outDir + "Image output\\";
		File fim = new File(imDir);
		if (!fim.exists() && writePictures) fim.mkdirs();
		String lattDir = outDir + "Lattices\\";
		File flatt = new File(lattDir);
		if (!flatt.exists()) flatt.mkdirs();
		File select = new File(outDir + "SelectedLayer.txt");
		try {
			select.createNewFile();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		int nlayers = mask.getNInLoop();
		
		writeBraggFile(outDir);
		
		String braggOut = "The bragg peaks were " + Printer.vectorP(bragg[0]) + " and " + Printer.vectorP(bragg[1]) + ".\r\n";
		braggOut += "In the image, that is " + Printer.vectorP(new int[] {bragg[0][0] + N/2, bragg[0][1] + N/2}) +" and " +Printer.vectorP(new int[] {bragg[1][0] + N/2, bragg[1][1] + N/2}) + ".\r\n";
		braggOut += mask.toString();
		ColumnIO.writeString(braggOut, outDir + "Bragg Info.txt");
		
		AtomicCoordinatesSet latt = this.getLattice();
		ColumnIO.writeString(latt.toString(), lattDir + "lattice.txt");
		ColumnIO.writeString(latt.toString(), outDir + "lattice.txt");
		ColumnIO.writeString(latt.getRt2Lattice().toString(), lattDir + "latticeRt2.txt");
		ColumnIO.writeString(latt.getRt2Lattice().getRt2Lattice().toString(), lattDir + "lattice2x2.txt");
		
		
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
			FieldOps.shift(fftz, shiftedFFTZ[i], bragg[i][0], bragg[i][1]);
		double temp;
		
		double naturalMin = FieldOps.magMin(fftz);
		//First loop: get expU for each bragg and mask, and write the output before application of the ufield
		Topomap phase;
		double[] indices = ArrayOps.generateArrayInclBoth(1, mask.getNInLoop(), mask.getNInLoop());
		for (int n = 0; n < 2; n++){
			for (int i = 0; i <mask.getNInLoop(); i++)
			{
				LayerUtil.putExpUFourier(shiftedFFTZ[n], mask.getMask(n, i), expU, expUfft);
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
			braggTrue[i][0] = bragg[i][0]*2*Math.PI/t.nx;
			braggTrue[i][1] = bragg[i][1]*2*Math.PI/t.nx;
		}
		double[][][][] u = new double [mask.getNInLoop()][t.nx][t.ny][2];
		BufferedImage uIM;
		double[][][] after = new double[mask.getNInLoop()][t.nx][t.ny];
		BicubicSplineInterpolatingFunction interp = null;
		BicubicSplineInterpolator erp = new BicubicSplineInterpolator();
		interp = erp.interpolate(ArrayOps.generateArrayInclBoth(0, t.nx-1, t.nx), ArrayOps.generateArrayInclBoth(0, t.ny-1, t.ny), t.data);
		double mean = FieldOps.mean(t.data);
		for (int k = 0; k < mask.getNInLoop(); k++)
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
			FieldOps.putU(phaseCont[0], phaseCont[1], braggTrue, 1, u[k]);
			
			//The next five lines are supposed to put an atom on top of the center pixel. 
			double[][] phaseContCent0 = {{((phaseCont[0][N/2][N/2] + Math.PI) % (2*Math.PI) - Math.PI)}};
			double[][] phaseContCent1 = {{((phaseCont[1][N/2][N/2] + Math.PI) % (2*Math.PI) - Math.PI)}};
			double[][][] uCent = {{{0, 0}}};
			FieldOps.putU(phaseContCent0, phaseContCent1, braggTrue, 1, uCent); //the -0.5 is so that the "center of the pixel" is on top of the atomic coordinate
			FieldOps.addVectorTo(u[k], new double[]{uCent[0][0][0] - u[k][N/2][N/2][0] - 0.5, uCent[0][0][1] - u[k][N/2][N/2][1] - 0.5});
			
			if (writePictures){
				uIM = ImageEditing.getBufferedImage(u[k], null);
				SRAW.writeImage(imDir + "u_" + MovieMaker.fromInt(k), uIM);
			}
			if (!regularize) FieldOps.applyUFieldBiCubic(interp, u[k], after[k], mean);
		}
		
		double[][][] ux = new double [mask.getNInLoop()][][];
		double[][][] uy = new double [mask.getNInLoop()][][];
		for (int i = 0; i < mask.getNInLoop(); i++)
		{
			ux[i] = FieldOps.copy(u[i], 0);
			uy[i] = FieldOps.copy(u[i], 1);
		}
		
		
		if (regularize)
		{
			BraggVectorRegularizer reg = new BraggVectorRegularizer(angle, N, JOptionPane.showConfirmDialog(null, "Make the Bragg peaks even?") == JOptionPane.YES_OPTION); 
			reg.setFirstBraggPeak(braggTrue[0]);
			reg.setSecondBraggPeak(braggTrue[1]);
			reg.regularize();
			AtomicCoordinatesSet coordReg = reg.getRegLattice();
			ColumnIO.writeString(coordReg.toString(), lattDir + "coordReg.txt");
			ColumnIO.writeString(coordReg.getRt2Lattice().toString(), lattDir + "coordRegRt2.txt");
			ColumnIO.writeString(coordReg.getRt2Lattice().getRt2Lattice().toString(), lattDir + "coordReg2x2.txt");
			double[][][] uReg = new double [N][N][2];
			FieldOps.putUField(getLattice(), coordReg, uReg);
			double[][] uRegx = FieldOps.getIndex(uReg, 0);
			double[][] uRegy = FieldOps.getIndex(uReg, 1);
			
			Layer.writeBIN(Layer.newLayer(t, uRegx), outDir + "uRegX.bin");
			Layer.writeBIN(Layer.newLayer(t, uRegy), outDir + "uRegY.bin");
			int[][] braggReg = reg.getFinalBragg();
			ColumnIO.writeTable(FieldOps.transpose(braggReg), outDir + "BraggReg.txt", "");
			for (int i = 0; i < nlayers; i++)
			{
				FieldOps.add(ux[i], uRegx, ux[i]);
				FieldOps.add(uy[i], uRegy, uy[i]);
				FieldOps.putIntoArray(ux[i], uy[i], u[i]);
				FieldOps.applyUFieldBiCubic(interp, u[i], after[i], mean);
			}
		}
		
		if (nlayers > 1)
		{
			Topomap uxt = new Topomap(ux, indices, t.x, t.y, null);
			Topomap uyt = new Topomap(uy, indices, t.x, t.y, null);
			Topomap.writeBIN(uxt, outDir + "ux.bin");
			Topomap.writeBIN(uyt, outDir + "uy.bin");
			
			Topomap done = new Topomap(after, indices, t.x, t.y, null);
			Topomap.writeBIN(done, outDir + "after.bin");
			new TopomapViewer(done, outDir, N*2);
		}
		else
		{
			Layer uxt = new Layer(ux[0], t.x, t.y, t.v, t.current);
			Layer uyt = new Layer(uy[0], t.x, t.y, t.v, t.current);
			Layer.writeBIN(uxt, outDir + "ux.bin");
			Layer.writeBIN(uyt, outDir + "uy.bin");
			
			Layer done = new Layer(after[0], t.x, t.y, t.v, t.current);
			Layer.writeBIN(done, outDir + "after.bin");
		}
//		System.exit(0);

	}

	
	public void setMasker(PeakEncircler enc)
	{
		String[] radiusMessage = new String [2];
		radiusMessage[0] = "The edge of the shape will smoothly decay to zero like as a Gaussian.\r\n" +
				"We will use a range of Gaussian widths. What should be the minimum? (in pixels)";
		radiusMessage[1] = "Now enter the maximum.";
		int n = 0;
		
		double rmin, rmax;
		rmin = Double.parseDouble(JOptionPane.showInputDialog(radiusMessage[0]));
		rmax = Double.parseDouble(JOptionPane.showInputDialog(radiusMessage[1]));
		n = Integer.parseInt(JOptionPane.showInputDialog("How many calculations shall we make (including the minimum and maximum)?"));
		
		masker = new CustomShapeMask(rmin, rmax, n, enc.isInPeak, bragg);
		RHKFileOps.doFitting(t, 12);
		askForLineMaskerAndExecute();
	}
	/**
	 * This class takes the Bragg peaks from the user.
	 * @author madhavanlab2011
	 *
	 */
	public static class BraggPeakTaker extends JFrame implements KeyListener, MouseListener, MouseMotionListener
	{
		int size, sizeRatio;
		public int ox = 20, oy = 40;
		BufferedImage image;
		public int N;
		Box[] boxes = new Box[1];
		double[][] ftmag;
		int xdown, xup, ydown, yup, currentx, currenty;
		int currentBragg = 0;
		double[][] bragg = new double [2][2];
		int[][] braggi = new int [2][2];
		DriftCorrectionMethods parent;
		static int WIDTH = 1200, HEIGHT = 1024;
		double angle;
		public BraggPeakTaker(Layer t, DriftCorrectionMethods parent)
		{
			BufferedImage fft;
			this.parent = parent;
			int size = Math.max(t.nx, 1024);
			N = t.nx;
			size -= size % t.nx;
			sizeRatio = size/t.nx;
			ftmag = new double[t.nx][t.nx];
			FFTOps.obtainFFTmagCent(t.data, ftmag);
			FieldOps.log(ftmag);
			
			image = new BufferedImage(size, size, BufferedImage.TYPE_INT_RGB);
			ImageEditing.enlargeBasic(ImageEditing.getBufferedImage(ftmag, 1), image, sizeRatio);
			addKeyListener(this);
			addMouseListener(this);
			addMouseMotionListener(this);
			setSize(WIDTH, HEIGHT);
			repaint();
			show();
			setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			int angDeg = Integer.parseInt(JOptionPane.showInputDialog(null, "Enter the expected angle between the peaks, in degrees.", 90));
			angle = Math.toRadians(angDeg);
			JOptionPane.showMessageDialog(null, "Please draw a box around the first Bragg peak. It shouldn't be too big.");
		}
		
		public void paint(Graphics g)
		{
			g.clearRect(0, 0, 2000, 2000);
			g.drawImage(image, ox, oy, null);
			g.setColor(java.awt.Color.BLUE);
			if (boxes[0] != null) boxes[0].draw(g, this);

		}

		public void mouseDragged(MouseEvent arg0) {
			// TODO Auto-generated method stub
			int x = arg0.getX();
			int y = arg0.getY();
			currentx = x; currenty = y;
		}
		@Override
		public void mouseMoved(MouseEvent arg0) {
			// TODO Auto-generated method stub
			currentx = arg0.getX();
			currenty = arg0.getY();
		}
		@Override
		public void mouseClicked(MouseEvent arg0) {
		}
		@Override
		public void mouseEntered(MouseEvent arg0) {
		}
		@Override
		public void mouseExited(MouseEvent arg0) {
		}
		@Override
		public void mousePressed(MouseEvent arg0) {
			// TODO Auto-generated method stub
			xdown = arg0.getX();
			ydown = arg0.getY();
		}
		@Override
		public void mouseReleased(MouseEvent arg0) {
			// TODO Auto-generated method stub
			xup = arg0.getX();
			yup = arg0.getY();
			if(boxes[0] == null)
			{
				int[] dr = {(xup-xdown)/sizeRatio, (yup-ydown)/sizeRatio};
				createNewBox(this.ftrCent(xdown, ydown), dr);
				repaint();
			}
		}
		@Override
		public void keyPressed(KeyEvent arg0) {
			// TODO Auto-generated method stub
			
		}
		@Override
		public void keyReleased(KeyEvent arg0) {
			// TODO Auto-generated method stub
			
		}
		@Override
		public void keyTyped(KeyEvent arg0) {
			// TODO Auto-generated method stub
			if (arg0.getKeyChar() == 'a')// || arg0.getKeyChar() == 'A')
				moveBox(-1, 0);
			if (arg0.getKeyChar() == 'd')// || arg0.getKeyChar() == 'D')
				moveBox(+1, 0);
			if (arg0.getKeyChar() == 'w')// || arg0.getKeyChar() == 'W')
				moveBox(0, -1);
			if (arg0.getKeyChar() == 's')// || arg0.getKeyChar() == 's')
				moveBox(0, +1);
			if (arg0.getKeyChar() == 'A')// || arg0.getKeyChar() == 'A')
				resizeBox(-1, 0);
			if (arg0.getKeyChar() == 'D')// || arg0.getKeyChar() == 'D')
				resizeBox(+1, 0);
			if (arg0.getKeyChar() == 'W')// || arg0.getKeyChar() == 'W')
				resizeBox(0, -1);
			if (arg0.getKeyChar() == 'S')// || arg0.getKeyChar() == 's')
				resizeBox(0, +1);
			if (arg0.getKeyChar() == ' ')
			{
				bragg[currentBragg] = CentroidField.centroid(ftmag, boxes[0].x, boxes[0].dx, boxes[0].y, boxes[0].dy, true);
				double[] nowBragg = {bragg[currentBragg][0]-N/2, bragg[currentBragg][1]-N/2};
				currentBragg++;
				double[] nextBragg = new double [2];
				Matrix.putProductWith(Matrix.getRotationMatrix(angle), nowBragg, nextBragg);
				nextBragg[0] += N/2; nextBragg[1] += N/2;
				boxes[0].move((int)nextBragg[0], (int)nextBragg[1], this);
				if (currentBragg == 2){
					boolean forceEven = JOptionPane.showConfirmDialog(null, "To force even-numbered Bragg Peaks, click YES.") == JOptionPane.YES_OPTION;
					parent.forceEven = forceEven;
					parent.evenDecided = true;
					if (forceEven){
						parent.regularizationDecided = true;
						parent.willRegularize = JOptionPane.showConfirmDialog(null, "Are you planning to regularize the Bragg peaks?") == JOptionPane.YES_OPTION;
					}
					for (int i = 0; i < 2; i++){
						bragg[i][0] -= N/2;
						bragg[i][1] -= N/2;
						braggi[i][0] = FieldOps.round(bragg[i][0]);
						braggi[i][1] = FieldOps.round(bragg[i][1]);
					}
				parent.braggd = bragg;
				parent.setBraggiAndProceed(braggi, angle);
				}
				else{
					
					JOptionPane.showMessageDialog(null, "Now the other one.");
					repaint();
				}
			}
		}
		int[] ftrCent(int x, int y)
		{
			return new int[] {(x-ox)/sizeRatio - N/2, (y-oy)/sizeRatio - N/2};
		}
		public void createNewBox(int[] r, int [] dr)
		{
			boxes[0] = new Box(r[0]+N/2, r[1]+N/2, dr[0], dr[1]);
		}

		public void moveBox(int dx, int dy)
		{
			boxes[0].x += dx;
			boxes[0].y += dy;
			repaint();
		}
		public void resizeBox(int dx, int dy)
		{
			boxes[0].dx += dx;
			boxes[0].dy += dy;
			repaint();
		}
		 private static class Box
		 {
			 int x, y, dx, dy;

			public Box(int x, int y, int dx, int dy) {
				this.x = x;
				this.y = y;
				this.dx = dx;
				this.dy = dy;
			}
			
			public void draw(Graphics g, BraggPeakTaker parent)
			{
				g.drawRect(parent.ox + x*parent.sizeRatio, parent.oy + y*parent.sizeRatio, dx*parent.sizeRatio, dy*parent.sizeRatio);
			}
			public void move(int x, int y, BraggPeakTaker parent)
			{
				this.x = x > 0 ? x+dx < parent.N ? x : parent.N - dx : 0; this.y = y > 0 ? y+dy < parent.N ? y : parent.N - dy : 0;
			}
			public void setSize(int sizex, int sizey)
			{
				this.dy = sizex;
				this.dy = sizey;
			}
			public int[] getR()
			{
				return new int[] {x, y};
			}
			public int[] getDR()
			{
				return new int[] {dx, dy};
			}
		 }
	
	}
	
	public static class PeakEncircler extends JFrame implements KeyListener, MouseListener, MouseMotionListener
	{
		int size, sizeRatio;
		public int ox = 20, oy = 40;
		BufferedImage image;
		public int N;
		Square[] boxes = new Square[1];
		double[][] ftmag, peakArea;
		int xdown, xup, ydown, yup, currentx, currenty;
		int currentBragg = 0;
		double[][] bragg = new double [2][2];
		int[][] braggi = new int [2][2];
		DriftCorrectionMethods parent;
		static int WIDTH = 1200, HEIGHT = 1024;
		double angle;
		private boolean zoomedIn = false;
		ColorScale scale;
		int scaleIndex = 9;
		int currenti, currentj;
		boolean[][] isBorder;
		
		public boolean[][][] isInPeak;
		public PeakEncircler(Layer t, DriftCorrectionMethods parent, double angle, int[][] braggi)
		{
			this.angle = angle;
			this.parent = parent;
			this.braggi = braggi;
			int size = 1024;
			N = t.nx;
			size -= size % t.nx;
			sizeRatio = size/t.nx;
			ftmag = new double[t.nx][t.nx];
			FFTOps.obtainFFTmagCent(t.data, ftmag);
			FieldOps.log(ftmag);
			
			isInPeak = new boolean [2][N][N];
			scale = ColorScales.getNew(ftmag, scaleIndex);
			image = new BufferedImage(size, size, BufferedImage.TYPE_INT_RGB);
			ImageEditing.enlargeBasic(ImageEditing.getBufferedImage(ftmag, scale), image, sizeRatio);
			addKeyListener(this);
			addMouseListener(this);
			addMouseMotionListener(this);
			setSize(WIDTH, HEIGHT);
			repaint();
			show();
			setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			JOptionPane.showMessageDialog(null, "Draw a box enclosing the ENTIRE Bragg peak with a good margin. Then press Spacebar.");
		}
		
		public void paint(Graphics g)
		{
			g.clearRect(0, 0, 2000, 2000);
			g.drawImage(image, ox, oy, null);
			g.setColor(ColorScales.getUnusedColor(scaleIndex));
			if (boxes[0] != null && !zoomedIn) boxes[0].draw(g, this);

		}

		private void resetZoomedInImage() {
			// TODO Auto-generated method stub
			BufferedImage temp = ImageEditing.getBufferedImage(peakArea, isBorder, scale, ColorScales.getUnusedColor(scaleIndex));
			ImageEditing.enlargeBasic(temp, image, sizeRatio);
		}
		private void resetZoomedOutImage() {
			// TODO Auto-generated method stub
			int size = 1024;
			size -= size % N;
			sizeRatio = size/N;
			BufferedImage temp = ImageEditing.getBufferedImage(ftmag, scale);
			image = new BufferedImage(size, size, BufferedImage.TYPE_INT_RGB);
			ImageEditing.enlargeBasic(temp, image, sizeRatio);
			
		}

		public void mouseDragged(MouseEvent arg0) {
			// TODO Auto-generated method stub
			int x = arg0.getX();
			int y = arg0.getY();
			currentx = x; currenty = y;
			if (zoomedIn)
			{
				setCurrentIJ(currentx, currenty);
				if (!isBorder[currenti][currentj]){
					isBorder[currenti][currentj] = true;
					resetZoomedInImage();
					repaint();
				}
			}
		}
		@Override
		public void mouseMoved(MouseEvent arg0) {
			// TODO Auto-generated method stub
			currentx = arg0.getX();
			currenty = arg0.getY();
			setCurrentIJ(currentx, currenty);
			setTitle("" + currenti + ", " + currentj);
		}
		@Override
		public void mouseClicked(MouseEvent arg0) {
		}
		@Override
		public void mouseEntered(MouseEvent arg0) {
		}
		@Override
		public void mouseExited(MouseEvent arg0) {
		}
		@Override
		public void mousePressed(MouseEvent arg0) {
			// TODO Auto-generated method stub
			xdown = arg0.getX();
			ydown = arg0.getY();
		}
		@Override
		public void mouseReleased(MouseEvent arg0) {
			// TODO Auto-generated method stub
			xup = arg0.getX();
			yup = arg0.getY();
			if(/*boxes[0] == null*/ !zoomedIn)
			{
				int[] dr = {(xup-xdown)/sizeRatio, (yup-ydown)/sizeRatio};
				createNewBox(this.ftrCent(xdown, ydown), Math.max(dr[0], dr[1]));
				repaint();
			}
		}
		@Override
		public void keyPressed(KeyEvent arg0) {
			// TODO Auto-generated method stub
			
		}
		@Override
		public void keyReleased(KeyEvent arg0) {
			// TODO Auto-generated method stub
			
		}
		@Override
		public void keyTyped(KeyEvent arg0) {
			// TODO Auto-generated method stub
			if (arg0.getKeyChar() == 'a')// || arg0.getKeyChar() == 'A')
				moveBox(-1, 0);
			if (arg0.getKeyChar() == 'd')// || arg0.getKeyChar() == 'D')
				moveBox(+1, 0);
			if (arg0.getKeyChar() == 'w')// || arg0.getKeyChar() == 'W')
				moveBox(0, -1);
			if (arg0.getKeyChar() == 's')// || arg0.getKeyChar() == 's')
				moveBox(0, +1);
			if (arg0.getKeyChar() == 'A')// || arg0.getKeyChar() == 'A')
				resizeBox(-1);
			if (arg0.getKeyChar() == 'D')// || arg0.getKeyChar() == 'D')
				resizeBox(+1);
			if (arg0.getKeyChar() == 'W')// || arg0.getKeyChar() == 'W')
				resizeBox(-1);
			if (arg0.getKeyChar() == 'S')// || arg0.getKeyChar() == 's')
				resizeBox(+1);
			if (arg0.getKeyChar() == ' ')
			{
				if (!zoomedIn)
				{
					zoomImageIn();
					if (currentBragg == 0)
						JOptionPane.showMessageDialog(null, "Draw a shape completely enclosing the Bragg peak.\r\n" +
							"The Bragg peak will not include the boundary, but you can add an additional smoothing distance later.");
					
				}
				else{
					boolean[][] braggPeak = FieldOps.getContiguousBlob(braggi[currentBragg][0]-boxes[0].x+N/2, braggi[currentBragg][1]-boxes[0].y+N/2, isBorder);
					isBorder = braggPeak;
					FieldOps.negate(isBorder);
					resetZoomedInImage();
					repaint();
					
					int response = JOptionPane.showConfirmDialog(null, "Is this Bragg peak OK?", "Confirmation", JOptionPane.YES_NO_OPTION);
					if (response == JOptionPane.NO_OPTION)
					{
						isBorder = new boolean[boxes[0].s][boxes[0].s];
						resetZoomedInImage();
						repaint();
					}
					else if (response == JOptionPane.YES_OPTION)
					{
						FieldOps.negate(braggPeak);
						for (int i = 0; i < boxes[0].s; i++)
							for (int j = 0; j < boxes[0].s; j++)
								isInPeak[currentBragg][i+boxes[0].x][j+boxes[0].y] = braggPeak[i][j]; 
						if (currentBragg == 0)
						{
							resetZoomedOutImage();
							zoomedIn = false;
							boxes[0].move(braggi[1][0] - boxes[0].s/2 + N/2, braggi[1][1] -boxes[0].s/2 + N/2, this);

							repaint();
							currentBragg++;
							JOptionPane.showMessageDialog(null, "Repeat for the other Bragg peak.");
						}
						else
						{
							dispose();
							parent.setMasker(this);
						}
					}
				}
			}
		}
		private void zoomImageIn() {
			// TODO Auto-generated method stub
			peakArea = boxes[0].getContents(ftmag);
			int size = 1024;
			size -= size % boxes[0].s;
			isBorder = new boolean[boxes[0].s][boxes[0].s];
			sizeRatio = size/boxes[0].s;
			image = new BufferedImage(size, size, BufferedImage.TYPE_INT_RGB);
			BufferedImage temp = ImageEditing.getBufferedImage(peakArea, scale);
			ImageEditing.enlargeBasic(temp, image, sizeRatio);
			repaint();
			zoomedIn = true;
		}

		int[] ftrCent(int x, int y)
		{
			return new int[] {(x-ox)/sizeRatio - N/2, (y-oy)/sizeRatio - N/2};
		}
		void setCurrentIJ(int screenX, int screenY)
		{
			currenti = (screenX-ox)/sizeRatio;
			currentj = (screenY-oy)/sizeRatio;
		}
		public void createNewBox(int[] r, int s)
		{
			boxes[0] = new Square(r[0]+N/2, r[1]+N/2, s);
		}

		public void moveBox(int dx, int dy)
		{
			boxes[0].x += dx;
			boxes[0].y += dy;
			repaint();
		}
		public void resizeBox(int change)
		{
			boxes[0].s += change;
			boxes[0].s += change;
			repaint();
		}
		 private static class Square
		 {
			 int x, y, s;

			public Square(int x, int y, int s) {
				this.x = x;
				this.y = y;
				this.s = s;
			}
			
			public void draw(Graphics g, PeakEncircler parent)
			{
				g.drawRect(parent.ox + x*parent.sizeRatio, parent.oy + y*parent.sizeRatio, s*parent.sizeRatio, s*parent.sizeRatio);
			}
			public void move(int x, int y, PeakEncircler parent)
			{
				this.x = x > 0 ? x+s < parent.N ? x : parent.N - s : 0; this.y = y > 0 ? y+s < parent.N ? y : parent.N - s : 0;
			}
			public void setSize(int size)
			{
				this.s = size;
			}
			public int[] getR()
			{
				return new int[] {x, y};
			}
			public int[] getDR()
			{
				return new int[] {s, s};
			}
			public double[][] getContents(double[][] data)
			{
				double[][] results = new double [s][s];
				for (int i = 0; i < s; i++)
					for (int j = 0; j < s; j++)
						results[i][j] = data[i+x][j+y];
				return results;
			}
		 }
	}
	
	public static interface FourierFilterMaskSource
	{
		public double[][] getMask(int braggIndex, int loopIndex);
		public double getLengthScale(int loopIndex);
		public int getNInLoop();
		public String toString();
	}
	public static class CircleMask implements FourierFilterMaskSource
	{
		double[] rValues;
		int n;
		int N;
		
		public CircleMask(double rmin, double rmax, int n, int N)
		{
			double sqmin = rmin*rmin;
			double sqmax = rmax*rmax;
			rValues = ArrayOps.generateArrayInclBoth(sqmin, sqmax, n);
			for (double r : rValues)
				r = Math.sqrt(r);
			this.n = n;
			this.N = N;
		}
		public double[][] getMask(int braggIndex, int loopIndex) {
			double[][] mask = new double [N][N];
			TopomapUtil.FourierFilterMethods.putCircleMaskOrigin(rValues[loopIndex], mask);
			return mask;
		}
		public int getNInLoop() {
			return n;
		}
		public String toString()
		{
			String ans = "Cirlce Mask " + n + " radii:\r\n";
			ans += "Image # \t Radius\r\n";
			for (int i = 0; i < n; i++)
				ans += "" + i + "\t" + rValues[i] + "\r\n";
			return ans;
		}
		@Override
		public double getLengthScale(int loopIndex) {
			// TODO Auto-generated method stub
			return rValues[loopIndex];
		}
	}
	public static class FermiMask implements FourierFilterMaskSource
	{
		double[] rValues;
		int n;
		int N;
		double T;
		
		public FermiMask(double rmin, double rmax, int n, int N, double T)
		{
			this.T = T;
			double sqmin = rmin*rmin;
			double sqmax = rmax*rmax;
			this.n = n;
			rValues = ArrayOps.generateArrayInclBoth(sqmin, sqmax, n);
			for (int i = 0; i < rValues.length; i++)
				rValues[i] = Math.sqrt(rValues[i]);
			this.N = N;
		}
		public double[][] getMask(int braggIndex, int loopIndex) {
			double[][] mask = new double [N][N];
			TopomapUtil.FourierFilterMethods.putFermiMaskOrigin(rValues[loopIndex], T, mask);
			return mask;
		}
		public int getNInLoop() {
			return n;
		}
		public String toString()
		{
			String ans = "Fermi Mask " + n + " radii:\r\n";
			ans += "Width of edge = " + T + "\r\n"; 
			ans += "Image # \t Radius\r\n";
			for (int i = 0; i < n; i++)
				ans += "" + i + "\t" + rValues[i] + "\r\n";
			return ans;
		}
		public double getLengthScale(int loopIndex) {
			// TODO Auto-generated method stub
			return rValues[loopIndex];
		}
	}
	public static class GaussMask implements FourierFilterMaskSource
	{
		double[] rValues;
		int n;
		int N;
		
		public GaussMask(double rmin, double rmax, int n, int N)
		{
			double sqmin = rmin*rmin;
			double sqmax = rmax*rmax;
			this.n = n;
			rValues = ArrayOps.generateArrayInclBoth(sqmin, sqmax, n);
			for (double r : rValues)
				r = Math.sqrt(r);
			this.N = N;
		}
		public double[][] getMask(int braggIndex, int loopIndex) {
			double[][] mask = new double [N][N];
			TopomapUtil.FourierFilterMethods.putGaussMaskOrigin(rValues[loopIndex], mask);
			return mask;
		}
		public int getNInLoop() {
			return n;
		}
		public String toString()
		{
			String ans = "Gaussian Mask " + n + " decay lengths:\r\n";
			ans += "Image # \t Decay length\r\n";
			for (int i = 0; i < n; i++)
				ans += "" + i + "\t" + rValues[i] + "\r\n";
			return ans;
		}
		public double getLengthScale(int loopIndex) {
			// TODO Auto-generated method stub
			return rValues[loopIndex];
		}
	}
	public static class CustomShapeMask implements FourierFilterMaskSource
	{
		boolean[][][] isBraggPeak; //The original unshifted designations. (As in the image).
		double[][][] distanceFromTrue; //The distances, unshifted as in the image.
		double[] width;
		int[][] braggi;
		int n;
		int N;
		
		public CustomShapeMask(double wmin, double wmax, int n, boolean[][][] isBraggPeak, int[][] braggi)
		{
			this.n = n;
			this.N = isBraggPeak[0].length;
			this.braggi = braggi;
			width = ArrayOps.generateArrayInclBoth(wmin, wmax, n);
			this.isBraggPeak = isBraggPeak;
			
			distanceFromTrue = new double [2][N][N];
			for (int i = 0; i < 2; i++)
				Distance.putDistanceFromTrueSparse(isBraggPeak[i], distanceFromTrue[i]);
//			new LayerViewer(Layer.getFreeLayer(distanceFromTrue[0]), "", 512);
			
		}
		@Override
		public double[][] getMask(int braggIndex, int loopIndex) {
			double[][] unshiftedmask = new double [N][N];
			double[][] shiftedmask = new double [N][N];
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					if (distanceFromTrue[braggIndex][i][j] > 6*width[loopIndex])
						unshiftedmask[i][j] = 0;
					else
						unshiftedmask[i][j] = gauss(distanceFromTrue[braggIndex][i][j], width[loopIndex]);
			FieldOps.shift(unshiftedmask, shiftedmask, (braggi[braggIndex][0]+N/2), (braggi[braggIndex][1]+N/2));
			return shiftedmask;
		}
		private static double gauss(double r, double width)
		{
			return Math.exp(-(r*r)/(2*width*width));
		}
		@Override
		public int getNInLoop() {
			return n;
		}
		
		public String toString()
		{
			String ans = "Custom Shape mask\r\n";
			ArrayList<int[]> truePoints = new ArrayList<int[]>();
			int[] temp = new int [2];
			truePoints = FieldOps.getTruePoints(isBraggPeak[0]);
			for (int i = 0; i < 2; i++)
			{
				truePoints = FieldOps.getTruePoints(isBraggPeak[i]);
				ans += "Bragg peak #" +i + " = " + Printer.vectorP(braggi[0]) + " used " + truePoints.size() + " pixels:\r\n";
				ans += "Pixel Coordinates\t W.r.t. Peak\r\n";
				for (int j = 0; j < truePoints.size(); j++)
				{
					temp[0] = truePoints.get(j)[0] - (braggi[i][0] + N/2);
					temp[1] = truePoints.get(j)[1] - (braggi[i][1] + N/2);
					ans += Printer.vectorP(truePoints.get(j))  + "\t" + Printer.vectorP(temp) + "\r\n";
				}
			}
			
			ans += "The following decay lenghts were used:\r\n";
			ans += "Image # \t Decay length\r\n";
			for (int i = 0; i < n; i++)
				ans += "" + i + "\t" + width[i] + "\r\n";
			return ans;
		}
		public double getLengthScale(int loopIndex) {
			// TODO Auto-generated method stub
			return width[loopIndex];
		}
	}
	/**
	 * This one is supposed to be a four-pointed shape, the outline of which is like a circle
	 * inscribed in a square, translated by one radius with periodic boundary conditions on the square.
	 * 
	 * @author madhavanlab2011
	 *
	 */
	public static class FourPointedShapeMask implements FourierFilterMaskSource
	{
		boolean[][][] mask; //The original unshifted designations. (As in the image).
		double[][][] distanceFromTrue; //The distances, unshifted as in the image.
		double[] r; 
		int n;
		int N;
		double T;
		
		public FourPointedShapeMask(double rmin, double rmax, int n, int N, double T)
		{
			this.N = N;
			this.T = T;
			this.n = n;
			if (n ==1 ) r = new double[] {(rmin+rmax)/2};
			else r = ArrayOps.generateArrayInclBoth(rmin, rmax, n);
			mask = new boolean [n][N][N];
			distanceFromTrue = new double [n][N][N];
			
			int x, y;
			//set up the mask:
			boolean[][] miniMask, miniMaskClone;
			for (int i = 0; i < n; i++)
			{
//				miniMask = new boolean[2*(int)(r[i]+1)][2*(int)(r[i]+1)];
//				miniMaskClone = FieldOps.copy(miniMask);
//				for (int p = 0; p < miniMask.length; p++)
//					for (int q = 0; q < miniMask.length; q++)
//						if (Distance.distance(p-miniMask.length/2, q-miniMask.length/2) < r[i])
//							miniMask[p][q] = false;
//						else miniMask[p][q] = true;
//				FieldOps.shift(miniMask, miniMaskClone, miniMask.length/2, miniMask.length/2);
//				for (int p = 0; p < miniMask.length; p++){System.out.println();
//					for (int q = 0; q < miniMask.length; q++)
//					{
//						x = ((p-miniMask.length/2)+N/2)%N;
//						y = ((q-miniMask.length/2)+N/2)%N;
//						if (miniMaskClone[p][q])
//							mask[i][x][y] = true;
//						if (miniMaskClone[p][q]) System.out.print("T ");
//						else System.out.print("F ");
//					}
//				}
//				Distance.putDistanceFromTrueSparseCenter(mask[i], distanceFromTrue[i], r[i]+10*T);
				Distance.putDistanceFromPointyShapeCenter(distanceFromTrue[i], r[i]);
				FieldOps.shift(FieldOps.copy(distanceFromTrue[i]), distanceFromTrue[i]);
			}
//			new TopomapViewer(Topomap.newTopomap(null, distanceFromTrue), "", 512);
//			new LayerViewer(Layer.getFreeLayer(distanceFromTrue[0]), "", 512);
		}
		@Override
		public double[][] getMask(int braggIndex, int loopIndex) {
			double[][] unshiftedmask = new double [N][N];
//			double[][] shiftedmask = new double [N][N];
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					unshiftedmask[i][j] = gauss(distanceFromTrue[loopIndex][i][j], T);
//			FieldOps.shift(unshiftedmask, shiftedmask, (braggi[braggIndex][0]+N/2), (braggi[braggIndex][1]+N/2));
			return unshiftedmask;
		}
		private static double gauss(double r, double width)
		{
			return Math.exp(-(r*r)/(2*width*width));
		}
		
		public int getNInLoop() {
			return n;
		}
		
		public String toString()
		{
			String ans = "Pointy shape Mask " + n + " radii:\r\n";
			ans += "Width of edge = " + T + "\r\n"; 
			ans += "Image # \t Radius\r\n";
			for (int i = 0; i < n; i++)
				ans += "" + i + "\t" + r[i] + "\r\n";
			return ans;
		}
		public double getLengthScale(int loopIndex) {
			// TODO Auto-generated method stub
			return r[loopIndex];
		}

	}
	
	/**
	 * This class takes any mask and adds to it a pair of long vertical and horizontal lines,
	 * in attempt to deal with what happens at the edge of the map.
	 * 
	 * The lines attempt to discriminate by not counting pixels as large as the Bragg peak.
	 * 
	 * @author madhavanlab2011
	 *
	 */
	public static class MaskWithExtendedLines implements FourierFilterMaskSource
	{
		int length;
		int N;
		int n;
		int[][] braggi;
		double[][][] lineMask;
		FourierFilterMaskSource masker;
		
		public MaskWithExtendedLines(FourierFilterMaskSource masker, int[][] braggi, int length, Layer t)
		{
			this.braggi = braggi;
			this.N = t.nx;
			this.masker = masker;
			n = masker.getNInLoop();
			lineMask = new double [2][N][N];
			int xb, yb;
			int xp, xm, yp, ym;
			double[][] ftmag = FFTOps.obtainFFTmagCent(t.data);
			double braggMag;
			for (int i = 0; i < 2; i++)
			{
				xb = braggi[i][0] + t.nx/2;
				yb = braggi[i][1] + t.ny/2;
				braggMag = ftmag[xb][yb];
				
				for (int k = 1; k < length; k++) //go out in all 4 directions, up to length.
				{
					xp = (xb+k+t.nx)%t.nx;
					yp = (yb+k+t.ny)%t.ny;
					xm = (xb-k+t.nx)%t.nx;
					ym = (yb-k+t.ny)%t.ny;
					if (ftmag[xp][yb] <= braggMag/k) lineMask[i][k][0] = 1.0/k;
					if (ftmag[xm][yb] <= braggMag/k) lineMask[i][lineMask[0].length-k][0] = 1.0/k;
					if (ftmag[xb][yp] <= braggMag/k) lineMask[i][0][k] = 1.0/k;
					if (ftmag[xb][ym] <= braggMag/k) lineMask[i][0][lineMask[0][0].length-k] = 1.0/k;
				}
			}
			
		}

		@Override
		public double[][] getMask(int braggIndex, int loopIndex) {
			// TODO Auto-generated method stub
			double[][] theMask = masker.getMask(braggIndex, loopIndex);
			for (int i = 0; i < theMask.length; i++)
				for (int j = 0; j < theMask.length; j++)
					theMask[i][j] = Math.max(theMask[i][j], lineMask[braggIndex][i][j]);
			return theMask;
		}
		@Override
		public int getNInLoop() {
			// TODO Auto-generated method stub
			return n;
		}
		public double getLengthScale(int loopIndex) {
			// TODO Auto-generated method stub
			return masker.getLengthScale(loopIndex);
		}

		
	}
	
	/**
	 * We try to emulate the fourier transforms produced by the real space method, by
	 * including wide strips along the x- and y-directions.
	 * The width of these strips is given by some fraction of the length scale of the masker.
	 * @author madhavanlab2011
	 *
	 */
	public static class MaskWithExtendedStrips implements FourierFilterMaskSource
	{
		int length;
		int N;
		int n;
		int[][] braggi;
		double[][][][] stripMask;
		FourierFilterMaskSource masker;
		
		public MaskWithExtendedStrips(FourierFilterMaskSource masker, int[][] braggi, int length, Layer t)
		{
			this.braggi = braggi;
			this.N = t.nx;
			this.masker = masker;
			n = masker.getNInLoop();
			stripMask = new double [2][n][N][N];
			int xb, yb;
			int xp, xm, yp, ym; 
			int xbw, ybw; //relative to the width
			int wx, wy;
			double width;
			int wmin, wmax;
			
			double[][] ftmag = FFTOps.obtainFFTmagCent(t.data);
			double braggMag;
			for (int i = 0; i < 2; i++)
			{
				xb = braggi[i][0] + t.nx/2;
				yb = braggi[i][1] + t.ny/2;
				braggMag = ftmag[xb][yb];
				for (int q = 0; q < n; q++)
				{
					width = masker.getLengthScale(q)*2;
					wmin = FieldOps.round(-width/2);
					wmax = FieldOps.round(+width/2);
					for (int w = wmin+1; w < wmax; w++)
						for (int k = 1; k < length; k++) //go out in all 4 directions, up to length.
						{
							xp = (xb+k+t.nx)%t.nx;
							yp = (yb+k+t.ny)%t.ny;
							xm = (xb-k+t.nx)%t.nx;
							ym = (yb-k+t.ny)%t.ny;
							xbw = (xb+w+t.nx)%t.nx;
							ybw = (yb+w+t.ny)%t.ny;
							wx = (w+t.nx)%t.nx;
							wy = (w+t.ny)%t.ny;
							if (ftmag[xp][ybw] <= braggMag/k) stripMask[i][q][k][wy] = 1.0;
							if (ftmag[xm][ybw] <= braggMag/k) stripMask[i][q][stripMask[0][0].length-k][wy] = 1.0;
							if (ftmag[xbw][yp] <= braggMag/k) stripMask[i][q][wx][k] = 1.0;
							if (ftmag[xbw][ym] <= braggMag/k) stripMask[i][q][wx][stripMask[0][0][0].length-k] = 1.0;
						}
				}
			}
			
		}

		@Override
		public double[][] getMask(int braggIndex, int loopIndex) {
			// TODO Auto-generated method stub
			double[][] theMask = masker.getMask(braggIndex, loopIndex);
			for (int i = 0; i < theMask.length; i++)
				for (int j = 0; j < theMask.length; j++)
					theMask[i][j] = Math.max(theMask[i][j], stripMask[braggIndex][loopIndex][i][j]);
			return theMask;
		}
		@Override
		public int getNInLoop() {
			// TODO Auto-generated method stub
			return n;
		}
		public double getLengthScale(int loopIndex) {
			// TODO Auto-generated method stub
			return masker.getLengthScale(loopIndex);
		}

		
	}
	
	public double[][][] convertToLatticeCoordinates(double[][][] u, AtomicCoordinatesSet latt)
	{
		double[][][] uLatt = new double[u.length][u[0].length][2];
		double[] temp = latt.getOrigin().clone();
		latt.setOrigin(0, 0);
		for (int i = 0; i < u.length; i++)
			for (int j = 0; j < u[0].length; j++)
			{
				latt.putAtomicCoords(u[i][j], uLatt[i][j]);
			}
		latt.setOrigin(temp[0], temp[1]);
		return uLatt;
	}
	public void putConvertedToLatticeCoordinates(double[][][] u, AtomicCoordinatesSet latt, double[][][] uLatt)
	{
		double[] temp = latt.getOrigin().clone();
		latt.setOrigin(0, 0);
		for (int i = 0; i < u.length; i++)
			for (int j = 0; j < u[0].length; j++)
			{
				latt.putAtomicCoordsVector(u[i][j], uLatt[i][j]);
			}
		latt.setOrigin(temp[0], temp[1]);
	}
	public AtomicCoordinatesSet loadFromDir(String dir)
	{
		return new AtomicCoordinatesSet(ColumnIO.readString(dir + "lattice.txt"));
	}
	public double[][][][] loadUFromTopomaps(Topomap[] us)
	{
		Topomap ux = us[0], uy = us[1];
		double[][][][] u = new double [ux.nlayers][ux.nx][ux.ny][2];
		for (int k = 0; k < ux.nlayers; k++)
			for (int i = 0; i < ux.nx; i++)
				for (int j = 0; j < ux.ny; j++)
				{
					u[k][i][j][0] = ux.data[k][i][j];
					u[k][i][j][1] = uy.data[k][i][j];
				}
		return u;
	}
	public double[][][] loadUFromLayers(Layer[] us)
	{
		Layer ux = us[0], uy = us[1];
		double[][][] u = new double [ux.nx][ux.ny][2];
		for (int i = 0; i < ux.nx; i++)
			for (int j = 0; j < ux.ny; j++)
			{
				u[i][j][0] = ux.data[i][j];
				u[i][j][1] = uy.data[i][j];
			}
		return u;
	}
	public double[][][] loadUFromTopomaps(Topomap[] us, int k)
	{
		Topomap ux = us[0], uy = us[1];
		double[][][] u = new double [ux.nx][ux.ny][2];
			for (int i = 0; i < ux.nx; i++)
				for (int j = 0; j < ux.ny; j++)
				{
					u[i][j][0] = ux.data[k][i][j];
					u[i][j][1] = uy.data[k][i][j];
				}
		return u;
	}
	public void putFromTopomaps (Topomap[] us, int k, double[][][] u)
	{
		Topomap ux = us[0], uy = us[1];
			for (int i = 0; i < ux.nx; i++)
				for (int j = 0; j < ux.ny; j++)
				{
					u[i][j][0] = ux.data[k][i][j];
					u[i][j][1] = uy.data[k][i][j];
				}
	}
	public Topomap[] getUTopomaps(String dir)
	{
		if (dir == null) dir = FileOps.selectDir(fc);
		Topomap ux, uy;
		ux = Topomap.readBIN(dir + "ux.bin");
		uy = Topomap.readBIN(dir + "uy.bin");
		return new Topomap[] {ux, uy};
	}
	public Layer[] getULayers(String dir)
	{
		if (dir == null) dir = FileOps.selectDir(fc);
		Layer ux, uy;
		ux = Layer.readBIN(dir + "ux.bin");
		uy = Layer.readBIN(dir + "uy.bin");
		return new Layer[] {ux, uy};
	}
	public Topomap[] getPhaseTopomaps(String dir)
	{
		if (dir == null) dir = FileOps.selectDir(fc);
		Topomap ux, uy;
		ux = Topomap.readBIN(dir + "phase_b0.bin");
		uy = Topomap.readBIN(dir + "phase_b1.bin");
		return new Topomap[] {ux, uy};
	}
	public Topomap[] getULattTopomaps(String dir)
	{
		if (dir == null) dir = FileOps.selectDir(fc);
		Topomap ux, uy;
		ux = Topomap.readBIN(dir + "uLatt_a.bin");
		uy = Topomap.readBIN(dir + "uLatt_b.bin");
		return new Topomap[] {ux, uy};
	}
	public void convertUFromPixToLatt(String dir)
	{
		Topomap[] us = getUTopomaps(dir);
		Topomap[] uLatts = new Topomap [2];
		AtomicCoordinatesSet latt = loadFromDir(dir);
		double[][][] u = new double [us[0].nx][us[1].ny][2];
		double[][][][] uLattTopo = new double[2][us[0].nlayers][us[0].nx][us[1].ny];
		double[][][] uLatt = new double [u.length][u[0].length][2];
		
		for (int k = 0; k < us[0].nlayers; k++)
		{
			putFromTopomaps(us, k, u);
			putConvertedToLatticeCoordinates(u, latt, uLatt);
			for (int i = 0; i < us[0].nx; i++)
				for (int j = 0; j < us[0].ny; j++)
				{
					uLattTopo[0][k][i][j] = uLatt[i][j][0];
					uLattTopo[1][k][i][j] = uLatt[i][j][1];
				}
		}
		uLatts[0] = Topomap.newTopomap(us[0], uLattTopo[0]);
		uLatts[1] = Topomap.newTopomap(us[1], uLattTopo[1]);
		writeLatticeTopos(uLatts, dir);
	}
	public void writeLatticeTopos(Topomap[] uLatts, String dir)
	{
		Topomap.writeBIN(uLatts[0], dir + "uLatt_a.bin");
		Topomap.writeBIN(uLatts[1], dir + "uLatt_b.bin");
	}
	public void writeUTopos(Topomap[] u, String dir)
	{
		Topomap.writeBIN(u[0], dir + "ux.bin");
		Topomap.writeBIN(u[1], dir + "uy.bin");
	}
	public void subtractParabolicFit(String dir)
	{
		File outputDir = new File(dir + "analysis\\");
		if (!outputDir.exists()) outputDir.mkdir();
		Topomap[] u = getUTopomaps(dir);
		TopomapUtil.subtractParabolicFit(u);
		writeUTopos(u, outputDir.toString() + "\\");
	}
	public void subtractParabolicFitLatt(String dir)
	{
		File outputDir = new File(dir + "analysis\\");
		if (!outputDir.exists()) outputDir.mkdir();
		Topomap[] uLatt = getULattTopomaps(dir);
		TopomapUtil.subtractParabolicFit(uLatt);
		writeLatticeTopos(uLatt, outputDir.toString() + "\\");
	}
	public int readSelectedLayer(String dir)
	{
		try{
			int x = Integer.parseInt(ColumnIO.readString(dir + "SelectedLayer.txt").trim());
			return x;
		}
		catch(Exception e)
		{
			int y = Integer.parseInt(JOptionPane.showInputDialog("File not found. What is the index of the layer you choose?"));
			ColumnIO.writeString("" + y, dir + "SelectedLayer.txt");
			return y;
		}
	}
	public void exportSelectedLayer(String dir)
	{
		int i = readSelectedLayer(dir);
		Topomap after = Topomap.readBIN(dir + "after.bin");
		Topomap[] u = getUTopomaps(dir);
		Topomap[] uLatt = getULattTopomaps(dir);
		Topomap[] phase = getPhaseTopomaps(dir);
		
		String outDir = dir + "Selected Layer\\";
		if (!new File(outDir).exists()) new File(outDir).mkdir();
		Layer.writeBIN(u[0].getLayer(i), outDir + "ux.bin");
		Layer.writeBIN(u[1].getLayer(i), outDir + "uy.bin");
		Layer.writeBIN(after.getLayer(i), outDir + "after.bin");
		Layer.writeBIN(uLatt[0].getLayer(i), outDir + "uLatt_a.bin");
		Layer.writeBIN(uLatt[1].getLayer(i), outDir + "uLatt_b.bin");
		
		Layer[] phaseBare = new Layer[] {phase[0].getLayer(i), phase[1].getLayer(i)};
		
		Layer.writeBIN(phaseBare[0], outDir + "phase_b0.bin");
		Layer.writeBIN(phaseBare[1], outDir + "phase_b1.bin");
		
		int[][] temp = new int [after.nx][after.ny];
		for (int j = 0; j < 2; j++)
		{
			FieldOps.putPhaseSteps(phaseBare[j].data, after.nx/2, after.ny/2, temp);
			FieldOps.putAddedPhaseSteps(phaseBare[j].data, temp, 2*Math.PI, phaseBare[j].data);
			Layer.writeBIN(phaseBare[j], outDir + "phaseCont_b" + j + ".bin");
		}
		
	}
	public void exportSelectedLayerAnalysis(String dir, int whatToDo)
	{
		
		int i = readSelectedLayer(dir);
		String adir = dir + "analysis\\";
		Topomap after = Topomap.readBIN(dir + "after.bin");
		Topomap[] u = null, uLatt = null;
		if (whatToDo == 0 || whatToDo == 2)
			u = getUTopomaps(adir);
		if (whatToDo == 1 || whatToDo == 2)
			uLatt = getULattTopomaps(adir);
		
		Topomap[] phase = getPhaseTopomaps(dir);
		
		String outDir = adir + "Selected Layer\\";
		if (!new File(outDir).exists()) new File(outDir).mkdir();
		if (whatToDo == 0 || whatToDo == 2)
		{
			Layer.writeBIN(u[0].getLayer(i), outDir + "ux.bin");
			Layer.writeBIN(u[1].getLayer(i), outDir + "uy.bin");
		}
		if (whatToDo == 1 || whatToDo == 2)
		{
			Layer.writeBIN(uLatt[0].getLayer(i), outDir + "uLatt_a.bin");
			Layer.writeBIN(uLatt[1].getLayer(i), outDir + "uLatt_b.bin");
		}
		Layer.writeBIN(after.getLayer(i), outDir + "after.bin");
		Layer.writeBIN(phase[0].getLayer(i), outDir + "phase_b0.bin");
		Layer.writeBIN(phase[1].getLayer(i), outDir + "phase_b1.bin");
	}
	public static void main(String[] args)
	{
		Topomap.setStdDir();
//		Layer t = Layer.open();
		DriftCorrectionMethods doer = new DriftCorrectionMethods();

		int choice = Integer.parseInt(JOptionPane.showInputDialog(null, "Enter your choice: \r\n"
				+ "0 - Do AUTOMATIC drift correction\r\n"
				+ "1 - Do drift correction with user-specified mask and detailed image output\r\n"
				+ "2 - Do file operations to extract the selected layer of output (required detailed output)\r\n"
				+ "3 - Load a u-field and apply it to something\r\n"
				+ "4 - Load an auto-U-field and apply it to a dI/dV map.\r\n"
				+ "5 - Apply a drift correction topomap to a topomap.\r\n"
				+ "6 - Load an auto-U-field and apply it to a single layer (FANCY).\r\n"
				+ "7 - Get the U-field as a function of length scale\r\n"
				+ "8 - Do something to measure the strain in real space\r\n"
				+ "         (for k-space use LayerViewerFourier_2)",
				"Drift Correction Options", JOptionPane.INFORMATION_MESSAGE));
//				(null, "To perform drift correction, click Yes.\r\n" +
//				"To do file operations on an output folder, click No.\r\n" +
//				"To open a u-field and apply it to something, click Cancel.");
		if (choice == 0)
		{
			doer.doingItAutomatic = 1;
			doer.getBraggPeaksFromUser();
		}
		else if (choice == 1)
			doer.getBraggPeaksFromUser();
		else if (choice == 2)
		{
			String dir = FileOps.selectDir(null).toString();
			doer.convertUFromPixToLatt(dir);
			doer.subtractParabolicFitLatt(dir);
			doer.subtractParabolicFit(dir);
			doer.exportSelectedLayer(dir);
			doer.exportSelectedLayerAnalysis(dir, 2);
		}
		else if (choice == 3)
		{
			double[][][] u = //doer.loadUFromLayers(doer.getULayers(null));
					DriftCorrectionMethods.getULayerFromFolder(null, doer.fc);
			int subchoice = JOptionPane.showConfirmDialog(null, "To apply U to a 3D map, click Yes. To apply to a single layer, click No.");
			if (subchoice == JOptionPane.YES_OPTION)
			{
				Topomap map = Topomap.open(doer.fc);
				File save = FileOps.selectSave(doer.fc);
				int secEstimate = (int)(Math.pow(map.nx/192.0, 4) * 1 * FANCY_TIME_192x192);
				boolean fancy = JOptionPane.showConfirmDialog(null, "Attempt fancy application of the u-field?\r\nFancy way estimated to take " + (secEstimate/60) + " minutes per layer.") == JOptionPane.YES_OPTION;
				Topomap done = TopomapUtil.applyUField(u, map);
				Topomap.writeBIN(done, save.toString());
				new TopomapViewer(done, doer.sourceDir, 512);
			}
			else if (subchoice == JOptionPane.NO_OPTION)
			{
				Layer l = Layer.open(doer.fc);
				double[][] newdat = l.data.clone();
				File save = FileOps.selectSave(doer.fc);
				int secEstimate = (int)(Math.pow(l.nx/192.0, 4) * 1 * FANCY_TIME_192x192);
				boolean fancy = JOptionPane.showConfirmDialog(null, "Attempt fancy application of the u-field?\r\nFancy way estimated to take " + (secEstimate/60) + " minutes per layer.") == JOptionPane.YES_OPTION;
				newdat = applyUField(fancy, newdat, u);
//				FieldOps.applyUFieldBiCubic(l.data, u, newdat);
				Layer done = Layer.newLayer(l, newdat);
				Layer.writeBIN(done, save.toString());
				new LayerViewer(done, doer.sourceDir, 512);
			}
		}
		else if (choice == 4)
		{
			JOptionPane.showMessageDialog(null, "Open the map");
			doer.applyAutoUFieldToTopomap(Topomap.open(doer.fc));
			System.exit(0);
		}
		else if (choice == 5)
		{
//			doer.applyAutoUFieldTopomapToTopomap(Topomap.open(doer.fc));
			doer.applyAutoUFieldTopomapToDirectory(FileOps.selectDir(doer.fc));
			System.exit(0);
		}
		else if (choice == 6)
		{		

			JOptionPane.showMessageDialog(null, "Open the drift correction folder.");
			double[][][] u = //DriftCorrectionAnalysis.getUFromAutoDriftCorrFolder(null, doer.fc);
					DriftCorrectionMethods.getULayerFromFolder(null, doer.fc);
			Layer l = Layer.open(doer.fc);
			double nxOver192 = l.nx/192.0;
			int secEstimate = (int)(Math.pow(nxOver192, 4) * FANCY_TIME_192x192);
			
			JOptionPane.showMessageDialog(null, "Calculation estimated to take " + (secEstimate/60) + " minutes.");
			double[][] result = applyUFieldSpecial(u, l.data);
			boolean[][] outsidePixels = FieldOps.getOutsidePixels(u);
			FieldOps.zero(result, outsidePixels);
			FieldOps.changeZeroToAverage(result);
			Layer answer = Layer.newLayer(l, result);
			
			String fn = doer.fc.getSelectedFile().getName();
			String dir = doer.fc.getCurrentDirectory().toString() + "\\";
			String fns = fn.substring(0, fn.length()-4);
			Layer.writeBIN(answer, dir + fns + "_afterFancy.bin");
		}
		if (choice == 7)
		{
			doer.doingItAutomatic = 2;
			doer.getBraggPeaksFromUser();
		}
		if (choice == 8)
		{
			doer.doingItAutomatic = 3;
			doer.getBraggPeaksFromUser();
		}
		if (choice == 9)
		{
			Topomap topography = Topomap.open(doer.fc);
			writeNewTranslationsToAFolder(topography, FileOps.selectDir(doer.fc));
			System.exit(0);
		}

//		
//		AtomicCoordinatesSet
		//		doer.getBraggPeaksFromUser(Layer.openFree(doer.fc));
		
//		Printer.printlnHorizontal(doer.bragg[0]);
//		Printer.printlnHorizontal(doer.bragg[1]);
	}
	 public static class BraggVectorRegularizer{
			
		public double angle = Math.PI/2;

		public int N;
		
		double[][] rotMat = new double [2][2];
		double[][] rotMatT = new double [2][2];
		
		double[][][] braggSetRot = new double [2][2][2];
		double[][][] unitBraggSetRot = new double [2][2][2];
		double[][] finalBragg = new double [2][2];
		
		boolean regEven = false; 
		
		//[0] is the first user-defined bragg peak and its +angle cousin. [1] = second user-defined Bragg peak and its -angle cousin 
		public BraggVectorRegularizer(double angle,
				int N, boolean regEven) {
			super();
			this.regEven = regEven;
			this.angle = angle;
			this.N = N;
			Matrix.putRotationMatrix(angle, rotMat);
			Matrix.putRotationMatrix(-angle, rotMatT);
		}
		
		public void setFirstBraggPeak(double[] bragg)
		{
				for (int j = 0; j < 2; j++)
					braggSetRot[0][0][j] = bragg[j];
				braggSetRot[0][1] = Matrix.getProductWith(rotMat, bragg);
		}
		public void setSecondBraggPeak(double[] bragg)
		{
				for (int j = 0; j < 2; j++)
					braggSetRot[1][1][j] = bragg[j];
				braggSetRot[1][0] = Matrix.getProductWith(rotMatT, bragg);
		}
		
		public void regularize()
		{
			double mag1 = Complex.mag(braggSetRot[0][0]);
			double mag2 = Complex.mag(braggSetRot[1][1]);
			double magAvg = (mag1+mag2)/2;
			
			System.out.println(mag1 + "\t" + mag2 + "\t" + magAvg);
			for (int i = 0; i < 2; i++)
				for (int j = 0; j < 2; j++)
					unitBraggSetRot[i][j] = Distance.unitVector(braggSetRot[i][j]);
			
			double angle1 = Complex.phase(unitBraggSetRot[0][0]);
			double angle2 = Complex.phase(unitBraggSetRot[1][0]);
			//In the case where angle1 or angle2 is close to 0, the other may be on the other side of 2pi. In that case special care must be taken:
			if (Math.abs(angle1-angle2) > Math.PI)
				if (angle1 > Math.PI) angle1 -= 2*Math.PI;
				else if (angle2 > Math.PI) angle2 -= 2*Math.PI;
			
			double angleAvg = (angle1+angle2)/2;
			System.out.println(angle1 +"\t"+ angle2 + "\t"+angleAvg);
			
			finalBragg[0][0] = magAvg*Math.cos(angleAvg);
			finalBragg[0][1] = magAvg*Math.sin(angleAvg);
			finalBragg[1][0] = magAvg*Math.cos(angleAvg+angle);
			finalBragg[1][1] = magAvg*Math.sin(angleAvg+angle);
			double[][] temp = new double [2][2];
			
			
			if (regEven){
				temp[0] = roundBraggIntEven(finalBragg[0], N);
				temp[1] = roundBraggIntEven(finalBragg[1], N);
			}
			else
			{
				temp[0] = roundBraggInt(finalBragg[0], N);
				temp[1] = roundBraggInt(finalBragg[1], N);
			}
//			parent.braggTrue = finalBragg;
		}
	 	public AtomicCoordinatesSet getRegLattice()
	 	{
			return AtomicCoordinatesSet.generateCentered(finalBragg, N);
	 	}
	 	public int[][] getFinalBragg()
	 	{
	 		int[][] ans = new int [2][2];
	 		for (int i = 0; i < ans.length; i++)
	 			for (int j = 0; j < 2; j++)
	 				ans[i][j] = FieldOps.round(finalBragg[i][j]);
	 		return ans;
	 	}
		//this returns the temp (pixel) variable for record-keeping purposes
		public static double[] roundBraggInt(double[] vec, int N)
		{
			double[] temp = new double[] {vec[0]/(2*Math.PI/N), vec[1]/(2*Math.PI/N)};
			temp[0] = FieldOps.round(temp[0]);
			temp[1] = FieldOps.round(temp[1]);
//			System.out.println(Printer.vectorP(temp));
			vec[0] = temp[0];//*(2*Math.PI/N);
			vec[1] = temp[1];//*(2*Math.PI/N);
			return temp; 
			
		}
		public static double[] roundBraggIntEven(double[] vec, int N)
		{
			double[] temp = new double[] {vec[0]/(2*Math.PI/N), vec[1]/(2*Math.PI/N)};
			temp[0] = FieldOps.roundEven(temp[0]);
			temp[1] = FieldOps.roundEven(temp[1]);
//			System.out.println(Printer.vectorP(temp));
			vec[0] = temp[0];//*(2*Math.PI/N);
			vec[1] = temp[1];//*(2*Math.PI/N);
			return temp; 
			
		}
	 }
	 
	 /**
	  * This program applies the u-field in the manner described in Ilija's IDL file from Harvard.
	  * That is, it takes the complex Fourier-transform of the data, but using a distorted position vector r-u.
	  * (Of course the distortion of space makes it impossible to use FFT in taking this, so that the method may take an extremely long time.)
	  * The inverse-FFT of that object is the final answer.
	  * 
	  * @param u
	  * @param data
	  * @return
	  */
	 public static double[][] applyUFieldSpecial(double[][][] u, double[][] data)
	 {
		 int nx = data.length, ny = data[0].length;
		 double[][][] fftz = new double [nx][ny][2];
		 
		 double[][] source = FieldOps.copy(data);
		 double mean = FieldOps.mean(source);
		 FieldOps.subtractAvg(source);
		 double qx, qy;
		 long start = System.currentTimeMillis();
		 double[][][] rMinusU = new double [nx][ny][2];
		 
		 for (int i = 0; i < nx; i++)
			 for (int j = 0; j < ny; j++)
			 {
				 rMinusU[i][j][0] = i + u[i][j][0];
				 rMinusU[i][j][1] = j + u[i][j][1];
			 }
		 
		 
		 //We will now try to save time by computing the sines and cosines outside of the for loops:
		 
		 long desiredPrecomputedDoubles = 2*((long)nx*(long)ny)*((long)nx + (long)ny); //This is the amount of space necessary to pre-compute the sin and cos of q dot rMinusU;
		 
		 System.out.println(maxPrecomputedDoubles + "\t" + desiredPrecomputedDoubles + "\t" + (desiredPrecomputedDoubles/maxPrecomputedDoubles));
		 
		 
		 //In this case the path is simple and we need only precalculate ALL the phases.
		 if (maxPrecomputedDoubles >= desiredPrecomputedDoubles) //pre-compute the entire pair of phase arrays and then apply trigonometric identities:
		 {
			 double[][][] cosXPart = new double [nx][nx][ny];
			 double[][][] sinXPart = new double [nx][nx][ny];
			 double[][][] cosYPart = new double [ny][nx][ny];
			 double[][][] sinYPart = new double [ny][nx][ny];
			 
			 for (int i = 0; i < nx; i++)
				 for (int j = 0; j < ny; j++)
				 {
					 for (int k = 0; k < nx; k++)
					 {
						qx = (k-nx/2)*2*Math.PI/nx;
						cosXPart[k][i][j] = Math.cos(qx*rMinusU[i][j][0]);
						sinXPart[k][i][j] = Math.sin(qx*rMinusU[i][j][0]);
					 }
					 for (int k = 0; k < ny; k++)
					 {
						qy = (k-ny/2)*2*Math.PI/ny;
						cosYPart[k][i][j] = Math.cos(qy*rMinusU[i][j][1]);
						sinYPart[k][i][j] = Math.sin(qy*rMinusU[i][j][1]);
					 }
				 }
			 
			 //now let us run through the loop as described in "old":
			 for (int i = 0; i < nx; i++)
				 for (int j = 0; j < ny; j++)
				 {
					 if (rMinusU[i][j][0] < 0 || rMinusU[i][j][0] > nx-1 || rMinusU[i][j][1] < 0 || rMinusU[i][j][1] > ny-1)
						 source[i][j] = 0;
				 }
				int im, jm; //the coordinates in the array of the vector -Q.

			 for (int i = 0; i <= nx/2; i++){System.out.print(" " + i);
				 for (int j = 0; j < ny; j++)
				 {
					 im = -i + nx;
					 jm = -j + ny;
					 if (i != 0 && j != 0 && i != im){
						 for (int p = 0; p < nx; p++)
							 for (int q = 0; q < ny; q++)
							 {
								 //The cosine of q dot R is cos (qxrx + qyry) and according to the sum of angles identity is below:
								 fftz[i][j][0] += source[p][q]*(cosXPart[i][p][q]*cosYPart[j][p][q] - sinXPart[i][p][q]*sinYPart[j][p][q]); 
								 fftz[i][j][1] += -source[p][q]*(sinXPart[i][p][q]*cosYPart[j][p][q] + cosXPart[i][p][q]*sinYPart[j][p][q]);
							 }
						 fftz[im][jm][0] = fftz[i][j][0];
						 fftz[im][jm][1] = -fftz[i][j][1];
					 }
					 else
						 for (int p = 0; p < nx; p++)
							 for (int q = 0; q < ny; q++)
							 {
								 fftz[i][j][0] += source[p][q]*(cosXPart[i][p][q]*cosYPart[j][p][q] - sinXPart[i][p][q]*sinYPart[j][p][q]); 
								 fftz[i][j][1] += -source[p][q]*(sinXPart[i][p][q]*cosYPart[j][p][q] + cosXPart[i][p][q]*sinYPart[j][p][q]);
							 }
				 }
			 }
				 double[][][] ansZ = new double [nx][ny][2];
				 
				 FFTOps.putIFFT(fftz, ansZ, true);
				 long stop = System.currentTimeMillis();
				 
				 System.out.println();
				 System.out.println("It took " + (stop-start)/1000 + " seconds.");
				 
				 double[][] ans = FieldOps.getIndex(ansZ, 0);
				 FieldOps.plusEquals(ans, mean);
				 return ans;
		 }
		 else //if on the other hand pre-computing everything would take too much memory, we must proceed by blocks.
		 {
			 //The block size must evenly divide the map size. If the map size is a prime number, screw you.
			 int maxBlockSize = (int)(maxPrecomputedDoubles/(nx*ny*4));
			 System.out.println("max size is " + maxBlockSize);
			 //Now, shrink the block size until it is an integer divisor of the map size:
			 while (maxBlockSize*(nx/maxBlockSize) != nx)
				 maxBlockSize--;
			 System.out.println("max size is " + maxBlockSize);
			 
			 //Now, define the doubles:
			 double[][][] cosXPart = new double [maxBlockSize][nx][ny];
			 double[][][] sinXPart = new double [maxBlockSize][nx][ny];
			 double[][][] cosYPart = new double [maxBlockSize][nx][ny];
			 double[][][] sinYPart = new double [maxBlockSize][nx][ny];
			 
			 for (int i = 0; i < nx; i++)
				 for (int j = 0; j < ny; j++)
				 {
					 if (rMinusU[i][j][0] < 0 || rMinusU[i][j][0] > nx-1 || rMinusU[i][j][1] < 0 || rMinusU[i][j][1] > ny-1)
						 source[i][j] = 0;
				 }
				int im, jm; //the coordinates in the array of the vector -Q.

			 int inb, jnb;
			 int nBlocks = nx/maxBlockSize;
			 for (int bi = 0; bi < nBlocks; bi++)
				 for (int bj = 0; bj < nBlocks; bj++)
				 {
					 System.out.println("\r\n" + bi + "\t" + bj);
					 //The following for loops will be USED only if bi*maxBlockSize <= nx/2, since the N^4 part only runs if i <= nx/2
					 if (bi*maxBlockSize <= nx/2) 
						 for (int p = 0; p < nx; p++)
							 for (int q = 0; q < ny; q++)
							 {
								 for (int i = bi*maxBlockSize; i < (bi+1)*maxBlockSize; i++)
								 {
									inb = i % maxBlockSize;
									qx = (i-nx/2)*2*Math.PI/nx;
									cosXPart[inb][p][q] = Math.cos(qx*rMinusU[p][q][0]);
									sinXPart[inb][p][q] = Math.sin(qx*rMinusU[p][q][0]);
								 }
								 for (int j = bj*maxBlockSize; j < (bj+1)*maxBlockSize; j++)
								 {
									 jnb = j % maxBlockSize;
									qy = (j-ny/2)*2*Math.PI/ny;
									cosYPart[jnb][p][q] = Math.cos(qy*rMinusU[p][q][1]);
									sinYPart[jnb][p][q] = Math.sin(qy*rMinusU[p][q][1]);
								 }
							 }
					 //having defined the block of pre-computed sines and cosines we proceed to
					 //compute block-wise the fftz as before:
						 for (int i = bi*maxBlockSize; i < (bi+1)*maxBlockSize && i <= nx/2; i++){//System.out.print(" " + i);
							 for (int j = bj*maxBlockSize; j < (bj+1)*maxBlockSize; j++)
							 {
								 im = -i + nx;
								 jm = -j + ny;
								 inb = i % maxBlockSize;
								 jnb = j % maxBlockSize;
								 if (i != 0 && j != 0 && i != im){
									 for (int p = 0; p < nx; p++)
										 for (int q = 0; q < ny; q++)
										 {
											 //The cosine of q dot R is cos (qxrx + qyry) and according to the sum of angles identity is below:
											 fftz[i][j][0] += source[p][q]*(cosXPart[inb][p][q]*cosYPart[jnb][p][q] - sinXPart[inb][p][q]*sinYPart[jnb][p][q]); 
											 fftz[i][j][1] += -source[p][q]*(sinXPart[inb][p][q]*cosYPart[jnb][p][q] + cosXPart[inb][p][q]*sinYPart[jnb][p][q]);
										 }
									 fftz[im][jm][0] = fftz[i][j][0];
									 fftz[im][jm][1] = -fftz[i][j][1];
								 }
								 else
									 for (int p = 0; p < nx; p++)
										 for (int q = 0; q < ny; q++)
										 {
											 fftz[i][j][0] += source[p][q]*(cosXPart[inb][p][q]*cosYPart[jnb][p][q] - sinXPart[inb][p][q]*sinYPart[jnb][p][q]); 
											 fftz[i][j][1] += -source[p][q]*(sinXPart[inb][p][q]*cosYPart[jnb][p][q] + cosXPart[inb][p][q]*sinYPart[jnb][p][q]);
										 }
							 }
						 }
				 }
			 double[][][] ansZ = new double [nx][ny][2];
			 
			 FFTOps.putIFFT(fftz, ansZ, true);
			 long stop = System.currentTimeMillis();
			 
			 System.out.println();
			 System.out.println("It took " + (stop-start)/1000 + " seconds.");
			 
			 double[][] ans = FieldOps.getIndex(ansZ, 0);
			 FieldOps.plusEquals(ans, mean);
			 return ans;
		 }
		 
		 
		 
		 //Because this method takes a Fourier transform, it gives "periodic" results in the sense that stuff which is pushed over one edge will reappear on the other edge.
		 //All pixels which are going to be pushed outside the edge have to be destroyed before the algorithm begins.
	 }
	 
	 public static double[][] applyUFieldSpecial_withShifting(double[][][] u, double[][] data, int[][] bragg, double[] placeToPutTranslation)
	 {
		 int nx = data.length, ny = data[0].length;
		 double[][][] fftz = new double [nx][ny][2];
//		 if (braggReserve == null) throw new NullPointerException();
		 if (bragg == null) throw new NullPointerException();
		 double[][] source = FieldOps.copy(data);
		 double mean = FieldOps.mean(source);
		 FieldOps.subtractAvg(source);
		 double qx, qy;
		 long start = System.currentTimeMillis();
		 double[][][] rMinusU = new double [nx][ny][2];
		 
		 for (int i = 0; i < nx; i++)
			 for (int j = 0; j < ny; j++)
			 {
				 rMinusU[i][j][0] = i + u[i][j][0];
				 rMinusU[i][j][1] = j + u[i][j][1];
			 }
		 
		 
		 //We will now try to save time by computing the sines and cosines outside of the for loops:
		 
		 long desiredPrecomputedDoubles = 2*((long)nx*(long)ny)*((long)nx + (long)ny); //This is the amount of space necessary to pre-compute the sin and cos of q dot rMinusU;
		 
		 System.out.println(maxPrecomputedDoubles + "\t" + desiredPrecomputedDoubles + "\t" + (desiredPrecomputedDoubles/maxPrecomputedDoubles));
		 
		 
		 //In this case the path is simple and we need only precalculate ALL the phases.
		 if (maxPrecomputedDoubles >= desiredPrecomputedDoubles) //pre-compute the entire pair of phase arrays and then apply trigonometric identities:
		 {
			 double[][][] cosXPart = new double [nx][nx][ny];
			 double[][][] sinXPart = new double [nx][nx][ny];
			 double[][][] cosYPart = new double [ny][nx][ny];
			 double[][][] sinYPart = new double [ny][nx][ny];
			 
			 for (int i = 0; i < nx; i++)
				 for (int j = 0; j < ny; j++)
				 {
					 for (int k = 0; k < nx; k++)
					 {
						qx = (k-nx/2)*2*Math.PI/nx;
						cosXPart[k][i][j] = Math.cos(qx*rMinusU[i][j][0]);
						sinXPart[k][i][j] = Math.sin(qx*rMinusU[i][j][0]);
					 }
					 for (int k = 0; k < ny; k++)
					 {
						qy = (k-ny/2)*2*Math.PI/ny;
						cosYPart[k][i][j] = Math.cos(qy*rMinusU[i][j][1]);
						sinYPart[k][i][j] = Math.sin(qy*rMinusU[i][j][1]);
					 }
				 }
			 
			 //now let us run through the loop as described in "old":
			 for (int i = 0; i < nx; i++)
				 for (int j = 0; j < ny; j++)
				 {
					 if (rMinusU[i][j][0] < 0 || rMinusU[i][j][0] > nx-1 || rMinusU[i][j][1] < 0 || rMinusU[i][j][1] > ny-1)
						 source[i][j] = 0;
				 }
				int im, jm; //the coordinates in the array of the vector -Q.

			 for (int i = 0; i <= nx/2; i++){System.out.print(" " + i);
				 for (int j = 0; j < ny; j++)
				 {
					 im = -i + nx;
					 jm = -j + ny;
					 if (i != 0 && j != 0 && i != im){
						 for (int p = 0; p < nx; p++)
							 for (int q = 0; q < ny; q++)
							 {
								 //The cosine of q dot R is cos (qxrx + qyry) and according to the sum of angles identity is below:
								 fftz[i][j][0] += source[p][q]*(cosXPart[i][p][q]*cosYPart[j][p][q] - sinXPart[i][p][q]*sinYPart[j][p][q]); 
								 fftz[i][j][1] += -source[p][q]*(sinXPart[i][p][q]*cosYPart[j][p][q] + cosXPart[i][p][q]*sinYPart[j][p][q]);
							 }
						 fftz[im][jm][0] = fftz[i][j][0];
						 fftz[im][jm][1] = -fftz[i][j][1];
					 }
					 else
						 for (int p = 0; p < nx; p++)
							 for (int q = 0; q < ny; q++)
							 {
								 fftz[i][j][0] += source[p][q]*(cosXPart[i][p][q]*cosYPart[j][p][q] - sinXPart[i][p][q]*sinYPart[j][p][q]); 
								 fftz[i][j][1] += -source[p][q]*(sinXPart[i][p][q]*cosYPart[j][p][q] + cosXPart[i][p][q]*sinYPart[j][p][q]);
							 }
				 }
			 }
		 }
		 else //if on the other hand pre-computing everything would take too much memory, we must proceed by blocks.
		 {
			 //The block size must evenly divide the map size. If the map size is a prime number, screw you.
			 int maxBlockSize = (int)(maxPrecomputedDoubles/(nx*ny*4));
			 System.out.println("max size is " + maxBlockSize);
			 //Now, shrink the block size until it is an integer divisor of the map size:
			 while (maxBlockSize*(nx/maxBlockSize) != nx)
				 maxBlockSize--;
			 System.out.println("max size is " + maxBlockSize);
			 
			 //Now, define the doubles:
			 double[][][] cosXPart = new double [maxBlockSize][nx][ny];
			 double[][][] sinXPart = new double [maxBlockSize][nx][ny];
			 double[][][] cosYPart = new double [maxBlockSize][nx][ny];
			 double[][][] sinYPart = new double [maxBlockSize][nx][ny];
			 
			 for (int i = 0; i < nx; i++)
				 for (int j = 0; j < ny; j++)
				 {
					 if (rMinusU[i][j][0] < 0 || rMinusU[i][j][0] > nx-1 || rMinusU[i][j][1] < 0 || rMinusU[i][j][1] > ny-1)
						 source[i][j] = 0;
				 }
				int im, jm; //the coordinates in the array of the vector -Q.

			 int inb, jnb;
			 int nBlocks = nx/maxBlockSize;
			 for (int bi = 0; bi < nBlocks; bi++)
				 for (int bj = 0; bj < nBlocks; bj++)
				 {
					 System.out.println("\r\n" + bi + "\t" + bj);
					 //The following for loops will be USED only if bi*maxBlockSize <= nx/2, since the N^4 part only runs if i <= nx/2
					 if (bi*maxBlockSize <= nx/2) 
						 for (int p = 0; p < nx; p++)
							 for (int q = 0; q < ny; q++)
							 {
								 for (int i = bi*maxBlockSize; i < (bi+1)*maxBlockSize; i++)
								 {
									inb = i % maxBlockSize;
									qx = (i-nx/2)*2*Math.PI/nx;
									cosXPart[inb][p][q] = Math.cos(qx*rMinusU[p][q][0]);
									sinXPart[inb][p][q] = Math.sin(qx*rMinusU[p][q][0]);
								 }
								 for (int j = bj*maxBlockSize; j < (bj+1)*maxBlockSize; j++)
								 {
									 jnb = j % maxBlockSize;
									qy = (j-ny/2)*2*Math.PI/ny;
									cosYPart[jnb][p][q] = Math.cos(qy*rMinusU[p][q][1]);
									sinYPart[jnb][p][q] = Math.sin(qy*rMinusU[p][q][1]);
								 }
							 }
					 //having defined the block of pre-computed sines and cosines we proceed to
					 //compute block-wise the fftz as before:
						 for (int i = bi*maxBlockSize; i < (bi+1)*maxBlockSize && i <= nx/2; i++){//System.out.print(" " + i);
							 for (int j = bj*maxBlockSize; j < (bj+1)*maxBlockSize; j++)
							 {
								 im = -i + nx;
								 jm = -j + ny;
								 inb = i % maxBlockSize;
								 jnb = j % maxBlockSize;
								 if (i != 0 && j != 0 && i != im){
									 for (int p = 0; p < nx; p++)
										 for (int q = 0; q < ny; q++)
										 {
											 //The cosine of q dot R is cos (qxrx + qyry) and according to the sum of angles identity is below:
											 fftz[i][j][0] += source[p][q]*(cosXPart[inb][p][q]*cosYPart[jnb][p][q] - sinXPart[inb][p][q]*sinYPart[jnb][p][q]); 
											 fftz[i][j][1] += -source[p][q]*(sinXPart[inb][p][q]*cosYPart[jnb][p][q] + cosXPart[inb][p][q]*sinYPart[jnb][p][q]);
										 }
									 fftz[im][jm][0] = fftz[i][j][0];
									 fftz[im][jm][1] = -fftz[i][j][1];
								 }
								 else
									 for (int p = 0; p < nx; p++)
										 for (int q = 0; q < ny; q++)
										 {
											 fftz[i][j][0] += source[p][q]*(cosXPart[inb][p][q]*cosYPart[jnb][p][q] - sinXPart[inb][p][q]*sinYPart[jnb][p][q]); 
											 fftz[i][j][1] += -source[p][q]*(sinXPart[inb][p][q]*cosYPart[jnb][p][q] + cosXPart[inb][p][q]*sinYPart[jnb][p][q]);
										 }
							 }
						 }
				 }
		 }
	 	//Now we have the fftz. Let us apply the proper translation operator, then ifft it.
	 	double[] r = FieldOps.getRForTopo(null, fftz, bragg);
	 	placeToPutTranslation[0] = r[0];
	 	placeToPutTranslation[1] = r[1];
//	 	r[0] = -r[0]; r[1] = -r[1];
	 	System.out.println("fftzO:\t" + ((bragg[0][0]+nx/2)%nx) + "\t" + ((bragg[0][1]+ny/2)%ny) + "\t"+ fftz[(bragg[0][0]+nx/2)%nx][(bragg[0][1]+ny/2)%ny][0] + "\t" + fftz[(bragg[0][0]+nx/2)%nx][(bragg[0][1]+ny/2)%ny][1]);
	 	System.out.println("fftzO:\t" + ((bragg[1][0]+nx/2)%nx) + "\t" + ((bragg[1][1]+ny/2)%ny) + "\t" + fftz[(bragg[1][0]+nx/2)%nx][(bragg[1][1]+ny/2)%ny][0] + "\t" + fftz[(bragg[1][0]+nx/2)%nx][(bragg[1][1]+ny/2)%ny][1]);
	 	double[][] ansre = FieldOps.translateFancy(null, fftz, r, false);
		FieldOps.plusEquals(ansre, mean);
	 	System.out.println("fftzN:\t" + fftz[(bragg[0][0]+nx/2)%nx][(bragg[0][1]+ny/2)%ny][0] + "\t" + fftz[(bragg[0][0]+nx/2)%nx][(bragg[0][1]+ny/2)%ny][1]);
	 	System.out.println("fftzN:\t" + fftz[(bragg[1][0]+nx/2)%nx][(bragg[1][1]+ny/2)%ny][0] + "\t" + fftz[(bragg[1][0]+nx/2)%nx][(bragg[1][1]+ny/2)%ny][1]);

		long stop = System.currentTimeMillis();
		 System.out.println();
		 System.out.println("It took " + (stop-start)/1000 + " seconds.");
		 return ansre;
		 //Because this method takes a Fourier transform, it gives "periodic" results in the sense that stuff which is pushed over one edge will reappear on the other edge.
		 //All pixels which are going to be pushed outside the edge have to be destroyed before the algorithm begins.
	 }
	 /**
	  * This is intended to define the translation that will be necessary when the U-field is fancily applied.
	  * The translation normally calculated from the fancy inverse FFT. But, actually,
	  * we only need two pixels of the fancy inverse FFT to be nonzero in order to make the calculation:
	  * Bragg1 and Bragg2. The method calcaultes the fancy inverse FFT at only those two pixels, then returns the translation.
	  * @param u
	  * @param data
	  * @param bragg
	  * @return
	  */
	 public static double[] getTranslationSpecialQuickDirty(double[][][] u, double[][] data, int[][] bragg)
	 {
		 int nx = data.length, ny = data[0].length;
		 double[][][] fftz = new double [nx][ny][2];
//		 if (braggReserve == null) throw new NullPointerException();
		 if (bragg == null) throw new NullPointerException();
		 double[][] source = FieldOps.copy(data);
		 double mean = FieldOps.mean(source);
		 FieldOps.subtractAvg(source);
		 double[] qx = new double [2], qy = new double [2];
		 long start = System.currentTimeMillis();
		 double[][][] rMinusU = new double [nx][ny][2];
		 
		 for (int i = 0; i < nx; i++)
			 for (int j = 0; j < ny; j++)
			 {
				 rMinusU[i][j][0] = i + u[i][j][0];
				 rMinusU[i][j][1] = j + u[i][j][1];
			 }
		 
		 
		 //We will now try to save time by computing the sines and cosines outside of the for loops:
		 
		 double[][][] cosXPart = new double [2][nx][ny];
		 double[][][] sinXPart = new double [2][nx][ny];
		 double[][][] cosYPart = new double [2][nx][ny];
		 double[][][] sinYPart = new double [2][nx][ny];
			 
			 for (int i = 0; i < nx; i++)
				 for (int j = 0; j < ny; j++)
				 {
					 for (int k = 0; k < 2; k++)
					 {
						qx[k] = bragg[k][0]*2*Math.PI/nx;
						cosXPart[k][i][j] = Math.cos(qx[k]*rMinusU[i][j][0]);
						sinXPart[k][i][j] = Math.sin(qx[k]*rMinusU[i][j][0]);
					 }
					 for (int k = 0; k < 2; k++)
					 {
						qy[k] = bragg[k][1]*2*Math.PI/ny;
						cosYPart[k][i][j] = Math.cos(qy[k]*rMinusU[i][j][1]);
						sinYPart[k][i][j] = Math.sin(qy[k]*rMinusU[i][j][1]);
					 }
				 }
			 
			 //now let us run through the loop as described in "old":
			 for (int i = 0; i < nx; i++)
				 for (int j = 0; j < ny; j++)
				 {
					 if (rMinusU[i][j][0] < 0 || rMinusU[i][j][0] > nx-1 || rMinusU[i][j][1] < 0 || rMinusU[i][j][1] > ny-1)
						 source[i][j] = 0;
				 }
				int im, jm; //the coordinates in the array of the vector -Q.

			for (int k = 0; k < 2; k++){
				int i = (bragg[k][0]+nx/2)%nx, j = (bragg[k][1]+ny/2)%ny;
				 im = -i + nx;
				 jm = -j + ny;
				 if (i != 0 && j != 0 && i != im){
					 for (int p = 0; p < nx; p++)
						 for (int q = 0; q < ny; q++)
						 {
							 //The cosine of q dot R is cos (qxrx + qyry) and according to the sum of angles identity is below:
							 fftz[i][j][0] += source[p][q]*(cosXPart[k][p][q]*cosYPart[k][p][q] - sinXPart[k][p][q]*sinYPart[k][p][q]); 
							 fftz[i][j][1] += -source[p][q]*(sinXPart[k][p][q]*cosYPart[k][p][q] + cosXPart[k][p][q]*sinYPart[k][p][q]);
						 }
					 fftz[im][jm][0] = fftz[i][j][0];
					 fftz[im][jm][1] = -fftz[i][j][1];
				 }
				 else
					 for (int p = 0; p < nx; p++)
						 for (int q = 0; q < ny; q++)
						 {
							 fftz[i][j][0] += source[p][q]*(cosXPart[k][p][q]*cosYPart[k][p][q] - sinXPart[k][p][q]*sinYPart[k][p][q]); 
							 fftz[i][j][1] += -source[p][q]*(sinXPart[k][p][q]*cosYPart[k][p][q] + cosXPart[k][p][q]*sinYPart[k][p][q]);
						 }
				 }
	 	//Now we have the fftz. Let us apply the proper translation operator, then ifft it.
	 	double[] r = FieldOps.getRForTopo(null, fftz, bragg);
	 	System.out.println("fftzO:\t" + ((bragg[0][0]+nx/2)%nx) + "\t" + ((bragg[0][1]+ny/2)%ny) + "\t"+ fftz[(bragg[0][0]+nx/2)%nx][(bragg[0][1]+ny/2)%ny][0] + "\t" + fftz[(bragg[0][0]+nx/2)%nx][(bragg[0][1]+ny/2)%ny][1]);
	 	System.out.println("fftzO:\t" + ((bragg[1][0]+nx/2)%nx) + "\t" + ((bragg[1][1]+ny/2)%ny) + "\t" + fftz[(bragg[1][0]+nx/2)%nx][(bragg[1][1]+ny/2)%ny][0] + "\t" + fftz[(bragg[1][0]+nx/2)%nx][(bragg[1][1]+ny/2)%ny][1]);
	 	double[][] ansre = FieldOps.translateFancy(null, fftz, r, false);
		FieldOps.plusEquals(ansre, mean);
	 	System.out.println("fftzN:\t" + fftz[(bragg[0][0]+nx/2)%nx][(bragg[0][1]+ny/2)%ny][0] + "\t" + fftz[(bragg[0][0]+nx/2)%nx][(bragg[0][1]+ny/2)%ny][1]);
	 	System.out.println("fftzN:\t" + fftz[(bragg[1][0]+nx/2)%nx][(bragg[1][1]+ny/2)%ny][0] + "\t" + fftz[(bragg[1][0]+nx/2)%nx][(bragg[1][1]+ny/2)%ny][1]);

		long stop = System.currentTimeMillis();
		 System.out.println();
		 System.out.println("It took " + (stop-start)/1000 + " seconds.");

	 	return r;//	 	r[0] = -r[0]; r[1] = -r[1];
//		 return ansre;
		 //Because this method takes a Fourier transform, it gives "periodic" results in the sense that stuff which is pushed over one edge will reappear on the other edge.
		 //All pixels which are going to be pushed outside the edge have to be destroyed before the algorithm begins.
	 }
	 public static double[][] applyUFieldSpecial_Old(double[][][] u, double[][] data)
	 {
		 int nx = data.length, ny = data[0].length;
		 double[][][] fftz = new double [nx][ny][2];
		 
		 double[][] source = data.clone();
		 double mean = FieldOps.mean(source);
		 FieldOps.subtractAvg(source);
		 double qx, qy;
		 long start = System.currentTimeMillis();
		 double[][][] rMinusU = new double [nx][ny][2];
		 
		 for (int i = 0; i < nx; i++)
			 for (int j = 0; j < ny; j++)
			 {
				 rMinusU[i][j][0] = i + u[i][j][0];
				 rMinusU[i][j][1] = j + u[i][j][1];
			 }
		 
		 //Because this method takes a Fourier transform, it gives "periodic" results in the sense that stuff which is pushed over one edge will reappear on the other edge.
		 //All pixels which are going to be pushed outside the edge have to be destroyed before the algorithm begins.
		 for (int i = 0; i < nx; i++)
			 for (int j = 0; j < ny; j++)
			 {
				 if (rMinusU[i][j][0] < 0 || rMinusU[i][j][0] > nx-1 || rMinusU[i][j][1] < 0 || rMinusU[i][j][1] > ny-1)
					 source[i][j] = 0;
			 }
//		 FieldOps.changeZeroToAverage(source);
//		 double qDotRMinusU;
//		 for (int i = 0; i < nx; i++){System.out.print(" " + i); //the old loop went through the entire range
//			 for (int j = 0; j < ny; j++)
//			 {
//				 qx = (i-nx/2)*2*Math.PI/nx;
//				 qy = (j-ny/2)*2*Math.PI/ny;
//				 for (int p = 0; p < nx; p++)
//					 for (int q = 0; q < ny; q++)
//					 {
//						 qDotRMinusU = (qx*rMinusU[p][q][0] + qy*rMinusU[p][q][1]);
//						 fftz[i][j][0] += data[p][q]*Math.cos(qDotRMinusU);
//						 fftz[i][j][1] += -data[p][q]*Math.sin(qDotRMinusU);
//					 }
//			 }
			 double qDotRMinusU;
			 int im, jm; //the coordinates in the array of the vector -Q.
			 
			 for (int i = 0; i <= nx/2; i++){System.out.print(" " + i);
				 for (int j = 0; j < ny; j++)
				 {
					 qx = (i-nx/2)*2*Math.PI/nx;
					 qy = (j-ny/2)*2*Math.PI/ny;
					 im = -i + nx;
					 jm = -j + ny;
					 if (i != 0 && j != 0 && i != im)
						 for (int p = 0; p < nx; p++)
							 for (int q = 0; q < ny; q++)
							 {
								 qDotRMinusU = (qx*rMinusU[p][q][0] + qy*rMinusU[p][q][1]);
								 fftz[i][j][0] += source[p][q]*Math.cos(qDotRMinusU);
								 fftz[i][j][1] += -source[p][q]*Math.sin(qDotRMinusU);
								 fftz[im][jm][0] = fftz[i][j][0];
								 fftz[im][jm][1] = -fftz[i][j][1];
							 }
					 else
						 for (int p = 0; p < nx; p++)
							 for (int q = 0; q < ny; q++)
							 {
								 qDotRMinusU = (qx*rMinusU[p][q][0] + qy*rMinusU[p][q][1]);
								 fftz[i][j][0] += source[p][q]*Math.cos(qDotRMinusU);
								 fftz[i][j][1] += -source[p][q]*Math.sin(qDotRMinusU);
							 }
						 
				 }
		 }
		 
		 double[][] fftmag = FieldOps.magnitude(fftz);
//		 FieldOps.log(fftmag);
//		 LayerViewer.show(Layer.getFreeLayer(fftmag), 512);
		 double[][][] ansZ = new double [nx][ny][2];
		 
		 FFTOps.putIFFT(fftz, ansZ, true);
		 long stop = System.currentTimeMillis();
		 
		 System.out.println();
		 System.out.println("It took " + (stop-start)/1000 + " seconds.");
		 
		 double[][] ans = FieldOps.getIndex(ansZ, 0);
		 FieldOps.plusEquals(ans, mean);
		 return ans;
	 }
	 
	 /**
	  * This is designed to apply the masks of a series, pick the best u-field automatically, and apply that best one.
	  * The best one is determined at the level of phaseCont, so I allocate nlayers of phaseCont despite the memory requirements.
	  * 
	  * @author madhavanlab2011
	  *
	  */
	 public static class UFieldCalculation extends DriftCorrectionCalculator
	 {
		 	double[][][] expUfft;
			double[][][][] shiftedFFTZ;
			FourierFilterMaskSource mask;
			public UFieldCalculation(FourierFilterMaskSource mask, Layer t, BicubicSplineInterpolatingFunction interp, int[][] bragg, double angle, boolean forceEven)
			{
				super(t, interp, bragg, angle, mask.getNInLoop(), forceEven);
				this.mask = mask;
				shiftedFFTZ = new double[mask.getNInLoop()][t.nx][t.ny][2];
				expUfft = new double [t.nx][t.ny][2];
				//Do preliminary steps
				FFTOps.putFFT(t.data, fftz, false);
				for (int i = 0; i < 2; i++)
					FieldOps.shift(fftz, shiftedFFTZ[i], bragg[i][0], bragg[i][1]);
			}
			
			public void doExpuCalc()
			{
				for (int n = 0; n < 2; n++)
					for (int i = 0; i <mask.getNInLoop(); i++)
					{
						LayerUtil.putExpUFourier(shiftedFFTZ[n], mask.getMask(n, i), expU, expUfft);
						FieldOps.phase(expU, expUPhase[n]);
						FieldOps.putPhaseSteps(expUPhase[n], N/2, N/2, phaseN[n]);
						FieldOps.putAddedPhaseSteps(expUPhase[n], phaseN[n], 2*Math.PI, phaseCont[i][n]);
						FieldOps.subtractAvg(phaseCont[i][n]);
						//an additional step will be used to place an atom in the center of the field of view, so that the default origin {N/2, N/2} has an atom on top of it.
//						FieldOps.shiftCenterModulo(phaseCont[i][n], Math.PI*2);
					}
			}

			public String[] getEnergyLines()
			{
				String[] ans = new String[nlayers+2];
				ans[0] = "Index\tLength Scale\tQuality";
				for (int i = 1; i < nlayers+1; i++)
					ans[i] = "" + (i-1) + "\t" + mask.getLengthScale(i-1) + "\t" + energy[i-1];
				ans[nlayers+1] = "The selected layer was " + selectedLayer + ".";
				return ans;
			}
			/**
			 * This embodies the previous methods getFinished_NonlinearOnly, getFinished_LinearOnly, getFinished_Both.
			 * The code will be 0 for nonlinear only, 1 for both, 2 for linear only. I expect only 0 and 1 will be used.
			 * @param mask
			 * @param t
			 * @param bragg
			 * @param fancy
			 * @param forceEven
			 * @param code
			 * @return
			 */
			public static UFieldCalculation getFinished(FourierFilterMaskSource mask, Layer t, int[][] bragg, double angleRad, boolean fancy, boolean forceEven, int code)
			{
				UFieldCalculation calc = new UFieldCalculation(mask, t, null, bragg, angleRad, forceEven);
				if (code == 2)
				{
					calc.doRegularization();
					calc.addRegularization();
					calc.applyUField(fancy);
					return calc;
				}
				else
				{
					calc.doExpuCalc();
					calc.pickBestLayer();
					calc.makeBestUField();
					if (code == 1){
						calc.doRegularization();
						calc.addRegularization();
					}
					calc.applyUField(fancy);
					return calc;
				}
			}
			
			public static FourierFilterMaskSource getDefaultMask(int[][] bragg, int N, int desiredMemUsage)
			{
				int memForExtraLayers = desiredMemUsage - getMinimumMemUsage(N);
				if (memForExtraLayers < 0){
					System.out.println("Error! You have not specified enough memory to do even one layer.");
					System.exit(0);
				}
				int nlayers = memForExtraLayers/getExtraMemUsagePerLayer(N);
				double braggMag = Distance.distance(bragg[0][0], bragg[0][1]);
				return new DriftCorrectionMethods.FourPointedShapeMask(1, braggMag/6, nlayers, N, 1);
			}
			
	 }
	 /**
	  * This class is designed to do the calculation in its real space version (takes longer but should in principle give slightly better results).
	  * It will also determine the best layer automatically. The real-space length scales will be determined automatically in atomic units.
	  * @author madhavanlab2011
	  *
	  */
	 public static class UFieldCalculationReal extends DriftCorrectionCalculator
	 {
		 	boolean stopEarly = false; //This one will stop the calculation as soon as the energy is < 0.
		 	//Since the length scales are in increasing order, the energy will increase so this is safe.
		 	
		 	double meanData;
		 	
			double[][] gaussMask;
			
			double[] lengths;
			
			public UFieldCalculationReal(Layer t, BicubicSplineInterpolatingFunction interp, int[][] bragg, double angle, int nlayers, boolean forceEven)
			{
				super(t, interp, bragg, angle, nlayers, forceEven);
				lengths = new double[nlayers];
				double unitLength = N/(Distance.distance(bragg[0][0], bragg[0][1])); //This is one atomic space.
				double maxLength = nlayers > 1 ? unitLength * Math.sqrt(4*nlayers) : unitLength;
				System.out.println(maxLength);
				
				lengths = ArrayOps.generateArrayNotInclLower(0, maxLength, nlayers);
				energy = new double [nlayers];
				//The preliminary steps consist of filling up the "phaseTopo" array.
				doSummation();
			}
			public UFieldCalculationReal(Layer t, BicubicSplineInterpolatingFunction interp, int[][] bragg, double angle, double[] lengths, boolean forceEven)
			{
				super (t, interp, bragg, angle, lengths.length, forceEven);
				this.lengths = lengths;
				//The preliminary steps consist of filling up the "phaseTopo" array.
				doSummation();
			}
			public void doSummation()
			{
				for (int i = 0; i < 2; i++)
					for (int j = 0; j < N; j++)
						for (int k = 0; k < N; k++){
							phaseTopo[i][j][k][0] = Math.cos(DataManip.dot(braggTrue[i], j, k))*source.data[j][k];
							phaseTopo[i][j][k][1] = -Math.sin(DataManip.dot(braggTrue[i], j, k))*source.data[j][k];
						}
			}				

			public void doExpUCalc()
			{
				for (int j = 0; j < nlayers; j++){
					for (int i = 0; i < 2; i++){
						System.out.println("On Bragg Vector # " + (i+1) + ", doing calculation " + (j+1) + " out of " + nlayers + ".");
						gaussMask = Mask.getGaussianMask(lengths[j]);
						DataManip.getDeviationGaussDefault(lengths[j], phaseTopo[i], expU, gaussMask);
						for (int m = 0; m < N; m++)
							for (int n = 0; n < N; n++)
									expUPhase[i][m][n] = FieldOps.atan(expU[m][n][0], expU[m][n][1]);
						//Once we put it in phaseCont, we'll be OK to reuse expU and expUPhase
						FieldOps.putPhaseSteps(expUPhase[i], N/2, N/2, phaseN[i]);
						FieldOps.putAddedPhaseSteps(expUPhase[i], phaseN[i], 2*Math.PI, phaseCont[j][i]);
						FieldOps.subtractAvg(phaseCont[j][i]);
						
					}
					if (stopEarly)
					{
						energy[j] = getEnergy(phaseCont[j]);
						if (energy[j] < 0)
						{
							selectedLayer = j;
							return;
						}
					}
	 			}
			}
			
			
			public String[] getEnergyLines()
			{
				String[] ans = new String[nlayers+2];
				ans[0] = "Index\tLength Scale\tQuality";
				for (int i = 1; i < nlayers+1; i++)
					ans[i] = "" + (i-1) + "\t" + lengths[i-1] + "\t" + energy[i-1];
				ans[nlayers+1] = "The selected layer was " + selectedLayer + ".";
				return ans;
			}
			
			
			
			
			/**
			 * This embodies the previous methods getFinished_NonlinearOnly, getFinished_LinearOnly, getFinished_Both.
			 * The code will be 0 for nonlinear only, 1 for both, 2 for linear only. I expect only 0 and 1 will be used.
			 * @param mask
			 * @param t
			 * @param bragg
			 * @param fancy
			 * @param forceEven
			 * @param code
			 * @return
			 */
			public static UFieldCalculationReal getFinished(int layers, Layer t, int[][] bragg, double angleRad, boolean fancy, boolean forceEven, int code)
			{
				UFieldCalculationReal calc = new UFieldCalculationReal(t, null, bragg, angleRad, layers, forceEven);
				if (code == 2)
				{
					calc.doRegularization();
					calc.addRegularization();
					calc.applyUField(fancy);
					return calc;
				}
				else
				{
					calc.stopEarly = true;
					calc.doExpUCalc();
					calc.pickBestLayer();
					calc.makeBestUField();
					if (code == 1){
						calc.doRegularization();
						calc.addRegularization();
					}
					calc.applyUField(fancy);
					return calc;
				}
			}
//			/**The default angle is 90 degrees since no regularization will take place.
//			 *
//			 * @param mask
//			 * @param t
//			 * @param bragg
//			 * @return
//			 */
//			public static UFieldCalculationReal getFinished_NonlinearOnly(int nlayers, Layer t, int[][] bragg, boolean fancy, boolean forceEven)
//			{
//				UFieldCalculationReal calc = new UFieldCalculationReal(t, null, bragg, Math.PI/2, nlayers, forceEven);
//				calc.stopEarly = true;
//				calc.doExpUCalc();
//				calc.pickBestLayer();
//				calc.makeBestUField();
//				calc.applyUField(fancy);
//				return calc;
//			}
//			public static UFieldCalculationReal getFinished_LinearOnly(int nlayers, Layer t, int[][] bragg, double angleRad, boolean fancy, boolean forceEven)
//			{
//				UFieldCalculationReal calc = new UFieldCalculationReal(t, null, bragg, angleRad, nlayers, forceEven);
//				calc.doRegularization();
//				calc.addRegularization();
//				calc.applyUField(fancy);
//				return calc;
//			}
//			public static UFieldCalculationReal getFinished_Both(int nlayers, Layer t, int[][] bragg, double angleRad, boolean fancy, boolean forceEven)
//			{
//				UFieldCalculationReal calc = new UFieldCalculationReal(t, null, bragg, angleRad, nlayers, forceEven);
//				calc.stopEarly = true;
//				calc.doExpUCalc();
//				calc.pickBestLayer();
//				calc.makeBestUField();
//				calc.doRegularization();
//				calc.addRegularization();
//				calc.applyUField(fancy);
//				return calc;
//			}
	 }

		
		/**
		 * This class is supposed to contain the functions common to UFieldCalculation and UFieldCalculationReal.
		 * @author Dan
		 *
		 */
		public static abstract class DriftCorrectionCalculator
		{
		 	int N;
		 	int nlayers;
		 	Layer source;
		 	int[][] bragg;
		 	double[][] braggTrue; // 
		 	boolean calcWasGood;
		 	double angle;
		 	BraggVectorRegularizer reg;
			double[][][] uReg;
			int[][] braggReg;
		 	AtomicCoordinatesSet latt = null, lattReg = null;
		 	
		 	boolean[][] outsidePixels;
		 	
		 	boolean stopEarly = false; //This one will stop the calculation as soon as the energy is < 0.
		 	//Since the length scales are in increasing order, the energy will increase so this is safe.
		 	
		 	double meanData;
		 	
			double[][][] expU;
			double[][][] expUPhase;
			int[][][] phaseN;
			double[][][][] phaseCont;
			double[][][][] phaseTopo;
			double[][][] u; //Here U is in the form [n][n][2]
			double[][] after;
			double[][][] fftz;
			BicubicSplineInterpolator erp = new BicubicSplineInterpolator(); // to get the values of the source at non-integer pixels.
			BicubicSplineInterpolatingFunction interp;
			double[][] gaussMask;
			
			double[] lengths;
			
			double[] extraTranslation = new double [2];
			
			double[] energy;
		 	int selectedLayer;
			public DriftCorrectionCalculator(Layer t, BicubicSplineInterpolatingFunction interp, int[][] bragg, double angle, int nlayers, boolean even)
			{
				source = t;
				if (interp == null && t.nx < 2048)
					this.interp = FieldOps.getBicubicInterpoation(source.data);//erp.interpolate(ArrayOps.generateArrayInclBoth(0, t.nx-1, t.nx), ArrayOps.generateArrayInclBoth(0, t.ny-1, t.ny), t.data);
				else if (interp != null)
					this.interp = interp;
				else interp = null;
				meanData = FieldOps.mean(t.data);
				N = t.nx;
				this.angle = angle;
				this.bragg = bragg;
				this.braggTrue = getBraggTrue(bragg, N);
				latt = new AtomicCoordinatesSet(braggTrue[0], braggTrue[1], new double [] {N/2, N/2}).getReciprocalLattice();
				
				this.nlayers = nlayers;
				expU = new double [t.nx][t.ny][2];
				expUPhase = new double [2][t.nx][t.ny];
				phaseN = new int [2][t.nx][t.ny];
				phaseCont = new double [nlayers][2][t.nx][t.ny];
				u = new double [t.nx][t.ny][2];
				after = new double [t.nx][t.ny];
				fftz = new double[t.nx][t.ny][2];
				FFTOps.putFFT(t.data, fftz, false);

				lengths = new double[nlayers];
				double unitLength = N/(Distance.distance(bragg[0][0], bragg[0][1])); //This is one atomic space.
				double maxLength = nlayers > 1 ? unitLength * Math.sqrt(4*nlayers) : unitLength;
				System.out.println(maxLength);
				
				lengths = ArrayOps.generateArrayNotInclLower(0, maxLength, nlayers);
				energy = new double [nlayers];

				reg = new BraggVectorRegularizer(angle, N, even);
				uReg = new double [N][N][2];
				phaseTopo = new double [2][N][N][2];
			}

			public static int getNLayers(int N, int desiredMemUsage)
			{
				int memForExtraLayers = desiredMemUsage - getMinimumMemUsage(N);
				if (memForExtraLayers < 0){
					System.out.println("Error! You have not specified enough memory to do even one layer.");
					System.exit(0);
				}
				return memForExtraLayers/getExtraMemUsagePerLayer(N);
			}
			
			public static int getMinimumMemUsage(int N)
			{
				return 8*N*N*24; //24 is approximately the number of NxN double fields we will have. 
			}
			public static int getExtraMemUsagePerLayer(int N)
			{
				return 8*N*N*3;
			}
			public abstract String[] getEnergyLines();
			public void doRegularization()
			{
				reg.setFirstBraggPeak(braggTrue[0]);
				reg.setSecondBraggPeak(braggTrue[1]);
				reg.regularize();
				lattReg = reg.getRegLattice();
				braggReg = reg.getFinalBragg();
				FieldOps.putUField(latt, lattReg, uReg);
			}

			/**
			 * In order to settle on the best u-field a definition of quality is necessary. We must have as few vortices as possible, but we also want 
			 * the shortest possible length scale (that is, phaseCont should look bubbly). So, the "energy" is defined as: # of vortices - sum(|grad phaseCont|)/(PI*N*N).
			 * This divisor will cause vortices to be prioritized over short length scale.
			 * @param phaseCont
			 * @return
			 */
			public static double getEnergy(double[][][] phaseCont)
			{
				double energy = 0;
				for (int i = 0; i < 2; i++)
				{
					energy += //FieldOps.getNVortices(phaseCont[i], Math.PI);
							FieldOps.getNLargeDifferences(phaseCont[i], Math.PI/4);
					
					energy -= ArrayOps.sum(FieldOps.gradMag(phaseCont[i]))/(Math.PI*phaseCont[0].length*phaseCont[0].length);
				}
				return energy;
			}
			public void pickBestLayer()
			{
				for (int i = 0; i < nlayers; i++)
					energy[i] = getEnergy(phaseCont[i]);
				
				calcWasGood = ArrayOps.min(energy) < 0;
				selectedLayer = ArrayOps.minIndex(energy);
			}

			public void makeBestUField()
			{			
				FieldOps.putU(phaseCont[selectedLayer][0], phaseCont[selectedLayer][1], braggTrue, 1, u);
			}

			/**
			 * This is expected to define the extra translation properly, regardless of whether the drift correction is applied fancily or not.
			 * @param fancy
			 */
			public void applyUField(boolean fancy)
			{
				if (fancy){
					int[][] shiftBragg = lattReg == null ? bragg : braggReg;
					after = DriftCorrectionMethods.applyUFieldSpecial_withShifting(u, source.data, shiftBragg, extraTranslation);
					outsidePixels = FieldOps.getOutsidePixels(u);
					FieldOps.zero(after, outsidePixels);
					FieldOps.changeZeroToAverage(after);
				}
				else{
					if (interp != null) FieldOps.applyUFieldBiCubic(interp, u, after, meanData);
					else	after = FieldOps.applyUField(source.data, u, 4, 4);
					int[][] shiftBragg = lattReg == null ? bragg : braggReg;
					extraTranslation = getTranslationSpecialQuickDirty(u, source.data, shiftBragg);
					//To encourage the user to use fancy drift correction which is the
					//only one that doesn't suck, I will calculate the ExtraTranslation but
					//NOT apply it.
					outsidePixels = FieldOps.getOutsidePixels(u);
				}
			}

			public void addRegularization()
			{
				FieldOps.add(u, uReg, u);
			}


			public void writeFileOutput(String dir, String layername, boolean didRegularization, Layer original)
			{
				String rk = this.getClass().equals(UFieldCalculation.class) ? "k" : "r";
				Layer.writeBIN(Layer.newLayer(source, after), dir + layername + "_after" + (didRegularization ? "Reg" : "") + rk + ".bin");
				
				File outdir = new File(dir + layername + "_Drift corr auto " + rk + (didRegularization ? " reg" : "") + "\\");
				if (!outdir.exists()) outdir.mkdir();
				String outputDir = outdir.toString() + "\\";
				
				Layer.writeBIN(Layer.newLayer(source, phaseCont[selectedLayer][0]), outputDir + "phaseCont_0.bin");
				Layer.writeBIN(Layer.newLayer(source, phaseCont[selectedLayer][1]), outputDir + "phaseCont_1.bin");
				//These two lines write the ACTUAL drift fields. This was added 2/1/2015 because
				//the actual drift fields may be shifted to make the Bragg peaks real, from 
				//what they would be calculated to be from phaseCont and uReg.
//				Layer.writeBIN(Layer.newLayer(source, FieldOps.getIndex(u, 0)), outputDir + "ux.bin");
//				Layer.writeBIN(Layer.newLayer(source, FieldOps.getIndex(u, 1)), outputDir + "uy.bin");

				SRAW.writeImage(outputDir + "outsidePixels.bmp", outsidePixels);
			
				ColumnIO.writeString("" + extraTranslation[0] + "\t" + extraTranslation[1], outputDir + "translation.txt");
				
				ColumnIO.writeString(latt.toString(), outputDir + "lattice.txt");
				ColumnIO.writeLines(getEnergyLines(), outputDir + "Selection info.txt");
				
				ColumnIO.writeTable(FieldOps.transpose(bragg), outputDir + "Bragg.txt", "");
				
				
				if (!calcWasGood) //write topomaps of the phase angles for reference
				{
					File outdirMaps = new File(outputDir + "phase maps\\");
					if (!outdirMaps.exists()) outdirMaps.mkdir();
					String outdirm = outdirMaps.toString() + "\\";
					double[][][][] pc = new double [2][nlayers][][];
					for (int i = 0; i < nlayers; i++)
					{
						pc[0][i] = phaseCont[i][0];
						pc[1][i] = phaseCont[i][1];
					}
					Topomap.writeBIN(new Topomap(pc[0], this.lengths, source.x, source.y, null), outdirm + "phaseCont_0.bin");
					Topomap.writeBIN(new Topomap(pc[1], this.lengths, source.x, source.y, null), outdirm + "phaseCont_1.bin");
				}
				if (didRegularization)
				{
					ColumnIO.writeString(lattReg.toString(), outputDir + "latticeReg.txt");
					Layer.writeBIN(Layer.newLayer(source, FieldOps.getIndex(uReg, 0)), outputDir + "uRegX.bin");
					Layer.writeBIN(Layer.newLayer(source, FieldOps.getIndex(uReg, 1)), outputDir + "uRegY.bin");
					ColumnIO.writeTable(FieldOps.transpose(braggReg), outputDir + "BraggReg.txt", "");
				}
				
				if (original != null)
				{
					Layer afterOrig = applyUField(original);
					Layer.writeBIN(afterOrig, dir + layername + "_orig_after" + (didRegularization ? "Reg" : "") + ".bin");
				}
			}

			/**
			 * This method DOES NOT attempt to write the after file. It ONLY writes the data to allow fancy application of the u-field later.
			 * @param outdir
			 * @param layername
			 * @param didRegularization
			 * @param original
			 * @param i
			 */
			public void writeFileOutputToDir(String outdirs, String layername, boolean didRegularization, Layer original)
			{
//				Layer.writeBIN(Layer.newLayer(source, after), dir + layername + "_after" + (didRegularization ? "Reg" : "") + "r " +  i + ".bin");
				
				File outdir = new File(outdirs);
				if (!outdir.exists()) outdir.mkdir();
				
				Layer.writeBIN(Layer.newLayer(source, phaseCont[selectedLayer][0]), outdirs + "phaseCont_0.bin");
				Layer.writeBIN(Layer.newLayer(source, phaseCont[selectedLayer][1]), outdirs + "phaseCont_1.bin");
				//These two lines write the ACTUAL drift fields. This was added 2/1/2015 because
				//the actual drift fields may be shifted to make the Bragg peaks real, from 
				//what they would be calculated to be from phaseCont and uReg.
//				Layer.writeBIN(Layer.newLayer(source, FieldOps.getIndex(u, 0)), outdirs + "ux.bin");
//				Layer.writeBIN(Layer.newLayer(source, FieldOps.getIndex(u, 1)), outdirs + "uy.bin");
				
				SRAW.writeImage(outdirs + "outsidePixels.bmp", outsidePixels);
				
				ColumnIO.writeString("" + extraTranslation[0] + "\t" + extraTranslation[1], outdirs + "translation.txt");

				ColumnIO.writeString(latt.toString(), outdirs + "lattice.txt");
				ColumnIO.writeLines(getEnergyLines(), outdirs + "Selection info.txt");
				
				ColumnIO.writeTable(FieldOps.transpose(bragg), outdirs + "Bragg.txt", "");

				
				if (didRegularization)
				{
					ColumnIO.writeString(lattReg.toString(), outdirs + "latticeReg.txt");
					Layer.writeBIN(Layer.newLayer(source, FieldOps.getIndex(uReg, 0)), outdirs + "uRegX.bin");
					Layer.writeBIN(Layer.newLayer(source, FieldOps.getIndex(uReg, 1)), outdirs + "uRegY.bin");
					ColumnIO.writeTable(FieldOps.transpose(braggReg), outdirs + "BraggReg.txt", "");
				}
			}

			public void writeFileOutput(String dir, String layername, boolean didRegularization, Layer original, int i, boolean realTrue)
			{
				String direlement = realTrue ? "r" : "k";
				Layer.writeBIN(Layer.newLayer(source, after), dir + layername + "_after" + (didRegularization ? "Reg" : "") + direlement + " " +  i + ".bin");
				
				File outdir = new File(dir + layername + "_Drift corr auto r" + (didRegularization ? " reg" : "") + "\\");
				if (!outdir.exists()) outdir.mkdir();
				String outputDir = outdir.toString() + "\\";
				
				Layer.writeBIN(Layer.newLayer(source, phaseCont[selectedLayer][0]), outputDir + "phaseCont_0_" + i + ".bin");
				Layer.writeBIN(Layer.newLayer(source, phaseCont[selectedLayer][1]), outputDir + "phaseCont_1_" + i + ".bin");
				//These two lines write the ACTUAL drift fields. This was added 2/1/2015 because
				//the actual drift fields may be shifted to make the Bragg peaks real, from 
				//what they would be calculated to be from phaseCont and uReg.
//				Layer.writeBIN(Layer.newLayer(source, FieldOps.getIndex(u, 0)), outputDir + "ux.bin");
//				Layer.writeBIN(Layer.newLayer(source, FieldOps.getIndex(u, 1)), outputDir + "uy.bin");
				
				SRAW.writeImage(outputDir + "outsidePixels.bmp", outsidePixels);
				
				ColumnIO.writeString("" + extraTranslation[0] + "\t" + extraTranslation[1], outputDir + "translation_" + i + ".txt");

				ColumnIO.writeString(latt.toString(), outputDir + "lattice.txt");
				ColumnIO.writeLines(getEnergyLines(), outputDir + "Selection info_" + i + ".txt");
				
				ColumnIO.writeTable(FieldOps.transpose(bragg), outputDir + "Bragg.txt", "");

				
				if (didRegularization && i == 0)
				{
					ColumnIO.writeString(lattReg.toString(), outputDir + "latticeReg.txt");
					Layer.writeBIN(Layer.newLayer(source, FieldOps.getIndex(uReg, 0)), outputDir + "uRegX.bin");
					Layer.writeBIN(Layer.newLayer(source, FieldOps.getIndex(uReg, 1)), outputDir + "uRegY.bin");
					ColumnIO.writeTable(FieldOps.transpose(braggReg), outputDir + "BraggReg.txt", "");
				}
				
				if (original != null)
				{
					Layer afterOrig = applyUField(original);
					Layer.writeBIN(afterOrig, dir + layername + "_orig_after" + (didRegularization ? "Reg" : "") + ".bin");
				}
			}
			public Layer applyUField(Layer t)
			{
				BicubicSplineInterpolatingFunction interp = FieldOps.getBicubicInterpoation(source.data);//erp.interpolate(ArrayOps.generateArrayInclBoth(0, t.nx-1, t.nx), ArrayOps.generateArrayInclBoth(0, t.ny-1, t.ny), t.data);

				double[][] after = FieldOps.applyUField(source.data, u, 4, 4);
				return Layer.newLayer(source, after);
			}


		}
		
		public static double[][] applyUField(boolean fancy, double[][] data, double[][][] u)
		{
			double meanData = FieldOps.mean(data);
			double[][] after = new double [data.length][data[0].length];
			
			if (fancy){
				after = DriftCorrectionMethods.applyUFieldSpecial(u, data);
				boolean[][] outsidePixels = FieldOps.getOutsidePixels(u);
				FieldOps.zero(after, outsidePixels);
				FieldOps.changeZeroToAverage(after);
			}
			else{
				BicubicSplineInterpolatingFunction interp = data.length >= 2048 ?  null : FieldOps.getBicubicInterpoation(data);
				if (interp != null) FieldOps.applyUFieldBiCubic(interp, u, after, meanData);
				else	after = FieldOps.applyUField(data, u, 4, 4);
				boolean[][] outsidePixels = FieldOps.getOutsidePixels(u);
			}
			return after;
		}
		public static Topomap applyTopomapOfUFieldsToLayer(Layer t, Topomap ux, Topomap uy)
		{
			double[][][] data = new double [ux.nlayers][t.nx][t.ny];
			for (int i = 0; i < data.length; i++)
			{
				data[i] = applyUField(false, t.data, FieldOps.getNM2Array(new double[][][] {ux.data[i], uy.data[i]}));
//				Layer.writeBIN(Layer.newLayer(t, data[i]));
			}
			return Topomap.newTopomap(ux, data);
		}
		
		public static void getUFieldTopomap_Real(Layer t, int[][] bragg, double[] L)
		{
			UFieldCalculationReal calc;
			Topomap ux, uy;
			double[][][] uxd = new double [L.length][t.nx][t.ny];
			double[][][] uyd = new double [L.length][t.nx][t.ny];
			
			for (int i = 0; i < L.length; i++)
			{
				calc = new UFieldCalculationReal(t, null, bragg, Math.PI/2, new double[] {L[i]}, false);
				calc.doExpUCalc();
				calc.pickBestLayer();
				calc.makeBestUField();
				uxd[i] = FieldOps.getIndex(calc.u, 0);
				uyd[i] = FieldOps.getIndex(calc.u, 1);
			}
			
			ux = new Topomap(uxd, L, t.x, t.y, null);
			uy = new Topomap(uyd, L, t.x, t.y, null);
			new TopomapViewer_complex2(ux, uy, Topomap.stddir, 512);
		}
		
		/**
		 * This returns seven objects: The four derivative elements in lattice directions:
		 * uaa, uab, uba, ubb in that order, followed by the shear element epsilon12, the rotational element w12, and finally
		 * half the trace of gradU, being the isotropic strain.
		 * 
		 * The calculation is of course done in real space.
		 * @param t
		 * @param braggi
		 * @param L
		 * @return
		 */
		public double[][][] measureTheStrain(Layer t, int[][] braggi, double L){
			UFieldCalculationReal calc;
//			double[][] bragg = new double [2][2]; //This is braggTrue.
			double[][] bragg = getBraggTrue(braggi, t.nx);
			calc = new UFieldCalculationReal(t, null, braggi, Math.PI/2, new double[] {L}, false);
			//Instead of doing calc.doExpUCalc we take each part of the phaseTopo and smooth it ourselves.
			double[][][][] ourPhaseTopo = FieldOps.copy(calc.phaseTopo);
			double[][][] realPart = new double [2][t.nx][t.ny];
			double[][][] imagPart = new double [2][t.nx][t.ny];
			double[][][] realSmoo = new double [2][t.nx][t.ny];
			double[][][] imagSmoo = new double [2][t.nx][t.ny];
			double mag;
			double[][][][] gradReal = new double [2][2][t.nx][t.ny];
			double[][][][] gradImag = new double [2][2][t.nx][t.ny];
			for (int i =0; i < 2; i++)
			{
				realPart[i] = FieldOps.getIndex(ourPhaseTopo[i], 0);
				imagPart[i] = FieldOps.getIndex(ourPhaseTopo[i], 1);
				realSmoo[i] = FieldOps.gaussSmooth(realPart[i], L);
				imagSmoo[i] = FieldOps.gaussSmooth(imagPart[i], L);
				//Now we divide by the magnitude of Smoo in order to make it really exp(-i kdotu) or whatever.
				for (int p = 0; p < t.nx; p++)
					for (int q = 0; q < t.ny; q++)
					{
						mag = Distance.distance(realSmoo[i][p][q], imagSmoo[i][p][q]);
						realSmoo[i][p][q] /= mag;
						imagSmoo[i][p][q] /= mag;
					}
				//Take the gradient
				gradReal[i] = FieldOps.gradient(realSmoo[i]);
				gradImag[i] = FieldOps.gradient(imagSmoo[i]);
				
				//Now we must multiply in the complex sense, the complex conjugate of exp(-i kdotu) by this thing, so that 
				//hopefully what remains is grad(k dot u) times i.
				double[] tempA = new double [2], tempB = new double [2], tempC = new double [2];
				for (int c = 0; c < 2; c++)
					for (int p = 0; p < t.nx; p++)
						for (int q = 0; q < t.ny; q++)
						{
							tempA[0] = gradReal[i][c][p][q];
							tempA[1] = gradImag[i][c][p][q];
							tempB[0] = realSmoo[i][p][q];
							tempB[1] = -imagSmoo[i][p][q];
							Complex.product(tempA, tempB, tempC);
							gradReal[i][c][p][q] = tempC[0];
							gradImag[i][c][p][q] = tempC[1];
						}
				//Now, gradReal should hopefully be zero, and gradImag should by +/- i * grad(k dot u)
			}
			
			//We must now take the imaginary gradient parts and convert them into the elements of the strain tensor
			//This uses the same math as Drift correction. We cannot assume that the Bragg vectors are orthogonal, and so
			//we choose the first Bragg vector and some other perpindicular vector to be a and b
			double[] mags = new double [2]; mags[0] = Complex.mag(bragg[0]); mags[1] = Complex.mag(bragg[1]);
			double[] k1Hat = new double [2], k2Hat = new double [2];
			k1Hat[0] = bragg[0][0]/mags[0]; k1Hat[1] = bragg[0][1]/mags[0];
			k2Hat[0] = bragg[1][0]/mags[1]; k2Hat[1] = bragg[1][1]/mags[1];

			double costh = k1Hat[0]*k2Hat[0] + k1Hat[1]*k2Hat[1];

			//use the gramm-schmidt procedure to form k2PrimeHat from k1Hat and k2Hat.
			double[] k2PrimeHat = new double [2];
			k2PrimeHat[0] = k2Hat[0] - costh*k1Hat[0];
			k2PrimeHat[1] = k2Hat[1] - costh*k1Hat[1];
			double magPrime = Complex.mag(k2PrimeHat);
			k2PrimeHat[0] /= magPrime; k2PrimeHat[1] /= magPrime;
			
			//The dot products of u with the unit vectors
			double[][] uaa = new double [t.nx][t.ny];
			double[][] uab = new double [t.nx][t.ny];
			double[][] uba = new double [t.nx][t.ny];
			double[][] ubb = new double [t.nx][t.ny];
			double[][][] gradImag1Prime = new double[2][t.nx][t.ny];
			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					uaa[i][j] = (gradImag[0][0][i][j]*k1Hat[0] + gradImag[0][1][i][j]*k1Hat[1])/mags[0];
					uab[i][j] = (gradImag[0][0][i][j]*k2PrimeHat[0] + gradImag[0][1][i][j]*k2PrimeHat[1])/mags[0];
					for (int k = 0; k < 2; k++){
						gradImag1Prime[k][i][j] = (gradImag[1][k][i][j] - costh*gradImag[0][k][i][j])/(magPrime);
					}
					uba[i][j] = (gradImag1Prime[0][i][j]*k1Hat[0] + gradImag1Prime[1][i][j]*k1Hat[1])/mags[1];
					ubb[i][j] = (gradImag1Prime[0][i][j]*k2PrimeHat[0] + gradImag1Prime[1][i][j]*k2PrimeHat[1])/mags[1];
				}
			
			double[][] e12 = new double [t.nx][t.ny];
			double[][] w12 = new double [t.nx][t.ny];
			double[][] tr_e = new double [t.nx][t.ny];
			double[][] anis = new double [t.nx][t.ny];
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
				{
					e12[i][j] = (uab[i][j] + uba[i][j])/2;
					w12[i][j] = (uab[i][j] - uba[i][j])/2;
					tr_e[i][j] = (uaa[i][j] + ubb[i][j])/2;
					anis[i][j] = (uaa[i][j] - ubb[i][j])/2;
				}
			
			//Now that I have the elements of the derivative tensor, I can make the other tensorial elements:
			//The cauchy strain tensor has its shear element e12 
		//Now I hopfully have the correct stuff of the strain tensor defined with respect to the lattice.
			//Let us check:
//			Layer u1x = Layer.newLayer(t, gradImag[0][0]);
//			Layer u1y = Layer.newLayer(t, gradImag[0][1]);
//			Layer u1xr = Layer.newLayer(t, gradReal[0][0]);
//			Layer u1yr = Layer.newLayer(t, gradReal[0][1]);
//			
//			String f = FileOps.selectSave(fc).toString();
//			Layer.writeBIN(u1x, f + "u1x.bin");
//			Layer.writeBIN(u1y, f + "u1y.bin");
//			Layer.writeBIN(u1xr, f + "u1xr.bin");
//			Layer.writeBIN(u1yr, f + "u1yr.bin");
			return new double[][][] {uaa, uab, uba, ubb, e12, w12, tr_e, anis};
		}

		public static double[][][] measureTheStrainFourier(Layer t, int[][] braggi, double[][] ffiltWeights){
			UFieldCalculationReal calc;
//			double[][] bragg = new double [2][2]; //This is braggTrue.
			double[][] bragg = getBraggTrue(braggi, t.nx);
			calc = new UFieldCalculationReal(t, null, braggi, Math.PI/2, new double[] {1}, false);
			//Instead of doing calc.doExpUCalc we take each part of the phaseTopo and smooth it ourselves.
			double[][][][] ourPhaseTopo = FieldOps.copy(calc.phaseTopo);
			
			double[][][] realPart = new double [2][t.nx][t.ny];
			double[][][] imagPart = new double [2][t.nx][t.ny];
			double[][][] realSmoo = new double [2][t.nx][t.ny];
			double[][][] imagSmoo = new double [2][t.nx][t.ny];
			double mag;
			double[][][][] gradReal = new double [2][2][t.nx][t.ny];
			double[][][][] gradImag = new double [2][2][t.nx][t.ny];
			
			for (int i = 0; i < 2; i++)
			{
				realPart[i] = FieldOps.getIndex(ourPhaseTopo[i], 0);
				imagPart[i] = FieldOps.getIndex(ourPhaseTopo[i], 1);
				//Now we apply the fourier-filter mask to the real and imaginary parts:
				FFT2DSmall fft;
				fft = FFTOps.obtainFFT(realPart[i]);
				FFTOps.supressModes(fft, ffiltWeights);
				fft.doIFFT();
				for (int j = 0; j < t.nx; j++)
					for (int k = 0; k < t.ny; k++)
						realSmoo[i][j][k] = fft.f[j][k][0];
				fft = FFTOps.obtainFFT(imagPart[i]);
				FFTOps.supressModes(fft, ffiltWeights);
				fft.doIFFT();
				for (int j = 0; j < t.nx; j++)
					for (int k = 0; k < t.ny; k++)
						imagSmoo[i][j][k] = fft.f[j][k][0];
				
				//Now we divide by the magnitude of Smoo in order to make it really exp(-i kdotu) or whatever.
				for (int p = 0; p < t.nx; p++)
					for (int q = 0; q < t.ny; q++)
					{
						mag = Distance.distance(realSmoo[i][p][q], imagSmoo[i][p][q]);
						realSmoo[i][p][q] /= mag;
						imagSmoo[i][p][q] /= mag;
					}
				//Take the gradient
				gradReal[i] = FieldOps.gradient(realSmoo[i]);
				gradImag[i] = FieldOps.gradient(imagSmoo[i]);
				
				//Now we must multiply in the complex sense, the complex conjugate of exp(-i kdotu) by this thing, so that 
				//hopefully what remains is grad(k dot u) times i.
				double[] tempA = new double [2], tempB = new double [2], tempC = new double [2];
				for (int c = 0; c < 2; c++)
					for (int p = 0; p < t.nx; p++)
						for (int q = 0; q < t.ny; q++)
						{
							tempA[0] = gradReal[i][c][p][q];
							tempA[1] = gradImag[i][c][p][q];
							tempB[0] = realSmoo[i][p][q];
							tempB[1] = -imagSmoo[i][p][q];
							Complex.product(tempA, tempB, tempC);
							gradReal[i][c][p][q] = tempC[0];
							gradImag[i][c][p][q] = tempC[1];
						}
				//Now, gradReal should hopefully be zero, and gradImag should by +/- i * grad(k dot u)
			}
			
			//We must now take the imaginary gradient parts and convert them into the elements of the strain tensor
			//This uses the same math as Drift correction. We cannot assume that the Bragg vectors are orthogonal, and so
			//we choose the first Bragg vector and some other perpindicular vector to be a and b
			double[] mags = new double [2]; mags[0] = Complex.mag(bragg[0]); mags[1] = Complex.mag(bragg[1]);
			double[] k1Hat = new double [2], k2Hat = new double [2];
			k1Hat[0] = bragg[0][0]/mags[0]; k1Hat[1] = bragg[0][1]/mags[0];
			k2Hat[0] = bragg[1][0]/mags[1]; k2Hat[1] = bragg[1][1]/mags[1];

			double costh = k1Hat[0]*k2Hat[0] + k1Hat[1]*k2Hat[1];

			//use the gramm-schmidt procedure to form k2PrimeHat from k1Hat and k2Hat.
			double[] k2PrimeHat = new double [2];
			k2PrimeHat[0] = k2Hat[0] - costh*k1Hat[0];
			k2PrimeHat[1] = k2Hat[1] - costh*k1Hat[1];
			double magPrime = Complex.mag(k2PrimeHat);
			k2PrimeHat[0] /= magPrime; k2PrimeHat[1] /= magPrime;
			
			//The dot products of u with the unit vectors
			double[][] uaa = new double [t.nx][t.ny];
			double[][] uab = new double [t.nx][t.ny];
			double[][] uba = new double [t.nx][t.ny];
			double[][] ubb = new double [t.nx][t.ny];
			double[][][] gradImag1Prime = new double[2][t.nx][t.ny];
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
				{
					uaa[i][j] = (gradImag[0][0][i][j]*k1Hat[0] + gradImag[0][1][i][j]*k1Hat[1])/mags[0];
					uab[i][j] = (gradImag[0][0][i][j]*k2PrimeHat[0] + gradImag[0][1][i][j]*k2PrimeHat[1])/mags[0];
					for (int k = 0; k < 2; k++){
						gradImag1Prime[k][i][j] = (gradImag[1][k][i][j] - costh*gradImag[0][k][i][j])/(magPrime);
					}
					uba[i][j] = (gradImag1Prime[0][i][j]*k1Hat[0] + gradImag1Prime[1][i][j]*k1Hat[1])/mags[1];
					ubb[i][j] = (gradImag1Prime[0][i][j]*k2PrimeHat[0] + gradImag1Prime[1][i][j]*k2PrimeHat[1])/mags[1];
				}
			
			double[][] e12 = new double [t.nx][t.ny];
			double[][] w12 = new double [t.nx][t.ny];
			double[][] tr_e = new double [t.nx][t.ny];
			double[][] anis = new double [t.nx][t.ny];
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
				{
					e12[i][j] = (uab[i][j] + uba[i][j])/2;
					w12[i][j] = (uab[i][j] - uba[i][j])/2;
					tr_e[i][j] = (uaa[i][j] + ubb[i][j])/2;
					anis[i][j] = (uaa[i][j] - ubb[i][j])/2;
				}
			
			//Now that I have the elements of the derivative tensor, I can make the other tensorial elements:
			//The cauchy strain tensor has its shear element e12 
		//Now I hopfully have the correct stuff of the strain tensor defined with respect to the lattice.
			//Let us check:
//			Layer u1x = Layer.newLayer(t, gradImag[0][0]);
//			Layer u1y = Layer.newLayer(t, gradImag[0][1]);
//			Layer u1xr = Layer.newLayer(t, gradReal[0][0]);
//			Layer u1yr = Layer.newLayer(t, gradReal[0][1]);
//			
//			String f = FileOps.selectSave(fc).toString();
//			Layer.writeBIN(u1x, f + "u1x.bin");
//			Layer.writeBIN(u1y, f + "u1y.bin");
//			Layer.writeBIN(u1xr, f + "u1xr.bin");
//			Layer.writeBIN(u1yr, f + "u1yr.bin");
			return new double[][][] {uaa, uab, uba, ubb, e12, w12, tr_e, anis};
		}

		public static void writeNewTranslationsToAFolder(Topomap t, String dir)
		{
			double[][][][] u = getUTopomapFromFolder(dir, null);
			for (int i = 0; i < t.nlayers; i++)
			{
				double[] dr = getTranslationSpecialQuickDirty(u[i], t.data[i], braggReserve);
				File translation = new File(dir + "translation_" + i + ".txt");	
				if (!translation.exists()) 	ColumnIO.writeString("" + dr[0] + "\t" + dr[1], dir + "translation_" + i + ".txt");

			}
		}

}
