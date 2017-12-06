package main;

import image.ImageEditing;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

import javax.swing.JFileChooser;

import schrodinger.MovieMaker;
import util.ArrayOps;
import util.Complex;
import util.FieldOps;
import util.Printer;
import util.color.ColorScales;
import util.fileops.ColumnIO;
import util.fileops.FileOps;
import util.fileops.Layer;
import util.fileops.Topomap;
import util.fourier.FFT2D;
import util.fourier.FFT2DSmall;
import util.fourier.FFT2D_Subset;
import util.fourier.FFTOps;
import util.fourier.ImpurityEmergenceEditor;
import util.fourier.SFT;
import util.geom.AtomicCoordinatesSet;
import util.geom.Mask;
import util.matrix.Matrix;
import util.robot.Robo;
import util.scalar.functions.TwoDArray;

public class DataManip {

	public static void main(String[] args)
	{
		String dir = "C:\\data\\lawlerhex\\";
		Topomap.setStdDir();
		JFileChooser fc = new JFileChooser(Topomap.stddir);
//		makeCentroidField();
//		subtractAverage();
//		doFFT();
//		doFFTIFFT();
//		smoothOut();
//		makeSubFFTMovie();
//		printSubFFTPeaks();
//		readRLatticeVectors();
//		makeSubFFTData();
		int npts = 1024;
		int nll = 15;
		double[] e = ArrayOps.generateArrayInclBoth(-1, 4, npts);
		double[] g = new double [npts];
		double de = 1.63;
		double[] maxima =  new double [nll];
		for (int i = 0; i < maxima.length; i++)
			maxima[i] = de*Math.sqrt(i);
		double width = 0.18;
		for (int i = 0; i < npts; i++)
			for (int j = 0; j < maxima.length; j++)
				g[i] += width/(Math.pow(e[i]-maxima[j],2) + width*width);
		
		String[] lines = new String[npts];
		for (int i = 0; i < npts; i++)
			lines[i] = "" + g[i] + "\t" + e[i];
		Printer.copyToClipboard(lines);
		System.out.println(Printer.arrayVertical(maxima));
//		suppressBPandWrite(new int[] {258, 269}, 1, "C:\\data\\lawlerhex\\kind.dat", "C:\\data\\lawlerhex\\kindmod");
//		int N = 512;
//		int[] bragg1 = {266, 248};
//		int[] bragg2 = {244, 251};
//		int[] r0 = {0, 0}, dr = {1, 1}, npts = {512, 512};
//		int L = 128;
//		int biggerL = 192;
//		int braggsmear = 1;
//		doDeviationJobGauss("C:\\data\\lawlerhex\\", "kind", new int[][] {bragg1, bragg2}, braggsmear, N, r0, dr, npts, L, biggerL);
//		dir = "C:\\data\\lintrans\\run136topo6_620\\cutoff_pt4a\\";
		dir = "C:\\data\\lawlerhex\\flat8022010\\subloc\\";
//		dir = "C:\\data\\lawlerhex\\run153topo15_830\\subloc\\";
//		dir = "C:\\data\\lawlerhex\\fakeuhex\\";
//		int[] L = new int[7];
//		L[0] = 1;
//		for (int i = 1; i < L.length; i++)
//			L[i] = L[i-1]*2;
//		doUFieldSummation(dir, "didv1", 1, L, 0, 0);
//		double[] L = new double[40];
//		double dl = 0.25;
//		L[0] = dl;
//		for (int i = 1; i < L.length; i++)
//			L[i] = L[i-1] + dl;

//		doUFieldSummationFT(dir, "topo1trans", new double[] {10.5}, 0.95, true);
//		doUFieldSummationFT(dir, "corr", L, 0.8, true);
//		doImpCorrAngleTest();
//		writeImpImagesvsEnergy();
//		getDifferences();
//		writeImpImagesvsEnergyRotated();
//		Layer t = Layer.openFree(fc);
//		AtomicCoordinatesSet latt = new AtomicCoordinatesSet(FileOps.openText(fc));
//		ArrayList<drawing.GaussSquareImpurityAngleSuite.Impurity> imps = drawing.GaussSquareImpurityAngleSuite.getImpurities(FileOps.selectOpen(fc));
//		
//		
//		
//		
//		writeImpImagesOneLayerWithRoot2Splitting(t, imps, latt, fc.getCurrentDirectory().toString() + "\\", false);
//		doUFieldSummationFTWriteBin(dir, "topo", new double[] {6.25, 0, 3.5}, 0.8, true);
//		FieldUtil.addPhaseWrite(dir, "topo", new double[] {6.25, 0, 3.5});
//		DataManip.finishUField(dir, "topo", new double[] {6.25, 0, 3.5}, "corr", true);
//		braggSmearTest("C:\\data\\lawlerhex\\stripe3\\", "topo");
		//		dir = "C:\\data\\lawlerhex\\todan\\";
//		testExpand(ColumnIO.readSquareTable(dir + "topoAxL64cont.dat"), dir);
//		int blocksize = 64;
//		testTranslate();
//		testAutocorrelate(dir);
//		applyTranslationEach(dir, "topo1trans.dat", "topo1transxL10.5cont.dat", "topo1transyL10.5cont.dat", "topo1t");
//		makeRescaleXYMov(dir, "5trans", 270.0811, 1000, 4.3717*Math.sqrt(3)/2);

		//		makeTranslationFFTMov("C:\\data\\lawlerhex\\fakeu\\", "topo.dat", "topoxL4cont.dat", "topoyL4cont.dat", "topocorrect", 2, 200);
//		dir += "bise1\\"; String name = "topo";
//		testTranslation(ColumnIO.readSquareTable(dir + name + ".dat"), dir, name);
//		makeCutoffMovieFFT(dir, "topocorrect", 64, 8, 0, 1, 200, "C:\\Users\\madhavanlab2011\\Documents\\BMP to AVI\\");
//		tryLinearTransBragg(dir, "topo", "bragg3.txt");
//		double[][] delta = new double [512][512];
//		delta[0][0] = 1;
//		delta[0][1] = 0.0001;
//		ColumnIO.writeBin(delta, "C:\\data\\lawlerhex\\delta\\topo.dat");
	}
	
	public static void devJ()
	{
		int N = 512;
		int[] bragg1 = {266, 248};
		int[] bragg2 = {244, 251};
		int[] r0 = {0, 0}, dr = {1, 1}, npts = {416, 416};
		int L = 96;
//		int biggerL = 192;
		int braggsmear = 1;
//		doDeviationJob("C:\\data\\lawlerhex\\", "kind", new int[][] {bragg1, bragEasyg2}, braggsmear, N, r0, dr, npts, L);
		
	}
	
	public static void makeDevJMovie()
	{
		
	}
	
	public static double atan(double x, double y)
	{
		if (x == 0) return y > 0 ? Math.PI/2 : -Math.PI/2;
		if (x < 0) return Math.atan(y/x) + Math.PI;
		if (y < 0) return Math.atan(y/x) + 2*Math.PI;
		return Math.atan(y/x);
	}
	
	public static double[][] magnitude(double[][] x, double[][] y)
	{
		int N = x.length, M = x[0].length;
		double[][] mag = new double[N][M];
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				mag[i][j] = Math.sqrt(x[i][j]*x[i][j] + y[i][j]*y[i][j]);
		return mag;
	}
	public static double[][] phase(double[][] x, double[][] y)
	{
		int N = x.length, M = x[0].length;
		double[][] phase = new double[N][M];
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
				phase[i][j] = atan(x[i][j], y[i][j]);
		return phase;
	}

	static void subtractAvg(double[][] data)
	{
		double avg = 0;
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++)
			{
				avg += data[i][j]; 
			}
		avg /= data.length*data[0].length;
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++)
			{
				data[i][j] -= avg; 
			}
		
	}
	static double[][] load(String filepath)
	{
		return SRAW.getData(filepath);
	}
	
	static void save(double[][] data, String filepath)
	{
		ColumnIO.writeTable(data, filepath);
	}
	static void copy(double[][] target, double[][][] source, int index)
	{
		int N = target.length, M = target[0].length;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
			{
				target[i][j] = source[i][j][index];
			}
	}
	
	
	public static void doSFT()
	{
		String in = "c:\\MinGW\\data\\lattice.txt";
		long t = System.currentTimeMillis();
		double[][] data = load(in);
		System.out.println("" + (System.currentTimeMillis() - t) + " to load.");
		t = System.currentTimeMillis();
		
		int N = data.length;
		
		SFT ft = new SFT(N/4);
		System.out.println("" + (System.currentTimeMillis() - t) + " to make ft.");
		t = System.currentTimeMillis();
		System.out.println(N);
		
		double[][] field = FieldOps.reduce(4, data);
		
		double[][][] datak = new double[2][N/4][N/4];
		ft.sft(field, datak);
		System.out.println();
		System.out.println("" + (System.currentTimeMillis() - t) + " to do ft.");
		t = System.currentTimeMillis();
		
		save(datak[0], "C:\\MinGW\\data\\re.txt");
		save(datak[1], "C:\\MinGW\\data\\im.txt");
	}
	public static void doFFT(){
		String in = "c:\\MinGW\\data\\1026008T.txt";
		long t = System.currentTimeMillis();
		double[][] data = load(in);
		System.out.println("" + (System.currentTimeMillis() - t) + " to load.");
		t = System.currentTimeMillis();
		
		int N = data.length;
		
		System.out.println(N);
		
		double[][] copy = new double [N][N];
		
		FFT2D fft = new FFT2D(data);
		System.out.println("" + (System.currentTimeMillis() - t) + " to create ft.");
		t = System.currentTimeMillis();
		
		fft.doFFT();

		System.out.println("" + (System.currentTimeMillis() - t) + " to do ft.");
		t = System.currentTimeMillis();
		
		double[][] field;
		
		copy(copy, fft.fHat2, 0);
		for (int i = 0; i < N; i++)
			copy[i][N/2] = 0;
		for (int i = 0; i < N; i++)
			copy[N/2][i] = 0;
		field = FieldOps.reduce(4, copy);
		save(field, "C:\\MinGW\\data\\1026008FTRE.txt");
		
		copy(copy, fft.fHat2, 1);
		for (int i = 0; i < N; i++)
			copy[i][N/2] = 0;
		for (int i = 0; i < N; i++)
			copy[N/2][i] = 0;
		field = FieldOps.reduce(4, copy);
		save(field, "C:\\MinGW\\data\\1026008FTIM.txt");
		
	}

	
	public static void doFFTIFFT(){
		String in = "c:\\MinGW\\data\\lattice.txt";
		long t = System.currentTimeMillis();
		double[][] data = load(in);
		System.out.println("" + (System.currentTimeMillis() - t) + " to load.");
		t = System.currentTimeMillis();
		double[][] field = FieldOps.reduce(4, data);
		int N = data.length;
		
		System.out.println(N);
		
		double[][] copyre = new double [N/4][N/4];
		double[][] copyim = new double [N/4][N/4];
		
		FFT2D fft = new FFT2D(field);
		System.out.println("" + (System.currentTimeMillis() - t) + " to create ft.");
		t = System.currentTimeMillis();
		
		fft.doFFT();

		System.out.println("" + (System.currentTimeMillis() - t) + " to do ft.");
		t = System.currentTimeMillis();
		
		copy(copyre, fft.fHat2, 0);
		save(copyre, "C:\\MinGW\\data\\1026008FTRE.txt");
		
		copy(copyim, fft.fHat2, 1);
		save(copyim, "C:\\MinGW\\data\\1026008FTIM.txt");
		System.out.println("" + (System.currentTimeMillis() - t) + " to save1.");
		t = System.currentTimeMillis();
		
		fft.doIFFT();
		System.out.println("" + (System.currentTimeMillis() - t) + " to ifft.");
		t = System.currentTimeMillis();
		
		copy(copyre, fft.f, 0);
//		field = reduce(4, copyre);
		save(copyre, "C:\\MinGW\\data\\1026008IFTRE.txt");
		copy(copyim, fft.f, 1);
//		field = reduce(4, copyim);
		save(copyim, "C:\\MinGW\\data\\1026008IFTIM.txt");
		
		System.out.println("" + (System.currentTimeMillis() - t) + " to save2.");
		t = System.currentTimeMillis();
		
		
	}

	//there are 3 bragg vectors, listed in order (clockwise or counterclockwise).
	public static void tryLinearTrans(int[][] bragg, int N)
	{
		double[][] braggi = new double [3][2];
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 2; j++)
			{
				braggi[i][j] = bragg[i][j] - N/2;
			}
		double th12, th23, th13;
		double th12o, th23o, th13o;
		
		double deform;
		double defmin = 0.8, ddef = 0.4;
		int npts = 1000;
		
		double[][] braggmod = new double [3][2];
		double[] mag = new double [3];
		for (int i = 0; i < 3; i++)
			mag[i] = Complex.mag(braggi[i]);
		
		for (int i = 0; i < npts; i++)
		{	
			deform = 0.8 + i*ddef/npts;
			for (int j = 0; j < 3; j++)
			{
					braggmod[j][0] = braggi[j][0]*deform; braggmod[j][1] = braggi[j][1];
					mag[j] = Complex.mag(braggmod[j]);
			}
			
			th12 = Math.acos(dot(braggmod[0], braggmod[1])/(mag[0]*mag[1]));
			th23 = Math.acos(dot(braggmod[1], braggmod[2])/(mag[1]*mag[2]));
			th13 = Math.acos(dot(braggmod[0], braggmod[2])/(mag[0]*mag[2]));
			th12o = (th12-Math.PI/3)*(th12-Math.PI/3);
			th23o = (th23-Math.PI/3)*(th23-Math.PI/3);
			th13o = (th13-2*Math.PI/3)*(th13-2*Math.PI/3);
			
			System.out.println(deform + "\t" + th12 + "\t" + th23 + "\t" + th13 + "\t" + th12o + "\t" + th23o + "\t" + th13o);
		}
	}
	public static void tryLinearTrans2(int[][] bragg, int N, String dir)
	{
		String subdir = "transphase\\";
		double[][] braggi = new double [3][2];
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < 2; j++)
			{
				braggi[i][j] = bragg[i][j] - N/2;
			}
		double th12, th23, th13;
		double th12o, th23o, th13o;
		
		double deform;
		int npts = 2048;
		double defmin = 0.8, defmax = 1.2;
		double ddef = (defmax - defmin)/(npts-1);
		double phi;
		double phimin = 0, phimax = Math.PI/4;
		double dphi = (phimax - phimin)/(npts-1);

		
		double[][] braggmod = new double [3][2];
		double[] mag = new double [3];
		for (int i = 0; i < 3; i++)
			mag[i] = Complex.mag(braggi[i]);
		
		double alpha, beta, gamma;
		double[][] th12op = new double [npts][npts];
		double[][] th23op = new double [npts][npts];
		double[][] th13op = new double [npts][npts];
		double[][] p1223 = new double [npts][npts];
		double[][] dmag12 = new double [npts][npts];
		double[][] dmag23 = new double [npts][npts];
		double[][] dmag1223 = new double [npts][npts];
		double[][] tout = new double [npts][npts];
		for (int i = 0; i < npts; i++)
			for (int k = 0; k < npts; k++)
			{	
				phi = phimin + i*dphi;
				deform = defmin + k*ddef;
				
				//The matric is (alpha gamma
//								gamma beta) where
				alpha = Math.cos(phi)*Math.cos(phi) + deform*Math.sin(phi)*Math.sin(phi);
				beta = deform*Math.cos(phi)*Math.cos(phi) + Math.sin(phi)*Math.sin(phi);
				gamma = (1 - deform)*Math.sin(phi)*Math.cos(phi);
				
				for (int j = 0; j < 3; j++)
				{
						braggmod[j][0] = braggi[j][0]*alpha + braggi[j][1]*gamma;
						braggmod[j][1] = braggi[j][0]*gamma + braggi[j][1]*beta;
						mag[j] = Complex.mag(braggmod[j]);
				}
				
				th12 = Math.acos(dot(braggmod[0], braggmod[1])/(mag[0]*mag[1]));
				th23 = Math.acos(dot(braggmod[1], braggmod[2])/(mag[1]*mag[2]));
				th13 = Math.acos(dot(braggmod[0], braggmod[2])/(mag[0]*mag[2]));
				th12o = (th12-Math.PI/3)*(th12-Math.PI/3);
				th23o = (th23-Math.PI/3)*(th23-Math.PI/3);
				th13o = (th13-2*Math.PI/3)*(th13-2*Math.PI/3);
				th12op[i][k] = 1/th12o;
				th13op[i][k] = 1/th13o;
				th23op[i][k] = 1/th23o;
				p1223[i][k] = th12op[i][k]*th23op[i][k];
				dmag12[i][k] = 1/Math.pow(mag[1] - mag[0], 2);
				dmag23[i][k] = 1/Math.pow(mag[2] - mag[1], 2);
				dmag1223[i][k] = dmag12[i][k]*dmag23[i][k];
				tout[i][k] = dmag1223[i][k]*p1223[i][k];
//				System.out.println(deform + "\t" + th12 + "\t" + th23 + "\t" + th13 + "\t" + th12o + "\t" + th23o + "\t" + th13o);
		}
		FieldOps.log(th12op);
		FieldOps.log(th13op);
		FieldOps.log(th23op);
		FieldOps.log(p1223);
		FieldOps.log(dmag12);
		FieldOps.log(dmag23);
		FieldOps.log(dmag1223);
		FieldOps.log(tout);
		if (!new File(dir + subdir).exists())
			new File(dir + subdir).mkdir();

		SRAW.writeImage(dir + subdir + "th12", th12op);
		SRAW.writeImage(dir + subdir + "th13", th13op);
		SRAW.writeImage(dir + subdir + "th23", th23op);
		SRAW.writeImage(dir + subdir + "p1223", p1223);
		SRAW.writeImage(dir + subdir + "mag12", dmag12);
		SRAW.writeImage(dir + subdir + "mag23", dmag23);
		SRAW.writeImage(dir + subdir + "mag1223", dmag1223);
		SRAW.writeImage(dir + subdir + "all", tout);
	}
	public static void tryLinearTransBragg(String dir, String name, String braggf)
	{
		int[][] braggi = readBragg3(dir + braggf);
		int braggsmear = readBraggSmear(dir + "bragg.txt");
		
		double[][] data = ColumnIO.readSquareTable(dir + name + ".dat" );
		int N = data.length;
		long t = System.currentTimeMillis();
		System.out.println((System.currentTimeMillis() - t)/1000 + " to load.");
		t = System.currentTimeMillis();
		
		double[][] fftmag = FFTOps.obtainFFTmagCent(data);
		t = printTime(" to fft", t);

		double[][] bragg = new double[3][2];
		for (int i = 0; i < 3; i++){
			bragg[i] = CentroidField.centroid(fftmag, braggi[i][0]-braggsmear, 2*braggsmear + 1, braggi[i][1]-braggsmear, 2*braggsmear + 1, false);
			System.out.println("(" + braggi[i][0] + ", " + braggi[i][1] + ")  ----> (" + bragg[i][0] + ", " + bragg[i][1] + ")" );
			bragg[i][0] -= N/2;
			bragg[i][1] -= N/2;
		}

		String subdir = "transphaseEx\\";
		double th12, th23, th13;
		double th12o, th23o, th13o;
		
		double deform;
		int npts = 2048;
		double defmin = 0.8, defmax = 1.2;
		double ddef = (defmax - defmin)/(npts-1);
		double phi;
		double phimin = 0, phimax = Math.PI/2;
		double dphi = (phimax - phimin)/(npts-1);

		
		double[][] braggmod = new double [3][2];
		double[] mag = new double [3];
		for (int i = 0; i < 3; i++)
			mag[i] = Complex.mag(bragg[i]);
		
		double alpha, beta, gamma;
		double[][] th12op = new double [npts][npts];
		double[][] th23op = new double [npts][npts];
		double[][] th13op = new double [npts][npts];
		double[][] p1223 = new double [npts][npts];
		double[][] dmag12 = new double [npts][npts];
		double[][] dmag23 = new double [npts][npts];
		double[][] dmag1223 = new double [npts][npts];
		double[][] tout = new double [npts][npts];
		for (int i = 0; i < npts; i++)
			for (int k = 0; k < npts; k++)
			{	
				phi = phimin + i*dphi;
				deform = defmin + k*ddef;
				
				//The matric is (alpha gamma
//								gamma beta) where
				alpha = Math.cos(phi)*Math.cos(phi) + deform*Math.sin(phi)*Math.sin(phi);
				beta = deform*Math.cos(phi)*Math.cos(phi) + Math.sin(phi)*Math.sin(phi);
				gamma = (1 - deform)*Math.sin(phi)*Math.cos(phi);
				
				for (int j = 0; j < 3; j++)
				{
						braggmod[j][0] = bragg[j][0]*alpha + bragg[j][1]*gamma;
						braggmod[j][1] = bragg[j][0]*gamma + bragg[j][1]*beta;
						mag[j] = Complex.mag(braggmod[j]);
				}
				
				th12 = Math.acos(dot(braggmod[0], braggmod[1])/(mag[0]*mag[1]));
				th23 = Math.acos(dot(braggmod[1], braggmod[2])/(mag[1]*mag[2]));
				th13 = Math.acos(dot(braggmod[0], braggmod[2])/(mag[0]*mag[2]));
				th12o = (th12-Math.PI/3)*(th12-Math.PI/3);
				th23o = (th23-Math.PI/3)*(th23-Math.PI/3);
				th13o = (th13-2*Math.PI/3)*(th13-2*Math.PI/3);
				th12op[i][k] = 1/th12o;
				th13op[i][k] = 1/th13o;
				th23op[i][k] = 1/th23o;
				p1223[i][k] = th12op[i][k]*th23op[i][k];
				dmag12[i][k] = 1/Math.pow(mag[1] - mag[0], 2);
				dmag23[i][k] = 1/Math.pow(mag[2] - mag[1], 2);
				dmag1223[i][k] = dmag12[i][k]*dmag23[i][k];
				tout[i][k] = dmag1223[i][k]*p1223[i][k];
//				System.out.println(deform + "\t" + th12 + "\t" + th23 + "\t" + th13 + "\t" + th12o + "\t" + th23o + "\t" + th13o);
		}
		FieldOps.log(th12op);
		FieldOps.log(th13op);
		FieldOps.log(th23op);
		FieldOps.log(p1223);
		FieldOps.log(dmag12);
		FieldOps.log(dmag23);
		FieldOps.log(dmag1223);
		FieldOps.log(tout);
		if (!new File(dir + subdir).exists())
			new File(dir + subdir).mkdir();

		SRAW.writeImage(dir + subdir + "th12", th12op);
		SRAW.writeImage(dir + subdir + "th13", th13op);
		SRAW.writeImage(dir + subdir + "th23", th23op);
		SRAW.writeImage(dir + subdir + "p1223", p1223);
		SRAW.writeImage(dir + subdir + "mag12", dmag12);
		SRAW.writeImage(dir + subdir + "mag23", dmag23);
		SRAW.writeImage(dir + subdir + "mag1223", dmag1223);
		SRAW.writeImage(dir + subdir + "all", tout);
	}
	public static void makeSubFFTMovie()
	{
		long t = System.currentTimeMillis();
		int L = 256;
		double[][] data = ColumnIO.readSquareTable("C:\\data\\lawlerhex\\8302010_0010\\topo.dat");
		int N = data.length;
		System.out.println("" + (System.currentTimeMillis() - t) + "to load");
		FFT2D_Subset fft = new FFT2D_Subset(data, L);
		double[][] out = new double[L][L];
		
		for (int i = 0; i < N-L; i++)
		{
			makeFFTPicture(data, (N-L)/2, i, fft, out, i);
		}
		
		Robo r = new Robo();
		r.flickMouse();
		Robo.wait(5000);
		r.typeString(MovieMaker.BMPtoAVICommand("fft", 0, 1535, 20));
		r.pressEnter();
		
	}
	public static void makeSubFFTData()
	{
		long t = System.currentTimeMillis();
		double[][] data = load("C:\\data\\todan.txt");
		System.out.println("" + (System.currentTimeMillis() - t) + "to load");
		t = System.currentTimeMillis();
		FFT2D_Subset fft = new FFT2D_Subset(data, 1024);
		int fftsize = 1024;
		
		int n = 6; //lattice vectors.
		
		double[][] out = new double[fftsize][fftsize];
		
		//places to search in the fft.
		int xmin[] = new int [] {282, 172, 140, 215, 326, 361};
		int xmax[] = new int [] {297, 184, 152, 230, 341, 372};
		int ymin[] = new int [] {189, 203, 266, 310, 297, 239};
		int ymax[] = new int [] {203, 214, 274, 321, 311, 246};
		
		for (int i = 0; i < n; i++)
		{
			xmin[i] += 256;
			xmax[i] += 256;
			ymin[i] += 256;
			ymax[i] += 256;
		}
		//How many fft's? try 65,536.
		int x, y;
		int npts = 256;
		int datasize = data.length;
		int dr = (datasize - fftsize)/npts;
		double[][][] veccomps = new double[2*n][npts][npts];
		double[] temp = new double[2];
		int m;
		for (int i = 0; i < npts; i++){
			System.out.println();

			for (int j = 0; j < npts; j++)
			{
				x = dr*i; y = dr*j;
				fft.doFFT(x, y);
				FieldOps.magnitude(fft.fHat2, out);
				
				System.out.print("(" + x + ", " + y + ") ");
				
				for (int k = 0; k < n; k++)
				{
					m = 2*k;
					temp = CentroidField.centroid(out, xmin[k], xmax[k]-xmin[k], ymin[k], ymax[k]-ymin[k], false);
					veccomps[m][i][j] = temp[0]; veccomps[m+1][i][j] = temp[1];
				}
				//			makeFFTPicture(data, 2*i, 1750, fft, out, i);
			}
		}
		System.out.println();
		System.out.println("" + (System.currentTimeMillis() - t) + " to do");
		t = System.currentTimeMillis();
		
		String[] names = new String[] {"X1", "Y1", "X2", "Y2", "X3", "Y3", "X4", "Y4", "X5", "Y5", "X6", "Y6"};
		for (int i = 0; i < 2*n; i++)
		{
			save(veccomps[i], "C:\\MinGW\\data\\veccomps\\" + names[i] + ".txt");
		}
		System.out.println("" + (System.currentTimeMillis() - t) + " to save.");
		t = System.currentTimeMillis();
	}
	public static void printSubFFTPeaks()
	{
		long t = System.currentTimeMillis();
		double[][] data = load("C:\\data\\todan.txt");
		System.out.println("" + (System.currentTimeMillis() - t) + "to load");
		FFT2D_Subset fft = new FFT2D_Subset(data, 1024);
		double[][] out = new double[1024][1024];
		
		int xmin[] = new int [] {282, 172, 140, 215, 326, 361};
		int xmax[] = new int [] {297, 184, 152, 230, 341, 372};
		int ymin[] = new int [] {189, 203, 266, 310, 297, 239};
		int ymax[] = new int [] {203, 214, 274, 321, 311, 246};

		for (int i = 0; i < xmin.length; i++)
		{
			xmin[i] += 256;
			xmax[i] += 256;
			ymin[i] += 256;
			ymax[i] += 256;
		}

		for (int i = 0; i < 1536; i++)
		{
			printPeakMovement(data, 2*i, 1750, fft, out, i, xmin, xmax, ymin, ymax);
		}
		
//		Robo r = new Robo();
//		r.flickMouse();
//		Robo.wait(5000);
//		r.typeString(MovieMaker.BMPtoAVICommand("fft", 0, 1535, 20));
//		r.pressEnter();
		
	}
	
	public static void readRLatticeVectors()
	{
		double[][] data = load("C:\\MinGW\\data\\latticecut.txt");
		System.out.println(data.length + "\t" + data[0].length);
		//theta[i][j] is the angle between the ith and the jth reciprocal lattice vectors.
		int n = data.length/2;
		int N = data[0].length;
		double[][] theta = new double[n][n];
		
		double mag1, mag2, dot;
		double[] mag = new double[n];
		double[] thetax = new double[n];
		for (int i = 0; i < data[0].length; i++)
		{
			for (int j = 0; j < n; j++)
			{
				mag[j] = Math.sqrt(data[j*2][i]*data[j*2][i] + data[j*2 + 1][i]*data[j*2 + 1][i]);
				thetax[j] = DataManip.atan(data[j*2][i], data[j*2 + 1][i])*180/Math.PI;
				for (int k = 0; k < n; k++)
				{
					dot = data[j*2][i]*data[k*2][i] + data[j*2 + 1][i]*data[k*2 + 1][i];
					theta[j][k] = Math.acos(dot/mag[j]*mag[k]);
					theta[j][k] *= 180/Math.PI;
				}
				System.out.print(mag[j] + "\t" + thetax[j] + "\t");
			}
			System.out.println();
		}
	}

	public static void makeFFTPicture(double[][] data, int x, int y, FFT2D_Subset fft, double[][] mag, int index)
	{
		fft.doFFT(x, y);
		FieldOps.magnitude(fft.fHat2, mag);
		FieldOps.log(mag);
//		SRAW.writeImage("C:\\Program Files\\BMP to AVI\\fft" + MovieMaker.fromInt(index) + ".bmp", mag, 256, 256, 512, 512);
		SRAW.writeImage("C:\\Program Files\\BMP to AVI\\fft" + MovieMaker.fromInt(index) + ".bmp", mag);
	}
	public static void printPeakMovement(double[][] data, int x, int y, FFT2D_Subset fft, double[][] mag, int index, int[] xmin, int[] xmax, int[] ymin, int[] ymax)
	{
		fft.doFFT(x, y);
		FieldOps.magnitude(fft.fHat2, mag);
		printPeaks(mag, xmin.length, xmin, xmax, ymin, ymax);
	}
	
	
	public static void printPeaks(double[][] field, int npeaks, int[] xmin, int[] xmax, int[] ymin, int[] ymax)
	{
		double[] mean = new double[2], sigma = new double[2], msq = new double[2];
		String print = "";
		for (int i = 0; i < npeaks; i++)
		{
			mean = CentroidField.centroid(field, xmin[i], xmax[i]-xmin[i], ymin[i], ymax[i]-ymin[i], false);
			msq = CentroidField.cenSq(field, xmin[i], xmax[i]-xmin[i], ymin[i], ymax[i]-ymin[i], false);
			sigma[0] = Math.sqrt(msq[0] - mean[0]*mean[0]);
			sigma[1] = Math.sqrt(msq[1] - mean[1]*mean[1]);
			print += mean[0] + "\t" + mean[1] + "\t" + sigma[0] + "\t" + sigma[1] + "\t\t";
		}
		System.out.println(print);
	}
	public static void makeCentroidField()
	{
		String txt = "C:\\MinGW\\data\\";
		long t = System.currentTimeMillis();

//		double[][] real = load(txt + "retap.txt");
//		double[][] imag = load(txt + "imtap.txt");
//		System.out.println("" + (System.currentTimeMillis() - t) + " to load.");
//		t = System.currentTimeMillis();
//		int N = 512;
//		double[][] data = new double[N][N];
//		for (int i = 0; i < N; i++)
//			for(int j = 0; j < N; j++)
//				data[i][j] = real[i][j]*real[i][j] + imag[i][j]*imag[i][j];

		double[][] data = load(txt + "magsquaretap3.txt");
		System.out.println("" + (System.currentTimeMillis() - t) + " to load.");
		t = System.currentTimeMillis();
		double[][] field = CentroidField.centroidField2(data, 0, 48, 12);

		System.out.println("" + (System.currentTimeMillis() - t) + " to do.");
		t = System.currentTimeMillis();

		save(field, txt + "mags3 64 6.txt");
		System.out.println("" + (System.currentTimeMillis() - t) + " to save.");
		t = System.currentTimeMillis();
	}
	public static void smoothOut()
	{
		String txt = "C:\\MinGW\\data\\";
		long t = System.currentTimeMillis();
		double[][] real = load(txt + "retap.txt");
		double[][] imag = load(txt + "imtap.txt");
		System.out.println("" + (System.currentTimeMillis() - t) + " to load.");
		t = System.currentTimeMillis();
		int N = 512;
		double[][] data = new double[N][N];
		for (int i = 0; i < N; i++)
			for(int j = 0; j < N; j++)
				data[i][j] = real[i][j]*real[i][j] + imag[i][j]*imag[i][j];
		
		CentroidField.smooth(data, 3);
		
		System.out.println("" + (System.currentTimeMillis() - t) + " to do.");
		t = System.currentTimeMillis();
		
		save(data, txt + "magsquaretap3.txt");
		System.out.println("" + (System.currentTimeMillis() - t) + " to save.");
		t = System.currentTimeMillis();
	}
	public static void subtractAverage()
	{
		String txt = "C:\\MinGW\\data\\";
		long t = System.currentTimeMillis();
		double[][] data = load(txt + "table.txt");
		System.out.println("" + (System.currentTimeMillis() - t) + " to load.");
		t = System.currentTimeMillis();

		subtractAvg(data);
		
		System.out.println("" + (System.currentTimeMillis() - t) + " to do.");
		t = System.currentTimeMillis();

		save(data, txt + "tableprime.txt");
		System.out.println("" + (System.currentTimeMillis() - t) + " to save.");
		t = System.currentTimeMillis();
	}

	//This returns a 2-dimensional complex array: [nx][ny][2].
	//In this method we sum over a rectangle of size L, instead of 
	//summing over a larger area with a decaying Gaussian factor of width L, 
	//as in the paper. We take all big arrays as parameters.
	public static double[][][] getDeviation(int N, int L, double[] q, int xmin, int dx, int nx, int ymin, int dy, int ny, double[][] topo, double[][][] expu, double[][][] Texpmr)
	{
		//q is the Bragg peak position, for the entire topograph. It is expressed
		//as a pixel displacement from the origin of the FFT.
	
		//Our unit system is "natural": The unit of length is 1 pixel; the size of the 
		//topograph is then NxN. The correct physical value of q is		
		double[] qright = {2*Math.PI*q[0]/N, 2*Math.PI*q[1]/N};
		
		//expd is e^(i*q*u(r)) according to the paper.
//		double[][] Tmag, Tphase;
		int x, y;
		double dot = 0;
		double[] cmplx = new double[2];
		int i, j, m, n;
		for (i = 0; i < nx; i++){System.out.print(i + " ");
			for (j = 0; j < ny; j++)
			{
				expu[i][j][0] = 0; expu[i][j][1] = 0;
				x = xmin + i*dx;
				y = ymin + j*dy;
				for (m = 0; m < L; m++)
					for (n = 0; n < L; n++)
					{
						dot = dot(qright, x+m, y+n);
						Complex.expmi(dot, cmplx);
//						System.out.println(dot + ", " + cmplx[0]);
						Complex.times(cmplx, topo[x+m][y+n], cmplx);
						Texpmr[m][n][0] = cmplx[0];
						Texpmr[m][n][1] = cmplx[1];
						Complex.sum(expu[i][j], cmplx, expu[i][j]);
					}
//				Tmag = FieldOps.magnitude(Texpmr);
//				Tphase = FieldOps.phase(Texpmr);
//				double max = ArrayOps.max(Tmag);
//				double min = ArrayOps.min(Tmag);
//				System.out.print("[" + min + ", " + max + "]\t");
//				SRAW.writeImage("C:\\data\\lawlerhex\\kd"+"_"+i+"_"+j+".bmp", Tmag, Tphase, new ColorScales.MYC2d(max, min, 2*Math.PI));
			}
		}
		System.out.println();
		return expu;
	}
	public static double[][][] getDeviation(int L, double[][][] cosmsintopo, int xmin, int dx, int nx, int ymin, int dy, int ny, double[][][] expu)
	{
		//q is the Bragg peak position, for the entire topograph. It is expressed
		//as a pixel displacement from the origin of the FFT.
	
		//Our unit system is "natural": The unit of length is 1 pixel; the size of the 
		//topograph is then NxN. The correct physical value of q is		
		
		//expd is e^(i*q*u(r)) according to the paper.
	//	double[][] Tmag, Tphase;
		int x, y;
		int i, j, m, n;
		for (i = 0; i < nx; i++){/*System.out.print(i + " ");*/
			for (j = 0; j < ny; j++)
			{
				expu[i][j][0] = 0; expu[i][j][1] = 0;
				x = xmin + i*dx;
				y = ymin + j*dy;
				for (m = 0; m < L; m++)
					for (n = 0; n < L; n++)
					{
						Complex.sum(expu[i][j], cosmsintopo[x+m][y+n], expu[i][j]);
					}
			}
		}
//		System.out.println();
		return expu;
	}
	public static double[][][] getDeviationGauss(int L, double[][][] cosmsintopo, int xmin, int dx, int nx, int ymin, int dy, int ny, double[][][] expu, double[][] gauss)
	{
		//L is now the actual width of the Gaussian. We will go from -3L to +3L. We assume
		//that the gaussian array is filled out at least that much.
		
		//q is the Bragg peak position, for the entire topograph. It is expressed
		//as a pixel displacement from the origin of the FFT.
	
		//Our unit system is "natural": The unit of length is 1 pixel; the size of the 
		//topograph is then NxN. The correct physical value of q is		
		
		//expd is e^(i*q*u(r)) according to the paper.
	//	double[][] Tmag, Tphase;
		int x, y;
		int N = cosmsintopo.length;
		double[] cmplx = new double[2];
		int i, j;
		int xprime, yprime;
		int xpmin, xpmax, ypmin, ypmax;
		double gsum = 0;
		for (i = 0; i < nx; i++)/*{System.out.print(i + " ");*/
			for (j = 0; j < ny; j++)
			{
				gsum = 0;
				expu[i][j][0] = 0; expu[i][j][1] = 0;
				x = xmin + i*dx;
				y = ymin + j*dy;
				//this is done so that xprime runs from 0 to N, or from x - biggerL/2 to x+biggerL/2
				xpmin = Math.max(0, x - gauss.length);
				xpmax = Math.min(N, x + gauss.length);
				ypmin = Math.max(0, y - gauss.length);
				ypmax = Math.min(N, y + gauss.length);
				for (xprime = xpmin; xprime < xpmax; xprime++)
					for (yprime = ypmin; yprime < ypmax; yprime++)
					{
						gsum += gauss[Math.abs(x-xprime)][Math.abs(y-yprime)];
						Complex.times(cosmsintopo[xprime][yprime], gauss[Math.abs(x-xprime)][Math.abs(y-yprime)], cmplx);
						Complex.sum(expu[i][j], cmplx, expu[i][j]);
					}
				expu[i][j][0] /= gsum; expu[i][j][1] /= gsum;
	
			}
//		}
//		System.out.println();
		return expu;
	}
	public static double[][][] getDeviationGaussDefault(double L, double[][][] cosmsintopo, double[][][] expu, double[][] gauss)
	{
		//L is now the actual width of the Gaussian. We will go from -3L to +3L. We assume
		//that the gaussian array is filled out at least that much.
		
		//q is the Bragg peak position, for the entire topograph. It is expressed
		//as a pixel displacement from the origin of the FFT.
	
		//Our unit system is "natural": The unit of length is 1 pixel; the size of the 
		//topograph is then NxN. The correct physical value of q is		
		
		//expd is e^(i*q*u(r)) according to the paper.
	//	double[][] Tmag, Tphase;
		int x, y;
		int N = cosmsintopo.length;
		double[] cmplx = new double[2];
		int i, j;
		int xprime, yprime;
		int xpmin, xpmax, ypmin, ypmax;
		double gsum = 0;
		int nx = cosmsintopo.length;
		int ny = cosmsintopo[0].length;
		for (i = 0; i < nx; i++)/*{System.out.print(i + " ");*/
			for (j = 0; j < ny; j++)
			{
				gsum = 0;
				expu[i][j][0] = 0; expu[i][j][1] = 0;
				x = i;
				y = j;
				//this is done so that xprime runs from 0 to N, or from x - biggerL/2 to x+biggerL/2
				xpmin = FieldOps.round(Math.max(0, x - 3*L - 1));
				xpmax = FieldOps.round(Math.min(N, x + 3*L + 1));
				ypmin = FieldOps.round(Math.max(0, y - 3*L - 1));
				ypmax = FieldOps.round(Math.min(N, y + 3*L + 1));
				for (xprime = xpmin; xprime < xpmax; xprime++)
					for (yprime = ypmin; yprime < ypmax; yprime++)
					{
						gsum += gauss[Math.abs(x-xprime)][Math.abs(y-yprime)];
						Complex.times(cosmsintopo[xprime][yprime], gauss[Math.abs(x-xprime)][Math.abs(y-yprime)], cmplx);
						Complex.sum(expu[i][j], cmplx, expu[i][j]);
					}
				expu[i][j][0] /= gsum; expu[i][j][1] /= gsum;
	
			}
//		}
//		System.out.println();
		return expu;
	}
	
	//This method does the entire fft and finds the bragg peaks, by taking
	//the center of mass of the FFT within n pixels of the indicated Bragg peaks.
	public static void doDeviationJob(String dir, String name, int[][] braggi, int n, int N, int[] r1, int[] dr, int[]npts, int L)
	{
		long t = System.currentTimeMillis();

		double[][] data = ColumnIO.readSquareTable(dir + name + ".dat" );
		System.out.println((System.currentTimeMillis() - t)/1000 + " to load.");
		t = System.currentTimeMillis();
		
		double[][] mag = FFTOps.obtainFFTmagCent(data);
		t = printTime(" to fft", t);

		double[][] bragg = new double[2][2];
		for (int i = 0; i < 2; i++){
			bragg[i] = CentroidField.centroid(mag, braggi[i][0]-n, 2*n + 1, braggi[i][1]-n, 2*n + 1, false);
			System.out.println("(" + braggi[i][0] + ", " + braggi[i][1] + ")  ----> (" + bragg[i][0] + ", " + bragg[i][1] + ")" );
			bragg[i][0] -= N/2;
			bragg[i][1] -= N/2;
		}
		
		double[][] braggTrue = new double[2][2];
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				braggTrue[i][j] = bragg[i][j]*2*Math.PI/N;
		
		//We store the trigonometric factors in memory to save computing time:
		double[][][] cosmsintopo = new double[N][N][2];
//		System.out.println(bragg[0] + ", " + bragg[1]);
//		System.out.println(qright[0] + ", " + qright[1]);
		double[][][] dev = new double[npts[0]][npts[1]][2];
		double[][] phase = new double[npts[0]][npts[1]];
		mag = new double[npts[0]][npts[1]];
		for (int i = 1; i < 3; i++)
		{
			//set the table values.
			for (int j = 0; j < N; j++)
				for (int k = 0; k < N; k++){
					cosmsintopo[j][k][0] = Math.cos(dot(braggTrue[i-1], j, k))*data[j][k];
					cosmsintopo[j][k][1] = -Math.sin(dot(braggTrue[i-1], j, k))*data[j][k];
				}
			SRAW.writeImage(dir + name + "wave" + i + ".bmp", cosmsintopo, false);
			getDeviation(L, cosmsintopo, r1[0], dr[0], npts[0], r1[1], dr[1], npts[1], dev);
			t = printTime(" to do.", t);
			
//			mag = FieldOps.magnitude(dev);
			FieldOps.magnitude(dev, mag);
			FieldOps.phase(dev, phase);
			double max = ArrayOps.max(mag);
//			System.out.println(max + ", " + avg);
			SRAW.writeImage(dir + name + "dev" + i + ".bmp", mag, phase, new ColorScales.MYC2d(max, 0, 2*Math.PI));
			ColumnIO.writeTable(mag, dir + name + "devmag" + i + ".txt");
			ColumnIO.writeTable(phase, dir + name + "devphase" + i + ".txt");
			//For visual purposes.
			
		}
		
	}
	public static void makeDevMovie(String dir, String name, int N, int[] r1, int[] dr, int[]npts, int[] L)
	{
		long t = System.currentTimeMillis();

		double[][] data = ColumnIO.readSquareTable(dir + name + ".dat" );
		System.out.println((System.currentTimeMillis() - t)/1000 + " to load.");
		t = System.currentTimeMillis();
		
		double[][] mag = FFTOps.obtainFFTmagCent(data);
		t = printTime(" to fft", t);

		int[][] braggi = readBragg(dir + "bragg.txt");
		int n = readBraggSmear(dir + "bragg.txt");
		double[][] bragg = new double[2][2];
		for (int i = 0; i < 2; i++){
			bragg[i] = CentroidField.centroid(mag, braggi[i][0]-n, 2*n + 1, braggi[i][1]-n, 2*n + 1, false);
			System.out.println("(" + braggi[i][0] + ", " + braggi[i][1] + ")  ----> (" + bragg[i][0] + ", " + bragg[i][1] + ")" );
			bragg[i][0] -= N/2;
			bragg[i][1] -= N/2;
		}
		
		double[][] braggTrue = new double[2][2];
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				braggTrue[i][j] = bragg[i][j]*2*Math.PI/N;
		
		//We store the trigonometric factors in memory to save computing time:
		double[][][] cosmsintopo = new double[N][N][2];
//		System.out.println(bragg[0] + ", " + bragg[1]);
//		System.out.println(qright[0] + ", " + qright[1]);
		double[][][] dev = new double[npts[0]][npts[1]][2];
		double[][] phase = new double[npts[0]][npts[1]];
		mag = new double[npts[0]][npts[1]];
		String marker;
		String dir2 = "C:\\Program Files\\BMP to AVI\\";
		for (int i = 1; i < 3; i++)
		{
			marker = i==1 ? "x" : "y";
			//set the table values.
			for (int j = 0; j < N; j++)
				for (int k = 0; k < N; k++){
					cosmsintopo[j][k][0] = Math.cos(dot(braggTrue[i-1], j, k))*data[j][k];
					cosmsintopo[j][k][1] = -Math.sin(dot(braggTrue[i-1], j, k))*data[j][k];
				}
			SRAW.writeImage(dir + name + "wave" + i + ".bmp", cosmsintopo, false);
			for (int j = 0; j < L.length; j++)
			{
				getDeviation(L[j], cosmsintopo, r1[0], dr[0], npts[0], r1[1], dr[1], npts[1], dev);
//				t = printTime(" to do.", t);
				
	//			mag = FieldOps.magnitude(dev);
				FieldOps.magnitude(dev, mag);
				FieldOps.phase(dev, phase);
				double max = ArrayOps.max(mag);
	//			System.out.println(max + ", " + avg);
				SRAW.writeImage(dir2 + name + "dev" + marker + MovieMaker.fromInt(j) + ".bmp", dev, true);
//				ColumnIO.writeTable(mag, dir + name + "devmag" + i + ".txt");
//				ColumnIO.writeTable(phase, dir + name + "devphase" + i + ".txt");
				//For visual purposes.
			}
			t = printTime(" to do.", t);
			MovieMaker.typeBMPtoAVICommand(name+ "dev" + marker, 0, L.length - 1, 20);
		}
		
	}
	public static double[][] getNormExactBragg(String dir, double[][] data)
	{
		int N = data.length;
		double[][] mag = FFTOps.obtainFFTmagCent(data);
		int[][] braggi = readBragg(dir + "bragg.txt");
		int n = readBraggSmear(dir + "bragg.txt");
		double[][] bragg = new double[2][2];
		for (int i = 0; i < 2; i++){
			bragg[i] = CentroidField.centroid(mag, braggi[i][0]-n, 2*n + 1, braggi[i][1]-n, 2*n + 1, false);
			System.out.println("(" + braggi[i][0] + ", " + braggi[i][1] + ")  ----> (" + bragg[i][0] + ", " + bragg[i][1] + ")" );
			bragg[i][0] -= N/2;
			bragg[i][1] -= N/2;
		}
		
		double[][] braggTrue = new double[2][2];
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				braggTrue[i][j] = bragg[i][j]*2*Math.PI/N;
		
		return braggTrue;
	}
	public static void makeDevMovieGauss(String dir, String name, int[] r1, int[] dr, int[]npts, int[] L, int Lcutoff)
	{
		long t = System.currentTimeMillis();
		
		double[][] data = ColumnIO.readSquareTable(dir + name + ".dat" );
		int N = data.length;
		System.out.println((System.currentTimeMillis() - t)/1000 + " to load.");
		t = System.currentTimeMillis();
		
		double[][] mag = FFTOps.obtainFFTmagCent(data);
		t = printTime(" to fft", t);

		int[][] braggi = readBragg(dir + "bragg.txt");
		int n = readBraggSmear(dir + "bragg.txt");
		double[][] bragg = new double[2][2];
		for (int i = 0; i < 2; i++){
			bragg[i] = CentroidField.centroid(mag, braggi[i][0]-n, 2*n + 1, braggi[i][1]-n, 2*n + 1, false);
			System.out.println("(" + braggi[i][0] + ", " + braggi[i][1] + ")  ----> (" + bragg[i][0] + ", " + bragg[i][1] + ")" );
			bragg[i][0] -= N/2;
			bragg[i][1] -= N/2;
		}
		
		double[][] braggTrue = new double[2][2];
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				braggTrue[i][j] = bragg[i][j]*2*Math.PI/N;
		
		writeBraggTrue(dir + "bragg.txt", braggTrue);
		
		int biggerL = L[L.length - 1]*4;
		double[][] gauss = new double[biggerL][biggerL];
		
		//We store the trigonometric factors in memory to save computing time:
		double[][][] cosmsintopo = new double[N][N][2];
//		System.out.println(bragg[0] + ", " + bragg[1]);
//		System.out.println(qright[0] + ", " + qright[1]);
		double[][][] dev = new double[npts[0]][npts[1]][2];
		String marker;
		String dir2 = "C:\\Users\\madhavanlab2011\\Documents\\BMP to AVI\\";
		for (int i = 2; i < 3; i++)
		{
			marker = i==1 ? "x" : "y";
			//set the table values.
			for (int j = 0; j < N; j++)
				for (int k = 0; k < N; k++){
					cosmsintopo[j][k][0] = Math.cos(dot(braggTrue[i-1], j, k))*data[j][k];
					cosmsintopo[j][k][1] = -Math.sin(dot(braggTrue[i-1], j, k))*data[j][k];
				}
			SRAW.writeImage(dir + name + "wave" + marker + ".bmp", cosmsintopo, false);
			for (int j = 0; j < L.length; j++)
			{
				for (int k = 0; k < biggerL; k++)
					for (int l = 0; l < biggerL; l++)
						gauss[k][l] = Math.exp((-((double)(k*k) + (l*l))/(2*L[j]*L[j])));
				getDeviationGauss(L[j], cosmsintopo, r1[0], dr[0], npts[0], r1[1], dr[1], npts[1], dev, gauss);
//				t = printTime(" to do.", t);
				
	//			mag = FieldOps.magnitude(dev);
	//			System.out.println(max + ", " + avg);
				SRAW.writeImage(dir + "devout\\" + name + "dev" + marker + "L" + L[j] + ".bmp", dev, true);
				if (L[j] > Lcutoff)
					ColumnIO.writeBinsPolar(dev, dir + name + marker + "L" + L[j]);
				//For visual purposes.
			}
			t = printTime(" to do.", t);
//			MovieMaker.typeBMPtoAVICommand(name+ "dev" + marker, 0, L.length - 1, 10);
		}
		
	}
	public static void doUFieldSummation(String dir, String name, int detail, int[] L, int Lcutoff, int xory)
	{
		long t = System.currentTimeMillis();
		
		double[][] data = ColumnIO.readSquareTable(dir + name + ".dat" );
		int N = data.length;
		int npts = N/detail;
		System.out.println((System.currentTimeMillis() - t)/1000 + " to load.");
		t = System.currentTimeMillis();
		
		double[][] mag = FFTOps.obtainFFTmagCent(data);
		t = printTime(" to fft", t);

		int[][] braggi = readBragg(dir + "bragg.txt");
		int n = readBraggSmear(dir + "bragg.txt");
		double[][] bragg = new double[2][2];
		for (int i = 0; i < 2; i++){
			bragg[i] = CentroidField.centroid(mag, braggi[i][0]-n, 2*n + 1, braggi[i][1]-n, 2*n + 1, false);
			System.out.println("(" + braggi[i][0] + ", " + braggi[i][1] + ")  ----> (" + bragg[i][0] + ", " + bragg[i][1] + ")" );
			bragg[i][0] -= N/2;
			bragg[i][1] -= N/2;
		}
		
		double[][] braggTrue = new double[2][2];
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				braggTrue[i][j] = bragg[i][j]*2*Math.PI/N;
		
		writeBraggTrue(dir + "bragg.txt", braggTrue);
		
		int biggerL = L[L.length - 1]*4;
		double[][] gauss = new double[biggerL][biggerL];
		
		//We store the trigonometric factors in memory to save computing time:
		double[][][] cosmsintopo = new double[N][N][2];
//		System.out.println(bragg[0] + ", " + bragg[1]);
//		System.out.println(qright[0] + ", " + qright[1]);
		double[][][] dev = new double[npts][npts][2];
		String marker;
		String dirdevout = dir + "devout\\";
		String dirdevoutbin = dir + "devoutbin\\";
		if (!new File(dirdevout).exists())
			new File(dirdevout).mkdir();
		if (!new File(dirdevoutbin).exists())
			new File(dirdevoutbin).mkdir();
		
		for (int i = 1+xory; i < 2+xory; i++)
		{
			marker = i==1 ? "x" : "y";
			//set the table values.
			for (int j = 0; j < N; j++)
				for (int k = 0; k < N; k++){
					cosmsintopo[j][k][0] = Math.cos(dot(braggTrue[i-1], j, k))*data[j][k];
					cosmsintopo[j][k][1] = -Math.sin(dot(braggTrue[i-1], j, k))*data[j][k];
				}
			ColumnIO.writeBinsPolar(cosmsintopo, dirdevoutbin + name + "wave" + marker);
			SRAW.writeImage(dirdevout + name + "wave" + marker + ".bmp", cosmsintopo, false);
			for (int j = 0; j < L.length; j++)
			{
				for (int k = 0; k < biggerL; k++)
					for (int l = 0; l < biggerL; l++)
						gauss[k][l] = Math.exp((-((double)(k*k) + (l*l))/(2*L[j]*L[j])));
				getDeviationGauss(L[j], cosmsintopo, 0, detail, npts, 0, detail, npts, dev, gauss);
//				t = printTime(" to do.", t);
				
	//			mag = FieldOps.magnitude(dev);
	//			System.out.println(max + ", " + avg);
				t = printTimeNL(" to write L=" + L[j] + "\t", t);
				SRAW.writeImage(dirdevout + name + marker + "L" + L[j] + ".bmp", dev, true);
				if (L[j] > Lcutoff)
					ColumnIO.writeBinsPolar(dev, dirdevoutbin + name + marker + "L" + L[j]);
				//For visual purposes.
			}
//			MovieMaker.typeBMPtoAVICommand(name+ "dev" + marker, 0, L.length - 1, 10);
		}
		
	}
	public static void doUFieldSummationFT(String dir, String name, double[] L, double xi, boolean do3bragg)
	{
		long t = System.currentTimeMillis();
		
		double[][] data = ColumnIO.readSquareTable(dir + name + ".dat" );
		int N = data.length;
		System.out.println((System.currentTimeMillis() - t)/1000 + " to load.");
		t = System.currentTimeMillis();
		
		double[][][] fftz = new double [N][N][2];
		double[][][] nfftz = new double [N][N][2];
		FFTOps.obtainFFTCent(data, fftz);
//		double[][] braggTrue = getBraggTrue(dir, data, do3bragg);
		double[][] braggTrue = getIntBraggTrue(dir, data, do3bragg);
//		double[][] braggTrue = getSymIntBraggTrue(dir, data, do3bragg, Math.PI/3);
		int nb = braggTrue.length;
		t = printTime(" to fft", t);
		//this time we gaussianly suppress the fft and see what happens.
		double[][] gauss = new double[N][N];
				//We store the trigonometric factors in memory to save computing time:
		double[][][] cosmsintopo = new double[N][N][2];
//		System.out.println(bragg[0] + ", " + bragg[1]);
//		System.out.println(qright[0] + ", " + qright[1]);
		double[][][] dev = new double[N][N][2];
		String marker;
		String dirdevout = dir + "devoutft\\";
		String dirdevoutbin = dir + "devoutbinft\\";
		if (!new File(dirdevout).exists())
			new File(dirdevout).mkdir();
		if (!new File(dirdevoutbin).exists())
			new File(dirdevoutbin).mkdir();
		
		for (int i = 1; i < nb + 1; i++)
		{
			if (i == 1) marker = "x";
			else if (i == 2) marker = "y";
			else marker = "z";
			//set the table values.
			for (int j = 0; j < N; j++)
				for (int k = 0; k < N; k++){
					cosmsintopo[j][k][0] = Math.cos(dot(braggTrue[i-1], j, k))*data[j][k];
					cosmsintopo[j][k][1] = -Math.sin(dot(braggTrue[i-1], j, k))*data[j][k];
				}
			FFTOps.obtainFFTCent(cosmsintopo, fftz);
//			ColumnIO.writeBinsPolar(cosmsintopo, dirdevoutbin + name + "wave" + marker);
			SRAW.writeImage(dirdevout + name + "wave" + marker + ".bmp", cosmsintopo, false);
			SRAW.writeImage(dirdevout + name + "wave" + marker + "fft.bmp", fftz, false, true);
			for (int j = 0; j < L.length; j++)
			{
//				TwoDArray.getGaussMask(L[j], N, gauss);
				TwoDArray.getCrossMask(L[j], xi, N, gauss);

				for (int k = 0; k < N; k++)
					for (int l = 0; l < N; l++)
					{
						nfftz[k][l][0] = fftz[k][l][0]*gauss[Math.abs(k-N/2)][Math.abs(l-N/2)];
						nfftz[k][l][1] = fftz[k][l][1]*gauss[Math.abs(k-N/2)][Math.abs(l-N/2)];
					}
				FFTOps.obtainIFFTCent(nfftz, dev);
				
				//				t = printTime(" to do.", t);
				
	//			mag = FieldOps.magnitude(dev);
	//			System.out.println(max + ", " + avg);
				t = printTimeNL(" to write L=" + L[j] + "\t", t);
				SRAW.writeImage(dirdevout + name + marker + "L" + L[j] + ".bmp", dev, true);
				SRAW.writeImagPhase(dirdevout + name + marker + "L" + L[j] + "phase.bmp", dev);
				//solely for image output:
				for (int m = 0;  m < N; m++)
					for (int p = 0; p < N; p++)
						if (Complex.mag(nfftz[m][p]) == 0) nfftz[m][p][0] = Double.MIN_VALUE;
				SRAW.writeImage(dirdevout + name + marker + "L" + L[j] + "fft.bmp", nfftz, true, true);
//				ColumnIO.writeBinsPolar(dev, dirdevoutbin + name + "dev" + marker + "L" + L[j]);
				//For visual purposes.
			}
//			MovieMaker.typeBMPtoAVICommand(name+ "dev" + marker, 0, L.length - 1, 10);
		}
		
	}
	public static double[][] getBraggTrue(String dir, double[][] data, boolean do3bragg)
	{
		int N = data.length;
		
		double[][] mag = FFTOps.obtainFFTmagCent(data);
		int[][] braggi;
		braggi = do3bragg ? readBragg3(dir + "bragg.txt") : readBragg(dir + "bragg.txt");
		int nb = do3bragg ? 3 : 2;
		int n = do3bragg ? readBraggSmear3(dir + "bragg.txt") : readBraggSmear(dir + "bragg.txt");
		double[][] bragg = new double[nb][2];
		for (int i = 0; i < nb; i++){
			bragg[i] = CentroidField.centroid(mag, braggi[i][0]-n, 2*n + 1, braggi[i][1]-n, 2*n + 1, false);
			System.out.println("(" + braggi[i][0] + ", " + braggi[i][1] + ")  ----> (" + bragg[i][0] + ", " + bragg[i][1] + ")" );
			bragg[i][0] -= N/2;
			bragg[i][1] -= N/2;
		}
		
		double[][] braggTrue = new double[nb][2];
		for (int i = 0; i < nb; i++)
			for (int j = 0; j < 2; j++)
				braggTrue[i][j] = bragg[i][j]*2*Math.PI/N;
		
		writeBraggTrue(dir + "bragg.txt", braggTrue);
		return braggTrue;
	}
	public static double[][] getIntBraggTrue(String dir, double[][] data, boolean do3bragg)
	{
		int N = data.length;
		
		double[][] mag = FFTOps.obtainFFTmagCent(data);
		int[][] braggi;
		braggi = do3bragg ? readBragg3(dir + "bragg.txt") : readBragg(dir + "bragg.txt");
		int nb = do3bragg ? 3 : 2;
		int n = do3bragg ? readBraggSmear3(dir + "bragg.txt") : readBraggSmear(dir + "bragg.txt");
		double[][] bragg = new double[nb][2];
		double[] temp = new double [2];
		for (int i = 0; i < nb; i++){
			temp = CentroidField.centroid(mag, braggi[i][0]-n, 2*n + 1, braggi[i][1]-n, 2*n + 1, false);
			bragg[i][0] = (int)(temp[0]+0.5);
			bragg[i][1] = (int)(temp[1]+0.5);
			System.out.println("(" + braggi[i][0] + ", " + braggi[i][1] + ")  ----> (" + temp[0] + ", " + temp[1] + ")  ----> (" + bragg[i][0] + ", " + bragg[i][1] + ")");
			bragg[i][0] -= N/2;
			bragg[i][1] -= N/2;
		}
		
		double[][] braggTrue = new double[nb][2];
		for (int i = 0; i < nb; i++)
			for (int j = 0; j < 2; j++)
				braggTrue[i][j] = bragg[i][j]*2*Math.PI/N;
		
		writeBraggTrue(dir + "bragg.txt", braggTrue);
		return braggTrue;
	}
	public static double[][] getSymIntBraggTrue(String dir, double[][] data, boolean do3bragg, double theta)
	{
		int N = data.length;
		
		double[][] mag = FFTOps.obtainFFTmagCent(data);
		int[][] braggi;
		braggi = do3bragg ? readBragg3(dir + "bragg.txt") : readBragg(dir + "bragg.txt");
		int nb = do3bragg ? 3 : 2;
		int n = do3bragg ? readBraggSmear3(dir + "bragg.txt") : readBraggSmear(dir + "bragg.txt");
		double[][] bragg = new double[nb][2];
		for (int i = 0; i < 1; i++){
			bragg[i] = CentroidField.centroid(mag, braggi[i][0]-n, 2*n + 1, braggi[i][1]-n, 2*n + 1, false);
			System.out.println("(" + braggi[i][0] + ", " + braggi[i][1] + ")  ----> (" + bragg[i][0] + ", " + bragg[i][1] + ")" );
			bragg[i][0] = FieldOps.round(bragg[i][0]);
			bragg[i][1] = FieldOps.round(bragg[i][1]);
			bragg[i][0] -= N/2;
			bragg[i][1] -= N/2;
		}
		
		
		double[][] braggTrue = new double[nb][2];
		double bmag = Complex.mag(bragg[0]);
		double btheta = FieldOps.atan(bragg[0][0], bragg[0][1]);
		for (int i = 1; i < nb; i++)
		{
			bragg[i][0] = bmag*Math.cos(btheta - i*theta);
			bragg[i][1] = bmag*Math.sin(btheta - i*theta);
			bragg[i][0] = FieldOps.round(bragg[i][0]);
			bragg[i][1] = FieldOps.round(bragg[i][1]);
		}
		for (int i = 0; i < nb; i++)
			for (int j = 0; j < 2; j++)
				braggTrue[i][j] = bragg[i][j]*2*Math.PI/N;
		
		writeBraggTrue(dir + "bragg.txt", braggTrue);
		return braggTrue;
		
	}
	public static void doUFieldSummationFTWriteBin(String dir, String name, double[] L, double xi, boolean do3bragg)
	{
		long t = System.currentTimeMillis();
		
		double[][] data = ColumnIO.readSquareTable(dir + name + ".dat" );
		int N = data.length;
		System.out.println((System.currentTimeMillis() - t)/1000 + " to load.");
		t = System.currentTimeMillis();
		
		double[][][] fftz = new double [N][N][2];
		double[][][] nfftz = new double [N][N][2];
		FFTOps.obtainFFTCent(data, fftz);
//		double[][] braggTrue = getBraggTrue(dir, data, do3bragg);
		double[][] braggTrue = getIntBraggTrue(dir, data, do3bragg);
//		double[][] braggTrue = getSymIntBraggTrue(dir, data, do3bragg, Math.PI/3);
		int nb = braggTrue.length;
		
		//this time we gaussianly suppress the fft and see what happens.
		double[][] gauss = new double[N][N];
				//We store the trigonometric factors in memory to save computing time:
		double[][][] cosmsintopo = new double[N][N][2];
		double[][][] dev = new double[N][N][2];
		String marker;
		String dirdevout = dir + "devoutft\\";
		String dirdevoutbin = dir + "devoutbinft\\";
		if (!new File(dirdevout).exists())
			new File(dirdevout).mkdir();
		if (!new File(dirdevoutbin).exists())
			new File(dirdevoutbin).mkdir();
		
		for (int i = 1; i < nb + 1; i++)
		{
			if (L[i-1] == 0) continue;
			if (i == 1) marker = "x";
			else if (i == 2) marker = "y";
			else marker = "z";
			//set the table values.
			for (int j = 0; j < N; j++)
				for (int k = 0; k < N; k++){
					cosmsintopo[j][k][0] = Math.cos(dot(braggTrue[i-1], j, k))*data[j][k];
					cosmsintopo[j][k][1] = -Math.sin(dot(braggTrue[i-1], j, k))*data[j][k];
				}
			FFTOps.obtainFFTCent(cosmsintopo, fftz);
//			TwoDArray.getGaussMask(L[j], N, gauss);
			TwoDArray.getCrossMask(L[i-1], xi, N, gauss);
			for (int k = 0; k < N; k++)
				for (int l = 0; l < N; l++)
					{
						nfftz[k][l][0] = fftz[k][l][0]*gauss[Math.abs(k-N/2)][Math.abs(l-N/2)];
						nfftz[k][l][1] = fftz[k][l][1]*gauss[Math.abs(k-N/2)][Math.abs(l-N/2)];
					}
			FFTOps.obtainIFFTCent(nfftz, dev);
			ColumnIO.writeBinsPolar(dev, dirdevoutbin + name + "dev" + marker + "L" + L[i-1]);
			t = printTimeNL(" to write L=" + L[i-1] + "\t", t);
		}
		
	}

	public static void braggSmearTest(String dir, String name)
	{
		int[][] braggi = readBragg(dir + "bragg.txt");
		int n = readBraggSmear(dir + "bragg.txt");
		double[] bragg = new double[2];
		double[][] data = ColumnIO.readSquareTable(dir + name + ".dat" );
		int N = data.length;
		double[][] mag = FFTOps.obtainFFTmagCent(data);
		System.out.println("fft");
		int x;
		int y;
		for (int j = 0; j < 3*n; j++)
			for (int k = 0; k < 3*n; k++){
					x = braggi[0][0]-n-(3*n/2)+j;
					y = braggi[0][1]-n-(3*n/2)+k;
					bragg = CentroidField.centroid(mag, x, 2*n + 1, y, 2*n + 1, false);
					System.out.print((x+n) + "\t" + (y+n) + "\t" + bragg[0] + "\t" + bragg[1] + "\t");
					x = braggi[1][0]-n-(3*n/2)+j;
					y = braggi[1][1]-n-(3*n/2)+k;
					bragg = CentroidField.centroid(mag, x, 2*n + 1, y, 2*n + 1, false);
					System.out.println((x+n) + "\t" + (y+n) + "\t" + bragg[0] + "\t" + bragg[1] + "\t");
				}

	}
	//This is the same as the above, except a Gaussian weight factor is used instead of a sharp cutoff.
	public static void doDeviationJobGauss(String dir, String name, int[][] braggi, int n, int N, int[] r1, int[] dr, int[]npts, int L, int biggerL)
	{
		long t = System.currentTimeMillis();

		double[][] data = ColumnIO.readSquareTable(dir + name + ".dat" );
		System.out.println((System.currentTimeMillis() - t)/1000 + " to load.");
		t = System.currentTimeMillis();
		
		double[][] mag = FFTOps.obtainFFTmagCent(data);
		t = printTime(" to fft", t);

		double[][] bragg = new double[2][2];
		for (int i = 0; i < 2; i++){
			bragg[i] = CentroidField.centroid(mag, braggi[i][0]-n, 2*n + 1, braggi[i][1]-n, 2*n + 1, false);
			System.out.println("(" + braggi[i][0] + ", " + braggi[i][1] + ")  ----> (" + bragg[i][0] + ", " + bragg[i][1] + ")" );
			bragg[i][0] -= N/2;
			bragg[i][1] -= N/2;
		}
		
		double[][] braggTrue = new double[2][2];
		for (int i = 0; i < 2; i++)
			for (int j = 0; j < 2; j++)
				braggTrue[i][j] = bragg[i][j]*2*Math.PI/N;
		
		double[][] gauss = new double [biggerL][biggerL];
		for (int i = 0; i < biggerL; i++)
			for (int j = 0; j < biggerL; j++)
				gauss[i][j] = Math.exp((-((double)(i*i) + (j*j))/(2*L*L)));
		SRAW.writeImage("C:\\data\\lawlerhex\\"+name+"gauss.bmp", gauss);
		//We store the trigonometric factors in memory to save computing time:
		double[][][] cosmsintopo = new double[N][N][2];
//		System.out.println(bragg[0] + ", " + bragg[1]);
//		System.out.println(qright[0] + ", " + qright[1]);
		double[][][] dev = new double[npts[0]][npts[1]][2];
		double[][] phase = new double[npts[0]][npts[1]];
		mag = new double[npts[0]][npts[1]];
		for (int i = 1; i < 3; i++)
		{
			//set the table values.
			for (int j = 0; j < N; j++)
				for (int k = 0; k < N; k++){
					cosmsintopo[j][k][0] = Math.cos(dot(braggTrue[i-1], j, k))*data[j][k];
					cosmsintopo[j][k][1] = -Math.sin(dot(braggTrue[i-1], j, k))*data[j][k];
				}
			SRAW.writeImage(dir + name + "wave" + i + ".bmp", cosmsintopo, false);
			getDeviationGauss(L, cosmsintopo, r1[0], dr[0], npts[0], r1[1], dr[1], npts[1], dev, gauss);
			t = printTime(" to do.", t);
			
//			mag = FieldOps.magnitude(dev);
			FieldOps.magnitude(dev, mag);
			FieldOps.phase(dev, phase);
			double max = ArrayOps.max(mag);
//			System.out.println(max + ", " + avg);
			SRAW.writeImage(dir + name + "dev" + i + ".bmp", mag, phase, new ColorScales.MYC2d(max, 0, 2*Math.PI));
			ColumnIO.writeTable(mag, dir + name + "devmag" + i + ".txt");
			ColumnIO.writeTable(phase, dir + name + "devphase" + i + ".txt");
			//For visual purposes.
			
		}
		
	}
	
	//This suppresses the third Bragg peak in the FFT and writes the modified
	//data to the disk
	public static void suppressBPandWrite(int[] bp, int k, String in, String out)
	{
		long t = System.currentTimeMillis();
		FFT2DSmall fft = FFTOps.obtainFFT(in);
		t = printTime(" to load and fft", t);
		int N = fft.falt.length;
		int x, y;
		int xm, ym;
		double[][] mag = FieldOps.magnitude(fft.fHat);
		double min = ArrayOps.min(mag), avg = ArrayOps.mean(mag);
		double value = min*100;
		double ratio;
		for (int i = -k; i <= k; i++)
			for(int j = -k; j <= k; j++)
			{
				x = util.fourier.Util.reltofHat((bp[0] + i - N/2), N);
				y = util.fourier.Util.reltofHat((bp[1] + j - N/2), N);
				xm = util.fourier.Util.reltofHat(-(bp[0] + i - N/2), N);
				ym = util.fourier.Util.reltofHat(-(bp[1] + j - N/2), N);
				//x, y are the Bragg peak. To keep the data real we must
				//kill (x, y) and (-x, -y) bragg peaks. 
				//clearly mag[x][y] = mag[-x][-y] x and y relative to 0, 0.
				ratio = mag[x][y]/value;
				fft.fHat[x][y][0] /= ratio;
				fft.fHat[x][y][1] /= ratio;
				fft.fHat[xm][ym][0] /= ratio;
				fft.fHat[xm][ym][1] /= ratio;
				mag[x][y] /= ratio;
				mag[xm][ym] /= ratio;
			}
		
//		FFTOps.writeFFTBMP(out+"fft", fft, true, true, mag);
		FFTOps.writeFFTBMPs(out+"fft", fft, true, true);
		t = printTime(" writing bmp", t);

		fft.doIFFT();
		t = printTime(" doing ifft", t);

		FieldOps.magnitude(fft.f, mag);
		ColumnIO.writeBin(mag, out + ".dat");
		t = printTime(" writing bin", t);

		SRAW.writeImage(out + ".bmp", mag);
		t = printTime(" writing last bmp", t);
	}
	
	public static double dot(double[] a, int b0, int b1)
	{
		return a[0]*b0 + a[1]*b1;
	}
	public static double dot(double[] a, double[] b)
	{
		return a[0]*b[0] + a[1]*b[1];
	}	
	
	public static long printTime(String phrase, long t)
	{
		System.out.println((System.currentTimeMillis() - t)/1000 + phrase);
		return System.currentTimeMillis();
	}
	public static long printTimeNL(String phrase, long t)
	{
		System.out.print((System.currentTimeMillis() - t)/1000 + phrase);
		return System.currentTimeMillis();
	}
	//This reads 2 lines: (x1, y1)
	//						(x2, y2)
	public static int[][] readBragg(String file)
	{
		Scanner s = null;
		try {
			s = new Scanner(new File(file));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		String line1 = s.nextLine();
		String line2 = s.nextLine();
		int[] b1 = {Integer.parseInt(line1.substring(0, line1.indexOf(","))), Integer.parseInt(line1.substring(line1.indexOf(",") + 2))};
		int[] b2 = {Integer.parseInt(line2.substring(0, line2.indexOf(","))), Integer.parseInt(line2.substring(line2.indexOf(",") + 2))};
		return new int[][] {b1, b2};
	}
	public static int[][] readBragg3(String file)
	{
		Scanner s = null;
		try {
			s = new Scanner(new File(file));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		String line1 = s.nextLine();
		String line2 = s.nextLine();
		String line3 = s.nextLine();
		int[] b1 = {Integer.parseInt(line1.substring(0, line1.indexOf(","))), Integer.parseInt(line1.substring(line1.indexOf(",") + 2))};
		int[] b2 = {Integer.parseInt(line2.substring(0, line2.indexOf(","))), Integer.parseInt(line2.substring(line2.indexOf(",") + 2))};
		int[] b3 = {Integer.parseInt(line3.substring(0, line3.indexOf(","))), Integer.parseInt(line3.substring(line3.indexOf(",") + 2))};
		return new int[][] {b1, b2, b3};
	}
	public static double[][] readBraggTrue(String file, boolean do3bragg)
	{
		Scanner s = null;
		try {
			s = new Scanner(new File(file));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		if (!do3bragg){
			s.nextLine();
			s.nextLine();
			s.nextLine();
			s.nextLine();
//			s.nextLine();
			String line1 = s.nextLine();
			String line2 = s.nextLine();
			double[] b1 = {Double.parseDouble(line1.substring(0, line1.indexOf(","))), Double.parseDouble(line1.substring(line1.indexOf(",") + 2))};
			double[] b2 = {Double.parseDouble(line2.substring(0, line2.indexOf(","))), Double.parseDouble(line2.substring(line2.indexOf(",") + 2))};
			return new double[][] {b1, b2};
		}
		else{
			s.nextLine();
			s.nextLine();
			s.nextLine();
			s.nextLine();
			s.nextLine();
			String line1 = s.nextLine();
			String line2 = s.nextLine();
			String line3 = s.nextLine();
			double[] b1 = {Double.parseDouble(line1.substring(0, line1.indexOf(","))), Double.parseDouble(line1.substring(line1.indexOf(",") + 2))};
			double[] b2 = {Double.parseDouble(line2.substring(0, line2.indexOf(","))), Double.parseDouble(line2.substring(line2.indexOf(",") + 2))};
			double[] b3 = {Double.parseDouble(line3.substring(0, line3.indexOf(","))), Double.parseDouble(line3.substring(line3.indexOf(",") + 2))};
			return new double[][] {b1, b2, b3};
		}
	}
	public static void writeBraggTrue(String file, double[][] bragg)
	{
		Scanner s = null;
		try {
			s = new Scanner(new File(file));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		String contents = "";
		while (s.hasNextLine())
			contents += s.nextLine() + "\r\n";
		contents += "\r\n";
		for (int i = 0; i < bragg.length; i++)
			contents += bragg[i][0] + ", " + bragg[i][1] + "\r\n";
		ColumnIO.writeString(contents, file);
	}
	//returns the third line of Bragg.
	public static int readBraggSmear(String file)
	{
		Scanner s = null;
		try {
			s = new Scanner(new File(file));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		s.nextLine();
		s.nextLine();
		return Integer.parseInt(s.nextLine());
	}
	public static int readBraggSmear3(String file)
	{
		Scanner s = null;
		try {
			s = new Scanner(new File(file));
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		s.nextLine();
		s.nextLine();
		s.nextLine();
		return Integer.parseInt(s.nextLine());
	}
	public static void applyTranslationEach(String dir, String topo, String phasex, String phasey, String out)
	{
		double[][] data = ColumnIO.readSquareTable(dir + topo);
		double[][] px = ColumnIO.readSquareTable(dir + phasex);
		double[][] py = ColumnIO.readSquareTable(dir + phasey);
		double[][] bragg = readBraggTrue(dir + "bragg.txt", false);
		double[][][] u = FieldOps.getU(px, py, bragg, 1);
		FieldOps.subtractAvg(u); //to make the "average translation" zero
		ColumnIO.writeBin(u, 0, dir + out + "ux.dat");
		ColumnIO.writeBin(u, 1, dir + out + "uy.dat");
		
		SRAW.writeImage(dir + out + "u.bmp", u, false);
		SRAW.writeImage(dir + out + "ux.bmp", u, 0);
		SRAW.writeImage(dir + out + "uy.bmp", u, 1);
		
//		double[][] answer = FieldOps.translateMap(data, u);
		while (u.length < data.length)
			u = FieldOps.expand(u);
		double[][] answer = FieldOps.applyUField(data, u, 8, 8);
		
		ColumnIO.writeBin(answer, dir + out + ".dat");
		SRAW.writeImage(dir + out + ".bmp", answer);
		FFT2DSmall f = new FFT2DSmall(answer); f.doFFT();
		FFTOps.writeFFTBMPCent(dir + out + "fft", f, true, true);
	
	}
	public static void finishUField(String dir, String topo, double[] L, String out, boolean did3bragg)
	{
		double[][] data = ColumnIO.readSquareTable(dir + topo + ".dat");
		String[] phase, mark = {"x", "y", "z"};
		phase = new String [L.length];
		for (int i = 0; i < L.length; i++)
			phase[i] = "devoutbinft\\" + topo +  mark[i] + "L" + L[i] + "cont.dat";
		String phasex = "", phasey = "";
		double[][] bragguse = new double[2][];
		double[][] bragg = readBraggTrue(dir + "bragg.txt", did3bragg);
		if (L.length == 3)
		{
			if (L[0] == 0)
			{
				phasex = phase[1]; bragguse[0] = bragg[1];
				phasey = phase[2]; bragguse[1] = bragg[2];
			}
			else if (L[1] == 0)
			{
				phasex = phase[0]; bragguse[0] = bragg[0];
				phasey = phase[2]; bragguse[1] = bragg[2];
			}
			else if (L[2] == 0)
			{
				phasex = phase[0]; bragguse[0] = bragg[0];
				phasey = phase[1]; bragguse[1] = bragg[1];
			}
		}
		else
		{
			phasex = phase[0];
			phasey = phase[1];
			bragguse[0] = bragg[0];
			bragguse[1] = bragg[1];
		}
		double[][] px = ColumnIO.readSquareTable(dir + phasex);
		double[][] py = ColumnIO.readSquareTable(dir + phasey);
		
		double[][][] u = FieldOps.getU(px, py, bragguse, 1);
		FieldOps.subtractAvg(u); //to make the "average translation" zero
		ColumnIO.writeBin(u, 0, dir + out + "ux.dat");
		ColumnIO.writeBin(u, 1, dir + out + "uy.dat");
		
		SRAW.writeImage(dir + out + "u.bmp", u, false);
		SRAW.writeImage(dir + out + "ux.bmp", u, 0);
		SRAW.writeImage(dir + out + "uy.bmp", u, 1);
		
//		double[][] answer = FieldOps.translateMap(data, u);
//		while (u.length < data.length)
//			u = FieldOps.expand(u);
		double[][] answer = FieldOps.applyUFieldSmooth(data, u, 8, 8);
		
		ColumnIO.writeBin(answer, dir + out + ".dat");
		SRAW.writeImage(dir + out + ".bmp", answer);
		FFT2DSmall f = new FFT2DSmall(answer); f.doFFT();
		FFTOps.writeFFTBMPCent(dir + out + "fft", f, true, true);
	
	}

	public static void applyTranslationBlock(String dir, String topo, String phasex, String phasey, String out, int blocksize)
	{
		double[][] data = ColumnIO.readSquareTable(dir + topo);
		double[][] px = ColumnIO.readSquareTable(dir + phasex);
		double[][] py = ColumnIO.readSquareTable(dir + phasey);
		double[][] bragg = readBraggTrue(dir + "bragg.txt", false);
		double[][][] u = FieldOps.getU(px, py, bragg, 1);
		FieldOps.subtractAvg(u); //to make the "average translation" zero
		ColumnIO.writeBin(u, 0, dir + out + "ux.dat");
		ColumnIO.writeBin(u, 1, dir + out + "uy.dat");
		
		SRAW.writeImage(dir + out + "u.bmp", u, false);
		SRAW.writeImage(dir + out + "ux.bmp", u, 0);
		SRAW.writeImage(dir + out + "uy.bmp", u, 1);
		
		double[][] answer = FieldOps.translateMapBlocks(data, u, blocksize);
		
		SRAW.writeImage(dir + out, answer);
		FFT2DSmall f = new FFT2DSmall(answer); f.doFFT();
		FFTOps.writeFFTBMPs(dir + out + "fft", f, true, true);
	
	}
	
	public static void makeCutoffMovie(String dir, String name, int L, int freq, double gammamin, double gammamax, int npics, String writeDir)
	{
		double[][] data = ColumnIO.readSquareTable(dir + name + ".dat");
		double[][] localAvg = FieldOps.gaussSmooth(data, L, freq);
		
		int N = data.length;
		double[][] answer = new double [N][N];
		double dgamma = (gammamax - gammamin)/(npics-1);
		double gamma;
		for (int i = 0; i < npics; i++)
		{
			gamma = gammamin + i*dgamma;
			
			FieldOps.cutOffExtremes(data, localAvg, gamma, answer);
			SRAW.writeImage(writeDir + name + MovieMaker.fromInt(i), answer);
		}
		MovieMaker.typeBMPtoAVICommand(name, 0, npics-1, 20);
	}
	public static void makeCutoffMovieFFT(String dir, String name, int L, int freq, double gammamin, double gammamax, int npics, String writeDir)
	{
		double[][] data = ColumnIO.readSquareTable(dir + name + ".dat");
		double[][] localAvg = FieldOps.gaussSmooth(data, L, freq);
		
		int N = data.length;
		double[][] answer = new double [N][N];
		FFT2DSmall fft = new FFT2DSmall(answer);
		double dgamma = (gammamax - gammamin)/(npics-1);
		double gamma;
		double[][] mag = new double[N][N];
		for (int i = 0; i < npics; i++)
		{
			gamma = gammamin + i*dgamma;
			FieldOps.cutOffExtremes(data, localAvg, gamma, answer);
			fft.setF(answer);
			fft.doFFT();
			FieldOps.magnitude(fft.fHat2(), mag);
			FieldOps.log(mag);
			SRAW.writeImage(writeDir + name +"fft" + MovieMaker.fromInt(i), mag);
		}
		MovieMaker.typeBMPtoAVICommand(name + "fft", 0, npics-1, 20);
	}
	
	
//	public static void makeTranslationMov(String dir, String topo, String phasex, String phasey, String out, double max, int npts)
//	{
//		String dir2 = "C:\\Users\\madhavanlab2011\\Documents\\BMP to AVI\\";
//		double[][] data = ColumnIO.readSquareTable(dir + topo);
//		double[][] px = ColumnIO.readSquareTable(dir + phasex);
//		double[][] py = ColumnIO.readSquareTable(dir + phasey);
//		double[][] bragg = readBraggTrue(dir + "bragg.txt");
//		
//		double[][] answer;
//		for (int i = 0; i < npts; i++){
//			answer = FieldOps.translateFromPhases(data, px, py, bragg, -max + (2*max*i)/npts);
////			ColumnIO.writeBin(answer, dir + out + ".dat");
//			SRAW.writeImage(dir2 + out + MovieMaker.fromInt(i) + ".bmp", answer);
//		}
//		MovieMaker.typeBMPtoAVICommand(out, 0, npts-1, 20);
//		
//	}
//	public static void makeTranslationFFTMov(String dir, String topo, String phasex, String phasey, String out, double max, int npts)
//	{
//		String dir2 = "C:\\Users\\madhavanlab2011\\Documents\\BMP to AVI\\";
//		double[][] data = ColumnIO.readSquareTable(dir + topo);
//		double[][] px = ColumnIO.readSquareTable(dir + phasex);
//		double[][] py = ColumnIO.readSquareTable(dir + phasey);
//		double[][] bragg = readBraggTrue(dir + "bragg.txt");
//		
//		double[][] answer;
//		double[][][] fHat2 = new double[px.length][py.length][2];
//		double[][] mag = new double[px.length][px.length];
//		FFT2DSmall fft = new FFT2DSmall(px);
//		for (int i = 0; i < npts; i++){
//			answer = FieldOps.translateFromPhases(data, px, py, bragg, -max + (			2*max*i)/npts);
////			ColumnIO.writeBin(answer, dir + out + ".dat");
//			fft.setF(answer);
//			fft.doFFT();
//			fft.fHat2(fHat2);
//			FieldOps.magnitude(fHat2, mag);
////			FieldOps.log(mag);
//			SRAW.writeImage(dir2 + out + "fftl" + MovieMaker.fromInt(i) + ".bmp", mag, false, -22);
//			
////			SRAW.writeImage(dir2 + out + MovieMaker.fromInt(i) + ".bmp", answer);
//		}
//		MovieMaker.typeBMPtoAVICommand(out + "fftl", 0, npts-1, 16);
//		
//	}
	public static void testTranslation(double[][] data, String dir, String name)
	{
		String dir2 = "C:\\Users\\madhavanlab2011\\Documents\\BMP to AVI\\";
		int N = data.length;
		double[][] data2 = new double[N][N];
		double nx = 1, ny = 0;
		double k = 20;
		for (int i = 0; i < 1000; i++)
		{
			FieldOps.translatePixelGroup(data, data2, 0, 0, N, N, (nx*i)/k, (ny*i)/k);
			SRAW.writeImage(dir2 + name + MovieMaker.fromInt(i) + ".bmp", data2);
			for (int m = 0; m < N; m++)
				for (int n = 0; n < N; n++)
					data2[m][n] = 0;
		}
		MovieMaker.typeBMPtoAVICommand(name, 0, 999, 20);
	}
	
	public static void makeRescaleXYMov(String dir, String name, double kmag, double topolength, double bragglength)
	{
		String dir2 = "C:\\Users\\madhavanlab2011\\Documents\\BMP to AVI\\";
		double klength = topolength/kmag;
		double factor = bragglength/klength;
		System.out.println(klength + "\t" + factor);
		double[][] data = ColumnIO.readSquareTable(dir + name + ".dat");
		int nframes = 1000;
		double ratio = 0.5, dratio = 1.0/nframes;
		double[][] result;
//		double[][] fftmag;
//		for (int i = 0; i < 1000; i++)
//		{
//			ratio = 0.5 + i*dratio;
//			result = FieldOps.rescaleXY(data, factor*ratio, 16);
////			ColumnIO.writeBin(result, dir + name + "_resc.dat");
//			SRAW.writeImage(dir2 + name + "_resc" + MovieMaker.fromInt(i), result);
//			fftmag = FFTOps.obtainFFTmagCent(result);
//			FieldOps.log(fftmag);
//			SRAW.writeImage(dir2 + name + "_rescfft" + MovieMaker.fromInt(i), fftmag);
//		}
		MovieMaker.typeBMPtoAVICommand(name + "_resc", 0, nframes-1, 20);
		Robo.wait(90000);
		MovieMaker.typeBMPtoAVICommand(name + "_rescfft", 0, nframes-1, 20);
	}
	public static void testExpand(double[][] data, String dir)
	{
		SRAW.writeImage(dir + "test.bmp", FieldOps.expand(data));
		SRAW.writeImage(dir + "test4.bmp", FieldOps.expand4(data));
	}
	public static void testTranslate()
	{
		double[][] source = {{1, 1}};
		double[][] target = new double [8][8];
		for (int i = 0; i < 1000; i++)
		{
			FieldOps.zero(target);
			FieldOps.translateOnePixelAlt(source, target, 0, 0, ((double)i)/120, 0);
			FieldOps.translateOnePixelAlt(source, target, 0, 1, ((double)i)/125, -1);
			for (int j = 0; j < 8; j++)
				System.out.print("\t" + target[j][0]);
//			System.out.println();
			System.out.println("\t" + ArrayOps.sum(target));
		}
	}
	public static void testAutocorrelate(String dir)
	{
		double[][] data = ColumnIO.readSquareTable(dir + "topo.dat");
		double[][] fftmag = FFTOps.obtainFFTmagCent(data);
		FieldOps.log(fftmag);
		double[][] answer = new double [data.length][data.length];
		FieldOps.autocorrelate(fftmag, answer);
		SRAW.writeImage(dir + "auto3", answer);
	}
	
	public static void doImpCorrAngleTest()
	{
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		int npts = 1000;
		ArrayList<ImpurityEmergenceEditor.Impurity> imps = ImpurityEmergenceEditor.getImpurities(FileOps.selectOpen(fc), 1);
		double[][] data = FileOps.openBin(fc);
		double[][] mask = Mask.rectMaskCent(128, 128, 64, 64);
		double[][] aroundImp = new double [128][128];
		double[][] corr = new double [imps.size()][npts];
		
		int ix, iy;
		
		for (int i = 0; i < imps.size(); i++)
		{
			ix = (int)Math.round(imps.get(i).position[0]);
			iy = (int)Math.round(imps.get(i).position[1]);
			FieldOps.expandBi(data, 8, ix-8, 16, iy-8, 16, aroundImp);
			FieldOps.changeZeroToAverage(aroundImp);
			SRAW.writeImage(fc.getCurrentDirectory().toString() + "\\impanglepics\\imp" + MovieMaker.fromInt(i), aroundImp);			
		}
	}
	public static void writeImpImagesvsEnergy()
	{
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		int npts = 1000;
		ArrayList<ImpurityEmergenceEditor.Impurity> imps = ImpurityEmergenceEditor.getImpurities(FileOps.selectOpen(fc), 1);
		AtomicCoordinatesSet latt = new AtomicCoordinatesSet(FileOps.openText());
		Topomap t = Topomap.open(fc);
		double[][][] latticeSites, latticeSitesBlowup;
		
//		double[][] data = FileOps.openBin(fc);
		double[][] mask = Mask.rectMaskCent(128, 128, 64, 64);
		double[][] aroundImp = new double [128][128];
		double[][] corr = new double [imps.size()][npts];
		
		File f;
		int[] bestLayers = new int []
				{
					395, 395, 395, 304, 258, 314, 395, 325, 271, 257,
					257, 328, 395, 395, 332, 332, 312, 312, 354, 305,
					283, 312, 304, 334, 292, 248, 318, 339, 291, 319,
					319, 376, 389, 368, 310, 269, 395, 395, 367, 391,
					340, 382, 356, 313, 340, 270, 391, 395, 374, 368,
					276, 361, 276, 380, 293, 332
				};

		BufferedImage image = new BufferedImage(t.nx, t.ny, BufferedImage.TYPE_INT_RGB);
		BufferedImage image2 = new BufferedImage(t.nx, t.ny, BufferedImage.TYPE_INT_RGB);
		double[][] pinkRing = new double [t.nx][t.ny];
		 double[] atomc1 = latt.getAtomicCoords(0, 0);
		 double[] atomc2 = latt.getAtomicCoords(0, t.nx);
		 double[] atomc3 = latt.getAtomicCoords(t.nx, 0);
		 double[] atomc4 = latt.getAtomicCoords(t.nx, t.nx);
		 double[] min = new double [2], max = new double [2];
		 min[0] = Math.min(Math.min(atomc1[0], atomc2[0]),Math.min(atomc3[0], atomc4[0]));
		 min[1] = Math.min(Math.min(atomc1[1], atomc2[1]),Math.min(atomc3[1], atomc4[1]));
		 max[0] = Math.max(Math.max(atomc1[0], atomc2[0]),Math.max(atomc3[0], atomc4[0]));
		 max[1] = Math.max(Math.max(atomc1[1], atomc2[1]),Math.max(atomc3[1], atomc4[1]));
		 
		 int xmn = FieldOps.round(min[0]);
		 int xmx = FieldOps.round(max[0])+1;
		 int ymn = FieldOps.round(min[1]);
		 int ymx = FieldOps.round(max[1])+1;
		 latticeSites = new double [xmx-xmn][ymx-ymn][2];
		 latticeSitesBlowup = new double [xmx-xmn][ymx-ymn][2];
		 for (int i = 0; i < xmx-xmn; i++)
			 for (int j = 0; j < ymx-ymn; j++)
				 latticeSites[i][j] = latt.getPixelCoords(i+xmn, j+ymn);
		 //1 pixel shift
		 for (int i = 0; i < xmx-xmn; i++)
			 for (int j = 0; j < ymx-ymn; j++)
				 latticeSites[i][j][0] += 1;
		 
		 
		int ix, iy;
		int zoom = 4;
		int initsize = 32;
		double ringrsmall = 10; double ringtsmall = 0.5;
		BufferedImage blowup = new BufferedImage(zoom*initsize, zoom*initsize, BufferedImage.TYPE_INT_RGB);
		BufferedImage blowup2 = new BufferedImage(zoom*initsize, zoom*initsize, BufferedImage.TYPE_INT_RGB);
		double[][] dots = new double [zoom*initsize][zoom*initsize];
		
		for (int i = 0; i < imps.size(); i++)
		{
			f = new File(fc.getCurrentDirectory().toString() + "\\images\\imp" + MovieMaker.fromInt(i) + "\\");
			if (!f.exists()) f.mkdir();
			ix = (int)Math.round(imps.get(i).position[0]);
			iy = (int)Math.round(imps.get(i).position[1]);
			 for (int p = 0; p < xmx-xmn; p++)	//determine the position of the lattice sites
				 for (int q = 0; q < ymx-ymn; q++){
					 latticeSitesBlowup[p][q][0] = latticeSites[p][q][0] - (ix-initsize/2);							 
					 latticeSitesBlowup[p][q][1] = latticeSites[p][q][1] - (iy-initsize/2);
					 latticeSitesBlowup[p][q][0] *= zoom;
					 latticeSitesBlowup[p][q][1] *= zoom;
				 }
			FieldOps.putGaussDots(dots, latticeSitesBlowup, 1, 1, true);
			for (int k = 0; k < t.nlayers; k++)
			{
				FieldOps.expandBi(t.data[k], zoom, ix-initsize/2, initsize, iy-initsize/2, initsize, aroundImp);
				FieldOps.changeZeroToAverage(aroundImp);
//				SRAW.writeImage(f.toString() + "\\imp" + MovieMaker.fromInt(k), aroundImp);
				SRAW.writeImage(blowup, aroundImp, null);
				
				
				ImageEditing.moveEachPixelTowardAColor(blowup, dots, Color.BLUE, blowup2);
				SRAW.writeImage(f.toString() + "\\imp" + MovieMaker.fromInt(k), blowup2);
				if (k == bestLayers[i])
				{
					double rangemax = ringrsmall + 4*ringtsmall;
					double g, r;
					for (int m = 0; m < t.nx; m++){
						for (int p = 0; p < t.ny; p++){ pinkRing[m][p] = 0;
							if (Math.abs(m - imps.get(i).position[0]) < rangemax && Math.abs(p - imps.get(i).position[1]) < rangemax)
							{
								r = Math.sqrt(Math.pow(m - imps.get(i).position[0], 2) + Math.pow(p - imps.get(i).position[1], 2));
								g = Math.exp(-(r-ringrsmall)*(r-ringrsmall)/(ringtsmall*ringtsmall));
								pinkRing[m][p] += g;
							}
							pinkRing[m][p] = Math.min(pinkRing[m][p], 1);
						}
					}
					SRAW.writeImage(image, t.data[k], null);
					ImageEditing.moveEachPixelTowardAColor(image, pinkRing, Color.MAGENTA, image2);
					SRAW.writeImage(f.toString() + "\\circle" + k, image2);

				}
					
			}
		}
	}
	public static void writeImpImagesOneLayer(Layer t, ArrayList<ImpurityEmergenceEditor.Impurity> imps, AtomicCoordinatesSet latt, String dir, boolean drawDots)
	{
//		JFileChooser fc = new JFileChooser(Topomap.stddir);
		int npts = 1000;
//		ArrayList<ImpurityEmergenceEditor.Impurity> imps = ImpurityEmergenceEditor.getImpurities(FileOps.selectOpen(fc), 1);
//		AtomicCoordinatesSet latt = new AtomicCoordinatesSet(FileOps.openText());
//		Topomap t = Topomap.open(fc);
		double[][][] latticeSites, latticeSitesBlowup;
		
//		double[][] data = FileOps.openBin(fc);
//		double[][] mask = Mask.rectMaskCent(128, 128, 64, 64);
		
		File f;

		BufferedImage image = new BufferedImage(t.nx, t.ny, BufferedImage.TYPE_INT_RGB);
		BufferedImage image2 = new BufferedImage(t.nx, t.ny, BufferedImage.TYPE_INT_RGB);
		double[][] pinkRing = new double [t.nx][t.ny];
		 double[] atomc1 = latt.getAtomicCoords(0, 0);
		 double[] atomc2 = latt.getAtomicCoords(0, t.nx);
		 double[] atomc3 = latt.getAtomicCoords(t.nx, 0);
		 double[] atomc4 = latt.getAtomicCoords(t.nx, t.nx);
		 double[] min = new double [2], max = new double [2];
		 min[0] = Math.min(Math.min(atomc1[0], atomc2[0]),Math.min(atomc3[0], atomc4[0]));
		 min[1] = Math.min(Math.min(atomc1[1], atomc2[1]),Math.min(atomc3[1], atomc4[1]));
		 max[0] = Math.max(Math.max(atomc1[0], atomc2[0]),Math.max(atomc3[0], atomc4[0]));
		 max[1] = Math.max(Math.max(atomc1[1], atomc2[1]),Math.max(atomc3[1], atomc4[1]));
		 
		 int xmn = FieldOps.round(min[0]);
		 int xmx = FieldOps.round(max[0])+1;
		 int ymn = FieldOps.round(min[1]);
		 int ymx = FieldOps.round(max[1])+1;
		 latticeSites = new double [xmx-xmn][ymx-ymn][2];
		 latticeSitesBlowup = new double [xmx-xmn][ymx-ymn][2];
		 for (int i = 0; i < xmx-xmn; i++)
			 for (int j = 0; j < ymx-ymn; j++)
				 latticeSites[i][j] = latt.getPixelCoords(i+xmn, j+ymn);
//		 //1 pixel shift
//		 for (int i = 0; i < xmx-xmn; i++)
//			 for (int j = 0; j < ymx-ymn; j++)
//				 latticeSites[i][j][0] += 1;
		 
		 
		int ix, iy;
		int zoom = 4;
		int initsize = 32;
		double ringrsmall = 10; double ringtsmall = 0.5;
		BufferedImage blowup = new BufferedImage(zoom*initsize, zoom*initsize, BufferedImage.TYPE_INT_RGB);
		BufferedImage blowup2 = new BufferedImage(zoom*initsize, zoom*initsize, BufferedImage.TYPE_INT_RGB);
		double[][] dots = new double [zoom*initsize][zoom*initsize];
		double[][] aroundImp = new double [zoom*initsize][zoom*initsize];
		f = new File(dir + "\\images\\imps\\");
		if (!f.exists()) f.mkdirs();

		for (int i = 0; i < imps.size(); i++)
		{
			ix = (int)Math.round(imps.get(i).position[0]);
			iy = (int)Math.round(imps.get(i).position[1]);
			 for (int p = 0; p < xmx-xmn; p++)	//determine the position of the lattice sites
				 for (int q = 0; q < ymx-ymn; q++){
					 latticeSitesBlowup[p][q][0] = latticeSites[p][q][0] - (ix-initsize/2);							 
					 latticeSitesBlowup[p][q][1] = latticeSites[p][q][1] - (iy-initsize/2);
					 latticeSitesBlowup[p][q][0] *= zoom;
					 latticeSitesBlowup[p][q][1] *= zoom;
				 }
			if (drawDots)	FieldOps.putGaussDots(dots, latticeSitesBlowup, 1, 1, true);
			FieldOps.expandBi(t.data, zoom, ix-initsize/2, initsize, iy-initsize/2, initsize, aroundImp);
			FieldOps.changeZeroToAverage(aroundImp);
//				SRAW.writeImage(f.toString() + "\\imp" + MovieMaker.fromInt(k), aroundImp);
			SRAW.writeImage(blowup, aroundImp, null);
			if (drawDots)
			{
				ImageEditing.moveEachPixelTowardAColor(blowup, dots, Color.BLUE, blowup2);
				SRAW.writeImage(f.toString() + "\\imp" + MovieMaker.fromInt(i), blowup2);
			}
			else
				SRAW.writeImage(f.toString() + "\\imp" + MovieMaker.fromInt(i), blowup);
				

			double rangemax = ringrsmall + 4*ringtsmall;
			double g, r;
			for (int m = 0; m < t.nx; m++){
				for (int p = 0; p < t.ny; p++){ pinkRing[m][p] = 0;
					if (Math.abs(m - imps.get(i).position[0]) < rangemax && Math.abs(p - imps.get(i).position[1]) < rangemax)
					{
						r = Math.sqrt(Math.pow(m - imps.get(i).position[0], 2) + Math.pow(p - imps.get(i).position[1], 2));
						g = Math.exp(-(r-ringrsmall)*(r-ringrsmall)/(ringtsmall*ringtsmall));
						pinkRing[m][p] += g;
					}
					pinkRing[m][p] = Math.min(pinkRing[m][p], 1);
				}
			}
			SRAW.writeImage(image, t.data, null);
			ImageEditing.moveEachPixelTowardAColor(image, pinkRing, Color.MAGENTA, image2);
			SRAW.writeImage(f.toString() + "\\circle" + i, image2);
		}
	}
	public static void writeImpImagesOneLayerWithRoot2Splitting(Layer t, ArrayList<drawing.GaussSquareImpurityAngleSuite.Impurity> imps, AtomicCoordinatesSet latt, String dir, boolean drawDots)
	{
//		JFileChooser fc = new JFileChooser(Topomap.stddir);
//		ArrayList<ImpurityEmergenceEditor.Impurity> imps = ImpurityEmergenceEditor.getImpurities(FileOps.selectOpen(fc), 1);
//		AtomicCoordinatesSet latt = new AtomicCoordinatesSet(FileOps.openText());
//		Topomap t = Topomap.open(fc);
		double[][][] latticeSites, latticeSitesBlowup;

		for (int i = 0; i < imps.size(); i++)
			if (imps.get(i).latticePos[0] == 0) imps.get(i).latticePos = latt.getAtomicCoords(imps.get(i).position);

		double x, y;
		int xm, ym;

//		double[][] data = FileOps.openBin(fc);
//		double[][] mask = Mask.rectMaskCent(128, 128, 64, 64);
		
		File[] f;

		BufferedImage image = new BufferedImage(t.nx, t.ny, BufferedImage.TYPE_INT_RGB);
		BufferedImage image2 = new BufferedImage(t.nx, t.ny, BufferedImage.TYPE_INT_RGB);
		double[][] pinkRing = new double [t.nx][t.ny];
		 double[] atomc1 = latt.getAtomicCoords(0, 0);
		 double[] atomc2 = latt.getAtomicCoords(0, t.nx);
		 double[] atomc3 = latt.getAtomicCoords(t.nx, 0);
		 double[] atomc4 = latt.getAtomicCoords(t.nx, t.nx);
		 double[] min = new double [2], max = new double [2];
		 min[0] = Math.min(Math.min(atomc1[0], atomc2[0]),Math.min(atomc3[0], atomc4[0]));
		 min[1] = Math.min(Math.min(atomc1[1], atomc2[1]),Math.min(atomc3[1], atomc4[1]));
		 max[0] = Math.max(Math.max(atomc1[0], atomc2[0]),Math.max(atomc3[0], atomc4[0]));
		 max[1] = Math.max(Math.max(atomc1[1], atomc2[1]),Math.max(atomc3[1], atomc4[1]));
		 
		 int xmn = FieldOps.round(min[0]);
		 int xmx = FieldOps.round(max[0])+1;
		 int ymn = FieldOps.round(min[1]);
		 int ymx = FieldOps.round(max[1])+1;
		 latticeSites = new double [xmx-xmn][ymx-ymn][2];
		 latticeSitesBlowup = new double [xmx-xmn][ymx-ymn][2];
		 for (int i = 0; i < xmx-xmn; i++)
			 for (int j = 0; j < ymx-ymn; j++)
				 latticeSites[i][j] = latt.getPixelCoords(i+xmn, j+ymn);
//		 //1 pixel shift
//		 for (int i = 0; i < xmx-xmn; i++)
//			 for (int j = 0; j < ymx-ymn; j++)
//				 latticeSites[i][j][0] += 1;

		int[] subA = new int[imps.size()];
		String[] output = new String[] {"A", "B"};

		 
		int ix, iy;
		int zoom = 4;
		int initsize = 64;
		double ringrsmall = 15; double ringtsmall = 0.5;
		BufferedImage blowup = new BufferedImage(zoom*initsize, zoom*initsize, BufferedImage.TYPE_INT_RGB);
		BufferedImage blowup2 = new BufferedImage(zoom*initsize, zoom*initsize, BufferedImage.TYPE_INT_RGB);
		double[][] dots = new double [zoom*initsize][zoom*initsize];
		double[][] aroundImp = new double [zoom*initsize][zoom*initsize];
		f = new File[] {new File(dir + "\\images\\imps\\" + output[0] + "\\"), new File(dir + "\\images\\imps\\" + output[1] + "\\")};
		if (!f[0].exists()) f[0].mkdirs();
		if (!f[1].exists()) f[1].mkdirs();

		String outdir;
		for (int i = 0; i < imps.size(); i++)
		{
			x = imps.get(i).latticePos[0];
			y = imps.get(i).latticePos[1];
			xm = FieldOps.roundDown(x);
			ym = FieldOps.roundDown(y);
			subA[i] = (xm+ym + 100000)%2;

			
			ix = (int)Math.round(imps.get(i).position[0]);
			iy = (int)Math.round(imps.get(i).position[1]);
			 for (int p = 0; p < xmx-xmn; p++)	//determine the position of the lattice sites
				 for (int q = 0; q < ymx-ymn; q++){
					 latticeSitesBlowup[p][q][0] = latticeSites[p][q][0] - (ix-initsize/2);							 
					 latticeSitesBlowup[p][q][1] = latticeSites[p][q][1] - (iy-initsize/2);
					 latticeSitesBlowup[p][q][0] *= zoom;
					 latticeSitesBlowup[p][q][1] *= zoom;
				 }
			if (drawDots)	FieldOps.putGaussDots(dots, latticeSitesBlowup, 1, 1, true);
			FieldOps.expandBi(t.data, zoom, ix-initsize/2, initsize, iy-initsize/2, initsize, aroundImp);
			FieldOps.changeZeroToAverage(aroundImp);
//				SRAW.writeImage(f.toString() + "\\imp" + MovieMaker.fromInt(k), aroundImp);
			SRAW.writeImage(blowup, aroundImp, 5);
			if (drawDots)
			{
				ImageEditing.moveEachPixelTowardAColor(blowup, dots, Color.MAGENTA, blowup2);
				SRAW.writeImage(f[subA[i]].toString() + "\\imp" + MovieMaker.fromInt(i), blowup2);
			}
			else
				SRAW.writeImage(f[subA[i]].toString() + "\\imp" + MovieMaker.fromInt(i), blowup);
				

			double rangemax = ringrsmall + 4*ringtsmall;
			double g, r;
			for (int m = 0; m < t.nx; m++){
				for (int p = 0; p < t.ny; p++){ pinkRing[m][p] = 0;
					if (Math.abs(m - imps.get(i).position[0]) < rangemax && Math.abs(p - imps.get(i).position[1]) < rangemax)
					{
						r = Math.sqrt(Math.pow(m - imps.get(i).position[0], 2) + Math.pow(p - imps.get(i).position[1], 2));
						g = Math.exp(-(r-ringrsmall)*(r-ringrsmall)/(ringtsmall*ringtsmall));
						pinkRing[m][p] += g;
					}
					pinkRing[m][p] = Math.min(pinkRing[m][p], 1);
				}
			}
			SRAW.writeImage(image, t.data, 5);
			ImageEditing.moveEachPixelTowardAColor(image, pinkRing, Color.MAGENTA, image2);
			SRAW.writeImage(f[subA[i]].toString() + "\\circle" + i, image2);
		}
	}
	public static void writeImpImagesvsEnergyRotated()
	{
		Topomap.setStdDir();
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		int npts = 1000;
		ArrayList<ImpurityEmergenceEditor.Impurity> imps = ImpurityEmergenceEditor.getImpurities(FileOps.selectOpen(fc), 1);
		AtomicCoordinatesSet latt = new AtomicCoordinatesSet(FileOps.openText(fc));
		
		double[][] matrix = new double [2][2];
		double[][] matTranspose = new double [2][2];
		double theta = Math.toRadians(24.5671713206013148);
		Matrix.putRotationMatrix(theta, matrix);
		Matrix.putRotationMatrix(-theta, matTranspose);
		AtomicCoordinatesSet lattRot = latt.getRotatedCopy(matrix, latt.getOrigin());
		
		Topomap t = Topomap.open(fc);
		double[][][] latticeSites, latticeSitesBlowup, latticeSites2;
		
//		double[][] data = FileOps.openBin(fc);
//		double[][] mask = Mask.rectMaskCent(128, 128, 64, 64);
		int zoom = 6;
		int initsize = 20;
		double[][] aroundImp = new double [zoom*initsize][zoom*initsize];
//		double[][] corr = new double [imps.size()][npts];
		
		File f;
		int[] bestLayers = new int []
				{
					395, 395, 395, 304, 258, 314, 395, 325, 271, 257,
					257, 328, 395, 395, 332, 332, 312, 312, 354, 305,
					283, 312, 304, 334, 292, 248, 318, 339, 291, 319,
					319, 376, 389, 368, 310, 269, 395, 395, 367, 391,
					340, 382, 356, 313, 340, 270, 391, 395, 374, 368,
					276, 361, 276, 380, 293, 332
				};

		BufferedImage image = new BufferedImage(t.nx, t.ny, BufferedImage.TYPE_INT_RGB);
		BufferedImage image2 = new BufferedImage(t.nx, t.ny, BufferedImage.TYPE_INT_RGB);
		double[][] pinkRing = new double [t.nx][t.ny];
		 double[] atomc1 = latt.getAtomicCoords(-t.nx, -t.ny);
		 double[] atomc2 = latt.getAtomicCoords(-t.nx, t.ny);
		 double[] atomc3 = latt.getAtomicCoords(t.nx, -t.ny);
		 double[] atomc4 = latt.getAtomicCoords(t.nx, t.ny);
		 double[] min = new double [2], max = new double [2];
		 min[0] = Math.min(Math.min(atomc1[0], atomc2[0]),Math.min(atomc3[0], atomc4[0]));
		 min[1] = Math.min(Math.min(atomc1[1], atomc2[1]),Math.min(atomc3[1], atomc4[1]));
		 max[0] = Math.max(Math.max(atomc1[0], atomc2[0]),Math.max(atomc3[0], atomc4[0]));
		 max[1] = Math.max(Math.max(atomc1[1], atomc2[1]),Math.max(atomc3[1], atomc4[1]));
		 
		 int xmn = FieldOps.round(min[0]);
		 int xmx = FieldOps.round(max[0])+1;
		 int ymn = FieldOps.round(min[1]);
		 int ymx = FieldOps.round(max[1])+1;
		 latticeSites = new double [xmx-xmn][ymx-ymn][2];
		 latticeSites2 = new double [xmx-xmn][ymx-ymn][2];
		 latticeSitesBlowup = new double [xmx-xmn][ymx-ymn][2];
		 for (int i = 0; i < xmx-xmn; i++)
			 for (int j = 0; j < ymx-ymn; j++)
				 latticeSites[i][j] = latt.getPixelCoords(i+xmn, j+ymn);
		 //1 pixel shift
//		 for (int i = 0; i < xmx-xmn; i++)
//			 for (int j = 0; j < ymx-ymn; j++)
//				 latticeSites[i][j][0] += 1;
		 
		 
		int ix, iy;
		double[] io = new double [2]; //the corner of the image
		double[] ioRot = new double [2];
		double ringrsmall = 10; double ringtsmall = 0.5;
		BufferedImage blowup = new BufferedImage(zoom*initsize, zoom*initsize, BufferedImage.TYPE_INT_RGB);
		BufferedImage blowup2 = new BufferedImage(zoom*initsize, zoom*initsize, BufferedImage.TYPE_INT_RGB);
		double[][] dots = new double [zoom*initsize][zoom*initsize];
		double[][][] dots2 = new double [4][t.nx][t.ny];
		boolean[][] isInPq = new boolean[xmx-xmn][ymx-ymn];
		
		for (int i = 0; i < imps.size(); i++)
		{
//			if (i != 29) continue;// && i != 31) continue;
			f = new File(fc.getCurrentDirectory().toString() + "\\images\\imp" + MovieMaker.fromInt(i) + "\\");
			if (!f.exists()) f.mkdir();
			ix = (int)Math.round(imps.get(i).position[0]);
			iy = (int)Math.round(imps.get(i).position[1]);
			io = new double [] {ix-initsize/2 + 1/(2.0*zoom), iy-initsize/2 + 1/(2.0*zoom)};
			Matrix.putProductWithOrigin(matTranspose, io, ioRot, imps.get(i).position);
			for (int p = 0; p < xmx-xmn; p++)	//determine the position of the lattice sites
				 for (int q = 0; q < ymx-ymn; q++){
					 Matrix.putProductWithOrigin(matrix, latticeSites[p][q], latticeSites2[p][q], imps.get(i).position);
					 latticeSitesBlowup[p][q][0] = latticeSites2[p][q][0] - io[0];							 
					 latticeSitesBlowup[p][q][1] = latticeSites2[p][q][1] - io[1];
//					 latticeSitesBlowup[p][q][0] = latticeSites2[p][q][0] - ioRot[0];							 
//					 latticeSitesBlowup[p][q][1] = latticeSites2[p][q][1] - ioRot[1];
					 latticeSitesBlowup[p][q][0] *= zoom;
					 latticeSitesBlowup[p][q][1] *= zoom;
					 if (latticeSitesBlowup[p][q][0] > 0 && latticeSitesBlowup[p][q][0] < 128 && latticeSitesBlowup[p][q][1] > 0 && latticeSitesBlowup[p][q][1] < 128)
					 {
						 isInPq[p][q] = true;
//						 System.out.println(Printer.vectorP(latticeSites[p][q]) + "\t" + Printer.vectorP(latticeSites2[p][q]) + "\t" + Printer.vectorP(latticeSitesBlowup[p][q]));
//						 System.out.println(Printer.vectorP(imps.get(i).position));
//						 System.out.println(Distance.distance(latticeSites[p][q], imps.get(i).position));
//						 System.out.println(Distance.distance(latticeSites2[p][q], imps.get(i).position));
//						 System.out.println(Distance.distance(latticeSitesBlowup[p][q][0] - (imps.get(i).position[0] - io[0]), latticeSitesBlowup[p][q][1] - (imps.get(i).position[1] - io[1])));
					 }
					 else
						 isInPq[p][q] = false;
				 }
//			FieldOps.putGaussDots(dots2[0], latticeSites, 0.5, isInPq);
//			FieldOps.putGaussDots(dots2[3], latticeSites2, 0.5, isInPq);
//			FieldOps.putGaussDots(dots2[1], new double[][] {io}, 0.5);
//			FieldOps.putGaussDots(dots2[2], new double[][] {ioRot}, 0.5);
//			FieldOps.putGaussDots(dots2[3], new double[][] {imps.get(i).position}, 0.5);
			
			FieldOps.putGaussDots(dots, latticeSitesBlowup, 1, isInPq);
			for (int k = 0; k < t.nlayers; k++)
			{
				FieldOps.expandBiRotated(t.data[k], -theta, imps.get(i).position, zoom, ix-initsize/2, initsize, iy-initsize/2, initsize, aroundImp);
				FieldOps.changeZeroToAverage(aroundImp);
//				SRAW.writeImage(f.toString() + "\\imp" + MovieMaker.fromInt(k), aroundImp);
				SRAW.writeImage(blowup, aroundImp, 9);
//				ColumnIO.writeBin(aroundImp, f.toString() + "\\imp" + MovieMaker.fromInt(k) + ".dat");
				
				
				ImageEditing.moveEachPixelTowardAColor(blowup, dots, Color.MAGENTA, blowup2);
				SRAW.writeImage(f.toString() + "\\imp" + MovieMaker.fromInt(k), blowup2);
				if (k == /*bestLayers[i]*/ 250)
				{
					double rangemax = ringrsmall + 4*ringtsmall;
					double g, r;
					for (int m = 0; m < t.nx; m++){
						for (int p = 0; p < t.ny; p++){ pinkRing[m][p] = 0;
							if (Math.abs(m - imps.get(i).position[0]) < rangemax && Math.abs(p - imps.get(i).position[1]) < rangemax)
							{
								r = Math.sqrt(Math.pow(m - imps.get(i).position[0], 2) + Math.pow(p - imps.get(i).position[1], 2));
								g = Math.exp(-(r-ringrsmall)*(r-ringrsmall)/(ringtsmall*ringtsmall));
								pinkRing[m][p] += g;
							}
							pinkRing[m][p] = Math.min(pinkRing[m][p], 1);
						}
					}
					SRAW.writeImage(image, t.data[k], null);
					ImageEditing.moveEachPixelTowardAColor(image, pinkRing, Color.MAGENTA, image2);
//					ImageEditing.moveEachPixelTowardAColor(image2, dots2[0], Color.BLUE, image2);
//					ImageEditing.moveEachPixelTowardAColor(image2, dots2[1], Color.CYAN, image2);
//					ImageEditing.moveEachPixelTowardAColor(image2, dots2[2], Color.GREEN, image2);
//					ImageEditing.moveEachPixelTowardAColor(image2, dots2[3], Color.ORANGE, image2);
					SRAW.writeImage(f.toString() + "\\circle" + k, image2);

				}
					
			}
		}
	}
	public static void getDifferences()
	{
		double[][] temp = new double [128][128];
		double[][] m1 = new double [128][128];
		double[][] m2 = new double [128][128];
		double[][] mat = new double [2][2];
		int[] c = new int[] {64, 64};
//		Mask.putRectMaskRotated(128, 128, 65, 65, c, 0, temp, m1, mat);
//		Mask.putRectMaskRotated(128, 128, 65, 65, c, Math.PI/2, temp, m2, mat);
		double[][] minus = FieldOps.minus(m1, m2);
		FileOps.writeTableASCII(minus);
		
	}
}

