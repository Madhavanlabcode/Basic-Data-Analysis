package util.fourier;

import image.ImageEditing;

import java.awt.image.BufferedImage;
import java.io.File;

import edu.emory.mathcs.jtransforms.fft.DoubleFFT_1D;

import util.ArrayOps;
import util.Complex;
import util.FieldOps;
import util.color.ColorScales;
import util.fileops.ColumnIO;
import main.SRAW;

public class FFTOps {

	public static void writeFFTBMP(String in, String out, boolean justmag, boolean logmag)
	{
		long t = System.currentTimeMillis(), tf;
		
		System.out.println("start");
		double[][] data = SRAW.getData(in);
		System.out.println((System.currentTimeMillis() - t)/1000 + " to load.");
		t = System.currentTimeMillis();
		
		FFT2D fft = new FFT2D(data);
		System.out.println((System.currentTimeMillis() - t)/1000 + " making fft.");
		t = System.currentTimeMillis();
		
		fft.doFFT();
		System.out.println((System.currentTimeMillis() - t)/1000 + " doing fft.");
		t = System.currentTimeMillis();
		
		double[][] mag = FieldOps.magnitude(fft.fHat2);

		if (logmag)
		FieldOps.log(mag);
		
		if (justmag)
		{
			SRAW.writeImage(out, data);
			return;
		}
	}

	public static double[][] loadSquareBin(String in)
	{
		return ColumnIO.readSquareTable(in);
	}
	public static FFT2DSmall obtainFFT(String in)
	{
		double[][] data = loadSquareBin(in);
		FFT2DSmall fft = new FFT2DSmall(data);
		fft.doFFT();
		return fft;
	}
	public static FFT2DSmall obtainFFT(double[][] data)
	{
		FFT2DSmall fft = new FFT2DSmall(data);
		fft.doFFT();
		return fft;
	}
	public static FFT2DSmall obtainFFT(double[][][] data)
	{
		FFT2DSmall fft = new FFT2DSmall(data);
		fft.doFFT();
		return fft;
	}
	
	public static void writeBMPDir(String dir, boolean log, String keyword)
	{
		File[] files = new File(dir).listFiles();
		FFT2DSmall fft;
		double[][] x;
		for (int k = 0; k < files.length; k++)
		{
			if (files[k].toString().endsWith(".dat") && files[k].toString().contains(keyword))
			{
				fft = obtainFFT(ColumnIO.readSquareTable(files[k].toString()));
				FFTOps.writeFFTBMPCent(files[k].toString().substring(0, files[k].toString().length() - 4) + "fft.bmp", fft, true, true);
				//				x = new double [fft.N][fft.N];
//				for (int i = 0; i < fft.N; i++)
//					for (int j = 0; j < fft.N; j++)
//					{
//						x[i][j] = Math.sqrt(Math.pow(fft.getFHat2Re(i, j), 2) + Math.pow(fft.getFHat2Im(i, j), 2));
//					}
//				if (log)
//					FieldOps.log(x);
//				SRAW.writeImage(files[k].toString().substring(0, files[k].toString().length() - 4) + "fft.bmp", x);
			}
		}
	}
	public static void writeBMPComplexDir(String dir, boolean log, String keyword)
	{
		File[] files = new File(dir).listFiles();
		FFT2DSmall fft = null;
		boolean found = false;
		for (int k = 0; k < files.length; k++)
		{
			found = false;
			if (files[k].toString().endsWith(".dat") && files[k].toString().contains(keyword) && files[k].toString().contains("mag"))
			{
				for (int i = k; i < files.length && !found; i++)
				{
					if (files[i].toString().endsWith(".dat") && files[i].toString().contains(keyword) && files[i].toString().contains("phase")){
						fft = obtainFFT(ColumnIO.readSquareTables(files[k].toString(), files[i].toString(), true));
						found = true;
					}
				}
				FFTOps.writeFFTBMPCent(files[k].toString().substring(0, files[k].toString().length() - 7), fft, true, false);
			}
		}
	}
	public static void writeBMPDirScale(String dir, boolean log, double min, double max)
	{
		File[] files = new File(dir).listFiles();
		FFT2DSmall fft;
		double[][] x;
		for (int k = 0; k < files.length; k++)
		{
			if (files[k].toString().endsWith(".dat"))
			{
				fft = obtainFFT(ColumnIO.readSquareTable(files[k].toString()));
				x = new double [fft.N][fft.M];
				for (int i = 0; i < fft.N; i++)
					for (int j = 0; j < fft.M; j++)
					{
						x[i][j] = Math.sqrt(Math.pow(fft.getFHat2Re(i, j), 2) + Math.pow(fft.getFHat2Im(i, j), 2));
					}
				if (log)
					FieldOps.log(x);
				SRAW.writeImage(files[k].toString().substring(0, files[k].toString().length() - 4) + "fft_s.bmp", x, min, max);
			}
		}
	}
	public static void writeBMPDirScaleMiddle(String dir, boolean log, double min, double max, String keyword, int factor)
	{
		File[] files = new File(dir).listFiles();
		FFT2DSmall fft;
		double[][] x;
		int dx, dy, minx, miny;
		for (int k = 0; k < files.length; k++)
		{
			if (files[k].toString().endsWith(".dat") && files[k].toString().contains(keyword))
			{
				fft = obtainFFT(ColumnIO.readSquareTable(files[k].toString()));
				x = new double [fft.N][fft.M];
				for (int i = 0; i < fft.N; i++)
					for (int j = 0; j < fft.M; j++)
					{
						x[i][j] = Math.sqrt(Math.pow(fft.getFHat2Re(i, j), 2) + Math.pow(fft.getFHat2Im(i, j), 2));
					}
				if (log)
					FieldOps.log(x);
				dx = x.length/factor;
				dy = x[0].length/factor;
				minx = (x.length-dx)/2;
				miny = (x[0].length-dy)/2;
		
				SRAW.writeImage(files[k].toString().substring(0, files[k].toString().length() - 4) + "fft_sm.bmp", x, min, max, minx, miny, dx, dy);
			}
		}
	}
	public static void writeBMPDirScale(String dir, boolean log, double min, double max, String keyword)
	{
		File[] files = new File(dir).listFiles();
		FFT2DSmall fft;
		double[][] x;
		for (int k = 0; k < files.length; k++)
		{
			if (files[k].toString().endsWith(".dat") && files[k].toString().contains(keyword))
			{
				fft = obtainFFT(ColumnIO.readSquareTable(files[k].toString()));
				x = new double [fft.N][fft.M];
				for (int i = 0; i < fft.N; i++)
					for (int j = 0; j < fft.M; j++)
					{
						x[i][j] = Math.sqrt(Math.pow(fft.getFHat2Re(i, j), 2) + Math.pow(fft.getFHat2Im(i, j), 2));
					}
				if (log)
					FieldOps.log(x);
				SRAW.writeImage(files[k].toString().substring(0, files[k].toString().length() - 4) + "fft_s.bmp", x, min, max);
			}
		}
	}
	public static void writeBMPDir2D(String dir, boolean log)
	{
		File[] files = new File(dir).listFiles();
		FFT2DSmall fft;
		for (int k = 0; k < files.length; k++)
		{
			if (files[k].toString().endsWith(".dat"))
			{
				fft = obtainFFT(ColumnIO.readSquareTable(files[k].toString()));
//				SRAW.writeImage(files[k].toString().substring(0, files[k].toString().length() - 4) + "fft.bmp", x);
				SRAW.writeImage(files[k].toString().substring(0, files[k].toString().length() - 4) + "fft_z.bmp", fft.fHat2(), false, log);
			}
		}

	}
	public static void supressModes(FFT2DSmall fft, int[] x, int[] y, double factor)
	{
		int N = fft.fHat.length;
		int M = fft.fHat[0].length;
		for (int i = 0; i < x.length; i++)
		{
			fft.fHat[(x[i]+N)%N][(y[i]+M)%M][0] /= factor;
			fft.fHat[(x[i]+N)%N][(y[i]+M)%M][1] /= factor;
		}
	}
	public static void supressModesTo(FFT2DSmall fft, int[] x, int[] y, double finalMag)
	{
		int N = fft.fHat.length;
		int M = fft.fHat[0].length;
		double nowmag;
		double factor;
		for (int i = 0; i < x.length; i++)
		{
			nowmag = Complex.mag(fft.fHat[(x[i]+N)%N][(y[i]+M)%M]);
			factor = nowmag/finalMag;
			fft.fHat[(x[i]+N)%N][(y[i]+M)%M][0] /= factor;
			fft.fHat[(x[i]+N)%N][(y[i]+M)%M][1] /= factor;
		}
	}
	public static void supressModesTo(FFT2DSmall fft, boolean[][] suppress, double finalMag)
	{
		int N = fft.fHat.length;
		int M = fft.fHat[0].length;
		double nowmag;
		double factor;
		for (int i = 0; i < fft.fHat.length; i++)
			for (int j = 0; j < fft.fHat[0].length; j++)
				if (suppress[i][j])
				{
					nowmag = Complex.mag(fft.fHat[i][j]);
					factor = nowmag/finalMag;
					fft.fHat[i][j][0] /= factor;
					fft.fHat[i][j][1] /= factor;
				}
	}
	public static void supressModes(FFT2DSmall fft, double[][] factors)
	{
		for (int i = 0; i < fft.fHat.length; i++)
			for (int j = 0; j < fft.fHat[0].length; j++)
				{
					fft.fHat[i][j][0] *= factors[i][j];
					fft.fHat[i][j][1] *= factors[i][j];
				}
	}
	public static double[][] getFourierFilteredIFFT(double[][] rawData, boolean[][] suppress)
	{
		FFT2DSmall fft = new FFT2DSmall(rawData);
		fft.doFFT();
		int N = fft.fHat.length;
		int M = fft.fHat[0].length;
		double nowmag;
		double phase;
//		double factor;
		double finalMag = FieldOps.magMin(fft.fHat);
		if (finalMag <= 0)
			finalMag = Double.MIN_VALUE;
		for (int i = 0; i < fft.fHat.length; i++)
			for (int j = 0; j < fft.fHat[0].length; j++)
				if (suppress[i][j])
				{
					phase = Complex.phase(fft.fHat[i][j]);
//					nowmag = Complex.mag(fft.fHat[i][j]);
//					factor = nowmag/finalMag;
					fft.fHat[i][j][0] = finalMag*Math.cos(phase);
					fft.fHat[i][j][1] = finalMag*Math.sin(phase);
				}
		fft.doIFFT();
		return FieldOps.getIndex(fft.f, 0);
//		return ArrayOps.toDouble(suppress);
	}
	public static void supressElseModesTo(FFT2DSmall fft, int[] x, int[] y, double finalMag)
	{
		int N = fft.fHat.length;
		int M = fft.fHat[0].length;
		double nowmag;
		double factor;
		boolean chosen = false;
		int xtrue, ytrue;
		for (int i = 0; i < N; i++)
			for (int j = 0; j < M; j++)
			{
				chosen = false;
				for (int k = 0; k < x.length; k++){
					xtrue = (x[k]+N)%N;
					ytrue = (y[k]+M)%M;
					if (xtrue == i && ytrue == j) chosen = true;
				}
				if (!chosen)
				{
					nowmag = Complex.mag(fft.fHat[(i+N)%N][(j+M)%M]);
					factor = nowmag/finalMag;
					if (factor != 0){
						fft.fHat[(i+N)%N][(j+M)%M][0] /= factor;
						fft.fHat[(i+N)%N][(j+M)%M][1] /= factor;
					}
				}
				else 
				{
//					System.out.println("Chosen");
				}
			}
	}

	public static double[][] obtainFFTmag(double[][] data)
	{
		FFT2DSmall fft = new FFT2DSmall(data);
		fft.doFFT();
		return FieldOps.magnitude(fft.fHat);
	}
	public static double[][] obtainFFTmagCent(double[][] data)
	{
		FFT2DSmall fft = new FFT2DSmall(data);
		fft.doFFT();
		double[][] ans = new double[fft.N][fft.M];
		for (int i = 0; i < fft.N; i++)
			for (int j = 0; j < fft.M; j++)
			{
				ans[i][j] = fft.getFHat2Mag(i, j);
			}
		return ans;
	}
	public static void putFFTmagCent(double[][] data, double[][] ans, boolean log)
	{
		FFT2DSmall fft = new FFT2DSmall(data);
		fft.doFFT();
		for (int i = 0; i < fft.N; i++)
			for (int j = 0; j < fft.M; j++)
			{
				ans[i][j] = fft.getFHat2Mag(i, j);
			}
		if (log)
			FieldOps.log(ans);
	}
	public static void obtainFFTCent(double[][] data, double[][][] ans)
	{
		FFT2DSmall fft = new FFT2DSmall(data);
		fft.doFFT();
		for (int i = 0; i < fft.N; i++)
			for (int j = 0; j < fft.M; j++)
			{
				ans[i][j][0] = fft.getFHat2Re(i, j);
				ans[i][j][1] = fft.getFHat2Im(i, j);
			}
	}
	public static void obtainFFTCent(double[][][] data, double[][][] ans)
	{
		FFT2DSmall fft = new FFT2DSmall(data);
		fft.doFFT();
		for (int i = 0; i < fft.N; i++)
			for (int j = 0; j < fft.M; j++)
			{
				ans[i][j][0] = fft.getFHat2Re(i, j);
				ans[i][j][1] = fft.getFHat2Im(i, j);
			}
	}
	//assume the result is real. use the real part
	public static void obtainIFFTCent(double[][][] ft, double[][] ans)
	{
		FFT2DSmall fft = new FFT2DSmall(ft, true);
		fft.doIFFT();
		for (int i = 0; i < fft.N; i++)
			for (int j = 0; j < fft.M; j++)
			{
				ans[i][j] = fft.f[i][j][0];
			}
	}
	public static void obtainIFFTCent(double[][][] ft, double[][][] ans)
	{
		FFT2DSmall fft = new FFT2DSmall(ft, true);
		fft.doIFFT();
		for (int i = 0; i < fft.N; i++)
			for (int j = 0; j < fft.M; j++)
			{
				ans[i][j][0] = fft.f[i][j][0];
				ans[i][j][1] = fft.f[i][j][1];
			}
	}
	public static void putIFFT(double[][][] ft, double[][][] ans, boolean centered)
	{
		FFT2DSmall fft = new FFT2DSmall(ft, centered);
		fft.doIFFT();
		for (int i = 0; i < fft.N; i++)
			for (int j = 0; j < fft.M; j++)
			{
				ans[i][j][0] = fft.f[i][j][0];
				ans[i][j][1] = fft.f[i][j][1];
			}
	}
	public static void putFFT(double[][][] data, double[][][] ans, boolean centered)
	{
		FFT2DSmall fft = new FFT2DSmall(data);
		fft.doFFT();
		if (centered)
			for (int i = 0; i < fft.N; i++)
				for (int j = 0; j < fft.M; j++)
				{
					ans[i][j][0] = fft.getFHat2Re(i, j);
					ans[i][j][1] = fft.getFHat2Im(i, j);
				}
		else
			for (int i = 0; i < fft.N; i++)
				for (int j = 0; j < fft.M; j++)
				{
					ans[i][j][0] = fft.fHat[i][j][0];
					ans[i][j][1] = fft.fHat[i][j][1];				
				}
	}
	public static void putFFT(double[][] data, double[][][] ans, boolean centered)
	{
		FFT2DSmall fft = new FFT2DSmall(data);
		fft.doFFT();
		if (centered)
			for (int i = 0; i < fft.N; i++)
				for (int j = 0; j < fft.M; j++)
				{
					ans[i][j][0] = fft.getFHat2Re(i, j);
					ans[i][j][1] = fft.getFHat2Im(i, j);
				}
		else
			for (int i = 0; i < fft.N; i++)
				for (int j = 0; j < fft.M; j++)
				{
					ans[i][j][0] = fft.fHat[i][j][0];
					ans[i][j][1] = fft.fHat[i][j][1];				
				}
	}
	public static void obtainFFTmagCent(double[][] data, double[][] ans)
	{
		FFT2DSmall fft = new FFT2DSmall(data);
		fft.doFFT();
		for (int i = 0; i < fft.N; i++)
			for (int j = 0; j < fft.M; j++)
			{
				ans[i][j] = fft.getFHat2Mag(i, j);
			}
	}
	public static void writeFFTReIm(String out, FFT2DSmall fft)
	{
		ColumnIO.writeBin(fft.fHat, 0, out + "re.dat");
		ColumnIO.writeBin(fft.fHat, 1, out + "im.dat");
	}
	public static void writeFFTBMP(String out, FFT2DSmall fft, boolean log, boolean justmag)
	{
		double[][] x = new double[fft.N][fft.M];
		FieldOps.magnitude(fft.fHat, x);
//		for (int i = 0; i < fft.N; i++)
//			for (int j = 0; j < fft.N; j++)
//			{
//				x[i][j] = Math.sqrt(Math.pow(fft.getFHat2Re(i, j), 2) + Math.pow(fft.getFHat2Im(i, j), 2));
//			}
		if (log)
			FieldOps.log(x);
		
		if (justmag)
		{
			SRAW.writeImage(out + ".bmp", x);
		}
		else
		{
			double[][] phase = new double[fft.N][fft.M];
			for (int i = 0; i < fft.N; i++)
				for (int j = 0; j < fft.M; j++)
				{
					phase[i][j] = FieldOps.atan(fft.getFHat2Re(i, j), fft.getFHat2Im(i, j));
				}
			double max = ArrayOps.max(x), min = ArrayOps.min(x);
			System.out.println("[" + min + ", " + max + "]");
			ColorScales.MYC2d scale = new ColorScales.MYC2d(max, min, 2*Math.PI);
			SRAW.writeImage(out + "fft.bmp", x, phase, scale);
		}
	}
	
	//This writes both the original and shifted bmp.
	public static void writeFFTBMPs(String out, FFT2DSmall fft, boolean log, boolean justmag)
	{
		double[][] x = new double[fft.N][fft.M];
		FieldOps.magnitude(fft.fHat, x);
		if (log)
			FieldOps.log(x);
		
		if (justmag)
		{
//			SRAW.writeImage(out + ".bmp", x, true, -22);
			SRAW.writeImage(out + ".bmp", x);
		}
		else
		{
			double[][] phase = FieldOps.phase(fft.fHat);
			double max = ArrayOps.max(x), min = ArrayOps.min(x);
			System.out.println("[" + min + ", " + max + "]");
			ColorScales.MYC2d scale = new ColorScales.MYC2d(max, min, 2*Math.PI);
			SRAW.writeImage(out + ".bmp", x, phase, scale);
		}
		
		//redo shifted
		for (int i = 0; i < fft.N; i++)
			for (int j = 0; j < fft.M; j++)
			{
				x[i][j] = Math.sqrt(Math.pow(fft.getFHat2Re(i, j), 2) + Math.pow(fft.getFHat2Im(i, j), 2));
			}
		if (log)
			FieldOps.log(x);
		
		if (justmag)
		{
//			SRAW.writeImage(out + "cent.bmp", x, true, -22);
			SRAW.writeImage(out + "cent.bmp", x);
		}
		else
		{
			double[][] phase = new double[fft.N][fft.M];
			for (int i = 0; i < fft.N; i++)
				for (int j = 0; j < fft.M; j++)
				{
					phase[i][j] = FieldOps.atan(fft.getFHat2Re(i, j), fft.getFHat2Im(i, j));
				}
			double max = ArrayOps.max(x), min = ArrayOps.min(x);
			System.out.println("[" + min + ", " + max + "]");
			ColorScales.MYC2d scale = new ColorScales.MYC2d(max, min, 2*Math.PI);
			SRAW.writeImage(out + "fft.bmp", x, phase, scale);
		}
	}
	public static void writeFFTBMPCent(String out, FFT2DSmall fft, boolean log, boolean justmag)
	{
		double[][] x = new double[fft.N][fft.M];
		for (int i = 0; i < fft.N; i++)
			for (int j = 0; j < fft.M; j++)
			{
				x[i][j] = Math.sqrt(Math.pow(fft.getFHat2Re(i, j), 2) + Math.pow(fft.getFHat2Im(i, j), 2));
			}
		if (log)
			FieldOps.log(x);
		
		if (justmag)
		{
//			SRAW.writeImage(out + "cent.bmp", x, true, -22);
			SRAW.writeImage(out, x);
		}
		else
		{
			double[][] phase = new double[fft.N][fft.M];
			for (int i = 0; i < fft.N; i++)
				for (int j = 0; j < fft.M; j++)
				{
					phase[i][j] = FieldOps.atan(fft.getFHat2Re(i, j), fft.getFHat2Im(i, j));
				}
			double max = ArrayOps.max(x), min = ArrayOps.min(x);
			System.out.println("[" + min + ", " + max + "]");
			ColorScales.MYC2d scale = new ColorScales.MYC2d(max, min, 2*Math.PI);
			SRAW.writeImage(out + "fft.bmp", x, phase, scale);
		}
	}
	public static BufferedImage getImageCent(double[][][] fHat, boolean log, boolean justmag, double hardMin)
	{
		double[][][] fHat2 = new double [fHat.length][fHat[0].length][2];
		FieldOps.shift(fHat, fHat2);
		double[][] x = new double[fHat.length][fHat[0].length];
		for (int i = 0; i < fHat.length; i++)
			for (int j = 0; j < fHat[0].length; j++)
			{
				x[i][j] = Complex.mag(fHat2[i][j]);
				if (x[i][j] < hardMin)
					x[i][j] = hardMin;
			}
		if (log)
			FieldOps.log(x);
		
		if (justmag)
			return ImageEditing.getBufferedImage(x);
		double[][] phase = new double[fHat.length][fHat[0].length];
		for (int i = 0; i < fHat.length; i++)
			for (int j = 0; j < fHat[0].length; j++)
			{
				phase[i][j] = Complex.phase(fHat2[i][j]);
			}
		double max = ArrayOps.max(x), min = ArrayOps.min(x);
		System.out.println("[" + min + ", " + max + "]");
		ColorScales.MYC2d scale = new ColorScales.MYC2d(max, min, 2*Math.PI);
		return ImageEditing.getBufferedImage(x, phase, scale);
	}
	public static void writeFFTBMPCent(BufferedImage image, FFT2DSmall fft, boolean log, boolean justmag)
	{
		double[][] x = new double[fft.N][fft.M];
		for (int i = 0; i < fft.N; i++)
			for (int j = 0; j < fft.M; j++)
			{
				x[i][j] = Math.sqrt(Math.pow(fft.getFHat2Re(i, j), 2) + Math.pow(fft.getFHat2Im(i, j), 2));
			}
		if (log)
			FieldOps.log(x);
		
		if (justmag)
		{
//			SRAW.writeImage(out + "cent.bmp", x, true, -22);
			SRAW.writeImage(image, x, null);
		}
		else
		{
			SRAW.writeImage(image, fft.fHat2(), null);
		}
	}

	public static void writeFFTBMP(String out, FFT2DSmall fft, boolean log, boolean justmag, double[][] x)
	{
		if (log)
			FieldOps.log(x);
		
		if (justmag)
		{
			SRAW.writeImage(out + ".bmp", x);
		}
		else
		{
			double[][] phase = new double[fft.N][fft.M];
			for (int i = 0; i < fft.N; i++)
				for (int j = 0; j < fft.M; j++)
				{
					phase[i][j] = FieldOps.atan(fft.getFHat2Re(i, j), fft.getFHat2Im(i, j));
				}
			double max = ArrayOps.max(x), min = ArrayOps.min(x);
			System.out.println("[" + min + ", " + max + "]");
			ColorScales.MYC2d scale = new ColorScales.MYC2d(max, min, 2*Math.PI);
			SRAW.writeImage(out + "fft.bmp", x, phase, scale);
		}
	}

	public static void writeFFTAll(String in, String out, boolean logscale, boolean justmag)
	{
		long t = System.currentTimeMillis();
		FFT2DSmall fft = obtainFFT(in);
		t = printTime(" to load and FFT", t);
		
		
//		writeFFTReIm(out, fft);
//		t = printTime(" to write dats", t);
		
		writeFFTBMPs(out, fft, logscale, justmag);
		t = printTime(" to write bmps", t);
	}
	public static void writeFFTBin(String in, String out)
	{
		long t = System.currentTimeMillis();
		
		double[][] data = ColumnIO.readSquareTable(in);
		System.out.println((System.currentTimeMillis() - t)/1000 + " to load.");
		t = System.currentTimeMillis();
		
		int N = data.length;
		
		SFT sft = new SFT(N);
		double[][][] ft = new double[2][N][N];
		System.out.println((System.currentTimeMillis() - t)/1000 + " making sft.");
		t = System.currentTimeMillis();
		
		
		sft.sft(data, ft);
		System.out.println((System.currentTimeMillis() - t)/1000 + " doing sft.");
		t = System.currentTimeMillis();
		
		ColumnIO.writeBin(ft[0], out + "re.dat");
		ColumnIO.writeBin(ft[1], out + "im.dat");

	}
	
	public static void main(String[] args)
	{
		String dir = "C:\\data\\lintrans\\run136topo6_620\\cutoff_pt4a\\devoutbin\\";
//		dir = "C:\\data\\lawlerhex\\superm\\devoutbin\\";
//		String in = dir + "topo.dat";
//		String out = dir + "topofft";
//		writeFFTAll(in, out, true, true);
//		writeBMPDir(dir, true, "resc");
		writeBMPComplexDir(dir, true, "16");
//		writeBMPDirScaleMiddle(dir, true, 3.17, 8.07, "didv", 2);
		
//		suppress5Write(in, out+".dat", 10000);
//		FieldUtil.makeImage(out+".dat", out+".bmp");
//		writeBMPDir("C:\\Users\\madhavanlab2011\\Documents\\PPTs\\1\\", true);
//		writeBMPDir2D("C:\\Users\\madhavanlab2011\\Documents\\PPTs\\1\\", true);
	}
	public static void writeSFTBMP(String in, String out, boolean justmag, boolean logmag)
	{
		long t = System.currentTimeMillis(), tf;
		
		System.out.println("start");
		double[][] data = SRAW.getData(in);
		System.out.println((System.currentTimeMillis() - t)/1000 + " to load.");
		t = System.currentTimeMillis();
		
		int N = data.length;
		
		SFT sft = new SFT(N);
		double[][][] ft = new double[N][N][2];
		System.out.println((System.currentTimeMillis() - t)/1000 + " making sft.");
		t = System.currentTimeMillis();
		
		
		sft.sft(data, ft);
		System.out.println((System.currentTimeMillis() - t)/1000 + " doing sft.");
		t = System.currentTimeMillis();
		
		double[][] mag = FieldOps.magnitude(ft);
		double[][] phase = FieldOps.phase(ft);
		
		if (logmag)
			FieldOps.log(mag);
		
		if (justmag)
		{
			SRAW.writeImage(out, mag);
			return;
		}
		else 
		{
			SRAW.writeImage(out, mag);
		}
	}
	public static void writeSFTFiles(String in, String out, boolean justmag, boolean logmag)
	{
		long t = System.currentTimeMillis();
		
		System.out.println("start");
//		double[][] data = SRAW.getData(in);
		double[][] data = ColumnIO.readSquareTable(in);
		System.out.println((System.currentTimeMillis() - t)/1000 + " to load.");
		t = System.currentTimeMillis();
		
		int N = data.length;
		
		SFT sft = new SFT(N);
		double[][][] ft = new double[2][N][N];
		System.out.println((System.currentTimeMillis() - t)/1000 + " making sft.");
		t = System.currentTimeMillis();
		
		
		sft.sft(data, ft);
		System.out.println((System.currentTimeMillis() - t)/1000 + " doing sft.");
		t = System.currentTimeMillis();
		
		ColumnIO.writeBin(ft[0], out + "re.dat");
		ColumnIO.writeBin(ft[1], out + "im.dat");
	}
	
	public static void suppressAvgWrite(String in, String out, double factor)
	{
		suppressModesWrite(in, new int[] {0}, new int[] {0}, factor, out);
	}
	public static void suppress5Write(String in, String out, double factor)
	{
		suppressModesWrite(in, new int[] {0, 1, 1, -1, -1}, new int[] {0, 1, -1, 1, -1}, factor, out);
	}
	public static void suppressModesWrite(String in, int[] x, int[] y, double factor, String out)
	{
		FFT2DSmall fft = obtainFFT(in);
		supressModes(fft, x, y, factor);
		fft.doIFFT();
		ColumnIO.writeBin(fft.f, 0, out);
	}
	
	public static long printTime(String phrase, long t)
	{
		System.out.println((System.currentTimeMillis() - t)/1000 + phrase);
		return System.currentTimeMillis();
	}
	public static double[] get1DFFTMag(double[] data)
	{
		DoubleFFT_1D fft = null;
		double[] temp = null;
		double[] ans = new double [data.length];
		fft = new DoubleFFT_1D(data.length);
		temp = new double [2*data.length];
		for (int i = 0; i < data.length; i++)
		{
			temp[2*i] = data[i];
			temp[2*i+1] = 0;
		}
		fft.complexForward(temp);
		for (int i = 0; i < data.length; i++)
		{
			ans[i] = Complex.mag(temp[2*i], temp[2*i+1]);
		}
		return ans;
	}
	/**
	 * The fmin, fmax constitute the array of bounds of the regions to be filtered.
	 * @param values
	 * @param fmin
	 * @param fmax
	 * @return
	 */
	public static double[] getFourierFiltered(double[] values, int[] fmin, int[] fmax) {
		int nlayers = values.length;
		double magmin = ArrayOps.min(FFTOps.get1DFFTMag(values));
		double[] fftz = FFTOps.get1DFFTComplex(values);
		double[] ifftReal = new double [nlayers];
		//We have to assume here that lastT is the spectrum whose average we desire to modify.
		for (int k = 0; k < fmin.length; k++){
			if (fmax[k] > nlayers/2) {
				System.out.println("Error. " + fmax[k] + " was outside of the range. Returning original values.");
				FFTOps.putIFFTReal(fftz.clone(), ifftReal);
	//				return;
			}
	
			
			int f2max = nlayers-fmax[k];
			int f2min = nlayers-fmin[k];
			double ratio;
			for (int i = 0; i < values.length; i++){
				if ((i >= fmin[k] && i <= fmax[k]) || (i >= f2max && i <= f2min))
				{
					ratio = Complex.mag(fftz[2*i], fftz[2*i+1])/magmin;
//					if (ratio < 100) ratio = 100;
					fftz[2*i] /= ratio;
					fftz[2*i+1] /= ratio;
	//					System.out.println(i + "\t" + ratio + "\t" + fmin + "\t" + fmax + "\t" + magmin + "\t" + Complex.mag(fftClone[2*i], fftClone[2*i+1]));
				}
			}
		}
		FFTOps.putIFFTReal(fftz, ifftReal);
		return ifftReal;
	}

	public static double[] get1DFFTComplex(double[] data) //this returns the temporary array
	{
		DoubleFFT_1D fft = null;
		double[] temp = null;
		fft = new DoubleFFT_1D(data.length);
		temp = new double [2*data.length];
		for (int i = 0; i < data.length; i++)
		{
			temp[2*i] = data[i];
			temp[2*i+1] = 0;
		}
		fft.complexForward(temp);
		return temp;
	}
	public static double[] getIFFTReal(double[] temp)
	{
		
		DoubleFFT_1D fft = null;
		fft = new DoubleFFT_1D(temp.length/2);
		double[] ans = new double [temp.length/2];
		fft.complexInverse(temp, true);
		for (int i = 0; i < temp.length/2; i++)
		{
			ans[i] = temp[2*i];
		}
		return ans;
	}
	public static void putIFFTReal(double[] temp, double[] ans)
	{
		
		DoubleFFT_1D fft = null;
		fft = new DoubleFFT_1D(temp.length/2);
		fft.complexInverse(temp, true);
		for (int i = 0; i < temp.length/2; i++)
		{
			ans[i] = temp[2*i];
		}
	}

	//This writes only the vicinity of the bragg peaks to a file.
	//hex so bragg3.
//	public static void writeBraggOnly(String dir, String namein, String nameout, String bragg)
//	{
//		FFT2DSmall f = new FFT2DSmall(ColumnIO.readSquareTable(dir + namein + ".dat"));
//		f.doFFT();
//		double[][][] fhat2 = f.fHat2();
//		int[][] braggpts = DataManip.readBragg(dir + bragg);
//		int[] bragg3 = {bragg[1][0] - bragg[0][0], bragg[1][1] - bragg[0][1]};
//	}
	
}
