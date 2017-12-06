package main;

import image.imageIO;

import java.awt.Color;
import java.awt.image.BufferedImage;

import util.ArrayOps;
import util.Complex;
import util.color.ColorScale;
import util.color.ColorScale2d;
import util.color.ColorScales;
import util.fileops.ColumnIO;
import util.FieldOps;

public class SRAW {

	//In this method note that j corresponds to x, i to y, since the 1-d array
	//is written in sequential order
	public static void writeBMP(String in, String out)
	{
		writeImage(out, getData(in));
	}
	
	public static double[][] getData(String in)
	{
		return ColumnIO.readAllColumns(new java.io.File(in), null);
	}
//	public static double[][] getData(String in, String delimiter)
//	{
//		return ColumnIO.readAllColumnsB(new java.io.File(in), null, "\t");
//	}

	public static void writeImage(String out, double[][] data)
	{
		double max = ArrayOps.max(data), min = ArrayOps.min(data);
		System.out.print("[" + min + ", " + max + "]\t");
		ColorScales.LinearBRYW scale = new ColorScales.LinearBRYW(max, min);
		int[] rgba = new int [data.length*data[0].length];
		java.awt.Color c;
		int n = 0;
		for (int i = 0; i < data[0].length; i++)
			for (int j = 0; j < data.length; j++)
			{
				c = scale.of(data[j][i]);
				rgba[n] = c.getRGB();
				n++;
			}
		imageIO.saveImage3(rgba, out, data.length, data[0].length);
	}
	//This assumes that the buffered image is in the same transposed form as the array would have been
	public static void writeImage(String out, boolean[][] data)
	{
		int[] rgba = new int [data.length*data[0].length];
		java.awt.Color c;
		int n = 0;
		for (int i = 0; i < data[0].length; i++)
			for (int j = 0; j < data.length; j++)
			{
				c = data[j][i] ? Color.white : Color.black;
				rgba[n] = c.getRGB();
				n++;
			}
		imageIO.saveImage3(rgba, out, data.length, data[0].length);
	}
	public static void writeInflatedImage(String out, double[][] data, int xpixperpoint, int ypixperpoint)
	{
		double max = ArrayOps.max(data), min = ArrayOps.min(data);
		System.out.print("[" + min + ", " + max + "]\t");
		ColorScales.LinearBRYW scale = new ColorScales.LinearBRYW(max, min);
		int[] rgba = new int [data.length*data[0].length*xpixperpoint*ypixperpoint];
		java.awt.Color c;
		int n = 0;
		for (int i = 0; i < data[0].length*ypixperpoint; i++)
			for (int j = 0; j < data.length*xpixperpoint; j++)
			{
				c = scale.of(data[j/xpixperpoint][i/ypixperpoint]);
				rgba[n] = c.getRGB();
				n++;
			}
		imageIO.saveImage3(rgba, out, data.length*xpixperpoint, data[0].length*ypixperpoint);
	}
	public static void writeImagPhase(String out, double[][][] data)
	{
		double max = Math.PI*2, min = 0;
		System.out.print("[" + min + ", " + max + "]\t");
		ColorScales.LinearBRYW scale = new ColorScales.LinearBRYW(max, min);
		int[] rgba = new int [data.length*data[0].length];
		java.awt.Color c;
		int n = 0;
		for (int i = 0; i < data[0].length; i++)
			for (int j = 0; j < data.length; j++)
			{
				c = scale.of(Complex.phase(data[j][i]));
				rgba[n] = c.getRGB();
				n++;
			}
		imageIO.saveImage3(rgba, out, data.length, data[0].length);
	}
	public static void writeImage(String out, double[][] data, double min, double max)
	{
		System.out.print("[" + min + ", " + max + "]\t");
		ColorScales.LinearBRYW scale = new ColorScales.LinearBRYW(max, min);
		int[] rgba = new int [data.length*data[0].length];
		java.awt.Color c;
		int n = 0;
		for (int i = 0; i < data[0].length; i++)
			for (int j = 0; j < data.length; j++)
			{
				c = scale.of(data[j][i]);
				rgba[n] = c.getRGB();
				n++;
			}
		imageIO.saveImage3(rgba, out, data.length, data[0].length);
	}
	public static void writeImage(String out, double[][] data, boolean minext, double minpar)
	{
		double max = ArrayOps.max(data), min = minext ? minpar : ArrayOps.min(data);
		System.out.print("[" + min + ", " + max + "]\t");
		ColorScales.LinearBRYW scale = new ColorScales.LinearBRYW(max, min);
		int[] rgba = new int [data.length*data[0].length];
		java.awt.Color c;
		int n = 0;
		for (int i = 0; i < data[0].length; i++)
			for (int j = 0; j < data.length; j++)
			{
				c = scale.of(data[j][i]);
				rgba[n] = c.getRGB();
				n++;
			}
		imageIO.saveImage3(rgba, out, data.length, data[0].length);
	}
	public static void writeImage(String out, int[][] data)
	{
		double max = ArrayOps.max(data), min = ArrayOps.min(data);
		System.out.print("[" + min + ", " + max + "]\t");
		ColorScales.LinearBRYW scale = new ColorScales.LinearBRYW(max, min);
		int[] rgba = new int [data.length*data[0].length];
		java.awt.Color c;
		int n = 0;
		for (int i = 0; i < data[0].length; i++)
			for (int j = 0; j < data.length; j++)
			{
				c = scale.of(data[j][i]);
				rgba[n] = c.getRGB();
				n++;
			}
		imageIO.saveImage3(rgba, out, data.length, data[0].length);
	}
	public static void writeImage(String out, double[][] data, int x, int y, int dx, int dy)
	{
		ColorScales.LinearBRYW scale = new ColorScales.LinearBRYW(ArrayOps.max(data), ArrayOps.min(data));
		int[] rgba = new int [dx*dy];
		java.awt.Color c;
		int n = 0;
		for (int i = y; i < y+dy; i++)
			for (int j = x; j < x+dx; j++)
			{
				c = scale.of(data[j][i]);
				rgba[n] = c.getRGB();
				n++;
			}
		imageIO.saveImage3(rgba, out, dx, dy);
	}
	public static void writeImage(String out, double[][] data, double min, double max, int x, int y, int dx, int dy)
	{
		ColorScales.LinearBRYW scale = new ColorScales.LinearBRYW(max, min);
		int[] rgba = new int [dx*dy];
		java.awt.Color c;
		int n = 0;
		for (int i = y; i < y+dy; i++)
			for (int j = x; j < x+dx; j++)
			{
				c = scale.of(data[j][i]);
				rgba[n] = c.getRGB();
				n++;
			}
		imageIO.saveImage3(rgba, out, dx, dy);
	}
//	public static void writeImage(String out, double[][] data)
//	{
//		ColorScales.LinearBRYW scale = new ColorScales.LinearBRYW(ArrayOps.max(data), ArrayOps.min(data));
//		int[][] rgba = new int [data.length*data[0].length][4];
//		java.awt.Color c;
//		int n = 0;
//		for (int i = 0; i < data[0].length; i++)
//			for (int j = 0; j < data.length; j++)
//			{
//				c = scale.of(data[j][i]);
//				rgba[n][0] = c.getRed();
//				rgba[n][1] = c.getGreen();
//				rgba[n][2] = c.getBlue();
//				rgba[n][3] = c.getAlpha();
//				n++;
//			}
//		imageIO.saveImage2(rgba, out, data.length, data[0].length);
//	}
	public static void writeImage(String out, double[][][] data, int thirdindex)
	{
		double max = ArrayOps.max(data, thirdindex), min = ArrayOps.min(data, thirdindex);
		System.out.print("[" + min + ", " + max + "]\t");
		ColorScales.LinearBRYW scale = new ColorScales.LinearBRYW(max, min);
		int[] rgba = new int [data.length*data[0].length];
		java.awt.Color c;
		int n = 0;
		for (int i = 0; i < data[0].length; i++)
			for (int j = 0; j < data.length; j++)
			{
				c = scale.of(data[j][i][thirdindex]);
				rgba[n] = c.getRGB();
				n++;
			}
		imageIO.saveImage3(rgba, out, data.length, data[0].length);
		
	}
	public static void writeImage(String out, double[][] mag, double[][] phase, ColorScale2d scale)
	{
		int[][] rgba = new int [mag.length*mag[0].length][4];
		java.awt.Color c;
		int n = 0;
		for (int i = 0; i < mag[0].length; i++)
			for (int j = 0; j < mag.length; j++)
			{
				c = scale.of(mag[j][i], phase[j][i]);
				rgba[n][0] = c.getRed();
				rgba[n][1] = c.getGreen();
				rgba[n][2] = c.getBlue();
				rgba[n][3] = c.getAlpha();
				n++;
			}
		imageIO.saveImage2(rgba, out, mag.length, mag[0].length);
	}
	public static void writeImage(String out, double[][][] complex, boolean minzero)
	{
		int[] rgba = new int [complex.length*complex[0].length];
		java.awt.Color c;
		int n = 0;
		double max = FieldOps.magMax(complex);
		double min = minzero ? 0 : FieldOps.magMin(complex);
		System.out.print("[" + min + ", " + max + "]\t");
		ColorScales.MYC2d scale = new ColorScales.MYC2d(max, min, 2*Math.PI);
		for (int i = 0; i < complex[0].length; i++)
			for (int j = 0; j < complex.length; j++)
			{
				c = scale.of(Complex.mag(complex[j][i]), FieldOps.atan(complex[j][i][0], complex[j][i][1]));
				rgba[n] = c.getRGB();
				n++;
			}
		System.out.print("[" + min + ", " + max + "]");
		imageIO.saveImage3(rgba, out, complex.length, complex[0].length);
	}
	public static void writeImage(String out, double[][][] complex, boolean minzero, boolean log)
	{
		int[] rgba = new int [complex.length*complex[0].length];
		java.awt.Color c;
		int n = 0;
		double altmin = !log ? 0 : Math.log(Math.sqrt(Double.MIN_VALUE));
		double max = log ? Math.log(FieldOps.magMax(complex)) : FieldOps.magMax(complex);
		double min = minzero ? altmin : log ? Math.log(FieldOps.magMin(complex)) : FieldOps.magMin(complex);
		System.out.print("[" + min + ", " + max + "]\t");
		ColorScales.MYC2d scale = new ColorScales.MYC2d(max, min, 2*Math.PI);
		double d;
		for (int i = 0; i < complex[0].length; i++)
			for (int j = 0; j < complex.length; j++)
			{
				d = log ? Math.log(Complex.mag(complex[j][i])) : Complex.mag(complex[j][i]);
				if (minzero && log && d == Double.NaN) d = altmin;
				c = scale.of(d, FieldOps.atan(complex[j][i][0], complex[j][i][1]));
				rgba[n] = c.getRGB();
				n++;
			}
		imageIO.saveImage3(rgba, out, complex.length, complex[0].length);
	}
	public static void writeImage(BufferedImage out, double[][][] complex, ColorScale2d scale)
	{
		java.awt.Color c;
		if (scale == null)
			scale = new ColorScales.MYC2d(FieldOps.magMax(complex), FieldOps.magMin(complex), 2*Math.PI);
		for (int i = 0; i < complex[0].length; i++)
			for (int j = 0; j < complex.length; j++)
			{
				c = scale.of(Complex.mag(complex[j][i]), FieldOps.atan(complex[j][i][0], complex[j][i][1]));
				out.setRGB(j, i, c.getRGB());
			}
	}
	public static void writeImage(BufferedImage out, double[][][] complex, ColorScale2d scale, int sizeratio)
	{
		java.awt.Color c;
		int m, n;
		if (scale == null)
			scale = new ColorScales.MYC2d(FieldOps.magMax(complex), FieldOps.magMin(complex), 2*Math.PI);
		for (int i = 0; i < complex[0].length; i++)
			for (int j = 0; j < complex.length; j++)
			{
				c = scale.of(Complex.mag(complex[j][i]), FieldOps.atan(complex[j][i][0], complex[j][i][1]));
				for (m = 0; m < sizeratio; m++)
					for (n = 0; n < sizeratio; n++)
						out.setRGB(sizeratio*j+m, sizeratio*i+n, c.getRGB());
//				out.setRGB(j, i, c.getRGB());
			}
	}
	public static void writeImage(BufferedImage out, double[][] data, ColorScale scale)
	{
		java.awt.Color c;
		if (scale == null) scale = ColorScales.getNew(data, 0);
		for (int i = 0; i < data[0].length; i++)
			for (int j = 0; j < data.length; j++)
			{
				c = scale.of(data[j][i]);
				out.setRGB(j, i, c.getRGB());
			}
	}
	public static void writeImage(BufferedImage out, double[][] data, int scaleIndex)
	{
		java.awt.Color c;
		ColorScale scale = ColorScales.getNew(data, scaleIndex);
		for (int i = 0; i < data[0].length; i++)
			for (int j = 0; j < data.length; j++)
			{
				c = scale.of(data[j][i]);
				out.setRGB(j, i, c.getRGB());
			}
	}
	public static void writeImage(BufferedImage out, double[][] data, int scaleIndex,  double min, double max)
	{
		java.awt.Color c;
		ColorScale scale = ColorScales.getNew(min, max, scaleIndex);
		for (int i = 0; i < data[0].length; i++)
			for (int j = 0; j < data.length; j++)
			{
				c = scale.of(data[j][i]);
				out.setRGB(j, i, c.getRGB());
			}
	}
	public static void writeImage(BufferedImage out, double[][] data, ColorScale scale, int sizeratio)
	{
		java.awt.Color c;
		int i, j, m, n;
		for (i = 0; i < data.length; i++)
			for (j = 0; j < data[0].length; j++)
			{
				c = scale.of(data[i][j]);
				for (m = 0; m < sizeratio; m++)
					for (n = 0; n < sizeratio; n++)
						out.setRGB(sizeratio*i+m, sizeratio*j+n, c.getRGB());
			}
	}
	public static void writeImage(String out, double[][] data, boolean log)
	{
		int[] rgba = new int [data.length*data[0].length];
		java.awt.Color c;
		int n = 0;
		double max = log ? Math.log(ArrayOps.max(data)) : ArrayOps.max(data);
		double min = log ? Math.log(ArrayOps.min(data)) : ArrayOps.min(data);
		System.out.print("[" + min + ", " + max + "]\t");
		ColorScales.LinearBRYW scale = new ColorScales.LinearBRYW(max, min);
		for (int i = 0; i < data[0].length; i++)
			for (int j = 0; j < data.length; j++)
			{
				c = scale.of(log ? Math.log(data[j][i]) : data[j][i]);
				rgba[n] = c.getRGB();
				n++;
			}
		imageIO.saveImage3(rgba, out, data.length, data[0].length);
	}
	public static void writeImage(String out, double[][] data, ColorScale scale)
	{
		int[][] rgba = new int [data.length*data[0].length][4];
		java.awt.Color c;
		int n = 0;
		for (int i = 0; i < data[0].length; i++)
			for (int j = 0; j < data.length; j++)
			{
				c = scale.of(data[j][i]);
				rgba[n][0] = c.getRed();
				rgba[n][1] = c.getGreen();
				rgba[n][2] = c.getBlue();
				rgba[n][3] = c.getAlpha();
				n++;
			}
		imageIO.saveImage2(rgba, out, data.length, data[0].length);
	}
	public static void writeImage(String out, double[][] data, ColorScale scale, boolean renormalize)
	{
		if(renormalize)
			scale.renormalize(ArrayOps.max(data), ArrayOps.min(data));
		int[][] rgba = new int [data.length*data[0].length][4];
		java.awt.Color c;
		int n = 0;
		for (int i = 0; i < data[0].length; i++)
			for (int j = 0; j < data.length; j++)
			{
				c = scale.of(data[j][i]);
				rgba[n][0] = c.getRed();
				rgba[n][1] = c.getGreen();
				rgba[n][2] = c.getBlue();
				rgba[n][3] = c.getAlpha();
				n++;
			}
		imageIO.saveImage2(rgba, out, data.length, data[0].length);
	}
	public static void writeImage(String out, BufferedImage image)
	{
		int[] rgb = new int[image.getWidth()*image.getHeight()];
		rgb = image.getRGB(0, 0, image.getWidth(), image.getHeight(), rgb, 0, image.getWidth());
		imageIO.saveImage3(rgb, out, image.getWidth(), image.getHeight());
	}
	
	public static void writeBMPSmall(String in, String out, int xmin, int ymin, int xmax, int ymax)
	{
		int dx = xmax - xmin, dy = ymax - ymin;
		double[][] data = ColumnIO.readAllColumns(new java.io.File(in), null);
		ColorScales.LinearBRYW scale = new ColorScales.LinearBRYW(ArrayOps.max(data), ArrayOps.min(data));
		int[][] rgba = new int [dx*dy][4];
		java.awt.Color c;
		int n = 0;
		for (int i = xmin; i < xmax; i++)
			for (int j = ymin; j < ymax; j++)
			{
				c = scale.of(data[i][j]);
				rgba[n][0] = c.getRed();
				rgba[n][1] = c.getGreen();
				rgba[n][2] = c.getBlue();
				rgba[n][3] = c.getAlpha();
				n++;
			}
		imageIO.saveImage2(rgba, out, dx, dy);
		//		BufferedImage image = new BufferedImage(data.length, data[0].length, BufferedImage.TYPE_4BYTE_ABGR);
//		for (int i = 0; i < data.length; i++)
//			for (int j = 0; j < data[0].length; j++)
//			{
//				image.setRGB(i, j, scale.of(data[i][j]).getRGB());
//			}
//		CS585Image imageSave = new CS585Image(image);
//		imageSave.save(out); 
	}
	
	public static void main(String[] args)
	{
		String name = "todan";
		String txt = "C:\\data\\";
		String bmp = "C:\\data\\";
		double[][] data = getData(txt + name + ".txt");
		System.out.println("got");
//		double[][][] grad = FieldOps.gradient(data);
//		System.out.println("got");
		writeImage(bmp + name + ".bmp", data);
//		System.out.println("got");
		
	}
	
	public static void printTot(String in)
	{
		double[][] data = ColumnIO.readAllColumns(new java.io.File(in), null);
		System.out.println("" + util.ArrayOps.sum(data));
	}
	public static double getR(String in, int x, int y)
	{
		double[][] data = ColumnIO.readAllColumns(new java.io.File(in), null);
		double r = 0;
		for (int i = 0; i < data.length; i++)
			for(int j = 0; j < data[i].length; j++)
				r += Math.sqrt((i-x)*(i-x) + (j-y)*(j-y))*data[i][j];
		return r/util.ArrayOps.sum(data);
	}
	public static void printInfo(String in, int k, int x, int y)
	{
		double[][] data = ColumnIO.readAllColumns(new java.io.File(in), null);
		double r = 0, rSq = 0, drSq = 0, sum = 0;
		double xcm = 0, ycm = 0;
		for (int i = 0; i < data.length; i++)
			for(int j = 0; j < data[i].length; j++)
			{
				sum += data[i][j];
				drSq = (i-x)*(i-x) + (j-y)*(j-y);
				rSq += drSq*data[i][j];
				r += Math.sqrt(drSq)*data[i][j];
				xcm += i*data[i][j];
				ycm += j*data[i][j];
			}
		rSq /= sum;
		r /= sum;
		xcm /= sum;
		ycm /= sum;
		
		System.out.println(k + "\t" + r + "\t" + rSq +"\t" + sum + "\t" + xcm + "\t" + ycm);
	}
	public static void buildRHist(String in, int x, int y)
	{
		double[][] data = ColumnIO.readAllColumns(new java.io.File(in), null);
		double[] occupancy = new double[500];
		for (int i = 0; i < 500; i++)
		{
			occupancy[i] = 0;
		}
		double r = 0, rSq = 0, drSq = 0;
		for (int i = 0; i < data.length; i++)
			for(int j = 0; j < data[i].length; j++)
			{
				drSq = (i-x)*(i-x) + (j-y)*(j-y);
				r= Math.sqrt(drSq);
				occupancy[(int)(r/2)] += r != 0 ? data[i][j]/r : 0;
			}
		for (int i =0 ; i < 500; i++)
			System.out.println(i*2 + "\t" + occupancy[i]);
	}
	public static void buildRHist2(String in, int x, int y)
	{
		double[][] data = ColumnIO.readAllColumns(new java.io.File(in), null);
		double[] occupancy = new double[500];
		for (int i = 0; i < 500; i++)
		{
			occupancy[i] = 0;
		}
		double r = 0, rSq = 0, drSq = 0;
		for (int i = 0; i < data.length; i++)
			for(int j = 0; j < data[i].length; j++)
			{
				drSq = (i-x)*(i-x) + (j-y)*(j-y);
				r = drSq;
				occupancy[(int)(r/2)] += r != 0 ? data[i][j]/r : 0;
			}
		for (int i =0 ; i < 500; i++)
			System.out.println(i*2 + "\t" + occupancy[i]);
	}
	public static void buildRHist3(String in, int x, int y)
	{
		double[][] data = ColumnIO.readAllColumns(new java.io.File(in), null);
		double[] occupancy = new double[500];
		int[] cells = new int[500];
		for (int i = 0; i < 500; i++)
		{
			occupancy[i] = 0;
			cells[i] = 0;
		}
		double r = 0, drSq = 0;
		for (int i = 0; i < data.length; i++)
			for(int j = 0; j < data[i].length; j++)
			{
				drSq = (i-x)*(i-x) + (j-y)*(j-y);
				r= Math.sqrt(drSq);
				cells[(int)(r/2)]++;
			}
		for (int i = 0; i < data.length; i++)
			for(int j = 0; j < data[i].length; j++)
			{
				drSq = (i-x)*(i-x) + (j-y)*(j-y);
				r= Math.sqrt(drSq);
				occupancy[(int)(r/2)] += 1/(double)cells[(int)(r/2)];
			}
		for (int i =0 ; i < 500; i++)
			System.out.println(i + "\t" + occupancy[i] + "\t" + cells[i]);
	}
	
	public static double[] readFile(String in, int lineNum)
	{
		return ColumnIO.readNthLine(new java.io.File(in), lineNum);
	}
	
	public static void renderPathsIntoFields(String txt)
	{
		double[] x, y;
		int N = 500, half = N/2;
		int[][] occ = new int [N][N];
		for (int i = 0; i < N; i++)
			for (int j = 0; j < N; j++)
				occ[i][j] = 0;
		
		double min, max;
		String name = "x data 150000";
		String name2 = "y data 150000";
		String out = "field at ";
		for (int i = 0; i < 1000; i++){
			x = readFile(txt + name + ".txt", i+1);
			System.out.print("x read; ");
			y = readFile(txt + name2+ ".txt", i+1);
			System.out.print("y read; ");
			for (int j = 0; j < x.length; j++)
			{
				occ[(int)(x[j] + half)][(int)(y[j] + half)]++;
			}
			
			System.out.print("field seeded; ");
			out = "field n=" + i;
			ColumnIO.writeTable(occ, out + ".txt", txt);

			for (int k = 0; k < N; k++)
				for (int j = 0; j < N; j++)
					occ[k][j] = 0;
			System.out.println("Field written " + i);
//			min = ArrayOps.min(data);
//			max = ArrayOps.max(data);
		}
	
	}
	
	public static void renameFiles(String dir, String newName)
	{
		java.io.File[] list = new java.io.File(dir).listFiles();
		int index, nsteps, no;
		for (int i = 0; i < list.length; i++)
			if(list[i].toString().contains("nsteps")){
			nsteps = list[i].toString().indexOf('=') + 2;
			no = list[i].toString().indexOf("psmall") - 1;
			index = Integer.parseInt(list[i].toString().substring(nsteps, no));
			System.out.println(nsteps + "\t" + no + "\t" + list[i].toString() + "\t" + index);
			list[i].renameTo(new java.io.File(dir + newName + index + ".bmp"));
			}
	}
	
	
}
