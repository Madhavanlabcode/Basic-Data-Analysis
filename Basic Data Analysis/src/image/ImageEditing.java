package image;

import impurity.PointImp;

import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.Toolkit;
import java.awt.datatransfer.DataFlavor;
import java.awt.datatransfer.Transferable;
import java.awt.datatransfer.UnsupportedFlavorException;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.ArrayList;

import javax.imageio.ImageIO;
import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import main.SRAW;

import org.apache.commons.math3.analysis.interpolation.BicubicSplineInterpolatingFunction;

import schrodinger.MovieMaker;
import util.ArrayOps;
import util.Complex;
import util.FieldOps;
import util.NumFormat;
import util.color.ColorScale;
import util.color.ColorScale2d;
import util.color.ColorScales;
import util.fileops.FileOps;
import util.fileops.Topomap;
import util.geom.AtomicCoordinatesSet;
import util.geom.Distance;
import util.matrix.Matrix;

public class ImageEditing {

	
	
	
	public static BufferedImage getBufferedImage(double[][] data)
	{
		BufferedImage image = new BufferedImage(data.length, data[0].length, BufferedImage.TYPE_INT_RGB);
		ColorScales.LinearBRYW scale = new ColorScales.LinearBRYW(ArrayOps.max(data), ArrayOps.min(data));
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++)
				image.setRGB(i, j, scale.of(data[i][j]).getRGB());
		return image;

	}
	public static BufferedImage getBufferedImage(int[][] data)
	{
		BufferedImage image = new BufferedImage(data.length, data[0].length, BufferedImage.TYPE_INT_RGB);
		ColorScales.LinearBRYW scale = new ColorScales.LinearBRYW(ArrayOps.max(data), ArrayOps.min(data));
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++)
				image.setRGB(i, j, scale.of(data[i][j]).getRGB());
		return image;

	}
	public static BufferedImage getBufferedImage(double[][] data, int scaleIndex)
	{
		BufferedImage image = new BufferedImage(data.length, data[0].length, BufferedImage.TYPE_INT_RGB);
		ColorScale scale = ColorScales.getNew(data, scaleIndex);
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++)
				image.setRGB(i, j, scale.of(data[i][j]).getRGB());
		return image;
	}
	public static BufferedImage getBufferedImage(double[][] data, ColorScale scale)
	{
		BufferedImage image = new BufferedImage(data.length, data[0].length, BufferedImage.TYPE_INT_RGB);
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++)
				image.setRGB(i, j, scale.of(data[i][j]).getRGB());
		return image;
	}
	public static BufferedImage getBufferedImage_NotLargerThan(double[][] data, ColorScale scale, int sx, int sy)
	{
		if (data.length <= sx && data[0].length <= sy)
			return getBufferedImage(data, scale);
		
		BufferedImage image = new BufferedImage(sx, sy, BufferedImage.TYPE_INT_RGB);
		double[] x = ArrayOps.generateArrayInclBoth(0, data.length, sx);
		double[] y = ArrayOps.generateArrayInclBoth(0, data[0].length, sy);
		double mean = FieldOps.mean(data);
		for (int i = 0; i < sx; i++)
			for (int j = 0; j < sy; j++)
				image.setRGB(i, j, scale.of(FieldOps.getValueAt(data, x[i], y[j], mean)).getRGB());
		return image;
	}
	public static BufferedImage getBufferedImage_FixedSizeInt(double[][] data, ColorScale scale, int sx, int sy)
	{
		BufferedImage image = new BufferedImage(sx, sy, BufferedImage.TYPE_INT_RGB);
		double[] x = ArrayOps.generateArrayInclBoth(0-0.5, data.length-1+0.5, sx);
		double[] y = ArrayOps.generateArrayInclBoth(0-0.5, data[0].length-1+0.5, sy);
		double mean = FieldOps.mean(data);
		for (int i = 0; i < sx; i++)
			for (int j = 0; j < sy; j++)
				image.setRGB(i, j, scale.of(FieldOps.evaluateAtInt(data, x[i], y[j], mean)).getRGB());
		return image;
	}
	public static void putBufferedImage_FixedSizeInt(double[][] data, ColorScale scale, BufferedImage image)
	{
		int sx = image.getWidth(), sy = image.getHeight();
		double[] x = ArrayOps.generateArrayInclBoth(0, data.length-1, sx);
		double[] y = ArrayOps.generateArrayInclBoth(0, data[0].length-1, sy);
		double mean = FieldOps.mean(data);
		for (int i = 0; i < sx; i++)
			for (int j = 0; j < sy; j++)
				image.setRGB(i, j, scale.of(FieldOps.evaluateAtInt(data, x[i], y[j], mean)).getRGB());
	}
	public static void putBufferedImage_FixedSizeIntTranspose(double[][] data, ColorScale scale, BufferedImage image)
	{
		int sx = image.getWidth(), sy = image.getHeight();
		double[] x = ArrayOps.generateArrayInclBoth(0, data.length-1, sy);
		double[] y = ArrayOps.generateArrayInclBoth(0, data[0].length-1, sx);
		double mean = FieldOps.mean(data);
		for (int i = 0; i < sy; i++)
			for (int j = 0; j < sx; j++)
				image.setRGB(j, i, scale.of(FieldOps.evaluateAtInt(data, x[i], y[j], mean)).getRGB());
	}
	public static BufferedImage getBufferedImage(double[][] data, boolean[][] marked, ColorScale scale, Color mark)
	{
		BufferedImage image = new BufferedImage(data.length, data[0].length, BufferedImage.TYPE_INT_RGB);
		if (scale == null) scale = ColorScales.getNew(data, 11);
		for (int i = 0; i < data.length; i++)
			for (int j = 0; j < data[0].length; j++)
				if (!marked[i][j])
					image.setRGB(i, j, scale.of(data[i][j]).getRGB());
				else
					image.setRGB(i, j, mark.getRGB());
		return image;
	}
	public static BufferedImage getBufferedImage(boolean[][] marked)
	{
		BufferedImage image = new BufferedImage(marked.length, marked[0].length, BufferedImage.TYPE_INT_RGB);
		for (int i = 0; i < marked.length; i++)
			for (int j = 0; j < marked[0].length; j++)
				if (!marked[i][j])
					image.setRGB(i, j, Color.BLACK.getRGB());
				else
					image.setRGB(i, j, Color.white.getRGB());
		return image;
	}
	public static void setTrueToColor(boolean[][] marked, BufferedImage image, Color mark)
	{
		for (int i = 0; i < marked.length; i++)
			for (int j = 0; j < marked[0].length; j++)
				if (marked[i][j])
					image.setRGB(i, j, mark.getRGB());
	}
	public static BufferedImage getBufferedImageTricolor(ArrayList<int[]> white, ArrayList<int[]> red, ArrayList<int[]> blue, int nx, int ny)
	{
		BufferedImage image = new BufferedImage(nx, ny, BufferedImage.TYPE_INT_RGB);
		ImageEditing.blacken(image);
		for (int i = 0; i < white.size(); i++)
			image.setRGB(white.get(i)[0], white.get(i)[1], Color.WHITE.getRGB());
		for (int i = 0; i < red.size(); i++)
			image.setRGB(red.get(i)[0], red.get(i)[1], Color.red.getRGB());
		for (int i = 0; i < blue.size(); i++)
			image.setRGB(blue.get(i)[0], blue.get(i)[1], Color.blue.getRGB());
		return image;
	}
	public static BufferedImage getMultiColorImage(ArrayList<int[]>[] markedPixels, Color[] colors, int nx, int ny)
	{
		BufferedImage image = new BufferedImage(nx, ny, BufferedImage.TYPE_INT_RGB);
		ImageEditing.blacken(image);
		for (int k = 0; k < markedPixels.length; k++)
			for (int i = 0; i < markedPixels[k].size(); i++)
				image.setRGB(markedPixels[k].get(i)[0], markedPixels[k].get(i)[1], colors[k].getRGB());
		return image;
	}
	/**
	 * Assumes complex is [N][M][2].
	 * @param complex
	 * @param scale
	 * @return
	 */
	public static BufferedImage getBufferedImage(double[][][] complex, ColorScale2d scale)
	{
		java.awt.Color c;
		BufferedImage out = new BufferedImage(complex.length, complex[0].length, BufferedImage.TYPE_INT_RGB);
		if (scale == null)
			scale = new ColorScales.MYC2d(FieldOps.magMax(complex), FieldOps.magMin(complex), 2*Math.PI);
		for (int i = 0; i < complex[0].length; i++)
			for (int j = 0; j < complex.length; j++)
			{
				c = scale.of(Complex.mag(complex[j][i]), FieldOps.atan(complex[j][i][0], complex[j][i][1]));
				out.setRGB(j, i, c.getRGB());
			}
		return out;
	}
	public static BufferedImage getBufferedImage(double[][] mag, double[][] phase, ColorScale2d scale)
	{
		if (scale == null)
			scale = new ColorScales.MYC2d(FieldOps.max(mag), FieldOps.min(mag), 2*Math.PI);
		java.awt.Color c;
		BufferedImage out = new BufferedImage(mag.length, mag[0].length, BufferedImage.TYPE_INT_RGB);
		int n = 0;
		for (int i = 0; i < mag[0].length; i++)
			for (int j = 0; j < mag.length; j++)
			{
				c = scale.of(mag[j][i], phase[j][i]);
				out.setRGB(j, i, c.getRGB());
			}
		return out;

	}
	public static BufferedImage getBufferedImageTranspose(double[][] data)
	{
		BufferedImage image = new BufferedImage(data[0].length, data.length, BufferedImage.TYPE_INT_RGB);
		ColorScales.LinearBRYW scale = new ColorScales.LinearBRYW(ArrayOps.max(data), ArrayOps.min(data));
		for (int i = 0; i < data[0].length; i++)
			for (int j = 0; j < data.length; j++)
				image.setRGB(i, j, scale.of(data[j][i]).getRGB());
		return image;

	}
	public static BufferedImage getBufferedImageTranspose(double[][] data, int scaleIndex)
	{
		BufferedImage image = new BufferedImage(data[0].length, data.length, BufferedImage.TYPE_INT_RGB);
		ColorScale scale = ColorScales.getNew(data, scaleIndex);
		for (int i = 0; i < data[0].length; i++)
			for (int j = 0; j < data.length; j++)
				image.setRGB(i, j, scale.of(data[j][i]).getRGB());
		return image;

	}
	
	//This will assign each pixel in the output image according to outputColor = inputColor + (outputColor-inputColor)*amount
	public static BufferedImage moveEachPixelTowardAColor(BufferedImage input, double[][] amount, Color color)
	{
		BufferedImage output = new BufferedImage(input.getWidth(), input.getHeight(), BufferedImage.TYPE_INT_RGB);
		moveEachPixelTowardAColor(input, amount, color, output);
		return output;
	}
	public static BufferedImage getSubset(BufferedImage input, int xm, int dx, int ym, int dy)
	{
		BufferedImage sub = new BufferedImage(dx, dy, BufferedImage.TYPE_INT_RGB);
		for (int i = 0; i < dx; i++)
			for (int j = 0; j < dy; j++)
				sub.setRGB(i, j, input.getRGB(xm+i, ym+j));
		return sub;
	}
	public static void moveEachPixelTowardAColor(BufferedImage input, double[][] amount, Color color, BufferedImage output)
	{
		int dr, dg, db;
		Color in, out;
		double a;
		for (int i = 0; i < amount.length; i++)
			for (int j = 0; j < amount[0].length; j++)
			{	
				a = amount[i][j];
				if (a >= 1.0/256)
				{
					in = new Color(input.getRGB(i, j));
					dr = color.getRed() - in.getRed();
					dg = color.getGreen() - in.getGreen();
					db = color.getBlue() - in.getBlue();
					out = new Color((int)(in.getRed() + dr*a), (int)(in.getGreen() + dg*a), (int)(in.getBlue() + db*a));
					output.setRGB(i, j, out.getRGB());
				}
				else
					output.setRGB(i, j, input.getRGB(i, j));
			}
	}
	
	//assumes that the source has even number dimensions
	public static void reduce2(BufferedImage source, BufferedImage target)
	{
		int r, g, b;
		Color c;
		for (int i = 0; i < source.getWidth()-2; i+=2)
			for (int j = 0; j < source.getHeight()-2; j+=2)
			{
				r = 0; g = 0; b = 0;
				c = new Color(source.getRGB(i, j));
				r += c.getRed();
				g += c.getGreen();
				b += c.getBlue();
				c = new Color(source.getRGB(i+1, j));
				r += c.getRed();
				g += c.getGreen();
				b += c.getBlue();
				c = new Color(source.getRGB(i, j+1));
				r += c.getRed();
				g += c.getGreen();
				b += c.getBlue();
				c = new Color(source.getRGB(i+1, j+1));
				r += c.getRed();
				g += c.getGreen();
				b += c.getBlue();
				c = new Color(r/4, g/4, b/4);
				target.setRGB(i/2, j/2, c.getRGB());
			}
	}
	public static void reduce(BufferedImage source, BufferedImage target, int fact)
	{
		int r, g, b;
		Color c;
		for (int i = 0; i < source.getWidth()-fact; i+=fact)
			for (int j = 0; j < source.getHeight()-fact; j+=fact)
			{
				
				r = 0; g = 0; b = 0;
				for (int m = 0; m < fact; m++)
					for (int n = 0; n < fact; n++)
					{
						c = new Color(source.getRGB(i+m, j+n));
						r += c.getRed();
						g += c.getGreen();
						b += c.getBlue();

					}
				c = new Color(r/(fact*fact), g/(fact*fact), b/(fact*fact));
				target.setRGB(i/fact, j/fact, c.getRGB());
			}
	}
	
	//This produces target as target = p0*(1-p) + p1*p, where p is between 0 and 1.
	public static void fuse(BufferedImage p0, BufferedImage p1, double p, BufferedImage target)
	{
//		Color c0, c1, c;
		int[] c0 = new int [3], c1 = new int[3];
		int r, g, b;
		if (p < 1/256.0) {
			copy(p0, target);
			return;
		}
		else if (p > 255.0/256.0)
		{
			copy(p1, target);
			return;
		}
		for (int i = 0; i < p0.getWidth(); i++)
			for (int j = 0; j < p1.getHeight(); j++){
//				c0 = new Color(p0.getRGB(i, j));
//				c1 = new Color(p1.getRGB(i, j));
				putColor(p0.getRGB(i, j), c0);
				putColor(p1.getRGB(i, j), c1);
				r = (int)(p*c1[0] + (1-p)*c0[0]);
				g = (int)(p*c1[1] + (1-p)*c0[1]);
				b = (int)(p*c1[2] + (1-p)*c0[2]);
//				c = new Color(r, g, b);
				target.setRGB(i, j, RGB(r, g, b));
			}
	}
	
	//Creates a combined image of two images horizontal to each other
	public static BufferedImage weldHorizontal(BufferedImage a, BufferedImage b)
	{
		BufferedImage ans = new BufferedImage(a.getWidth() + b.getWidth(), Math.max(a.getHeight(), b.getHeight()), BufferedImage.TYPE_INT_RGB);
		for (int i = 0; i < a.getWidth(); i++)
			for (int j = 0; j < a.getHeight(); j++)
				ans.setRGB(i, j, a.getRGB(i,j));
		for (int i = 0; i < b.getWidth(); i++)
			for (int j = 0; j < b.getHeight(); j++)
				ans.setRGB(i+a.getWidth(), j, b.getRGB(i,j));
		return ans;
	}
	//assumes that the source has even number dimensions
	public static void enlargeBasic(BufferedImage source, BufferedImage target, int factor)
	{
		for (int i = 0; i < source.getWidth(); i++)
			for (int j = 0; j < source.getHeight(); j++)
			{
				for (int m = 0; m < factor; m++)
					for (int n = 0; n < factor; n++)
						target.setRGB(factor*i+m, factor*j+n, source.getRGB(i, j));
			}
	}
	public static BufferedImage getEnlarged(BufferedImage source, int factor)
	{
		BufferedImage target = new BufferedImage (source.getWidth()*factor, source.getHeight()*factor, BufferedImage.TYPE_INT_RGB);
		for (int i = 0; i < source.getWidth(); i++)
			for (int j = 0; j < source.getHeight(); j++)
			{
				for (int m = 0; m < factor; m++)
					for (int n = 0; n < factor; n++)
						target.setRGB(factor*i+m, factor*j+n, source.getRGB(i, j));
			}
		return target;
	}
	public static BufferedImage enlargeBasicStretch(BufferedImage source, int fx, int fy)
	{
		BufferedImage target = new BufferedImage(fx*source.getWidth(), fy*source.getHeight(), source.getType());
		enlargeBasicStretch(source, target, fx, fy);
		return target;
	}
	public static void enlargeBasicStretch(BufferedImage source, BufferedImage target, int fx, int fy)
	{
		for (int i = 0; i < source.getWidth(); i++)
			for (int j = 0; j < source.getHeight(); j++)
			{
				for (int m = 0; m < fx; m++)
					for (int n = 0; n < fy; n++)
						target.setRGB(fx*i+m, fy*j+n, source.getRGB(i, j));
			}
	}
	
	public static void copy(BufferedImage source, BufferedImage target)
	{
		for (int i = 0; i < source.getWidth(); i++)
			for (int j = 0; j < source.getHeight(); j++)
				target.setRGB(i, j, source.getRGB(i, j));
	}
	public static void copyInto(BufferedImage source, BufferedImage target, int x0, int y0)
	{
		for (int i = 0; i < source.getWidth(); i++)
			for (int j = 0; j < source.getHeight(); j++)
				target.setRGB(i+x0, j+y0, source.getRGB(i, j));
	}
	
	
	private static void putColor(int rgb, int[] color)
	{
        color[2] = ((rgb>>16)&0xff);
        color[1] = ((rgb>>8)&0xff);
        color[0] = (rgb&0xff);
	}
	private static int RGB(int r, int g, int b)
	{
		return   ((r&0xff) | ((g&0xff)<<8)
                | ((b&0xff)<<16) | ((Integer.MAX_VALUE)<<24));

	}
	private static int RGB(int[] c)
	{
		return   ((c[0]&0xff) | ((c[1]&0xff)<<8)
                | ((c[2]&0xff)<<16) | ((Integer.MAX_VALUE)<<24));
	}
	public static void writeLine(BufferedImage image, String line, int x, int y, Color c)
	{
		java.awt.Graphics g = image.getGraphics();
		g.setColor(c);
		Font f = new Font(Font.SERIF, Font.PLAIN, 20);
		g.setFont(f);
		g.drawString(line, x, y);
	}
	public static void writeLine(BufferedImage image, String line, int x, int y, Color c, String fonttype, int fontstyle, int fontsize)
	{
		java.awt.Graphics g = image.getGraphics();
		g.setColor(c);
		Font f = new Font(fonttype, fontstyle, fontsize);
		g.setFont(f);
		g.drawString(line, x, y);
	}
	public static void drawLine(BufferedImage image, int x1, int y1, int x2, int y2, Color c)
	{
		java.awt.Graphics g = image.getGraphics();
		g.setColor(c);
		g.drawLine(x1, y1, x2, y2);
	}
	public static void drawPlus(BufferedImage image, int x, int y, int halflength, Color c)
	{
		java.awt.Graphics g = image.getGraphics();
		g.setColor(c);
		g.drawLine(x-halflength, y, x+halflength, y);
		g.drawLine(x, y-halflength, x, y+halflength);
	}
	public static void drawThickLine(BufferedImage image, int x1, int y1, int x2, int y2, Color c, int thickness)
	{
		java.awt.Graphics g = image.getGraphics();
		g.setColor(c);
		
		g.drawLine(x1, y1, x2, y2);
		for (int i = 0; i < thickness; i++)
		{
			g.drawLine(x1+i, y1+i, x2-i, y2-i);
			g.drawLine(x1-i, y1-i, x2-i, y2-i);
			g.drawLine(x1-i, y1-i, x2+i, y2+i);
			g.drawLine(x1+i, y1+i, x2+i, y2+i);
		}
	}
	public static void flipX(BufferedImage image) {
		int[][] y = new int[image.getWidth()][image.getHeight()];
		for (int i = 0; i <y.length; i++)
			for (int j = 0; j < y[i].length; j++)
				y[i][j] = image.getRGB(y.length-1-i,j);
		for (int i = 0; i < y.length; i++)
			for (int j = 0; j < y[i].length; j++)
				image.setRGB(i, j, y[i][j]);

		
	}
	public static void blacken(BufferedImage ans) {
		for (int i = 0; i < ans.getWidth(); i++)
			for (int j = 0; j < ans.getHeight(); j++)
				ans.setRGB(i, j, Color.black.getRGB());
	}
	public static void setColor(BufferedImage ans, Color c) {
		for (int i = 0; i < ans.getWidth(); i++)
			for (int j = 0; j < ans.getHeight(); j++)
				ans.setRGB(i, j, c.getRGB());
	}
	
	public static BufferedImage open(String filepath)
	{
		BufferedImage got = null;
		try {
			got = ImageIO.read(new java.io.File(filepath));
			return got;
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}
	public static boolean[][] loadFromImage(String path, Color trueColor)
	{
		BufferedImage image = ImageEditing.open(path);
		boolean[][] ans = new boolean [image.getWidth()][image.getHeight()];
		for (int i = 0; i < ans.length; i++)
			for (int j = 0; j < ans[0].length; j++)
				ans[i][j] = new Color(image.getRGB(i, j)).equals(trueColor);
		return ans;
	}
	public static double[][] getPhase(BufferedImage myc){
		double[][] ans = new double [myc.getWidth()][myc.getHeight()];
		Color c;
		int r, g, b;
		int max;
		for (int i = 0; i < ans.length; i++)
			for (int j = 0; j < ans[0].length; j++)
			{
				c = new Color(myc.getRGB(i, j));
				r = c.getRed();
				g = c.getGreen();
				b = c.getBlue();
				max = Math.max(Math.max(r, g), b);
				if (max == r) //first third: magenta to yellow
					ans[i][j] = 0 + (g/(double)max)*2*Math.PI/3;
				if (max == g) //yellow to cyan
					ans[i][j] = 2*Math.PI/3 + (b/(double)max)*2*Math.PI/3;
				if (max == b)
					ans[i][j] = 4*Math.PI/3 + (r/(double)max)*2*Math.PI/3;
			}
		
	return ans;
		
	}
	public static double[][] getFromImage(BufferedImage im){
		double[][] ans = new double [im.getWidth()][im.getHeight()];
		Color c;
		int r, g, b;
		for (int i = 0; i < ans.length; i++)
			for (int j = 0; j < ans[0].length; j++)
			{
				c = new Color(im.getRGB(i, j));
				r = c.getRed();
				g = c.getGreen();
				b = c.getBlue();
				ans[i][j] = r+g+b;
			}
		System.out.println("loaded image data");
	return ans;
		
	}
	
	/**
	 * The magnitude will range from 0 to 1 
	 * @param myc
	 * @return
	 */
	public static double[][][] getComplex(BufferedImage myc){
		double[][][] ans = new double [myc.getWidth()][myc.getHeight()][2];
		double phase=0;
		Color c;
		int r, g, b;
		int max;
		for (int i = 0; i < ans.length; i++)
			for (int j = 0; j < ans[0].length; j++)
			{
				c = new Color(myc.getRGB(i, j));
				r = c.getRed();
				g = c.getGreen();
				b = c.getBlue();
				max = Math.max(Math.max(r, g), b);
				if (max == r) //first third: magenta to yellow
					phase = 0 + (g/(double)max)*2*Math.PI/3;
				if (max == g) //yellow to cyan
					phase = 2*Math.PI/3 + (b/(double)max)*2*Math.PI/3;
				if (max == b)
					phase = 4*Math.PI/3 + (r/(double)max)*2*Math.PI/3;
				ans[i][j][0] = (max/255.0)*Math.cos(phase);
				ans[i][j][1] = (max/255.0)*Math.sin(phase);
			}
		
	return ans;
		
	}
	public static BufferedImage createTileImage(double[][][] ds, ColorScale scale, int nPerLine, double[] v, int scaleIndex) {
		
		double nlines = ((double)ds.length/nPerLine);
		int ny;
		if (nlines != (int)nlines) ny = (int)nlines+1;
		else ny = (int) nlines;
		
		int remainder = ny*nPerLine - ds.length;
		
		int tnx = ds[0].length, tny = ds[0][0].length;
		
		BufferedImage image = new BufferedImage(ds[0].length*nPerLine, ds[0][0].length*ny, BufferedImage.TYPE_INT_RGB);
		ImageEditing.setColor(image, Color.WHITE);
		
		int x0, y0, index;
		BufferedImage temp; Graphics g;
		boolean noscale = scale == null;
		for (int i = 0; i < ds.length; i++)
		{
			if (noscale)
				scale = ColorScales.getNew(ds[i], scaleIndex);
			index = remainder + i;
			x0 = ds[0].length*(index % nPerLine);
			y0 = ds[0][0].length*(index/nPerLine);
			temp = ImageEditing.getBufferedImage(ds[i], scale);
			g = temp.getGraphics();
			g.setColor(ColorScales.getUnusedColor(scale));
			g.setFont(new Font("Arial", Font.PLAIN, tnx/10));
			g.drawString(NumFormat.voltage(v[i]), 5, tny-5);
			ImageEditing.copyInto(temp, image, x0, y0);
		}
		return image;
	}
	
	public static BufferedImage createTileImage(BufferedImage[] parts, int nPerLine, Color background)
	{
		
		double nlines = ((double)parts.length/nPerLine);
		int ny;
		if (nlines != (int)nlines) ny = (int)nlines+1;
		else ny = (int) nlines;
		
		int remainder = ny*nPerLine - parts.length;
		int x0, y0, index;
		BufferedImage image = new BufferedImage(parts[0].getWidth()*nPerLine, parts[0].getHeight()*ny, BufferedImage.TYPE_INT_RGB);
		ImageEditing.setColor(image, background);
		for (int i = 0; i < parts.length; i++)
		{
			
			index = remainder + i;
			x0 = parts[0].getWidth()*(index % nPerLine);
			y0 = parts[0].getHeight()*(index/nPerLine);
			ImageEditing.copyInto(parts[i], image, x0, y0);
		}
		return image;
	}
	public static void makeGIF(String basename, JFileChooser fc, double[] kCalcFinal, boolean[] spectraSkipped)
	{
		Topomap t = Topomap.open(fc);
		AtomicCoordinatesSet latt = new AtomicCoordinatesSet(FileOps.openText(fc));
		double[] onDir = latt.getA();
		double braggAngleDeg = Math.toDegrees(latt.getAngleBetweenVectors());
		double[] offDir = Matrix.getProductWith(Matrix.getRotationMatrix(Math.toRadians(braggAngleDeg/2)), onDir);
		
		BufferedImage[] gifstack = new BufferedImage[t.nlayers];

		double[] ondR = new double [2], offdR = new double [2];
		ondR[0] = onDir[0]/Complex.mag(onDir);
		ondR[1] = onDir[1]/Complex.mag(onDir);
		offdR[0] = offDir[0]/Complex.mag(offDir);
		offdR[1] = offDir[1]/Complex.mag(offDir);
		double[] dr;
		int ox = 0, oy = 0;
		int choice = JOptionPane.showConfirmDialog(null, "Is this the on-direction, or off-direction linecut? (Yes=On, Rt2=cancel)");
		if (choice == JOptionPane.YES_OPTION)
			dr = ondR;
		else if (choice == JOptionPane.NO_OPTION)
			dr = offdR;
		else{
			int[][] bragg = AtomicCoordinatesSet.generateBragg(latt.getRt2Lattice(), t.nx);
			dr = new double [] {bragg[1][0] - bragg[0][0], bragg[1][1]-bragg[0][1]};
			dr = Distance.unitVector(dr);
			ox = bragg[0][0];
			oy = bragg[0][1];
		}
		
		double kmag;
		
		
		Graphics g;
		for (int i = 0; i < t.nlayers; i++)
		{
			gifstack[i] = ImageEditing.getBufferedImage(t.data[i]);
			g = gifstack[i].getGraphics();
			g.setColor(java.awt.Color.BLUE);
			g.drawString(NumFormat.voltage(t.v[i]), 5, t.ny-5);
			if (!spectraSkipped[i]){
				kmag = kCalcFinal[i]*t.nx/t.xLength; //Assuming that t is a symmetrized FFT, the length scale is 1/m.
				double[] r = new double[] {t.nx/2 + ox + kmag*dr[0], t.ny/2 + oy + kmag*dr[1]};
				g.setColor(Color.BLUE);
				drawPlus(g, FieldOps.round(r[0]), FieldOps.round(r[1]), 3);
			}
			
		}
		GifSequenceWriter.writeGifSequence(gifstack, 500, basename + "Modes Movie.gif");
	}
	public static void drawPlus(Graphics g, int x, int y, int halflength)
	{
		g.drawLine(x-halflength, y, x+halflength, y);
		g.drawLine(x, y-halflength, x, y+halflength);
	}

	/**
	 * This is to create an image with an appropriate-sized label for the bias voltage.
	 * label codes: 0 = conventional label
	 * 				1 = Integer label index: (i+1)
	 * 				2 = No label.
	 * @param ds
	 * @param scaleindex
	 * @param v
	 * @return
	 */
	public static BufferedImage[] createImageArray(double[][][] ds, int scaleindex, double[] v, int label, double minperc, double maxperc, int resizeFactor, ColorScale scale) {
		
		BufferedImage[] output = new BufferedImage [v.length];
		BufferedImage temp;
		Graphics g;
		double[] bounds;
		int tnx = ds[0].length*resizeFactor, tny = ds[0][0].length*resizeFactor;
		for (int i = 0; i < ds.length; i++)
		{
			bounds = FieldOps.getPercentiles(ds[i], minperc, maxperc);
			output[i] = new BufferedImage(tnx, tny, BufferedImage.TYPE_INT_RGB);
			temp = ImageEditing.getBufferedImage(ds[i], scale == null ? ColorScales.getNew(bounds[0], bounds[1], scaleindex) : scale);
			enlargeBasic(temp, output[i], resizeFactor);
			g = output[i].getGraphics();
			g.setColor(ColorScales.getUnusedColor(scaleindex));
			g.setFont(new Font("Arial", Font.PLAIN, (int) Math.sqrt(tnx)));
			if (label == 0)
				g.drawString(NumFormat.voltage(v[i]), 5, tny-5);
			else if (label == 1)
				g.drawString("" + (int)(i+1),5, tny-5);
			else;
		}
		return output;
	}
	public static BufferedImage[] createImageArray(double[][][] ds, int scaleindex, double[] v, int label, double minperc, double maxperc, int[] resizeFactor, ColorScale scale) {
		
		BufferedImage[] output = new BufferedImage [v.length];
		BufferedImage temp;
		Graphics g;
		double[] bounds;
		int tnx = ds[0].length*resizeFactor[0], tny = ds[0][0].length*resizeFactor[1];
		for (int i = 0; i < ds.length; i++)
		{
			bounds = FieldOps.getPercentiles(ds[i], minperc, maxperc);
			output[i] = new BufferedImage(tnx, tny, BufferedImage.TYPE_INT_RGB);
			temp = ImageEditing.getBufferedImage(ds[i], scale == null ? ColorScales.getNew(bounds[0], bounds[1], scaleindex) : scale);
			enlargeBasicStretch(temp, output[i], resizeFactor[0], resizeFactor[1]);
			g = output[i].getGraphics();
			g.setColor(ColorScales.getUnusedColor(scaleindex));
			g.setFont(new Font("Arial", Font.PLAIN, (int) Math.sqrt(tnx)));
			if (label == 0)
				g.drawString(NumFormat.voltage(v[i]), 5, tny-5);
			else if (label == 1)
				g.drawString("" + (int)(i+1),5, tny-5);
			else;
		}
		return output;
	}
	/**
	 * Returns the ith image that would have been returned by createImageArray
	 * @param ds
	 * @param scaleindex
	 * @param v
	 * @param label
	 * @param minperc
	 * @param maxperc
	 * @param resizeFactor
	 * @param scale
	 * @param i
	 * @return
	 */
	public static BufferedImage createSingleImage(double[][][] ds, int scaleindex, double[] v, int label, double minperc, double maxperc, int resizeFactor, ColorScale scale, int i) {
		
		BufferedImage output;
		BufferedImage temp;
		Graphics g;
		double[] bounds;
		int tnx = ds[0].length*resizeFactor, tny = ds[0][0].length*resizeFactor;
		bounds = FieldOps.getPercentiles(ds[i], minperc, maxperc);
		output = new BufferedImage(tnx, tny, BufferedImage.TYPE_INT_RGB);
		temp = ImageEditing.getBufferedImage(ds[i], scale == null ? ColorScales.getNew(bounds[0], bounds[1], scaleindex) : scale);
		enlargeBasic(temp, output, resizeFactor);
		g = output.getGraphics();
		g.setColor(ColorScales.getUnusedColor(scaleindex));
		g.setFont(new Font("Arial", Font.PLAIN, (int) Math.sqrt(tnx)));
		if (label == 0)
			g.drawString(NumFormat.voltage(v[i]), 5, tny-5);
		else if (label == 1)
			g.drawString("" + (int)(i+1),5, tny-5);
		else;

		return output;
	}
	//Here I assume the image is a centered Fourier transform.
	public static void drawBrillouinZone(BufferedImage kspace, AtomicCoordinatesSet latt, Color lineColor)
	{
		boolean squareLattice = latt.getAngleBetweenVectors() > 5*Math.PI/12 && latt.getAngleBetweenVectors() <= 7*Math.PI/12;
		int ox = kspace.getWidth()/2;
		int oy = kspace.getHeight()/2;
		int nx = kspace.getWidth();
		int ny = kspace.getHeight();
		int[][] bragg = AtomicCoordinatesSet.generateBragg(latt, nx);
		Graphics g = kspace.getGraphics();
		g.setColor(lineColor);
		int p1x, p2x, p1y, p2y;
		if (squareLattice)
		{
			//This position being PI, PI;
			p1x = (bragg[0][0]+bragg[1][0])/2; p1y = (bragg[0][1]+bragg[1][1])/2;
			p2x = (bragg[0][0]-bragg[1][0])/2; p2y = (bragg[0][1]-bragg[1][1])/2; //This PI, -PI;
			g.drawLine(ox+p1x, oy+p1y, ox+p2x, oy+p2y);
			p1x = p2x; p1y = p2y;
			p2x = (-bragg[0][0]-bragg[1][0])/2; p2y = (-bragg[0][1]-bragg[1][1])/2; //This PI, -PI;
			g.drawLine(ox+p1x, oy+p1y, ox+p2x, oy+p2y);
			p1x = p2x; p1y = p2y;
			p2x = (-bragg[0][0]+bragg[1][0])/2; p2y = (-bragg[0][1]+bragg[1][1])/2; //This PI, -PI;
			g.drawLine(ox+p1x, oy+p1y, ox+p2x, oy+p2y);
			p1x = p2x; p1y = p2y;
			p2x = (+bragg[0][0]+bragg[1][0])/2; p2y = (+bragg[0][1]+bragg[1][1])/2; //This PI, -PI;
			g.drawLine(ox+p1x, oy+p1y, ox+p2x, oy+p2y);
		}
		else{
			//Now we have to draw six lines dammit... we require 3 bragg peaks and to wrap around in the manner of the previous
			//thing, we must ensure that they are neighborly ordered...
			int[] bragg3 = new int [2];
			int[] bragg1 = new int [2];
			int[] bragg2 = new int [2];
			//If the existing bragg peaks are at 60 degrees the third peak is simply their difference...
			System.out.println(latt.getAngleBetweenVectors());
			if (Distance.getAngleBetween(bragg[0], bragg[1]) < 5*Math.PI/12){
				bragg3[0] = bragg[1][0]-bragg[0][0];
				bragg3[1] = bragg[1][1]-bragg[0][1];
				bragg1 = bragg[0];
				bragg2 = bragg[1];
				//In neighborly order the order must be 1, 2, 3. Therefore if bragg3 is closer to bragg1, it must be moved:
				if (Distance.getAngleBetween(bragg3, bragg2) > Distance.getAngleBetween(bragg3, bragg1))
				{
					bragg3[0] = -bragg3[0];
					bragg3[1] = -bragg3[1];
				}
			}
			//If on the other hand the original vectors are at 120 degrees, it is their sum which is in between them:
			else if (Distance.getAngleBetween(bragg[0], bragg[1]) > 7*Math.PI/12){
				bragg3[0] = bragg[1][0]+bragg[0][0];
				bragg3[1] = bragg[1][1]+bragg[0][1];
				bragg1 = bragg[0];
				bragg2 = bragg[1];
				//HOwever at present 3 lies in between 1 and 2. We must switch 2 and 3:
				int[] temp = bragg3.clone();
				bragg3 = bragg2.clone();
				bragg2 = temp.clone();
			}
			//Now we proceed around the hexagon:
			p1x = (bragg1[0]+bragg2[0])/3; p1y = (bragg1[1]+bragg2[1])/3;
			p2x = (bragg2[0]+bragg3[0])/3; p2y = (bragg2[1]+bragg3[1])/3;
			g.drawLine(ox+p1x, oy+p1y, ox+p2x, oy+p2y);
			p1x = p2x; p1y = p2y;
			p2x = (bragg3[0]-bragg1[0])/3; p2y = (bragg3[1]-bragg1[1])/3;
			g.drawLine(ox+p1x, oy+p1y, ox+p2x, oy+p2y);
			p1x = p2x; p1y = p2y;
			p2x = (-bragg2[0]-bragg1[0])/3; p2y = (-bragg2[1]-bragg1[1])/3;
			g.drawLine(ox+p1x, oy+p1y, ox+p2x, oy+p2y);
			p1x = p2x; p1y = p2y;
			p2x = (-bragg2[0]-bragg3[0])/3; p2y = (-bragg2[1]-bragg3[1])/3;
			g.drawLine(ox+p1x, oy+p1y, ox+p2x, oy+p2y);
			p1x = p2x; p1y = p2y;
			p2x = (-bragg3[0]+bragg1[0])/3; p2y = (-bragg3[1]+bragg1[1])/3;
			g.drawLine(ox+p1x, oy+p1y, ox+p2x, oy+p2y);
			p1x = p2x; p1y = p2y;
			p2x = (bragg2[0]+bragg1[0])/3; p2y = (bragg2[1]+bragg1[1])/3;
			g.drawLine(ox+p1x, oy+p1y, ox+p2x, oy+p2y);
		}
		
	}
	public static void writeBiasVoltage(BufferedImage image, int scaleindex, int label, double v)
	{
		int tnx = image.getWidth(), tny = image.getHeight();
		Graphics g = image.getGraphics();
		g.setColor(ColorScales.getUnusedColor(scaleindex));
		g.setFont(new Font("Arial", Font.PLAIN, (int) Math.sqrt(tnx)));
		if (label == 0)
			g.drawString(NumFormat.voltage(v), 5, tny-5);
		else;
	
	}
	public static BufferedImage getCopy(BufferedImage source) {
		BufferedImage target = new BufferedImage (source.getWidth(), source.getHeight(), BufferedImage.TYPE_INT_RGB);
		for (int i = 0; i < source.getWidth(); i++)
			for (int j = 0; j < source.getHeight(); j++)
				target.setRGB(i, j, source.getRGB(i, j));
		return target;
	}
	
	public static class ImageTransferable implements Transferable
    {
        private Image image;

        public ImageTransferable (Image image)
        {
            this.image = image;
        }

        public Object getTransferData(DataFlavor flavor)
            throws UnsupportedFlavorException
        {
            if (isDataFlavorSupported(flavor))
            {
                return image;
            }
            else
            {
                throw new UnsupportedFlavorException(flavor);
            }
        }

        public boolean isDataFlavorSupported (DataFlavor flavor)
        {
            return flavor == DataFlavor.imageFlavor;
        }

        public DataFlavor[] getTransferDataFlavors ()
        {
            return new DataFlavor[] { DataFlavor.imageFlavor };
        }
    }
	
	public static void copyToClipboard(BufferedImage im)
	{
        ImageTransferable transferable = new ImageTransferable(im);
        Toolkit.getDefaultToolkit().getSystemClipboard().setContents(transferable, null);
	}
	public static BufferedImage getImpurityModuloMap(double[][] impLattPos,
			PointImp[] imps, int size, int mapSize) {
		BufferedImage image = new BufferedImage(size, size, BufferedImage.TYPE_INT_RGB);
		setColor(image, Color.GREEN);
		Graphics g = image.getGraphics();
		int r, b;
		int radius = (int)(size/(2*Math.sqrt(imps.length)));
		for (int i = 0; i < imps.length; i++)
		{
			r = (int)(imps[i].pixelPos[0]*256/mapSize);
			b = (int)(imps[i].pixelPos[1]*256/mapSize);
			r = Math.max(r, 0);
			r = Math.min(r, 255);
			b = Math.max(b, 0);
			b = Math.min(b, 255);
			
			g.setColor(new Color(r, 0, b));
			g.fillOval((int)(impLattPos[i][0]*size)-radius, (int)(impLattPos[i][1]*size)-radius, 2*radius, 2*radius);
		}
		return image;
	}
	public static BufferedImage[] getImpurityModuloMapMovie(double[][] impLattPos,
			PointImp[] imps, Topomap t, int[] maxi) {
		
		BufferedImage[] ans = new BufferedImage [t.nlayers];
		int fontsize = (int)Math.sqrt(t.nx);
		Graphics g;
		int r, b;
		int radius = (int)(t.nx/(2*Math.sqrt(imps.length)));
		for (int j = 0; j < t.nlayers; j++){
			ans[j] = new BufferedImage(t.nx, t.ny+(int)Math.sqrt(2*fontsize*fontsize), BufferedImage.TYPE_INT_RGB);
			setColor(ans[j], Color.GREEN);
			g = ans[j].getGraphics();
			
			for (int i = 0; i < imps.length; i++)
			{
				r = (int)(imps[i].pixelPos[0]*256/t.nx);
				b = (int)(imps[i].pixelPos[1]*256/t.ny);
				r = Math.max(r, 0);
				r = Math.min(r, 255);
				b = Math.max(b, 0);
				b = Math.min(b, 255);
				if (maxi[i] == j){
					g.setColor(new Color(r, 0, b));
					g.fillOval((int)(impLattPos[i][0]*t.nx)-radius, (int)(impLattPos[i][1]*t.ny)-radius, 2*radius, 2*radius);
				}
				g.setFont(new Font("Arial", Font.PLAIN, (int) Math.sqrt(t.nx)));
				g.setColor(Color.BLACK);
				g.drawString(NumFormat.voltage(t.v[j]), 5, t.ny+fontsize);
				g.drawLine(0, t.ny-1, t.nx, t.ny-1);
			}
		}
		return ans;
	}
	public static void writeForAVIMovie(BufferedImage[] stack, String basename)
	{
		for (int i = 0; i < stack.length; i++)
			SRAW.writeImage(MovieMaker.avidir + basename + MovieMaker.fromInt(i), stack[i]);
		System.out.println(MovieMaker.BMPtoAVICommand(basename, 0, stack.length-1, 5));

	}
	public static BufferedImage getBufferedImageBicubic(double[][] data, double[] x, double [] y, ColorScale cs, int csi, BicubicSplineInterpolatingFunction interp)
	{
		int nx = x.length, ny = y.length;
		double mean = FieldOps.mean(data);
		BufferedImage ans = new BufferedImage(nx, nx, BufferedImage.TYPE_INT_RGB);
		double[][] localData = new double [nx][ny];

		for (int i = 0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				if ((x[i] < 0 || x[i] > data.length-1) || (y[j] < 0 || y[j] > data[0].length-1))
					localData[i][j] = mean;
				else
					localData[i][j] = interp.value(x[i], y[j]);
			}
		return ImageEditing.getBufferedImage(localData, (cs == null ? ColorScales.getNew(localData, csi) : cs)  );
			
	}
	public static void color_ARGB(BufferedImage ks, Color color) {
		for (int i = 0; i < ks.getWidth(); i++)
			for (int j = 0; j < ks.getHeight(); j++)
			{
				ks.setRGB(i, j, color.getRGB());
			}
	}
	public static void renderTransparent(BufferedImage im){
		color_ARGB(im, new java.awt.Color(0, 0, 0, 0));
	}

}
