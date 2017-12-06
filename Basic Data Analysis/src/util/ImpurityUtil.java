package util;

import image.ImageEditing;
import impurity.PointImp;

import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import main.SRAW;

import org.apache.commons.math3.analysis.interpolation.BicubicSplineInterpolatingFunction;
import org.apache.commons.math3.analysis.interpolation.BicubicSplineInterpolator;

import schrodinger.MovieMaker;
import util.fileops.FileOps;
import util.fileops.Layer;
import util.fileops.PointSpectra;
import util.fileops.Topomap;
import util.geom.AtomicCoordinatesSet;
import util.geom.Distance;

public class ImpurityUtil {

	public static void main(String[] args)
	{
		Topomap.setStdDir();
		JFileChooser fc = new JFileChooser(Topomap.stddir);
//		Layer t = Layer.open(fc);
		AtomicCoordinatesSet latt = new AtomicCoordinatesSet(FileOps.openText(fc));
		PointImp[] imps = PointImp.readFromGaussSquareFile(FileOps.selectOpen(fc));
//		boolean[] nema = splitImpsWRTLayerValue(imps, latt, t, 0);
//		PointImp[][] splitNema = splitImps(imps, nema);
		boolean[] rt = splitImpsWRTRoot2(imps, latt, 0.5, 0.5);
		PointImp[][] splitNema = splitImps(imps, rt);
		File f = FileOps.selectSave(fc);
		PointImp.writeToFile(splitNema[0], new File(f.toString() + "_A.txt"));
		PointImp.writeToFile(splitNema[1], new File(f.toString() + "_B.txt"));
		
	}
	public static double[][] getMaskingLayer(Layer t, double gaussradius, JFileChooser fc)
	{
		PointImp[] imps = PointImp.readFromGaussSquareFile(FileOps.selectOpen(fc));
		return getMaskingLayer(t, gaussradius, imps);
	}
	public static double[][] getMaskingLayer(Layer t, double gaussradius, PointImp[] imps)
	{
//		PointImp[] imps = PointImp.readFromGaussSquareFile(FileOps.selectOpen(fc));
		double[][] dist = getDistanceTo(imps, t.nx, t.ny);
		double[][] map = new double [t.nx][t.ny];
		for (int i = 0; i < t.nx; i++)
			for (int j = 0; j < t.ny; j++)
				map[i][j] = Math.exp(-dist[i][j]*dist[i][j]/(gaussradius*gaussradius));
		
		boolean onEqualsOne = JOptionPane.showConfirmDialog(null, "Backround = 0?") == JOptionPane.YES_OPTION;
		if (!onEqualsOne){
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
					map[i][j] = 1 - map[i][j];
		}
		return map;
	}
	
	/**
	 * Slow, stupid brute-force calculation of the distance of each pixel between 0 and nx-1, ny-1 to the nearest impurity in the list.
	 * @param imps
	 * @param nx
	 * @param ny
	 * @return
	 */
	public static double[][] getDistanceTo(PointImp[] imps, int nx, int ny)
	{
		double[][] ans = new double [nx][ny];
		double[] distances = new double[imps.length];
		for (int i = 0 ; i < nx; i++){System.out.print(" " + i );
			for (int j = 0; j < ny; j++)
				for (int k = 0; k < imps.length; k++)
				{
					distances[k] = Distance.distance(new int [] {i, j}, imps[k].pixelPos);
					ans[i][j] = ArrayOps.min(distances);
				}
		}
		return ans;
	}
	public static double[][] getHeightOneGaussians(PointImp[] imps, int nx, int ny, double r)
	{
		double[][] dist = getDistanceTo(imps, nx, ny);
		double[][] ans = new double [nx][ny];
		for (int i = 0 ; i < nx; i++){
			for (int j = 0; j < ny; j++)
				ans[i][j] = Math.exp(-dist[i][j]/r*r);
		}
		return ans;
	}
	public static double[][] getOneOverR_Power_FromAll(PointImp[] imps, int nx, int ny, double epsilon, int power)
	{
		double[][] ans = new double [nx][ny];
		double[] distances = new double[imps.length];
		double max = 1/Math.pow(epsilon, power);
		for (int i = 0 ; i < nx; i++){System.out.print(" " + i );
			for (int j = 0; j < ny; j++)
				for (int k = 0; k < imps.length; k++)
				{
					distances[k] = Distance.distance(new int [] {i, j}, imps[k].pixelPos);
					
					ans[i][j] -= distances[k] > epsilon ? 1/Math.pow(distances[k], power) : max;
				}
		}
		return ans;
	}
	public static double[][] getYukawa_Power_FromAll(PointImp[] imps, int nx, int ny, double epsilon, int power, double explength)
	{
		double[][] ans = new double [nx][ny];
		double[] distances = new double[imps.length];
		double max = 1/Math.pow(epsilon, power);
		for (int i = 0 ; i < nx; i++){System.out.print(" " + i );
			for (int j = 0; j < ny; j++)
				for (int k = 0; k < imps.length; k++)
				{
					distances[k] = Distance.distance(new int [] {i, j}, imps[k].pixelPos);
					
					ans[i][j] -= distances[k] > epsilon ? Math.exp(-distances[k]/explength)/Math.pow(distances[k], power) : max;
				}
		}
		return ans;
	}
	
	public static double[][] countImpuritiesDensity(PointImp[] imps, int nx, int ny, double radius)
	{
		double[][] ans = new double [nx][ny];
		double[] distances = new double[imps.length];
		int pmin = FieldOps.roundDown(-radius-1);
		int pmax = FieldOps.round(radius+1);
		int n;
		for (int i = 0 ; i < nx; i++){System.out.print(" " + i );
			for (int j = 0; j < ny; j++){
				n = 0;
				for (int k = 0; k < imps.length; k++){
					distances[k] = Distance.distance(new int [] {i, j}, imps[k].pixelPos);
					if (distances[k] < radius) ans[i][j]++;
				}
				for (int p = FieldOps.roundDown(-radius-1); p < FieldOps.round(radius+1); p++)
					for (int q = FieldOps.roundDown(-radius-1); q < FieldOps.round(radius+1); q++)
					{
						if (FieldOps.withinBounds(i+p, j+q, ans) && Distance.distance(p, q) <= radius)
							n++;
					}
				ans[i][j] /= n/(Math.PI*radius*radius);
			}
		}
		return ans;
	}
	public static double[][] countImpuritiesGauss(PointImp[] imps, int nx, int ny, double radius)
	{
		double[][] ans = new double [nx][ny];
		double[] distances = new double[imps.length];
		util.scalar.functions.Functions.Gaussian2D gauss = new util.scalar.functions.Functions.Gaussian2D(0, radius);
		
		double[][] pixelWeight = new double [nx][ny];
		for (int i = 0 ; i < nx/2; i++){System.out.print(" " + i );
			for (int j = 0; j < ny/2; j++){
				for (int p = FieldOps.roundDown(i-4*radius-1); p < FieldOps.round(i + 4*radius+1); p++)
					for (int q = FieldOps.roundDown(j-4*radius-1); q < FieldOps.round(j + 4*radius+1); q++)
					{
						if (FieldOps.withinBounds(p, q, ans))
							pixelWeight[i][j] += gauss.of(Distance.distance(new int [] {p, q}, new int[] {i, j}));
					}
			}
		}
		System.out.println("\r\n" + FieldOps.max(pixelWeight));
		//only do 1/4th of the area with the stupid weights. Then symmetrize it
		for (int i = 0 ; i < nx; i++)
			for (int j = 0; j < ny; j++)
				if (i >= nx/2 && j >= ny/2)
					pixelWeight[i][j] = pixelWeight[nx-i-1][ny-j-1];
				else if (i >= nx/2)
					pixelWeight[i][j] = pixelWeight[nx-i-1][j];
				else if (j >= ny/2)
					pixelWeight[i][j] = pixelWeight[i][ny-j-1];
		System.out.println();
		for (int i = 0 ; i < nx; i++){System.out.print(" " + i );
			for (int j = 0; j < ny; j++){
				for (int k = 0; k < imps.length; k++){
					distances[k] = Distance.distance(new int [] {i, j}, imps[k].pixelPos);
					ans[i][j] += gauss.of(distances[k]);
				}
				ans[i][j] /= pixelWeight[i][j];
			}
		}
		return ans;
	}
	public static double[][][] countImpuritiesGauss(PointImp[] imps, int nx, int ny, double[] radii)
	{
		double[][][] ans = new double [radii.length][nx][ny];
		for (int i = 0; i < radii.length; i++){System.out.println("\r\n" + i + "\t");
			ans[i] = countImpuritiesGauss(imps, nx, ny, radii[i]);
		}
		return ans;
	}
	/**
	 * Here u is [nx][ny][2];
	 * @param imps
	 * @param u
	 * @return
	 */
	public static PointImp[] shiftWithUField(PointImp[] imps, double[][][] u)
	{
		BicubicSplineInterpolatingFunction interp0 = null;
		BicubicSplineInterpolator erp0 = new BicubicSplineInterpolator();
		interp0 = erp0.interpolate(ArrayOps.generateArrayInclBoth(0, u.length, u.length), ArrayOps.generateArrayInclBoth(0, u.length, u.length), FieldOps.getIndex(u, 0));
		BicubicSplineInterpolatingFunction interp1 = null;
		BicubicSplineInterpolator erp1 = new BicubicSplineInterpolator();
		interp1 = erp1.interpolate(ArrayOps.generateArrayInclBoth(0, u.length, u.length), ArrayOps.generateArrayInclBoth(0, u.length, u.length), FieldOps.getIndex(u, 1));
		
		ArrayList<PointImp> impList = new ArrayList<PointImp>();
		
		for (int i = 0; i < imps.length; i++)
		{
			try
			{
				impList.add(new PointImp(new double [] {imps[i].pixelPos[0] - interp0.value(imps[i].pixelPos[0], imps[i].pixelPos[1]) + 0.5, imps[i].pixelPos[1] - interp1.value(imps[i].pixelPos[0], imps[i].pixelPos[1])+0.5}) );
			}
			catch(org.apache.commons.math3.exception.OutOfRangeException e)
			{
				System.out.println("Impurity " + i + " was out of bounds.");
				e.printStackTrace();
			}
		}
		
		
		return getFromList(impList);
	}

	public static PointImp[] getFromList(ArrayList<PointImp> imps)
	{
		PointImp[] ans = new PointImp[imps.size()];
		for (int i = 0; i < imps.size(); i++)
			ans[i] = imps.get(i);
		return ans;
	}
	public static double[][] getPixelWeight(int nx, int ny, double radius)
	{
		util.scalar.functions.Functions.Gaussian gauss = new util.scalar.functions.Functions.Gaussian(0, radius);
		double[][] pixelWeight = new double [nx][ny];
		for (int i = 0 ; i < nx/2; i++){System.out.print(" " + i );
			for (int j = 0; j < ny/2; j++){
				for (int p = FieldOps.roundDown(i-4*radius-1); p < FieldOps.round(i + 4*radius+1); p++)
					for (int q = FieldOps.roundDown(j-4*radius-1); q < FieldOps.round(j + 4*radius+1); q++)
					{
						if (FieldOps.withinBounds(p, q, pixelWeight))
							pixelWeight[i][j] += gauss.of(Distance.distance(new int [] {p, q}, new int[] {i, j}));
					}
			}
		}
		
		//only do 1/4th of the area with the stupid weights. Then symmetrize it
		for (int i = 0 ; i < nx; i++)
			for (int j = 0; j < ny; j++)
				if (i >= nx/2 && j >= ny/2)
					pixelWeight[i][j] = pixelWeight[nx-i-1][ny-j-1];
				else if (i >= nx/2)
					pixelWeight[i][j] = pixelWeight[nx-i-1][j];
				else if (j >= ny/2)
					pixelWeight[i][j] = pixelWeight[i][ny-j-1];
		return pixelWeight;
	}
	
	/**
	 * This one returns a data field consisting of a normalized gaussian of radius radius centered at the defect site. The gaussian is also edge-normalized. 
	 * @param imps
	 * @param nx
	 * @param ny
	 * @param radius
	 * @return
	 */
	public static double[][] makeImpurityGaussMaskNormal(PointImp[] imps, int nx, int ny, double radius, double[][] pixelWeight, int k)
	{
		double[][] ans = new double [nx][ny];
		double distance;
		util.scalar.functions.Functions.Gaussian gauss = new util.scalar.functions.Functions.Gaussian(0, radius);
		
		System.out.println();
		for (int i = 0 ; i < nx; i++){System.out.print(" " + i );
			for (int j = 0; j < ny; j++){
					distance = Distance.distance(new int [] {i, j}, imps[k].pixelPos);
					ans[i][j] += gauss.of(distance);
					ans[i][j] /= pixelWeight[i][j];
			}
		}
		double sum = ArrayOps.sum(ans);
		ArrayOps.multiply(ans, 1/sum);
		
		return ans;
	}
	
	public static double[][] getLatticePositionsModuloN(PointImp[] imps, int N, double offsetX, double offsetY, AtomicCoordinatesSet latt)
	{
		double[][] ans = new double [imps.length][2];
		for (int i = 0; i < imps.length; i++)
		{
			ans[i] = latt.getAtomicCoords(imps[i].pixelPos);
			ans[i][0] %= N;
			ans[i][1] %= N;
			//These are to make sure that impurities with negative coordinates wind up within (0, N).
			ans[i][0] += N+offsetX;
			ans[i][1] += N+offsetY;
			ans[i][0] %= N;
			ans[i][1] %= N;
		}
		return ans;
	}
	public static boolean[] splitImpsWRTRoot2(PointImp[] imps, AtomicCoordinatesSet latt, boolean onAtom)
	{
		double[][] latticePos = new double [imps.length][2];
		for (int i = 0; i < imps.length; i++)
			latticePos[i] = latt.getAtomicCoords(imps[i].pixelPos);

		double x, y;
		int xm, ym;
		boolean[] subA = new boolean[imps.length];
		for (int i = 0; i < imps.length; i++)
		{
			x = latticePos[i][0];
			y = latticePos[i][1];
			//If the impurity is expected to the be the middle of the cell whose corner is (xm, ym), we should
			//round its position DOWN to find the integer atom associated with it.
			//But, if it is AT the integer position (xm, ym) we should completely round so that everything in
			//a square centered at (xm, ym) is associated with (xm, ym)
			if (!onAtom){
				xm = FieldOps.roundDown(x);
				ym = FieldOps.roundDown(y);
			}
			else
			{
				xm = FieldOps.round(x);
				ym = FieldOps.round(y);
				
			}
			subA[i] = (xm+ym)%2 == 0;
		}
		return subA;
	}
	
	/**
	 * This subtracts the offest in lattice units from each imp and then rounds to the nearest integer. We then see whether the sum is even or odd.
	 * 
	 * @param imps
	 * @param latt
	 * @param offsetX
	 * @param offsetY
	 * @return
	 */
	public static boolean[] splitImpsWRTRoot2(PointImp[] imps, AtomicCoordinatesSet latt, double offsetX, double offsetY)
	{
		double[][] latticePos = new double [imps.length][2];
		for (int i = 0; i < imps.length; i++)
			latticePos[i] = latt.getAtomicCoords(imps[i].pixelPos);

		double x, y;
		int xm, ym;
		boolean[] subA = new boolean[imps.length];
		for (int i = 0; i < imps.length; i++)
		{
			x = latticePos[i][0];
			y = latticePos[i][1];
			//If the impurity is expected to the be the middle of the cell whose corner is (xm, ym), we should
			//round its position DOWN to find the integer atom associated with it.
			//But, if it is AT the integer position (xm, ym) we should completely round so that everything in
			//a square centered at (xm, ym) is associated with (xm, ym)
			xm = FieldOps.round(x - offsetX);
			ym = FieldOps.round(y - offsetY);
			subA[i] = (xm+ym)%2 == 0;
		}
		return subA;
	}
	
	public static boolean[] splitImpsWRTLayerValue(PointImp[] imps, AtomicCoordinatesSet latt, Layer t, double threshold)
	{
		boolean[] ans = new boolean [imps.length];
		for (int i = 0; i < imps.length; i++)
			ans[i] = t.evaluateAt(imps[i].pixelPos) > threshold;
		return ans;
	}
	
	public static void shiftList(PointImp[] imps, double[] r)
	{
		for (int i = 0; i < imps.length; i++){
			imps[i].pixelPos[0] += r[0];
			imps[i].pixelPos[1] += r[1];
		}
	}
	public static PointImp[] removeOutOfBounds(PointImp[] imps, int nx, int ny)
	{
		ArrayList<PointImp> guys = new ArrayList<PointImp>();
		for (int i = 0; i < imps.length; i++){
			if (imps[i].pixelPos[0] > 0 && imps[i].pixelPos[0] < nx && imps[i].pixelPos[1] > 0 && imps[i].pixelPos[1] < ny)
				guys.add(imps[i]);
		}
		return getFromList(guys);
	}
	public static PointImp[][] splitImps(PointImp[] imps, boolean[] bools)
	{
		int na = 0, nb = 0;
		for (int i = 0; i < imps.length; i++)
			if (bools[i]) na++;
			else nb++;
		
		PointImp[] a = new PointImp[na];
		PointImp[] b = new PointImp[nb];
		
		int ia = 0, ib = 0;
		for (int i = 0; i < imps.length; i++)
			if (bools[i]) a[ia++] = imps[i];
			else b[ib++] = imps[i];
		
		return new PointImp[][] {a, b};
	}
	
	
	/**
	 * returns a list containing all the members of a, plus all the members of b which are not within epsilon of a member of a.
	 * @param a
	 * @param b
	 * @param epsilon
	 * @return
	 */
	public static PointImp[] combineTwoListsAndRemoveDuplicates(PointImp[] a, PointImp[] b, double epsilon)
	{
		ArrayList<PointImp> imps = new ArrayList<PointImp>();
		for (int i = 0; i < a.length; i++)
			imps.add(a[i]);
		
		boolean duplicate = false;
		for (int i = 0; i < b.length; i++){
			duplicate = false;
			for (int j = 0; j < a.length; j++){
				if (Distance.distance(b[i].pixelPos, a[j].pixelPos) < epsilon){
					duplicate = true;
					break;
				}
			}
			if (!duplicate)
				imps.add(b[i]);
		}
		
		return getFromList(imps);
	}
	
	public static PointSpectra getSpectraAt(PointImp[] imps, Topomap t, boolean smooth, double smoothingLength)
	{
		double[] x = new double [imps.length];
		double[] y = new double [imps.length];
		double[][] data = new double [imps.length][t.nlayers];
		
		for (int i = 0; i < imps.length; i++)
		{
			x[i] = imps[i].pixelPos[0];
			y[i] = imps[i].pixelPos[1];
		}
		
		if (!smooth)
		{
			for (int i = 0; i < imps.length; i++)
				for (int j = 0; j < t.nlayers; j++)
					data[i][j] = t.getValueAt(imps[i].pixelPos, j);
		}
		else
		{	
			double[][] temp;
			for (int j = 0; j < t.nlayers; j++){
				temp = FieldOps.gaussSmooth(t.data[j], smoothingLength);
				for (int i = 0; i < imps.length; i++)
					data[i][j] = FieldOps.getValueAt(temp, imps[i].pixelPos[0], imps[i].pixelPos[1], t.mean[j]);
			}
		}	
		return new PointSpectra(data, t.v, x, y);
	}
	public static PointSpectra getLocalDeviationSpectra(PointImp[] imps, Topomap t, double smoothingLengthSmall, double smoothingLengthLarge)
	{
		double[] x = new double [imps.length];
		double[] y = new double [imps.length];
		double[][] data = new double [imps.length][t.nlayers];
		
		for (int i = 0; i < imps.length; i++)
		{
			x[i] = imps[i].pixelPos[0];
			y[i] = imps[i].pixelPos[1];
		}
		double[][] tempS, tempL;
		for (int j = 0; j < t.nlayers; j++){
			tempS = FieldOps.gaussSmooth(t.data[j], smoothingLengthSmall);
			tempL = FieldOps.gaussSmooth(t.data[j], smoothingLengthLarge);
			for (int i = 0; i < imps.length; i++)
				data[i][j] = FieldOps.getValueAt(tempS, imps[i].pixelPos[0], imps[i].pixelPos[1], t.mean[j]) - FieldOps.getValueAt(tempL, imps[i].pixelPos[0], imps[i].pixelPos[1], t.mean[j]);
		}
		return new PointSpectra(data, t.v, x, y);
	}
	public static PointSpectra getLocalDeviationSpectraNorm(PointImp[] imps, Topomap t, double smoothingLengthSmall, double smoothingLengthLarge)
	{
		double[] x = new double [imps.length];
		double[] y = new double [imps.length];
		double[][] data = new double [imps.length][t.nlayers];
		
		for (int i = 0; i < imps.length; i++)
		{
			x[i] = imps[i].pixelPos[0];
			y[i] = imps[i].pixelPos[1];
		}
		double[][] tempS, tempL;
		double sigma;
		for (int j = 0; j < t.nlayers; j++){
			tempS = FieldOps.gaussSmooth(t.data[j], smoothingLengthSmall);
			tempL = FieldOps.gaussSmooth(t.data[j], smoothingLengthLarge);
			sigma = FieldOps.sigma(t.data[j]);
			for (int i = 0; i < imps.length; i++)
				data[i][j] = (FieldOps.getValueAt(tempS, imps[i].pixelPos[0], imps[i].pixelPos[1], t.mean[j]) - FieldOps.getValueAt(tempL, imps[i].pixelPos[0], imps[i].pixelPos[1], t.mean[j]))/sigma;
		}
		return new PointSpectra(data, t.v, x, y);
	}
	public static PointSpectra getLocalPercentiles(PointImp[] imps, Topomap t, double smoothingLengthSmall)
	{
		double[] x = new double [imps.length];
		double[] y = new double [imps.length];
		double[][] data = new double [imps.length][t.nlayers];
		
		for (int i = 0; i < imps.length; i++)
		{
			x[i] = imps[i].pixelPos[0];
			y[i] = imps[i].pixelPos[1];
		}
		double[] sort;
		double[][] tempS;
		for (int j = 0; j < t.nlayers; j++){
			tempS = FieldOps.gaussSmooth(t.data[j], smoothingLengthSmall);
			sort = FieldOps.getSortedDump(t.data[j]);
			System.out.print(j + " ");
			for (int i = 0; i < imps.length; i++)
				data[i][j] = ArrayOps.indexOf(sort, FieldOps.getValueAt(tempS, imps[i].pixelPos[0], imps[i].pixelPos[1], t.mean[j]), true)/((double)sort.length);
		}
		return new PointSpectra(data, t.v, x, y);
	}
		
	
	public static PointImp[][] writeImagesOfImps_SplitRt2(Layer t, PointImp[] imps, boolean[] subLatt, double length, int npix, String dir)
	{
		java.io.File[] subdirs = new java.io.File [2];
		subdirs[0] = new java.io.File(dir + "A\\");
		subdirs[1] = new java.io.File(dir + "B\\");
		for (int i = 0; i < 2; i++)
			if (!subdirs[i].exists()) subdirs[i].mkdirs();
		
		double pixLength = length*t.nx/t.xLength; //We use the x-direction in calculating the pixelLength.
		double[] x, y;
		BufferedImage im;
//		boolean onTrue = JOptionPane.showConfirmDialog(null, "Do we expect the impurity to be ON the lattice site (click YES) or shifted from it by half the unit cell (click NO)?") == JOptionPane.YES_OPTION; 
//		boolean[] subLatt = splitImpsWRTRoot2(imps, latt, onTrue);
		BicubicSplineInterpolator erp = new BicubicSplineInterpolator();
		BicubicSplineInterpolatingFunction interp = erp.interpolate(ArrayOps.generateArrayNotInclUpper(0, t.nx, t.nx), ArrayOps.generateArrayNotInclUpper(0, t.ny, t.ny), t.data);
		for (int i = 0 ; i < imps.length; i++)
		{
			x = ArrayOps.generateArrayInclBoth(imps[i].pixelPos[0]-pixLength/2, imps[i].pixelPos[0]+pixLength/2, npix);
			y = ArrayOps.generateArrayInclBoth(imps[i].pixelPos[1]-pixLength/2, imps[i].pixelPos[1]+pixLength/2, npix);
			im = ImageEditing.getBufferedImageBicubic(t.data, x, y, null, 17, interp);
			SRAW.writeImage(subdirs[(subLatt[i] ? 0 : 1)].toString() + "\\" + MovieMaker.fromInt(i), im);
		}
		
		return splitImps(imps, subLatt);
	}
	public static PointImp[] resizeList(PointImp[] imps, int nx, int nxNew) {
		PointImp[] ans = new PointImp[imps.length];
		double fraction = 1/(((double)nx)/((double)nxNew));
		for (int i = 0; i < imps.length; i++)
		{
			ans[i] = new PointImp(new double[] {imps[i].pixelPos[0]*fraction, imps[i].pixelPos[1]*fraction});
		}
		return ans;
	}
	
}