package util;

import image.ImageEditing;

import java.awt.image.BufferedImage;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;

import javax.swing.JFileChooser;

import drawing.LayerViewer;
import main.SRAW;
import schrodinger.MovieMaker;
import util.color.ColorScale;
import util.color.ColorScales;
import util.fileops.ColumnIO;
import util.fileops.FileOps;
import util.fileops.Layer;
import util.fileops.Topomap;
import util.geom.AtomicCoordinatesSet;
import util.geom.MapRectCoordSystem;
import util.regression.ACM_CustomFunctions;
import util.regression.ACM_NonLinearFitter;

/**
 * This is intended to contain utilities for looking at individual unit cells
 * @author madhavanlab2011
 *
 */
public class UnitCellUtil {

	public static void main(String[] args)
	{
		Topomap.setStdDir();
		JFileChooser fc = new JFileChooser(Topomap.stddir);
		
		Layer t = Layer.openFree(fc);
		AtomicCoordinatesSet latt = new AtomicCoordinatesSet(FileOps.openText(fc));
		
		int scalingFactor = 0; //every 2 scalingFactor multiplies the size of the viewed cell by 2.
		
		for(int i=0;i<scalingFactor;i++){
			latt = latt.getRt2Lattice();
		}
				
		Layer avg = getAverageUnitCell(t.data, latt, 0.25, 0.25, 30);
		LayerViewer.show(avg, 512, true);
		
//		writeLatticeSitePictures(latt, t, 64, true);
		
//		writeLatticeSitePicturesFit(latt, t, 64, true, "CosSqDimple");
//		rectifyOrigin_LowMemory(latt, t, 32, "CosSqDimple", 8, fc.getSelectedFile().getParent() + "\\");
		//		writeLatticeSiteLiteralPictures(latt, t, true, 4);
	}
	
	public static double[] getLatticeSiteMeans(AtomicCoordinatesSet latt, Layer t, int celldetail)
	{
		double[] xa = ArrayOps.generateArrayInclBoth(-0.5, 0.5, celldetail);
		double[] ya = xa.clone();
		
		UnitCell2D[] cells = getUnitCells(xa, ya, t, latt);
		double[] mean = new double[cells.length];
		for (int i = 0; i < mean.length; i++)
		{
			mean[i] = cells[i].getMean();
		}
		return mean;
	}
	
	/**
	 * Returns the list of gaussian parameters for each atom.
	 * @param latt
	 * @param t
	 * @param celldetail
	 * @param biCubic 
	 * @return the return array is [nparam][natoms] so that returnarray[i] is the list for a given parameter.
	 */
	public static double[][] fitAtomsToPeak(AtomicCoordinatesSet latt, Layer t, int celldetail, String fname, int margin)
	{
		double[] xa = ArrayOps.generateArrayInclBoth(-0.5, 0.5, celldetail);
		double[] ya = xa.clone();

		UnitCell2D[] cells = getUnitCells(xa, ya, t, latt, margin);
		
//		UnitCell2D.writeBIN(cells, "C:\\data\\analysis\\PbSnSe\\nice topos\\x03\\atom fitting\\" + "UC bicubic.bin");
		
		int nparam = ACM_CustomFunctions.getNParameters(fname);
		double [][] parameters = new double [nparam][cells.length];
		double[] temp;
		for (int i = 0; i < cells.length; i++)
		{
			System.out.print(i + " out of " + cells.length + "   \t");
			temp = ACM_NonLinearFitter.fitToFunction(cells[i].data, fname);
			System.out.print(temp[0] + "\t" + temp[1] + "\t" + temp[2] + "\t" + temp[3] + "\t" + temp[4] + "\t");
			for (int j = 0; j < temp.length; j++)
				parameters[j][i] = temp[j];
		}
		return parameters;
	}
	
	/**
	 * This fits all the atoms to peaked functions. Then, it finds out where the center is, on average. Then, it figures out the overall
	 * translation which, when applied to the lattice, will place the average center at 0, 0 in the unit cell.
	 * @param latt
	 * @param t
	 * @param celldetail
	 * @param fname
	 * @param margin
	 * @return
	 */
	public static void rectifyOrigin(AtomicCoordinatesSet latt, Layer t, int celldetail, String fname, int margin, String lattDir)
	{
		double[] xa = ArrayOps.generateArrayInclBoth(-0.5, 0.5, celldetail);
		double[] ya = xa.clone();

		UnitCell2D[] cells = getUnitCells(xa, ya, t, latt, margin);
		
//		UnitCell2D.writeBIN(cells, "C:\\data\\analysis\\PbSnSe\\nice topos\\x03\\atom fitting\\" + "UC bicubic.bin");
		
		int nparam = ACM_CustomFunctions.getNParameters(fname);
		double [][] parameters = new double [nparam][cells.length];
		double[] temp;
		for (int i = 0; i < cells.length; i++)
		{
			System.out.print(i + " out of " + cells.length + "   \t");
			temp = ACM_NonLinearFitter.fitToFunction(cells[i].data, fname);
			System.out.print(temp[0] + "\t" + temp[1] + "\t" + temp[2] + "\t" + temp[3] + "\t" + temp[4] + "\t");
			for (int j = 0; j < temp.length; j++)
				parameters[j][i] = temp[j];
		}
		
		double[][] purgedParams = LayerUtil.purgeUnwantedAtoms(LayerUtil.chooseWantedAtoms(parameters, celldetail), parameters);
		
		double meanx = ArrayOps.mean(purgedParams[0]);
		double meany = ArrayOps.mean(purgedParams[1]);
		
		double[] latticeUnitsPt = cells[0].getMetricCoords(meanx, meany);
		double[] pixelTrans = latt.getPixelCoords(latticeUnitsPt);
		pixelTrans[0] -= latt.getOrigin()[0];
		pixelTrans[1] -= latt.getOrigin()[1];
		
		latt.moveOrigin(pixelTrans[0], pixelTrans[1]);
		FileOps.writeString(latt.toString());
		ColumnIO.writeString(latt.toString(), lattDir + "1x1 translated.txt");
	}
	
	/**
	 * This one does not create the full array of unit cells, but rather, creates and fits them one at a time.
	 * @param latt
	 * @param t
	 * @param celldetail
	 * @param fname
	 * @param margin
	 * @param lattDir
	 */
	public static void rectifyOrigin_LowMemory(AtomicCoordinatesSet latt, Layer t, int celldetail, String fname, int margin, String lattDir)
	{
		double[] xa = ArrayOps.generateArrayInclBoth(-0.5, 0.5, celldetail);
		double[] ya = xa.clone();

//		UnitCell2D[] cells = getUnitCells(xa, ya, t, latt, margin);
		ArrayList<int[]> sites = getLatticeSites(latt, t.nx, margin);
		UnitCell2D cell = null;
//		UnitCell2D.writeBIN(cells, "C:\\data\\analysis\\PbSnSe\\nice topos\\x03\\atom fitting\\" + "UC bicubic.bin");
		
		int nparam = ACM_CustomFunctions.getNParameters(fname);
		double [][] parameters = new double [nparam][sites.size()];
		double[] temp;
		for (int i = 0; i < sites.size(); i++)
		{
			cell = new UnitCell2D(xa, ya, sites.get(i)[0], sites.get(i)[1], t, latt);
			
			System.out.print(i + " out of " + sites.size() + "   \t");
			temp = ACM_NonLinearFitter.fitToFunction(cell.data, fname);
			System.out.print(temp[0] + "\t" + temp[1] + "\t" + temp[2] + "\t" + temp[3] + "\t" + temp[4] + "\t");
			for (int j = 0; j < temp.length; j++)
				parameters[j][i] = temp[j];
		}
		
		double[][] purgedParams = LayerUtil.purgeUnwantedAtoms(LayerUtil.chooseWantedAtoms(parameters, celldetail), parameters);
		
		double meanx = ArrayOps.mean(purgedParams[0]);
		double meany = ArrayOps.mean(purgedParams[1]);
		
		double[] latticeUnitsPt = cell.getMetricCoords(meanx, meany);
		double[] pixelTrans = latt.getPixelCoords(latticeUnitsPt);
		pixelTrans[0] -= latt.getOrigin()[0];
		pixelTrans[1] -= latt.getOrigin()[1];
		
		latt.moveOrigin(pixelTrans[0], pixelTrans[1]);
//		FileOps.writeString(latt.toString());
		ColumnIO.writeString(latt.toString(), lattDir + "1x1 translated.txt");
	}
	
	public static double[][] fitAtomsToPeakWithPicture(AtomicCoordinatesSet latt, Layer t, String fname, String picDir, boolean absoluteScale, UnitCell2D[] cells)
	{
//		UnitCell2D.writeBIN(cells, "C:\\data\\analysis\\PbSnSe\\nice topos\\x03\\atom fitting\\" + "UC bicubic.bin");
		
		ColorScale c = null;
		if (absoluteScale) c = ColorScales.getNew(t.data, 1);
		
		if (!new File(picDir).exists()) new File(picDir).mkdir();
		
		int blowUpFactor = 2;
		
		int nparam = ACM_CustomFunctions.getNParameters(fname);
		double [][] parameters = new double [nparam][cells.length];
		double[][] fitData;
		double[] temp;
		for (int i = 0; i < cells.length; i++)
		{
			System.out.print(i + " out of " + cells.length + "   \t");
			temp = ACM_NonLinearFitter.fitToFunction(cells[i].data, fname);
			System.out.print(temp[0] + "\t" + temp[1] + "\t" + temp[2] + "\t" + temp[3] + "\t" + temp[4] + "\t");
			for (int j = 0; j < temp.length; j++)
				parameters[j][i] = temp[j];

			if (!absoluteScale) c = ColorScales.getNew(cells[i].data, 1);
			SRAW.writeImage(picDir + "z" + MovieMaker.fromInt(i) + "cell", ImageEditing.getEnlarged(ImageEditing.getBufferedImage(cells[i].data, c), blowUpFactor));
			fitData = ACM_NonLinearFitter.getExpectedData(cells[i].nx, cells[i].ny, temp, fname);
			SRAW.writeImage(picDir + "z" + MovieMaker.fromInt(i) + "fit", ImageEditing.getEnlarged(ImageEditing.getBufferedImage(fitData, c), blowUpFactor));
		}
		return parameters;
	}
	public static double[] getLatticeSiteMeans(AtomicCoordinatesSet latt, Layer t, int celldetail, double middleFrac)
	{
		double[] xa = ArrayOps.generateArrayInclBoth(-middleFrac/2, middleFrac/2, celldetail);
		double[] ya = xa.clone();
		
		UnitCell2D[] cells = getUnitCells(xa, ya, t, latt);
		double[] mean = new double[cells.length];
		for (int i = 0; i < mean.length; i++)
		{
			mean[i] = cells[i].getMean();
		}
		return mean;
	}
	public static double[] getLatticeSiteMaxima(AtomicCoordinatesSet latt, Layer t, double fraction)
	{
		double[] xa = ArrayOps.generateArrayInclBoth(-fraction/2, fraction/2, 2);
		double[] ya = xa.clone();
		
		UnitCell2D[] cells = getUnitCells(xa, ya, t, latt);
		double[] max = new double[cells.length];
		for (int i = 0; i < max.length; i++)
		{
			max[i] = ArrayOps.max(cells[i].literalData);
		}
		return max;
	}
	public static UnitCell2D[] getUnitCells(AtomicCoordinatesSet latt, Layer t, double fraction, int npts)
	{
		double[] xa = ArrayOps.generateArrayInclBoth(-fraction/2, fraction/2, npts);
		double[] ya = xa.clone();
		
		UnitCell2D[] cells = getUnitCells(xa, ya, t, latt);
		return cells;
		
	}
	public static double[] getLatticeSiteMeansLiteral(AtomicCoordinatesSet latt, Layer t, double middleFrac)
	{
		double[] xa = ArrayOps.generateArrayInclBoth(-middleFrac/2, middleFrac/2, 4);
		double[] ya = xa.clone();
		
		UnitCell2D[] cells = getUnitCells(xa, ya, t, latt);
		double[] mean = new double[cells.length];
		for (int i = 0; i < mean.length; i++)
		{
			mean[i] = cells[i].getMeanLiteral(t, latt);
		}
		return mean;
	}
	public static void writeLatticeSitePictures(AtomicCoordinatesSet latt, Layer t, int celldetail, boolean useOneScale)
	{
		double[] xa = ArrayOps.generateArrayInclBoth(-0.5, 0.5, celldetail);
		double[] ya = xa.clone();
		
		UnitCell2D[] cells = getUnitCells(xa, ya, t, latt);
		File base = FileOps.selectSave(null);
		ColorScale c = ColorScales.getNew(t.data, 1);
		for (int i = 0; i < cells.length; i++)
		{
			if (useOneScale)
				SRAW.writeImage(base.toString() + MovieMaker.fromInt(i), cells[i].getImage(c));
			else
				SRAW.writeImage(base.toString() + MovieMaker.fromInt(i), cells[i].getImage(ColorScales.getNew(cells[i].data, 1)));
		}
		
	}
	public static void writeLatticeSitePicturesFit(AtomicCoordinatesSet latt, Layer t, int celldetail, boolean useOneScale, String fname)
	{
		double[] xa = ArrayOps.generateArrayInclBoth(-0.5, 0.5, celldetail);
		double[] ya = xa.clone();
		
		UnitCell2D[] cells = getUnitCells(xa, ya, t, latt);
		File base = FileOps.selectSave(null);
		ColorScale c = ColorScales.getNew(t.data, 1);
		double[] temp;
		double[][] fit;
		for (int i = 0; i < cells.length; i++)
		{
//			if (i != 0) System.exit(0);
			temp = ACM_NonLinearFitter.fitToFunction(cells[i].data, fname);
			fit = ACM_NonLinearFitter.getExpectedData(celldetail, celldetail, temp, fname);
			
			if (useOneScale){
				SRAW.writeImage(base.toString() + MovieMaker.fromInt(i), cells[i].getImage(c));
//				SRAW.writeImage(base.toString() + "fit" + MovieMaker.fromInt(i), ImageEditing.getBufferedImage(fit, c));
			}
			else{
				SRAW.writeImage(base.toString() + MovieMaker.fromInt(i), cells[i].getImage(ColorScales.getNew(cells[i].data, 1)));
			}
		}
		
	}
	public static void writeLatticeSiteLiteralPictures(AtomicCoordinatesSet latt, Layer t, boolean useOneScale, int blowUpFactor)
	{
		double[] xa = ArrayOps.generateArrayInclBoth(-0.5, 0.5, 2);
		double[] ya = xa.clone();
		
		UnitCell2D[] cells = getUnitCells(xa, ya, t, latt);
		File base = FileOps.selectSave(null);
		ColorScale c = ColorScales.getNew(t.data, 1);
		for (int i = 0; i < cells.length; i++)
		{
			if (useOneScale)
				SRAW.writeImage(base.toString() + MovieMaker.fromInt(i), cells[i].getLiteralImage(c, t, latt, blowUpFactor));
			else
				SRAW.writeImage(base.toString() + MovieMaker.fromInt(i), cells[i].getLiteralImage(ColorScales.getNew(cells[i].data, 1), t, latt, blowUpFactor));
		}
		
	}
	
	/**
	 * Returns the integer lattice coordinates of all unit cells located within the range 0, N in x and y.
	 * @param a
	 * @param N
	 * @return
	 */
	public static ArrayList<int[]> getLatticeSites(AtomicCoordinatesSet a, int N)
	{
		 double[] atomc1 = a.getAtomicCoords(0, 0);
		 double[] atomc2 = a.getAtomicCoords(0, N);
		 double[] atomc3 = a.getAtomicCoords(N, 0);
		 double[] atomc4 = a.getAtomicCoords(N, N);
		 double[] min = new double [2], max = new double [2];
		 min[0] = Math.min(Math.min(atomc1[0], atomc2[0]),Math.min(atomc3[0], atomc4[0]));
		 min[1] = Math.min(Math.min(atomc1[1], atomc2[1]),Math.min(atomc3[1], atomc4[1]));
		 max[0] = Math.max(Math.max(atomc1[0], atomc2[0]),Math.max(atomc3[0], atomc4[0]));
		 max[1] = Math.max(Math.max(atomc1[1], atomc2[1]),Math.max(atomc3[1], atomc4[1]));
		 
		 
		 int xmn = FieldOps.round(min[0]);
		 int xmx = FieldOps.round(max[0])+1;
		 int ymn = FieldOps.round(min[1]);
		 int ymx = FieldOps.round(max[1])+1;
		 double[] pix = new double [2];
		 ArrayList<int[]> sites = new ArrayList<int[]>();
		 
		 for (int i = xmn; i < xmx; i++)
			 for (int j = ymn; j < ymx; j++)
			 {
				 pix = a.getPixelCoords(i, j);
				 if (pix[0] > 0 && pix[1] > 0 && pix[0] < N-1 && pix[1] < N-1)
					 sites.add(new int[] {i, j});
			 }
		 return sites;
	}
	public static void purgeLatticeSites(Layer t, AtomicCoordinatesSet a, ArrayList<int[]> sites)
	{
		 
		 double[] pix = new double [2];
		 double mean = //FieldOps.mean(t.data);
		t.data[0][t.nx-1]; 
		for (int i = 0; i < sites.size(); i++){
				 pix = a.getPixelCoords(sites.get(i)[0], sites.get(i)[1]);
				 double d = t.data[FieldOps.round(pix[0])][FieldOps.round(pix[1])];
//				 System.out.println(d);
				 if (d == mean){
					 sites.remove(i);
					 i--;
				 }
		 }
	}
	
	/**
	 * Returns all lattice sites a distance greater than margin from the edge.
	 * @param a
	 * @param N
	 * @param margin
	 * @return
	 */
	public static ArrayList<int[]> getLatticeSites(AtomicCoordinatesSet a, int N, int margin)
	{
		 double[] atomc1 = a.getAtomicCoords(0, 0);
		 double[] atomc2 = a.getAtomicCoords(0, N);
		 double[] atomc3 = a.getAtomicCoords(N, 0);
		 double[] atomc4 = a.getAtomicCoords(N, N);
		 double[] min = new double [2], max = new double [2];
		 min[0] = Math.min(Math.min(atomc1[0], atomc2[0]),Math.min(atomc3[0], atomc4[0]));
		 min[1] = Math.min(Math.min(atomc1[1], atomc2[1]),Math.min(atomc3[1], atomc4[1]));
		 max[0] = Math.max(Math.max(atomc1[0], atomc2[0]),Math.max(atomc3[0], atomc4[0]));
		 max[1] = Math.max(Math.max(atomc1[1], atomc2[1]),Math.max(atomc3[1], atomc4[1]));
		 
		 
		 int xmn = FieldOps.round(min[0]);
		 int xmx = FieldOps.round(max[0])+1;
		 int ymn = FieldOps.round(min[1]);
		 int ymx = FieldOps.round(max[1])+1;
		 double[] pix = new double [2];
		 ArrayList<int[]> sites = new ArrayList<int[]>();
		 
		 for (int i = xmn; i < xmx; i++)
			 for (int j = ymn; j < ymx; j++)
			 {
				 pix = a.getPixelCoords(i, j);
				 if (pix[0] > margin && pix[1] > margin && pix[0] < N-1-margin && pix[1] < N-1-margin)
					 sites.add(new int[] {i, j});
			 }
		 return sites;
	}
	/**
	 * This method updated 9/13/2014 to exclude data where the layer value is equal to the mean of the layer. This is intended to get rid of "edge space.";
	 * @param xa
	 * @param ya
	 * @param t
	 * @param latt
	 * @return
	 */
	public static UnitCell2D[] getUnitCells(double[] xa, double[] ya, Layer t, AtomicCoordinatesSet latt)
	{
		int N = t.nx;
		ArrayList<int[]> sites = getLatticeSites(latt, N);
		purgeLatticeSites(t, latt, sites);
		UnitCell2D[] cells = new UnitCell2D[sites.size()];
		for (int i = 0; i < sites.size(); i++)
		{
			cells[i] = new UnitCell2D(xa, ya, sites.get(i)[0], sites.get(i)[1], t, latt);
			if (i % 1000 == 0) System.out.println("" + i + " out of " + sites.size());
		}
		return cells;
	}
	public static UnitCell2D[] getUnitCells(double fraction, Layer t, AtomicCoordinatesSet latt, int celldetail, int margin)
	{
		double[] xa = ArrayOps.generateArrayInclBoth(-fraction/2, fraction/2, celldetail);
		double[] ya = xa.clone();
	
		UnitCell2D[] cells = getUnitCells(xa, ya, t, latt, margin);
		return cells;
	}
	public static UnitCell2D[] getUnitCells(double[] xa, double[] ya, Layer t, AtomicCoordinatesSet latt, int margin)
	{
		int N = t.nx;
		ArrayList<int[]> sites = getLatticeSites(latt, N, margin);
		UnitCell2D[] cells = new UnitCell2D[sites.size()];
		for (int i = 0; i < sites.size(); i++)
		{
			cells[i] = new UnitCell2D(xa, ya, sites.get(i)[0], sites.get(i)[1], t, latt);
			if (i % 1000 == 0) System.out.println("" + i + " out of " + sites.size());
		}
		return cells;
	}
	/**
	 * This class contains the data for a unit cell of the topography. 
	 * @author madhavanlab2011
	 *
	 */
	public static class UnitCell2D extends MapRectCoordSystem{
		

		private int a; //the coordinates of the unit cell which this is.
		private int b;
		double[][] data;
		public double min, max;
		double valueAtInt;
		double[] literalData;
		int[][] literalDataPositions;
		//The fields in MapRectCoordSystem will be in lattice units.
		public UnitCell2D(double[] x, double[] y, int a, int b, Layer t, AtomicCoordinatesSet latt) {
			super(x, y);
			this.a = a;
			this.b = b;
			
			data = new double [nx][ny];
			setDataBicubic(latt, t);
			setLiteralData1D(latt, t);
		}
		/**
		 * This one reads the data directly. The user is responsible for making sure that Layer t is the same layer as the "data" was originally taken via bicubic interpolation.
		 * (This is for loading from a file without repeating the bicubic interpolation.)
		 * @param x
		 * @param y
		 * @param a
		 * @param b
		 * @param t
		 * @param latt
		 */
		public UnitCell2D(double[] x, double[] y, int a, int b, double[][] data, Layer t, AtomicCoordinatesSet latt) {
			super(x, y);
			this.a = a;
			this.b = b;
			this.data = data;
				
			setLiteralData1D(latt, t);
		}
		
		/**
		 * This does not include the literal data which can be extracted without rebuilding the bicubic interpolator, but it does include "data".
		 * @return
		 */
		public int getByteLength()
		{
			return 8*x.length + 8*y.length + 8 + 8 + 8*data.length*data[0].length;
		}

		public double[] getLiteralData() {
			return literalData;
		}
		public double[][] getData()
		{
			return data;
		}

		public double getMean()
		{
			return FieldOps.mean(data);
		}
		
		/**
		 * This returns the mean of all ORIGINAL pixels which are inside the unit cell.
		 * @return
		 */
		public double getMeanLiteral(Layer l, AtomicCoordinatesSet latt)
		{
			int npix = 0;
			double sum = 0;
			int[] bounds = getLiteralBounds(latt, l);
			for (int i = bounds[0]; i < bounds[1]; i++)
				for (int j = bounds[2]; j < bounds[3]; j++)
					if (isCoordInside(latt.getAtomicCoords(i, j)))
					{
						sum += l.data[i][j];
						npix++;
					}
			System.out.println(npix + "\t" + bounds[0] + "\t" + bounds[1] + "\t" + bounds[2] + "\t" + bounds[3]);
			return sum/npix;
		}
		
		/**
		 * returns the bounds of the square entirely containing the unit cell. If the edge is near, the bounds are limited thereto.
		 * returns {xmin, xmax, ymin, ymax}.
		 * The upper bound is meant to be used in a for loop with < (so the array length is included)
		 * @param latt
		 * @param t
		 * @return
		 */
		public int[] getLiteralBounds(AtomicCoordinatesSet latt, Layer t)
		{
			double[] xa = new double [4];
			double[] ya = new double [4];
			double[] temp = new double [2];
			
			temp = latt.getPixelCoords(getA() + x[0], getB() + y[0]);
			xa[0] = temp[0]; ya[0] = temp[1];
			temp = latt.getPixelCoords(getA() + x[x.length-1], getB() + y[0]);
			xa[1] = temp[0]; ya[1] = temp[1];
			temp = latt.getPixelCoords(getA() + x[0], getB() + y[y.length-1]);
			xa[2] = temp[0]; ya[2] = temp[1];
			temp = latt.getPixelCoords(getA() + x[x.length-1], getB() + y[y.length-1]);
			xa[3] = temp[0]; ya[3] = temp[1];
			
			int[] ans = new int[] {
							(int)(ArrayOps.min(xa)),
							(int)(ArrayOps.max(xa)+2),
							(int)(ArrayOps.min(ya)),
							(int)(ArrayOps.max(ya)+2)
					};
			if (ans[0] < 0) ans[0] = 0;
			if (ans[2] < 0) ans[2] = 0;
			if (ans[1] > t.nx) ans[1] = t.nx;
			if (ans[3] > t.ny) ans[3] = t.ny;
			return ans;
		}
		
		public double[][] getLiteralData(AtomicCoordinatesSet latt, Layer t, boolean flatOutside)
		{
			int[] bounds = getLiteralBounds(latt, t);
			double[][] lit = new double [bounds[1]-bounds[0]][bounds[3]-bounds[2]];
			for (int i = bounds[0]; i < bounds[1]; i++)
				for (int j = bounds[2]; j < bounds[3]; j++)
					if (isCoordInside(latt.getAtomicCoords(i, j)) || !flatOutside)
						lit[i-bounds[0]][j-bounds[2]] = t.data[i][j];
					else
						lit[i-bounds[0]][j-bounds[2]] = 0;
			
			FieldOps.changeZeroToAverage(lit);
			return lit;
		}
		/**
		 * This method sets the data array with bicubic interpolation, without interpolating the entire layer. Instead, a
		 * "literal data" array is interpolated.
		 * @param latt
		 * @param t
		 */
		public void setDataBicubic(AtomicCoordinatesSet latt, Layer t) {
			// TODO Auto-generated method stub
			int[] bounds = this.getLiteralBounds(latt, t);
			double[][] lit = this.getLiteralData(latt, t, false);
			
			//copy the "axis label" arrays
			double[] xn = new double [lit.length];
			double[] yn = new double [lit[0].length];
			for (int i = 0; i < xn.length; i++)
				xn[i] = i+bounds[0];
			for (int i = 0; i < yn.length; i++)
				yn[i] = i+bounds[2];
			
			Layer temp = new Layer(lit, xn, yn, 0, 0);
			for (int i = 0; i < nx; i++)
				for (int j = 0; j < ny; j++)
				{
					data[i][j] = temp.evaluateBiCubicMetric(latt.getPixelCoords(x[i] + a, y[j] + b));
				}
			min = FieldOps.min(data);
			max = FieldOps.max(data);
			valueAtInt = temp.evaluateBiCubicMetric(latt.getPixelCoords(a, b));
		}

		public void setLiteralData1D(AtomicCoordinatesSet latt, Layer t)
		{
			int[] bounds = getLiteralBounds(latt, t);
			int n = 0;
			for (int i = bounds[0]; i < bounds[1]; i++)
				for (int j = bounds[2]; j < bounds[3]; j++)
					if (isCoordInside(latt.getAtomicCoords(i, j)))
						n++; else;
			literalData = new double [n];
			literalDataPositions = new int [n][2];
			n = 0;		
			for (int i = bounds[0]; i < bounds[1]; i++)
				for (int j = bounds[2]; j < bounds[3]; j++)
					if (isCoordInside(latt.getAtomicCoords(i, j)))
					{
							literalData[n] = t.data[i][j];
							literalDataPositions[n] = new int[] {i, j};
							n++;
					}
		}
		public boolean isCoordInside(double[] latticeCoord)
		{
			return latticeCoord[0] >= getA() + x[0] && latticeCoord[0] <= getA() + x[x.length-1] && latticeCoord[1] >= getB() + y[0] && latticeCoord[1] <= getB() + y[y.length-1];
		}
		public BufferedImage getImage(ColorScale s)
		{
			return ImageEditing.getBufferedImage(data, s);
		}
		
		public BufferedImage getLiteralImage(ColorScale s, Layer t, AtomicCoordinatesSet latt, int blowUpFactor)
		{
			if (blowUpFactor == 1)
				return ImageEditing.getBufferedImage(getLiteralData(latt, t, true), s);
			else return ImageEditing.getEnlarged(ImageEditing.getBufferedImage(getLiteralData(latt, t, true), s), blowUpFactor);
		}

		public int getA() {
			return a;
		}
		public int getB() {
			return b;
		}
		public static void writeBIN(UnitCell2D[] cells, String filepath)
		{
			int byteLength = 4;
			for (int i = 0; i < cells.length; i++)
				byteLength += cells[i].getByteLength();
			File file = new File(filepath);
			
			FileOutputStream outf = null;
			BufferedOutputStream outbuff = null;
			DataOutputStream outd = null;
			
			try {
				outf = new FileOutputStream(file);
				outbuff = new BufferedOutputStream(outf, byteLength/4);
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
			outd = new DataOutputStream(outbuff);
			try {
				outd.writeInt(cells.length);
				for (int i = 0; i < cells.length; i++)
				{
					outd.writeInt(cells[i].a);
					outd.writeInt(cells[i].b);
					outd.writeInt(cells[i].nx);
					outd.writeInt(cells[i].ny);
					for (int j = 0; j < cells[i].nx; j++)
						outd.writeDouble(cells[i].x[j]);
					for (int j = 0; j < cells[i].nx; j++)
						outd.writeDouble(cells[i].y[j]);
					for (int j = 0; j < cells[i].nx; j++)
						for (int k = 0; k < cells[i].ny; k++)
							outd.writeDouble(cells[i].data[j][k]);
				}
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			try {
//				outbuff.flush();
//				outf.flush();
				outd.close();
				outbuff.close();
				outf.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		public double getValueAtInt()
		{
			return valueAtInt;
		}
		public static UnitCell2D[] readBin(String filepath, AtomicCoordinatesSet latt, Layer t)
		{
			int nx, ny, a, b, ncells;
			double[] x = null, y = null;
			double[][] data = null;
			File file = new File(filepath);
			FileInputStream inf = null;
			BufferedInputStream inbuff = null;
			DataInputStream ind = null;
			try {
				inf = new FileInputStream(file);
				inbuff = new BufferedInputStream(inf, (int)file.length()/4);
			} catch (FileNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
				return null;
			}
			ind = new DataInputStream(inbuff);
			UnitCell2D[] cells = null;
			try {
				ncells = ind.readInt();
				cells = new UnitCell2D[ncells];
				for (int i = 0; i < ncells; i++)
				{
					a = ind.readInt();
					b = ind.readInt();
					nx = ind.readInt();
					ny = ind.readInt();
					x = new double [nx];
					y = new double [ny];
					data = new double [nx][ny];
					for (int j = 0; j < nx; j++)
						x[j] = ind.readInt();
					for (int j = 0; j < ny; j++)
						y[j] = ind.readInt();
					for (int j = 0; j < nx; j++)
						for (int k = 0; k < ny; k++)
							data[j][k] = ind.readDouble();
					cells[i] = new UnitCell2D(x, y, a, b, data, t, latt);
				}
				ind.close();
				inbuff.close();
				inf.close();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		return cells;
	}
		
		public static class UnitCellReferenceTable
		{
			public UnitCell2D[] cells;
			
			public int[][] table;
			public int amin, bmin;
			
			public UnitCellReferenceTable(UnitCell2D[] cells)
			{
				this.cells = cells;
				amin = Integer.MAX_VALUE;
				bmin = Integer.MAX_VALUE;
				int amax = Integer.MIN_VALUE, bmax = Integer.MIN_VALUE;
				
				for (int i = 0; i < cells.length; i++)
				{
					amin = Math.min(cells[i].a, amin);
					bmin = Math.min(cells[i].b, bmin);
					amax = Math.max(cells[i].a, amax);
					bmax = Math.max(cells[i].b, bmax);
				}
				
				table = new int [amax-amin + 1][bmax - bmin + 1];
				for (int i = 0; i < table.length; i++)
					for (int j = 0; j < table[0].length; j++)
						table[i][j] = -1;
				
				for (int i = 0; i < cells.length; i++)
					table[cells[i].a - amin][cells[i].b - bmin] = i;
			}
			
			public int getIndex(int a, int b)
			{
				return table[a-amin][b-bmin];
			}
		}
	}
	public static double[] getValuesAtInt(AtomicCoordinatesSet latt, Layer t)
	{
		double[] xa = ArrayOps.generateArrayInclBoth(-0.5, 0.5, 2);
		double[] ya = xa.clone();
		
		UnitCell2D[] cells = getUnitCells(xa, ya, t, latt);
		double[] val = new double[cells.length];
		for (int i = 0; i < cells.length; i++)
		{
			val[i] = cells[i].valueAtInt;
		}
		return val;
	}
	
	public static Layer getAverageUnitCell(double[][] source, AtomicCoordinatesSet latt, double aOffset, double bOffset, int detail)
	{
		double[][] nentries = new double [detail][detail];
		double[][] avg = new double [detail][detail];
		
		for(int i=0;i<detail;i++){
			for(int j=0;j<detail;j++){
				avg[i][j]=0;
				nentries[i][j]=0;
			}
		}
		
		int offsetA = (int)(aOffset*detail);
		int offsetB = (int)(bOffset*detail);
		double[] unity = ArrayOps.generateArrayNotInclLower(0, 1, detail);
		
		int nx = source.length, ny = source[0].length;
		for (int i =0; i < nx; i++)
			for (int j = 0; j < ny; j++)
			{
				double[] lp = latt.getAtomicCoords(i, j);
//				System.out.print(Printer.arrayLnHorizontal(lp));
				lp[0] %= 1;
				lp[1] %= 1;
				lp[0] += 1;
				lp[1] += 1;
				lp[0] %= 1;
				lp[1] %= 1;
//				System.out.println(Printer.arrayLnHorizontal(lp));
				
				
				int la = (int) Math.floor(lp[0]*detail);
				int lb = (int) Math.floor(lp[1]*detail);
				
				double lx = la-lp[0]*detail;
				double ly = lb-lp[1]*detail;
				
				avg[la][lb] += source[i][j]*(1-lx)*(1-ly);
				nentries[la][lb] += (1-lx)*(1-ly);
				if(la==0){
					avg[detail-1][lb] += source[i][j]*(lx)*(1-ly);
					nentries[detail-1][lb] += (lx)*(1-ly);
					
					if(lb==0){
						avg[detail-1][detail-1] += source[i][j]*(lx)*(ly);
						nentries[detail-1][detail-1] += (lx)*(ly);
						
						avg[la][detail-1] += source[i][j]*(1-lx)*(ly);
						nentries[la][detail-1] += (1-lx)*(ly);
					}else{
						avg[detail-1][lb-1] += source[i][j]*(lx)*(ly);
						nentries[detail-1][lb-1] += (lx)*(ly);
						
						avg[la][lb-1] += source[i][j]*(1-lx)*(ly);
						nentries[la][lb-1] += (1-lx)*(ly);
					}
				}else{
					avg[la-1][lb] += source[i][j]*(lx)*(1-ly);
					nentries[la-1][lb] += (lx)*(1-ly);
					
					if(lb==0){
						avg[la-1][detail-1] += source[i][j]*(lx)*(ly);
						nentries[la-1][detail-1] += (lx)*(ly);
						
						avg[la][detail-1] += source[i][j]*(1-lx)*(ly);
						nentries[la][detail-1] += (1-lx)*(ly);
					}else{
						avg[la-1][lb-1] += source[i][j]*(lx)*(ly);
						nentries[la-1][lb-1] += (lx)*(ly);
						
						avg[la][lb-1] += source[i][j]*(1-lx)*(ly);
						nentries[la][lb-1] += (1-lx)*(ly);
					}
				}

			}
		LayerViewer.show(Layer.getFreeLayer(nentries), 512, true);
		for (int i = 0; i < detail; i++)
			for (int j = 0; j < detail; j++)
			{
				if (nentries[i][j] != 0) avg[i][j] /= nentries[i][j];
			}
		
		double[][] ans = new double [detail][detail];
		FieldOps.shift(avg, ans, offsetA, offsetB);
		return new Layer(ans, unity, unity, 1, 1);
	}
}
