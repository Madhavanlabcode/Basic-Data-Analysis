package util;

import impurity.PointImp;

import java.io.File;
import java.util.ArrayList;

import javax.swing.JFileChooser;

import drawing.GraphDrawerCart;
import drawing.LineCutDrawer;
import main.SRAW;
import flanagan.analysis.Regression;
import util.fileops.ColumnIO;
import util.fileops.FileOps;
import util.fileops.Layer;
import util.fileops.PointSpectra;
import util.fileops.Topomap;
import util.fourier.FFTOps;
import util.landau.Eq48Zhang;
import util.landau.LandauLevelModel;
import util.regression.ACM_NonLinearFitter;
import util.regression.ACM_NonLinearFitter.FittingResult;
import util.scalar.functions.Functions.CubicSplineWrapper;

public class SpectraUtil {
	
	
	static JFileChooser fc;
	
	public static PointSpectra deleteFlatSpectra(PointSpectra t){
		//deletes spectra whose endpoints and center point are all the same
		boolean[] isFlat = new boolean[t.nspec];
		int nFlat = 0;
		
		for(int x=t.nspec;x>0;x--){
			if(t.data[x-1][0]==t.data[x-1][t.nlayers-1] && t.data[x-1][0]==t.data[x-1][t.nlayers/2]){
				isFlat[x-1]=true;
				nFlat++;
			}else{
				isFlat[x-1]=false;
			}
		}
		System.out.println(nFlat + " flat spectra.");
		
		double[] x = new double [t.nspec-nFlat];
		double[] y = new double [t.nspec-nFlat];
		double[][] data = new double[t.nspec-nFlat][t.nlayers];
		
		int nSkipped=0;
		for (int j = 0; j < t.nspec; j++){
			if (!isFlat[j]){
				x[j-nSkipped] = t.x[j];
				y[j-nSkipped] = t.y[j];
				for (int k = 0; k < t.nlayers; k++)
					data[j-nSkipped][k] = t.data[j][k];
			}else{
				nSkipped++;
			}
		}

		return new PointSpectra(data, t.v, x, y);
	}
	
	//cuts off the extreme voltages from a PointSpectra object.
	public static PointSpectra truncate(PointSpectra data, double min, double max, boolean specFlipped)
	{
		int imin = ArrayOps.indexOf(data.v, min, !specFlipped);
		int imax = ArrayOps.indexOf(data.v, max, !specFlipped);
		if (!specFlipped)
			return PointSpectra.truncateTo(data, imin, imax);
		else
			return PointSpectra.truncateTo(data, imax, imin);
	}
	public static PointSpectra truncate(PointSpectra data, int imin, int imax, boolean specFlipped)
	{
		if (!specFlipped)
			return PointSpectra.truncateTo(data, imin, imax);
		else
			return PointSpectra.truncateTo(data, imax, imin);
	}
	public static PointSpectra deleteOneSpectrum(PointSpectra t, int i)
	{
		double[] x = new double [t.nspec-1];
		double[] y = new double [t.nspec-1];
		double[][] data = new double[t.nspec-1][t.nlayers];
		for (int j = 0; j < t.nspec; j++)
			if (j < i){
				x[j] = t.x[j];
				y[j] = t.y[j];
				for (int k = 0; k < t.nlayers; k++)
					data[j][k] = t.data[j][k];
			}
			else if (j > i)
			{
				x[j-1] = t.x[j];
				y[j-1] = t.y[j];
				for (int k = 0; k < t.nlayers; k++)
					data[j-1][k] = t.data[j][k];
			}
		return new PointSpectra(data, t.v, x, y);
	}
	public static PointSpectra getSubset(PointSpectra t, int imin, int imax)
	{
		int di = imax - imin;
		double[] x = new double [di];
		double[] y = new double [di];
		double[][] data = new double[di][t.nlayers];
		for (int j = imin; j < imax; j++)
		{
			x[j-imin] = t.x[j];
			y[j-imin] = t.y[j];
			for (int k = 0; k < t.nlayers; k++)
				data[j-imin][k] = t.data[j][k];
		}
		return new PointSpectra(data, t.v, x, y);
	}
	public static PointSpectra collapseTimeWise(PointSpectra t, int factor)
	{
		int di = t.nspec/factor;
		double[] x = new double [di];
		double[] y = new double [di];
		double[][] data = new double[di][t.nlayers];
		for (int j = 0; j < di; j++)
		{
			x[j] = t.x[j];
			y[j] = t.y[j];
			for (int k = 0; k < t.nlayers; k++)
				for (int i = 0; i < factor; i++)
					data[j][k] += t.data[j*factor + i][k]/factor;
		}
		return new PointSpectra(data, t.v, x, y);
	}
	
	//The first is the FFT of the individual spectra. The second is the FFT of the average of them all.
	public static PointSpectra[] getFFT(PointSpectra t)
	{
		double[][] data = new double [t.nspec][t.nlayers];

		double[] fftAvg;
		
		double[] newV = ArrayOps.generateArrayInclBoth(0, t.nlayers-1, t.nlayers);
		
		for (int i = 0; i < t.nspec; i++)
			data[i] = FFTOps.get1DFFTMag(t.data[i]);
		
		PointSpectra all = new PointSpectra(data, newV, t.x, t.y);
		
		fftAvg = FFTOps.get1DFFTMag(t.average);
		
		PointSpectra avg = new PointSpectra(new double[][] {fftAvg}, newV, new double[] {0}, new double [] {0});
		
		return new PointSpectra[] {all, avg};
	}
	public static PointSpectra getFourierFiltered(PointSpectra t, int[] fmin, int[] fmax)
	{
		double[][] data = new double [t.nspec][t.nlayers];

		double[] fftAvg;
		
		double[] newV = ArrayOps.generateArrayInclBoth(0, t.nlayers-1, t.nlayers);
		
		for (int i = 0; i < t.nspec; i++)
			data[i] = FFTOps.getFourierFiltered(t.data[i], fmin, fmax);
		
		PointSpectra all = new PointSpectra(data, newV, t.x, t.y);
		return all;
	}

	public static PointSpectra[] subtBack(PointSpectra topo, int npoly)
	{
		Regression reg;
		double[] g = new double[topo.nlayers];
		double[] a;
		double[][] poly = new double [topo.nspec][topo.nlayers];
		double[][] subt = new double [topo.nspec][topo.nlayers];
		double[][] powers = new double[topo.nlayers][npoly+1];
		for (int k = 0; k < topo.nlayers; k++)
			for (int i = 0; i < npoly+1; i++)
				powers[k][i] = Math.pow(topo.v[k], i);
			
		
			for (int j = 0; j < topo.nspec; j++)
			{
				for (int k = 0; k < topo.nlayers; k++)
					g[k] = topo.data[j][k];
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
						poly[j][k] += a[m]*powers[k][m];
					subt[j][k] = topo.data[j][k]-poly[j][k];
				}
			}
		return new PointSpectra[] {new PointSpectra(subt, topo.v, topo.x, topo.y), new PointSpectra(poly, topo.v, topo.x, topo.y)};
	}
	
	public static PointSpectra getGaussSmoothed(PointSpectra t, double L)
	{
		return PointSpectra.newPointSpectra(t, gaussSmooth(t.data, L));
	}
	
	public static PointSpectra flipArrayOfSpectra(PointSpectra t)
	{
		double[][] ansdat = FieldOps.copy(t.data);
		ArrayOps.flipX(ansdat);
		double[] ansx = t.x.clone(), ansy = t.y.clone();
		
		ArrayOps.flip(ansx);
		ArrayOps.flip(ansy);
		
		return new PointSpectra(ansdat, t.v, ansx, ansy);
	}
	public static void doFlipArrayOfSpectra()
	{
		PointSpectra.writeBIN(flipArrayOfSpectra(PointSpectra.open(fc)), fc);
	}
	
	public static double[][] gaussSmooth(double[][] spec, double L)
	{
		int Lp = (int)(L+1);
		int biggerL = FieldOps.round(8*Lp);
		int nspec = spec.length, nlayers = spec[0].length;
		double[] gauss = new double [nlayers];
		double[][] smooth = new double[nspec][nlayers];
			for (int l = 0; l < nlayers; l++)
				gauss[l] = Math.exp(-(l*l)/(2*L*L));
		
		double gsum;
		int x, y;
		int xprime, yprime, xpmin, xpmax, ypmin, ypmax;
		for (int i = 0; i < nspec; i++){
			for (int j = 0; j < nlayers; j++)
			{
				gsum = 0;
				smooth[i][j] = 0;
				x = j;
				//this is done so that xprime runs from 0 to N, or from x - biggerL/2 to x+biggerL/2
				xpmin = Math.max(0, x - 6*Lp);
				xpmax = Math.min(nlayers, x + 6*Lp);
				for (xprime = xpmin; xprime < xpmax; xprime++)
				{
					gsum += gauss[Math.abs(x-xprime)];
					smooth[i][j] += gauss[Math.abs(x-xprime)]*spec[i][xprime];
				}
				smooth[i][j] /= gsum;
			}
		}
		return smooth;
	}

	
	public static void makePic(PointSpectra t)
	{
		SRAW.writeImage(FileOps.selectSave(fc).toString(), t.data);
	}
	public static void makePic(PointSpectra t, double minV, double maxV, int npts)
	{
		SRAW.writeImage(FileOps.selectSave(fc).toString(), getValuesBetweenIncl(t, minV, maxV, npts));
	}
	public static Layer getDataBetween(PointSpectra ivsK, double kmin, double kmax, int nk, int ne)
	{
		Layer t = ivsK.toLayer();
		double[] k = ArrayOps.generateArrayInclBoth(kmin, kmax, nk);
		double[] e = ArrayOps.generateArrayInclBoth(t.y[0], t.y[t.ny-1], ne);
		
		double[][] ans = new double [nk][ne];
		for (int i = 0; i < nk; i++)
			for (int j = 0; j < ne; j++)
				ans[i][j] = t.evaluateAtMetric(k[i], e[j]);
		return new Layer (ans, k, e, 1, 1);
	}
	public static double[][] getValuesBetweenIncl(PointSpectra t, double min, double max, int npts)
	{
		double[] v = ArrayOps.generateArrayInclBoth(min, max, npts);
		double[][] ans = new double [npts][t.nspec];
		
		for (int i = 0; i < npts; i++)
			for (int j = 0; j < t.nspec; j++)
				ans[i][j] = t.evaluateAt(v[i], j);
		
		return ans;
	}
	public static void doPeakFitting_MaximumCentroid(String dir, String ending, PointSpectra t, int radius, int nameIndex)//, boolean normalizeEachLine)
	{
		String outdir = dir + ending + "peakMaxCent_" + radius + "_" + nameIndex + "\\";// + (normalizeEachLine ? "norm_layer" : "") + "\\";
		String suboutdir = outdir + "Curves\\";
		File outd = new File (outdir);
		if (!outd.exists()) outd.mkdirs();
		File subd = new File (suboutdir);
		if (!subd.exists()) subd.mkdirs();
		
//		int layerpixheight = 6;
		
		double[][] data = t.data;
	
		//obtain E vs K.
		double[] k = t.v;
		int bmin = 2, bmax = data[0].length-1;
		double [] kmax = new double [t.y.length];
		for (int i = 0; i < t.y.length; i++)
		{
			kmax[i] = k[ArrayOps.maxIndex(data[i], bmin, bmax)];
		}
		int peakmin, peakmax;
		double[] peakRegion, peakRegionK;
		double[] centI = new double [t.y.length];
		for (int i = 0; i < t.y.length; i++)
		{
			peakmin = FieldOps.round(ArrayOps.maxIndex(data[i], bmin, bmax))-radius;
			peakmax = peakmin+2*radius+1;
			peakRegion = new double [peakmax-peakmin];
			peakRegionK = new double [peakmax-peakmin];
			
			for (int j = 0; j < peakRegion.length; j++)
			{
				peakRegion[j] = data[i][j+peakmin];
				peakRegionK[j] = k[j+peakmin];
			}
			
//			reg = new Regression(peakRegionK, peakRegion);
//			reg.linear();
			centI[i] = ArrayOps.indexToValue(peakRegionK, ArrayOps.centIndex(peakRegion));
			ColumnIO.writeTwoColumns(peakRegion, peakRegionK, suboutdir + "Curve " + i + "_" + t.y[i] + ".txt", "");
		}
		ColumnIO.writeTwoColumns(centI, t.y, "EvsK_scatter.txt", outdir);
		//Note re above: the Lorenztian is (Height/pi)*(Width/2)/((x-mean)^2 + (width/2)^2)
//		layerpixheight = (bmax-bmin)/(t.y.length+1);
//		System.out.println(layerpixheight);
	}
	public static PointSpectra averageEachPoint(PointSpectra[] splitByPosition)
	{
		int nspec = splitByPosition.length, nlayers = splitByPosition[0].nlayers;
		double[][] data = new double [nspec][nlayers];
		double[] x = new double [nspec];
		double[] y = new double [nspec];
		double[] v = splitByPosition[0].v;
		
		for (int i = 0; i < nspec; i++)
		{
			data[i] = splitByPosition[i].average;
			x[i] = splitByPosition[i].x[0];
			y[i] = splitByPosition[i].y[0];
		}
		return new PointSpectra(data, v, x, y);
	}
	public static PointSpectra cubicSpline(PointSpectra t, int fnpts)
	{
		double[][] fdata = new double [t.nspec][fnpts];
		CubicSplineWrapper wrap;
		double[] fv = ArrayOps.generateArrayInclBoth(t.v[0], t.v[t.v.length-1], fnpts);
		for (int i = 0; i < t.nspec; i++)
		{
			wrap = new CubicSplineWrapper(t.v, t.data[i]);
			for (int j = 0; j < fnpts; j++)
			{
				fdata[i][j] = wrap.of(fv[j]);
			}
		}
		return new PointSpectra(fdata, fv, t.x, t.y);
	}
	public static PointSpectra[] takeOnly(PointSpectra[] t, int imin, int imax)
	{
		PointSpectra[] result = new PointSpectra[imax-imin];
		for (int i = 0; i < result.length; i++)
			result[i] = t[i+imin];
		return result;
	}
	public static Layer convertToLayer(PointSpectra t)
	{
		return new Layer(FieldOps.transpose(t.data), t.v, t.x, 1, 1);
	}

	public static PointSpectra[] splitByMaxIntensityIndex(PointSpectra t)
	{
		int[] maxIndex = new int [t.nspec];
		for (int i = 0; i < t.nspec; i++)
			maxIndex[i] = ArrayOps.maxIndex(t.data[i]);
		int[][] hist = ArrayOps.getHistogram(maxIndex);
		
		for (int i = 0; i < hist[0].length; i++)
			System.out.println("" + t.v[hist[0][i]] + "\t" + hist[1][i]);
		
		int nbins = hist[0].length;
		PointSpectra[] ans = new PointSpectra[nbins];
		for (int i = 0; i < nbins; i++)
		{
			double[][] data = new double [hist[1][i]][t.nlayers];
			double[] x = new double [hist[1][i]];
			double[] y = new double [hist[1][i]];
			int n = 0;
			for (int j = 0; j < t.nspec; j++){
				if (maxIndex[j] == hist[0][i]){
					x[n] = t.x[j];
					y[n] = t.y[j];
					for (int k = 0; k < t.nlayers; k++)
						data[n][k] = t.data[j][k];
					n++;
				}
			}
			ans[i] = new PointSpectra(data, t.v, x, y);
			
		}
		return ans;
	}
	
	/**
	 * This assumes unique locations.
	 * @param t
	 * @return
	 */
	public static PointImp[] getImpList(PointSpectra t)
	{
		PointImp[] imps = new PointImp[t.nspec];
		for (int i = 0; i < t.nspec; i++)
		{
			imps[i] = new PointImp(new double[] {t.x[i], t.y[i]});
		}
		return imps;
	}
	
	public static PointSpectra getParabolaCentervsV(PointSpectra t, int di)
	{
		int npts = t.nlayers - di;
		double[][] data = new double [t.nspec][npts];
		double[][] yCalc = new double[npts][];
		double[][] xCalc = new double[npts][];
		for (int i = 0; i < t.nspec; i++)
			for (int j = 0; j < npts; j++)
			{
				FittingResult parab = new SpectraFittingRange(j, di, t.v, t.data[i]).fitToParabola();
				data[i][j] = parab.fitParams[0];
				yCalc[j] = parab.yCalc.clone();
				xCalc[j] = parab.x.clone();
			}
		double[] v = new double [npts];
		for (int i = 0; i < npts; i++){
			v[i] = t.v[i+di/2];
//			System.out.print(Printer.arrayLnHorizontal(yCalc[i]));
		}
		GraphDrawerCart.plotGraphs(xCalc, yCalc);
		
		
		
		return new PointSpectra(data, v, t.x, t.y);
	}
	
	public static PointSpectra replaceEachPixelGroupWithAverage(PointSpectra t)
	{
		PointSpectra[] split = t.splitByPosition();
		double[][] data = new double [split.length][t.nlayers];
		double[] x = new double [split.length];
		double[] y = new double [split.length];
		for (int i = 0; i < split.length; i++)
		{
			x[i] = split[i].x[0];
			y[i] = split[i].y[0];
			for (int j = 0; j < t.nlayers; j++)
				data[i][j] = split[i].average[j];
		}
		return new PointSpectra(data, t.v, x, y);
	}
	/***
	 * Here the line cuts are quasi-assumed to be taken in ascending order of fields. Each line cut must have the same number of position-split points, and should be taken in the same direction.
	 * @param lineCuts
	 * @return
	 */
	public static PointSpectra[] makeFieldDependenceDataset(PointSpectra[] lineCuts, double[] B)
	{
		ArrayList<PointSpectra[]> list = new ArrayList<PointSpectra[]>();
		for (int i = 0; i < lineCuts.length; i++)
			list.add(lineCuts[i].splitByPosition());
		int nlayers = lineCuts[0].nlayers;
		int npoints = list.get(0).length;
		double[] y = new double [B.length];
		double[] v = lineCuts[0].v;
		PointSpectra[] results = new PointSpectra[npoints+1];
		for (int i = 0; i < npoints; i++)
		{
			double[][] data = new double [lineCuts.length][nlayers];
			for (int j = 0; j < nlayers; j++)
				for (int k = 0; k < lineCuts.length; k++)
					data[k][j] = list.get(k)[i].average[j];
			results[i] = new PointSpectra(data, v, B, y);
		}
		
		double[][] average = new double [lineCuts.length][nlayers];
		for (int i = 0; i < lineCuts.length; i++)
			for (int j = 0; j < nlayers; j++)
			{
				for (int k = 0; k < npoints; k++)
					average[i][j] += list.get(i)[k].average[j];
				average[i][j] /= npoints;
			}
		results[npoints] = new PointSpectra(average, v, B, y);
		return results;
	}
	public static PointSpectra[] makeFieldDependenceDataset(LandauLevelFittingFolder f)
	{
		String[] paths = new String [f.names.length];
		for (int i = 0; i < paths.length; i++)
			paths[i] = f.dir + f.names[i] + ".bin";
		
		PointSpectra[] lineCuts = autoOpen(paths);
		return makeFieldDependenceDataset(lineCuts, f.B);
	}
	
	
	public static PointSpectra subtractFirstSpectrum(PointSpectra t)
	{
		int nlayers = t.nlayers;
		int nspec = t.nspec;
		double[][] data = new double [nspec][nlayers];
		for (int i = 1; i < nspec; i++)
			for (int j = 0; j < nlayers; j++)
				data[i][j] = t.data[i][j] - t.data[0][j];
		return PointSpectra.newPointSpectra(t, data);
	}
	public static double[][] adaptivelyFitAllSpectraToAParabola(PointSpectra t, int jstart, int jend, double initCent, double hWidth)
	{
		int nused = Math.abs(jstart-jend) + 1;
		double[] maxPos = new double [nused];
		double[] maxCurv = new double [nused];
		double[] maxRaw = new double [nused];
		double[] maxOffset = new double [nused];
		SpectraUtil.SpectraFittingRange sfr = new SpectraUtil.SpectraFittingRange(initCent-hWidth, initCent+hWidth, t.v, t.data[jstart]);
		FittingResult[] r = new FittingResult[nused];
		r[0] = sfr.fitToParabola();
		maxPos[0] = r[0].fitParams[0];
		maxCurv[0] = r[0].fitParams[1];
		maxRaw[0] = r[0].x[ArrayOps.maxIndex(r[0].y)];
		maxOffset[0] = r[0].fitParams[2];
		boolean tryingToGetMinimum = maxCurv[0] > 0;
		
		int delta = jstart > jend ? -1 : 1;
		
		double[][] xs = new double [2*nused][], ys = new double [2*nused][];
		
		double dh = (ArrayOps.max(r[0].y) - ArrayOps.min(r[0].y))/2;
		System.out.println(dh);
		xs[0] = r[0].x;
		xs[1] = r[0].x;
		ys[0] = r[0].y;
		ys[1] = r[0].yCalc;
		for (int j = jstart+delta; j != jend + delta; j += delta)
		{
//			if (counter % 10000 == 0) System.out.println(counter + "\t" + queue.size());
			int i = Math.abs(j - jstart);
			double guess = maxPos[j-delta];
			sfr = new SpectraUtil.SpectraFittingRange(guess-hWidth, guess+hWidth, t.v, t.data[j]);
			r[i] = sfr.fitToParabola();
			maxPos[i] = r[i].fitParams[0];
			maxCurv[i] = r[i].fitParams[1];
			maxRaw[i] = r[i].x[tryingToGetMinimum ? ArrayOps.maxIndex(r[i].y) : ArrayOps.maxIndex(r[i].y)];
			maxOffset[i] = r[i].fitParams[2];
//			System.out.print("" + points[2] + "\t" + points[3] + "\t" + Printer.arrayLnHorizontal(sfr.fitY));
			xs[2*i] = r[i].x;
			xs[2*i+1] = r[i].x;
			ys[2*i] = r[i].y;
			ys[2*i+1] = r[i].yCalc;
			ys[2*i] = ArrayOps.add(ys[2*i], i*dh);
			ys[2*i+1] = ArrayOps.add(ys[2*i+1], i*dh);
			System.out.println(i + " of " + nused);
		}
		
		GraphDrawerCart.plotGraphs(xs, ys);
		return new double[][] {maxPos, maxCurv, maxRaw, maxOffset};
		
	}

	
	public static PointSpectra[] autoOpen(String[] files)
	{
		PointSpectra[] ans = new PointSpectra[files.length];
		for (int i = 0; i < files.length; i++){
			ans[i] = PointSpectra.open(files[i]);
			if (ans[i].v[0] > ans[i].v[ans[i].v.length-1]) ans[i].flipV();
		}
		return ans;
	}
	public static PointSpectra[] autoOpen_map(String[] files, double smoothlength)
	{
		PointSpectra[] ans = new PointSpectra[files.length];
		for (int i = 0; i < files.length; i++){
			ans[i] = TopomapUtil.getTimeOrderedPointSpectra(Topomap.readBIN(files[i].toString()), 0);
			if (smoothlength > 0)
				for (int j = 0; j< ans[i].nspec; j++)
				{
					ans[i].data[j] = ArrayOps.gaussSmooth(ans[i].data[j], smoothlength);
					
				}
			if (ans[i].v[0] > ans[i].v[ans[i].v.length-1]) ans[i].flipV();
		}
		return ans;
	}
	
	public static void showLineCutAverages(PointSpectra t)
	{
		PointSpectra average  = replaceEachPixelGroupWithAverage(t);
		new LineCutDrawer(average);
	}
	
	public static PointSpectra getFromFieldDependenceDataset(LandauLevelFittingFolder f, int degree)
	{
		String[] paths = new String [f.names.length];
		for (int i = 0; i < paths.length; i++)
			paths[i] = f.dir + f.names[i] + (f.textDir ? ".txt" : ".bin");
		
		PointSpectra[] lineCuts = autoOpen(paths);
		PointSpectra[] fieldDep = makeFieldDependenceDataset(lineCuts, f.B);
		LineCutDrawer lc0 = new LineCutDrawer(fieldDep[fieldDep.length-1]);
		lc0.fc = fc;
		PointSpectra mZero = subtractFirstSpectrum(fieldDep[fieldDep.length-1]);
		
		PointSpectra[] subBack;
		if (degree > 0)
			subBack = subtBack(mZero, degree);
		else
			subBack = new PointSpectra[] {mZero, mZero};
		System.out.println(subBack[0].v[1] - subBack[0].v[0]);
		LineCutDrawer lc1 = new LineCutDrawer(subBack[0]);
		lc1.fc = fc;
		LineCutDrawer lc2 = new LineCutDrawer(subBack[1]);
		lc2.fc = fc;
		
//		PointSpectra parabFit = getParabolaCentervsV(subBack[0], 12);
//		new LineCutDrawer(parabFit);
		return subBack[0];
	}
	public static PointSpectra getFromFieldDependenceDataset_Map(LandauLevelFittingFolder f, int degree)
	{
		String[] paths = new String [f.names.length];
		for (int i = 0; i < paths.length; i++)
			paths[i] = f.dir + f.names[i] + (f.textDir ? ".txt" : ".bin");
		
		PointSpectra[] lineCuts = autoOpen_map(paths, 2);
		PointSpectra[] fieldDep = makeFieldDependenceDataset(lineCuts, f.B);
		PointSpectra mZero = subtractFirstSpectrum(fieldDep[fieldDep.length-1]);
		
		PointSpectra[] subBack;
		if (degree > 0)
			subBack = subtBack(mZero, degree);
		else
			subBack = new PointSpectra[] {mZero, mZero};
		System.out.println(subBack[0].v[1] - subBack[0].v[0]);
		new LineCutDrawer(subBack[0]);
		new LineCutDrawer(subBack[1]);
		
//		PointSpectra parabFit = getParabolaCentervsV(subBack[0], 12);
//		new LineCutDrawer(parabFit);
		return subBack[0];
	}
	
	public static LandauLevelPeak[][] getFromACertainFolder(LandauLevelFittingFolder f, int polyDegree, double[] excludedIndices){
		
		int nUnusedIndices = f.B.length-f.Btext.length;
		int offsetLL = 0;
		PointSpectra fieldDep = getFromFieldDependenceDataset(f, polyDegree);
		double[][] usedData = new double [f.namesT.length][fieldDep.nlayers];
		for (int i = 0; i < usedData.length; i++)
			for (int j = 0; j < fieldDep.nlayers; j++)
				usedData[i][j] = fieldDep.data[i+nUnusedIndices][j];
		
		PointSpectra usedF = new PointSpectra(usedData, fieldDep.v, f.Btext, new double [f.Btext.length]);
		
		
		String[] paths = new String [f.namesT.length];
		for (int i = 0; i < paths.length; i++){
			paths[i] = f.maxDir + f.namesT[i] + ".txt";
			System.out.println(paths[i]);
		}
		
		double[][] firstMaxima = new double[paths.length][];
		double[][] adjustedMaxima = new double [paths.length][];
		double[][] preIndices = new double[paths.length][];
		
		ArrayList<LandauLevelPeak> peaks = new ArrayList<LandauLevelPeak>();
		double hwidth = f.window;
		if (!f.preIndexed){
			for (int i = 0; i < paths.length; i++){
				firstMaxima[i] = ColumnIO.readNColumns(paths[i], 1, 0)[0];
				adjustedMaxima[i] = firstMaxima[i].clone();
			}
		}
		else
		{
			for (int i = 0; i < paths.length; i++){
				double[][] temp =  ColumnIO.readNColumns(paths[i], 2, 0);
				preIndices[i] = temp[0];
				firstMaxima[i] = temp[1];
				adjustedMaxima[i] = firstMaxima[i].clone();
			}
			
		}

		int nmaxima = 0;
		for (int i = 0; i < paths.length; i++)
			for (int j = 0; j < adjustedMaxima[i].length; j++)
				nmaxima++;

		double[][] xSpec = new double [usedF.nspec+2*nmaxima][];
		double[][] ySpec = new double [usedF.nspec+2*nmaxima][];
		int index = 0;
		double dh = 0.05;
		double smallNum = dh/20;
		double specMax  = ArrayOps.max(usedF.data);
		double[] offset = ArrayOps.generateArray(0, dh, paths.length);
		for (int i = 0; i < paths.length; i++){
			xSpec[index] = usedF.v.clone();
			ySpec[index] = usedF.data[i].clone();
			ArrayOps.addEquals(ySpec[index], offset[i]);
			index++;
			for (int j = 0; j < adjustedMaxima[i].length; j++){
				SpectraFittingRange sfr = new SpectraFittingRange(firstMaxima[i][j]-hwidth, firstMaxima[i][j]+hwidth,usedF.v,usedF.data[i]);
				FittingResult parab = sfr.fitToParabola();
				adjustedMaxima[i][j] = parab.fitParams[0];
				peaks.add(new LandauLevelPeak(adjustedMaxima[i][j], f.preIndexed ? preIndices[i][j] : (j-offsetLL), f.B[i+nUnusedIndices]));
//				if (preIndices[i][j] == 0.5) //handle the mystery peak
//					peaks.get(peaks.size()-1).mysteryOffset = 0.5;
					
				peaks.get(peaks.size()-1).fittingRange = sfr;
				System.out.println(peaks.get(peaks.size()-1));
		
				xSpec[index] = sfr.fitX.clone();
				ySpec[index] = parab.yCalc.clone();
				ArrayOps.addEquals(ySpec[index], offset[i]);
				index++;
				xSpec[index] = new double [] {adjustedMaxima[i][j], adjustedMaxima[i][j]};
				ySpec[index] = new double [] {specMax+smallNum, specMax+3*smallNum};
				ArrayOps.addEquals(ySpec[index], offset[i]);
				index++;
			}
			System.out.print(Printer.arrayLnHorizontal(firstMaxima[i]));
			System.out.print(Printer.arrayLnHorizontal(adjustedMaxima[i]));
		}
		
		LandauLevelPeak[] pk = new LandauLevelPeak[peaks.size()];
		for (int i = 0; i < pk.length; i++)
		{
			pk[i] = peaks.get(i);
		}
		//Now we have the alleged adjusted maxima
		LandauLevelPeak[][] split = splitByIndex(pk, excludedIndices);
		String[] lines = Printer.getColumnSeries(xSpec, ySpec);
		String fff = f.maxDir;//FileOps.selectSave(fc);
		ColumnIO.writeLines(lines,fff.toString() + "table for Origin.txt");
		double[][][] fand = LandauLevelPeak.getOrderedPairsEvsB(split);
		String[] fanDiagramLines = Printer.getColumnSeries(fand[0], fand[1]);
		ColumnIO.writeLines(fanDiagramLines,fff.toString() + "EvsB.txt");
		double[][][] evssqrtnb = LandauLevelPeak.getOrderedPairs(split);
		String[] nblines = Printer.getColumnSeries(evssqrtnb[0], evssqrtnb[1]);
		ColumnIO.writeLines(nblines,fff.toString() + "EvssqrtnB.txt");
		double[][][] evsnb = LandauLevelPeak.getOrderedPairsEvsNB(split);
		String[] nblines2 = Printer.getColumnSeries(evsnb[0], evsnb[1]);
		ColumnIO.writeLines(nblines2,fff.toString() + "EvsnB.txt");
		double[][][] evssqrtnPlusHalfb = LandauLevelPeak.getOrderedPairs(split, 0.5);
		String[] nphblines = Printer.getColumnSeries(evssqrtnPlusHalfb[0], evssqrtnPlusHalfb[1]);
		ColumnIO.writeLines(nphblines,fff.toString() + "EvssqrtnPlusHalfb.txt");
		String[] nbelines = new String[pk.length];
		for (int i = 0; i < pk.length; i++)
			nbelines[i] = "" + pk[i].n + "\t" + pk[i].B + "\t" + pk[i].max;
		ColumnIO.writeLines(nbelines, fff.toString() + "Landau level list.txt");
		//Let us also save the differences in energy between adjacent LL's as a function of field:
		double[][][] differences_B = new double [2][split.length-1][];
		for (int i = 0; i < split.length-1; i++){
			int nl = Math.max(split[i].length, split[i+1].length);
			differences_B[0][i] = new double [nl];
			differences_B[1][i] = new double [nl];
			ArrayList<Double> tempX = new ArrayList<Double>();
			ArrayList<Double> tempY = new ArrayList<Double>();
			for (int j = 0; j < split[i].length; j++){
				double B = split[i][j].B;
				for (int k = 0; k < split[i+1].length; k++)
					if (B == split[i+1][k].B){
						tempX.add(B);
						tempY.add(split[i+1][k].max - split[i][j].max);
					}
			}
			differences_B[0][i] = new double [tempX.size()];
			differences_B[0][i] = new double [tempY.size()];
			for (int j = 0; j < tempX.size(); j++){
				differences_B[0][i][j] = tempX.get(j);
				differences_B[1][i][j] = tempY.get(j);
			}
		}
		lines = Printer.getColumnSeries(differences_B[0], differences_B[1]);
		ColumnIO.writeLines(lines,fff.toString() + "differences_vs_Field.txt");
		
		//		for (int i = 0; i < lines.length; i++)
//			System.out.println(lines[i]);
		
		//We must write the spectra

		return split;
	}
	public static LandauLevelPeak[][] getFromACertainFolder_Map(LandauLevelFittingFolder f, int polyDegree, double[] excludedIndices){
		
		int nUnusedIndices = f.B.length-f.Btext.length;
		int offsetLL = 0;
		PointSpectra fieldDep = getFromFieldDependenceDataset_Map(f, polyDegree);
		double[][] usedData = new double [f.namesT.length][fieldDep.nlayers];
		for (int i = 0; i < usedData.length; i++)
			for (int j = 0; j < fieldDep.nlayers; j++)
				usedData[i][j] = fieldDep.data[i+nUnusedIndices][j];
		
		PointSpectra usedF = new PointSpectra(usedData, fieldDep.v, f.Btext, new double [f.Btext.length]);
		
		
		String[] paths = new String [f.namesT.length];
		for (int i = 0; i < paths.length; i++)
			paths[i] = f.maxDir + f.namesT[i] + ".txt";
		
		double[][] firstMaxima = new double[paths.length][];
		double[][] adjustedMaxima = new double [paths.length][];
		double[][] preIndices = new double[paths.length][];
		
		ArrayList<LandauLevelPeak> peaks = new ArrayList<LandauLevelPeak>();
		double hwidth = 0.004;
		if (!f.preIndexed){
			for (int i = 0; i < paths.length; i++){
				firstMaxima[i] = ColumnIO.readNColumns(paths[i], 1, 0)[0];
				adjustedMaxima[i] = firstMaxima[i].clone();
			}
		}
		else
		{
			for (int i = 0; i < paths.length; i++){
				double[][] temp =  ColumnIO.readNColumns(paths[i], 2, 0);
				preIndices[i] = temp[0];
				firstMaxima[i] = temp[1];
				adjustedMaxima[i] = firstMaxima[i].clone();
			}
			
		}

		int nmaxima = 0;
		for (int i = 0; i < paths.length; i++)
			for (int j = 0; j < adjustedMaxima[i].length; j++)
				nmaxima++;

		double[][] xSpec = new double [usedF.nspec+2*nmaxima][];
		double[][] ySpec = new double [usedF.nspec+2*nmaxima][];
		int index = 0;
		double dh = 0.2;
		double smallNum = dh/20;
		double specMax  = ArrayOps.max(usedF.data);
		double[] offset = ArrayOps.generateArray(0, dh, paths.length);
		for (int i = 0; i < paths.length; i++){
			xSpec[index] = usedF.v.clone();
			ySpec[index] = usedF.data[i].clone();
			ArrayOps.addEquals(ySpec[index], offset[i]);
			index++;
			for (int j = 0; j < adjustedMaxima[i].length; j++){
				SpectraFittingRange sfr = new SpectraFittingRange(firstMaxima[i][j]-hwidth, firstMaxima[i][j]+hwidth,usedF.v,usedF.data[i]);
				FittingResult parab = sfr.fitToParabola();
				adjustedMaxima[i][j] = parab.fitParams[0];
				peaks.add(new LandauLevelPeak(adjustedMaxima[i][j], f.preIndexed ? (int)preIndices[i][j] : (j-offsetLL), f.B[i+nUnusedIndices]));
				peaks.get(peaks.size()-1).fittingRange = sfr;
				System.out.println(peaks.get(peaks.size()-1));
		
				xSpec[index] = sfr.fitX.clone();
				ySpec[index] = parab.yCalc.clone();
				ArrayOps.addEquals(ySpec[index], offset[i]);
				index++;
				xSpec[index] = new double [] {adjustedMaxima[i][j], adjustedMaxima[i][j]};
				ySpec[index] = new double [] {specMax+smallNum, specMax+5*smallNum};
				ArrayOps.addEquals(ySpec[index], offset[i]);
				index++;
			}
			System.out.print(Printer.arrayLnHorizontal(firstMaxima[i]));
			System.out.print(Printer.arrayLnHorizontal(adjustedMaxima[i]));
		}
		
		LandauLevelPeak[] pk = new LandauLevelPeak[peaks.size()];
		for (int i = 0; i < pk.length; i++)
		{
			pk[i] = peaks.get(i);
		}
		//Now we have the alleged adjusted maxima
		LandauLevelPeak[][] split = splitByIndex(pk, excludedIndices);
		String[] lines = Printer.getColumnSeries(xSpec, ySpec);
		for (int i = 0; i < lines.length; i++)
			System.out.println(lines[i]);

		return split;
	}
	
	public static void main(String[] args)
	{
		Topomap.setStdDir();
//		H:\Work computer\data\analysis\In doped Bi2Se3\3.75 percent\mar16th
		String input = "H:\\Work computer\\data\\analysis\\In doped Bi2Se3\\3.75 percent\\mar16th\\";
		String output = input + "sm4\\";
//		processAllInFolder(input, output);
		
		LandauLevelFittingFolder f2p5 = LandauLevelFittingFolder.getFeb3rd_Folder_2p5Percent();
//		LandauLevelFittingFolder f0 = LandauLevelFittingFolder.getMar3rd_Folder_0Percent();
//		LandauLevelFittingFolder f0p = LandauLevelFittingFolder.getMar13th_Folder_0Percent();
//		LandauLevelFittingFolder f0p2 = LandauLevelFittingFolder.getMar13th_Folder_0Percen_no2();
//		LandauLevelFittingFolder f1 = LandauLevelFittingFolder.get1PercentFolder();
//		LandauLevelFittingFolder f5 = LandauLevelFittingFolder.getRun377Folder_5Percent();
//		LandauLevelFittingFolder fh = 
//		LandauLevelFittingFolder f8 = LandauLevelFittingFolder.run_of_Feb_25_2015_8Percent();
//		LandauLevelFittingFolder f3p75 = LandauLevelFittingFolder.getMar16th_Folder_375Percent();
//		
//		PointSpectra ps = PointSpectra.open(fc);
//		PointSpectra.writeBIN(getSubset(ps, 3*ps.nspec/4, ps.nspec), fc);
//		fc = new JFileChooser(Topomap.stddir); 
		
//		PointSpectra[] ans = makeFieldDependenceDataset(f5);
//		LandauLevelPeak[][] split = getFromACertainFolder_Map(f8, 3, null);
//		LandauLevelPeak[][] split = getFromACertainFolder(f5, 6, null);
		LandauLevelPeak[][] split = getFromACertainFolder(f2p5, 5, null);
//		LandauLevelPeak[][] split = getFromACertainFolder(f1, 5, null);
//		LandauLevelPeak[][] split = getFromACertainFolder(f5, 3, null);
//		LandauLevelPeak[][] split = getFromACertainFolder(f0, 3, null);
//		LandauLevelPeak[][] split = getFromACertainFolder(f0p, 1, null);
//		LandauLevelPeak[][] split = getFromACertainFolder(f0p2, 1, null);
//		LandauLevelPeak[][] split = getFromACertainFolder(f3p75, 1, null);
		
//		new LineCutDrawer(ans[ans.length-1]);
			
//		PointSpectra fieldDep = getFromFieldDependenceDataset_Map(f5, 0 );
//		doAFolder(LandauLevelFittingFolder.getRun377Folder_5Percent(), 0, null);
//		doAFolder(LandauLevelFittingFolder.getRun377Folder_5Percent(), 0, null);
//		double[] x = data[1], y = data[0];
		
//		
//		LandauLevelPeak[][] split = splitByIndex(getFromASCII(FileOps.selectOpen(fc)), null);//new int [] {4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24});
//		Eq48Zhang zh = new Eq48Zhang();x
//		zh.fit(split);
//		new Eq48Zhang_SliderPanel_2(split, zh);

		//		ArrayOps.flip(x);
//		ArrayOps.flip(y);
//		GraphDrawerCart.plotGraph(x, y);
//		double[] fft = FFTOps.get1DFFTMag(y);
//		double[] newY = FFTOps.getFourierFiltered(y, new int[] {0}, new int [] {15});
//		double[] indices = ArrayOps.generateArrayNotInclUpper(0, x.length, x.length);
//		GraphDrawerCart.plotGraph(indices, fft);
//		GraphDrawerCart.plotGraph(x, newY);
//		GraphDrawerCart.plotGraph(indices, FFTOps.get1DFFTMag(newY));
		
		
//		doFlipArrayOfSpectra();
//		PointSpectra t = PointSpectra.open(fc);
//		PointSpectra[] split = splitByMaxIntensityIndex(t);
//		File f = FileOps.selectSave(fc);
//		for (int i = 0; i < split.length; i++){
//			GraphDrawerCart.plotGraphs(t.v, split[i].data);
//			PointSpectra.writeBIN(split[i], f.toString() + MovieMaker.fromInt(i));
//			PointImp.writeToFile(SpectraUtil.getImpList(split[i]), new File(f.toString() + "_imps_" + MovieMaker.fromInt(i)));
//		}
		
//		Layer.writeBIN(getDataBetween(t, 0.38e10, 0.49e10, 512, 512), fc);
//		PointSpectra[] split = t.splitByPosition();
//		PointSpectra[] lineCutFull = takeOnly(split, 64, split.length);
//		PointSpectra cut = averageEachPoint(lineCutFull);
//		PointSpectra t = PointSpectra.open("C:\\data\\analysis\\Cu Data\\267\\Intensity spectra I vs k cut.bin");
//		PointSpectra backsub= subtBack(cut, 5)[0];
//		SpectraDrawer s = new SpectraDrawer(cut, null);
//		Layer.writeBIN(convertToLayer(t), fc);
//		makePic(backsub);
		//		PointSpectra tc = cubicSpline(t, 2048);
//		System.out.println(FileOps.addExtraBackslashes(fc.getSelectedFile()));
////		FileOps.writeTableASCII(t.data);
//		makePic(truncate(t, 3.8e9, 4.8e9, false));
//		doPeakFitting_MaximumCentroid(fc.getCurrentDirectory().toString(), "\\centroid fits\\", tc, 64, 1);
//		Printer.printlnVertical(t.v);
//		Printer.printlnVertical(t.y);
//		JFileChooser fc = new JFileChooser();
//		fc.showSaveDialog(null);
//		SRAW.writeImage(fc.getSelectedFile().toString(), t.data);
		
		Topomap t = Topomap.open();
		double[] spec = t.getAverageSpectrum();
//				t.getSpectrum(300, 300);
		double[] is = ArrayOps.generateArray(0, 1, t.nlayers);
//		spec = ArrayOps.subtractPolynomialFit(is, spec, 1)[0];
//		spec = ArrayOps.subtractLineForPeriodicness(spec)[0];
		double[] isd = ArrayOps.generateArrayInclBoth(-10, is.length+10, 1024);
		double[] interp = ArrayOps.getInterpolatedFourier(spec, isd);
		GraphDrawerCart.plotGraphs(new double[][] {is,  isd}, new double[][] {spec, interp});
	}
	public static LandauLevelFittingFolder getBi2Te3Kane()
	{
//		PointSpectra t = PointSpectra.open(fc);
		
		String dir = "C:\\data\\analysis\\Bi2Te3\\LL eph coupling\\peak files\\";
		
		String[]  names =	{
				"3.5T", "4T", "4.5T", "5T", "5.5T","6T","6.5T","7T","7.4T"
			};
		double[] B = {3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.4};
		
		//Now, take only the spectra to be used:
		String[]  namesT=	{
				"3.5T", "4T", "4.5T", "5T", "5.5T","6T","6.5T","7T","7.4T"
			};
		String textDir = "C:\\data\\analysis\\Bi2Te3\\LL eph coupling\\peak files\\";
		double[] Btext = {3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7, 7.4};
		return new LandauLevelFittingFolder(dir, names, B, namesT,  Btext, false);
	}
	public static void doAFolder(LandauLevelFittingFolder f, int polyDegree, double[] excludedIndices)
	{
		LandauLevelPeak[][] split = getFromACertainFolder(f, polyDegree, excludedIndices);
		
		for (int i = 0; i < split.length; i++)
			for (int j = 0; j < split[i].length; j++)
				System.out.println(split[i][j]);
		
		double[][][] table = LandauLevelPeak.getOrderedPairs(split);
		GraphDrawerCart.plotGraphs(table[0], table[1]);
		LandauLevelPeak.printOrderedPairs(split);
		String[] lines = Printer.getColumnSeries(table[0], table[1]);
		for (int i = 0; i < lines.length; i++)
			System.out.println(lines[i]);
		
		System.out.println("Now gamma = 1/2:");
		double[][][] table2 = LandauLevelPeak.getOrderedPairs(split, 0.5);
		String[] lines2 = Printer.getColumnSeries(table2[0], table2[1]);
		for (int i = 0; i < lines.length; i++)
			System.out.println(lines2[i]);
		
		//		double[][] data = FileOps.openTable(fc);
		
		Eq48Zhang zh = new Eq48Zhang();
		zh.fit(split);
		
	}
	public static LandauLevelPeak[] getFromASCII(File f)
	{
		double[][] table = ColumnIO.readNColumns(f, 3, 1);
		LandauLevelPeak[] ans = new LandauLevelPeak[table[0].length];
		for (int i = 0; i < table[0].length; i++)
			ans[i] = new LandauLevelPeak(table[2][i], (int)table[1][i], table[0][i]);
		return ans;
	}
	public static LandauLevelPeak[][] splitByIndex(LandauLevelPeak[] pk)
	{
		ArrayList<Double> uniqueIndices = new ArrayList<Double>();
		for (int i = 0; i < pk.length; i++)
		{
			boolean old = false;
			for (int j = 0; j < uniqueIndices.size() && !old; j++)
				if (pk[i].n == uniqueIndices.get(j)) old = true;
			
			if (!old) uniqueIndices.add(pk[i].n);
		}
		for (int i = 0; i < uniqueIndices.size(); i++)
			System.out.println(uniqueIndices.get(i));
		
		int[] numbers = new int[uniqueIndices.size()];
		for (int i = 0; i < pk.length; i++)
			for (int j = 0; j < uniqueIndices.size(); j++)
				if (pk[i].n == uniqueIndices.get(j)) numbers[j]++;
		
		LandauLevelPeak[][] ans = new LandauLevelPeak[uniqueIndices.size()][];
		for (int i = 0; i < uniqueIndices.size(); i++)
		{
			int n = 0;
			ans[i] = new LandauLevelPeak[numbers[i]];
			for (int j = 0; j < pk.length; j++)
				if (pk[j].n == uniqueIndices.get(i))
					ans[i][n++] = pk[j];
		}
		return ans;
	}
	public static LandauLevelPeak[][] splitByIndex(LandauLevelPeak[] pk, double[] excludedIndices)
	{
		if (excludedIndices == null)
			excludedIndices = new double [] {};
		ArrayList<Double> uniqueIndices = new ArrayList<Double>();
		for (int i = 0; i < pk.length; i++)
		{
			boolean old = false;
			for (int j = 0; j < uniqueIndices.size() && !old; j++)
				if (pk[i].n == uniqueIndices.get(j)) old = true;
			
			if (!old && !ArrayOps.contains(excludedIndices, pk[i].n)) uniqueIndices.add(pk[i].n);
		}
		for (int i = 0; i < uniqueIndices.size(); i++)
			System.out.println(uniqueIndices.get(i));
		
		int[] numbers = new int[uniqueIndices.size()];
		for (int i = 0; i < pk.length; i++)
			for (int j = 0; j < uniqueIndices.size(); j++)
				if (pk[i].n == uniqueIndices.get(j)) numbers[j]++;
		
		LandauLevelPeak[][] ans = new LandauLevelPeak[uniqueIndices.size()][];
		for (int i = 0; i < uniqueIndices.size(); i++)
		{
			int n = 0;
			ans[i] = new LandauLevelPeak[numbers[i]];
			for (int j = 0; j < pk.length; j++)
				if (pk[j].n == uniqueIndices.get(i))
					ans[i][n++] = pk[j];
		}
		return ans;
	}
	public static void processAllInFolder(String dirIn, String dirOut)
	{
		if (!new File(dirOut).exists()) new File (dirOut).mkdir();
		
		File[] inhere = new File(dirIn).listFiles();
		for (int i = 0; i < inhere.length; i++)
		{
			if (inhere[i].toString().endsWith(".bin")){
				PointSpectra t = PointSpectra.open(inhere[i].toString());
				t.data = SpectraUtil.gaussSmooth(t.data, 4);
				PointSpectra.writeBIN(t, dirOut + inhere[i].getName());
			}
		}
	}
	public static class SpectraFittingRange
	{
		public double center, min, max;
		public double[] x, y;
		
		double[] fitY;
		double[] fitX;
		
		public SpectraFittingRange(double min, double max, double[] x,
				double[] y) {
			super();
			this.min = min;
			this.max = max;
			this.center = (min+max)/2;
			
			this.x = x;
			this.y = y;
		
			int npts = 0;
			for (int i = 0; i < x.length; i++)
				if (x[i] >= min && x[i] <= max) npts++;
			
			fitY = new double [npts];
			fitX = new double [npts];
			int n = 0;
			for (int i = 0; i < x.length; i++)
				if (x[i] >= min && x[i] <= max)
				{
					fitX[n] = x[i];
					fitY[n] = y[i];
					n++;
				}
		}
		/**
		 * Assumes imin + di <= x.length
		 * @param imin
		 * @param di
		 * @param x
		 * @param y
		 */
		public SpectraFittingRange(int imin, int di, double[] x,
				double[] y) {
			super();
			this.min = x[imin];
			this.max = x[imin+di-1];
			this.center = (min+max)/2;
			
			this.x = x;
			this.y = y;
		
			fitY = new double [di];
			fitX = new double [di];
			int n = 0;
			for (int i = 0; i < di; i++){
				fitX[i] = x[i+imin];
				fitY[i] = y[i+imin];
			}
		}
		
		public FittingResult fitToParabola()
		{
			return ACM_NonLinearFitter.fitToParabola(fitX, fitY);
		}

	}
	
	public static class LandauLevelPeak{
		public SpectraFittingRange fittingRange;
		public double max;
		public double n;
		public double B;
		//This is alsways zero except for the mystery peak where we make it 0.5
		public double mysteryOffset = 0;
		public LandauLevelPeak(double max, double n, double B) {
			super();
			this.max = max;
			this.n = n;
			this.B = B;
		}
		public String toString()
		{
			return "" + max + "\t" + n + "\t" + B;
		}
		
		/**
		 * returns {sign n * sqrt(nB), max};
		 * @return
		 */
		public double[] getOrderedPair()
		{
			return new double[] {(n >= 0 ? 1 : -1)*Math.sqrt(Math.abs(n+mysteryOffset)*B), max};
		}
		public double[] getOrderedPair(double gamma)
		{
			return new double[] {(n >= 0 ? 1 : -1)*Math.sqrt((Math.abs(n + gamma))*B), max};
		}
		
		public double[] getOrderedPairEvsB()
		{
			return new double[] {B, max};
		}
		public double[] getOrderedPairEvsNB()
		{
			return new double[] {n*B, max};
		}
		public static double[][][] getOrderedPairs(LandauLevelPeak[][] split)
		{
			double[][][] ans = new double [2][split.length][];
			for (int i = 0; i < split.length; i++)
			{
				ans[0][i] = new double[split[i].length];
				ans[1][i] = new double[split[i].length];
				for (int j = 0; j < split[i].length; j++)
				{
					double[] r = split[i][j].getOrderedPair();
					ans[0][i][j] = r[0];
					ans[1][i][j] = r[1];
				}
			}
			return ans;
			
		}
		public static double[][][] getOrderedPairsEvsB(LandauLevelPeak[][] split)
		{
			double[][][] ans = new double [2][split.length][];
			for (int i = 0; i < split.length; i++)
			{
				ans[0][i] = new double[split[i].length];
				ans[1][i] = new double[split[i].length];
				for (int j = 0; j < split[i].length; j++)
				{
					double[] r = split[i][j].getOrderedPairEvsB();
					ans[0][i][j] = r[0];
					ans[1][i][j] = r[1];
				}
			}
			return ans;
			
		}
		public static double[][][] getOrderedPairsEvsNB(LandauLevelPeak[][] split)
		{
			double[][][] ans = new double [2][split.length][];
			for (int i = 0; i < split.length; i++)
			{
				ans[0][i] = new double[split[i].length];
				ans[1][i] = new double[split[i].length];
				for (int j = 0; j < split[i].length; j++)
				{
					double[] r = split[i][j].getOrderedPairEvsNB();
					ans[0][i][j] = r[0];
					ans[1][i][j] = r[1];
				}
			}
			return ans;
			
		}
		public static double[][][] getOrderedPairs(LandauLevelPeak[][] split, double gamma)
		{
			double[][][] ans = new double [2][split.length][];
			for (int i = 0; i < split.length; i++)
			{
				ans[0][i] = new double[split[i].length];
				ans[1][i] = new double[split[i].length];
				for (int j = 0; j < split[i].length; j++)
				{
					double[] r = split[i][j].getOrderedPair(gamma);
					ans[0][i][j] = r[0];
					ans[1][i][j] = r[1];
				}
			}
			return ans;
			
		}
		public static double[][][] getOrderedPairsModel(LandauLevelPeak[][] split, LandauLevelModel m)
		{
			double[][][] ans = new double [2][split.length][];
			for (int i = 0; i < split.length; i++)
			{
				ans[0][i] = new double[split[i].length];
				ans[1][i] = new double[split[i].length];
				for (int j = 0; j < split[i].length; j++)
				{
					double[] r = split[i][j].getOrderedPair();
					ans[0][i][j] = r[0];
					ans[1][i][j] = m.energy(split[i][j].n, split[i][j].B);
				}
			}
			return ans;
			
		}
		public static void printOrderedPairs(LandauLevelPeak[][] split)
		{
			double[][][] ans = new double [2][split.length][];
			for (int i = 0; i < split.length; i++)
			{
				ans[0][i] = new double[split[i].length];
				ans[1][i] = new double[split[i].length];
				for (int j = 0; j < split[i].length; j++)
				{
					double[] r = split[i][j].getOrderedPair();
					System.out.print("" + split[i][j].n + "\t" + Printer.arrayLnHorizontal(r));
					ans[0][i][j] = r[0];
					ans[1][i][j] = r[1];
				}
			}
//			return ans;
			
		}
		
		public static LandauLevelPeak[] readFromTxt(File f){
			double[][] columns = ColumnIO.readNColumns(f, 3, 0);
			LandauLevelPeak[] pk = new LandauLevelPeak[columns[0].length];
			for (int i = 0; i < pk.length; i++)
				pk[i] = new LandauLevelPeak(columns[2][i], columns[0][i], columns[1][i]);
			return pk;
		}
	}
	
	public static class LandauLevelFittingFolder{
		public String dir, maxDir;
		public String[] names;
		public double[] B;
		//Now, take only the spectra to be used:
		public String[] namesT;
		
		public boolean textDir;
		public double[] Btext;
		public double window;
		boolean preIndexed;
		public boolean text = false;
		public boolean map = false;
		public LandauLevelFittingFolder(String dir, String[] names, double[] B,
				String[] namesT, double[] Btext, boolean preIndexed) {
			super();
			this.dir = dir;
			this.names = names;
			this.B = B;
			this.namesT = namesT;
			this.Btext = Btext;
			this.preIndexed = preIndexed;
			Topomap.stddir = dir;
		}
		public static LandauLevelFittingFolder run_of_Feb_25_2015_8Percent()
		{
			String dir = "H:\\Work computer\\data\\analysis\\In doped Bi2Se3\\nanonis rhk\\run 475 8 percent\\";
			String[] names = {"0.5T","3T","4T","5T","5.5T","6T","6.5T","7T","7.5T"};
			double[] B = {0.5, 3, 4, 5, 5.5, 6, 6.5, 7, 7.5};
			String[] namesT = {"3T","4T","5T","5.5T","6T","6.5T","7T","7.5T"};
			double[] Btext = {3, 4, 5, 5.5, 6, 6.5, 7, 7.5};
			boolean preIndexed = false;
			LandauLevelFittingFolder f = new LandauLevelFittingFolder(dir, names, B, namesT, Btext, preIndexed);
			f.maxDir = dir + "\\provisional maxima\\";
			f.textDir = false;
			f.map = true;
			f.window = 0.004;
			return f;
			
		}
		public static LandauLevelFittingFolder run_of_Oct_25_2014()
		{
			String dir = "C:\\data\\analysis\\In doped Bi2Se3\\new data on rhk\\set number 2 10 25 2014\\sm2\\";
			String[] names = {"0.5T","3.5T","4.5T","5T","5.5T","6T","6.5T","7T","7.5T"};
			double[] B = {0.5, 3.5, 4.5, 5, 5.5, 6, 6.5, 7, 7.5};
			String[] namesT = {"3.5T","4.5T","5T","5.5T","6T","6.5T","7T","7.5T"};
			double[] Btext = {3.5, 4.5, 5, 5.5, 6, 6.5, 7, 7.5};
			boolean preIndexed = true;
			LandauLevelFittingFolder f = new LandauLevelFittingFolder(dir, names, B, namesT, Btext, preIndexed);
			f.maxDir = "C:\\data\\analysis\\In doped Bi2Se3\\new data on rhk\\set number 2 10 25 2014\\sm2\\provisional maxima\\";
			f.textDir = false;
			f.window = 0.004;
			return f;
		}
		public static LandauLevelFittingFolder run_of_Nov_2_2014()
		{
			String dir = "G:\\Work Computer\\data\\analysis\\BTS\\LL sets\\set 2 11 2 2014\\sm1\\";
			String[] names = {"0.0T","0.5T","3.0T", "4.0T", "4.5T","5.0T","5.5T","6.0T", "6.5T","7.0T","7.5T"};
			double[] B = {0,0.5, 3.0, 4.0, 4.5, 5, 5.5, 6, 6.5,7, 7.5};
			String[] namesT = {"3.0T", "4.0T", "4.5T","5.0T","5.5T","6.0T", "6.5T","7.0T","7.5T"};
			double[] Btext = {3.0, 4.0, 4.5,5, 5.5,6, 6.5,7, 7.5};
			boolean preIndexed = true;
			LandauLevelFittingFolder f = new LandauLevelFittingFolder(dir, names, B, namesT, Btext, preIndexed);
			f.maxDir = "G:\\Work Computer\\data\\analysis\\BTS\\LL sets\\set 2 11 2 2014\\sm1\\provisional maxima\\full mystery 0 m 123\\";
			f.textDir = false;
			f.window = 0.004;
			return f;
		}
		public static LandauLevelFittingFolder getRun377Folder_5Percent()
		{
//			PointSpectra t = PointSpectra.open(fc);
			
			String dir = "H:\\Work Computer\\data\\analysis\\In doped Bi2Se3\\5 percent\\run 377\\sm2\\";
			
			String[] names = 
				{
					"0T","3T","4T","5T","5.5T","6T","6.5T","7T","7.5T"
				};
			double[] B = {0, 3, 4, 5, 5.5, 6, 6.5, 7, 7.5};
			
			//Now, take only the spectra to be used:
			String[]  namesT=	{
					"3T", "4T", "5T","5.5T","6T","6.5T","7T","7.5T"
				};
			String textDir = "H:\\Work Computer\\data\\analysis\\In doped Bi2Se3\\5 percent\\run 377\\human maxima simple\\";
			double[] Btext = {3, 4, 5, 5.5, 6, 6.5, 7, 7.5};
			LandauLevelFittingFolder f = new LandauLevelFittingFolder(dir, names, B, namesT, Btext, false);
//			return new LandauLevelFittingFolder(dir, names, B, namesT,  Btext, false);
			f.maxDir = textDir;
			f.textDir = false;
			f.window = 0.004;
			return f;
		}
		public static LandauLevelFittingFolder get1PercentFolder()
		{
//			PointSpectra t = PointSpectra.open(fc);
			
			String dir = "H:\\Work Computer\\data\\analysis\\In doped Bi2Se3\\Field dependence set\\sm2\\";
			
			String[] names = 
				{
					"0T","3T","5T","5.5T","6T","6.5T","7T","7.5T"
				};
			double[] B = {0, 3, 5, 5.5, 6, 6.5, 7, 7.5};
			
			//Now, take only the spectra to be used:
			String[] namesT = 			{
					"5T","5.5T","6T","6.5T","7T","7.5T"
				};
//			String textDir = "H:\\Work Computer\\data\\analysis\\In doped Bi2Se3\\Field dependence set\\provisional maxima excplicit 00\\";
//			String textDir = "H:\\Work Computer\\data\\analysis\\In doped Bi2Se3\\Field dependence set\\provisional maxima excplicit 00p5\\";
//			String textDir = "H:\\Work Computer\\data\\analysis\\In doped Bi2Se3\\Field dependence set\\provisional maxima explicit m2m10p1p2\\";
			String textDir = dir + "pm minus1\\";
			double[] Btext = {5, 5.5, 6, 6.5, 7, 7.5};
			LandauLevelFittingFolder f = new LandauLevelFittingFolder(dir, names, B, namesT, Btext, true);
//			return new LandauLevelFittingFolder(dir, names, B, namesT,  Btext, false);
			f.maxDir = textDir;
			f.textDir = false;
			f.window = 0.005;
			return f;
		}
		public static LandauLevelFittingFolder getFeb3rd_Folder_2p5Percent()
		{
//			PointSpectra t = PointSpectra.open(fc);
			
			String dir = "H:\\Work Computer\\data\\analysis\\In doped Bi2Se3\\2.5 percent\\Feb 3\\sm2\\";
			
			String[] names = 
				{
					"0T","3T","4.5T","5T","5.5T","6T","6.5T","7T","7.5T"
				};
			double[] B = {0, 3, 4.5, 5, 5.5, 6, 6.5, 7, 7.5};
			
			//Now, take only the spectra to be used:
			String[]  namesT=	{
					"3T", "4.5T", "5T","5.5T","6T","6.5T","7T","7.5T"
				};
//			String textDir = "H:\\Work Computer\\data\\analysis\\In doped Bi2Se3\\2.5 percent\\Feb 3\\sm2\\human maxima explicit indices 00\\";
//			String textDir = "H:\\Work Computer\\data\\analysis\\In doped Bi2Se3\\2.5 percent\\Feb 3\\sm2\\human maxima explicit m10p1p2\\";
			String textDir = "H:\\Work Computer\\data\\analysis\\In doped Bi2Se3\\2.5 percent\\Feb 3\\sm2\\human maxima naive\\";
			double[] Btext = {3, 4.5, 5, 5.5, 6, 6.5, 7, 7.5};
			LandauLevelFittingFolder f = new LandauLevelFittingFolder(dir, names, B, namesT, Btext, true);
			f.maxDir = textDir;
			f.textDir = false;
			f.window = 0.005;
			return f;
		}
		public static LandauLevelFittingFolder getMar3rd_Folder_0Percent()
		{
//			PointSpectra t = PointSpectra.open(fc);
			
			String dir = "H:\\Work computer\\data\\analysis\\In doped Bi2Se3\\0 percent r481\\linecuts\\sm4\\";
			
			String[] names = 
				{
					"0T","3T","4.5T","5.5T","6.5T","7.5T"
				};
			double[] B = {0, 3, 4.5, 5.5, 6.5, 7.5};
			
			//Now, take only the spectra to be used:
			String[]  namesT=	{
					"3T","4.5T","5.5T","6.5T","7.5T"
				};
//			String textDir = "H:\\Work Computer\\data\\analysis\\In doped Bi2Se3\\2.5 percent\\Feb 3\\sm2\\human maxima explicit indices 00\\";
//			String textDir = dir + "\\pm skipmystery\\";
//			String textDir = dir + "\\pm WEIRD m0p5\\";
			String textDir = dir + "\\pm naive\\";
//			String textDir = dir + "\\pm minus1\\";
			
			double[] Btext = {3, 4.5, 5.5, 6.5, 7.5};
			LandauLevelFittingFolder f = new LandauLevelFittingFolder(dir, names, B, namesT, Btext, true);
			f.maxDir = textDir;
			f.textDir = false;
			f.window = 0.004;
			return f;
		}
		public static LandauLevelFittingFolder getMar13th_Folder_0Percent()
		{
//			PointSpectra t = PointSpectra.open(fc);
			
			String dir = "H:\\Work computer\\data\\analysis\\In doped Bi2Se3\\0 percent r481\\second dataset\\sm4\\";
			
			String[] names = 
				{
					"0T","3T","4T","5T","6T","7T","7.5T"
				};
			double[] B = {0, 3, 4, 5, 6, 7,7.5};
			
			//Now, take only the spectra to be used:
			String[]  namesT=	{
					"3T","4T","5T","6T","7T","7.5T"
				};
//			String textDir = "H:\\Work Computer\\data\\analysis\\In doped Bi2Se3\\2.5 percent\\Feb 3\\sm2\\human maxima explicit indices 00\\";
			String textDir = dir + "pm naive\\";
			double[] Btext = {3, 4, 5, 6, 7, 7.5};
			LandauLevelFittingFolder f = new LandauLevelFittingFolder(dir, names, B, namesT, Btext, true);
			f.maxDir = textDir;
			f.textDir = false;
			f.window = 0.004;
			return f;
		}
		public static LandauLevelFittingFolder getMar13th_Folder_0Percen_no2()
		{
//			PointSpectra t = PointSpectra.open(fc);
			
			String dir = "H:\\Work computer\\data\\analysis\\In doped Bi2Se3\\0 percent r481\\third dataset\\sm4\\";
			
			String[] names = 
				{
					"0.5T","3T","6T","7.5T"
				};
			double[] B = {0.5, 3, 6,7.5};
			
			//Now, take only the spectra to be used:
			String[]  namesT=	{
					"3T","6T","7.5T"
				};
//			String textDir = "H:\\Work Computer\\data\\analysis\\In doped Bi2Se3\\2.5 percent\\Feb 3\\sm2\\human maxima explicit indices 00\\";
			String textDir = dir + "pm naive\\";
			double[] Btext = {3, 6, 7.5};
			LandauLevelFittingFolder f = new LandauLevelFittingFolder(dir, names, B, namesT, Btext, true);
			f.maxDir = textDir;
			f.textDir = false;
			f.window = 0.004;
			return f;
		}
		public static LandauLevelFittingFolder getMar16th_Folder_375Percent()
		{
//			PointSpectra t = PointSpectra.open(fc);
			
			String dir = "H:\\Work computer\\data\\analysis\\In doped Bi2Se3\\3.75 percent\\mar16th\\sm4\\";
			
			String[] names = 
				{
					"0.5T","7.5T"
				};
			double[] B = {0.5, 7.5};
			
			//Now, take only the spectra to be used:
			String[]  namesT=	{
					"7.5T"
				};
//			String textDir = "H:\\Work Computer\\data\\analysis\\In doped Bi2Se3\\2.5 percent\\Feb 3\\sm2\\human maxima explicit indices 00\\";
			String textDir = dir + "pm naive\\";
			double[] Btext = {7.5};
			LandauLevelFittingFolder f = new LandauLevelFittingFolder(dir, names, B, namesT, Btext, true);
			f.maxDir = textDir;
			f.textDir = false;
			f.window = 0.004;
			return f;
		}
	}

	
	
	
}
