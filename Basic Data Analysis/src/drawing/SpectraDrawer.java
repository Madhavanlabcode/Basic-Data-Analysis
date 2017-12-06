package drawing;

import image.ImageEditing;

import java.awt.Color;
import java.awt.image.BufferedImage;

import javax.swing.JFileChooser;
import javax.swing.JOptionPane;

import analysis_gui.SpectraViewerTool;

import util.ArrayOps;
import util.Complex;
import util.FieldOps;
import util.Printer;
import util.SpectraUtil;
import util.color.ColorScales;
import util.fileops.PointSpectra;
import util.fileops.Topomap;
import util.fourier.FFTOps;

public class SpectraDrawer implements KeyStrokeProccessor{
	
	
	public GraphDrawerCart g;
	PointSpectra t, oldT;
	PointSpectra lastT;
	
	PointSpectra[] fft = null;
	boolean drawFFT = false;
	
	int nspec;
	
	double dh; //height offset.
	double[] x; //energy numbers.
	public double[][] spec; //[nspec][t.nlayers]. This is nominally a copy of t.data, but modifiable.;
	boolean drawThin = true;
	
	boolean specFlipped = false;
	
	boolean drawingLines = true; //else, render the line cut as an image
	BufferedImage image = null;
	BufferedImage sizedImage = null;
	int fx, fy; //image size ratios

	
	Color def = Color.RED;
	
	public SpectraViewerTool svt = null;
	public static JFileChooser fc;
	
	int selectedSpec = 0;
	
	boolean paintingJustOne = false;
	private boolean showAverage;
	
	int colorScaleIndex = 0;
	boolean colorscaleInverted = false;
	
	public SpectraDrawer(PointSpectra t, SpectraViewerTool svt)
	{
		fc = new JFileChooser(Topomap.stddir);
		this.oldT = t.copy();
		this.lastT = t.copy();
		this.t = t;
		x = t.v;
		this.svt = svt;
		this.nspec = t.nspec;
//		nspec = points.length;
		spec = new double [nspec][t.nlayers];
		
		specFlipped = t.v[0] > t.v[t.nlayers-1];
		
		copySpec();
		setDefaultDH();
		g = new GraphDrawerCart("Spectra Viewer", true);
		g.extra = this;
		setUpDrawer();
		image = new BufferedImage(t.nlayers, nspec, BufferedImage.TYPE_INT_RGB);
		resizeImage();
		g.showWindow();
	}
	public void resetTotally(PointSpectra t)
	{
		
		this.t = t;
		x = t.v;
//		this.svt = svt;
		this.nspec = t.nspec;
//		nspec = points.length;
		spec = new double [nspec][t.nlayers];
		
		specFlipped = t.v[0] > t.v[t.nlayers-1];
		
		copySpec();
		setDefaultDH();
		setUpDrawer();
	}
	
	public void changeNspec(int nspec)
	{
		this.nspec = nspec;
		spec = new double [nspec][t.nlayers];
		copySpec();
		this.setDefaultDH();
		setUpDrawer();
	}
	public void setDefaultDH()
	{
//		dh = (ArrayOps.max(spec) - ArrayOps.min(spec))/nspec;
		dh = 0;
	}
	public void copySpec()
	{
		for (int i = 0; i < nspec; i++)
			for (int j = 0; j < t.nlayers; j++)
				spec[i][j] = t.data[i][j] + i*dh;
	}
	public void setUpDrawer()
	{
		double max = ArrayOps.max(spec), min = ArrayOps.min(spec);
		g.setNumPlots(nspec+1);
		for (int i = 0; i < nspec; i++)
		{
			g.setXY(new GraphDrawerCart.GraphObject(x, spec[i]), false, false, i);
			g.setColor(def, i);
			g.setThin(drawThin, i);
			g.plot[i].isDrawing = !showAverage;
		}
		if (!drawFFT)
			g.setXY(new GraphDrawerCart.GraphObject(x, t.average), false, false, nspec);
		else
			g.setXY(new GraphDrawerCart.GraphObject(x, fft[1].average), false, false, nspec);
		g.setColor(def, nspec);
		g.setThin(drawThin, nspec);
		g.plot[nspec].isDrawing = showAverage;
		
		if (!specFlipped)
			g.setDomain(t.v[0], t.v[t.nlayers-1]);
		else
			g.setDomain(t.v[t.nlayers-1], t.v[0]);
			
		g.setRange(min, max);
		g.repaint();
	}
	
	public void refresh()
	{
		copySpec();
		if (!drawingLines)
			refreshImage();
		g.repaint();
	}
	public void dispose()
	{
		g.dispose();
	}
	public void setCoordinates(int x, int y, int w, int h)
	{
		g.setLocation(x, y);
		g.setSize(w, h);
	}
	public void processKeyStroke(char ch) {
		// TODO Auto-generated method stub
//		if (ch == 'n') {
//			this.changeNspec(Integer.parseInt(JOptionPane.showInputDialog("Enter the desired number of spectra")));
//		}
		if (ch == 'w')
		{
			this.dh = Double.parseDouble(JOptionPane.showInputDialog("Enter the desired offset in 'arbitrary units'."));
			copySpec();
			double max = ArrayOps.max(spec), min = ArrayOps.min(spec);
			g.setRange(min, max);
			g.repaint();
		}
		if (ch == 't'){
//			this.dh = Double.parseDouble(JOptionPane.showInputDialog("Enter the desired offset in 'arbitrary units'."));
//			double max = ArrayOps.max(spec), min = ArrayOps.min(spec);
//			g.setRange(min, max);
//			g.repaint();
			double min = Double.parseDouble(JOptionPane.showInputDialog("Enter the new minimum. (Old minimum " + ArrayOps.min(x) + ")"));
			double max = Double.parseDouble(JOptionPane.showInputDialog("Enter the new maximum. (Old maximum " + ArrayOps.max(x) + ")"));

			lastT = t.copy();
			PointSpectra t = SpectraUtil.truncate(this.t, min, max, specFlipped);
			this.resetTotally(t);
		}
		if (ch == 'T'){//Truncate by index
//			this.dh = Double.parseDouble(JOptionPane.showInputDialog("Enter the desired offset in 'arbitrary units'."));
//			double max = ArrayOps.max(spec), min = ArrayOps.min(spec);
//			g.setRange(min, max);
//			g.repaint();
			int min =Integer.parseInt(JOptionPane.showInputDialog("Enter the new minimum index"));
			int max = Integer.parseInt(JOptionPane.showInputDialog("Enter the number of data points to take off of the high end."));

			lastT = t.copy();
			PointSpectra t = SpectraUtil.truncate(this.t, min, (this.t.nlayers-max), specFlipped);
			this.resetTotally(t);
		}
		if (ch == 'i'){
			drawThin = !drawThin;
			for (int i = 0; i < nspec; i++)
			{
				g.setThin(drawThin, i);
				g.repaint();
			}
		}
		if (ch == 'f')
		{
			drawFFT = !drawFFT;
			if (drawFFT)
			{
				lastT = t.copy();
				fft = SpectraUtil.getFFT(t);
				this.resetTotally(fft[0]);
			}
			else
			{
				this.resetTotally(lastT);
			}
		}
		if (ch == 'L')
		{
			lastT = t.copy();
			int[][] bounds = this.getFourierFilterBounds();
			this.resetTotally(SpectraUtil.getFourierFiltered(t, bounds[0], bounds[1]));
		}
		if (ch == 'I')
		{
			if (!drawFFT)
				switchLineDrawing();
			else
			{
				int[][] boundArray = this.getFourierFilterBounds();
				if (boundArray == null) return;
				this.doIFFTAveragePrintFiltered(boundArray[0], boundArray[1]);
				this.doIFFTAveragePrintFiltered(boundArray[0], boundArray[1]);
 			}
		}
		if (ch == 's') //save everything
		{
			PointSpectra.writeBIN(t, fc);
		}
		if (ch == 'u')
		{
			String input = JOptionPane.showInputDialog("Enter the limits, comma-separated");
			int imin = Integer.parseInt(input.split(",")[0]);
			int imax = Integer.parseInt(input.split(",")[1]);
			this.resetTotally(SpectraUtil.getSubset(t, imin, imax));
		}
		if (ch == 'o')
		{
			PointSpectra p = PointSpectra.open(fc);
			this.oldT = p;
			this.lastT = p;
			this.resetTotally(p);
		}
//		if (ch)
		if (ch == 'a')
		{
			selectedSpec--;
			selectedSpec += nspec;
			selectedSpec %= nspec;
			resetSelection();
		}
		if (ch == 'd')
		{
			selectedSpec++;
			selectedSpec += nspec;
			selectedSpec %= nspec;
			resetSelection();
		}
		if (ch == 'k')
		{
			lastT = t.copy();
			PointSpectra p = SpectraUtil.deleteOneSpectrum(t, selectedSpec);
			this.resetTotally(p);
		}
		if (ch == 'e')
		{
			double y = g.currentMouseData[1];
			int index = 0;
			double d = Double.MAX_VALUE;
			
			System.out.println(Printer.vectorP(g.currentMouseData));
			for (int i = 0; i < t.nspec; i++)
			{
				if (Math.abs(y - t.evaluateAt(g.currentMouseData[0], i)) < d)
				{
					d = Math.abs(y - t.evaluateAt(g.currentMouseData[0], i));
					System.out.println("" + d +"\t" + i);
					index = i;
				}
			}
			selectedSpec = index;
			resetSelection();
		}
		if (ch == 'p')
		{
			int n = Integer.parseInt(JOptionPane.showInputDialog(null, "Enter the degree of polynomial to subtract"));
			if (n > 0)
			{
				PointSpectra[] fit = SpectraUtil.subtBack(t, n);
				lastT = t.copy();
				this.resetTotally(fit[0]);
			}
		}
		if (ch == 'P')
		{
//			g.printAll(0, t.nspec);
			if (showAverage)
				for (int i = 0; i < t.nlayers; i++)
				{
					System.out.println("" + t.v[i] + "\t" + t.average[i]);
				}
			else
				PointSpectra.writeASCII(t, fc);
			
		}
		if (ch == 'R')
		{
			resetTotally(SpectraUtil.collapseTimeWise(t, Integer.parseInt(JOptionPane.showInputDialog("Enter the averaging factor"))));
		}
//		if (ch == 'P')
//		{
//			lastT = t.copy();
//			for (int i = 0; i < t.nspec; i++)
//				for (int j = 0; j < t.nlayers; j++)
//					t.data[i][j] -= t.v[j];
//			
//			for (int j = 0; j < t.nlayers; j++)
//				t.average[j] -= t.v[j];
//			this.resetTotally(t);
//			
//		}
		if (ch == 'z')
		{
			this.resetTotally(oldT.copy());
		}
		if (ch == 'Z')
		{
			this.resetTotally(lastT);
		}
		if (ch == 'n')
		{
			paintingJustOne = !paintingJustOne;
			resetSelection();
		}
		if (ch == 'h')
		{
			double smoothness = Double.parseDouble(JOptionPane.showInputDialog("Enter smoothing length scale (pixels)."));
			lastT = t.copy();
			resetTotally(SpectraUtil.getGaussSmoothed(t, smoothness));
		}
		if (ch == 'M')
		{
			String[] tok = JOptionPane.showInputDialog("Enter the new bounds of the dataset.").split(",");
			double min = Double.parseDouble(tok[0]);
			double max = Double.parseDouble(tok[1]);
			lastT = t.copy();
			FieldOps.cutoffExtremesValue(t.data, min, max);
			t.resetAverage();
			this.resetTotally(t);
		}
		if (ch == 'c')
		{
			showAverage = !showAverage;
			System.out.println();
			System.out.println();
			System.out.println(Printer.arrayVertical(t.average));
			setUpDrawer();
			
		}
		if (ch == '.')
		{
			colorScaleIndex++;
			colorScaleIndex += ColorScales.NSCALES;
			colorScaleIndex %= ColorScales.NSCALES;
			if (!drawingLines)
			{
				g.imageScale = ColorScales.getNew(spec, colorScaleIndex, colorscaleInverted);
				g.refreshImage();
				refresh();
			
			}
		}
		if (ch == ',')
		{
			colorScaleIndex--;
			colorScaleIndex += ColorScales.NSCALES;
			colorScaleIndex %= ColorScales.NSCALES;
			if (!drawingLines)
			{
				g.imageScale = ColorScales.getNew(spec, colorScaleIndex, colorscaleInverted);
				g.refreshImage();
				refresh();
			
			}
		}
		if (ch == '-')
		{
			colorscaleInverted = !colorscaleInverted;
		}
		if (ch == ' ')
		{
			//print the average spectrum. 
			System.out.println();
			for (int i = 0; i < t.nlayers; i++)
				System.out.println(t.average[i]);
		}
		System.out.print("" + ch + ' ');
		
	}
	private int[][] getFourierFilterBounds()
	{
		String input = JOptionPane.showInputDialog("Enter the frequencies you want to filter in the following form:\r\n" +
				"A series of integer ranges a-b, separated by commas.\r\n" +
				"if there is only one range, omit any commma.\r\n" + 
				"a must be less than OR EQUAL TO b. If the range is to consist of just one frequency, do not omit the dash but simply let a and b be the same.\r\n" +
				"Both a and be will be included in the excluded frequencies.\r\n" +
				"You do NOT (repeat NOT) have to include the complementary range (i.e. the mirror reflection of a-b about the halfway line).\r\n");
		if (input.trim().equals(""))
			return null;
		else if (!input.contains(","))
		{
			String[] tok = input.split("-");
			int fmin = Integer.parseInt(tok[0]);
			int fmax = Integer.parseInt(tok[1]);
			return new int[][] {new int[] {fmin}, new int[] {fmax}};
		}
		else
		{
			String[] ranges = input.split(",");
			int[] fmin = new int [ranges.length];
			int[] fmax = new int [ranges.length];
			for (int i = 0; i < ranges.length; i++)
			{
				String[] tok = ranges[i].split("-");
				fmin[i] = Integer.parseInt(tok[0]);
				fmax[i] = Integer.parseInt(tok[1]);
			}
			return new int[][] {fmin, fmax};
		}

	}
	private void doIFFTAveragePrintFiltered(int[] fmin, int[] fmax) {
		double[] values = lastT.average;
		double magmin = ArrayOps.min(FFTOps.get1DFFTMag(values));
		double[] fftz = FFTOps.get1DFFTComplex(values);
		double[] ifftReal = new double [t.nlayers];
		//We have to assume here that lastT is the spectrum whose average we desire to modify.
		for (int k = 0; k < fmin.length; k++){
			if (fmax[k] >= t.nlayers/2) {
				System.out.println("Error. " + fmax + " was outside of the range. Returning original values.");
				FFTOps.putIFFTReal(fftz.clone(), ifftReal);
	//				return;
			}
	
			
			int f2max = t.nlayers-1-fmax[k];
			int f2min = t.nlayers-1-fmin[k];
			double ratio;
			for (int i = 0; i < values.length; i++){
				if ((i >= fmin[k] && i <= fmax[k]) || (i >= f2max && i <= f2min))
				{
					ratio = Complex.mag(fftz[2*i], fftz[2*i+1])/magmin;
					fftz[2*i] /= ratio;
					fftz[2*i+1] /= ratio;
	//					System.out.println(i + "\t" + ratio + "\t" + fmin + "\t" + fmax + "\t" + magmin + "\t" + Complex.mag(fftClone[2*i], fftClone[2*i+1]));
				}
			}
		}
		FFTOps.putIFFTReal(fftz, ifftReal);
		double[] fftmag = FFTOps.get1DFFTMag(ifftReal);
		double[] oldfftmag = FFTOps.get1DFFTMag(values);
		
		//Now, we print everthing: column 1 voltage. Column 2 IFFT. Column 3 FFT with suppressed shit.
		System.out.println("\r\nVoltage\tFiltered Data\tOld average\tFiltered FFT");
		for (int i = 0; i < t.nlayers; i++)
			System.out.println(lastT.v[i] + "\t" + ifftReal[i] + "\t" + values[i] + "\t "+ fftmag[i] + "\t" + oldfftmag[i]);
		
		//Just for the hell of it:
		GraphDrawerCart.plotGraph(lastT.v, ifftReal);

	}
	private void switchLineDrawing() {
		// TODO Auto-generated method stub
		drawingLines = !drawingLines;
		
		g.paintingArrays = drawingLines;
		g.imageScale = ColorScales.getNew(spec, colorScaleIndex);
		g.refreshImage();
		refresh();
	}
	private void refreshImage() {
		// TODO Auto-generated method stub
		image = ImageEditing.getBufferedImageTranspose(spec, 0);
		if (specFlipped) ImageEditing.flipX(image);
		try{
			ImageEditing.enlargeBasicStretch(image, sizedImage, fx, fy);
			
		}
		catch (Exception e)
		{
			e.printStackTrace();
		}
	}
	public void resizeImage() {
		int[] dim = g.getPlaneSize();
		fy = dim[1]/nspec;
		fx = dim[0]/t.nlayers;
		if (fx == 0) fx=1; if (fy==0) fy=1;
		sizedImage = new BufferedImage(fx*t.nlayers, fy*nspec, BufferedImage.TYPE_INT_RGB);
		g.image = sizedImage;
	}
	private void resetSelection() {
		// TODO Auto-generated method stub
		g.setTitle("Selected Spectrum: " + selectedSpec + " at (" + t.x[selectedSpec] + ", " + t.y[selectedSpec] + ")");
		resetColors();
		g.repaint();
	}
	public void resetColors()
	{
		for (int i = 0; i < nspec; i++)
			if (i == selectedSpec)
				g.setColor(Color.GREEN, i);
			else
				g.setColor(Color.RED, i);
		
		for (int i = 0; i < nspec; i++)
			g.plot[i].isDrawing = (i == selectedSpec || !paintingJustOne);
				
	}
	public static void main(String[] args)
	{
		Topomap.setStdDir();
		fc = new JFileChooser(Topomap.stddir);
		PointSpectra t = PointSpectra.open(fc);
		new SpectraDrawer(t, null);
	}
}
