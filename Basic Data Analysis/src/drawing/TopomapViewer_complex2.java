package drawing;

import image.ImageEditing;

import java.awt.Color;
import java.awt.Graphics;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;

import javax.swing.BoxLayout;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.border.LineBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import main.SRAW;
import schrodinger.MovieMaker;
import util.ArrayOps;
import util.Complex;
import util.ExportUtil;
import util.FieldOps;
import util.Printer;
import util.TopomapUtil;
import util.color.ColorScale1D;
import util.color.ColorScale2d;
import util.color.ColorScaleHolder;
import util.color.ColorScales;
import util.fileops.FileOps;
import util.fileops.Layer;
import util.fileops.RHKFileOps;
import util.fileops.Topomap;
import util.fourier.FFT3D_Wrapper;
import util.fourier.FFTOps;
import util.fourier.ImpurityListEditor;
import util.fourier.ImpurityListEditor.Impurity;
import util.geom.AtomicCoordinatesSet;
import util.geom.Distance;
import util.regression.FFTCutFitter;
import util.regression.FFTCutFitter.EQPair;


public class TopomapViewer_complex2 extends JFrame implements MouseListener, MouseMotionListener, MouseWheelListener, KeyListener{
	
	
	Topomap t, u;
	int N;
	double fmax, fmin, fdelta;
	double[][][] fftmag = null;
	boolean showFFT = false;
	double[][][] drawFieldC;
	
	GraphDrawerCart spec;
	double[] spectrum;
	double[] twoV;
	
	static JFileChooser fc;
	
	BufferedImage image;
	boolean refresh = true;
	
//	ColorScale scale;
	ColorScale2d cscale;
	double scalemin, scalemax;
	
	int defaultSize = 512;
	
	int sx = 512, sy = 512;
	int ox = 20, oy = 40;
	
	int[] writepoint = {ox + sx + 50, oy + 10};
	int linesize = 15;
	
	int zoomLevel = 0, zoomFactor = 1;
	int sizeratio = 1;
	
	int WIDTH = 1600, HEIGHT = 1080;
	
	//para is the current index
	int para = 0;
	double imageOverField = 1;

	int currentx, currenty;
	int currenti, currentj;
	double currentid, currentjd;
	double currentr, currentphi;
	int calcx = 0, calcy = 0;
	String dir;
	String name = null;
	double[] pbounds;
	SliderPanel s;
	int snpts;
	
	//line cut feature
	LineCutDrawer lc = null;
	StripCutTool sct = null;
	double[][] line = null; //the two points
	public int nearestEndIndex = 0;
	double d1, d2;
	
	//color scale change feater.
//	int currentCScale = 0;
	ColorScaleHolder csh;
	private boolean changeScaleWithLayer;
	private boolean ftlog = true;
	private boolean refreshFFT = false;
	
	double[] min, max, delta;
	double[] realmin, realmax, realdelta;
	double[] fftmin, fftmax, fftdelta;
	
	ArrayList<Impurity> imps = null;

	
	FFTCutFitter fitter;
	LayerViewer lv;
	public TopomapViewer_complex2(Topomap t, Topomap u, String dir,  int size)
	{
		snpts = t.nlayers-1;
		this.u = u;
		this.t = t;
		this.dir = dir;
		this.defaultSize = size;
		
		if (fc != null && fc.getSelectedFile() != null) name = fc.getSelectedFile().getName().substring(0, fc.getSelectedFile().getName().length()-4);
		fc = new JFileChooser(dir);
		drawFieldC = new double [t.nx][t.ny][2];
		spectrum = new double [t.nlayers];
		t.putSpectrum(0, 0, spectrum);
		spec = new GraphDrawerCart("Spectrum", t.v, spectrum);
		twoV = new double[] {t.v[0], t.v[0]};
		double min = ArrayOps.min(t.data);
		double max = ArrayOps.max(t.data);
		this.realmin = new double [t.nlayers];
		this.realmax = new double [t.nlayers];
		realdelta = new double[t.nlayers];
		this.fftmin = new double [t.nlayers];
		this.fftmax = new double [t.nlayers];
		fftdelta = new double[t.nlayers];

		this.min = realmin;
		this.max = realmax;
		this.delta = realdelta;
		recalculateBounds();
		cscale = new ColorScales.MYC2d(this.min[para] + delta[para]*0, this.min[para] + delta[para]*1, 2*Math.PI);
		delta = realdelta;
		csh = new ColorScaleHolder(min, max);
		spec.setXY(new GraphDrawerCart.GraphObject(twoV, new double[] {min, max}), false, false, 1);
		spec.setRange(min, max);
		spec.showWindow();
		
		N = t.nx;
		sx = N;
		resetDrawField();
		writepoint[0] = ox + sizeratio*sx + 50;
		WIDTH = sizeratio*t.nx + 450;
		HEIGHT = sizeratio*t.ny + 50;
		
		sy = this.drawFieldC[0].length;
		setFieldInfo();
		image = new BufferedImage(sx*sizeratio, sy*sizeratio, BufferedImage.TYPE_INT_RGB);
//		resetColorScale(ArrayOps.min(t.data), ArrayOps.max(t.data), true);
		formImage();
//		formFTImage();
		//gradcalc.activate();
		showWindow();
		setTitle("Title");

		s = new SliderPanel(this, new JFrame());
		s.show();
		addKeyListener(this);
		addMouseWheelListener(this);
		addMouseMotionListener(this);
		addMouseListener(this);
	}
	public void setFieldInfo()
	{
//		double[] bounds = calc.getMinMax(real);
		double[] bounds = null;
		bounds = new double[] {FieldOps.magMin(drawFieldC), FieldOps.magMax(drawFieldC)};
			
		fmax = bounds[1]; fmin = bounds[0];
		fdelta = fmax - fmin;
		setTitle("Range: [" + fmin + ", " + fmax + "]" + "     " + "p = " + para);
	}
	public void resetCSHFromSlider(double downnum, double upnum)
	{
		scalemax = min[para] + delta[para]*upnum;
		scalemin = min[para] + delta[para]*downnum;
		
		setTitle("Range: [" + scalemin + ", " + scalemax + "]" + "     " + "p = " + para);
//		csh.reBoundColorScale(scalemin, scalemax);
		cscale.renormalize(scalemax, scalemin);
	}
//	public void resetColorScale()
//	{
//		if (real) scale = ColorScales.getNew(scalemin, scalemax, currentCScale);
//		else cscale = new ColorScales.MYC2d(scalemax, scalemin, 2*Math.PI);
//	}

	public void formImage()
	{
		SRAW.writeImage(image, drawFieldC, cscale, sizeratio);
	}
	public void showWindow()
	{
		setSize(WIDTH, HEIGHT);
		repaint();
		show();
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}	
	public void paint(Graphics g)
	{
		if (refresh){
			g.clearRect(0, 0, 2000, 2000);
			g.drawImage(image, ox, oy, null);
			if (lc != null || sct != null) drawLine(g, csh.getUnusedColor());
			if (imps != null && !showFFT)
			{
				drawImpMarks(g, csh.getUnusedColor());
			}
			if (fitter != null)
			{
				double[] rOn = fitter.onLatt.getPixelCoords(new double[] {N/2, 0});
				double[] rOff = fitter.offLatt.getPixelCoords(new double[] {N/2, 0});
				g.setColor(csh.getUnusedColor());
				g.drawLine(screenX(N/2), screenY(N/2), screenX(rOn[0]), screenY(rOn[1]));
				g.setColor(Color.BLACK);
				g.drawLine(screenX(N/2), screenY(N/2), screenX(rOff[0]), screenY(rOff[1]));
				Printer.printlnHorizontal(rOn);
			
			}
			
			refresh = false;
			
		}
		g.setColor(Color.BLACK);
		drawText(g);
	}
	public void drawText(Graphics g)
	{
		g.clearRect(writepoint[0], writepoint[1]-linesize, 500, 800);
		String it = "";
		it += "Current mouse position in map: (" + currenti + ", " + currentj + ")" + "\r\n";
		it += "With respect to the center, that is (" + (currenti-t.nx/2) + ", " + (currentj-t.ny/2) + ") \r\n";
		if (!showFFT)
			it += "In metric units that is " + Printer.vectorP(t.getMetricCoords(currenti, currentj)) + "\r\n"; 
		else
			it += "In metric units that is " + Printer.vectorP(t.getFourierMetricCoords(currenti, currentj)) + "\r\n"; 
			
		it += "This vector has magnitude " + currentr + " \r\n";
		if (!showFFT)
			it += "(metric " + Distance.distance(t.getMetricCoords(t.nx/2, t.ny/2), t.getMetricCoords(currenti, currentj)) + ")\r\n";
		else
			it += "(metric " + Complex.mag(t.getFourierMetricCoords(currenti, currentj)) + ")\r\n";
			
		it += "And corresponding angle " + Math.toDegrees(FieldOps.atan((currenti-t.nx/2), (currentj-t.ny/2))) + "degrees. \r\n";
		it += "The value of the drawn field there is " + Complex.mag(drawFieldC[currenti][currentj]) + ".\r\n";
		it += "With phase " + Complex.phase(drawFieldC[currenti][currentj]);
		if (lc != null)
			it += "The line is from " + Printer.vectorP(line[0]) + " to " + Printer.vectorP(line[1]) + "\r\n";
		String[] lines = it.split("\r\n");
		for (int i = 0; i < lines.length; i++)
			g.drawString(lines[i], writepoint[0], writepoint[1] + i*linesize);
	}
	public void drawImpMarks(Graphics g, Color c)
	{
		double[] r = new double[2];
		g.setColor(c);
		for (int i = 0; i < imps.size(); i++)
		{
			r = imps.get(i).getPosition();
			drawPlus(g, (int)(r[0]*sizeratio)+ox, (int)(r[1]*sizeratio)+oy, 5);
		}
	}
	public void drawLine(Graphics g, Color c)
	{
		g.setColor(c);
		g.drawString("A", (int)(line[0][0]*sizeratio)+ox, (int)(line[0][1]*sizeratio)+oy);
		g.drawString("B", (int)(line[1][0]*sizeratio)+ox, (int)(line[1][1]*sizeratio)+oy);
		g.drawLine((int)(line[0][0]*sizeratio)+ox, (int)(line[0][1]*sizeratio)+oy, (int)(line[1][0]*sizeratio)+ox, (int)(line[1][1]*sizeratio)+oy);
	}
	public void drawPlus(Graphics g, int x, int y, int r)
	{
		g.drawLine(x-r, y, x+r, y);
		g.drawLine(x, y-r, x, y+r);
	}

	public void putFFTSpectrum()
	{
		for (int i = 0; i < t.nlayers; i++)
		{
			spectrum[i] = fftmag[i][currenti][currentj];
		}
	}
	public void moveLine(int dx, int dy)
	{
		line[nearestEndIndex][0] += dx/(double)sizeratio;
		line[nearestEndIndex][1] += dy/(double)sizeratio;
	}
	public void moveStripEnd(int dx, int dy)
	{
	
		
	}
	
	public void resetDrawField()
	{
		if (!showFFT) FieldOps.putNM2Array(new double[][][] {t.data[para], u.data[para]}, drawFieldC);
//		else drawField = fftmag[para];
		imageOverField = 1;
		while (drawFieldC.length > defaultSize && drawFieldC.length/2 >= defaultSize)
		{
			this.drawFieldC = FieldOps.reduce(2, this.drawFieldC);
			sx = this.drawFieldC.length;
			zoomLevel++;
			zoomFactor *= 2;
			imageOverField /= 2;
		}
		while (sizeratio*sx < defaultSize)
		{
			sizeratio*=2;
			imageOverField *= 2;
		}
		
	}
	
	public void resetGraphics(boolean redoColorScale, double downNum, double upNum)
	{
		resetDrawField();
		//		setFieldInfo();
		if (redoColorScale)
			cscale = new ColorScales.MYC2d(min[para] + delta[para]*upNum, min[para] + delta[para]*downNum, 2*Math.PI);
	    this.formImage();
	    refresh = true;
	    repaint();
	}
	
	
	public void keyPressed(KeyEvent arg0) {
	}
	public void keyReleased(KeyEvent arg0) {
	}
	public void keyTyped(KeyEvent arg0) {
		System.out.println(arg0.getKeyChar());
		if (arg0.getKeyChar() == 'f')
		{
			//take Fourier transform.
			if (fftmag == null)
			{
				refreshFFT = false;
				fftmag = new double [t.nlayers][t.nx][t.ny];
				for (int i = 0; i < t.nlayers; i++)
				{
					setTitle("Doing FFT: " + i);
					FFTOps.obtainFFTmagCent(t.data[i], fftmag[i]);
					if (ftlog)
						FieldOps.log(fftmag[i]);
				}
				recalculateFFTBounds();
				setTitle("Done");
			}
			else if (refreshFFT)
			{
				refreshFFT = false;
				for (int i = 0; i < t.nlayers; i++)
				{
					setTitle("Doing FFT: " + i);
					FFTOps.obtainFFTmagCent(t.data[i], fftmag[i]);
					if (ftlog)
						FieldOps.log(fftmag[i]);
				}
				recalculateFFTBounds();
			}
			showFFT = !showFFT;
			if (showFFT) {
				min = fftmin;
				max = fftmax;
				delta = fftdelta;
				spec.setRange(ArrayOps.min(fftmag), ArrayOps.max(fftmag));
			}
			else{
				min = realmin;
				max = realmax;
				delta = realdelta;
				spec.setRange(ArrayOps.min(t.data), ArrayOps.max(t.data));
			}
			
			if (lc != null)
			{
				if (showFFT)
					lc.setDataset(fftmag);
				else
					lc.setDataset(t.data);
			}
			else if (sct != null)
			{
				if (showFFT)
					sct.setData(fftmag);
				else
					sct.setData(t.data);
			}
			resetGraphics(true, 0, 1);
		}
		if (arg0.getKeyChar() == 'g')
		{
			ftlog = !ftlog;
			refreshFFT = true;
		}
		if (arg0.getKeyChar() == 'd')
			doAnything();
		if (arg0.getKeyChar() == ' '){
			if (this.showFFT)
				csh.reBoundColorScale(ArrayOps.min(fftmag), ArrayOps.max(fftmag));
			else
				csh.reBoundColorScale(ArrayOps.min(t.data), ArrayOps.max(t.data));
		}
		if (arg0.getKeyChar() == 'p')
		{
			File f = FileOps.selectSave(fc);
			if (f != null)
				SRAW.writeImage(f.toString(), image);
		}
		if (arg0.getKeyChar() == 'G')
		{
				int typeOfGIF = Integer.parseInt(JOptionPane.showInputDialog("Select the output type:\r\n" +
						" 0 - GIF with topo and FFT inset\r\n" +
						" 1 - GIF with FFT only (and voltage)\r\n" + 
						" 2 - Large BMP with layers as tiles\r\n" + 
						" 3 - Movie from BMP to AVI\r\n"));
				File f;
				switch(typeOfGIF){
					case 0:
						f = FileOps.selectSave(fc);
						Layer topo = Layer.openFree(fc);
						RHKFileOps.doFitting(topo, RHKFileOps.getUserFittingChoice());
						if (showFFT)
							ExportUtil.exportGIFForPPTSlide(t, topo, csh.getScale(), null, false, true, true, 500, f.toString(), csh.getCurrentCS());
						else
							ExportUtil.exportGIFForPPTSlide(t, topo, null, csh.getScale(), false, true, true, 500, f.toString(), csh.getCurrentCS());
						break;
					case 1:
						f = FileOps.selectSave(fc);
						TopomapUtil.writeGIFWithVoltage(t, f.toString(), Integer.parseInt(JOptionPane.showInputDialog("How much to blow it up by?")), csh.getCurrentCS());
						break;
					case 2:
						f = FileOps.selectSave(fc);
						int nPerLine = Integer.parseInt(JOptionPane.showInputDialog("How many tiles per line?"));
//						BufferedImage image = ImageEditing.createTileImage(showFFT ? fftmag : t.data, csh.getScale(), nPerLine, t.v);
						BufferedImage image = ImageEditing.createTileImage(showFFT ? fftmag : t.data, null, nPerLine, t.v, csh.getCurrentCS());
						SRAW.writeImage(f.toString(), image);
						ImageEditing.copyToClipboard(image);
						break;
					case 3:
						String shortname = name.contains(" ") ? name.split(" ")[0] : name;
						double minperc = Double.parseDouble(JOptionPane.showInputDialog("Minimum distribution cutoff?"));
						double maxperc = Double.parseDouble(JOptionPane.showInputDialog("Maximum distribution cutoff?"));
						int resizeFactor = Integer.parseInt(JOptionPane.showInputDialog("Blow up the image by what factor?", "" + 1));
						int imin, imax;
						String input = JOptionPane.showInputDialog("Enter the index-limits of the movie.", "" + 0 + "," + t.nlayers);
						imin = Integer.parseInt(input.split(",")[0].trim());
						imax = Integer.parseInt(input.split(",")[1].trim());

						BufferedImage[] stack = ImageEditing.createImageArray(showFFT ? fftmag : t.data, csh.getCurrentCS(), t.v, 0, minperc, maxperc, resizeFactor, null);
						BufferedImage[] substack = new BufferedImage[imax-imin];
						for (int i = imin; i < imax; i++)
							substack[i-imin] = stack[i];
						if (imps != null)
							for (int i = 0; i < substack.length; i++)
							{
								Graphics g = substack[i].getGraphics();
								g.setColor(csh.getUnusedColor());
								double[] r;
								for (int j = 0; j < imps.size(); j++)
								{
									r = imps.get(j).getPosition();
									drawPlus(g, (int)(r[0]*resizeFactor), (int)(r[1]*resizeFactor), 5);
								}
							}
						
						for (int i = 0; i < substack.length; i++)
							SRAW.writeImage(MovieMaker.avidir + shortname + MovieMaker.fromInt(i), substack[i]);
						System.out.println(MovieMaker.BMPtoAVICommand(shortname, 0, substack.length-1, 5));
						break;
				}	
		}
		if (arg0.getKeyChar() == 'l') //toggle the line cut tool
		{
			if (lc == null)
			{	
				if (line == null) line = new double[][] {{0, 0}, {t.nx, t.ny}};
				lc = new LineCutDrawer(t, 2, line);
				lc.setFC(fc);
			}
			else {
				lc.dispose();
				lc = null;
			}
			refresh = true;
			repaint();
		}
		if (arg0.getKeyChar() == 'L') //put a layer onto the spectra window
		{
			if (!spec.drawingLayer){
				Layer layer = Layer.open(fc);
				ColorScale1D scale = ColorScales.getNew(layer.data, csh.getCurrentCS());
				spec.drawingLayer = true;
				spec.layer = layer;
				spec.cscale = scale;
				spec.layerImage = ImageEditing.getBufferedImage(layer.data, scale);
				spec.plot[0].c = ColorScales.getUnusedColor(scale);
				spec.repaint();
			}
			else
			{
				spec.drawingLayer = false;
			}
		}
		if (arg0.getKeyChar() == 'c')
		{
			if (sct == null)
			{
				if (line == null) line = new double[][] {{t.nx/2, t.ny/2}, {t.nx, t.ny}};
				sct = new StripCutTool(t, 2, 0, t.nx/2, 1);
				sct.putLine(line);
				refresh = true;
				repaint();

			}
		}
		if (arg0.getKeyChar() == 'C') //copy the current image to clipboard
		{
//			ImageEditing.copyToClipboard(image);
			BufferedImage export = new BufferedImage(WIDTH, HEIGHT, BufferedImage.TYPE_INT_RGB);
			Graphics g = export.getGraphics();
			refresh = true;
			paint(g);
			ImageEditing.copyToClipboard(ImageEditing.getSubset(export, ox, sx*sizeratio, oy, sy*sizeratio));
		}
		if (arg0.getKeyChar() == 'V')
		{
			ImageEditing.copyToClipboard(csh.getScale().getScaleImage(8));
		}
		if (arg0.getKeyChar() == 'q')
		{
			csh.incrementScaleIndex();
			formImage();
			refresh = true;
			repaint();
		}
		if (arg0.getKeyChar() == 'e')
		{
			csh.decrementScaleIndex();
			formImage();
			refresh = true;
			repaint();
		}
		if (arg0.getKeyChar() == 'r')
		{
			changeScaleWithLayer = !changeScaleWithLayer;
		}
		if (arg0.getKeyChar() == 'x'){
			int pixcoord = JOptionPane.showConfirmDialog(null, " YES: Use original length coordinates in layer \r\n NO: Use pixel length coordinates in layer", "Length Units", JOptionPane.YES_NO_OPTION);
			if (pixcoord == JOptionPane.YES_OPTION)
				Layer.writeBIN(t.getLayer(para));
			else Layer.writeBIN(t.getLayerPixels(para));
		}
		if (arg0.getKeyChar() == 's')
		{
			for (int i = 0; i < t.nlayers; i++)
			{
				System.out.println(t.v[i] + "\t" + spectrum[i]);
			}
		}
		if (arg0.getKeyChar() == 'S')
		{
			int choice = Integer.parseInt(JOptionPane.showInputDialog("Select what to save:\r\n"+
					"0 - Real part\r\n"+
					"1 - Imaginary part\r\n"+
					"2 - Magnitude\r\n"+
					"3 - phase\r\n" +
					"4 - Divergence\r\n"+
					"5 - Curl\r\n"));
			if (choice == 0)
				Topomap.writeBIN(t, fc);
			else if (choice == 1)
				Topomap.writeBIN(u, fc);
			else if (choice == 2)
			{
				double[][][] mag = new double [t.nlayers][t.nx][t.ny];
				for (int i = 0; i < t.nlayers; i++)
					mag[i] = FieldOps.magnitude(t.data[i], u.data[i]);
				Topomap.writeBIN(Topomap.newTopomap(t, mag), fc);
			}
			else if (choice == 3)
			{
				double[][][] phase = new double [t.nlayers][t.nx][t.ny];
				for (int i = 0; i < t.nlayers; i++)
					phase[i] = FieldOps.phase(t.data[i], u.data[i]);
				Topomap.writeBIN(Topomap.newTopomap(t, phase), fc);
			}
			else if (choice == 4)
			{
				double[][][] phase = new double [t.nlayers][t.nx][t.ny];
				for (int i = 0; i < t.nlayers; i++)
					phase[i] = FieldOps.divergence_2NM(new double[][][] {t.data[i], u.data[i]});
				Topomap.writeBIN(Topomap.newTopomap(t, phase), fc);
			}
			else if (choice == 5)
			{
				double[][][] phase = new double [t.nlayers][t.nx][t.ny];
				for (int i = 0; i < t.nlayers; i++)
					phase[i] = FieldOps.curl(FieldOps.getNM2Array(new double[][][] {t.data[i], u.data[i]}));
				Topomap.writeBIN(Topomap.newTopomap(t, phase), fc);
			}
		}		
		if (arg0.getKeyChar() == 'M')
		{
			File f = FileOps.selectSave(fc);
			Layer.writeBIN(TopomapUtil.getSpectralDistributionFancy(128, 128, t), f.toString());
		}
		if (arg0.getKeyChar() == 'P') //print entire stack of images with current color scale
		{
			File f = FileOps.selectSave(fc);
			String name = f.toString();
//			BufferedImage im = new BufferedImage(image.getWidth(), image.getHeight(), BufferedImage.TYPE_INT_RGB);
			BufferedImage im = new BufferedImage(t.nx, t.ny, BufferedImage.TYPE_INT_RGB);
			for (int i = 0; i < t.nlayers; i++)
			{
//				SRAW.writeImage(name + MovieMaker.fromInt(i), t.data[i], csh.getScale());  //To print the bare image:
				if (!showFFT)
					SRAW.writeImage(im, t.data[i], csh.getScale());
				else
					SRAW.writeImage(im, fftmag[i], csh.getScale());
//				int wx = im.getWidth() - 150;
//				int wy = im.getWidth() -20;
//				ImageEditing.writeLine(im, "V = " + t.v[i], wx, wy, csh.getUnusedColor(), Font.SERIF, Font.PLAIN, 32);
				SRAW.writeImage(name + MovieMaker.fromInt(i), im);
			}
		}
		if (arg0.getKeyChar() == 'I' && lc == null)
		{
			if (imps == null) {imps = ImpurityListEditor.getImpurities(fc, 1); refresh = false; repaint();}
			else imps = null;
		}
		if (arg0.getKeyChar() == '+')
		{
			if (para < t.nlayers-1){ para++;
				twoV[0] = t.v[para];
				twoV[1] = t.v[para];
				spec.repaint();
				
				resetGraphics(changeScaleWithLayer, ((double)s.min.getValue()/1000), ((double)s.max.getValue()/1000));
				s.frame.setTitle("" + t.v[para]);
			}
		}
		if (arg0.getKeyChar() == '-')
		{
			if (para > 0){ 
				para--;
				twoV[0] = t.v[para];
				twoV[1] = t.v[para];
				spec.repaint();
				
				resetGraphics(changeScaleWithLayer, ((double)s.min.getValue()/1000), ((double)s.max.getValue()/1000));
				s.frame.setTitle("" + t.v[para]);
			}
		}
		if (arg0.getKeyChar() == 'T' && fitter == null)
		{
			fitter = new FFTCutFitter(t.data, new AtomicCoordinatesSet(FileOps.openText(fc)));
		}
		if (lc != null) lc.processKeyStroke(arg0.getKeyChar());
		else if (sct != null) sct.processKeyStroke(arg0.getKeyChar());		
		
		if (fitter != null){
			if (arg0.getKeyChar() == 'j') //construct the e, q point at the present mouse position and add it to the list
			{
				//first, determine whether it is on- or off- lattice.
				double[] r =new double[] {currentid-t.nx/2, currentjd-t.nx/2};
				double angon = fitter.onLatt.getAngleBetween(r, 0);
				double angoff = fitter.offLatt.getAngleBetween(r, 0);
				double dist = showFFT ? Complex.mag(t.getFourierMetricCoords(currentid, currentjd)) : Distance.distance(t.getMetricCoords(t.nx/2, t.ny/2), t.getMetricCoords(currentid, currentjd));
				boolean onit = angon < angoff;
				System.out.println("" + onit + "\t" + angon + "\t" + angoff);
				double distproj = onit ? dist*Math.cos(angon) : dist*Math.cos(angoff);
				fitter.user.addEQPair(new EQPair(t.v[para], distproj, onit));
			}
			if (arg0.getKeyChar() == 'J')
			{
				fitter.user.removeLastEQPair();
			}
			fitter.processKeyStroke(arg0.getKeyChar());
		}
		
	}
	public void mouseWheelMoved(MouseWheelEvent arg0) {
	}
	public void mouseDragged(MouseEvent arg0) {
		int x = arg0.getX();
		int y = arg0.getY();
		
		double r, phi;
		double[] pos = new double [] {((x-ox)/(double)sizeratio) - t.nx/2, ((y-oy)/(double)sizeratio) - t.ny/2};
		r = Complex.mag(pos);
		phi = Complex.phase(pos);
		
		if (lc != null)
		{
			moveLine((x - currentx), (y - currenty));
			lc.putPoints(line);
			lc.refresh();
			refresh = true;
			repaint();
		}
		else if (sct != null)
		{
			
			double dr, dphi;
			dr = r - currentr;
			dphi = phi-currentphi;
			System.out.println(dr + "\t" + dphi + "\t" + nearestEndIndex);
			if (nearestEndIndex == 0 && sct.cut[0].r0 + dr > 0)
				sct.setR0(sct.cut[0].r0 + dr);
			else if (sct.cut[0].r1 + dr > 0)
				sct.setR1(sct.cut[0].r1 + dr);
			sct.setPhi(sct.cut[0].phi + dphi);
			sct.putLine(line);
			sct.refresh();
			refresh = true;
			repaint();
		}
//		System.out.println(x);
		currentx = x; currenty = y;
		currentr = r;
		currentphi = phi;
		}
	public void mouseMoved(MouseEvent arg0) {
		// TODO Auto-generated method stub
		currentx = arg0.getX();
		currenty = arg0.getY();
		calcx = currentx - ox;
		calcy = currenty - oy;
		currentid = calcx/(double)sizeratio;
		currentjd = calcy/(double)sizeratio;
		currenti = (int)currentid;
		currentj = (int)currentjd;
		currentr = Distance.distance((currenti-t.nx/2), (currentj-t.ny/2));
		currentphi = FieldOps.atan(currentid-t.nx/2, currentjd-t.ny/2);
		if (currenti < 0) currenti = 0;
		if (currenti >= t.nx) currenti = t.nx-1;
		if (currentj < 0) currentj = 0;
		if (currentj >= t.ny) currentj = t.ny-1;
		refresh = false;
		repaint();
		
		spec.setTitle("Spectrum [" + currenti + ", " + currentj + "]");
		if (withinBounds(currenti, currentj, t.data[0]) && !showFFT)
			t.putSpectrum(currenti, currentj, spectrum);
		else if (showFFT){
			putFFTSpectrum();
			//			spec.setRange(ArrayOps.min(spectrum), ArrayOps.max(spectrum));
		}
		spec.repaint();
		
		if (lc != null || sct != null)
		{
			d1 = Distance.distance(currentid - line[0][0], currentjd - line[0][1]);
			d2 = Distance.distance(currentid - line[1][0], currentjd - line[1][1]);
			if (d1 < d2) nearestEndIndex = 0;
			else nearestEndIndex = 1;
//			System.out.println(nearestEndIndex);
		}		
//		refresh = false;
//		repaint();
	}
	public int screenX(double currentid)
	{
		return (int)(currentid*sizeratio + ox);
	}
	public int screenY(double currentjd)
	{
		return (int)(currentjd*sizeratio + oy);
	}
	public void mouseClicked(MouseEvent arg0) {
		
		if (fitter != null)
		{
			
		}
		
	}
	public void mouseEntered(MouseEvent arg0) {
	}
	public void mouseExited(MouseEvent arg0) {
	}
	public void mousePressed(MouseEvent arg0) {
	}
	public void mouseReleased(MouseEvent arg0) {
	}
	
	public static class SliderPanel extends JPanel implements ChangeListener
	{
		TopomapViewer_complex2 parent;
		JFrame frame;
		public JSlider s;
		public JSlider min, max;
		int parentmember;
		int oldvalue = 0;
		int oldminv = 0, oldmaxv = 999;
		static int npts = 1001;
		
		public SliderPanel(TopomapViewer_complex2 parent, JFrame frame)
		{
			npts = parent.snpts;
			this.frame = frame;
			this.parent = parent;
			this.setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
			this.setBorder(new LineBorder(Color.GRAY));
			s = new JSlider(0, npts, 0);
			s.setValue(0);
			
			s.setSnapToTicks(true);
			s.setMinorTickSpacing(1);
			s.addChangeListener(this);
			
			min = new JSlider(0, 1000, oldminv);
			max = new JSlider(0, 1000, oldmaxv);
			min.addChangeListener(this);
			max.addChangeListener(this);
			add(s);
			add(min);
			add(max);
			frame.setSize(3*npts+40, 80);
			frame.add(this);
		}
		public void show(){
			frame.show();
		}
		@Override
		public void stateChanged(ChangeEvent arg0) {
//			if (s.getValue() == oldvalue && min.getValue() == oldminv && max.getValue() == oldmaxv) return;
			if (min.getValue() != oldminv || max.getValue() != oldmaxv)
			{
				parent.resetCSHFromSlider(((double)min.getValue()/1000), ((double)max.getValue()/1000));
				parent.formImage();
				parent.refresh = true;
				parent.repaint();
				oldminv = min.getValue();
				oldmaxv = max.getValue();
				return;
			}
//			else if (s.getValue() == oldvalue) return;
//			parent.resetColorScale(((double)min.getValue()/1000), ((double)max.getValue()/1000));
			else
			{
				parent.para = s.getValue();
				parent.twoV[0] = parent.t.v[parent.para];
				parent.twoV[1] = parent.t.v[parent.para];
				parent.spec.repaint();
				
				parent.resetGraphics(parent.changeScaleWithLayer, ((double)min.getValue()/1000), ((double)max.getValue()/1000));
				frame.setTitle("" + parent.t.v[parent.para]);
				System.out.println(s.getValue());
				oldvalue = s.getValue();
			}
		}
	}
	
	public static boolean withinBounds(int i, int j, double[][] array)
	{
		return i >= 0 && j >= 0 && i < array.length && j < array[i].length;
	}
	public void recalculateBounds()
	{
		for (int i = 0; i < t.nlayers; i++)
		{
			this.realmin[i] = FieldOps.magMin(FieldOps.getNM2Array(new double[][][] {t.data[i], u.data[i]}));
			this.realmax[i] = FieldOps.magMax(FieldOps.getNM2Array(new double[][][] {t.data[i], u.data[i]}));
			this.realdelta[i] = this.max[i]-this.min[i];
		}
	}
	public void recalculateFFTBounds()
	{
		for (int i = 0; i < t.nlayers; i++)
		{
			this.fftmin[i] = FieldOps.min(fftmag[i]);
			this.fftmax[i] = FieldOps.max(fftmag[i]);
			this.fftdelta[i] = this.fftmax[i]-this.fftmin[i];
		}
	}
	public void doAnything()
	{
		int o = Integer.parseInt(JOptionPane.showInputDialog(null, "What would you like to do? (Note 3 will be applied to both real and imaginary parts)\r\n"+
				"1 - Do a work-function fitting\r\n" +
				"2 - Symmetrize\r\n" +
				"3 - Do some data processing technique to each layer\r\n"+
				"4 - Rescale the bias voltages by a constant factor (for Z-map)\r\n"+ 
				"5 - print a spectrum of the correlation with a layer (opens a layer)\r\n"+
				"6 - Split the map into bins based on an existing layer\r\n" +
				"7 - Subtract a polynomial fit from each spectrum\r\n" + 
				"8 - Get locally-averaged spectra at impurities\r\n" +
				"9 - Get a simple histogram of the spectra\r\n" +
				"10 - Smoothly suppress features to the else average\r\n" +
				"11 - Replace the map with its full 3D FFT\r\n" + 
				"12 - Sum over one of the dimensions\r\n" +
				"13 - Smooth each spectrum with a gaussian\r\n" + 
				"14 - Subtract a polynomial from each spectrum\r\n" +
				"15 - Graph the average FFT of all the spectra\r\n" + 
				"16 - Try to filter out the noise\r\n" + 
				"17 - Replace each spectrum by its FFT\r\n" +
				"18 - Replace each spectrum by its derivative\r\n"
				));
		
		
		switch(o)
		{
		case 1:
			TopomapUtil.fitEachSpectrumExponential(t, "");
			break;
		case 2:
			int o2 = JOptionPane.showConfirmDialog(null, "Yes for 6-fold\r\nNo for 4-fold");
			t = TopomapUtil.makeSymmetrizedFFTs(t, new AtomicCoordinatesSet(FileOps.openText(fc)), o2 == JOptionPane.NO_OPTION);
			resetGraphics(true, 0, 1);
			break;
		case 3:
			int o3 = RHKFileOps.getUserFittingChoice();
			RHKFileOps.doFittingTopo(t, o3);
			RHKFileOps.doFittingTopo(u, o3);
			recalculateBounds();
			break;
		case 4:
			double fz = Double.parseDouble(JOptionPane.showInputDialog("Enter the factor to multiply Z by."));
			for (int i = 0; i < t.nlayers; i++)
				t.v[i] *= fz;
		break;
		case 5:
			Layer l = Layer.open(fc);
			System.out.println();
			for (int i = 0; i < t.nlayers; i++)
				System.out.println("" + t.v[i] + "\t" + FieldOps.correlation(l.data, t.data[i]));
		break;
		case 6:
			TopomapUtil.splitTopomapByBinsOfLayer(t, dir, name, Integer.parseInt(JOptionPane.showInputDialog("Enter the number of bins")), fc);
		break;
		case 7:
			TopomapUtil.subtractPolynomialFitEachSpectrum(t, Integer.parseInt(JOptionPane.showInputDialog("Enter the degree of polynomial")));
			break;
		case 8:
			TopomapUtil.writeAverageSpectraAroundImps(t, Double.parseDouble(JOptionPane.showInputDialog("Enter the gaussian radius in pixels")), fc);
			break;
		case 9:
			Layer hist = TopomapUtil.getSpectralDistributionBasic(1, Integer.parseInt(JOptionPane.showInputDialog("Enter the number of bins.")), t);
			int size = t.nlayers;
			while (size*((double)hist.ny/t.nlayers) < 768)
				size*=2;
			
			lv = new LayerViewer(hist, dir, size);
			lv.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			break;
		case 10:
			Layer mask = Layer.open(fc);
			for (int i = 0; i < t.nlayers; i++)
				FieldOps.changeToWithoutMask(t.data[i], mask.data);
			break;
		case 11:
			t.data = FFT3D_Wrapper.getFFTMagCentXY(t.data, JOptionPane.showConfirmDialog(null, "Take the log?") == JOptionPane.YES_OPTION);
			break;
		case 12:
			lv = new LayerViewer(TopomapUtil.sumOneDimension(t, Integer.parseInt(JOptionPane.showInputDialog("Enter the code:\r\n\t 0 - Energy \r\n\t 1 - X \r\n\t 2 - Y"))), dir, 1024);
			lv.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			break;
		case 13:
			double[] tempSpec;
			double L = Double.parseDouble(JOptionPane.showInputDialog("Enter the smoothing length in pixels."));
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
				{
					tempSpec = ArrayOps.gaussSmooth(t.getSpectrum(i, j), L);
					for (int k = 0; k < t.nlayers; k++)
						t.data[k][i][j] = tempSpec[k];
				}
			break;
		case 14:
			int degree = Integer.parseInt(JOptionPane.showInputDialog("Enter the degree of polynomial."));
			if (JOptionPane.showConfirmDialog(null, "Save the fits before subtracting?") == JOptionPane.YES_OPTION)
			{
				double[][][] fits = new double[t.nlayers][t.nx][t.ny];
				tempSpec = new double[t.nlayers];
				double[][] tempans;
				for (int i = 0; i < t.nx; i++)
					for (int j = 0; j < t.ny; j++)
					{
						t.putSpectrum(i, j, tempSpec);
						tempans = ArrayOps.subtractPolynomialFit(t.v, tempSpec, degree);
						for (int k = 0; k < t.nlayers; k++)
						{
							t.data[k][i][j] = tempans[0][k];
							fits[k][i][j] = tempans[1][k];
						}
					}
				
				Topomap.writeBIN(Topomap.newTopomap(t, fits), fc);
			}
			else
			{
				TopomapUtil.subtractPolynomialFitEachSpectrum(t, degree);
			}
			break;
		case 15:
			if (JOptionPane.showConfirmDialog(null, "Take the log?") == JOptionPane.NO_OPTION)
				GraphDrawerCart.plotGraph(ArrayOps.generateArrayNotInclUpper(0, t.nlayers, t.nlayers), TopomapUtil.getAverageFFTOfSpectra(t));
			else{
				double[] temp = TopomapUtil.getAverageFFTOfSpectra(t);
				FieldOps.log(temp);
				GraphDrawerCart.plotGraph(ArrayOps.generateArrayNotInclUpper(0, t.nlayers, t.nlayers), temp);
				
			}
			break;
		case 16:
			degree = Integer.parseInt(JOptionPane.showInputDialog("Enter the degree of polynomial."));
			double[][][] fits = new double[t.nlayers][t.nx][t.ny];
			tempSpec = new double[t.nlayers];
			double[][] tempans = null;
			for (int i = 0; i < t.nx; i++){System.out.print("" + i + " ");
				for (int j = 0; j < t.ny; j++)
				{
					t.putSpectrum(i, j, tempSpec);
					tempans = ArrayOps.subtractPolynomialFit(t.v, tempSpec, degree);
					for (int k = 0; k < t.nlayers; k++)
					{
						t.data[k][i][j] = tempans[0][k];
						fits[k][i][j] = tempans[1][k];
					}
				}
			}
			if (JOptionPane.showConfirmDialog(null, "Save the fits?") == JOptionPane.YES_OPTION)
				Topomap.writeBIN(Topomap.newTopomap(t, fits), fc);
				
			String input = JOptionPane.showInputDialog("Enter upper left corner of the box, spearated by commas.");
			int imin = Integer.parseInt(input.split(",")[0].trim());
			int jmin = Integer.parseInt(input.split(",")[1].trim());
//			input = JOptionPane.showInputDialog("Enter width and height of the box, separated by commas.");
//			int width = Integer.parseInt(input.split(",")[0].trim());
//			int height = Integer.parseInt(input.split(",")[1].trim());
//			ArrayList<TopomapUtil.FilterBox> box = new ArrayList<TopomapUtil.FilterBox>();
//			box.add(new TopomapUtil.FilterBox(imin, jmin, width, height));
//			Topomap ans = TopomapUtil.fourierFilterTheMap(t, box);
//			for (int i = 0; i < t.nx; i++)
//				for (int j = 0; j < t.ny; j++)
//					for (int k = 0; k < t.nlayers; k++)
//						t.data[k][i][j] = ans.data[k][i][j] + fits[k][i][j];
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++){
					tempans[0] = FFTOps.getFourierFiltered(t.getSpectrum(i, j), new int[] {imin}, new int[] {jmin});
					for (int k = 0; k < t.nlayers; k++)
						t.data[k][i][j] = tempans[0][k] + fits[k][i][j];
				}
				break;
		case 17:
			TopomapUtil.replaceWithFFTOfSpectra(t);
			break;
		case 18:
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
				{
					tempSpec = ArrayOps.getDerivative(t.getSpectrum(i, j));
					for (int k = 0; k < t.nlayers; k++)
						t.data[k][i][j] = tempSpec[k];
				}
			break;
		}
		refreshFFT = true;
		recalculateBounds();
	}
	public static void main(String[] args)
	{
		Topomap.setStdDir();
		fc = new JFileChooser(Topomap.stddir);
		boolean doLayers = JOptionPane.showConfirmDialog(null, "Open two layers?") == JOptionPane.YES_OPTION;
		Topomap t, u;
		if (doLayers){
			t = Topomap.newTopomap(new Layer[] {Layer.openFree(fc)});
			u = Topomap.newTopomap(new Layer[] {Layer.openFree(fc)});
		}
		else{
			t = Topomap.open(fc);
			u = Topomap.open(fc);
		}
		int sizechoice = t.nx;
		while (sizechoice < 600) sizechoice *= 2;
		new TopomapViewer_complex2(t, u, fc.getCurrentDirectory().toString() + "\\", sizechoice);
	}

}
