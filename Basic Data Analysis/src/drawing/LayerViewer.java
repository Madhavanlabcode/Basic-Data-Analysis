package drawing;

import image.GifSequenceWriter;
import image.ImageEditing;
import impurity.PointImp;


//Added a library to get the screen dimensions
import java.awt.Toolkit;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Point;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionListener;
import java.awt.event.MouseWheelEvent;
import java.awt.event.MouseWheelListener;
import java.awt.image.BufferedImage;
import java.util.ArrayList;

import javax.swing.BoxLayout;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.border.LineBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import main.SRAW;
import util.ArrayOps;
import util.Complex;
import util.FieldOps;
import util.LayerUtil;
import util.NumFormat;
import util.Printer;
import util.SpectraUtil;
import util.TopomapUtil;
import util.color.ColorScale;
import util.color.ColorScale1D;
import util.color.ColorScale2d;
import util.color.ColorScales;
import util.fileops.ColumnIO;
import util.fileops.FileOps;
import util.fileops.Layer;
import util.fileops.PointSpectra;
import util.fileops.RHKFileOps;
import util.fileops.Topomap;
import util.fourier.FFTOps;
import util.fourier.ImpurityListEditor.Impurity;
import util.fourier.ImpurityListEditor;
import util.geom.AtomicCoordinatesSet;
import util.geom.Distance;
import util.regression.TwoDGaussianFreeFitter;
import util.robot.Robo;


public class LayerViewer extends JFrame implements MouseListener, MouseMotionListener, MouseWheelListener, KeyListener{
	private static final long serialVersionUID = 101516031674408198L;
	Layer t, backup = null;
	int N;
	double fmax, fmin, fdelta;
	double[][] fftmag = null;
	double[][][] fftz = null;
	boolean showFFT = false;
	double[][] drawField;
	double[][][] drawFieldC;
	boolean real = true;
	
	static JFileChooser fc;
	Robo rob=new Robo();
	
	BufferedImage image;
	boolean refresh = true;
	
	ColorScale scale;
	ColorScale2d cscale;
	double scalemin, scalemax;
	double csdown, csup;
		
	int sx = 512, sy = 512;
	int ox = 20, oy = 60;
	
	int[] writepoint = {ox + sx + 50, oy + 10};
	int linesize = 15;
	
	int zoomLevel = 0, zoomFactor = 1;
	int sizeratio = 1;
	
	int WIDTH = 1600, HEIGHT = 1080;
	
	//para is the current index
	int para = 0;
	
	int currentx, currenty;
	int currenti, currentj;
	double currentid, currentjd;
	int calcx = 0, calcy = 0;
	String dir;
	double[] pbounds;
	SliderPanel s;

	int currentCScale = 0;
	boolean expandedScaleRange = false;
	
	boolean ftlog = true;
	boolean switchLog = false;
	
	//for PointSpectra.
	PointSpectra spec = null;
	PointSpectra[] splitPos = null;
	boolean showingSpec = false;
	SpectraDrawer specView = null;
	int selectedSpectra = 0;
	
	LayerCutDrawer sectionDrawer = null;
	double[][] line = null; //the two points
	public int nearestEndIndex = 0;
	double d1, d2;
	double sOfCursor;

	int defaultSize = 1024;
	double imageOverField = 1;
	
	//for the booleans
	boolean[][] truePix = null;
	boolean[][] truePixImage = null;
	boolean drawingTP = false;
	
//	double lineTime = 29.999; //to find in k-space the temporal frequency
	double mouseFreq = 0;
	private ArrayList<Impurity> imps, imps2;
	
	Point mouseDownCompCoords=null;
	
	JMenuBar menuBar;
	
	TwoDGaussianFreeFitter fitter = null;
	public LayerViewer(Layer inLayer, String inDir, int size)
	{
		
		t = inLayer;
//		this.dir = dir;
		if (fc != null)
			dir = fc.getCurrentDirectory().toString() + "\\";
		else
		{
			dir = inDir;
			Topomap.setStdDir();
			fc = new JFileChooser(Topomap.stddir);
		}	
			defaultSize = size;
		
		if (defaultSize == 512)
		{
			WIDTH = 1000;
			HEIGHT = 562;
		}
		
		
		N = t.nx;
		sx = N;
		drawField = t.data;
		
		while (drawField.length > defaultSize && drawField.length/2 >= defaultSize)
		{
			this.drawField = FieldOps.reduce(2, this.drawField);
			sx = this.drawField.length;
			zoomLevel++;
			zoomFactor *= 2;
			imageOverField /= 2;
		}
		//Change this if you want size ratios other than 2^n
		while (sizeratio*sx < defaultSize)
		{
			sizeratio++;
			imageOverField++;
		}
		
		writepoint[0] = ox + sizeratio*sx + 50;
		WIDTH = sizeratio*sx + 450;
		HEIGHT = sizeratio*sx + 50;

		sy = this.drawField[0].length;
		setFieldInfo();
		image = new BufferedImage(sx*sizeratio, sy*sizeratio, BufferedImage.TYPE_INT_RGB);
		csup = 1;
		csdown = 0;
		resetColorScale(ArrayOps.min(t.data), ArrayOps.max(t.data), true);
		formImage();
//		formFTImage();
		//gradcalc.activate();
		
		//Setting up menu bar
		menuBar = new JMenuBar();
		setJMenuBar(menuBar);
		
	    // File Menu, F - Mnemonic
	    JMenu fileMenu = new JMenu("File");
	    fileMenu.setMnemonic(KeyEvent.VK_F);
	    menuBar.add(fileMenu);
	    
	    //File->Save as bin
	    JMenuItem saveAsBinMI = new JMenuItem("Save as a .bin");
	    fileMenu.add(saveAsBinMI);
	    class saveAsBinAL implements ActionListener{
	    	Layer t;
	    	saveAsBinAL(Layer tPrime){
	    		t=tPrime;
	    	}
	    	public void actionPerformed(ActionEvent e){
	    		Layer.writeBIN(t, fc);
	    	}
	    }
	    saveAsBinMI.addActionListener(new saveAsBinAL(t));
	   
	    
	    //File->Save layer as csv
	    JMenuItem saveAsCsvMI = new JMenuItem("Save as comma-separated text file");
	    fileMenu.add(saveAsCsvMI);
	    class saveAsCsvAL implements ActionListener{
	    	Layer t;
	    	saveAsCsvAL(Layer tPrime){
	    		t=tPrime;
	    	}
	    	public void actionPerformed(ActionEvent e){
	    		ColumnIO.writeTableCSV(t.data, FileOps.selectSave(fc).toString());
	    	}
	    }
	    saveAsCsvMI.addActionListener(new saveAsCsvAL(t));
	    
	    //File->Save layer as table
	    JMenuItem saveAsTableMI = new JMenuItem("Save as tab-separated table");
	    fileMenu.add(saveAsTableMI);
	    class saveAsTableAL implements ActionListener{
	    	Layer t;
	    	saveAsTableAL(Layer tPrime){
	    		t=tPrime;
	    	}
	    	public void actionPerformed(ActionEvent e){
	    		ColumnIO.writeTable(t.data, FileOps.selectSave(fc).toString());
	    	}
	    }
	    saveAsTableMI.addActionListener(new saveAsTableAL(t));
	    
	    //File->Save layer as text table
	    JMenuItem saveAsTextTableMI = new JMenuItem("Save as a text table");
	    fileMenu.add(saveAsTextTableMI);
	    class saveAsTextTableAL implements ActionListener{
	    	Layer t;
	    	saveAsTextTableAL(Layer tPrime){
	    		t=tPrime;
	    	}
	    	public void actionPerformed(ActionEvent e){
				boolean transpose = JOptionPane.showConfirmDialog(null, "Transpose the data?") == JOptionPane.YES_OPTION;
				double[][] table = new double [t.ny+1][t.nx+1];
				for (int i = 0; i < t.nx; i++)
					table[0][i+1] = t.x[i];
				for (int j = 0; j < t.ny; j++)
					table[j+1][0] = t.y[j];
				for (int i = 0; i < t.nx; i++)
					for (int j = 0; j < t.ny; j++)
						table[j+1][i+1] = t.data[i][j];
				if (!transpose)
					table = FieldOps.transpose(table);
				
				FileOps.writeTableASCII(fc, table);
	    	}
	    }
	    saveAsTextTableMI.addActionListener(new saveAsTextTableAL(t));

	    //File->merge layers to topomap
	    JMenuItem mergeToMapMI = new JMenuItem("Merge layers into a topomap");
	    fileMenu.add(mergeToMapMI);
	    class mergeToMapAL implements ActionListener{
	    	mergeToMapAL(){}
	    	public void actionPerformed(ActionEvent e) {
	    		String input = JOptionPane.showInputDialog("Enter the number of layers");
				int numLay = Integer.parseInt(input);
				Layer[] lays = new Layer[numLay];
				for(int x=0;x<numLay;x++){
		    		lays[x]=Layer.openFree(fc);
				}
				new TopomapViewer(Topomap.newTopomap(lays),dir,size);
	    	}
	    }
	    mergeToMapMI.addActionListener(new mergeToMapAL());
	    
	    //File->Copy current image to clipboard
        JMenuItem toClipMI = new JMenuItem("Copy current image to clipboard");
        fileMenu.add(toClipMI);
	    class toClipAL implements ActionListener{
	    	public void actionPerformed(ActionEvent e) {
				BufferedImage export = new BufferedImage(WIDTH, HEIGHT, BufferedImage.TYPE_INT_RGB);
				Graphics g = export.getGraphics();
				refresh = true;
				paint(g);
				try
				{
					ImageEditing.copyToClipboard(ImageEditing.getSubset(export, ox, sx*sizeratio, oy, sy*sizeratio));
				}
				catch(Exception arg0){
					ImageEditing.copyToClipboard(image);
				}
	    	}
	    }
	    toClipMI.addActionListener(new toClipAL());
	    
        //File->Copy image to clipboard with custom settings
        JMenuItem imageToClipMI = new JMenuItem("Paste image to clipboard with custom settings");
        fileMenu.add(imageToClipMI);
	    class imageToClipAL implements ActionListener{
	    	Layer t;
	    	imageToClipAL(Layer tPrime){
	    		t=tPrime;
	    	}
	    	public void actionPerformed(ActionEvent e) {
				String input = JOptionPane.showInputDialog("Enter the color scale range (comma separated).", "" + scale.getMin() + "," + scale.getMax());
				double min = Double.parseDouble(input.split(",")[0].trim());
				double max = Double.parseDouble(input.split(",")[1].trim());
				ColorScale1D cs = ColorScales.getNew(min, max, currentCScale);

				int[] sizeRatio = new int [2];
				sizeRatio[0] = Integer.parseInt(JOptionPane.showInputDialog("By how much would you like to blow up the image? (X)", "" + sizeratio));
				sizeRatio[1] = Integer.parseInt(JOptionPane.showInputDialog("By how much would you like to blow up the image? (Y)", "" + sizeratio));
				BufferedImage fullImage = ImageEditing.getBufferedImage(showFFT ? fftmag : t.data, cs);
				BufferedImage im = ImageEditing.enlargeBasicStretch(fullImage, sizeRatio[0], sizeRatio[1]);
				ImageEditing.copyToClipboard(im);
	    	}
	    }
	    imageToClipMI.addActionListener(new imageToClipAL(t));
	    
	    //File->Blow up image and copy to clipboard with data on top
        JMenuItem imagePlusDataMI = new JMenuItem("Blow up image and copy to clipboard with data on top");
        fileMenu.add(imagePlusDataMI);
	    class imagePlusDataAL implements ActionListener{
	    	Layer t;
	    	imagePlusDataAL(Layer tPrime){
	    		t=tPrime;
	    	}
	    	public void actionPerformed(ActionEvent e) {
				int hugeRatio = 18;
				
				BufferedImage huge = ImageEditing.getBufferedImage(t.data, scale);
				huge = ImageEditing.getEnlarged(huge, hugeRatio);
				Graphics g = huge.getGraphics();
				g.setColor(ColorScales.getUnusedColor(scale));
				g.setFont(new Font("Lucida Console", Font.PLAIN, hugeRatio/3));
				for (int i = 0; i < t.nx; i++)
					for (int j = 0; j < t.ny; j++)
					{
						double offset = (i % 2) == 0 ? 0.7 : 0.4;
						g.drawString(NumFormat.scientific(t.data[i][j], 2), i*hugeRatio, (int)((j+offset)*hugeRatio));
					}
				ImageEditing.copyToClipboard(huge);
	    	}
	    }
	    imagePlusDataMI.addActionListener(new imagePlusDataAL(t));
	    
	    //File->open new layer
        JMenuItem openLayerMI = new JMenuItem("Open a new layer");
        fileMenu.add(openLayerMI);
	    class openLayerAL implements ActionListener{
	    	Layer t;
	    	openLayerAL(Layer tPrime){
	    		t=tPrime;
	    	}
	    	public void actionPerformed(ActionEvent e) {
	    		FieldOps.copy(Layer.openFree(fc).data, t.data);
				resetGraphics(true);
	    	}
	    }
	    openLayerMI.addActionListener(new openLayerAL(t));
	    
        //Display menu
	    JMenu dispMenu = new JMenu("Display");
	    menuBar.add(dispMenu);
	    
	    //Display->refresh graphics
        JMenuItem refreshMI = new JMenuItem("Refresh picture (shortcut spacebar)");
        dispMenu.add(refreshMI);
        refreshMI.addActionListener(new ActionListener(){
      	  	public void actionPerformed(ActionEvent e) {
    			refresh = true;
    			repaint();
      	  	}
        });

	    //Display->forward in color scale
        JMenuItem forwardScaleMI = new JMenuItem("Increment color scale (shortcut 'q')");
        dispMenu.add(forwardScaleMI);
        forwardScaleMI.addActionListener(new ActionListener(){
      	  	public void actionPerformed(ActionEvent e) {
    			currentCScale++;
    			currentCScale += ColorScales.NSCALES;
    			currentCScale %= ColorScales.NSCALES;
    			resetColorScale(csdown, csup, false);
    			formImage();
    			setTitle("Color scale: " + (currentCScale+1) + " out of " + ColorScales.NSCALES);
    			refresh = true;
    			repaint();
      	  	}
        });
	    
	    //Display->backward in color scale
        JMenuItem backwardScaleMI = new JMenuItem("Decrement color scale (shortcut 'e')");
        dispMenu.add(backwardScaleMI);
        backwardScaleMI.addActionListener(new ActionListener(){
      	  	public void actionPerformed(ActionEvent e) {
	      	  	currentCScale--;
				currentCScale += ColorScales.NSCALES;
				currentCScale %= ColorScales.NSCALES;
				resetColorScale(csdown, csup, false);
				formImage();
				setTitle("Color scale: " + (currentCScale+1) + " out of " + ColorScales.NSCALES);
				refresh = true;
				repaint();
      	  	}
        });
	    
	    //Display->set color scale to extremes
	    JMenuItem setColorExtremesMI = new JMenuItem("Set color scale by extremes");
	    dispMenu.add(setColorExtremesMI);
	    class setColorExtremesAL implements ActionListener{
	    	Layer t;
	    	setColorExtremesAL(Layer tPrime){
	    		t=tPrime;
	    	}
	    	public void actionPerformed(ActionEvent e) {
	    		if(showFFT)
					resetColorScale(ArrayOps.min(fftmag), ArrayOps.max(fftmag), true);
				else
					resetColorScale(ArrayOps.min(t.data), ArrayOps.max(t.data), true);
	    		
	    		refresh=true;
	    		repaint();
	    	}
	    }
	    setColorExtremesMI.addActionListener(new setColorExtremesAL(t));
	    
	    //Display->toggle extended color scale
        JMenuItem extendColorMI = new JMenuItem("Toggle extended color scale");
        dispMenu.add(extendColorMI);
        extendColorMI.addActionListener(new ActionListener(){
      	  	public void actionPerformed(ActionEvent e) {
      	  		expandedScaleRange = !expandedScaleRange;
      	  	}
        });
        
        //Display->Copy color scale gradient to clipboard
        JMenuItem copyScaleMI = new JMenuItem("Copy color scale gradient to clipboard");
        dispMenu.add(copyScaleMI);
        copyScaleMI.addActionListener(new ActionListener(){
      	  	public void actionPerformed(ActionEvent e) {
      	  		ImageEditing.copyToClipboard(((ColorScale1D)scale).getScaleImage(8));
      	  	}
        });
	    
	    // Fourier Transform menu
	    JMenu fourierMenu = new JMenu("Fourier Transforms");
	    menuBar.add(fourierMenu);
		
        //Fourier->take transform
        JMenuItem fftMI = new JMenuItem("Take/revert FFT (shortcut 'f')");
        fourierMenu.add(fftMI);
    	class fftAL implements ActionListener{
    		private Layer t;
    		fftAL(Layer tPrime){
    			t=tPrime;
    		}
    	    public void actionPerformed(ActionEvent e) {
    	    	//take Fourier transform.
    			if (fftmag == null || switchLog == true)
    			{
    				fftmag = new double [t.nx][t.ny];
    				FFTOps.obtainFFTmagCent(t.data, fftmag);
    				if (ftlog)
    					FieldOps.log(fftmag);
    				setTitle("Done");
    				switchLog = false;
    			}
    			showFFT = !showFFT;
    			resetGraphics(true);
           }
        }
        fftMI.addActionListener(new fftAL(t));
        
        //Fourier->toggle log scale in transform
        JMenuItem toggleLogMI = new JMenuItem("Toggle log scale in transform");
        fourierMenu.add(toggleLogMI);
        toggleLogMI.addActionListener(new ActionListener(){
        	public void actionPerformed(ActionEvent e){
    			ftlog = !ftlog;
    			switchLog = true;
        	}
        });
        
        //Selection menu
        JMenu selectionMenu = new JMenu("Selection");
        menuBar.add(selectionMenu);
        
        //Selection->crop
        JMenuItem cropMI = new JMenuItem("Crop layer");
        selectionMenu.add(cropMI);
    	class cropAL implements ActionListener{
    		private Layer t;
    		cropAL(Layer tPrime){
    			t=tPrime;
    		}
    	    public void actionPerformed(ActionEvent e) {
    			String imput2 = JOptionPane.showInputDialog("Enter coordinates of the upper-left corner of the box, comma separated.");
    			int xi = Integer.parseInt(imput2.split(",")[0].trim());
    			int yi = Integer.parseInt(imput2.split(",")[1].trim());
    			imput2 = JOptionPane.showInputDialog("Enter the width and height of the box, comma separated.");
    			int xf = xi + Integer.parseInt(imput2.split(",")[0].trim());
    			int yf = yi  + Integer.parseInt(imput2.split(",")[1].trim());
    			LayerViewer.show(LayerUtil.cropLayer(t, xi, xf, yi, yf), sx, false);
           }
        }
    	cropMI.addActionListener(new cropAL(t));
    	
    	//Selection->Write points above the cutoff to a point imps file
        JMenuItem highToImpsMI = new JMenuItem("Save points above cutoff as a points file");
        selectionMenu.add(highToImpsMI);
    	class highToImpsAL implements ActionListener{
    		private Layer t;
    		highToImpsAL(Layer tPrime){
    			t=tPrime;
    		}
    	    public void actionPerformed(ActionEvent e) {
    			PointImp[] imps = PointImp.getAboveCutoffContiguous(t, scale.getMax());
    			PointImp.writeToFile(imps, FileOps.selectSave(fc));
           }
        }
    	highToImpsMI.addActionListener(new highToImpsAL(t));
    	
    	//Selection->toggle viewing a points file
        JMenuItem viewImpsMI = new JMenuItem("Toggle viewing a points file");
        selectionMenu.add(viewImpsMI);
    	class viewImpsAL implements ActionListener{
    	    public void actionPerformed(ActionEvent e) {
    			if (imps == null) {imps = ImpurityListEditor.getImpurities(fc, 1); refresh = false; repaint();}
    			else{
    				imps = null;
    				refresh = true;
    				repaint();
    			}
           }
        }
    	viewImpsMI.addActionListener(new viewImpsAL());
    	
    	//Selection->toggle viewing two points files
        JMenuItem viewTwoImpsMI = new JMenuItem("Toggle viewing two points files");
        selectionMenu.add(viewTwoImpsMI);
    	class viewTwoImpsAL implements ActionListener{
    	    public void actionPerformed(ActionEvent e) {
    			if (imps == null) {
    				imps = ImpurityListEditor.getImpurities(fc, 1);
    				imps2 = ImpurityListEditor.getImpurities(fc, 1);
    				refresh = false; repaint();
    				}
    			else{
    				imps = null;
    				imps2 = null;
    				refresh = true;
    				repaint();
    			}
           }
        }
    	viewTwoImpsMI.addActionListener(new viewTwoImpsAL());
    	
    	//Selection->Use mask bitmap to split layer into two
        JMenuItem maskLayerMI = new JMenuItem("Use mask (bitmap or Layer) to split layer into two");
        selectionMenu.add(maskLayerMI);
        class maskLayerAL implements ActionListener{
        	Layer t;
        	maskLayerAL(Layer tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
        		Layer spf = Layer.openFree(fc);
    			new LayerViewer(Layer.newLayer(t, FieldOps.spatialFilter(t.data, spf.data)), dir, defaultSize);
    			new LayerViewer(Layer.newLayer(t, FieldOps.spatialFilter(t.data, ArrayOps.subtract(ArrayOps.max(spf.data),spf.data))), dir, defaultSize);
        	}
        }
        maskLayerMI.addActionListener(new maskLayerAL(t));
	
        //Analysis menu
        JMenu analysisMenu = new JMenu("Analysis");
        menuBar.add(analysisMenu);
        
        //Analysis->Fit a bump
        JMenuItem bumpFitMI = new JMenuItem("Fit a bump");
        analysisMenu.add(bumpFitMI);
        class bumpFitAL implements ActionListener{
        	Layer t;
        	String dir;
        	bumpFitAL(Layer tPrime, String dirPrime){
        		t=tPrime;
        		dir=dirPrime;
        	}
        	public void actionPerformed(ActionEvent e){
        		if (fitter == null)
        		{
        			fitter = new TwoDGaussianFreeFitter(showFFT ? fftmag : t.data);
        			JOptionPane.showMessageDialog(null, "Mouse over the peak center and press F.");
        		}
        		else
        		{
        			switch(fitter.userStage)
        			{
        			case 0:
        				fitter.initGuesses[0] = currentid;
        				fitter.initGuesses[1] = currentjd;
        				fitter.userStage++;
        				JOptionPane.showMessageDialog(null, "Mouse over the right side of the peak and press F.");
        				break;
        			case 1:
        				fitter.initGuesses[2] = Math.abs(currentid - fitter.initGuesses[0]);
        				fitter.userStage++;
        				JOptionPane.showMessageDialog(null, "Mouse over the top side of the peak and press F.");
        				break;
        			case 2:
        				fitter.initGuesses[3] = Math.abs(currentjd - fitter.initGuesses[1]);
        				fitter.userStage++;
        				fitter.initGuesses[5] = FieldOps.min(fitter.fftmag);
        				fitter.initGuesses[4] = fitter.fftmag[FieldOps.round(fitter.initGuesses[0])][FieldOps.round(fitter.initGuesses[1])];
        				System.out.println(Printer.arrayVertical(fitter.initGuesses));
        				double[] param = fitter.fitToFunction(false);
        				System.out.println(Printer.arrayVertical(param));
        				new LayerViewer(Layer.newLayer(t, fitter.values), dir, defaultSize);
        				new LayerViewer(Layer.newLayer(t, FieldOps.minus(t.data, fitter.values)), dir, defaultSize);
        				break;
        			}
        		}
        	}
        }
        bumpFitMI.addActionListener(new bumpFitAL(t,dir));
        
        //Analysis->Print correlation
        JMenuItem correlateMI = new JMenuItem("Print correlation");
        analysisMenu.add(correlateMI);
        class correlateAL implements ActionListener{
        	Layer t;
        	correlateAL(Layer tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
        		System.out.println(FieldOps.correlation(t.data, Layer.openFree(fc).data));
        	}
        }
        correlateMI.addActionListener(new correlateAL(t));
        
        
        //Analysis->Histogram
        JMenuItem histMI = new JMenuItem("Histogram of values");
        analysisMenu.add(histMI);
        class histAL implements ActionListener{
        	Layer t;
        	histAL(Layer tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
    			double[][] field = showFFT ? fftmag : t.data;
    			double[] layerData = FieldOps.getArray(field);
    			if (JOptionPane.showConfirmDialog(null, "Automatic histogram?") == JOptionPane.YES_OPTION){
    				double[] bins = ArrayOps.generateArrayInclBoth(FieldOps.min(field), FieldOps.max(field), (int)Math.sqrt(2*t.nx*t.ny));
    				int[] pixelHist = ArrayOps.getHistogram(layerData, bins);
    				GraphDrawerCart.plotGraph(bins, pixelHist);
    				ArrayOps.quicksort(layerData);
    				for (int i = 0; i < layerData.length; i++)
    					if (i % 100 == 0)
    						System.out.println("" + ((double)i)/layerData.length + "\t" + layerData[i]);
    			}
    			else {
    				double binwidth = Double.parseDouble(JOptionPane.showInputDialog("Enter desired bin width"));
    				double range = FieldOps.max(t.data) - FieldOps.min(t.data);
    				int nbins = (int)(range/binwidth)+1;
    				double[] bins = ArrayOps.generateArray(FieldOps.min(t.data), binwidth, nbins);
    				int[] pixelHist = ArrayOps.getHistogram(layerData, bins);
    				double[] frequency = new double [nbins];
    				GraphDrawerCart.plotGraph(bins, pixelHist);
    				for (int i = 0; i < bins.length; i++)
    				{
    					frequency[i] = ((double)pixelHist[i])/(t.nx*t.ny);
    					System.out.println("" + (bins[i]+binwidth/2) + "\t" + pixelHist[i] + "\t" + frequency[i]);
    				}
    			}
        	}
        }
        histMI.addActionListener(new histAL(t));
        
        //Analysis->enter line time
        JMenuItem lineTimeMI = new JMenuItem("Enter line time");
        analysisMenu.add(lineTimeMI);
        class lineTimeAL implements ActionListener{
        	Layer t;
        	lineTimeAL(Layer tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
				int choice = JOptionPane.showConfirmDialog(null, "To enter the line time (Nanonis) click Yes.\r\nTo enter ms/point (RHK) click No.");
				if (choice == JOptionPane.CANCEL_OPTION)
					return;
				if (choice == JOptionPane.YES_OPTION)
					t.lineTime = Double.parseDouble(JOptionPane.showInputDialog("Enter the time in s."));
				else t.lineTime = Double.parseDouble(JOptionPane.showInputDialog("Enter the ms/point in ms."))*t.nx/1000;
        	}
        }
        lineTimeMI.addActionListener(new lineTimeAL(t));
        
        //Analysis->toggle line cut tool
        JMenuItem lineCutMI = new JMenuItem("Toggle line cut tool");
        analysisMenu.add(lineCutMI);
        class lineCutAL implements ActionListener{
        	Layer t;
        	LayerViewer thisViewer;
        	lineCutAL(Layer tPrime, LayerViewer thisPrime){
        		t=tPrime;
        		thisViewer=thisPrime;
        	}
        	public void actionPerformed(ActionEvent e){
    			if (sectionDrawer == null)
    			{	
    				if (line == null) line = new double[][] {{0, 0}, {t.nx, t.ny}};
    				sectionDrawer = new LayerCutDrawer(thisViewer, 64, line);
    				sectionDrawer.setFC(fc);
    			}
    			else {
    				sectionDrawer.dispose();
    				sectionDrawer = null;
    			}
    			refresh = true;
    			repaint();
        	}
        }
        lineCutMI.addActionListener(new lineCutAL(t,this));
        
        ///////////////////////////////////////////////////////////////////////////
        //File-> Save line cut to a txt file
        JMenuItem saveCutasTxtMI = new JMenuItem("Save line cut as tab separated text file");
        analysisMenu.add(saveCutasTxtMI);
        class saveCutasTxtAL implements  ActionListener{
        	LayerViewer lv;
        	saveCutasTxtAL(LayerViewer lvPrime){
        		lv=lvPrime;
        	}
        	public void actionPerformed(ActionEvent e){
        		//if we have created a sectionDrawer for the line cut and the spectra have the correct
        		//length, the create a 2xN array of the x and y values of the line cut
        		if(line !=null){
        			//sectionDrawer= new LayerCutDrawer(lv, 128, line);
        			ColumnIO.writeTwoColumnsTxt(lv.sectionDrawer.s, lv.sectionDrawer.spec, FileOps.selectSave(fc).toString());
        		}
        	}
        }
        saveCutasTxtMI.addActionListener(new saveCutasTxtAL(this));
        ///////////////////////////////////////////////////////////////////////////////////
        
        
        //Analysis->bin the layer and get a histogram of impurity concentration
        JMenuItem impConcentrationMI = new JMenuItem("Histogram of impurity concentration in bins");
        analysisMenu.add(impConcentrationMI);
        class impConcentrationAL implements ActionListener{
        	Layer t;
        	impConcentrationAL(Layer tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
    			if (imps != null)
    			{	
    				int nbins = Integer.parseInt(JOptionPane.showInputDialog("How many bins?"));
    				int[] nimps = new int[nbins];
    				int[] impbin = new int[imps.size()];
    				int[][] bins = FieldOps.getPercentileBinsForField(t.data, nbins);
    				
    				int xi, yi;
    				for (int i = 0; i < imps.size(); i++)
    				{
    					xi = FieldOps.round(imps.get(i).getPosition()[0]);
    					yi = FieldOps.round(imps.get(i).getPosition()[1]);
    					xi = Math.min(xi, t.nx-1);
    					xi = Math.max(xi, 0);
    					yi = Math.min(yi, t.ny-1);
    					yi = Math.max(yi, 0);
    					nimps[bins[xi][yi]]++;
    					impbin[i] = bins[xi][yi]; 
    				}
    				
    				System.out.println("Bin\tN_imps");
    				for (int i = 0; i < nbins; i++)
    				{
    					System.out.println("" + i + "\t" + nimps[i] + "");
    				}
    				
    				System.out.println("imp\tbin\tx\ty");
    				for (int i = 0; i < imps.size(); i++)
    					System.out.println("" + i + "\t" + impbin[i] + "\t" + imps.get(i).getPosition()[0] + "\t" + imps.get(i).getPosition()[1]);
    			}
        	}
        }
        impConcentrationMI.addActionListener(new impConcentrationAL(t));
        
        //Analysis->Take the radial average
        JMenuItem radialAvgMI = new JMenuItem("Save the radial average ");
        analysisMenu.add(radialAvgMI);
        class radialAvgAL implements ActionListener{
        	Layer t;
        	radialAvgAL(Layer tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
        		int nBins = Math.min(Integer.parseInt(JOptionPane.showInputDialog("Enter number of radial bins to put the data in")),t.x.length);
        		double[] total= new double[nBins];
        		int[] count = new int[nBins];
        		for(int i=0;i<nBins;i++){
        			count[i]=0;
        			total[i]=0.0;
        		}
        		
        		double binSize = Math.sqrt(t.x.length * t.x.length + t.y.length * t.y.length)/nBins/2;
        		System.out.println("binSize: " + binSize);
        		
        		if(showFFT){
        			for(int j=0;j<t.x.length;j++){
        				for(int k=0;k<t.y.length;k++){
        					total[(int)(Math.sqrt((j-t.x.length/2) * (j-t.x.length/2)+(k-t.y.length/2)*(k-t.y.length/2))/binSize)]+=fftmag[j][k];
        					count[(int)(Math.sqrt((j-t.x.length/2) * (j-t.x.length/2)+(k-t.y.length/2)*(k-t.y.length/2))/binSize)]++;
        				}
        			}
        		}
        		else{
        			for(int j=0;j<t.x.length;j++){
        				for(int k=0;k<t.y.length;k++){
        					total[(int)(Math.sqrt((j-t.x.length/2)*(j-t.x.length/2)+(k-t.y.length/2)*(k-t.y.length/2))/binSize)]+=t.data[j][k];
        					count[(int)(Math.sqrt((j-t.x.length/2)*(j-t.x.length/2)+(k-t.y.length/2)*(k-t.y.length/2))/binSize)]++;
        				}
        			}
        		}
        		
        		double[][] radialAvg = new double[1][nBins];
        		for(int i=0;i<nBins;i++){
        			radialAvg[0][i]=total[i]/count[i];
        		}
        		System.out.println("Took radial average");
        		ColumnIO.writeTableCSV(radialAvg, FileOps.selectSave(fc).toString());
        	}
        }
        radialAvgMI.addActionListener(new radialAvgAL(t));
        
        //Manipulation menu
        JMenu manipMenu = new JMenu("Manipulation");
        menuBar.add(manipMenu);
        
        //Manip->undo to last backup
        JMenuItem undoMI = new JMenuItem("Undo to latest backup");
        manipMenu.add(undoMI);
    	class undoAL implements ActionListener{
    		Robo rob;
    		undoAL(Robo robPrime){
    			rob=robPrime;
    		}
    	    public void actionPerformed(ActionEvent e) {
    			rob.typeChar('Z');
           }
        }
    	undoMI.addActionListener(new undoAL(rob));
        
        //Manip->Noise or Background subtracting tool
        JMenuItem subtractToolMI = new JMenuItem("Noise or background subtracting tool");
        manipMenu.add(subtractToolMI);
    	class subtractToolAL implements ActionListener{
    		private Layer t;
    		subtractToolAL(Layer tPrime){
    			t=tPrime;
    		}
    	    public void actionPerformed(ActionEvent e) {
    			backup = Layer.newLayer(t, FieldOps.copy(t.data));
    			int o = RHKFileOps.getUserFittingChoice();
    			RHKFileOps.doFitting(t, o);
    			resetGraphics(true);
           }
        }
    	subtractToolMI.addActionListener(new subtractToolAL(t));
    	
    	//Manip->cut data outside current color scale
        JMenuItem cutScaleMI = new JMenuItem("Flatten datapoints outside current color scale");
        manipMenu.add(cutScaleMI);
    	class cutScaleAL implements ActionListener{
    		private Layer t;
    		cutScaleAL(Layer tPrime){
    			t=tPrime;
    		}
    	    public void actionPerformed(ActionEvent e) {
    	    	if(spec==null){
	    			backup = Layer.newLayer(t, FieldOps.copy(t.data));
	    			FieldOps.cutoffExtremesValue(t.data, scale.getMin(), scale.getMax());
	    			resetGraphics(true);
    	    	}
           }
        }
    	cutScaleMI.addActionListener(new cutScaleAL(t));
    	
    	//Manip->take second derivative
        JMenuItem takeSecDerivMI = new JMenuItem("Take 2nd derivative of data along lattice directions");
        manipMenu.add(takeSecDerivMI);
    	class takeSecDerivAL implements ActionListener{
    		private Layer t;
    		takeSecDerivAL(Layer tPrime){
    			t=tPrime;
    		}
    	    public void actionPerformed(ActionEvent e) {
    			AtomicCoordinatesSet latt = new AtomicCoordinatesSet(ColumnIO.getASCII(FileOps.selectOpen(fc).toString()));
//    			t = TopomapUtil.makeSymmetrizedFFT(t, latt, JOptionPane.showConfirmDialog(null, "Lattice is square?") == JOptionPane.YES_OPTION);
    			int index = Integer.parseInt(JOptionPane.showInputDialog("Enter your choice for the second derivative:\r\n"+
    					"0 - faa\r\n"+
    					"1 - fab\r\n"+
    					"2 - fba\r\n"+
    					"3 - fbb\r\n"+
    					"4 - 1/2(fab+fba)\r\n"+
    					"5 - 1/2(fab-fba) which is zero identically\r\n"+
    					"6 - 1/2(faa+fbb)\r\n"+
    					"7 - 1/2(faa-fbb)\r\n"));
    			
    			backup = Layer.newLayer(t, FieldOps.copy(t.data));
    			
    			AtomicCoordinatesSet recip = latt.getReciprocalLattice();
    			double[][] units = FieldOps.getGrahamSchmidtUnitVectors(recip.getA(), recip.getB());
    			FieldOps.copy(FieldOps.getDirectional2ndDerivativeStuff(t.data, units[0], units[1])[index],t.data);
    			resetColorScale();
    			formImage();
    			repaint();
           }
        }
    	takeSecDerivMI.addActionListener(new takeSecDerivAL(t));
    	
        //Manipulation-> take autocorrelation
        JMenuItem autocorrMI = new JMenuItem("Take autocorrelation");
        manipMenu.add(autocorrMI);
        class autocorrAL implements ActionListener{
        	Layer t;
        	autocorrAL(Layer tPrime){
        		t=tPrime;
        	}
        	public void actionPerformed(ActionEvent e){
    	    	if(spec==null){
	    			backup = Layer.newLayer(t, FieldOps.copy(t.data));
	    			t.data = FieldOps.getAutocorrelationFourier(t.data);
	    			resetGraphics(true);
    	    	}
        	}
        }
        autocorrMI.addActionListener(new autocorrAL(t));
    	
    	//Point spectra menu
    	JMenu specMenu = new JMenu("Point Spectra");
    	menuBar.add(specMenu);
    	
    	//Spec-> open point spectra
        JMenuItem openSpecMI = new JMenuItem("Open point spectra");
        specMenu.add(openSpecMI);
        openSpecMI.addActionListener(new ActionListener(){
      	  	public void actionPerformed(ActionEvent e) {
    			//open point spectra;
    			spec = PointSpectra.open(fc);
    			if (spec != null)
    			{
    				splitPos = spec.splitByPosition();
    				Robo.wait(50);
    				if (!showingSpec) switchSpectraDrawing();
    				else{
    					refresh = true;
    					repaint();
    				}
    			}
    		}
        });
        
        //Spec-> show point spectra
        JMenuItem showSpecMI = new JMenuItem("Toggle showing point spectra");
        specMenu.add(showSpecMI);
        showSpecMI.addActionListener(new ActionListener(){
      	  	public void actionPerformed(ActionEvent e) {
      	  		if(spec!=null)
      	  			switchSpectraDrawing();
    		}
        });
        
        //Spec-> toggle spectrum viewer
        JMenuItem specViewerMI = new JMenuItem("Toggle spectrum viewer");
        specMenu.add(specViewerMI);
        specViewerMI.addActionListener(new ActionListener(){
      	  	public void actionPerformed(ActionEvent e) {
      	  		if(spec!=null){
      				if (specView == null)
      					specView = new SpectraDrawer(splitPos[selectedSpectra], null);
      				else
      				{
      					specView.dispose();
      					specView = null;
      				}
      	  		}
    		}
        });
        
    	//Spec->Print selected spectrum
        JMenuItem printSpecMI = new JMenuItem("Print selected spectrum");
        specMenu.add(printSpecMI);
        printSpecMI.addActionListener(new ActionListener(){
      	  	public void actionPerformed(ActionEvent e) {
      	  		if(spec!=null)
      	  			System.out.println(Printer.getTable(new double[][] {spec.v, splitPos[selectedSpectra].average}));
      	  	}
        });
        
        //Spec-> create gif
        JMenuItem specToGifMI = new JMenuItem("Create GIF from spectra");
        specMenu.add(specToGifMI);
        specToGifMI.addActionListener(new ActionListener(){
      	  	public void actionPerformed(ActionEvent e) {
      	  		if(spec!=null){
      				BufferedImage[] gifstack = new BufferedImage [splitPos.length];
	    			BufferedImage[] imstack = new BufferedImage [splitPos.length];
	    			BufferedImage[] graphStack = new BufferedImage [splitPos.length];
	    			int sizeFactor = GraphDrawerCart.DEFAULT_WIDTH/image.getWidth();
	    			if (sizeFactor == 0)
	    				sizeFactor = 1;
	    			
	    			BufferedImage baseImage, bigImage, graphImage;
	    			for (int i = 0; i < splitPos.length; i++)
	    			{
	    				baseImage = ImageEditing.getCopy(image);
	    				imstack[i] = baseImage;
	    				drawSpecSites(baseImage, i);
	    				graphImage = GraphDrawerCart.drawPlot(spec.v, splitPos[i].average, ArrayOps.min(spec.v), ArrayOps.max(spec.v), ArrayOps.min(spec.data), ArrayOps.max(spec.data));
	    				graphStack[i] = graphImage;
	    				bigImage = ImageEditing.getEnlarged(baseImage,sizeFactor);
	    				gifstack[i] = new BufferedImage(GraphDrawerCart.DEFAULT_WIDTH + sizeFactor*image.getWidth(), GraphDrawerCart.DEFAULT_HEIGHT, BufferedImage.TYPE_INT_RGB);
	    				ImageEditing.setColor(gifstack[i], Color.WHITE);
	    				ImageEditing.copyInto(graphImage, gifstack[i], 0, 0);
	    				ImageEditing.copyInto(bigImage, gifstack[i], graphImage.getWidth(), graphImage.getHeight() - bigImage.getHeight());
	    			}
	    			String path = FileOps.selectSave(fc).toString();
	    			String wholepath = path + "_both";
	    			String graphpath = path + "_graph";
	    			String impath = path + "_image";
	    			GifSequenceWriter.writeGifSequence(gifstack, 500, wholepath.endsWith(".gif") ? wholepath : wholepath + ".gif");
	    			GifSequenceWriter.writeGifSequence(graphStack, 500, graphpath.endsWith(".gif") ? graphpath : graphpath + ".gif");
	    			GifSequenceWriter.writeGifSequence(imstack, 500, impath.endsWith(".gif") ? impath : impath + ".gif");
      	  		}
      	  	}
        });
        
        //Spec->increment spectrum position
        JMenuItem plusSpecMI = new JMenuItem("Increment spectrum position");
        specMenu.add(plusSpecMI);
        plusSpecMI.addActionListener(new ActionListener(){
      	  	public void actionPerformed(ActionEvent e) {
      	  		if(spec!=null){
      				selectedSpectra--;
	    			selectedSpectra += splitPos.length;
	    			selectedSpectra %= splitPos.length;
	    		    refresh = true;
	    		    repaint();
	    		    setTitle("" + selectedSpectra);
	    			if (specView != null)
	    				specView.resetTotally(splitPos[selectedSpectra]);
      	  		}
    		}
        });
        
        //Spec->dencrement spectrum position
        JMenuItem minusSpecMI = new JMenuItem("Decrement spectrum position");
        specMenu.add(minusSpecMI);
        minusSpecMI.addActionListener(new ActionListener(){
      	  	public void actionPerformed(ActionEvent e) {
      	  		if(spec!=null){
      				selectedSpectra++;
	    			selectedSpectra += splitPos.length;
	    			selectedSpectra %= splitPos.length;
	    		    refresh = true;
	    		    repaint();
	    		    setTitle("" + selectedSpectra);
	    			if (specView != null)
	    				specView.resetTotally(splitPos[selectedSpectra]);
      	  		}
    		}
        });
        
        //Spec->treat rows as separate spectra
        JMenuItem layerSpecMI = new JMenuItem("Manipulations assuming each row is a spectrum");
        specMenu.add(layerSpecMI);
    	class layerSpecAL implements ActionListener{
    		Robo rob;
    		layerSpecAL(Robo robPrime){
    			rob=robPrime;
    		}
    	    public void actionPerformed(ActionEvent e) {
    	    	rob.typeChar('R');
    	    }
        }
    	layerSpecMI.addActionListener(new layerSpecAL(rob));
    	
//    	////////////////////////////////////
//    	//Zoom-> Gives ability to zoom in and then zoom out
//		JMenu zoomMenu = new JMenu("Zoom");
//		menuBar.add(zoomMenu);
//    	
//		//Zoom in-> zoom in by selecting a point and square box dimension;
//		JMenuItem zoomInMI = new JMenuItem("Zoom In");
//		zoomMenu.add(zoomInMI);
//		class zoomInAL implements ActionListener{
//			Layer t;
//			zoomInAL(Layer tPrime){
//				t=tPrime;
//			}
//			public void actionPerformed(ActionEvent e){
//				t.zoomIn();
//				java.awt.Dimension screenSize=Toolkit.getDefaultToolkit().getScreenSize();
//				sizeratio=(int)(screenSize.getHeight()/t.nx);
//			}
//		}
//		zoomInMI.addActionListener(new zoomInAL(t));
//		
//		//Zoom out -> zoom out to the previous size
//		JMenuItem zoomOutMI = new JMenuItem("Zoom Out");
//		zoomMenu.add(zoomOutMI);
//		class zoomOutAL implements ActionListener{
//			Layer t;
//			zoomOutAL(Layer tPrime){
//				t=tPrime;
//			}
//			public void actionPerformed(ActionEvent e){
//				//t.zoomOut();
//				//sizeratio=change ratio
//			}
//		}
//		zoomOutMI.addActionListener(new zoomOutAL(t));
//		
//		/////////////////////////////////////////////////////////
		
    	
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
		bounds = new double[] {ArrayOps.min(drawField), ArrayOps.max(drawField)};
			
		fmax = bounds[1]; fmin = bounds[0];
		fdelta = fmax - fmin;
		setTitle("Range: [" + fmin + ", " + fmax + "]" + "     " + "p = " + para);
	}
	public void resetColorScale(double downnum, double upnum, boolean hardLimits)
	{	
		if (!hardLimits){
			csdown = downnum;
			csup = upnum;
		}
		if (!expandedScaleRange)
		{
			scalemax = fmin + fdelta*upnum;
			scalemin = fmin + fdelta*downnum;
		}
		else //in this case make the color scale the middle 0.5 of the range 0 to 1.
		{
			scalemax = fmin + fdelta*(upnum-0.25)*2;
			scalemin = fmin + fdelta*(downnum-0.25)*2;			
		}
		
		if (hardLimits) scale = ColorScales.getNew(downnum, upnum, currentCScale);
		else scale = ColorScales.getNew(scalemin, scalemax, currentCScale);
	}
	public void resetColorScale()
	{
		if (real) scale = ColorScales.getNew(scalemin, scalemax, currentCScale);
		else cscale = new ColorScales.MYC2d(scalemax, scalemin, 2*Math.PI);
	}

	public void formImage()
	{
		if (real) SRAW.writeImage(image, drawField, scale, sizeratio);
		else SRAW.writeImage(image, drawFieldC, cscale);
		
		if (drawingTP && !showFFT) ImageEditing.setTrueToColor(truePixImage, image, ColorScales.getUnusedColor(currentCScale)) ;
	}
	public void showWindow()
	{
		setSize(WIDTH, HEIGHT);
		repaint();
		setVisible(true);
		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}	
	public void paint(Graphics g)
	{
		if (refresh){
			g.clearRect(0, 0, 2000, 2000);
			g.drawImage(image, ox, oy, null);
			refresh = false;
		}
		if (showingSpec && spec != null)
		{
			drawSpecSites(g);
		}
		if (imps != null && !showFFT)
		{
			drawImpMarks(g, ColorScales.getUnusedColor(currentCScale), 0);
		}
		if (imps2 != null && !showFFT)
		{
			drawImpMarks(g, ColorScales.getUnusedColor2(currentCScale), 1);
		}

		if (sectionDrawer != null) drawLine(g, ColorScales.getUnusedColor(currentCScale));
		
		menuBar.paintImmediately(0,0,WIDTH,50);
		drawText(g);
	}
	public void drawLine(Graphics g, Color c)
	{
		g.setColor(c);
		g.drawString("A", (int)(line[0][0]*imageOverField)+ox, (int)(line[0][1]*imageOverField)+oy);
		g.drawString("B", (int)(line[1][0]*imageOverField)+ox, (int)(line[1][1]*imageOverField)+oy);
		g.drawLine((int)(line[0][0]*imageOverField)+ox, (int)(line[0][1]*imageOverField)+oy, (int)(line[1][0]*imageOverField)+ox, (int)(line[1][1]*imageOverField)+oy);
	}
	public void drawImpMarks(Graphics g, Color c, int k)
	{
		ArrayList<Impurity> imps = k == 0 ? this.imps : imps2;
		double[] r = new double[2];
		g.setColor(c);
		for (int i = 0; i < imps.size(); i++)
		{
			r = imps.get(i).getPosition();
			drawPlus(g, (int)(r[0]*sizeratio)+ox, (int)(r[1]*sizeratio)+oy, 5);
		}
	}

	public void drawSpecSites(Graphics g)
	{
		double[] r = new double [2];
		int x, y;
		for (int i = 0; i < splitPos.length; i++)
		{
			t.putPixelCoords(splitPos[i].x[0], splitPos[i].y[0], r);
//			r = new double[] {splitPos[i].x[0], splitPos[i].y[0]};
			x = screenX(r[0]);
			y = screenY(r[1]);
			System.out.println("" + i + "\t" + Printer.vectorP(new int [] {x, y}) + "\t" + Printer.vectorP(new double[] {splitPos[i].x[0], splitPos[i].y[0]}));;
			if (i == selectedSpectra) g.setColor(Color.MAGENTA);
			else g.setColor(ColorScales.getUnusedColor(currentCScale));
			drawPlus(g, x, y, 5);
			if (splitPos[i].nspec > 1)
				g.drawString("" + splitPos[i].nspec, x+10, y);
		}
	}
	public void drawSpecSites(BufferedImage image, int spec)
	{
		double[] r = new double [2];
		int x, y;
		Graphics g = image.getGraphics();
		for (int i = 0; i < splitPos.length; i++)
		{
			t.putPixelCoords(splitPos[i].x[0], splitPos[i].y[0], r);
			x = FieldOps.round(r[0]*sizeratio);
			y = FieldOps.round(r[1]*sizeratio);
			if (i == spec) g.setColor(Color.MAGENTA);
			else g.setColor(ColorScales.getUnusedColor(currentCScale));
			drawPlus(g, x, y, 5);
		}
	}
	public static void drawPlus(Graphics g, int x, int y, int halflength)
	{
		g.drawLine(x-halflength, y, x+halflength, y);
		g.drawLine(x, y-halflength, x, y+halflength);
	}
	public void drawText(Graphics g)
	{
		g.setColor(Color.BLACK);
		g.clearRect(writepoint[0], writepoint[1]-linesize, 500, 800);
		String it = "";
		it += "Current mouse position in map: (" + currenti + ", " + currentj + ")" + "\r\n";
		it += "With respect to the center, that is (" + (currenti-t.nx/2) + ", " + (currentj-t.ny/2) + ") \r\n";
		if (!showFFT)
			it += "In metric units that is " + Printer.vectorPFormat(t.getMetricCoords(currenti, currentj)) + "\r\n"; 
		else
			it += "In metric units that is " + Printer.vectorPFormat(t.getFourierMetricCoords(currenti, currentj)) + "\r\n"; 
		it += "This vector has magnitude " + NumFormat.scientific(Distance.distance((currenti-t.nx/2), (currentj-t.ny/2)),3) + " \r\n";
		if (!showFFT)
			it += "(metric " + NumFormat.scientific(Distance.distance(t.getMetricCoords(t.nx/2, t.ny/2), t.getMetricCoords(currenti, currentj)),3) + ")\r\n";
		else
			it += "(metric " + NumFormat.scientific(Complex.mag(t.getFourierMetricCoords(currenti, currentj)),3) + ")\r\n";
		it += "And corresponding angle " + NumFormat.scientific(Math.toDegrees(FieldOps.atan((currenti-t.nx/2), (currentj-t.ny/2))),3) + "degrees. \r\n";
		it += "The value of the drawn field there is " + NumFormat.scientific((showFFT ? fftmag[currenti][currentj] : t.data[currenti][currentj]),3) + ".\r\n";
		if (showFFT) it += "Wavelength = " + NumFormat.scientific(t.xLength/(Distance.distance((currenti-t.nx/2), (currentj-t.ny/2))),3) + " m\r\n";
		if (showFFT && fftz != null) it += "Phase is " + NumFormat.scientific(FieldOps.atan(fftz[currenti][currentj][0], fftz[currenti][currentj][1]),3);
		String[] lines = it.split("\r\n");
		for (int i = 0; i < lines.length; i++)
			g.drawString(lines[i], writepoint[0], writepoint[1] + i*linesize);
	}
	public void resetGraphics(boolean redoColorScale)
	{
		if (!showFFT) drawField = t.data;
		else drawField = fftmag;
		zoomLevel = 0;
		zoomFactor = 1;
		sizeratio = 1;
		imageOverField = 1;
		
		while (drawField.length > defaultSize && drawField.length/2 >= defaultSize)
		{
			this.drawField = FieldOps.reduce(2, this.drawField);
			sx = this.drawField.length;
			zoomLevel++;
			zoomFactor *= 2;
			imageOverField /= 2;
		}
		//Can change this if you want size ratios other than 2^n 
		while (sizeratio*sx < defaultSize)
		{
			sizeratio++;
			imageOverField++;
		}
		setFieldInfo();
		if (redoColorScale)
			resetColorScale(0, 1, false);
	    this.formImage();
	    refresh = true;
	    repaint();
	}
	public void switchSpectraDrawing()
	{
		showingSpec = !showingSpec;
		System.out.println(showingSpec);
		refresh = true;
		repaint();
	}
	//if 
	public void keyPressed(KeyEvent arg0) {
	}
	public void keyReleased(KeyEvent arg0) {
	}
	public void keyTyped(KeyEvent arg0) {
		System.out.println(arg0.getKeyChar());
		if (arg0.getKeyChar() == 'f')
		{
			//take Fourier transform.
			if (fftmag == null || switchLog || true)
			{
				fftmag = new double [t.nx][t.ny];
				FFTOps.obtainFFTmagCent(t.data, fftmag);
				if (ftlog)
					FieldOps.log(fftmag);
				setTitle("Done");
				switchLog = false;
			}
			showFFT = !showFFT;
			resetGraphics(true);
		}
		if (arg0.getKeyChar() == 'F' && fitter == null)
		{
			fitter = new TwoDGaussianFreeFitter(showFFT ? fftmag : t.data);
			JOptionPane.showMessageDialog(null, "Mouse over the peak center and press F.");
		}
		else if (arg0.getKeyChar() == 'F' && fitter != null)
		{
			switch(fitter.userStage)
			{
			case 0:
				fitter.initGuesses[0] = currentid;
				fitter.initGuesses[1] = currentjd;
				fitter.userStage++;
				JOptionPane.showMessageDialog(null, "Mouse over the right side of the peak and press F.");
				break;
			case 1:
				fitter.initGuesses[2] = Math.abs(currentid - fitter.initGuesses[0]);
				fitter.userStage++;
				JOptionPane.showMessageDialog(null, "Mouse over the top side of the peak and press F.");
				break;
			case 2:
				fitter.initGuesses[3] = Math.abs(currentjd - fitter.initGuesses[1]);
				fitter.userStage++;
				fitter.initGuesses[5] = FieldOps.min(fitter.fftmag);
				fitter.initGuesses[4] = fitter.fftmag[FieldOps.round(fitter.initGuesses[0])][FieldOps.round(fitter.initGuesses[1])];
				System.out.println(Printer.arrayVertical(fitter.initGuesses));
				double[] param = fitter.fitToFunction(false);
				System.out.println(Printer.arrayVertical(param));
				new LayerViewer(Layer.newLayer(t, fitter.values), dir, defaultSize);
				new LayerViewer(Layer.newLayer(t, FieldOps.minus(t.data, fitter.values)), dir, defaultSize);
				break;
			}
		}
		if (arg0.getKeyChar() == ' '){
			refresh=true;
    		repaint();
		}
		
		if (arg0.getKeyChar() == 'p')
		{
			String input = JOptionPane.showInputDialog("Enter the color scale range (comma separated).", "" + scale.getMin() + "," + scale.getMax());
			double min = Double.parseDouble(input.split(",")[0].trim());
			double max = Double.parseDouble(input.split(",")[1].trim());
			ColorScale1D cs = ColorScales.getNew(min, max, currentCScale);

			int[] sizeRatio = new int [2];
			sizeRatio[0] = Integer.parseInt(JOptionPane.showInputDialog("By how much would you like to blow up the image? (X)", "" + sizeratio));
			sizeRatio[1] = Integer.parseInt(JOptionPane.showInputDialog("By how much would you like to blow up the image? (Y)", "" + sizeratio));
			BufferedImage fullImage = ImageEditing.getBufferedImage(showFFT ? fftmag : t.data, cs);
			BufferedImage im = ImageEditing.enlargeBasicStretch(fullImage, sizeRatio[0], sizeRatio[1]);
			ImageEditing.copyToClipboard(im);
		}
		if (arg0.getKeyChar() == 'P')
		{
			backup = Layer.newLayer(t, FieldOps.copy(t.data));
			int o = RHKFileOps.getUserFittingChoice();
			RHKFileOps.doFitting(t, o);
			resetGraphics(true);
		}
		if (arg0.getKeyChar() == 'x' && spec == null)
		{
			backup = Layer.newLayer(t, FieldOps.copy(t.data));
			FieldOps.cutoffExtremesValue(t.data, scale.getMin(), scale.getMax());
			resetGraphics(true);
		}
		if (arg0.getKeyChar() == 'X') //crop this
		{
			String imput2 = JOptionPane.showInputDialog("Enter coordinates of the upper-left corner of the box, comma separated.");
			int xi = Integer.parseInt(imput2.split(",")[0].trim());
			int yi = Integer.parseInt(imput2.split(",")[1].trim());
			imput2 = JOptionPane.showInputDialog("Enter the width and height of the box, comma separated.");
			int xf = xi + Integer.parseInt(imput2.split(",")[0].trim());
			int yf = yi  + Integer.parseInt(imput2.split(",")[1].trim());
			LayerViewer.show(LayerUtil.cropLayer(t, xi, xf, yi, yf), this.sx, false);
		}
		if (arg0.getKeyChar() == 'S')
		{
			Layer.writeBIN(t, fc);
		}
		if (arg0.getKeyChar() == 'A')
		{
			ColumnIO.writeTable(t.data, FileOps.selectSave(fc).toString());
		}
		if (arg0.getKeyChar() == 's' && spec != null)
		{
			System.out.println(Printer.getTable(new double[][] {spec.v, splitPos[selectedSpectra].average}));
		}
		//Create spectra GIF
		if (arg0.getKeyChar() == 'G' && spec != null)
		{
			BufferedImage[] gifstack = new BufferedImage [splitPos.length];
			BufferedImage[] imstack = new BufferedImage [splitPos.length];
			BufferedImage[] graphStack = new BufferedImage [splitPos.length];
			int sizeFactor = GraphDrawerCart.DEFAULT_WIDTH/image.getWidth();
			if (sizeFactor == 0)
				sizeFactor = 1;
			
			BufferedImage baseImage, bigImage, graphImage;
			for (int i = 0; i < splitPos.length; i++)
			{
				baseImage = ImageEditing.getCopy(image);
				imstack[i] = baseImage;
				drawSpecSites(baseImage, i);
				graphImage = GraphDrawerCart.drawPlot(spec.v, splitPos[i].average, ArrayOps.min(spec.v), ArrayOps.max(spec.v), ArrayOps.min(spec.data), ArrayOps.max(spec.data));
				graphStack[i] = graphImage;
				bigImage = ImageEditing.getEnlarged(baseImage,sizeFactor);
				gifstack[i] = new BufferedImage(GraphDrawerCart.DEFAULT_WIDTH + sizeFactor*image.getWidth(), GraphDrawerCart.DEFAULT_HEIGHT, BufferedImage.TYPE_INT_RGB);
				ImageEditing.setColor(gifstack[i], Color.WHITE);
				ImageEditing.copyInto(graphImage, gifstack[i], 0, 0);
				ImageEditing.copyInto(bigImage, gifstack[i], graphImage.getWidth(), graphImage.getHeight() - bigImage.getHeight());
			}
			String path = FileOps.selectSave(fc).toString();
			String wholepath = path + "_both";
			String graphpath = path + "_graph";
			String impath = path + "_image";
			GifSequenceWriter.writeGifSequence(gifstack, 500, wholepath.endsWith(".gif") ? wholepath : wholepath + ".gif");
			GifSequenceWriter.writeGifSequence(graphStack, 500, graphpath.endsWith(".gif") ? graphpath : graphpath + ".gif");
			GifSequenceWriter.writeGifSequence(imstack, 500, impath.endsWith(".gif") ? impath : impath + ".gif");
			
		}
		if (arg0.getKeyChar() == 'G' && spec == null)
		{
			System.out.println(FieldOps.correlation(t.data, Layer.openFree(fc).data));
		}
		if (arg0.getKeyChar() == 'j')
			expandedScaleRange = !expandedScaleRange;
		if (arg0.getKeyChar() == 'q')
		{
			currentCScale++;
			currentCScale += ColorScales.NSCALES;
			currentCScale %= ColorScales.NSCALES;
			resetColorScale(csdown, csup, false);
			formImage();
			setTitle("Color scale: " + (currentCScale+1) + " out of " + ColorScales.NSCALES);
			refresh = true;
			repaint();
		}
		if (arg0.getKeyChar() == 'e')
		{
			currentCScale--;
			currentCScale += ColorScales.NSCALES;
			currentCScale %= ColorScales.NSCALES;
			resetColorScale(csdown, csup, false);
			formImage();
			setTitle("Color scale: " + (currentCScale+1) + " out of " + ColorScales.NSCALES);
			refresh = true;
			repaint();
		}
		if (arg0.getKeyChar() == 'E') //export csv
		{
			ColumnIO.writeTableCSV(t.data, FileOps.selectSave(fc).toString());
		}
		if (arg0.getKeyChar() == 'g'){
			ftlog = !ftlog;
			switchLog = true;
		}
		if (arg0.getKeyChar() == 'o')
		{
			//open point spectra;
			spec = PointSpectra.open(fc);
			if (spec != null)
			{
				splitPos = spec.splitByPosition();
				Robo.wait(50);
				if (!showingSpec) switchSpectraDrawing();
				else{
					refresh = true;
					repaint();
				}
			}
		}
		if (arg0.getKeyChar() == 'h' && spec != null)
		{
			//show point spectra
			switchSpectraDrawing();
		}
		if (arg0.getKeyChar() == 'h' && spec == null)
		{
			double[][] field = showFFT ? fftmag : t.data;
			double[] layerData = FieldOps.getArray(field);
			if (JOptionPane.showConfirmDialog(null, "Automatic histogram?") == JOptionPane.YES_OPTION){
				double[] bins = ArrayOps.generateArrayInclBoth(FieldOps.min(field), FieldOps.max(field), (int)Math.sqrt(2*t.nx*t.ny));
				int[] pixelHist = ArrayOps.getHistogram(layerData, bins);
				GraphDrawerCart.plotGraph(bins, pixelHist);
				ArrayOps.quicksort(layerData);
				for (int i = 0; i < layerData.length; i++)
					if (i % 100 == 0)
						System.out.println("" + ((double)i)/layerData.length + "\t" + layerData[i]);
			}
			else {
				double binwidth = Double.parseDouble(JOptionPane.showInputDialog("Enter desired bin width"));
				double range = FieldOps.max(t.data) - FieldOps.min(t.data);
				int nbins = (int)(range/binwidth)+1;
				double[] bins = ArrayOps.generateArray(FieldOps.min(t.data), binwidth, nbins);
				int[] pixelHist = ArrayOps.getHistogram(layerData, bins);
				double[] frequency = new double [nbins];
				GraphDrawerCart.plotGraph(bins, pixelHist);
				for (int i = 0; i < bins.length; i++)
				{
					frequency[i] = ((double)pixelHist[i])/(t.nx*t.ny);
					System.out.println("" + (bins[i]+binwidth/2) + "\t" + pixelHist[i] + "\t" + frequency[i]);
				}
			}
		}
		//switching spectra position
		if (arg0.getKeyChar() == 'z' && spec != null && showingSpec)
		{
			selectedSpectra--;
			selectedSpectra += splitPos.length;
			selectedSpectra %= splitPos.length;
		    refresh = true;
		    repaint();
		    setTitle("" + selectedSpectra);
			if (specView != null)
				specView.resetTotally(splitPos[selectedSpectra]);
		}
		if (arg0.getKeyChar() == 'x' && spec != null && showingSpec)
		{
			selectedSpectra++;
			selectedSpectra += splitPos.length;
			selectedSpectra %= splitPos.length;
		    refresh = true;
		    repaint();
		    setTitle("" + selectedSpectra);
			if (specView != null)
				specView.resetTotally(splitPos[selectedSpectra]);
		}
		if (arg0.getKeyChar() == 't' && spec != null && showingSpec)
		{
			//turn on and off the spectra viewer
			if (specView == null)
				specView = new SpectraDrawer(splitPos[selectedSpectra], null);
			else
			{
				specView.dispose();
				specView = null;
			}
		}
		if (arg0.getKeyChar() == 'T')
		{
//			System.out.println(t.lineTime);
			if (t.lineTime == 0)
			{
				int choice = JOptionPane.showConfirmDialog(null, "To enter the line time (Nanonis) click Yes.\r\nTo enter ms/point (RHK) click No.");
				if (choice == JOptionPane.CANCEL_OPTION)
					return;
				if (choice == JOptionPane.YES_OPTION)
					t.lineTime = Double.parseDouble(JOptionPane.showInputDialog("Enter the time in s."));
				else t.lineTime = Double.parseDouble(JOptionPane.showInputDialog("Enter the ms/point in ms."))*t.nx/1000;
			}
		}
		if (arg0.getKeyChar() == 'Z')
		{
			if (backup != null)
			{
				t = Layer.newLayer(backup, backup.data);
				resetGraphics(true);
			}
		}
		if (arg0.getKeyChar() == 'l') //toggle the line cut tool
		{
			if (sectionDrawer == null)
			{	
				if (line == null) line = new double[][] {{0, 0}, {t.nx, t.ny}};
				//number of points in a line cut is second argument
				sectionDrawer = new LayerCutDrawer(this, 64, line);
				sectionDrawer.setFC(fc);
			}
			else {
				sectionDrawer.dispose();
				sectionDrawer = null;
			}
			refresh = true;
			repaint();
		}
		if (arg0.getKeyChar() == 'B')
		{
			if (truePix == null)
			{
				truePix = ImageEditing.loadFromImage(FileOps.selectOpen(fc).toString(), Color.WHITE);
				truePixImage = FieldOps.resizeByFactor(truePix, imageOverField);
				drawingTP = true;
				refresh = true;
				formImage();
				repaint();
			}
			else
			{
				drawingTP = !drawingTP;
				refresh = true;
				formImage();
				repaint();
			}
		}
		//print to the cliboard a blown-up crop of the current image
		if (arg0.getKeyChar() == 'N')
		{
//			System.out.println(sizeratio);
			String input = JOptionPane.showInputDialog("Enter upper left corner of the box, spearated by commas, with respect to the center.");
			int x = Integer.parseInt(input.split(",")[0].trim()) + t.nx/2;
			int y = Integer.parseInt(input.split(",")[1].trim()) + t.ny/2;

			int boxSize = Integer.parseInt(JOptionPane.showInputDialog("Size of the box?"));
			BufferedImage temp = ImageEditing.getSubset(image, sizeratio*x, boxSize, sizeratio*y, boxSize);
			int blowupfactor = Integer.parseInt(JOptionPane.showInputDialog("Blow-up factor?", "" + 1));
			ImageEditing.copyToClipboard(ImageEditing.getEnlarged(temp, blowupfactor));
			
		}
		
		if (arg0.getKeyChar() == 'R') //Treat the layer as if it is a point spectra object:
		{
//			showCustomLayers();
			String menu = "What would you like to do?\r\n" + 
					"0 - smooth each line with a Gaussian\r\n" + 
					"1 - replace each line with -1* its second derivative\r\n" + 
					"2 - Replace each line with its derivative\r\n" + 
					"3 - adaptively fit some lines to parabolas and print the results\r\n";
			int choice = Printer.getAnInt(menu);
			backup = Layer.newLayer(t, FieldOps.copy(t.data));
			PointSpectra ps = t.toSpectra();
			if (choice == 0){
				ps = SpectraUtil.getGaussSmoothed(ps, Printer.getADouble("Smoothing length?"));
				t = ps.toLayer();
				resetGraphics(true);
			}
			if (choice == 1){
				for (int i = 0; i < ps.nspec; i++)
					ps.data[i] = ArrayOps.getDerivative(ArrayOps.getDerivative(ps.data[i]));
				FieldOps.negate(ps.data);
				t = ps.toLayer();
				resetGraphics(true);
			}
			if (choice == 2){
				for (int i = 0; i < ps.nspec; i++)
					ps.data[i] = ArrayOps.getDerivative(ps.data[i]);
//				FieldOps.negate(ps.data);
				t = ps.toLayer();
				resetGraphics(true);
			}
			if (choice == 3){
				int jstart = Printer.getAnInt("First line?", 0);
				int jend = Printer.getAnInt("Last line inclusive?", t.ny-1);
				int delta = jstart > jend ? -1 : 1;
				double fcent = Printer.getADouble("Fitting center in metric units?");
				double hwidth = Printer.getADouble("Half width in metric units?");
				double[][] results = SpectraUtil.adaptivelyFitAllSpectraToAParabola(ps, jstart, jend, fcent, hwidth);
				String[] lines = new String [results[0].length+2];
				lines[0] = "jstart = " + jstart + "\tjend = " + jend + "\tfcent = " + fcent + "\thwidth = " + hwidth;
				lines[1] = "Index\tFit center\t";
				for (int i = 0; i < results[0].length; i++)
				{
					lines[i+2] = "" + i + "\t" + results[0][i];
				}
				ColumnIO.writeLines(lines, FileOps.selectSave(fc).toString());
				BufferedImage crosses = new BufferedImage(image.getWidth(), image.getHeight(), BufferedImage.TYPE_INT_ARGB);
				ImageEditing.renderTransparent(crosses);
				Graphics g = crosses.getGraphics();
				g.setColor(ColorScales.getUnusedColor(scale));
				for (int i = 0; i < results[0].length; i++)
				{
					int y = screenY(jstart + (i*delta) + 0.5) - oy;
					int x = screenX(t.getPixelCoords(results[0][i], y)[0]) - ox;
					drawPlus(g, x, y, 2);
				}
				ImageEditing.copyToClipboard(crosses);

			}
			
		}
		if (arg0.getKeyChar() == 'r')
		{
			//write to a text file.
			boolean transpose = JOptionPane.showConfirmDialog(null, "Transpose the data?") == JOptionPane.YES_OPTION;
			double[][] table = new double [t.ny+1][t.nx+1];
			for (int i = 0; i < t.nx; i++)
				table[0][i+1] = t.x[i];
			for (int j = 0; j < t.ny; j++)
				table[j+1][0] = t.y[j];
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
					table[j+1][i+1] = t.data[i][j];
			if (!transpose)
				table = FieldOps.transpose(table);
			
			FileOps.writeTableASCII(fc, table);
			
		}
		if (arg0.getKeyChar() == 'M')
		{
			int hugeRatio = 18;
			
			BufferedImage huge = ImageEditing.getBufferedImage(t.data, scale);
			huge = ImageEditing.getEnlarged(huge, hugeRatio);
			Graphics g = huge.getGraphics();
			g.setColor(ColorScales.getUnusedColor(scale));
			g.setFont(new Font("Lucida Console", Font.PLAIN, hugeRatio/3));
			for (int i = 0; i < t.nx; i++)
				for (int j = 0; j < t.ny; j++)
				{
					double offset = (i % 2) == 0 ? 0.7 : 0.4;
					g.drawString(NumFormat.scientific(t.data[i][j], 2), i*hugeRatio, (int)((j+offset)*hugeRatio));
				}
			ImageEditing.copyToClipboard(huge);
		}
		//Bin the layer and output a histogram of how many impurities are in each bin
		if (arg0.getKeyChar() == 'b' && imps != null)
		{
			int nbins = Integer.parseInt(JOptionPane.showInputDialog("How many bins?"));
			int[] nimps = new int[nbins];
			int[] impbin = new int[imps.size()];
			int[][] bins = FieldOps.getPercentileBinsForField(t.data, nbins);
			
			int xi, yi;
			for (int i = 0; i < imps.size(); i++)
			{
				xi = FieldOps.round(imps.get(i).getPosition()[0]);
				yi = FieldOps.round(imps.get(i).getPosition()[1]);
				xi = Math.min(xi, t.nx-1);
				xi = Math.max(xi, 0);
				yi = Math.min(yi, t.ny-1);
				yi = Math.max(yi, 0);
				nimps[bins[xi][yi]]++;
				impbin[i] = bins[xi][yi]; 
			}
			
			System.out.println("Bin\tN_imps");
			for (int i = 0; i < nbins; i++)
			{
				System.out.println("" + i + "\t" + nimps[i] + "");
			}
			
			System.out.println("imp\tbin\tx\ty");
			for (int i = 0; i < imps.size(); i++)
				System.out.println("" + i + "\t" + impbin[i] + "\t" + imps.get(i).getPosition()[0] + "\t" + imps.get(i).getPosition()[1]);;
			
			
		}
		if (arg0.getKeyChar()=='O')
		{
			FieldOps.copy(Layer.openFree(fc).data, t.data);
			resetGraphics(true);
		}
		if (arg0.getKeyChar() == 'c' && showFFT) //make an image of the complex fft
		{
			if (fftz == null)
				fftz = new double [t.nx][t.ny][2];
			
			FFTOps.putFFT(t.data, fftz, false);
			boolean save = JOptionPane.showConfirmDialog(null, "Save real and imaginary parts?") == JOptionPane.YES_OPTION;
			if (save)
			{
				double[][] real = FieldOps.getIndex(fftz, 0);
				double[][] real2 = new double [t.nx][t.ny];
				FieldOps.shift(real, real2);
				double[][] imag = FieldOps.getIndex(fftz, 1);
				double[][] imag2 = new double [t.nx][t.ny];
				FieldOps.shift(imag, imag2);
				String s = FileOps.selectSave(fc).toString();
				Layer.writeBIN(Layer.newLayer(t, real2), s + "_re.bin");
				Layer.writeBIN(Layer.newLayer(t, imag2), s + "_im.bin");
			}
			ImageEditing.enlargeBasic(FFTOps.getImageCent(fftz, ftlog, false, FieldOps.min(fftmag)), image, sizeratio);
			refresh = true;
			repaint();
		}
		if (arg0.getKeyChar() == 'C') //copy the current image to clipboard
		{
//			ImageEditing.copyToClipboard(image);
			BufferedImage export = new BufferedImage(WIDTH, HEIGHT, BufferedImage.TYPE_INT_RGB);
			Graphics g = export.getGraphics();
			refresh = true;
			paint(g);
			try
			{
				ImageEditing.copyToClipboard(ImageEditing.getSubset(export, ox, sx*sizeratio, oy, sy*sizeratio));
			}
			catch(Exception e){
				ImageEditing.copyToClipboard(image);
			}
		}
		if (arg0.getKeyChar() == 'W')
		{
			String imput2 = JOptionPane.showInputDialog("Enter minimum and maximum of the color scale.", "" + scale.getMin() + "," + scale.getMax());
			double min = Double.parseDouble(imput2.split(",")[0].trim());
			double max = Double.parseDouble(imput2.split(",")[1].trim());
			BufferedImage output = ImageEditing.getBufferedImage(showFFT ? fftmag : t.data, ColorScales.getNew(min, max, currentCScale));
			ImageEditing.copyToClipboard(output);
			SRAW.writeImage(FileOps.selectSave(fc).toString(), output);
		}
		if (arg0.getKeyChar() == 'i')
		{
			PointImp[] imps = PointImp.getAboveCutoffContiguous(t, scale.getMax());
			PointImp.writeToFile(imps, FileOps.selectSave(fc));
		}
		if (arg0.getKeyChar() == 'I')
		{
			if (imps == null) {imps = ImpurityListEditor.getImpurities(fc, 1); refresh = false; repaint();}
			else{
				imps = null;
				refresh = true;
				repaint();
			}
		}
		if (arg0.getKeyChar() == 'J')
		{
			if (imps == null) {
				imps = ImpurityListEditor.getImpurities(fc, 1);
				imps2 = ImpurityListEditor.getImpurities(fc, 1);
				refresh = false; repaint();
				}
			else{
				imps = null;
				imps2 = null;
				refresh = true;
				repaint();
			}
		}
		if (arg0.getKeyChar() == 'V')
		{
			ImageEditing.copyToClipboard(((ColorScale1D)scale).getScaleImage(8));
		}

		if (sectionDrawer != null) sectionDrawer.processKeyStroke(arg0.getKeyChar());
		if (specView != null) specView.processKeyStroke(arg0.getKeyChar());
		
		if (arg0.getKeyChar() == 'Y')
		{
			AtomicCoordinatesSet latt = new AtomicCoordinatesSet(ColumnIO.getASCII(FileOps.selectOpen(fc).toString()));
			t = TopomapUtil.makeSymmetrizedFFT(t, latt, JOptionPane.showConfirmDialog(null, "Lattice is square?") == JOptionPane.YES_OPTION);
			Layer.writeBIN(t,fc);
			this.resetColorScale();
			this.formImage();
			this.repaint();
		}
		if (arg0.getKeyChar() == 'U')
		{
			AtomicCoordinatesSet latt = new AtomicCoordinatesSet(ColumnIO.getASCII(FileOps.selectOpen(fc).toString()));
//			t = TopomapUtil.makeSymmetrizedFFT(t, latt, JOptionPane.showConfirmDialog(null, "Lattice is square?") == JOptionPane.YES_OPTION);
			int index = Integer.parseInt(JOptionPane.showInputDialog("Enter your choice for the second derivative:\r\n"+
					"0 - faa\r\n"+
					"1 - fab\r\n"+
					"2 - fba\r\n"+
					"3 - fbb\r\n"+
					"4 - 1/2(fab+fba)\r\n"+
					"5 - 1/2(fab-fba) which is zero identically\r\n"+
					"6 - 1/2(faa+fbb)\r\n"+
					"7 - 1/2(faa-fbb)\r\n"));
			
			backup = Layer.newLayer(t, FieldOps.copy(t.data));
			
			AtomicCoordinatesSet recip = latt.getReciprocalLattice();
			double[][] units = FieldOps.getGrahamSchmidtUnitVectors(recip.getA(), recip.getB());
					t.data = FieldOps.getDirectional2ndDerivativeStuff(t.data, units[0], units[1])[index];
			this.resetColorScale();
			this.formImage();
			this.repaint();
		}

	}
	public void mouseWheelMoved(MouseWheelEvent arg0) {
	}
	public void mouseDragged(MouseEvent arg0) {
		int x = arg0.getX();
		int y = arg0.getY();
		System.out.println(arg0.getModifiers());
		if (sectionDrawer != null)
		{
			if (arg0.getModifiers() == 17)//MouseEvent.SHIFT_DOWN_MASK)
				moveLineBoth((x - currentx), (y - currenty));
			else
				moveLine((x - currentx), (y - currenty));
			sectionDrawer.putPoints(line);
			sectionDrawer.refresh();
			refresh = true;
			repaint();
		}
		else{
            Point currCoords = arg0.getLocationOnScreen();
            this.setLocation(currCoords.x - mouseDownCompCoords.x, currCoords.y - mouseDownCompCoords.y);
		}
		//fixed the mouse drag by updating the current position of the mouse every MouseEvent
		currentx=x;
		currenty=y;
	}
	
	public void mouseMoved(MouseEvent arg0) {
		currentx = arg0.getX();
		currenty = arg0.getY();
		calcx = currentx - ox;
		calcy = currenty - oy;
		currentid = calcx/(double)imageOverField;
		currentjd = calcy/(double)imageOverField;
		currenti = (int)(calcx/imageOverField);
		currentj = (int)(calcy/imageOverField);
		if (currenti < 0) currenti = 0;
		if (currenti >= t.nx) currenti = t.nx-1;
		if (currentj < 0) currentj = 0;
		if (currentj >= t.ny) currentj = t.ny-1;
		
		if (sectionDrawer != null)
		{
			d1 = Distance.distance(currenti - line[0][0], currentj - line[0][1]);
			d2 = Distance.distance(currenti - line[1][0], currentj - line[1][1]);
			if (d1 < d2) nearestEndIndex = 0;
			else nearestEndIndex = 1;
			
			sOfCursor = ((currenti-line[0][0])*(line[1][0]-line[0][0]) + (currentj-line[0][1])*(line[1][1]-line[0][1]))/Distance.distance(line[0], line[1]);
			//also calculate the horitontal distance along the line:
			sectionDrawer.setSVert(sOfCursor);
//			System.out.println(sOfCursor);
		}
		if (showFFT && t.lineTime != 0)
		{
			mouseFreq = (currenti-t.nx/2)/t.lineTime + (currentj - t.ny/2)/(t.ny*t.lineTime);
			if (mouseFreq < 0) mouseFreq += t.nx/t.lineTime;
			setTitle(""+mouseFreq);
		}
		refresh = false;
		repaint();
//		refresh = false;
//		repaint();
	}
	public void moveLine(int dx, int dy)
	{
		line[nearestEndIndex][0] += ((double)dx)/((double)sizeratio);
		line[nearestEndIndex][1] += ((double)dy)/((double)sizeratio);
		System.out.println(Printer.getTable(line) + "a\r\n");
	}
	public void moveLineBoth(int dx, int dy)
	{
		line[0][0] += dx/imageOverField;
		line[0][1] += dy/imageOverField;
		line[1][0] += dx/imageOverField;
		line[1][1] += dy/imageOverField;
		System.out.println(Printer.getTable(line) + "b\r\n");
	}
//	private void showCustomLayers() {
//		// TODO Auto-generated method stub
//		int sizeUnit = t.nx;
//		int miniSizeUnit = 32;
//		int b1x = 216, b1y = -218;
//		int b2x =  218, b2y = 216;
//		int blowUpFactor = (2*sizeUnit)/miniSizeUnit;
//		
//		double[][] bragg1Area = FieldOps.subset(fftmag, b1x - miniSizeUnit/2 + t.nx/2, b1x + miniSizeUnit/2 + t.nx/2, b1y - miniSizeUnit/2 + t.ny/2, b1y + miniSizeUnit/2 + t.ny/2);
//		double[][] bragg2Area = FieldOps.subset(fftmag, b2x - miniSizeUnit/2 + t.nx/2, b2x + miniSizeUnit/2 + t.nx/2, b2y - miniSizeUnit/2 + t.ny/2, b2y + miniSizeUnit/2 + t.ny/2);
//		
//		double[] x1 = ArrayOps.generateArrayNotInclUpper(b1x, b1x+miniSizeUnit, miniSizeUnit);
//		double[] x2 = ArrayOps.generateArrayNotInclUpper(b2x, b2x+miniSizeUnit, miniSizeUnit);
//		double[] y1 = ArrayOps.generateArrayNotInclUpper(b1y, b1y+miniSizeUnit, miniSizeUnit);
//		double[] y2 = ArrayOps.generateArrayNotInclUpper(b2y, b2y+miniSizeUnit, miniSizeUnit);
//		
//		LayerViewer.show(new Layer(bragg1Area, x1, y1, t.v, t.current), 1024, false);
//		LayerViewer.show(new Layer(bragg2Area, x2, y2, t.v, t.current), 1024, false);
//	}

	public double currenti(double currentx)
	{
		return (currentx-ox)/sizeratio;
	}
	public double currentj(double currenty)
	{
		return (currenty-ox)/sizeratio;
	}
	public int screenX(double currenti)
	{
		return FieldOps.round(currenti*sizeratio)+ox;
	}
	public int screenY(double currentj)
	{
		return FieldOps.round(currentj*sizeratio)+oy;
	}
	public void mouseClicked(MouseEvent arg0) {
	}
	public void mouseEntered(MouseEvent arg0) {
	}
	public void mouseExited(MouseEvent arg0) {
	}
	public void mousePressed(MouseEvent arg0) {
		mouseDownCompCoords = arg0.getPoint();
	}
	public void mouseReleased(MouseEvent arg0) {
		mouseDownCompCoords = null;
	}
	
	public static class SliderPanel extends JPanel implements ChangeListener
	{
		private static final long serialVersionUID = 3021416954431474662L;
		LayerViewer parent;
		JFrame frame;
		public JSlider min, max;
		int parentmember;
		int oldvalue = 0;
		int oldminv = 0, oldmaxv = 999;
		static int npts = 1001;
		
		public SliderPanel(LayerViewer parent, JFrame frame)
		{
			this.frame = frame;
			this.parent = parent;
			this.setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
			this.setBorder(new LineBorder(Color.GRAY));
			min = new JSlider(0, 1000, oldminv);
			max = new JSlider(0, 1000, oldmaxv);
			min.addChangeListener(this);
			max.addChangeListener(this);
			add(min);
			add(max);
			frame.setSize(3+40, 80);
			frame.add(this);
		}
		public void show(){
			frame.setVisible(true);
		}
		@Override
		public void stateChanged(ChangeEvent arg0) {
//			if (s.getValue() == oldvalue && min.getValue() == oldminv && max.getValue() == oldmaxv) return;
			if (min.getValue() != oldminv || max.getValue() != oldmaxv)
			{
				parent.resetColorScale(((double)min.getValue()/1000), ((double)max.getValue()/1000), false);
				parent.formImage();
				parent.refresh = true;
				parent.repaint();
				oldminv = min.getValue();
				oldmaxv = max.getValue();
				frame.setTitle("Min: " + parent.scalemin +", Max: " +parent.scalemax);
				return;
			}
//			else if (s.getValue() == oldvalue) return;
//			parent.resetColorScale(((double)min.getValue()/1000), ((double)max.getValue()/1000));
		}
	}
	
	public static boolean withinBounds(int i, int j, double[][] array)
	{
		return i > 0 && j > 0 && i < array.length && j < array[i].length;
	}
	public static void main(String[] args)
	{
		Topomap.setStdDir();
		fc = new JFileChooser(Topomap.stddir);
		Layer t = Layer.openFree(fc);
		//Changed to make the window display as large as possible while fitting in your display
		//If you want to make it smaller, increase the value of makeWindowSmaller from 0
		java.awt.Dimension screenSize = Toolkit.getDefaultToolkit().getScreenSize();
		int makeWindowSmaller=0;
		int largerDim= Math.max(t.nx,t.ny);
		float scalefactor = (float)screenSize.getHeight()/largerDim - makeWindowSmaller;
		int sizechoice = (int)(largerDim*scalefactor);
		new LayerViewer(t, Topomap.stddir, sizechoice);
	}
	
	public static void show(Layer t, int size, boolean closeDefault)
	{
		Topomap.setStdDir();
		fc = new JFileChooser(Topomap.stddir);
		LayerViewer lv = new LayerViewer(t, Topomap.stddir, size);
		if (!closeDefault)
			lv.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		else 
			lv.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}
}
