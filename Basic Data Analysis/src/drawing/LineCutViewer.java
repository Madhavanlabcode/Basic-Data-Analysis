package drawing;

import image.ImageEditing;

import java.awt.Color;
import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;
import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.border.LineBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import schrodinger.MovieMaker;
import util.ArrayOps;
import util.fileops.FileOps;
import util.fileops.LineOfSpectra;
import util.fileops.Topomap;
import util.geom.Distance;

/**
 * This class allows the viewing of single, defined Line cut. There will be two graph windows: one for the z-line, and one for the spectrum. The spectrum will be selected 
 * by adjusting a slider.
 * @author madhavanlab2011
 *
 */
public class LineCutViewer extends GraphDrawer_User{

	LineOfSpectra data;
	double[] s; //distance along the graph, from 0 to the end.
	double[] zBounds;
	double[] twoValuesS;
	int selected;
	GraphDrawerCart zGraph;
	GraphDrawerCart specGraph;
	
	SliderPanel sp;
	
	public LineCutViewer(LineOfSpectra data)
	{
		this.data = data;
		s = new double [data.nspec];
		for (int i = 0; i < s.length; i++)
			s[i] = Distance.distance(new double[] {data.x[0], data.y[0]}, new double[] {data.x[i], data.y[i]});
		
		zBounds = new double [] {ArrayOps.min(data.z), ArrayOps.max(data.z)};
		selected = 0;
		twoValuesS = new double [] {0, 0};
		zGraph = new GraphDrawerCart("Topography", s, data.z);
		zGraph.setXY(new GraphDrawerCart.GraphObject(twoValuesS, zBounds), false, false, 1);
		zGraph.setRange(zBounds[0], zBounds[1]);
		zGraph.showWindow();
		zGraph.user = this;

		specGraph = new GraphDrawerCart("Spectra", data.v, data.getSpectrum(0));
		specGraph.setNumPlots(data.nspec);
		for (int i = 0; i < data.nspec; i++)
		{
			specGraph.setXY(new GraphDrawerCart.GraphObject(data.v, data.getSpectrum(i)), false, false, i);
			specGraph.plot[i].isDrawing = i == 0;
			specGraph.plot[i].c = Color.RED;
		}
		specGraph.setRange(ArrayOps.min(data.data), ArrayOps.max(data.data));
		specGraph.showWindow();
		
		sp = new SliderPanel(this, new JFrame());
		sp.show();

	}
	
	public void setSelected(int value) {
		
		selected = value;
		twoValuesS[0] = s[selected];
		twoValuesS[1] = s[selected];
		zGraph.repaint();
		for (int i = 0; i < data.nspec; i++)
			specGraph.plot[i].isDrawing = i == selected;
		specGraph.repaint();
	}
	public static class SliderPanel extends JPanel implements ChangeListener
	{
		LineCutViewer parent;
		JFrame frame;
		public JSlider s;
		int npts;
		
		public SliderPanel(LineCutViewer parent, JFrame frame)
		{
			npts = parent.data.nspec;
			this.frame = frame;
			this.parent = parent;
			this.setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
			this.setBorder(new LineBorder(Color.GRAY));
			s = new JSlider(0, npts, 0);
			s.setValue(0);
			
			s.setSnapToTicks(true);
			s.setMinorTickSpacing(1);
			s.addChangeListener(this);
			
			add(s);
			frame.setSize(3*npts+40, 80);
			frame.add(this);
		}
		public void show(){
			frame.show();
		}
		@Override
		public void stateChanged(ChangeEvent arg0) {
//			if (s.getValue() == oldvalue && min.getValue() == oldminv && max.getValue() == oldmaxv) return;
//			else if (s.getValue() == oldvalue) return;
//			parent.resetColorScale(((double)min.getValue()/1000), ((double)max.getValue()/1000));
			parent.setSelected(s.getValue());
		}

	}
	@Override
	public void processKeyStroke(KeyEvent e) {
		// TODO Auto-generated method stub
		switch(e.getKeyChar())
		{
		case 'm': //make a movie with the two tiles stacked vertically. The total dimensions will be width, height in ratio 4:3. 
			int width = 800, height = 3*width/4;
			String basename = FileOps.selectSave(null).toString();
			BufferedImage upper, lower, total = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
			for (int i = 0; i < data.nspec; i++)
			{
				upper = GraphDrawerCart.drawPlots(new double[][] {s, new double[] {s[i], s[i]}}, new double[][] {data.z, zBounds}, width, height/2);
				lower = GraphDrawerCart.drawPlot(data.v, data.data[i], ArrayOps.min(data.v), ArrayOps.max(data.v), ArrayOps.min(data.data), ArrayOps.max(data.data), width, height/2);
				ImageEditing.copyInto(upper, total, 0, 0);
				ImageEditing.copyInto(lower, total, 0, height/2);
				try {
					ImageIO.write(total, "png", new File(basename + MovieMaker.fromInt(i) + ".png"));
				} catch (IOException e1) {
					// TODO Auto-generated catch block
					e1.printStackTrace();
				}
//				SRAW.writeImage(basename + MovieMaker.fromInt(i), total);
			}
			
			
		}
	}

	@Override
	public void processMouseClick(MouseEvent e) {
		// TODO Auto-generated method stub
		
	}
	
	public static void main(String[] args)
	{
		Topomap.setStdDir();
		File f = FileOps.selectOpen(null);
		new LineCutViewer(LineOfSpectra.readBIN(f.toString()));
	}

}
