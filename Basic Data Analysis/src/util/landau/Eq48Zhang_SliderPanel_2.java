package util.landau;

import java.awt.Color;

import javax.swing.BoxLayout;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.border.LineBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import util.ArrayOps;
import util.Printer;
import util.SpectraUtil.LandauLevelPeak;

import drawing.GraphDrawerCart;

public class Eq48Zhang_SliderPanel_2 extends JPanel implements ChangeListener{

	GraphDrawerCart g;
	LandauLevelModel m;
	
	LandauLevelPeak[][] peaks;
	
	int[] npts;
	double[][] values;
	
	double[] min, max;
	//These are model-specific and should probably not be inlcuded in a "general" constructor, but we do it anyway.
	
	JFrame frame;
	public JSlider[] para;
	int np = 4;
	
	double[][] x, y;
	
	public Eq48Zhang_SliderPanel_2(LandauLevelPeak[][] peaks, LandauLevelModel m)
	{
		npts = new int [np];
		for (int i = 0; i < np; i++)
			npts[i] = 1001;
		
		min = new double [np];
		max = new double [np];
		double edp_temp = 0;
		int countZero = 0;
		for (int i = 0; i < peaks.length; i++)
			for (int j = 0; j < peaks[i].length; j++)
				if (peaks[i][j].n == 0)
				{
					edp_temp += peaks[i][j].max;
					countZero++;
				}
		
		if (countZero != 0){
			edp_temp /= countZero;
			min = new double[] {edp_temp-0.1, -50, -300, 0.1};
			max = new double[] {edp_temp+0.1, 150, 300, 10};
		}
		else{
			min = new double[] {-0.3, -50, -300, 0.1};
			max = new double[] {0.0, 150, 300, 10};
			
		}
		min = new double[] {0, 30, 0, 0};
		max = new double[] {0, 30, 0, 4};
		values = new double[np][];
		for (int i = 0; i < np; i++)
		{
			values[i] = ArrayOps.generateArrayInclBoth(min[i], max[i], npts[i]+1);
		}
		this.frame = new JFrame();
		this.m = m;
		this.peaks = peaks;
		double[][][] xy = LandauLevelPeak.getOrderedPairsModel(peaks, m);
		double[][][] xy2 = LandauLevelPeak.getOrderedPairs(peaks);
		x = new double [xy[0].length+xy2[0].length][];
		y = new double [xy[1].length+xy2[1].length][];
		for (int i = 0; i < xy[0].length; i++)
		{
			x[i] = new double [xy[0][i].length];
			x[i+xy[0].length] = new double [xy[0][i].length];
			y[i] = new double [xy[0][i].length];
			y[i+xy[0].length] = new double [xy[0][i].length];
			for (int j = 0; j < xy[0][i].length; j++)
			{
				x[i][j] = xy[0][i][j];
				y[i][j] = xy[1][i][j];
				x[i+xy[0].length][j] = xy2[0][i][j];
				y[i+xy[0].length][j] = xy2[1][i][j];
			}
		}
		g = GraphDrawerCart.getNew(x, y);
		g.showWindow();
		this.setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
		this.setBorder(new LineBorder(Color.GRAY));
		para= new JSlider[np];
		for (int i = 0; i < np; i++){
			para[i] = new JSlider(0, npts[i], 0);
			para[i].setValue(0);
			para[i].setSnapToTicks(true);
			para[i].setMinorTickSpacing(1);
			para[i].addChangeListener(this);
			add(para[i]);
		}
		
		
		frame.setSize(npts[0]+40, 80);
		frame.add(this);
		frame.setVisible(true);
	}
	
	public void resetModel(double[] param){
		m.setParam(param);
		
		double[][][] xy = LandauLevelPeak.getOrderedPairsModel(peaks, m);
		for (int i = 0; i < xy[0].length; i++)
			for (int j = 0; j < xy[0][i].length; j++){
				x[i][j] = xy[0][i][j];
				y[i][j] = xy[1][i][j];
			}
		g.repaint();
	}

	
	@Override
	public void stateChanged(ChangeEvent arg0) {
		// TODO Auto-generated method stub
		double[] param = new double [np];
		for (int i = 0; i < np; i++)
			param[i] = values[i][para[i].getValue()];
		resetModel(param);
		frame.setTitle(Printer.vectorP(param));
	}
	
	public static void main (String[] args){
		int[] N = {0, 1, 2, 3, 4, 5, 6, 7, 8};
		double[] B = {1, 2, 3, 4, 5, 7, 8};
		LandauLevelPeak[][] peaks = new LandauLevelPeak[N.length][B.length];
		for (int i = 0; i < N.length; i++)
			for (int j = 0; j < B.length; j++)
				peaks[i][j] = new LandauLevelPeak(0, N[i], B[j]);
		new Eq48Zhang_SliderPanel_2(peaks, new Eq48Zhang());
	}

	
	
}
