package drawing;

import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.*;
import javax.swing.border.LineBorder;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

public class FieldPanel extends JFrame implements ActionListener{

	private static final long serialVersionUID = -8888590328274505702L;
	FieldDrawer drawer;
	FieldDrawer2 drawer2;
	
	BoundsPanel bounds;
	AnglePanel angle;
	public FieldPanel(FieldDrawer drawer)
	{
		this.drawer = drawer;
		JPanel main = new JPanel();
		setTitle("Controls");
		main.setLayout(new BoxLayout(main, BoxLayout.PAGE_AXIS));
		this.bounds = new BoundsPanel(this);
		main.add(bounds);
		this.angle = new AnglePanel(this);
		main.add(angle);
		add(main);
		setSize(200, 100);
		setVisible(true);

	}
	public FieldPanel(FieldDrawer2 drawer2)
	{
		this.drawer2 = drawer2;
		JPanel main = new JPanel();
		setTitle("Controls");
		main.setLayout(new BoxLayout(main, BoxLayout.PAGE_AXIS));
		this.bounds = new BoundsPanel(this);
		main.add(bounds);
		this.angle = new AnglePanel(this);
		main.add(angle);
		add(main);
		setSize(200, 100);
		setVisible(true);

	}
	
	@Override
	public void actionPerformed(ActionEvent arg0) {
		// TODO Auto-generated method stub
		
	}
	public static class BoundsPanel extends JPanel implements ChangeListener
	{
		private static final long serialVersionUID = 369430464950527803L;
		FieldPanel parent;
		public JSlider min, max;
		public BoundsPanel(FieldPanel parent)
		{
			this.parent = parent;
			this.setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
			this.setBorder(new LineBorder(Color.GRAY));
			min = new JSlider(0, 1000, 0);
			min.addChangeListener(this);
			max = new JSlider(0, 1000, 1000);
			max.addChangeListener(this);

			add(min);
			add(max);
		}
		@Override
		public void stateChanged(ChangeEvent arg0) {
			parent.drawer2.redrawColorScale(((double)min.getValue()/1000), ((double)max.getValue()/1000));
		}
	}
	
	public static class AnglePanel extends JPanel implements ChangeListener
	{
		private static final long serialVersionUID = 269926310445678633L;
		FieldPanel parent;
		public JSlider theta;
		public AnglePanel(FieldPanel parent)
		{
			this.parent = parent;
			this.setLayout(new BoxLayout(this, BoxLayout.Y_AXIS));
			this.setBorder(new LineBorder(Color.GRAY));
			theta = new JSlider(0, 1000, 0);
			theta.addChangeListener(this);
			
			add(theta);
		}
		@Override
		public void stateChanged(ChangeEvent arg0) {
			if (parent.drawer2.gradcalc.active)
			{
				parent.drawer2.gradcalc.resetTheta(theta.getValue()*Math.PI/500);
			}
		}
	}
}
