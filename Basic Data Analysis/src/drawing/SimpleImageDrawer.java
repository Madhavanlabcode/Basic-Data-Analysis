package drawing;

import java.awt.Graphics;
import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.image.BufferedImage;

import javax.swing.JFileChooser;
import javax.swing.JFrame;

import main.SRAW;
import util.fileops.FileOps;
import util.fileops.Topomap;

/**Draws an image an lets the user save it.
 * 
 * @author madhavanlab2011
 *
 */
public class SimpleImageDrawer extends JFrame implements KeyListener, MouseListener
{
	public BufferedImage image;
	JFileChooser fc;
	
	public SimpleImageDrawer(BufferedImage image)
	{
		fc = new JFileChooser(Topomap.stddir);
		this.image = image;
		this.setSize(image.getWidth() + 20, image.getHeight() + 60);
		this.show();
		this.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
		this.repaint();
		addKeyListener(this);
	}
	
	public void paint(Graphics g)
	{
		g.clearRect(0, 0, this.getWidth(), this.getHeight());
		g.drawImage(image, 10, 40, null);
	}

	@Override
	public void keyPressed(KeyEvent arg0) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void keyReleased(KeyEvent arg0) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void keyTyped(KeyEvent arg0) {
		// TODO Auto-generated method stub
		switch (arg0.getKeyChar())
		{
		case 'S':
			SRAW.writeImage(FileOps.selectSave(fc).toString(), image);
			break;
		}
	}

	@Override
	public void mouseClicked(MouseEvent e) {
	}
	@Override
	public void mouseEntered(MouseEvent e) {
	}
	@Override
	public void mouseExited(MouseEvent e) {
	}
	@Override
	public void mousePressed(MouseEvent e) {
	}
	@Override
	public void mouseReleased(MouseEvent e) {
	}
}
