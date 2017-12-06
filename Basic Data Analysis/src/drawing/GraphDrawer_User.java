package drawing;

import java.awt.event.KeyEvent;
import java.awt.event.MouseEvent;

public abstract class GraphDrawer_User {

	protected GraphDrawerCart g;
	
	public abstract void processKeyStroke(KeyEvent e);
	public abstract void processMouseClick(MouseEvent e);
	
}
