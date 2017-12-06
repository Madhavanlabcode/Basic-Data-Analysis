package drawing;

import main.FieldMain;

//We assume all input data is 2^n X 2^n square table, and n >= 9.
public class DataZoomPanel {

	int scalex, scaley;
	int recipx  = 0, recipy = 0;
	int sx, sy;
	int px = 512, py = 512;
	
	int zoom;
	//These are the relative coordinates. They run from [0, scalex/zoom).
	int fx, fy;
	
	FieldDrawer drawer;
	FieldMain parent;
	
	public double[][] data, field;
	
	public DataZoomPanel(FieldMain parent) 
	{
		this.parent = parent; this.drawer = parent.drawer;
		data = parent.data;
		
		//Smoothes the data down into 512x512
		sx = data.length;
		sy = data[0].length;
		scalex = sx/px; scaley = sy/py;
		if(scalex == 0 || scaley == 0)
		{
//			px = sx; py = sy; scalex = 1; scaley = 1;
			recipx = px/sx; recipy = py/sy;
		}
		else;
			
		field = new double[px][py];
		//This generates the field as the average over the whole data
		int i = 0, j = 0;
		for (i = 0; i < px; i++)
			for(j = 0; j < py; j++)
			{
				if (scalex != 0){
					field[i][j] = 0;
					for (int k = 0; k < scalex; k++)
						for (int m = 0; m < scaley; m++)
						{
							field[i][j] += data[i*scalex + k][j*scaley + m];
						}
					field[i][j] /= scalex*scaley;
					}
				else
				{
					field[i][j] = data[i/recipx][j/recipy];
				}
			}
		
		zoom = scalex;
		fx = 0; fy = 0;
	}
	
	public void zoomInOnce(int x, int y)
	{
		if (zoom==0) return;
		//x, y can be 0,1 0,1.
		zoom/=2;
		fx = fx*2 + x;
		fy = fy*2 + y;
		
		int i, j;
		for (i = 0; i < px; i++)
			for(j = 0; j < py; j++)
			{
				field[i][j] = 0;
				for (int k = 0; k < zoom; k++)
					for (int m = 0; m < zoom; m++)
					{
						field[i][j] += data[fx*px*zoom + i*zoom + k][fy*py*zoom + j*zoom + m];
					}
				field[i][j] /= zoom*zoom;
			}
		}
	
}
