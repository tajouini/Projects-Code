package guimodule;

import processing.core.PApplet;

public class MyDisplay extends PApplet {
	
	public void setup() {
		size(400, 400);
		background(200,200,200);
		
	}
	
	public void draw() {
		fill(255,255,0); // yellow
		ellipse(200,200,390,390);
		fill(0,0,0);	// yellow
		ellipse(120,130,50,50);
		ellipse(280,130,50,50);
		noFill();
		arc(200,280,90,90,0,PI);


	}
	

}
