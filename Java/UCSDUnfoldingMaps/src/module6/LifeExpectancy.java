package module6;

import processing.core.PApplet;
import de.fhpotsdam.unfolding.UnfoldingMap;
import de.fhpotsdam.unfolding.utils.MapUtils;
import parsing.ParseFeed;
import de.fhpotsdam.unfolding.providers.*;
import de.fhpotsdam.unfolding.providers.Google.*;

import java.util.List;
import de.fhpotsdam.unfolding.data.Feature;
import de.fhpotsdam.unfolding.data.GeoJSONReader;

import java.util.HashMap;


import de.fhpotsdam.unfolding.marker.Marker;

/**
 * Visualizes life expectancy in different countries. 
 * 
 * It loads the country shapes from a GeoJSON file via a data reader, and loads the population density values from
 * another CSV file (provided by the World Bank). The data value is encoded to transparency via a simplistic linear
 * mapping.
 */
public class LifeExpectancy extends PApplet {

	UnfoldingMap map;
	HashMap<String, Float> lifeExpMap;
	List<Feature> countries;
	List<Marker> countryMarkers;
	int color_var;

	public void setup() {
		size(800, 600, OPENGL);
		map = new UnfoldingMap(this, 50, 50, 700, 500, new Google.GoogleMapProvider());
		MapUtils.createDefaultEventDispatcher(this, map);

		// Load lifeExpectancy data
		lifeExpMap = ParseFeed.loadLifeExpectancyFromCSV(this,"LifeExpectancyWorldBank.csv");
		

		// Load country polygons and adds them as markers
		countries = GeoJSONReader.loadData(this, "countries.geo.json");
		countryMarkers = MapUtils.createSimpleMarkers(countries);
		map.addMarkers(countryMarkers);
		System.out.println(countryMarkers.get(0).getId());
		
		// Country markers are shaded according to life expectancy (only once)
		shadeCountries();
	}

	public void draw() {
		// Draw map tiles and country markers
		map.draw();
		addKey();
	}

	//Helper method to color each country based on life expectancy
	//Red-orange indicates low (near 40)
	//Blue indicates high (near 100)
	private void shadeCountries() {
		for (Marker marker : countryMarkers) {
			// Find data for country of the current marker
			String countryId = marker.getId();
			System.out.println(lifeExpMap.containsKey(countryId));
			if (lifeExpMap.containsKey(countryId)) {
				float lifeExp = lifeExpMap.get(countryId);
				// Encode value as brightness (values range: 40-90)
				int colorLevel = (int) map(lifeExp, 40, 90, 10, 255);
				 color_var = color(255-colorLevel, 100, colorLevel);
				marker.setColor(color_var);
			}
			else {
				marker.setColor(color(150,150,150));
			}
		}
	}
	
	private void addKey() {	
		// Remember you can use Processing's graphics methods here
		
		fill(255, 255, 240);
		
		int xbase = 25;
		int ybase = 50;
		
		rect(xbase, ybase, 150, 250);
		
		fill(0);
		textAlign(LEFT, CENTER);
		textSize(12);
		text("Life Expectancy Key", xbase+25, ybase+25);
		
		fill(color(0,102,204));
		int tri_xbase = xbase + 35-5;
		int tri_ybase = ybase + 50-5;
		rect(tri_xbase, tri_ybase,10,10);

		fill(0);
		textAlign(LEFT, CENTER);
		text(">80", tri_xbase + 20, tri_ybase);
		
		text("<80", xbase+50, ybase+70);
		text("<60", xbase+50, ybase+90);
		text("<50", xbase+50, ybase+110);		
		fill(98,83,137);
		rect(xbase+35-5, 
				ybase+70-5, 
				10, 
				10);
		fill(color_var);
		rect(xbase+35-5, ybase+90-5, 10, 10);
		fill(color(255,102,0));
		rect(xbase+35-5, ybase+110-5, 10, 10);
		
		
		
	}



}
