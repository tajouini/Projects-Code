
/**
 * Write a description of CharactersInPlay here.
 * 
 * @author (your name) 
 * @version (a version number or a date)
 */
import java.util.ArrayList; // import the ArrayList class
import edu.duke.FileResource;

public class CharactersInPlay {
    
    private ArrayList<String> names;
    private ArrayList<Integer> counts;
    
    public CharactersInPlay(){
        names = new ArrayList<String>();
        counts = new ArrayList<Integer>();
        
    }
    
    
    public void update(String person){
        if (names.indexOf(person)==-1){
            names.add(person);
            counts.add(1);
        }
        else {
            int val = counts.get(names.indexOf(person));
            counts.set(names.indexOf(person), val+1);
        }
    }
    
    public void findAllCharacters(){
        names.clear();
        counts.clear();
        FileResource file = new FileResource();
        for (String k: file.lines()){
            if (k.indexOf(".")>=0){
                String person=k.substring(0,k.indexOf("."));
                update(person);
            }
        
        }
    }
    
    public void charactersWithNumParts(int num1, int num2){
        for (int k=0; k<names.size(); k++){

            if (num1<=counts.get(k) && counts.get(k)<=num2){
                System.out.println(names.get(k)+ " with counts "+ counts.get(k));
            }
        }
    }
    
    public void tester(){
        findAllCharacters();
        int k_max=0;
        for (int k=0; k<names.size();k++){
            if (counts.get(k)>counts.get(k_max)){
                k_max = k;
            }
        }
        System.out.println(names.get(k_max)+" with counts "+counts.get(k_max));
        int num1 = 10;
        int num2 = 15;
        System.out.println("******Names with lines between****"+num1+ " and "+ num2);
        charactersWithNumParts(num1, num2);
        
    }

}
