
/**
 * Write a description of WordFrequencies here.
 * 
 * @author (your name) 
 * @version (a version number or a date)
 */
import edu.duke.FileResource;
import java.util.ArrayList; // import the ArrayList class
public class WordFrequencies {
    
    private ArrayList<String> myWords; 
    private ArrayList<Integer> myFreqs; 
    private int maxId = 0;
    public WordFrequencies(){
        myWords = new ArrayList<String>();
        myFreqs = new ArrayList<Integer>();

    }
    public void findUnique(){
        myWords.clear();
        myFreqs.clear();
        FileResource file = new FileResource();
        for (String word: file.words()){
            word = word.toLowerCase();
            if (myWords.indexOf(word)==-1){
                myWords.add(word);
                myFreqs.add(1);
            }
            else{
                int val = myFreqs.get(myWords.indexOf(word));
                myFreqs.set(myWords.indexOf(word), val+1);
                
            }
        }
    }
    public void findIndexMax(){
   
        for (int j=0; j< myFreqs.size(); j++){
            if (myFreqs.get(maxId)<myFreqs.get(j)){
                maxId  = j;
            }
        }
    }

    public void tester(){
        findUnique();
        System.out.println("Number of unique words is "+myWords.size());
        //for (int j=0; j<myWords.size(); j++){
        //    System.out.println(myWords.get(j)+ " appears "+myFreqs.get(j)+
        //    " times");
        //}
        findIndexMax();
         System.out.println("The word "+myWords.get(maxId)+" occurs most often and its count are: "
         + myFreqs.get(maxId));
    }
}
