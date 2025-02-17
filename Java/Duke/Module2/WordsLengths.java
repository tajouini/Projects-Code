
/**
 * Write a description of WordsLengths here.
 * 
 * @author (your name) 
 * @version (a version number or a date)
 */
import edu.duke.FileResource;

public class WordsLengths {
    
    public void countWordLengths(FileResource resource, int[] counts){
        for (String s: resource.words()){
            if (!Character.isLetter(s.charAt(0))||!Character.isLetter(s.charAt(s.length()-1))){
                counts[s.length()-1]+=1;
            }
            else{
            if (s.length()<counts.length){
                counts[s.length()]+=1;

            }
            else{
                counts[counts.length-1]+=1;
            }
            }
        }
        
    }
    
    public int indexOfMax(int[] vals){
        int ind_max = 0;
        for (int k=0; k<vals.length; k++){
            if (vals[k]>vals[ind_max]){
             ind_max =k;   
            }  
        }
        return ind_max;
    }
    
    
    public  void testCountWordLengths(){
        int[] counts = new int[31];
        FileResource file = new FileResource("manywords.txt");
        countWordLengths(file,counts);
        for (int k=0; k<counts.length; k++){
            System.out.println(counts[k]+" word of length "+k);

        }
            System.out.println(indexOfMax(counts));

    }
}
