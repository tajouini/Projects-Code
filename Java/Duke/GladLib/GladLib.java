import edu.duke.*;
import java.util.*;

public class GladLib {
    private ArrayList<String> adjectiveList;
    private ArrayList<String> nounList;
    private ArrayList<String> colorList;
    private ArrayList<String> countryList;
    private ArrayList<String> nameList;
    private ArrayList<String> animalList;
    private ArrayList<String> timeList;
    private ArrayList<String> verbList;
    private ArrayList<String> fruitList;
    private ArrayList<String> usedWordsList;
    private int numWords; 
    private Random myRandom;
    
    private static String dataSourceURL = "http://dukelearntoprogram.com/course3/data";
    private static String dataSourceDirectory = "data";
    
    public GladLib(){
        initializeFromSource(dataSourceDirectory);
        myRandom = new Random();
    }
    
    public GladLib(String source){
        initializeFromSource(source);
        myRandom = new Random();
    }
    
    private void initializeFromSource(String source) {
        adjectiveList= readIt(source+"/adjective.txt");    
        nounList = readIt(source+"/noun.txt");
        colorList = readIt(source+"/color.txt");
        countryList = readIt(source+"/country.txt");
        nameList = readIt(source+"/name.txt");        
        animalList = readIt(source+"/animal.txt");
        timeList = readIt(source+"/timeframe.txt");
        verbList = readIt(source+"/verb.txt");
        fruitList = readIt(source+"/fruit.txt");
        usedWordsList = new ArrayList<String>();
    }
    
    private String randomFrom(ArrayList<String> source){
        int index = myRandom.nextInt(source.size());
        return source.get(index);
    }
    
    private String getSubstitute(String label) {
       // usedWordsList.add(randomFrom(countryList));
        // usedWordsList.add(randomFrom(colorList));
        // usedWordsList.add(randomFrom(nounList));
        // usedWordsList.add(randomFrom(nameList));
        // usedWordsList.add(randomFrom(adjectiveList));
        // usedWordsList.add(randomFrom(animalList));
        // usedWordsList.add(randomFrom(timeList));
        // usedWordsList.add(randomFrom(verbList));
        // usedWordsList.add(randomFrom(fruitList));

        if (label.equals("country")) {
            usedWordsList.add(randomFrom(countryList));

            return randomFrom(countryList);
        }
        if (label.equals("color")){
            usedWordsList.add(randomFrom(colorList));

            return randomFrom(colorList);
        }
        if (label.equals("noun")){
            usedWordsList.add(randomFrom(nounList));

            return randomFrom(nounList);
        }
        if (label.equals("name")){
            usedWordsList.add(randomFrom(nameList));

            return randomFrom(nameList);
        }
        if (label.equals("adjective")){
            usedWordsList.add(randomFrom(adjectiveList));

            return randomFrom(adjectiveList);
        }
        if (label.equals("animal")){
            usedWordsList.add(randomFrom(animalList));

            return randomFrom(animalList);
        }
        if (label.equals("timeframe")){
            usedWordsList.add(randomFrom(timeList));

            return randomFrom(timeList);
        }
        if (label.equals("number")){
            return ""+myRandom.nextInt(50)+5;
        }
        if (label.equals("verb")){
            usedWordsList.add(randomFrom(verbList));

            return randomFrom(verbList);
        }
        if (label.equals("fruit")){
            usedWordsList.add(randomFrom(fruitList));

            return randomFrom(fruitList);
        }
        return "**UNKNOWN**";
    }
    
    private String processWord(String w){
        int first = w.indexOf("<");
        int last = w.indexOf(">",first);
        if (first == -1 || last == -1){
            return w;
        }
        numWords +=1;
        String prefix = w.substring(0,first);
        String suffix = w.substring(last+1);
        String sub = getSubstitute(w.substring(first+1,last));
        // move on until you choose a different word 
        for (int k = 0; k<usedWordsList.size(); k++){
            if (usedWordsList.get(k).equals(sub)){
            sub = getSubstitute(w.substring(first+1,last));
            }
        }
        return prefix+sub+suffix;
        
    }
    
    private void printOut(String s, int lineWidth){
        int charsWritten = 0;
        for(String w : s.split("\\s+")){
            if (charsWritten + w.length() > lineWidth){
                System.out.println();
                charsWritten = 0;
            }
            System.out.print(w+" ");
            charsWritten += w.length() + 1;
        }
    }
    
    private String fromTemplate(String source){
        String story = "";
        if (source.startsWith("http")) {
            URLResource resource = new URLResource(source);
            for(String word : resource.words()){
                story = story + processWord(word) + " ";
            }
        }
        else {
            FileResource resource = new FileResource(source);
            for(String word : resource.words()){
                story = story + processWord(word) + " ";
            }
        }
        return story;
    }
    
    private ArrayList<String> readIt(String source){
        ArrayList<String> list = new ArrayList<String>();
        if (source.startsWith("http")) {
            URLResource resource = new URLResource(source);
            for(String line : resource.lines()){
                list.add(line);
            }
        }
        else {
            FileResource resource = new FileResource(source);
            for(String line : resource.lines()){
                list.add(line);
            }
        }
        return list;
    }
    
    public void makeStory(){
        usedWordsList.clear();
        numWords = 0;
        System.out.println("\n");
        String story = fromTemplate("data/madtemplate2.txt");
        printOut(story, 60);
        System.out.println("\n");
        System.out.println("***We replaced "+ numWords+ " words***");
    }
    


}
