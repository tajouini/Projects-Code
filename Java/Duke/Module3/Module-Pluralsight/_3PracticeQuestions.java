 

public class _3PracticeQuestions {

    public static void main(String[] args) {

        // question 1
        System.out.println("Question 1");
        int x = 5;
        if(x < 5)
            System.out.println("a");
        x++;
        System.out.println("b");
        System.out.println(x);

        // question 2
        System.out.println("Question 2");
        int y = 2;
        if(y == 2) { // assignment vs. comparison
            System.out.println("a");
            y++;
        } else {
            System.out.println("b");
            System.out.println("c");
        }
        System.out.println(y);

        // question 3
        System.out.println("Question 3");
        int a = 3;
        int b = --a == 2 ? a++ : a--;
        System.out.println(a);
        System.out.println(b);
    }
}
