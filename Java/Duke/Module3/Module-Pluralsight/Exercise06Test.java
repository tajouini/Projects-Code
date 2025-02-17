 

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;

class Exercise06Test {

    @Test
    @DisplayName("Exercise 6: Insert an element into an array")
    void insertElement() {
        String[] array = {"one", "two", "three"};

        String[] result1 = new Exercise06().insertElement(array, "four", 0);
        assertArrayEquals(new String[]{"four", "one", "two", "three"}, result1);
        assertArrayEquals(new String[]{"one", "two", "three"}, array, "The method modified the argument array.");

        String[] result2 = new Exercise06().insertElement(array, "four", 1);
        assertArrayEquals(new String[]{"one", "four", "two", "three"}, result2);
        assertArrayEquals(new String[]{"one", "two", "three"}, array, "The method modified the argument array.");

        String[] result3 = new Exercise06().insertElement(array, "four", 2);
        assertArrayEquals(new String[]{"one", "two", "four", "three"}, result3);
        assertArrayEquals(new String[]{"one", "two", "three"}, array, "The method modified the argument array.");

        String[] result4 = new Exercise06().insertElement(array, "four", 3);
        assertArrayEquals(new String[]{"one", "two", "three", "four"}, result4);
        assertArrayEquals(new String[]{"one", "two", "three"}, array, "The method modified the argument array.");
    }
}
