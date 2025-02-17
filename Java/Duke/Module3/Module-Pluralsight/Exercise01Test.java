 

import org.junit.jupiter.api.DisplayName;
import org.junit.jupiter.api.Test;

import java.lang.reflect.Method;

import static org.junit.jupiter.api.Assertions.*;

class Exercise01Test {

    @Test
    @DisplayName("Exercise 1: Write a varargs method")
    void writeVarargsMethod() {
        // Note: The test cannot actually check if the second parameter of the method is a varargs parameter.
        Method concat = assertDoesNotThrow(() -> Exercise01.class.getDeclaredMethod("concat", String.class, Object[].class),
                "The class Exercise01 does not contain a method named 'concat' with the expected parameters.");

        // The method must return a String
        assertTrue(String.class.isAssignableFrom(concat.getReturnType()), "The return type of the method 'concat' is not String.");

        Exercise01 exercise = new Exercise01();

        assertEquals("", assertDoesNotThrow(() -> concat.invoke(exercise, ":", new Object[0]), "The 'concat' method threw an exception."));
        assertEquals("one", assertDoesNotThrow(() -> concat.invoke(exercise, ":", new Object[]{"one"}), "The 'concat' method threw an exception."));
        assertEquals("one:two", assertDoesNotThrow(() -> concat.invoke(exercise, ":", new Object[]{"one", "two"}), "The 'concat' method threw an exception."));
        assertEquals("one:two:three", assertDoesNotThrow(() -> concat.invoke(exercise, ":", new Object[]{"one", "two", "three"}), "The 'concat' method threw an exception."));
    }
}
