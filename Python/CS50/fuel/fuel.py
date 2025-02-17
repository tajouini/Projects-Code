# In a file called fuel.py, implement a program that prompts the user for a fraction, formatted as X/Y,
# wherein each of X and Y is an integer, and then outputs, as a percentage rounded to the nearest integer, how much fuel is in the tank.
# If, though, 1% or less remains, output E instead to indicate that the tank is essentially empty.
# And if 99% or more remains, output F instead to indicate that the tank is essentially full.
# If, though, X or Y is not an integer, X is greater than Y, or Y is 0, instead prompt the user again.
# (It is not necessary for Y to be 4.) Be sure to catch any exceptions like ValueError or ZeroDivisionError.


def tank_state():
    X = []
    Y = []
    while True:
        try:
            fraction_str = input("Fraction: ")
            for i in range(len(fraction_str)):
                if fraction_str[i] == "/":
                    X = int(fraction_str[0:i])
                    Y = int(fraction_str[i+1:])
                    break

            fraction = float(X/Y)

        except (ValueError, ZeroDivisionError, TypeError):
            pass

        else:
                if (X>Y):
                    pass
                else:
                    if fraction <= 0.01:
                        return "E"
                    elif fraction >= 0.99:
                        return "F"
                    else:
                        return str(round(fraction*100))+"%"

print(tank_state())
