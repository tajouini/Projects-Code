#In a file called camel.py,
# implement a program that prompts the user for the name of a variable in camel case and outputs the corresponding name in snake case.
# Assume that the user’s input will indeed be in camel case.

def fromCamelToSnake(input):
    output = ""
    for char in input:
        if char != char.lower(): # character in upper case
            output = output+ "_"+char.lower()
        else:
            output = output + char
    return output

def main():
    word = input("Enter a word in Camel case ")
    print(fromCamelToSnake(word))

main()
