# In a file called adieu.py, implement a program that prompts the user for names, one per line,
# until the user inputs control-d. Assume that the user will input at least one name.
# Then bid adieu to those names, separating two names with one and, three names with two commas and one and, and names with commas and one and
import inflect

list = []
p = inflect.engine()

while True:
    try:
        a = input("Name: ")
        list.append(a)

    except EOFError:
        mylist = p.join(list, final_sep=",")
        print("\nAdieu, adieu, to "+mylist)
        break
