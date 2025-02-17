# In a file called grocery.py, implement a program that prompts the user for items,
# one per line, until the user inputs control-d (which is a common way of ending one’s input to a program).
# Then output the user’s grocery list in all uppercase, sorted alphabetically by item,
# prefixing each line with the number of times the user inputted that item. No need to pluralize the items.
# Treat the user’s input case-insensitively.

grocery_list = {}

while True:
    try:
            item = input("")
            if item.upper() in grocery_list:
                grocery_list[item.upper()] += 1
            elif not item.upper()== "" :
                grocery_list[item.upper()] = 1
    except EOFError:
        break

for key in sorted(grocery_list):
   # print("\n")
    print(str(grocery_list[key])+" "+key.upper())



