word = input("Enter a word ")
new = ""
for i in word:
 if i == " ":
    i = "..."
 new = new +i
print(new)
