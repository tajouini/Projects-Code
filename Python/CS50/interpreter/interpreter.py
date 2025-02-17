# x is an integer
# y is +, -, *, or /
# z is an integer
def calculate (word):
 # xyz = word.replace(" ", "") # remove all white spaces
   # word = word.strip()
    count = 0;  # count how many white spaces
    for i in range(len(word)):
        if word[i] == " " and count == 0:
            count = count + 1
            x = int(word[0:i])
            y = word[i+1]
            z = int(word[i+3:])

    match y:

        case "+":
            answer = x+z
        case "-":
            answer = x-z
        case "*":
            answer = x*z
        case "/":
            answer = x/z
    return answer

word = input("input x y z ")
#print(int(word[0:1]))
print("{:.1f}".format(calculate(word)))

