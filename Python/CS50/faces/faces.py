def convert(word):
    new_word = ""
    j = ""
    for i in word:
        if (i == ":"):
                new_word = new_word+""
        elif (j==":"):
            if i == ")":
                new_word = new_word + "🙂"
            elif i == "(":
                new_word = new_word + "🙁"
        else:
           new_word = new_word+i
        j = i
    return new_word

def main():
    word_user = input("Enter a word ")
    print(convert(word_user))


main()
