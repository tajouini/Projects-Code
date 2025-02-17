# When texting or tweeting, it’s not uncommon to shorten words to save time or space,
# as by omitting vowels, much like Twitter was originally called twttr.
# In a file called twttr.py, implement a program that prompts the user for a str of text and
# then outputs that same text but with all vowels (A, E, I, O, and U) omitted, whether inputted in uppercase or lowercase.

def omit_vowels(str):
    str_without_vowels = ""
    for char in str:
        if char.capitalize() in ["A", "E", "I", "O", "U"]:
            str_without_vowels = str_without_vowels
        else:
            str_without_vowels = str_without_vowels + char

    return str_without_vowels


str_user = input ("Enter a string of text ")
print(omit_vowels(str_user))



