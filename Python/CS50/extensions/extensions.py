#In a file called extensions.py, implement a program that prompts the user for the name of a file and then outputs that file’s media type if the file’s name ends, case-insensitively, in any of these suffixes:
#.gif
#.jpg
#.jpeg
#.png
#.pdf
#.txt
#.zip
#If the file’s name ends with some other suffix or has no suffix at all, output application/octet-stream instead, which is a common default.


def output_fcn (word):
    match word.lower().strip()[-3:]:
        case "gif":
            output = "image/"+word.lower().strip()[-3:]
        case "jpg":
            output = "image/"+"jpeg"
        case "png":
            output = "image/"+word.lower().strip()[-3:]
        case "pdf":
            output = "application/"+word.lower().strip()[-3:]
        case "txt":
            name = word.lower().strip()[0:-4]
            output = "text/"+name
        case "zip":
            output = "application/"+word.lower().strip()[-3:]
        case _:
            output = "application/octet-stream"

    if word.lower().strip()[-4:]== "jpeg":
        output = "image/"+word.lower().strip()[-4:]
    return output

name_and_extension = input("Enter a file name and extension ")
print(output_fcn(name_and_extension))

