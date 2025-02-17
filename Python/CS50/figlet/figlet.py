from pyfiglet import Figlet
import sys

# from pyfiglet import Figlet
#f = pyfiglet.figlet_format("text to render", font="slant")
# print(f)
figlet = Figlet()



if (len(sys.argv) == 3):
    if sys.argv[1] not in ["-f", "--font"]:
        sys.exit("the first arg is not valid")
    if sys.argv[2] in figlet.getFonts():
        figlet.setFont(font=sys.argv[2])
        print(figlet.renderText(input("Input")))
    else:
        sys.exit("the font is not valid")
elif (len(sys.argv) == 1):
    print(figlet.renderText(input("Input")))
else:
    sys.exit("wrong number of args")


