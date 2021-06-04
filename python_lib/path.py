import os
import sys

directory = os.getcwd()

while directory != "/" and os.path.basename(directory) != "app" and os.path.basename(directory) != "app_new":
	directory = os.path.dirname(directory)

sys.path.append(directory)
