# Script for checking the score of the output file

if __name__ == '__main__'
    if (len(sys.argv) < 2):
        print ("Please provide 2 command line arguments")
        exit(0)

    inputfile = sys.argv[1]
    outputfile = sys.argv[2]
