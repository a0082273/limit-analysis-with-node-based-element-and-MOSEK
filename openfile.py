class InputFile():
    def __init__(self):
        self.f=None
        self.nline=None



    def open_file(self):
        infilename=raw_input() #for python2.x
        # infilename=input() #for python3.x
        i=len(infilename)
        if infilename[i-3:i]!='dat':
            print('The last 3 letters of a input file is ',infilename[i-3:i])
            print('A input file name must end .dat')
            sys.exit()
        self.f=open(infilename,'r')




    def get_number_of_line(self):
        self.nline=sum(1 for line in self.f)
