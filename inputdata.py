class InputData():
    def __init__(self,infile):
        self.f=infile.f
        self.nline=infile.nline
        self.problem_type=None
        self.element_type=None
        self.material_type=None


    def get_title(self):
        self.f.seek(0)
        for iline in range(self.nline):
            line=self.f.readline()
            if line.find('title')!=-1: break
        print('Title :',self.f.readline().strip())


    def get_problem(self):
        self.f.seek(0)
        for iline in range(self.nline):
            line=self.f.readline()
            if line.find('problem')!=-1: break
        self.problem_type=self.f.readline().strip()
        print('Problem type :',self.problem_type)


    def get_output_file_type(self):
        self.f.seek(0)
        for iline in range(self.nline):
            line=self.f.readline()
            io=line.find('output')
            if io !=-1: break
        noutfile=int(line.strip().split()[1])
        outfile_type=[]
        for ioutfile in range(noutfile):
            line=self.f.readline()
            outfile_type.append(line.strip())
        for ioutfile in range(noutfile):
            print('Type of output file',ioutfile+1,':',outfile_type[ioutfile])











if(__name__=="__main__"):
    openfile()
    output_data()
