import os, glob

def do_one_term( infilename, outfilename ):

    template = open("unimageit.com_template",mode='r').readlines()

    final_file = []

    for line in template:
        newline = line
        if line.find('INPUTSOURCEFILENAME')!=-1: newline = line.replace('INPUTSOURCEFILENAME',infilename)
        if line.find('OUTPUTCORRELATIONFILENAME')!=-1: newline = line.replace('OUTPUTCORRELATIONFILENAME',outfilename)
        final_file.append( newline )
    
    open("unimageit.com",mode='w').writelines(final_file)
    
    command = "./converts2c -histogram unimageit.com"
#    print command
    for line in os.popen(command).readlines(): print line.replace('\n','')
    
    
if __name__=="__main__":    
    
    for infilename in glob.glob("3d_cart_expanded_into_3d_sphr/source/term_*.coral"):
        outfilename = infilename.replace("/source/","/correlation/")
        do_one_term( infilename, outfilename )
#        break
