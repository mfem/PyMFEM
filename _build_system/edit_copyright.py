import os

txt1 = "//\n"
txt2 = "// Copyright (c) 2020-2025, Princeton Plasma Physics Laboratory, All rights reserved.\n"

def find_i_file(path):

    for x in os.listdir(path):
        if not x.endswith(".i"):
            continue
               
        fname = os.path.join(path, x)
        print("processing: ", x)
        fid = open(fname)
        lines = fid.readlines()
        fid.close()
        
        if lines[1].startswith("// Copyright"):
            lines[1] = txt2
        else:
            lines = [txt1, txt2, txt1]+lines


        fid = open(fname, "w")
        fid.write("".join(lines))
        fid.close()
        
    

if __name__ == '__main__':
    find_i_file("../mfem/_ser")
    find_i_file("../mfem/_par")
    find_i_file("../mfem/common")        
    
