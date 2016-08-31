def read_mat(file):
    fid = open(file, 'r')
    r = []
    c = []
    realdata = []
    imagdata = []

    for l in fid.readlines():
        xx = l.split()
        r.append(int(xx[0]))
        c.append(int(xx[1]))
        realdata.append(float(xx[2]))
        if len(xx) > 3:
           imagdata.append(float(xx[3]))

    return r, c, realdata, imagdata
