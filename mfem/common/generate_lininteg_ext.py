file = "../_ser/lininteg.py"

out = ["namespace mfem {"]
fid = open(file, 'r')
for line in fid.readlines():
    if line.startswith("class"):
        cname = (line.split(' ')[1]).split('(')[0]
         
    if line.startswith("    def __init__"):
        pp = ""
        if line.find("*args") != -1:
            pp = "    self._coeff = args"
        elif line.find(", QG") != -1:
            pp = "    self._coeff = QG"
        elif line.find(", QF)") != -1:
            pp = "    self._coeff = QF"
        elif line.find(", F)") != -1:
            pp = "    self._coeff = F"
        elif line.find(", f, s=1.0, ir=None)") != -1:
            pp = "    self._coeff = (f, ir)"
        elif line.find(", uD_, lambda_, mu_, alpha_, kappa_)") != -1:
            pp = "    self._coeff = uD_"
        elif line.find("(self)") != -1:
            pass
        else:
            print(cname)
            print(line)
            assert False, "No recipt for this pattern "
        if pp != "":
            out.append("%pythonappend " + cname + "::" + cname + " %{")
            out.append(pp)
            out.append("%}")
fid.close()
out.append("}")

fid = open("lininteg_ext.i", "w")
fid.write("\n".join(out))
fid.close()

