file = "../_ser/lininteg.py"

out = ["namespace mfem {"]
fid = open(file, 'r')
for line in fid.readlines():
    if line.startswith("class"):
        cname = (line.split(' ')[1]).split('(')[0]

    if line.startswith("    def __init__"):
        pp = []
        if line.find(", ir=None") != -1:
            line = line.replace(", ir=None", "")
            pp.append("    self._ir=ir")
        if line.find("*args") != -1:
            pp.append("    self._coeff = args")
        elif line.find("self, seed_=0") != -1:
            pp.append("    self._coeff = QG")
        elif line.find("self, vqfc") != -1:
            pp.append("    self._coeff = vqfc")
        elif line.find("self, qfc") != -1:
            pp.append("    self._coeff = qfc")
        elif line.find(", QG") != -1:
            pp.append("    self._coeff = QG")
        elif line.find(", QF)") != -1:
            pp.append("    self._coeff = QF")
        elif line.find(", F)") != -1:
            pp.append("    self._coeff = F")
        elif line.find(", f)") != -1:
            pp.append("    self._coeff = f")
        elif line.find(", f, s=1.0)") != -1:
            pp.append("    self._coeff = f")
        elif line.find(", uD_, lambda_, mu_, alpha_, kappa_)") != -1:
            pp.append("    self._coeff = uD_")
        elif line.find("(self)") != -1:
            pass
        else:
            print(cname)
            print(line)
            assert False, "No recipt for this pattern "
        if len(pp) > 0:
            out.append("%pythonappend " + cname + "::" + cname + " %{")
            for x in pp:
                out.append(x)
            out.append("%}")
fid.close()
out.append("}")
fid = open("lininteg_ext.i", "w")
fid.write("\n".join(out))
fid.close()

