file = "../_ser/bilininteg.py"

out = ["namespace mfem {"]
fid = open(file, 'r')
for line in fid.readlines():
    if line.startswith("class"):
        cname = (line.split(' ')[1]).split('(')[0]
         
    if line.startswith("    def __init__"):
        pp = ""
        if line.find("*args") != -1:
            pp = "    self._coeff = args"
        elif line.find("own_bfi_") != -1:
            pp = "    if own_bfi_ == 1:  bfi_.thisown = 0"
        elif line.find("integ, own_integ=1") != -1:
            pp = "    if own_integ == 1:  integ.thisown = 0"
        elif line.find("own_integs=1") != -1:
            pp = "    self.own_integs = own_integs"
        elif line.find(", vq)") != -1:
            pp = "    self._coeff = vq"
        elif line.find(", q)") != -1:
            pp = "    self._coeff = q"
        elif line.find(", q, a=1.0)") != -1:
            pp = "    self._coeff = q"
        elif line.find(", q, i") != -1:
            pp = "    self._coeff = q"         
        elif line.find(", sc)") != -1:
            pp = "    self._coeff = sc"
        elif line.find(", vc)") != -1:
            pp = "    self._coeff = vc"
        elif line.find("(self)") != -1: 
            pass
        elif line.find("(self, ir=None)") != -1: 
            pass
        elif line.find("(self, fes, e=1.0)") != -1:
            pass
        else:
            print(cname)
            print(line)
            assert False, "No recipt for this pattern "
        if pp != "":
            out.append("%pythonappend " + cname + "::" + cname + " %{")
            out.append(pp)
            out.append("%}")
            out.append("")            
fid.close()

out.append("%pythonappend SumIntegrator::AddIntegrator %{")
out.append("   if self.own_integs == 1: integ.thisown = 0")
out.append("%}")
out.append("}")            
out.append("")            

fid = open("bilininteg_ext.i", "w")
fid.write("\n".join(out))
fid.close()
