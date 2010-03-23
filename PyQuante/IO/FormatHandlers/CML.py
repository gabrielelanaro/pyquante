from PyQuante.IO.Data import Data
import xml.dom.minidom as minidom

from PyQuante.Element import sym2no,symbol

class Handler(object):
    key = "cml"
    ext = ".cml"
    description = "cml - Chemical File Format an XML-like file format for chemical file"
    def read(self, string):
        doc = minidom.parseString(string)
        elements = doc.getElementsByTagName("atom")
        atomlist=[]
        for atom in elements:
            sym = atom.getAttribute("elementType")
            atno = sym2no[sym]
            x = float(atom.getAttribute("x3"))
            y = float(atom.getAttribute("y3"))
            z = float(atom.getAttribute("z3"))
            atomlist.append((atno,(x,y,z)))

        data = Data()
        data.build_molecule(atomlist = atomlist)
        return data
    def write(self, data):
        impl = minidom.getDOMImplementation()
        doc = impl.createDocument(None, "molecule", None)
        if data.has("molecule"):
            atomarray = doc.createElement("atomArray")
            doc.documentElement.appendChild(atomarray)
            for atom in data.molecule.atoms:
                atomx = doc.createElement("atom")
                atomx.setAttribute("elementType",symbol[atom.atno])
                atomx.setAttribute("x3",str(atom.r[0]))
                atomx.setAttribute("y3",str(atom.r[1]))
                atomx.setAttribute("z3",str(atom.r[2]))
                atomarray.appendChild(atomx)
        return doc.toprettyxml()
