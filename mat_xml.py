import xml.etree.ElementTree as ET

def materials_to_xml(material_sections):
    """
    material_sections: list of tuples (combo, fields)
        combo: QComboBox with selected material
        fields: dict of QLineEdit fields keyed by "Name", "Type", "E_T", "G_TT", "E_L", "G_LT", "V_LT"
    Returns: XML string
    """
    root = ET.Element("Materials")

    # Reference material (hardcoded example, can be modified)
    ref = ET.SubElement(root, "Reference_Material", Lambda0="5e3", Mu0="20e3")

    for i, (combo, fields) in enumerate(material_sections, 1):
        name = combo.currentText()
        if name == "Select material...":
            continue

        law = fields["Type"].text()
        if law.lower().startswith("iso"):
            law_attr = "elasisoPF"
        else:
            law_attr = "elasanisoT"

        mat = ET.SubElement(root, "Material", numM=str(i),
                            Lib="beha/libUmatAmitex.so", Law=law_attr)

        # Read properties as floats
        def getf(key):
            try:
                return float(fields[key].text())
            except ValueError:
                return 0.0

        E_T = getf("E_T")
        G_TT = getf("G_TT")
        E_L = getf("E_L")
        G_LT = getf("G_LT")
        V_LT = getf("V_LT")

        if law_attr == "elasanisoT":
            # Compute 21 coefficients for transverse isotropic material
            # This is a common mapping, adjust if needed for your UMAT
            C = [0]*21
            # Example mapping (Voigt notation):
            C[0] = E_L/(1-V_LT**2*E_L/E_T)       # C11
            C[1] = E_T*V_LT/(1-V_LT**2*E_L/E_T)  # C12
            C[2] = C[1]                           # C13
            C[3] = E_T/(1)                        # C22 simplified
            C[4] = C[2]                           # C23
            C[5] = C[0]                           # C33
            C[6] = G_LT                            # C44
            C[7] = G_TT                            # C55
            C[8] = G_TT                            # C66
            C[9] = 0.1e-6                          # Thermal expansion
            C[10] = 22e-6
            C[11] = 22e-6
            C[12] = 0.0
            C[13] = 0.0
            C[14] = 1.0
            C[15] = 0.0
            C[16] = 1.0
            C[17] = 0.0
            C[18] = 0.0
            C[19] = 0.5e-3
            C[20] = 1.0
            # Add Coeff elements
            for idx, val in enumerate(C, 1):
                ET.SubElement(mat, "Coeff", Index=str(idx), Type="Constant", Value=str(val))

        elif law_attr == "elasisoPF":
            # Example: 6 coefficients for isotropic material
            C = [E_T, G_TT, V_LT, 0.0, 0.5e-3, 1.1e-2]
            for idx, val in enumerate(C, 1):
                ET.SubElement(mat, "Coeff", Index=str(idx), Type="Constant", Value=str(val))

        # Add IntVar elements (default to 0)
        for idx in range(1, 3):
            ET.SubElement(mat, "IntVar", Index=str(idx), Type="Constant", Value="0.")

    # Convert to pretty XML string
    def prettify(elem, level=0):
        indent = "    "
        txt = level * indent + f"<{elem.tag}"
        # Attributes
        if elem.attrib:
            txt += " " + " ".join(f'{k}="{v}"' for k, v in elem.attrib.items())
        children = list(elem)
        if children:
            txt += ">\n"
            for child in children:
                txt += prettify(child, level + 1) + "\n"
            txt += level * indent + f"</{elem.tag}>"
        else:
            txt += " />"
        return txt

    return '<?xml version="1.0" encoding="UTF-8"?>\n' + prettify(root)

