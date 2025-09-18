import xml.etree.ElementTree as ET

class MaterialManager:
    def __init__(self, filepath=None):
        self.filepath = filepath
        self.tree = None
        self.root = None
        if filepath:
            self.load_file(filepath)

    def load_file(self, filepath):
        """Load an XML material file."""
        self.tree = ET.parse(filepath)
        self.root = self.tree.getroot()
        self.filepath = filepath

    def save_file(self, filepath=None):
        """Save the XML material file."""
        if filepath is None:
            filepath = self.filepath
        if self.tree:
            self.tree.write(filepath, encoding="UTF-8", xml_declaration=True)

    def get_materials(self):
        """Return a list of material elements."""
        if self.root is not None:
            return self.root.findall("Material")
        return []

    def get_material_coeffs(self, material_index):
        """Get all Coeff elements for a material by index (1-based)."""
        materials = self.get_materials()
        if 0 < material_index <= len(materials):
            return materials[material_index-1].findall("Coeff")
        return []

    def set_coeff_value(self, material_index, coeff_index, value):
        """Set the value of a specific Coeff by indices (both 1-based)."""
        coeffs = self.get_material_coeffs(material_index)
        if 0 < coeff_index <= len(coeffs):
            coeffs[coeff_index-1].set("Value", str(value))

    def pretty_print(self):
        """Print XML tree nicely for debugging."""
        for mat in self.get_materials():
            print(f"Material numM={mat.get('numM')}")
            for coeff in mat.findall("Coeff"):
                print(f"  Coeff {coeff.get('Index')} = {coeff.get('Value')}")

