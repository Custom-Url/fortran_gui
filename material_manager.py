import csv
import os
from PyQt6.QtWidgets import (
    QDialog, QVBoxLayout, QLabel, QPushButton,
    QLineEdit, QComboBox, QFormLayout, QMessageBox
)

MATERIALS_FILE = "materials.csv"

FIELDNAMES = [
    "Name", "Type",
    "E_T", "G_TT",           # isotropic
    "E_L", "G_LT", "V_LT"    # transverse isotropic
]


def ensure_csv_exists():
    if not os.path.exists(MATERIALS_FILE):
        with open(MATERIALS_FILE, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
            writer.writeheader()


def load_materials():
    ensure_csv_exists()
    with open(MATERIALS_FILE, "r", newline="") as f:
        reader = csv.DictReader(f)
        return list(reader)


def get_material(name):
    materials = load_materials()
    for m in materials:
        if m["Name"] == name:
            return m
    return None


def save_all_materials(materials):
    with open(MATERIALS_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(materials)


def create_material(data):
    ensure_csv_exists()
    materials = load_materials()
    if any(m["Name"] == data["Name"] for m in materials):
        raise ValueError(f"Material '{data['Name']}' already exists!")
    materials.append(data)
    save_all_materials(materials)


def edit_material(data):
    materials = load_materials()
    for i, m in enumerate(materials):
        if m["Name"] == data["Name"]:
            materials[i] = data
            break
    save_all_materials(materials)


def delete_material(name):
    materials = [m for m in load_materials() if m["Name"] != name]
    save_all_materials(materials)


# -----------------------------
# Material Manager Dialog
# -----------------------------
class MaterialManager(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Material Manager")
        self.resize(500, 600)

        self.layout = QVBoxLayout()
        self.setLayout(self.layout)

        # Material selector
        self.material_combo = QComboBox()
        self.material_combo.addItem("Select material...")
        self.material_combo.currentTextChanged.connect(self.load_selected_material)
        self.layout.addWidget(QLabel("Select Material to Edit:"))
        self.layout.addWidget(self.material_combo)

        # Type selector
        self.type_combo = QComboBox()
        self.type_combo.addItems(["Isotropic", "Transversely Isotropic"])
        self.type_combo.currentTextChanged.connect(self.toggle_fields)
        self.layout.addWidget(QLabel("Material Type:"))
        self.layout.addWidget(self.type_combo)

        # Form fields
        self.fields = {}
        form_layout = QFormLayout()

        self.fields["Name"] = QLineEdit()
        self.fields["E_T"] = QLineEdit()
        self.fields["G_TT"] = QLineEdit()
        self.fields["E_L"] = QLineEdit()
        self.fields["G_LT"] = QLineEdit()
        self.fields["V_LT"] = QLineEdit()

        form_layout.addRow("Name:", self.fields["Name"])
        form_layout.addRow("Young's Modulus (Tranverse Plane):", self.fields["E_T"])
        form_layout.addRow("Shear Modulus (Tranverse Plane):", self.fields["G_TT"])
        form_layout.addRow("Young's Modulus (Longitudinal):", self.fields["E_L"])
        form_layout.addRow("Shear Modulus (Longitudinal-Tranverse):", self.fields["G_LT"])
        form_layout.addRow("Poissons Ratio:", self.fields["V_LT"])
        self.layout.addLayout(form_layout)

        # Connect signals
        self.fields["E_T"].textChanged.connect(self.autofill_iso)
        self.fields["G_TT"].textChanged.connect(self.autofill_iso)

        # Initialize
        self.update_fields(self.type_box.currentText())

        layout = QVBoxLayout()
        layout.addLayout(form_layout)
        self.setLayout(layout)

        # Buttons
        self.btn_create = QPushButton("Create New Material")
        self.btn_edit = QPushButton("Save Changes")
        self.btn_delete = QPushButton("Delete Material")
        self.btn_close = QPushButton("Close")

        self.layout.addWidget(self.btn_create)
        self.layout.addWidget(self.btn_edit)
        self.layout.addWidget(self.btn_delete)
        self.layout.addWidget(self.btn_close)

        # Connect buttons
        self.btn_create.clicked.connect(self.create_material_gui)
        self.btn_edit.clicked.connect(self.edit_material_gui)
        self.btn_delete.clicked.connect(self.delete_material_gui)
        self.btn_close.clicked.connect(self.close)

        # Load initial materials into combo
        self.refresh_material_list()
        self.toggle_fields()

    # -----------------------------
    # GUI Logic
    # -----------------------------
    def refresh_material_list(self):
        self.material_combo.blockSignals(True)
        self.material_combo.clear()
        self.material_combo.addItem("Select material...")
        for m in load_materials():
            self.material_combo.addItem(m["Name"])
        self.material_combo.blockSignals(False)

    def toggle_fields(self):
        """Grey out fields depending on type"""
        is_isotropic = (self.type_combo.currentText() == "Isotropic")

        # Grey out anisotropic fields if isotropic
        for key in ["E_L", "G_LT", "V_LT"]:
            self.fields[key].setDisabled(is_isotropic)

    def update_fields(self, text):
        """Enable/disable depending on isotropy"""
        is_iso = (text == "Isotropic")
        self.fields["E_L"].setDisabled(is_iso)
        self.fields["G_LT"].setDisabled(is_iso)
        self.fields["V_LT"].setDisabled(is_iso)

        if is_iso:
            self.autofill_iso()

    def autofill_iso(self):
        """Autofill isotropic fields when disabled"""
        if self.type_box.currentText() != "Isotropic":
            return

        try:
            E = float(self.fields["E_T"].text())
            G = float(self.fields["G_TT"].text())
            if E > 0 and G > 0:
                nu = E / (2 * G) - 1
                self.fields["E_L"].setText(f"{E:.4g}")
                self.fields["G_LT"].setText(f"{G:.4g}")
                self.fields["V_LT"].setText(f"{nu:.4g}")
        except ValueError:
            # ignore if inputs not valid yet
            pass
    def load_selected_material(self):
        name = self.material_combo.currentText()
        if name == "Select material...":
            for f in self.fields.values():
                f.clear()
            return

        mat = get_material(name)
        if mat:
            self.fields["Name"].setText(mat["Name"])
            self.type_combo.setCurrentText(mat["Type"])
            for key in ["E_T", "G_TT", "E_L", "G_LT", "V_LT"]:
                self.fields[key].setText(mat.get(key, ""))

    def collect_field_data(self):
        data = {}
        for key, widget in self.fields.items():
            data[key] = widget.text()
        data["Type"] = self.type_combo.currentText()
        return data

    def create_material_gui(self):
        data = self.collect_field_data()
        try:
            create_material(data)
            QMessageBox.information(self, "Material Created", f"Material '{data['Name']}' added.")
            self.refresh_material_list()
        except ValueError as e:
            QMessageBox.warning(self, "Error", str(e))

    def edit_material_gui(self):
        data = self.collect_field_data()
        edit_material(data)
        QMessageBox.information(self, "Material Updated", f"Material '{data['Name']}' updated.")
        self.refresh_material_list()

    def delete_material_gui(self):
        name = self.fields["Name"].text()
        if not name:
            return
        delete_material(name)
        QMessageBox.information(self, "Material Deleted", f"Material '{name}' deleted.")
        self.refresh_material_list()
        for f in self.fields.values():
            f.clear()
        self.material_combo.setCurrentIndex(0)

