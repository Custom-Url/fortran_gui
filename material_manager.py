import csv
import os
from PyQt6.QtWidgets import QDialog, QVBoxLayout, QLabel, QPushButton, QMessageBox

MATERIALS_FILE = "materials.csv"


def ensure_csv_exists():
    """Create the CSV file if it doesn't exist yet."""
    if not os.path.exists(MATERIALS_FILE):
        with open(MATERIALS_FILE, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["Name", "Type", "Lambda", "Mu"])
            writer.writeheader()


def load_materials():
    """Return a list of materials from the CSV."""
    ensure_csv_exists()
    with open(MATERIALS_FILE, "r", newline="") as f:
        reader = csv.DictReader(f)
        return list(reader)


def get_material(name):
    """Return a single material dict by name."""
    materials = load_materials()
    for m in materials:
        if m["Name"] == name:
            return m
    return None


def create_material(data):
    """Add a new material (dict) to the CSV."""
    ensure_csv_exists()
    materials = load_materials()
    if any(m["Name"] == data["Name"] for m in materials):
        raise ValueError(f"Material '{data['Name']}' already exists!")

    with open(MATERIALS_FILE, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["Name", "Type", "Lambda", "Mu"])
        writer.writerow(data)


def edit_material(updated_data):
    """Update an existing material by name."""
    ensure_csv_exists()
    materials = load_materials()
    for i, m in enumerate(materials):
        if m["Name"] == updated_data["Name"]:
            materials[i] = updated_data
            break

    with open(MATERIALS_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["Name", "Type", "Lambda", "Mu"])
        writer.writeheader()
        writer.writerows(materials)


def delete_material(name):
    """Delete a material from the CSV by name."""
    ensure_csv_exists()
    materials = [m for m in load_materials() if m["Name"] != name]

    with open(MATERIALS_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["Name", "Type", "Lambda", "Mu"])
        writer.writeheader()
        writer.writerows(materials)


# --- GUI Dialog (for later editing/creation/deletion in a window) ---
class MaterialManager(QDialog):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.setWindowTitle("Material Manager")
        self.resize(400, 300)

        layout = QVBoxLayout()
        layout.addWidget(QLabel("This is the Material Manager window.\n\n"
                                "Here you will be able to create, edit, and delete materials."))

        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        layout.addWidget(btn_close)

        self.setLayout(layout)

