import csv
import os

CSV_FILE = "materials.csv"
FIELDNAMES = ["Name", "Lambda", "Mu"]


def load_materials():
    """Load all materials from CSV into a list of dicts"""
    if not os.path.exists(CSV_FILE):
        return []
    with open(CSV_FILE, newline="") as f:
        reader = csv.DictReader(f)
        return list(reader)


def get_material(name):
    """Get a single material dict by name"""
    materials = load_materials()
    for mat in materials:
        if mat["Name"] == name:
            return mat
    return None


def create_material(new_data):
    """Add new material to CSV"""
    file_exists = os.path.exists(CSV_FILE)
    with open(CSV_FILE, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        if not file_exists:
            writer.writeheader()
        writer.writerow(new_data)


def edit_material(updated_data):
    """Edit an existing material in CSV"""
    materials = load_materials()
    for mat in materials:
        if mat["Name"] == updated_data["Name"]:
            mat.update(updated_data)
            break
    with open(CSV_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(materials)


def delete_material(name_to_delete):
    """Delete material from CSV by name"""
    materials = [m for m in load_materials() if m["Name"] != name_to_delete]
    with open(CSV_FILE, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(materials)

